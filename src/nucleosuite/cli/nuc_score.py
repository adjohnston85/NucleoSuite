"""CLI registration for generic nucleosome scoring."""

from __future__ import annotations

import argparse

from nucleosuite.cli.common import (
    add_mode_estimation_arguments,
    add_bam_fragment_arguments,
    add_randomization_arguments,
    add_interval_output_arguments,
    normalise_track_list,
    auto_or_integer_mode,
    resolve_fragment_mode,
    validate_bam_arguments,
)
from nucleosuite.scoring.basic_tracks import BASIC_TRACKS
from nucleosuite.scoring.pns import BNS_TRACKS, PNS_TRACKS, SNS_TRACKS, TNS_TRACKS, SCORING_METHODS


def register(subparsers):
    parser = subparsers.add_parser(
        "nuc-score",
        help="Build nucleosome-position score tracks from paired-end fragment geometry.",
        description=(
            "Build sinusoidal nucleosome score (SNS), probabilistic nucleosome score "
            "(PNS), boxcar nucleosome score (BNS), or triangular nucleosome score "
            "(TNS) kernels, optional auxiliary tracks, and shared peak calls."
        ),
    )
    add_bam_fragment_arguments(
        parser, default_lower=None, default_upper=None,
        default_max_duplicates=1, default_even_dyad="split",
    )
    add_randomization_arguments(parser)
    add_interval_output_arguments(parser)
    parser.add_argument(
        "--scoring-method", choices=SCORING_METHODS, default="sns",
        help=(
            "Nucleosome scoring kernel: sns uses the length-adaptive sinusoidal "
            "kernel; pns uses endpoint-derived probability triangles; bns uses a "
            "centred balanced boxcar; tns uses a centred unit-mass triangle "
            "(default: sns)."
        ),
    )
    parser.add_argument(
        "--mode", "--mode-length", dest="mode_length",
        type=auto_or_integer_mode, default="auto",
        help=(
            "Modal protected-DNA length defining the scoring geometry. auto "
            "estimates it from accepted fragments; an integer bypasses estimation "
            "(default: auto)."
        ),
    )
    add_mode_estimation_arguments(
        parser, default_search_lower=137, default_search_upper=197
    )
    parser.add_argument(
        "--frag-mode-padding",
        type=int,
        default=30,
        help=(
            "When --frag-lower and/or --frag-upper are omitted, derive each missing "
            "scoring bound from the resolved mode plus or minus this many bp "
            "(default: 30). Explicit fragment bounds override the corresponding "
            "automatic bound independently."
        ),
    )
    parser.add_argument(
        "--smooth-window", type=int, default=0,
        help="Savitzky-Golay window for optional score smoothing; 0 disables (default: 0).",
    )
    parser.add_argument(
        "--smooth-order", type=int, default=2,
        help="Savitzky-Golay polynomial order used when score smoothing is enabled (default: 2).",
    )
    parser.add_argument(
        "--dinuc-profile", action="store_true",
        help="Write observed dyad-aligned 16-dinucleotide profiles.",
    )
    parser.add_argument(
        "--dinuc-fraction", action="store_true",
        help="Report dinucleotide frequencies as fractions from 0 to 1 instead of percentages.",
    )
    parser.add_argument(
        "--split-ww-types", action="store_true",
        help="Write all-fragment and type1-type4 score, track, and peak outputs.",
    )
    parser.add_argument(
        "--score-mode", dest="score_mode", choices=("on", "off"), default="on",
        help="Enable or disable nucleosome-score track and peak generation (default: on).",
    )
    parser.add_argument(
        "--score-format", dest="score_format", choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="File format for nucleosome-score signal tracks (default: bigwig).",
    )
    parser.add_argument(
        "--other-format", choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="File format for coverage, dyad, and fragment-end tracks (default: bigwig).",
    )
    parser.add_argument(
        "--score-tracks", dest="score_tracks", nargs="*",
        default=["sns", "posSNS"],
        help=(
            "Signal tracks for the selected scoring method. SNS uses sns, posSNS, "
            "and optional sns_smoothed; PNS uses pns, posPNS, and optional "
            "pns_smoothed; BNS uses bns, posBNS, and optional bns_smoothed; TNS "
            "uses tns, posTNS, and optional tns_smoothed. The default writes the "
            "raw score and matching non-negative reference tracks."
        ),
    )
    parser.add_argument(
        "--other-tracks", nargs="*", default=list(BASIC_TRACKS),
        help=(
            "Other tracks: coverage dyad fragment_ends fragment_left_ends "
            "fragment_right_ends, or none (default: all five)."
        ),
    )
    parser.add_argument(
        "--peak-score-scale", type=float, default=1.0,
        help=(
            "Multiplier applied to nucleosome-score peak values written to BED as six-decimal "
            "floats (default: 1)."
        ),
    )
    parser.add_argument(
        "--bigbed-score-scale", type=float, default=None,
        help=(
            "Multiplier applied to floating BED peak scores when converting "
            "the bigBed score field to an integer in the 0-1000 range. "
            "By default SNS is not rescaled (1); PNS, BNS and TNS retain the "
            "1000-fold conversion used for their fractional score ranges."
        ),
    )
    parser.add_argument(
        "--min-region-length", type=int, default=50,
        help="Minimum retained positive score-region length in bp (default: 50).",
    )
    parser.add_argument(
        "--max-neg-run", type=int, default=0,
        help=(
            "Maximum consecutive zero-or-negative bases allowed inside a positive "
            "score region (default: 0)."
        ),
    )
    parser.add_argument(
        "--peak-coverage-threshold",
        type=float,
        help=(
            "Optional minimum fragment coverage at the nucleosome peak position "
            "written in BED column 7. Peaks below the threshold are not written; "
            "breakpoint peaks are unchanged (default: off)."
        ),
    )
    parser.add_argument(
        "--peak-calling",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Call nucleosome and breakpoint intervals while producing score tracks "
            "(default: enabled; --no-peak-calling writes tracks without either callset)."
        ),
    )
    parser.set_defaults(command_runner=run)


def _resolve_mode_and_fragment_range(args):
    """Resolve the scoring mode, then fill omitted fragment bounds from mode +/- padding."""

    import copy

    if args.frag_mode_padding < 0:
        raise ValueError("--frag-mode-padding must be non-negative")
    if args.frag_lower is not None and args.frag_lower < 1:
        raise ValueError("--frag-lower must be at least 1")
    if args.frag_upper is not None and args.frag_upper < 1:
        raise ValueError("--frag-upper must be at least 1")

    if args.mode_length == "auto":
        # Mode estimation is deliberately decoupled from the eventual scoring bounds.
        # This lets the mode be estimated from a stable search interval first, after
        # which the default scoring interval follows the observed mode.
        estimate_args = copy.copy(args)
        estimate_args.frag_lower = int(args.mode_search_lower)
        estimate_args.frag_upper = int(args.mode_search_upper)
        mode, estimate, mode_source, mode_seed = resolve_fragment_mode(
            estimate_args, args.mode_length, command="nuc-score"
        )
    else:
        mode, estimate, mode_source, mode_seed = resolve_fragment_mode(
            args, args.mode_length, command="nuc-score"
        )

    args.mode_length = int(mode)
    if args.frag_lower is None:
        args.frag_lower = max(1, args.mode_length - args.frag_mode_padding)
    if args.frag_upper is None:
        args.frag_upper = args.mode_length + args.frag_mode_padding
    return args.mode_length, estimate, mode_source, mode_seed


def run(args):
    mode, estimate, mode_source, mode_seed = _resolve_mode_and_fragment_range(args)
    validate_bam_arguments(args)
    if args.smooth_window < 0:
        raise ValueError("--smooth-window must be 0 or greater")
    if args.smooth_window and (args.smooth_window < 3 or args.smooth_window % 2 == 0):
        raise ValueError(
            "--smooth-window must be 0 or an odd integer of at least 3"
        )
    if args.smooth_order < 0:
        raise ValueError("--smooth-order must be non-negative")
    if args.min_region_length < 1 or args.max_neg_run < 0:
        raise ValueError("Peak-region lengths must be valid positive values")
    if args.bigbed_score_scale is None:
        args.bigbed_score_scale = 1.0 if args.scoring_method == "sns" else 1000.0
    if args.bigbed_score_scale < 0:
        raise ValueError("--bigbed-score-scale must be 0 or greater")
    if args.peak_coverage_threshold is not None:
        import math
        if not math.isfinite(args.peak_coverage_threshold) or args.peak_coverage_threshold < 0:
            raise ValueError(
                "--peak-coverage-threshold must be a finite value of 0 or greater"
            )
        if args.score_mode == "off":
            raise ValueError("--peak-coverage-threshold requires --score-mode on")
        if not args.peak_calling:
            raise ValueError("--peak-coverage-threshold requires peak calling")
    all_score_tracks = set(PNS_TRACKS) | set(SNS_TRACKS) | set(BNS_TRACKS) | set(TNS_TRACKS)
    args.score_tracks = normalise_track_list(
        args.score_tracks, all_score_tracks, "--score-tracks"
    )
    track_sets = {
        "sns": SNS_TRACKS,
        "pns": PNS_TRACKS,
        "bns": BNS_TRACKS,
        "tns": TNS_TRACKS,
    }
    default_track_map = {
        "pns": {
            "sns": "pns",
            "posSNS": "posPNS",
            "sns_smoothed": "pns_smoothed",
        },
        "bns": {
            "sns": "bns",
            "posSNS": "posBNS",
            "sns_smoothed": "bns_smoothed",
        },
        "tns": {
            "sns": "tns",
            "posSNS": "posTNS",
            "sns_smoothed": "tns_smoothed",
        },
    }
    if args.scoring_method != "sns":
        mapping = default_track_map[args.scoring_method]
        args.score_tracks = [mapping.get(track, track) for track in args.score_tracks]
    allowed_tracks = track_sets[args.scoring_method]
    invalid = [track for track in args.score_tracks if track not in allowed_tracks]
    if invalid:
        method = args.scoring_method.upper()
        raise ValueError(
            f"{method} scoring accepts {', '.join(allowed_tracks)} score tracks"
        )
    args.other_tracks = normalise_track_list(
        args.other_tracks, set(BASIC_TRACKS), "--other-tracks"
    )
    if args.score_mode == "off":
        args.score_tracks = []
        args.score_format = "none"
    elif any(track.endswith("_smoothed") for track in args.score_tracks) and args.smooth_window == 0:
        raise ValueError(
            "A smoothed score track requires --smooth-window to enable smoothing"
        )
    if args.dinuc_fraction and not args.dinuc_profile:
        raise ValueError("--dinuc-fraction requires --dinuc-profile")
    if (args.dinuc_profile or args.split_ww_types) and not args.fasta:
        raise ValueError("--dinuc-profile and --split-ww-types require --fasta")
    if not getattr(args, "_per_contig_worker", False):
        from nucleosuite.mode_estimation import write_single_mode_report
        from nucleosuite.workflows.common import default_output_prefix, input_paths_from_args

        requested_prefix = args.out_prefix or default_output_prefix(
            input_paths_from_args(args), args.contigs
        )
        report_prefix = (
            f"{requested_prefix}_method{args.scoring_method}_mode{args.mode_length}"
            f"_lower{args.frag_lower}_upper{args.frag_upper}"
            f"_smooth{args.smooth_window}x{args.smooth_order}"
        )
        report = write_single_mode_report(
            f"{report_prefix}_fragment_mode_estimation.tsv",
            estimate=estimate,
            resolved_mode=args.mode_length,
            mode_source=mode_source,
            seed=mode_seed,
        )
        print(f"[nuc-score] Fragment mode report: {report}", flush=True)
    from nucleosuite.workflows import pns as workflow
    from nucleosuite.parallel import run_native_per_contig
    return run_native_per_contig("nuc-score", args, workflow.run)
