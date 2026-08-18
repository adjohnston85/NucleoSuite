"""CLI registration for the full PNS workflow."""

from __future__ import annotations

from nucleosuite.cli.common import (
    add_bam_fragment_arguments,
    add_randomization_arguments,
    add_interval_output_arguments,
    normalise_track_list,
    validate_bam_arguments,
)
from nucleosuite.scoring.basic_tracks import BASIC_TRACKS
from nucleosuite.scoring.pns import BNS_TRACKS, PNS_TRACKS, TNS_TRACKS, SCORING_METHODS


def register(subparsers):
    parser = subparsers.add_parser(
        "pns",
        help="Build a nucleosome-position signal from paired-end fragment geometry.",
        description=(
            "Build probabilistic nucleosome score (PNS), boxcar nucleosome score "
            "(BNS), or triangular nucleosome score (TNS) kernels, optional auxiliary "
            "tracks, and shared peak calls."
        ),
    )
    add_bam_fragment_arguments(
        parser, default_lower=137, default_upper=197,
        default_max_duplicates=1, default_even_dyad="split",
    )
    add_randomization_arguments(parser)
    add_interval_output_arguments(parser)
    parser.add_argument(
        "--scoring-method", choices=SCORING_METHODS, default="pns",
        help=(
            "Nucleosome scoring kernel: pns uses endpoint-derived probability "
            "triangles; bns uses a centred balanced boxcar; tns uses a centred "
            "unit-mass triangle spanning the scoring support (default: pns)."
        ),
    )
    parser.add_argument(
        "--mode-length", type=int, default=167,
        help=(
            "Modal protected-DNA length that defines scoring support for accepted "
            "fragment lengths (default: 167 bp)."
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
        "--pns-mode", choices=("on", "off"), default="on",
        help="Enable or disable nucleosome-score track and peak generation (default: on).",
    )
    parser.add_argument(
        "--pns-format", choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="File format for nucleosome-score signal tracks (default: bigwig).",
    )
    parser.add_argument(
        "--other-format", choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="File format for coverage, dyad, and fragment-end tracks (default: bigwig).",
    )
    parser.add_argument(
        "--score-tracks", "--pns-tracks", dest="pns_tracks", nargs="*",
        default=["pns", "posPNS"],
        help=(
            "Signal tracks for the selected scoring method. PNS uses pns, posPNS, "
            "and optional pns_smoothed; BNS uses bns, posBNS, and optional "
            "bns_smoothed; TNS uses tns, posTNS, and optional tns_smoothed. The "
            "default writes the raw centred and uncentred tracks."
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
        "--bigbed-score-scale", type=float, default=1000.0,
        help=(
            "Multiplier applied to floating BED peak scores when converting "
            "the bigBed score field to an integer in the 0-1000 range "
            "(default: 1000)."
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
    parser.set_defaults(command_runner=run)


def run(args):
    validate_bam_arguments(args)
    if args.mode_length < 3:
        raise ValueError("--mode-length must be at least 3")
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
    if args.bigbed_score_scale < 0:
        raise ValueError("--bigbed-score-scale must be 0 or greater")
    if args.peak_coverage_threshold is not None:
        import math
        if not math.isfinite(args.peak_coverage_threshold) or args.peak_coverage_threshold < 0:
            raise ValueError(
                "--peak-coverage-threshold must be a finite value of 0 or greater"
            )
        if args.pns_mode == "off":
            raise ValueError("--peak-coverage-threshold requires --pns-mode on")
    all_score_tracks = set(PNS_TRACKS) | set(BNS_TRACKS) | set(TNS_TRACKS)
    args.pns_tracks = normalise_track_list(
        args.pns_tracks, all_score_tracks, "--score-tracks"
    )
    track_sets = {
        "pns": PNS_TRACKS,
        "bns": BNS_TRACKS,
        "tns": TNS_TRACKS,
    }
    default_track_map = {
        "bns": {
            "pns": "bns",
            "posPNS": "posBNS",
            "pns_smoothed": "bns_smoothed",
        },
        "tns": {
            "pns": "tns",
            "posPNS": "posTNS",
            "pns_smoothed": "tns_smoothed",
        },
    }
    if args.scoring_method != "pns":
        mapping = default_track_map[args.scoring_method]
        args.pns_tracks = [mapping.get(track, track) for track in args.pns_tracks]
    allowed_tracks = track_sets[args.scoring_method]
    invalid = [track for track in args.pns_tracks if track not in allowed_tracks]
    if invalid:
        method = args.scoring_method.upper()
        raise ValueError(
            f"{method} scoring accepts {', '.join(allowed_tracks)} score tracks"
        )
    args.other_tracks = normalise_track_list(
        args.other_tracks, set(BASIC_TRACKS), "--other-tracks"
    )
    if args.pns_mode == "off":
        args.pns_tracks = []
        args.pns_format = "none"
    elif any(track.endswith("_smoothed") for track in args.pns_tracks) and args.smooth_window == 0:
        raise ValueError(
            "A smoothed score track requires --smooth-window to enable smoothing"
        )
    if args.dinuc_fraction and not args.dinuc_profile:
        raise ValueError("--dinuc-fraction requires --dinuc-profile")
    if (args.dinuc_profile or args.split_ww_types) and not args.fasta:
        raise ValueError("--dinuc-profile and --split-ww-types require --fasta")
    from nucleosuite.workflows import pns as workflow
    from nucleosuite.parallel import run_native_per_contig
    return run_native_per_contig("pns", args, workflow.run)
