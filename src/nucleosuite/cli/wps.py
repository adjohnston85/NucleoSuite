"""CLI registration for the full WPS workflow."""

from __future__ import annotations

from nucleosuite.cli.common import (
    add_bam_fragment_arguments,
    add_randomization_arguments,
    add_interval_output_arguments,
    normalise_track_list,
    validate_bam_arguments,
)
from nucleosuite.scoring.basic_tracks import BASIC_TRACKS
from nucleosuite.scoring.wps import WPS_TRACKS


def register(subparsers):
    parser = subparsers.add_parser(
        "wps",
        help="Build a window-protection signal showing protected and exposed DNA.",
        description=(
            "Calculate raw WPS and its smoothed and running-median-adjusted tracks, "
            "then apply the selected PNS- or WPS-style peak caller."
        ),
    )
    add_bam_fragment_arguments(
        parser, default_lower=120, default_upper=180,
        default_max_duplicates=0, default_even_dyad="left",
    )
    add_randomization_arguments(parser)
    add_interval_output_arguments(parser)
    parser.add_argument(
        "--protection", type=int, default=120,
        help="Centred WPS protection-window width (default: 120 bp).",
    )
    parser.add_argument(
        "--baseline-window", type=int, default=1000,
        help="Running-median window subtracted from raw and smoothed WPS tracks (default: 1000 bp).",
    )
    parser.add_argument(
        "--sg-window", type=int, default=21,
        help="Savitzky-Golay smoothing window; use 0 for no smoothing (default: 21 bp).",
    )
    parser.add_argument(
        "--sg-order", type=int, default=2,
        help="Savitzky-Golay polynomial order (default: 2).",
    )
    parser.add_argument(
        "--score-format", choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="File format for requested WPS and auxiliary signal tracks (default: bigwig).",
    )
    parser.add_argument(
        "--score-tracks", nargs="*",
        default=["coverage", "sm_mWPS", "wps", "wps_smoothed", "mWPS", "dyad"],
        help=(
            "Tracks to write: wps, wps_smoothed, mWPS, sm_mWPS, coverage, dyad, "
            "fragment_ends, fragment_left_ends, fragment_right_ends, or none. "
            "Default: coverage, sm_mWPS, wps, wps_smoothed, mWPS, dyad."
        ),
    )
    parser.add_argument(
        "--peak-caller", choices=("wps", "pns", "none"), default="wps",
        help="Peak segmentation applied to --peak-track, or none to omit peak calls (default: wps).",
    )
    parser.add_argument(
        "--peak-track", choices=WPS_TRACKS, default="sm_mWPS",
        help="WPS-derived signal on which peaks are called (default: sm_mWPS).",
    )
    parser.add_argument(
        "--peak-score-scale", type=float, default=1.0,
        help="Multiplier applied before peak scores are rounded into the BED score field (default: 1).",
    )
    parser.add_argument("--peak-minlen", type=int, default=50, help="Minimum WPS candidate or long-region subrun length in bp (default: 50).")
    parser.add_argument("--peak-maxlen", type=int, default=150, help="Maximum selected WPS peak subrun length in bp (default: 150).")
    parser.add_argument("--peak-maxregion", type=int, default=450, help="Maximum merged positive WPS region length evaluated for peaks (default: 450 bp).")
    parser.add_argument("--peak-merge-gap", type=int, default=5, help="Maximum coordinate difference joining successive positive WPS positions (default: 5 bp).")
    parser.add_argument("--peak-varicutoff", type=float, default=5.0, help="Required strict lower bound for the maximum adjusted WPS value in a called peak (default: 5).")
    parser.set_defaults(command_runner=run)


def run(args):
    validate_bam_arguments(args)
    if args.protection < 2:
        raise ValueError("--protection must be at least 2")
    if args.baseline_window < 1:
        raise ValueError("--baseline-window must be positive")
    if args.sg_window < 0 or (args.sg_window and args.sg_window % 2 == 0):
        raise ValueError("--sg-window must be odd, or 0 to disable smoothing")
    valid_tracks = set(WPS_TRACKS) | set(BASIC_TRACKS)
    args.score_tracks = normalise_track_list(
        args.score_tracks, valid_tracks, "--score-tracks"
    )
    if args.peak_minlen < 1 or args.peak_maxlen < args.peak_minlen:
        raise ValueError("Require 1 <= --peak-minlen <= --peak-maxlen")
    if args.peak_maxregion < args.peak_maxlen:
        raise ValueError("--peak-maxregion must be at least --peak-maxlen")
    if args.peak_caller == "wps" and args.peak_track == "sm_mWPS":
        required_overlap = (
            args.peak_maxregion
            + args.baseline_window // 2
            + (args.sg_window // 2 if args.sg_window else 0)
        )
        if args.overlap_bp < required_overlap:
            raise ValueError(
                "WPS peak calling requires --overlap-bp >= "
                f"{required_overlap} for the selected peak and preprocessing windows"
            )
    from nucleosuite.workflows import wps as workflow
    from nucleosuite.parallel import run_native_per_contig
    return run_native_per_contig("wps", args, workflow.run)
