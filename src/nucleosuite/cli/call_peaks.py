"""CLI registration for BigWig peak calling."""

from __future__ import annotations

from nucleosuite.cli.common import add_interval_output_arguments


def resolve_smooth_window(method: str, wps_input_mode: str, value: int | None) -> int:
    """Return the method-specific smoothing window used by call-peaks."""
    if value is not None:
        return value
    return 21 if method == "wps" and wps_input_mode == "raw" else 0


def register(subparsers):
    parser = subparsers.add_parser(
        "call-peaks",
        aliases=["peak-call", "peakcall"],
        help="Convert a continuous nucleosome-score or WPS BigWig into discrete peak positions.",
    )
    parser.add_argument("-i", "--input-bigwig", required=True, help="Nucleosome-score or WPS BigWig signal to segment.")
    parser.add_argument(
        "--blacklist-bed",
        help="BED blacklist; blacklisted signal is missing and overlapping calls are discarded.",
    )
    parser.add_argument("-o", "--out-prefix", dest="out_prefix", default=None, help="Output prefix; default: derived from the input BigWig and caller.")
    add_interval_output_arguments(parser)
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser)
    parser.add_argument(
        "--peak-caller", "--method", dest="method",
        choices=("pns", "wps"), required=True,
        help="Peak-calling algorithm: PNS positive-region caller or Kircher-style WPS caller.",
    )
    parser.add_argument(
        "--call-type", "--signal", dest="signal",
        choices=("both", "nucleosome", "breakpoint"), default="both",
        help="Write nucleosome calls, breakpoint calls, or both (default: both).",
    )
    parser.add_argument(
        "-r", "--regions", nargs="+", default=None,
        help="Optional contigs/regions; supports all, autosomes and numeric ranges.",
    )
    parser.add_argument("--chunk-bp", type=int, default=5_000_000, help="Core BigWig read chunk length in bp (default: 5000000).")
    parser.add_argument(
        "--overlap-bp", type=int, default=1_000,
        help="Signal padding around chunks to avoid boundary-split peaks (default: 1000 bp).",
    )
    parser.add_argument(
        "--smooth-window", type=int, default=None,
        help=(
            "Savitzky-Golay window. Defaults to 0 for the nucleosome-score caller and 21 for raw WPS; "
            "use 0 to disable smoothing explicitly."
        ),
    )
    parser.add_argument("--smooth-order", type=int, default=2, help="Savitzky-Golay polynomial order when smoothing is active (default: 2).")
    parser.add_argument(
        "--score-scale", type=float, default=1.0,
        help="Multiplier applied to called peak scores before BED output (default: 1).",
    )
    parser.add_argument(
        "--bigbed-score-scale", type=float, default=None,
        help=(
            "Multiplier applied to floating PNS BED peak scores during bigBed conversion. "
            "PNS defaults to 1 (no additional rescaling). WPS bigBed scores use their "
            "standard BED integer scale."
        ),
    )

    pns = parser.add_argument_group("PNS-style caller")
    pns.add_argument("--min-region-length", type=int, default=50, help="Minimum positive nucleosome-score region length retained as a call (default: 50 bp).")
    pns.add_argument(
        "--max-neg-run", type=int, default=0,
        help=(
            "Maximum consecutive zero-or-negative bases allowed inside a nucleosome-score "
            "positive region (default: 0)."
        ),
    )

    wps = parser.add_argument_group("WPS-style caller")
    wps.add_argument(
        "--wps-input-mode",
        choices=("raw", "adjusted"),
        default="raw",
        help=(
            "Interpret the input as raw WPS and apply the 1-kb median "
            "adjustment plus Savitzky-Golay smoothing, or as an already adjusted "
            "sm_mWPS signal to evaluate directly (default: raw)."
        ),
    )
    wps.add_argument(
        "--wps-baseline-window",
        type=int,
        default=1000,
        help="Running-median window used when --wps-input-mode raw (default: 1000 bp).",
    )
    wps.add_argument("--wps-min-length", type=int, default=50, help="Minimum WPS candidate or long-region subrun length (default: 50 bp).")
    wps.add_argument("--wps-max-length", type=int, default=150, help="Maximum selected WPS peak subrun length (default: 150 bp).")
    wps.add_argument("--wps-max-region", type=int, default=450, help="Maximum merged positive WPS region length evaluated (default: 450 bp).")
    wps.add_argument("--wps-merge-gap", type=int, default=5, help="Maximum coordinate difference joining successive positive WPS positions (default: 5 bp).")
    wps.add_argument("--wps-score-cutoff", type=float, default=5.0, help="Required strict lower bound for the maximum adjusted WPS value in a call (default: 5).")
    parser.set_defaults(command_runner=run)


def run(args):
    args.smooth_window = resolve_smooth_window(
        args.method, args.wps_input_mode, args.smooth_window
    )
    if args.bigbed_score_scale is None:
        args.bigbed_score_scale = 1.0
    if args.chunk_bp < 1 or args.overlap_bp < 0:
        raise ValueError("Chunk size must be positive and overlap non-negative")
    if args.smooth_window < 0 or (
        args.smooth_window and (args.smooth_window < 3 or args.smooth_window % 2 == 0)
    ):
        raise ValueError(
            "--smooth-window must be 0 or an odd integer of at least 3"
        )
    if args.smooth_order < 0:
        raise ValueError("--smooth-order must be non-negative")
    if args.bigbed_score_scale < 0:
        raise ValueError("--bigbed-score-scale must be 0 or greater")
    if args.min_region_length < 1 or args.max_neg_run < 0:
        raise ValueError("Invalid PNS peak-region settings")
    if args.wps_baseline_window < 1:
        raise ValueError("--wps-baseline-window must be positive")
    if args.wps_min_length < 1 or args.wps_max_length < args.wps_min_length:
        raise ValueError("Require 1 <= --wps-min-length <= --wps-max-length")
    if args.wps_max_region < args.wps_max_length:
        raise ValueError("--wps-max-region must be at least --wps-max-length")
    if args.method == "wps" and args.wps_input_mode == "raw":
        required_overlap = (
            args.wps_max_region
            + args.wps_baseline_window // 2
            + (args.smooth_window // 2 if args.smooth_window else 0)
        )
        if args.overlap_bp < required_overlap:
            raise ValueError(
                "Raw-WPS peak calling requires --overlap-bp >= "
                f"{required_overlap} for the selected peak and smoothing windows"
            )
    from nucleosuite.workflows import call_peaks as workflow
    from nucleosuite.parallel import run_bigwig_per_contig
    return run_bigwig_per_contig(
        "call-peaks", args, workflow.run,
        bigwig_attr="input_bigwig", selector_attr="regions", prefix_attr="out_prefix"
    )
