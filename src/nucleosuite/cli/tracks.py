"""CLI for the combined multi-range track generator."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from nucleosuite.cli.common import add_interval_output_arguments


TRACK_HELP = (
    "Track names may include pns, posPNS, pns_smoothed, wps, wps_smoothed, mWPS, sm_mWPS, coverage, dyad, fragment_ends, "
    "fragment_left_ends, fragment_right_ends, pns_peaks, wps_peaks, dinuc_profile, ww_types and type_dyads."
)


def register(subparsers):
    parser = subparsers.add_parser(
        "tracks",
        help="Generate multiple fragment-derived track sets in one input pass.",
        description=(
            "Read each fragment once per genomic chunk and update all requested "
            "nucleosome scores, WPS, coverage, dyad, fragment-end and sequence-derived outputs for multiple "
            "fragment-length ranges."
        ),
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-b", "--bamfiles", "--bam", dest="bamfiles", nargs="+",
        help="Input coordinate-sorted paired-end BAM file(s).",
    )
    inputs.add_argument(
        "--fragments", "--fragment-bed", dest="fragment_files", nargs="+",
        help="Input fragment BED, BED.gz, bigBed or .bb file(s).",
    )
    parser.add_argument("--fasta", help="Optional indexed reference FASTA.")
    parser.add_argument(
        "--chrom-sizes",
        help="Chromosome sizes, BAM or CRAM for fragment interval input.",
    )
    parser.add_argument(
        "--blacklist-bed",
        help=(
            "Optional BED; complete fragments overlapping excluded regions are "
            "removed before any requested track or sequence profile is calculated."
        ),
    )
    parser.add_argument(
        "-c", "--contigs", nargs="+", default=None,
        help="Contigs or intervals; supports lists, ranges, autosomes, all, and chrom:START-END selectors.",
    )
    parser.add_argument(
        "--fragment-range",
        action="append",
        default=[],
        metavar="RANGE=TRACKS",
        help=(
            "Repeated range specification, for example "
            "137-197=pns,posPNS,wps,coverage,pns_peaks,wps_peaks or "
            "145=dyad,fragment_left_ends,fragment_right_ends. " + TRACK_HELP
        ),
    )
    parser.add_argument(
        "--spec-file",
        action="append",
        default=[],
        help=(
            "Tab-separated specification file with fragment_range, output_prefix, "
            "tracks and optional basic_scope columns. May be repeated."
        ),
    )
    parser.add_argument(
        "--output-dir", default="tracks",
        help="Root directory for range-specific output trees (default: tracks).",
    )
    parser.add_argument(
        "--output-prefix", default=None,
        help="Prefix used for inline range specifications; default: derived from the input collection.",
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Optional TSV report written only after every requested output completes.",
    )
    parser.add_argument(
        "--max-duplicates",
        type=int,
        default=1,
        help="Maximum retained copies of an identical complete fragment; 0 disables (default: 1).",
    )
    parser.add_argument(
        "--max-per-coordinate",
        type=int,
        default=0,
        help=(
            "Optional cap on dyad and fragment-end signal at one coordinate; "
            "0 disables the cap (default: 0)."
        ),
    )
    parser.add_argument(
        "--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams",
        help="Apply identical-fragment coordinate limits across all inputs or within each input (default: all_bams).",
    )
    parser.add_argument("--chunk-bp", type=int, default=100_000, help="Core genomic chunk length in bp (default: 100000).")
    parser.add_argument("--overlap-bp", type=int, default=1_000, help="Context added to each side of a chunk before core-owned output is retained (default: 1000 bp).")
    parser.add_argument("--subsample", type=float, default=None, help="Optional independent fragment-retention probability from 0 to 1.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed used for reproducible fragment subsampling.")
    parser.add_argument(
        "--even-dyad", choices=("split", "left", "right"), default="split",
        help="Represent an even-length centre as 0.5 on both central bases or 1 on the left or right base (default: split).",
    )
    parser.add_argument(
        "--output-format", choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="Signal-track format for requested outputs (default: bigwig).",
    )
    parser.add_argument(
        "--staged-bedgraph-root",
        default=os.environ.get("NUCLEOSUITE_STAGED_BEDGRAPH_ROOT"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--staged-bedgraph-source-root",
        default=os.environ.get("NUCLEOSUITE_STAGED_BEDGRAPH_SOURCE_ROOT"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--staged-bedgraph-source-id",
        default=os.environ.get("NUCLEOSUITE_STAGED_BEDGRAPH_SOURCE_ID"),
        help=argparse.SUPPRESS,
    )
    add_interval_output_arguments(parser, default="bed")
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser, combine_resources=True, resumable=True)

    # Nucleosome-score settings shared by every score-enabled range.
    parser.add_argument("--score-mode-length", type=int, default=167, help="Modal protected-DNA length defining nucleosome-score geometry for every scoring range (default: 167 bp).")
    parser.add_argument("--score-smooth-window", type=int, default=0, help="Savitzky-Golay score smoothing window; 0 disables smoothing (default: 0).")
    parser.add_argument("--score-smooth-order", type=int, default=2, help="Score Savitzky-Golay polynomial order (default: 2).")
    parser.add_argument("--score-min-region-length", type=int, default=50, help="Minimum positive score region length retained as a nucleosome call (default: 50 bp).")
    parser.add_argument("--score-max-neg-run", type=int, default=0, help="Maximum zero-or-negative run bridged inside a score region (default: 0 bp).")
    parser.add_argument(
        "--score-peak-score-scale", type=float, default=1.0,
        help="Multiplier applied to score-derived peak scores written as six-decimal BED floats (default: 1).",
    )
    parser.add_argument(
        "--bigbed-score-scale", type=float, default=None,
        help=(
            "Multiplier applied to floating score-derived peak BED scores during bigBed "
            "conversion before integer rounding/clamping. PNS defaults to 1 (no score "
            "rescaling). WPS peak conversion is unchanged."
        ),
    )

    # WPS settings shared by every WPS-enabled range.
    parser.add_argument("--wps-protection", type=int, default=120, help="WPS protection-window width used for every WPS range (default: 120 bp).")
    parser.add_argument("--wps-baseline-window", type=int, default=1000, help="Running-median WPS baseline window (default: 1000 bp).")
    parser.add_argument("--wps-sg-window", type=int, default=21, help="WPS Savitzky-Golay window; 0 disables smoothing (default: 21 bp).")
    parser.add_argument("--wps-sg-order", type=int, default=2, help="WPS Savitzky-Golay polynomial order (default: 2).")
    parser.add_argument(
        "--wps-peak-track",
        choices=("wps", "wps_smoothed", "mWPS", "sm_mWPS"),
        default="sm_mWPS",
        help="WPS-family signal evaluated by the WPS peak caller (default: sm_mWPS).",
    )
    parser.add_argument("--wps-peak-score-scale", type=float, default=1.0, help="Multiplier used when WPS peak scores are converted to BED scores (default: 1).")
    parser.add_argument("--wps-peak-minlen", type=int, default=50, help="Minimum WPS candidate or long-region subrun length (default: 50 bp).")
    parser.add_argument("--wps-peak-maxlen", type=int, default=150, help="Maximum selected WPS peak subrun length (default: 150 bp).")
    parser.add_argument("--wps-peak-maxregion", type=int, default=450, help="Maximum merged positive WPS region length evaluated (default: 450 bp).")
    parser.add_argument("--wps-peak-merge-gap", type=int, default=5, help="Maximum coordinate difference joining successive positive WPS positions (default: 5 bp).")
    parser.add_argument("--wps-peak-varicutoff", type=float, default=5.0, help="Required strict lower bound for a called peak's maximum adjusted WPS value (default: 5).")
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    parser.set_defaults(command_runner=run)
    return parser


def _validate(args) -> None:
    if bool(args.bamfiles) == bool(args.fragment_files):
        raise ValueError("Provide exactly one of --bam or --fragments")
    if not args.fragment_range and not args.spec_file:
        raise ValueError("Provide at least one --fragment-range or --spec-file")
    if args.max_duplicates < 0 or args.max_per_coordinate < 0:
        raise ValueError("Duplicate and coordinate limits must be 0 or greater")
    if args.chunk_bp < 1 or args.overlap_bp < 0:
        raise ValueError("Chunk size must be positive and overlap non-negative")
    if args.subsample is not None and not 0.0 <= args.subsample <= 1.0:
        raise ValueError("--subsample must be between 0 and 1")
    if args.score_mode_length < 3:
        raise ValueError("--score-mode-length must be at least 3")
    if args.score_smooth_window < 0 or (
        args.score_smooth_window and (
            args.score_smooth_window < 3 or args.score_smooth_window % 2 == 0
        )
    ):
        raise ValueError("--score-smooth-window must be 0 or an odd integer >= 3")
    if args.score_smooth_order < 0 or (
        args.score_smooth_window and args.score_smooth_order >= args.score_smooth_window
    ):
        raise ValueError("Invalid score smoothing order")
    if args.score_min_region_length < 1 or args.score_max_neg_run < 0:
        raise ValueError("Invalid score peak-region settings")
    if args.bigbed_score_scale < 0:
        raise ValueError("--bigbed-score-scale must be 0 or greater")
    if args.wps_protection < 2 or args.wps_baseline_window < 1:
        raise ValueError("Invalid WPS protection or baseline window")
    if args.wps_sg_window < 0 or (
        args.wps_sg_window and args.wps_sg_window % 2 == 0
    ):
        raise ValueError("--wps-sg-window must be odd, or 0 to disable")
    if args.wps_sg_order < 0 or (
        args.wps_sg_window and args.wps_sg_order >= args.wps_sg_window
    ):
        raise ValueError("Invalid WPS smoothing order")
    if args.wps_peak_minlen < 1 or args.wps_peak_maxlen < args.wps_peak_minlen:
        raise ValueError("Invalid WPS peak lengths")
    if args.wps_peak_maxregion < args.wps_peak_maxlen:
        raise ValueError("--wps-peak-maxregion must be at least --wps-peak-maxlen")
    staging_values = (
        args.staged_bedgraph_root,
        args.staged_bedgraph_source_root,
        args.staged_bedgraph_source_id,
    )
    if any(staging_values) and not all(staging_values):
        raise ValueError(
            "Staged bedGraph output requires root, source root and source identifier"
        )
    for path in args.spec_file:
        if not Path(path).is_file():
            raise FileNotFoundError(path)


def run(args):
    if args.bigbed_score_scale is None:
        args.bigbed_score_scale = 1.0
    args.scoring_method = "pns"
    _validate(args)
    from nucleosuite.workflows import tracks as workflow
    from nucleosuite.parallel import run_tracks_per_contig

    return run_tracks_per_contig(args, workflow.run)
