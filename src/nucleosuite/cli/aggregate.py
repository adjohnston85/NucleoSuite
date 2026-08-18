"""Command-line interface for ``nucleosuite aggregate``."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.align import AlignmentConfig, run_alignment
from nucleosuite.progress import ProgressReporter


def add_aggregate_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "aggregate",
        help="Aggregate per-base BigWig signal around regions or relative nucleosomes.",
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.set_defaults(command_function=run_from_args, command_runner=run_from_args)
    parser.add_argument("-w", "--bigwig", type=Path, required=True, help="BigWig signal to extract around each accepted feature.")
    parser.add_argument("-r", "--region-bed", dest="region_bed", type=Path, required=True, help="BED3+ reference features; the midpoint is used unless --point-col selects another coordinate.")
    parser.add_argument(
        "--blacklist-bed",
        type=Path,
        help=(
            "BED blacklist. Anchors overlapping it are excluded and overlapping "
            "bases elsewhere in each window are retained as missing values."
        ),
    )
    parser.add_argument("-n", "--nucleosome-bed", type=Path, help="Optional BED of nucleosome positions used to replace each feature midpoint with a strand-relative nucleosome centre.")
    parser.add_argument("--nucleosome-offset", type=int, default=1, help="Non-zero strand-relative nucleosome rank selected from --nucleosome-bed; positive is downstream and negative upstream (default: 1).")
    parser.add_argument("--state-bed", type=Path, help="Optional BED mask retaining reference features that overlap at least one state interval.")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory for outputs (default: current directory).")
    parser.add_argument("--output-prefix", help="Output filename prefix; default: derived from the region, signal, and optional selector inputs.")
    parser.add_argument("--heatmap-output", "--output", dest="heatmap_output", type=Path, help="Explicit heatmap image stem/path; final extension follows --plot-format (default stem: <output-dir>/<prefix>_heatmap).")
    parser.add_argument("--heatmap-matrix-output", type=Path, help="Explicit sorted/plotted matrix path; default: <output-dir>/<prefix>_heatmap_matrix.tsv.gz.")
    parser.add_argument("--aggregate-output", type=Path, help="Explicit complete mean-profile path; default: <output-dir>/<prefix>_aggregate_all.tsv.")
    parser.add_argument("--plotted-mean-output", type=Path, help="Explicit mean-of-plotted-rows path; default: <output-dir>/<prefix>_heatmap_mean.tsv.")
    parser.add_argument("--mean-plot-output", type=Path, help="Explicit plotted-mean image stem/path; final extension follows --plot-format (default stem: <output-dir>/<prefix>_heatmap_mean).")
    parser.add_argument("--summary-output", type=Path, help="Explicit processing-summary path; default: <output-dir>/<prefix>_summary.tsv.")
    parser.add_argument("--window-half", type=int, default=2500, help="Bases extracted on each side of the aggregation centre.")
    parser.add_argument("--chrom-col", type=int, default=1, help="One-based chromosome column in --region-bed (default: 1).")
    parser.add_argument("--start-col", type=int, default=2, help="One-based start-coordinate column in --region-bed (default: 2).")
    parser.add_argument("--end-col", type=int, default=3, help="One-based end-coordinate column in --region-bed (default: 3).")
    parser.add_argument("--strand-col", type=int, default=6, help="One-based strand column in --region-bed (default: 6).")
    parser.add_argument("--point-col", type=int, default=0, help="One-based absolute genomic centre column; 0 uses the interval midpoint (default: 0).")
    parser.add_argument("--skip-header", action="store_true", help="Skip the first physical line of --region-bed.")
    parser.add_argument(
        "--missing-strand", choices=["forward", "random", "error"], default="forward",
        help="How to handle missing or invalid strand values.",
    )
    parser.add_argument("--zero-thresh", "--zero_thresh", dest="zero_thresh", type=int, default=5, help="Reject vectors containing this many consecutive zeros; 0 disables the filter (default: 5).")
    parser.add_argument("--max-score", "--max_score", dest="max_score", type=float, default=300.0, help="Reject vectors containing a value greater than this threshold (default: 300).")
    parser.add_argument(
        "--nan-to-zero", "--nan_to_zero",
        dest="nan_to_zero", action=argparse.BooleanOptionalAction, default=True,
        help=(
            "Convert missing BigWig values (NaN) to zero before filtering and aggregation. "
            "Use --no-nan-to-zero to reject windows containing missing values."
        ),
    )
    parser.add_argument("--max-heatmap-rows", "--max-lines", dest="max_heatmap_rows", type=int, help="Maximum accepted rows included in the plotted/exported matrix; the complete aggregate still uses all accepted rows.")
    parser.add_argument("--subsample-mode", choices=["first", "random"], default="first", help="Choose the first rows or a reservoir sample when --max-heatmap-rows is set (default: first; use --seed for reproducible random sampling).")
    parser.add_argument("--stop-after-valid", type=int, help="Stop reading after this many accepted features; omit to scan the complete region BED.")
    parser.add_argument("--seed", type=int, help="Seed for random missing-strand orientation and random heatmap-row sampling.")
    parser.add_argument("--breadth", type=float, default=1.0, help="Central fraction of matrix columns retained for heatmap plotting/export, in (0,1] (default: 1).")
    parser.add_argument("--vmin", type=float, help="Explicit lower heatmap colour limit; default: symmetric automatic range.")
    parser.add_argument("--vmax", type=float, help="Explicit upper heatmap colour limit; default: symmetric automatic range.")
    parser.add_argument(
        "--sort-mode", choices=["center", "rise_after_min", "mean_absolute", "absmean", "max", "unsorted"],
        default="mean_absolute",
        help="Heatmap-row ordering: centre value, rise after upstream minimum, mean absolute signal, maximum, or input order (default: mean_absolute).",
    )
    parser.add_argument("--axis-label", "--label-replace", dest="axis_label", default="+1 NPS peak", help="Feature label used on relative-position axes (default: +1 NPS peak).")
    parser.add_argument("--mean-ylim", "--mean_ylim", dest="mean_ylim", type=float, help="Optional symmetric positive y-axis limit for the plotted mean.")
    parser.add_argument("--colorbar-label", default="Nucleosome Protection Score (NPS)", help="Heatmap colour-bar label.")
    parser.add_argument("--mean-ylabel", default="Mean NPS", help="Y-axis label for the mean-profile plot (default: Mean NPS).")
    parser.add_argument(
        "--nrl",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Call peaks across the complete aggregate alignment and fit separate "
            "positive- and negative-direction repeat lengths (default: enabled)."
        ),
    )
    parser.add_argument(
        "--nrl-peak-resolution",
        type=float,
        default=160.0,
        help=(
            "Minimum separation for unified aggregate peak calling. The default "
            "160 bp resolution gives 51 bp detection smoothing and 21 bp summit "
            "refinement. The called peak nearest 0 within half this resolution is "
            "shared as peak number 0 by both directional regressions."
        ),
    )
    parser.add_argument(
        "--nrl-regression-min",
        type=float,
        default=0.0,
        help=(
            "Minimum absolute distance from position 0 included in each directional "
            "regression; this does not restrict smoothing or peak calling (default: 0)."
        ),
    )
    parser.add_argument(
        "--nrl-regression-max",
        type=float,
        default=None,
        help=(
            "Maximum absolute distance from position 0 included in each directional "
            "regression; this does not restrict smoothing or peak calling "
            "(default: --window-half)."
        ),
    )
    parser.add_argument(
        "--nrl-regression-exclusion-start",
        "--nrl-exclusion-start",
        type=float,
        default=None,
        help=(
            "Inclusive signed start position of an optional interval excluded from "
            "both directional regressions. Supply together with "
            "--nrl-regression-exclusion-end. Peak calling and the unified profile "
            "are unchanged."
        ),
    )
    parser.add_argument(
        "--nrl-regression-exclusion-end",
        "--nrl-exclusion-end",
        type=float,
        default=None,
        help=(
            "Inclusive signed end position of an optional interval excluded from "
            "both directional regressions. Supply together with "
            "--nrl-regression-exclusion-start."
        ),
    )
    parser.add_argument("--dpi", type=int, default=300, help="Heatmap raster resolution in dots per inch (default: 300).")
    parser.add_argument(
        "-c", "--contigs", nargs="+", default=["all"],
        help="Contigs to process. With --cores > 1, each contig is written separately before combination.",
    )
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser)
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def _run_serial(args: argparse.Namespace) -> int:
    config = AlignmentConfig(**{
        key: value for key, value in vars(args).items()
        if key in AlignmentConfig.__dataclass_fields__
    })
    run_alignment(config, progress=ProgressReporter("aggregate"))
    return 0


def run_from_args(args: argparse.Namespace) -> int:
    from nucleosuite.aggregate_parallel import run_aggregate_per_contig
    return run_aggregate_per_contig(args, _run_serial)


def main(argv=None) -> int:
    """Standalone dispatcher for ``nucleosuite aggregate``."""
    root = argparse.ArgumentParser(prog="nucleosuite")
    subparsers = root.add_subparsers(dest="command", required=True)
    add_aggregate_parser(subparsers)
    args = root.parse_args(["aggregate", *(argv or [])])
    return run_from_args(args)
