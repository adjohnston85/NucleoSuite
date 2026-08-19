#!/usr/bin/env python3
"""Compare one main nucleosome callset with one or more comparison callsets."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy import stats

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.progress import ProgressReporter
from nucleosuite.compare_positions_legacy import (
    PositionRecord,
    MatchedPair,
    MatchResult,
    _CandidateCursor,
    _records_by_chrom,
    match_many_to_one,
    match_positions,
    match_unique,
    read_positions,
)


@dataclass(frozen=True)
class ComparisonSpec:
    label: str
    path: Path


@dataclass
class ComparisonArrays:
    label: str
    path: Path
    main_count: int
    compare_count: int
    query_source: str
    query_count: int
    target_count: int
    unmatched_no_target_chrom: int
    unmatched_distance: int
    unmatched_unique: int
    main_line: np.ndarray
    main_scores: np.ndarray
    compare_scores: np.ndarray
    signed_distance: np.ndarray
    absolute_distance: np.ndarray
    percentile: np.ndarray
    group_labels: np.ndarray

    @property
    def pair_count(self) -> int:
        return int(self.main_scores.size)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite compare-positions",
        description=(
            "Compare one main nucleosome BED with one or more comparison BEDs. "
            "Each comparison is matched once using the smaller callset as the query, "
            "with one-to-one unique matching. Matched pairs are ranked by the main "
            "BED score and divided into percentile groups (quartiles by default)."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument(
        "--main-bed",
        help="Main nucleosome BED, BED.gz, or bigBed file.",
    )
    parser.add_argument(
        "--compare-bed",
        dest="compare_beds",
        action="append",
        default=[],
        metavar="[LABEL=]BED",
        help=(
            "Comparison BED. Repeat for multiple callsets. Prefix a path with LABEL= "
            "to set its plot/legend label; otherwise the basename is used."
        ),
    )
    # Backward-compatible two-file aliases. They are intentionally omitted from
    # the documentation; the redesigned interface is --main-bed/--compare-bed.
    parser.add_argument("--bed-a", dest="legacy_bed_a", help=argparse.SUPPRESS)
    parser.add_argument("--bed-b", dest="legacy_bed_b", help=argparse.SUPPRESS)
    parser.add_argument("--main-label", default=None, help="Display label for the main BED.")
    parser.add_argument(
        "--main-summit-column", type=int, default=None,
        help="One-based absolute summit column for the main BED; default uses the BED midpoint.",
    )
    parser.add_argument(
        "--compare-summit-column", type=int, default=None,
        help="One-based absolute summit column used for all comparison BEDs; default uses BED midpoints.",
    )
    parser.add_argument(
        "--main-score-column", type=int, default=5,
        help="One-based numeric score column for the main BED (default: 5).",
    )
    parser.add_argument(
        "--compare-score-column", type=int, default=5,
        help="One-based numeric score column used for comparison BEDs (default: 5).",
    )
    parser.add_argument(
        "--blacklist-bed",
        help="BED blacklist; complete overlapping records are excluded from all inputs.",
    )
    parser.add_argument(
        "--max-distance", type=float, default=None,
        help="Maximum allowed absolute summit distance in bp; omit for no limit.",
    )
    parser.add_argument(
        "--percentile-interval", type=int, default=25,
        help=(
            "Width of main-score percentile groups among matched pairs. The default "
            "25 gives quartiles: 0-25, 25-50, 50-75, and 75-100."
        ),
    )
    parser.add_argument(
        "--distance-bins", default="5,10,20,50,100",
        help="Comma-separated upper bounds for score-correlation-by-distance summaries.",
    )
    parser.add_argument(
        "--histogram-bin-width", type=float, default=1.0,
        help="Distance-distribution bin width in bp (default: 1).",
    )
    parser.add_argument(
        "--histogram-x-min", type=float, default=0.0,
        help="Displayed lower x-axis limit for distance distributions (default: 0).",
    )
    parser.add_argument(
        "--histogram-x-max", type=float, default=300.0,
        help="Displayed upper x-axis limit for distance distributions (default: 300).",
    )
    parser.add_argument(
        "--distance-x-major-tick", type=float, default=None,
        help="Major x-axis tick interval in bp for numeric distance plots; default is automatic.",
    )
    parser.add_argument(
        "--distance-x-minor-tick", type=float, default=None,
        help="Minor x-axis tick interval in bp for numeric distance plots; default is automatic.",
    )
    parser.add_argument(
        "--score-normalization", choices=("raw", "zscore", "percentile"), default="zscore",
        help="Score representation for main-vs-comparison score-agreement plots (default: zscore).",
    )
    parser.add_argument(
        "--score-correlation", choices=("spearman", "pearson", "both"), default="spearman",
        help="Correlation displayed for main-vs-comparison score agreement (default: spearman).",
    )
    parser.add_argument(
        "--score-z-limit", type=float, default=10.0,
        help="Symmetric score z-axis limit for score-agreement plots; 0 disables (default: 10).",
    )
    parser.add_argument(
        "--score-distance-type", choices=("absolute", "signed"), default="absolute",
        help="Distance used for main-score-versus-distance plots (default: absolute).",
    )
    parser.add_argument(
        "--score-distance-correlation", choices=("spearman", "pearson", "both"), default="spearman",
        help="Correlation displayed for main peak score versus matched distance (default: spearman).",
    )
    parser.add_argument(
        "--score-distance-plot", choices=("hexbin", "scatter"), default="hexbin",
        help="Rendering for main-score-versus-distance plots (default: hexbin).",
    )
    parser.add_argument(
        "--plot-max-points", type=int, default=200000,
        help="Maximum pairs drawn in scatter/hexbin plots; 0 draws all pairs (default: 200000).",
    )
    parser.add_argument("--plot-seed", type=int, default=1, help="Plot subsampling seed (default: 1).")
    parser.add_argument(
        "--percentile-boxplot-y-max", type=float, default=500.0,
        help="Displayed upper y-axis limit for percentile distance boxplots; 0 disables (default: 500).",
    )
    parser.add_argument(
        "--stats", action="store_true",
        help="Perform pairwise comparison tests within each percentile group and annotate the boxplot.",
    )
    parser.add_argument(
        "--stats-test", choices=("nonparametric", "parametric"), default="nonparametric",
        help=(
            "Statistical family for within-percentile pairwise tests. Nonparametric "
            "uses paired Wilcoxon where shared main calls permit pairing, otherwise "
            "Mann-Whitney U; parametric uses paired t-tests or Welch t-tests (default: nonparametric)."
        ),
    )
    parser.add_argument(
        "--p-adjust", choices=("holm", "none"), default="holm",
        help="Multiple-testing correction applied independently within each percentile group (default: holm).",
    )
    parser.add_argument(
        "--p-display", choices=("value", "stars"), default="value",
        help="Display adjusted/raw p-values or significance stars above boxplots (default: value).",
    )
    parser.add_argument(
        "--skip-pairs-tsv", action="store_true",
        help="Do not write the detailed matched-pair TSV for each comparison.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Figure DPI; --plot-dpi takes precedence.")
    parser.add_argument(
        "-o", "--output-prefix", default=None,
        help="Output prefix. Default uses the main BED basename followed by _compare_positions.",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages.")
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive analysis. This "
            "budget is independent of suite --cores (default: 1)."
        ),
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def _default_label(path: str | Path) -> str:
    name = Path(path).name
    for suffix in (".bed.gz", ".bed.bgz", ".bigBed", ".bigbed", ".bed", ".bb", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("_.-")
    return token or "comparison"


def _parse_compare_spec(value: str) -> ComparisonSpec:
    if "=" in value:
        label, path_text = value.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError("--compare-bed LABEL=BED requires a non-empty label.")
    else:
        path_text = value
        label = _default_label(path_text)
    path = Path(path_text)
    return ComparisonSpec(label=label, path=path)


def _resolve_inputs(args: argparse.Namespace) -> tuple[Path, list[ComparisonSpec]]:
    main_value = getattr(args, "main_bed", None)
    compare_values = list(getattr(args, "compare_beds", None) or [])
    legacy_a = getattr(args, "legacy_bed_a", None)
    legacy_b = getattr(args, "legacy_bed_b", None)
    if main_value is None and legacy_a:
        main_value = legacy_a
    if not compare_values and legacy_b:
        compare_values = [legacy_b]
    if not main_value:
        raise ValueError("--main-bed is required.")
    if not compare_values:
        raise ValueError("At least one --compare-bed is required.")
    specs = [_parse_compare_spec(str(value)) for value in compare_values]
    labels = [spec.label for spec in specs]
    if len(labels) != len(set(labels)):
        duplicates = sorted({label for label in labels if labels.count(label) > 1})
        raise ValueError("Comparison labels must be unique: " + ", ".join(duplicates))
    return Path(main_value), specs


def _validate_args(args: argparse.Namespace) -> None:
    if not 1 <= int(args.percentile_interval) <= 100:
        raise ValueError("--percentile-interval must be between 1 and 100.")
    if args.max_distance is not None and (not math.isfinite(float(args.max_distance)) or float(args.max_distance) < 0):
        raise ValueError("--max-distance must be finite and zero or greater.")
    if not math.isfinite(float(args.histogram_bin_width)) or float(args.histogram_bin_width) <= 0:
        raise ValueError("--histogram-bin-width must be greater than zero.")
    if float(args.histogram_x_min) < 0 or float(args.histogram_x_max) <= float(args.histogram_x_min):
        raise ValueError("Histogram limits require 0 <= x-min < x-max.")
    if int(args.plot_max_points) < 0:
        raise ValueError("--plot-max-points must be zero or greater.")
    if int(args.dpi) < 1:
        raise ValueError("--dpi must be at least 1.")
    if not math.isfinite(float(args.score_z_limit)) or float(args.score_z_limit) < 0:
        raise ValueError("--score-z-limit must be finite and zero or greater.")
    if not math.isfinite(float(args.percentile_boxplot_y_max)) or float(args.percentile_boxplot_y_max) < 0:
        raise ValueError("--percentile-boxplot-y-max must be finite and zero or greater.")
    for value, option in (
        (args.main_score_column, "--main-score-column"),
        (args.compare_score_column, "--compare-score-column"),
        (args.main_summit_column, "--main-summit-column"),
        (args.compare_summit_column, "--compare-summit-column"),
    ):
        if value is not None and int(value) < 1:
            raise ValueError(f"{option} must be a one-based column number of at least 1.")


def _percentile_group_bounds(interval: int) -> list[tuple[int, int, str]]:
    groups: list[tuple[int, int, str]] = []
    lower = 0
    while lower < 100:
        upper = min(100, lower + int(interval))
        groups.append((lower, upper, f"{lower}-{upper}"))
        lower = upper
    return groups


def _assign_matched_percentiles(main_scores: np.ndarray, interval: int) -> tuple[np.ndarray, np.ndarray]:
    """Rank matched pairs by main score and return percentile + group label arrays."""
    n = int(main_scores.size)
    if n == 0:
        return np.asarray([], dtype=np.uint8), np.asarray([], dtype=object)
    # Stable sorting keeps ties deterministic without changing the user's explicit
    # requirement that matched pairs are sorted by main BED score.
    order = np.argsort(main_scores, kind="stable")
    percentile = np.empty(n, dtype=np.uint8)
    ranks = np.arange(1, n + 1, dtype=float)
    percentile[order] = np.ceil(ranks * 100.0 / n).astype(np.uint8)
    labels = np.empty(n, dtype=object)
    for lower, upper, label in _percentile_group_bounds(interval):
        mask = (percentile > lower) & (percentile <= upper)
        labels[mask] = label
    return percentile, labels


def _arrays_from_match(
    label: str,
    path: Path,
    result: MatchResult,
    main_count: int,
    compare_count: int,
    interval: int,
) -> ComparisonArrays:
    pairs = result.pairs
    n = len(pairs)
    main_line = np.fromiter((pair.a.line_number for pair in pairs), dtype=np.int64, count=n)
    main_scores = np.fromiter((pair.a.score for pair in pairs), dtype=np.float64, count=n)
    compare_scores = np.fromiter((pair.b.score for pair in pairs), dtype=np.float64, count=n)
    signed = np.fromiter((pair.signed_distance for pair in pairs), dtype=np.int64, count=n)
    absolute = np.abs(signed).astype(np.float64)
    percentile, group_labels = _assign_matched_percentiles(main_scores, interval)
    return ComparisonArrays(
        label=label,
        path=path,
        main_count=main_count,
        compare_count=compare_count,
        query_source="main" if result.query_source == "A" else "comparison",
        query_count=result.query_count,
        target_count=result.target_count,
        unmatched_no_target_chrom=result.unmatched_no_target_chrom,
        unmatched_distance=result.unmatched_distance,
        unmatched_unique=result.unmatched_unique,
        main_line=main_line,
        main_scores=main_scores,
        compare_scores=compare_scores,
        signed_distance=signed,
        absolute_distance=absolute,
        percentile=percentile,
        group_labels=group_labels,
    )


def _zscores(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        return values.astype(float)
    std = float(np.std(values))
    if not math.isfinite(std) or std == 0:
        return np.zeros(values.size, dtype=float)
    return (values - float(np.mean(values))) / std


def _percentile_ranks(values: np.ndarray) -> np.ndarray:
    if values.size == 0:
        return values.astype(float)
    ranks = stats.rankdata(values, method="average")
    if values.size == 1:
        return np.asarray([100.0])
    return 100.0 * (ranks - 1.0) / (values.size - 1.0)


def _selected_scores(result: ComparisonArrays, normalization: str) -> tuple[np.ndarray, np.ndarray, str]:
    if normalization == "raw":
        return result.main_scores, result.compare_scores, "raw score"
    if normalization == "percentile":
        return _percentile_ranks(result.main_scores), _percentile_ranks(result.compare_scores), "score percentile rank"
    return _zscores(result.main_scores), _zscores(result.compare_scores), "score z-score"


def _safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 2 or np.all(x == x[0]) or np.all(y == y[0]):
        return math.nan, math.nan
    try:
        if method == "spearman":
            output = stats.spearmanr(x, y)
        elif method == "pearson":
            output = stats.pearsonr(x, y)
        else:
            raise ValueError(method)
        return float(output.statistic), float(output.pvalue)
    except Exception:
        return math.nan, math.nan


def _linear_stats(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 2 or np.all(x == x[0]):
        return math.nan, math.nan, math.nan, math.nan
    try:
        regression = stats.linregress(x, y)
        return float(regression.slope), float(regression.intercept), float(regression.rvalue ** 2), float(regression.pvalue)
    except Exception:
        return math.nan, math.nan, math.nan, math.nan


def _plot_indices(count: int, maximum: int, seed: int) -> np.ndarray:
    if maximum <= 0 or count <= maximum:
        return np.arange(count, dtype=np.int64)
    rng = np.random.default_rng(int(seed))
    return np.sort(rng.choice(count, size=int(maximum), replace=False))


def _parse_distance_bins(text: str) -> list[float]:
    try:
        values = sorted({float(value.strip()) for value in str(text).split(",") if value.strip()})
    except ValueError as exc:
        raise ValueError("--distance-bins must be comma-separated numeric upper bounds.") from exc
    if not values or any(not math.isfinite(value) or value < 0 for value in values):
        raise ValueError("--distance-bins requires finite non-negative upper bounds.")
    return values


def _distance_label(lower: float | None, upper: float | None) -> str:
    def fmt(value: float) -> str:
        return str(int(value)) if float(value).is_integer() else f"{value:g}"
    if lower is None:
        return f"0-{fmt(upper or 0)}"
    if upper is None:
        return f">{fmt(lower)}"
    return f"{fmt(lower)}-{fmt(upper)}"


def _distance_bin_rows(result: ComparisonArrays, normalization: str, bounds: Sequence[float]) -> list[dict[str, object]]:
    main_scores, compare_scores, score_label = _selected_scores(result, normalization)
    rows: list[dict[str, object]] = []
    lower: float | None = None
    for upper in [*bounds, None]:
        if lower is None:
            mask = result.absolute_distance <= float(upper) if upper is not None else np.ones(result.pair_count, dtype=bool)
        elif upper is None:
            mask = result.absolute_distance > float(lower)
        else:
            mask = (result.absolute_distance > float(lower)) & (result.absolute_distance <= float(upper))
        x = main_scores[mask]
        y = compare_scores[mask]
        spearman, _ = _safe_corr(x, y, "spearman")
        pearson, _ = _safe_corr(x, y, "pearson")
        rows.append({
            "comparison": result.label,
            "distance_bin": _distance_label(lower, upper),
            "lower_exclusive": "" if lower is None else lower,
            "upper_inclusive": "" if upper is None else upper,
            "pair_count": int(mask.sum()),
            "spearman_score_correlation": spearman,
            "pearson_score_correlation": pearson,
            "plot_score_axis_label": score_label,
        })
        lower = upper
    return rows


def _histogram_rows(result: ComparisonArrays, x_min: float, x_max: float, width: float) -> list[dict[str, object]]:
    edges = np.arange(float(x_min), float(x_max) + float(width), float(width), dtype=float)
    if edges[-1] < float(x_max):
        edges = np.append(edges, float(x_max))
    counts, edges = np.histogram(result.absolute_distance, bins=edges)
    rows: list[dict[str, object]] = []
    total = max(1, int(result.pair_count))
    for index, count in enumerate(counts):
        rows.append({
            "comparison": result.label,
            "bin_start_inclusive": edges[index],
            "bin_end_exclusive": edges[index + 1],
            "pair_count": int(count),
            "density": float(count) / (total * float(edges[index + 1] - edges[index])),
        })
    return rows


def _write_tsv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str] | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def _write_pairs(path: Path, result: ComparisonArrays, pairs: Sequence[MatchedPair], main_label: str, args: argparse.Namespace) -> Path:
    selected_main, selected_compare, _ = _selected_scores(result, args.score_normalization)
    indices = set(_plot_indices(result.pair_count, args.plot_max_points, args.plot_seed).tolist())
    fields = [
        "comparison", "pair_id", "query_source", "chrom",
        "main_start", "main_end", "main_name", "main_summit", "main_score", "main_score_normalized", "main_line_number",
        "compare_start", "compare_end", "compare_name", "compare_summit", "compare_score", "compare_score_normalized", "compare_line_number",
        "signed_distance_compare_minus_main", "absolute_distance", "main_score_percentile", "percentile_group",
        "plot_selected", "plot_score_normalization", "plot_score_z_limit", "plot_correlation_method", "plot_label_a", "plot_label_b",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index, pair in enumerate(pairs):
            writer.writerow({
                "comparison": result.label,
                "pair_id": f"pair_{index + 1:09d}",
                "query_source": result.query_source,
                "chrom": pair.a.chrom,
                "main_start": pair.a.start,
                "main_end": pair.a.end,
                "main_name": pair.a.name,
                "main_summit": pair.a.summit,
                "main_score": pair.a.score,
                "main_score_normalized": selected_main[index],
                "main_line_number": pair.a.line_number,
                "compare_start": pair.b.start,
                "compare_end": pair.b.end,
                "compare_name": pair.b.name,
                "compare_summit": pair.b.summit,
                "compare_score": pair.b.score,
                "compare_score_normalized": selected_compare[index],
                "compare_line_number": pair.b.line_number,
                "signed_distance_compare_minus_main": pair.signed_distance,
                "absolute_distance": pair.absolute_distance,
                "main_score_percentile": int(result.percentile[index]),
                "percentile_group": str(result.group_labels[index]),
                "plot_selected": 1 if index in indices else 0,
                "plot_score_normalization": args.score_normalization,
                "plot_score_z_limit": args.score_z_limit,
                "plot_correlation_method": args.score_correlation,
                "plot_label_a": main_label,
                "plot_label_b": result.label,
            })
    return path


def _percentile_rows(result: ComparisonArrays) -> list[dict[str, object]]:
    return [
        {
            "comparison": result.label,
            "main_line_number": int(result.main_line[index]),
            "main_score": float(result.main_scores[index]),
            "main_score_percentile": int(result.percentile[index]),
            "percentile_group": str(result.group_labels[index]),
            "signed_distance_compare_minus_main": int(result.signed_distance[index]),
            "absolute_distance": float(result.absolute_distance[index]),
        }
        for index in range(result.pair_count)
    ]


def _percentile_summary_rows(results: Sequence[ComparisonArrays], interval: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for lower, upper, label in _percentile_group_bounds(interval):
            mask = result.group_labels == label
            values = result.absolute_distance[mask]
            if values.size:
                q1, median, q3 = np.percentile(values, [25, 50, 75])
                mean = float(np.mean(values))
                minimum = float(np.min(values))
                maximum = float(np.max(values))
            else:
                minimum = q1 = median = mean = q3 = maximum = math.nan
            rows.append({
                "comparison": result.label,
                "percentile_group": label,
                "group_lower_percentile": lower,
                "group_upper_percentile": upper,
                "matched_pair_count": int(values.size),
                "minimum_absolute_distance": minimum,
                "q1_absolute_distance": q1,
                "median_absolute_distance": median,
                "mean_absolute_distance": mean,
                "q3_absolute_distance": q3,
                "maximum_absolute_distance": maximum,
            })
    return rows


def _holm_adjust(p_values: Sequence[float]) -> list[float]:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.size, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if finite_indices.size == 0:
        return adjusted.tolist()
    order = finite_indices[np.argsort(values[finite_indices])]
    m = int(order.size)
    running = 0.0
    for rank, index in enumerate(order, start=1):
        candidate = min(1.0, float(values[index]) * (m - rank + 1))
        running = max(running, candidate)
        adjusted[index] = running
    return adjusted.tolist()


def _significance(p: float) -> str:
    if not math.isfinite(p):
        return "NA"
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _paired_group_values(a: ComparisonArrays, b: ComparisonArrays, group: str) -> tuple[np.ndarray, np.ndarray]:
    mask_a = a.group_labels == group
    mask_b = b.group_labels == group
    map_a = {int(line): float(value) for line, value in zip(a.main_line[mask_a], a.absolute_distance[mask_a])}
    map_b = {int(line): float(value) for line, value in zip(b.main_line[mask_b], b.absolute_distance[mask_b])}
    shared = sorted(set(map_a).intersection(map_b))
    return (
        np.asarray([map_a[line] for line in shared], dtype=float),
        np.asarray([map_b[line] for line in shared], dtype=float),
    )


def _pairwise_statistics(results: Sequence[ComparisonArrays], interval: int, family: str, adjustment: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for _lower, _upper, group in _percentile_group_bounds(interval):
        group_rows: list[dict[str, object]] = []
        for first, second in itertools.combinations(results, 2):
            values_first = first.absolute_distance[first.group_labels == group]
            values_second = second.absolute_distance[second.group_labels == group]
            paired_first, paired_second = _paired_group_values(first, second, group)
            paired = (
                paired_first.size >= 2
                and paired_first.size == values_first.size
                and paired_second.size == values_second.size
            )
            statistic = p_value = math.nan
            test_name = ""
            try:
                if family == "nonparametric":
                    if paired:
                        test_name = "Wilcoxon signed-rank"
                        if np.allclose(paired_first, paired_second):
                            statistic, p_value = 0.0, 1.0
                        else:
                            output = stats.wilcoxon(paired_first, paired_second, alternative="two-sided")
                            statistic, p_value = float(output.statistic), float(output.pvalue)
                    elif values_first.size and values_second.size:
                        test_name = "Mann-Whitney U"
                        output = stats.mannwhitneyu(values_first, values_second, alternative="two-sided")
                        statistic, p_value = float(output.statistic), float(output.pvalue)
                else:
                    if paired:
                        test_name = "paired t-test"
                        output = stats.ttest_rel(paired_first, paired_second, nan_policy="omit")
                        statistic, p_value = float(output.statistic), float(output.pvalue)
                    elif values_first.size >= 2 and values_second.size >= 2:
                        test_name = "Welch t-test"
                        output = stats.ttest_ind(values_first, values_second, equal_var=False, nan_policy="omit")
                        statistic, p_value = float(output.statistic), float(output.pvalue)
            except Exception:
                statistic = p_value = math.nan
            group_rows.append({
                "percentile_group": group,
                "comparison_1": first.label,
                "comparison_2": second.label,
                "test": test_name or ("not_tested" if not paired else "failed"),
                "paired": bool(paired),
                "n_1": int(values_first.size),
                "n_2": int(values_second.size),
                "n_paired": int(paired_first.size),
                "statistic": statistic,
                "p_value": p_value,
            })
        adjusted = _holm_adjust([float(row["p_value"]) for row in group_rows]) if adjustment == "holm" else [float(row["p_value"]) for row in group_rows]
        for row, p_adj in zip(group_rows, adjusted):
            row["p_adjusted"] = p_adj
            row["p_adjustment"] = adjustment
            row["significance"] = _significance(float(p_adj))
        rows.extend(group_rows)
    return rows


def _score_distance_stats(result: ComparisonArrays, distance_type: str) -> dict[str, object]:
    y = result.absolute_distance if distance_type == "absolute" else result.signed_distance.astype(float)
    spearman, spearman_p = _safe_corr(result.main_scores, y, "spearman")
    pearson, pearson_p = _safe_corr(result.main_scores, y, "pearson")
    slope, intercept, r2, regression_p = _linear_stats(result.main_scores, y)
    return {
        "comparison": result.label,
        "matched_pair_count": result.pair_count,
        "distance_type": distance_type,
        "spearman_rho": spearman,
        "spearman_p_value": spearman_p,
        "pearson_r": pearson,
        "pearson_p_value": pearson_p,
        "linear_slope_bp_per_score_unit": slope,
        "linear_intercept": intercept,
        "linear_r_squared": r2,
        "linear_slope_p_value": regression_p,
        "median_absolute_distance": float(np.median(result.absolute_distance)) if result.pair_count else math.nan,
        "mean_absolute_distance": float(np.mean(result.absolute_distance)) if result.pair_count else math.nan,
    }


def _summary_row(result: ComparisonArrays, main_path: Path, main_label: str) -> dict[str, object]:
    score_spearman, score_spearman_p = _safe_corr(result.main_scores, result.compare_scores, "spearman")
    score_pearson, score_pearson_p = _safe_corr(result.main_scores, result.compare_scores, "pearson")
    return {
        "main_label": main_label,
        "main_bed": str(main_path),
        "comparison": result.label,
        "compare_bed": str(result.path),
        "main_position_count": result.main_count,
        "compare_position_count": result.compare_count,
        "query_source": result.query_source,
        "query_position_count": result.query_count,
        "target_position_count": result.target_count,
        "matched_pair_count": result.pair_count,
        "unmatched_query_no_target_chromosome": result.unmatched_no_target_chrom,
        "unmatched_query_beyond_maximum_distance": result.unmatched_distance,
        "unmatched_query_unique_assignment": result.unmatched_unique,
        "median_absolute_distance": float(np.median(result.absolute_distance)) if result.pair_count else math.nan,
        "mean_absolute_distance": float(np.mean(result.absolute_distance)) if result.pair_count else math.nan,
        "spearman_main_vs_compare_score": score_spearman,
        "spearman_main_vs_compare_score_p_value": score_spearman_p,
        "pearson_main_vs_compare_score": score_pearson,
        "pearson_main_vs_compare_score_p_value": score_pearson_p,
    }


def _comparison_colors(count: int):
    from nucleosuite.plotting import category_colors
    return category_colors(count)


def _plot_score_agreement(prefix: Path, result: ComparisonArrays, main_label: str, args: argparse.Namespace) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import plot_path, save_figure

    x, y, axis_label = _selected_scores(result, args.score_normalization)
    indices = _plot_indices(result.pair_count, args.plot_max_points, args.plot_seed)
    figure, axis = plt.subplots(figsize=(7.5, 6.5))
    scatter = axis.scatter(
        x[indices], y[indices], c=result.absolute_distance[indices], s=9,
        alpha=0.55, linewidths=0, cmap="viridis", rasterized=True,
    )
    colorbar = figure.colorbar(scatter, ax=axis)
    colorbar.set_label("Absolute summit distance (bp)")
    if args.score_normalization == "zscore" and float(args.score_z_limit) > 0:
        axis.set_xlim(-float(args.score_z_limit), float(args.score_z_limit))
        axis.set_ylim(-float(args.score_z_limit), float(args.score_z_limit))
    axis.set_xlabel(f"{main_label} {axis_label}")
    axis.set_ylabel(f"{result.label} {axis_label}")
    axis.set_title("Score agreement coloured by summit distance")
    annotation: list[str] = []
    if args.score_correlation in {"spearman", "both"}:
        coefficient, p_value = _safe_corr(x, y, "spearman")
        annotation.append(f"Spearman ρ = {coefficient:.3f} (p={p_value:.3g})")
    if args.score_correlation in {"pearson", "both"}:
        coefficient, p_value = _safe_corr(x, y, "pearson")
        annotation.append(f"Pearson r = {coefficient:.3f} (p={p_value:.3g})")
    _slope, _intercept, r2, _p = _linear_stats(x, y)
    annotation.append(f"R² = {r2:.3f}")
    annotation.append(f"n = {result.pair_count:,}")
    axis.text(0.02, 0.98, "\n".join(annotation), transform=axis.transAxes, va="top", ha="left")
    figure.tight_layout()
    output = plot_path(Path(f"{prefix}_{_safe_token(result.label)}_score_agreement.png"))
    output = save_figure(figure, output, default_dpi=args.dpi)
    plt.close(figure)
    return output


def _plot_score_distance(prefix: Path, result: ComparisonArrays, args: argparse.Namespace, stats_row: dict[str, object]) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import plot_path, save_figure

    y_all = result.absolute_distance if args.score_distance_type == "absolute" else result.signed_distance.astype(float)
    indices = _plot_indices(result.pair_count, args.plot_max_points, args.plot_seed)
    x = result.main_scores[indices]
    y = y_all[indices]
    figure, axis = plt.subplots(figsize=(8.0, 6.0))
    if args.score_distance_plot == "hexbin":
        artist = axis.hexbin(x, y, gridsize=60, mincnt=1, bins="log", rasterized=True)
        colorbar = figure.colorbar(artist, ax=axis)
        colorbar.set_label("log10 plotted pair count")
    else:
        axis.scatter(x, y, s=7, alpha=0.35, linewidths=0, rasterized=True)
    slope = float(stats_row["linear_slope_bp_per_score_unit"])
    intercept = float(stats_row["linear_intercept"])
    finite_x = result.main_scores[np.isfinite(result.main_scores)]
    if finite_x.size and math.isfinite(slope) and math.isfinite(intercept):
        endpoints = np.asarray([float(np.min(finite_x)), float(np.max(finite_x))])
        axis.plot(endpoints, intercept + slope * endpoints, linestyle=":", linewidth=1.2, label="Linear fit")
    axis.set_xlabel("Main peak score")
    axis.set_ylabel("Absolute matched distance (bp)" if args.score_distance_type == "absolute" else "Signed compare − main distance (bp)")
    axis.set_title(f"Main peak score versus distance: {result.label}")
    annotation: list[str] = []
    if args.score_distance_correlation in {"spearman", "both"}:
        annotation.append(f"Spearman ρ = {float(stats_row['spearman_rho']):.3f} (p={float(stats_row['spearman_p_value']):.3g})")
    if args.score_distance_correlation in {"pearson", "both"}:
        annotation.append(f"Pearson r = {float(stats_row['pearson_r']):.3f} (p={float(stats_row['pearson_p_value']):.3g})")
    annotation.append(f"Linear R² = {float(stats_row['linear_r_squared']):.3f}")
    annotation.append(f"n = {result.pair_count:,}")
    axis.text(0.02, 0.98, "\n".join(annotation), transform=axis.transAxes, va="top", ha="left")
    figure.tight_layout()
    output = plot_path(Path(f"{prefix}_{_safe_token(result.label)}_main_score_vs_distance.png"))
    output = save_figure(figure, output, default_dpi=args.dpi)
    plt.close(figure)
    return output


def _plot_distance_histograms(prefix: Path, results: Sequence[ComparisonArrays], args: argparse.Namespace) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import apply_distance_x_axis, apply_integer_y_axis, plot_path, save_figure

    edges = np.arange(float(args.histogram_x_min), float(args.histogram_x_max) + float(args.histogram_bin_width), float(args.histogram_bin_width))
    colors = _comparison_colors(len(results))
    figure, axis = plt.subplots(figsize=(9.0, 5.8))
    for color, result in zip(colors, results):
        counts, resolved_edges = np.histogram(result.absolute_distance, bins=edges)
        centres = (resolved_edges[:-1] + resolved_edges[1:]) / 2.0
        axis.plot(centres, counts, label=result.label, color=color, linewidth=1.5)
    axis.set_xlim(float(args.histogram_x_min), float(args.histogram_x_max))
    axis.set_xlabel("Absolute summit distance (bp)")
    axis.set_ylabel("Matched pairs")
    axis.set_title("Matched-position distance distributions")
    apply_distance_x_axis(axis, major_interval=args.distance_x_major_tick, minor_interval=args.distance_x_minor_tick)
    apply_integer_y_axis(axis)
    axis.legend(frameon=False)
    figure.tight_layout()
    output = plot_path(Path(f"{prefix}_distance_histogram.png"))
    output = save_figure(figure, output, default_dpi=args.dpi)
    plt.close(figure)
    return output


def _plot_correlation_by_distance(prefix: Path, rows: Sequence[dict[str, object]], results: Sequence[ComparisonArrays], args: argparse.Namespace) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import plot_path, save_figure

    figure, axis = plt.subplots(figsize=(9.0, 5.8))
    colors = _comparison_colors(len(results))
    for color, result in zip(colors, results):
        subset = [row for row in rows if row["comparison"] == result.label]
        labels = [str(row["distance_bin"]) for row in subset]
        x = np.arange(len(labels), dtype=float)
        if args.score_correlation in {"spearman", "both"}:
            axis.plot(x, [float(row["spearman_score_correlation"]) for row in subset], marker="o", color=color, label=(result.label if args.score_correlation == "spearman" else f"{result.label} Spearman"))
        if args.score_correlation in {"pearson", "both"}:
            axis.plot(x, [float(row["pearson_score_correlation"]) for row in subset], marker="s", linestyle="--", color=color, label=(result.label if args.score_correlation == "pearson" else f"{result.label} Pearson"))
    axis.axhline(0, linewidth=0.8, color="black")
    if rows:
        labels = [str(row["distance_bin"]) for row in rows if row["comparison"] == results[0].label]
        axis.set_xticks(np.arange(len(labels)), labels, rotation=35, ha="right")
    axis.set_ylim(-1.05, 1.05)
    axis.set_xlabel("Absolute summit-distance bin (bp)")
    axis.set_ylabel("Main/comparison score correlation")
    axis.set_title("Score correlation by summit-distance bin")
    axis.legend(frameon=False)
    figure.tight_layout()
    output = plot_path(Path(f"{prefix}_correlation_by_distance.png"))
    output = save_figure(figure, output, default_dpi=args.dpi)
    plt.close(figure)
    return output


def _format_p(p: float, adjusted: bool) -> str:
    if not math.isfinite(p):
        return "NA"
    label = "p_adj" if adjusted else "p"
    if p < 0.0001:
        return f"{label}<1e-4"
    return f"{label}={p:.3g}"


def _plot_percentile_boxplot(
    prefix: Path,
    results: Sequence[ComparisonArrays],
    interval: int,
    stats_rows: Sequence[dict[str, object]],
    args: argparse.Namespace,
) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from nucleosuite.plotting import plot_path, save_figure

    groups = [label for _lower, _upper, label in _percentile_group_bounds(interval)]
    comparisons = [result.label for result in results]
    colors = _comparison_colors(len(results))
    figure, axis = plt.subplots(figsize=(max(9.0, 1.8 * len(groups) + 3.0), 6.5))
    group_centres = np.arange(1, len(groups) + 1, dtype=float)
    if len(results) == 1:
        offsets = np.asarray([0.0])
        width = 0.48
    else:
        spread = min(0.70, 0.14 * len(results))
        offsets = np.linspace(-spread / 2.0, spread / 2.0, len(results))
        width = min(0.16, 0.70 / len(results))
    positions: dict[tuple[str, str], float] = {}
    for comp_index, (result, color) in enumerate(zip(results, colors)):
        for group_index, group in enumerate(groups):
            values = result.absolute_distance[result.group_labels == group]
            pos = float(group_centres[group_index] + offsets[comp_index])
            positions[(group, result.label)] = pos
            if not values.size:
                continue
            box = axis.boxplot(
                [values], positions=[pos], widths=width, patch_artist=True,
                showfliers=False, whis=1.5,
            )
            box["boxes"][0].set_facecolor(color)
            box["boxes"][0].set_edgecolor(color)
            for key in ("whiskers", "caps", "medians"):
                for artist in box[key]:
                    artist.set_color(color if key != "medians" else "black")
    axis.set_xticks(group_centres, groups)
    axis.set_xlim(0.45, len(groups) + 0.55)
    if float(args.percentile_boxplot_y_max) > 0:
        axis.set_ylim(0.0, float(args.percentile_boxplot_y_max))
    axis.set_xlabel("Main peak score percentile group")
    axis.set_ylabel("Absolute matched distance (bp)")
    axis.legend([Patch(facecolor=color, edgecolor=color) for color in colors], comparisons, frameon=False)

    max_levels = 0
    if args.stats and stats_rows:
        by_group: dict[str, list[dict[str, object]]] = {group: [] for group in groups}
        for row in stats_rows:
            by_group.setdefault(str(row["percentile_group"]), []).append(row)
        max_levels = 0
        for group in groups:
            group_stats = by_group.get(group, [])
            max_levels = max(max_levels, len(group_stats))
            for level, row in enumerate(group_stats):
                x1 = positions.get((group, str(row["comparison_1"])))
                x2 = positions.get((group, str(row["comparison_2"])))
                if x1 is None or x2 is None:
                    continue
                y = 1.03 + level * 0.065
                transform = axis.get_xaxis_transform()
                axis.plot([x1, x1, x2, x2], [y - 0.015, y, y, y - 0.015], transform=transform, clip_on=False, color="black", linewidth=0.8)
                p_value = float(row["p_adjusted"] if args.p_adjust == "holm" else row["p_value"])
                text = str(row["significance"]) if args.p_display == "stars" else _format_p(p_value, args.p_adjust == "holm")
                axis.text((x1 + x2) / 2.0, y + 0.006, text, transform=transform, ha="center", va="bottom", fontsize=7, clip_on=False)
    title_y = 1.02 if max_levels == 0 else 1.09 + max_levels * 0.065
    axis.set_title("Matched distance by main-score percentile", y=title_y)
    figure.tight_layout()
    output = plot_path(Path(f"{prefix}_percentile_distance_boxplot.png"))
    output = save_figure(figure, output, default_dpi=args.dpi)
    plt.close(figure)
    return output


def run_comparison(args: argparse.Namespace) -> dict[str, Path]:
    _validate_args(args)
    main_path, specs = _resolve_inputs(args)
    reporter = ProgressReporter("compare-positions", quiet=bool(getattr(args, "quiet", False)))
    reporter.stage("Reading main nucleosome positions")

    blacklist = None
    if getattr(args, "blacklist_bed", None):
        from nucleosuite.core.blacklist import load_blacklist_unbounded
        blacklist = load_blacklist_unbounded(args.blacklist_bed)
    excluded_main = [0]
    main_records = read_positions(
        main_path, "A", args.main_summit_column, args.main_score_column,
        blacklist=blacklist, excluded_counter=excluded_main, progress=reporter,
    )
    if not main_records:
        raise ValueError("The main BED contains no usable positions after filtering.")
    main_label = args.main_label or _default_label(main_path)
    prefix = Path(args.output_prefix or f"{_default_label(main_path)}_compare_positions")
    prefix.parent.mkdir(parents=True, exist_ok=True)

    results: list[ComparisonArrays] = []
    all_pairs_for_writing: list[tuple[ComparisonArrays, list[MatchedPair]]] = []
    summary_rows: list[dict[str, object]] = []
    histogram_rows: list[dict[str, object]] = []
    correlation_rows: list[dict[str, object]] = []
    percentile_rows: list[dict[str, object]] = []
    score_distance_rows: list[dict[str, object]] = []
    outputs: dict[str, Path] = {}
    bounds = _parse_distance_bins(args.distance_bins)

    for index, spec in enumerate(specs, start=1):
        reporter.stage(f"Comparison {index}/{len(specs)}: {spec.label}")
        excluded_compare = [0]
        compare_records = read_positions(
            spec.path, "B", args.compare_summit_column, args.compare_score_column,
            blacklist=blacklist, excluded_counter=excluded_compare, progress=reporter,
        )
        if not compare_records:
            raise ValueError(f"Comparison BED {spec.path} contains no usable positions after filtering.")
        reporter.stage(
            f"One-to-one matching {spec.label}: main={len(main_records):,}; comparison={len(compare_records):,}"
        )
        match = match_positions(main_records, compare_records, "one-to-one", args.max_distance, progress=reporter)
        if not match.pairs:
            raise ValueError(f"No matched pairs were found for comparison {spec.label!r}.")
        arrays = _arrays_from_match(spec.label, spec.path, match, len(main_records), len(compare_records), args.percentile_interval)
        results.append(arrays)
        summary = _summary_row(arrays, main_path, main_label)
        summary["blacklist_overlapping_main_records_excluded"] = excluded_main[0]
        summary["blacklist_overlapping_compare_records_excluded"] = excluded_compare[0]
        summary_rows.append(summary)
        histogram_rows.extend(_histogram_rows(arrays, args.histogram_x_min, args.histogram_x_max, args.histogram_bin_width))
        correlation_rows.extend(_distance_bin_rows(arrays, args.score_normalization, bounds))
        percentile_rows.extend(_percentile_rows(arrays))
        sd_stats = _score_distance_stats(arrays, args.score_distance_type)
        score_distance_rows.append(sd_stats)

        if not args.skip_pairs_tsv:
            pair_path = Path(f"{prefix}_{_safe_token(spec.label)}_pairs.tsv")
            _write_pairs(pair_path, arrays, match.pairs, main_label, args)
            outputs[f"pairs_{_safe_token(spec.label)}"] = pair_path
        # Comparison-specific plots remain separate where overlaying millions of
        # points would obscure rather than clarify the relationship.
        outputs[f"score_agreement_plot_{_safe_token(spec.label)}"] = _plot_score_agreement(prefix, arrays, main_label, args)
        outputs[f"score_distance_plot_{_safe_token(spec.label)}"] = _plot_score_distance(prefix, arrays, args, sd_stats)

    summary_path = Path(f"{prefix}_summary.tsv")
    histogram_path = Path(f"{prefix}_distance_histogram.tsv")
    correlation_path = Path(f"{prefix}_correlation_by_distance.tsv")
    percentile_path = Path(f"{prefix}_percentile_distances.tsv")
    percentile_summary_path = Path(f"{prefix}_percentile_summary.tsv")
    score_distance_path = Path(f"{prefix}_main_score_vs_distance_statistics.tsv")
    _write_tsv(summary_path, summary_rows)
    _write_tsv(histogram_path, histogram_rows)
    _write_tsv(correlation_path, correlation_rows)
    _write_tsv(percentile_path, percentile_rows)
    _write_tsv(percentile_summary_path, _percentile_summary_rows(results, args.percentile_interval))
    _write_tsv(score_distance_path, score_distance_rows)
    outputs.update({
        "summary": summary_path,
        "distance_histogram": histogram_path,
        "correlation_by_distance": correlation_path,
        "percentile_distances": percentile_path,
        "percentile_summary": percentile_summary_path,
        "score_distance_statistics": score_distance_path,
    })

    statistics_rows: list[dict[str, object]] = []
    if args.stats and len(results) >= 2:
        statistics_rows = _pairwise_statistics(results, args.percentile_interval, args.stats_test, args.p_adjust)
        statistics_path = Path(f"{prefix}_percentile_statistics.tsv")
        _write_tsv(statistics_path, statistics_rows)
        outputs["percentile_statistics"] = statistics_path

    outputs["distance_histogram_plot"] = _plot_distance_histograms(prefix, results, args)
    outputs["correlation_by_distance_plot"] = _plot_correlation_by_distance(prefix, correlation_rows, results, args)
    outputs["percentile_boxplot"] = _plot_percentile_boxplot(prefix, results, args.percentile_interval, statistics_rows, args)

    if not args.quiet:
        reporter.emit(f"Main positions: {len(main_records):,}")
        for result in results:
            reporter.emit(f"{result.label}: {result.compare_count:,} positions; {result.pair_count:,} one-to-one pairs")
        for name, path in outputs.items():
            reporter.emit(f"{name}: {path}")
    reporter.complete()
    return outputs


def _run_serial(args: argparse.Namespace) -> int:
    try:
        run_comparison(args)
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


def run(args: argparse.Namespace) -> int:
    from nucleosuite.output_naming import parameterized_prefix
    from nucleosuite.partitioned import run_partitioned_command

    main_path, _specs = _resolve_inputs(args)
    requested = args.output_prefix or f"{_default_label(main_path)}_compare_positions"
    args.output_prefix = str(
        parameterized_prefix(
            requested,
            (
                ("pct", args.percentile_interval),
                ("maxdist", args.max_distance),
                ("dist", args.score_distance_type),
            ),
        )
    )
    base = Path(args.output_prefix).name
    return run_partitioned_command(
        "compare-positions",
        args,
        _run_serial,
        runner_module="nucleosuite.compare_positions",
        runner_function="_run_serial",
        primary_attr="main_bed",
        output_prefix_attr="output_prefix",
        path_attrs=("blacklist_bed",),
        named_path_list_attrs=("compare_beds",),
        base_name=base,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
