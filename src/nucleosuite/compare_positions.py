#!/usr/bin/env python3
"""Compare one main nucleosome callset with one or more comparison callsets."""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import itertools
import math
import re
import sys
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
from scipy import stats

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.regions import canonical_contig_key
from nucleosuite.io import open_text as open_interval_text
from nucleosuite.progress import ProgressReporter


@dataclass(frozen=True)
class PeakChrom:
    """Compact summit/score storage for one chromosome.

    Arrays are sorted by summit and then source line number.  Chromosome names
    are stored once at the dictionary level rather than repeated per peak.
    """

    name: str
    summits: np.ndarray
    scores: np.ndarray
    indices: np.ndarray

    @property
    def count(self) -> int:
        return int(self.summits.size)


@dataclass(frozen=True)
class CompactPeakSet:
    """One BED represented as chromosome-keyed numeric arrays."""

    path: Path
    source: str
    by_chrom: dict[str, PeakChrom]
    count: int
    excluded_blacklist: int = 0


@dataclass
class _CompactNode:
    """One unmatched coordinate group in the one-to-one matching frontier."""

    kind: str
    coordinate: int
    low: int
    high: int  # exclusive
    previous: int | None = None
    following: int | None = None
    version: int = 0
    alive: bool = True

    @property
    def empty(self) -> bool:
        return self.low >= self.high


@dataclass(frozen=True)
class CompactMatch:
    """Compact matched-pair arrays and one-to-one matching diagnostics."""

    query_source: str
    query_count: int
    target_count: int
    unmatched_no_target_chrom: int
    unmatched_distance: int
    unmatched_unique: int
    chrom_names: tuple[str, ...]
    chrom_code: np.ndarray
    main_index: np.ndarray
    compare_index: np.ndarray
    main_summit: np.ndarray
    main_scores: np.ndarray
    compare_scores: np.ndarray
    signed_distance: np.ndarray

    @property
    def pair_count(self) -> int:
        return int(self.main_scores.size)


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
    chrom_names: tuple[str, ...]
    chrom_code: np.ndarray
    main_line: np.ndarray
    compare_line: np.ndarray
    main_summit: np.ndarray
    main_scores: np.ndarray
    compare_scores: np.ndarray
    signed_distance: np.ndarray
    absolute_distance: np.ndarray
    percentile: np.ndarray
    group_index: np.ndarray
    group_names: tuple[str, ...]

    @property
    def pair_count(self) -> int:
        return int(self.main_scores.size)

    def group_mask(self, label: str) -> np.ndarray:
        try:
            index = self.group_names.index(label)
        except ValueError:
            return np.zeros(self.pair_count, dtype=bool)
        return self.group_index == index

    def group_label(self, pair_index: int) -> str:
        return self.group_names[int(self.group_index[pair_index])]



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
        help="Main nucleosome BED, BED.gz, or bigBed file. Use LABEL=BED to set the main callset label.",
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
        "--histogram-x-min", type=float, default=-250.0,
        help="Displayed lower x-axis limit for signed distance distributions (default: -250).",
    )
    parser.add_argument(
        "--histogram-x-max", type=float, default=250.0,
        help="Displayed upper x-axis limit for signed distance distributions (default: 250).",
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
        "--percentile-boxplot-y-max", type=float, default=200.0,
        help="Displayed upper y-axis limit for percentile distance boxplots; 0 disables (default: 200).",
    )
    outliers = parser.add_mutually_exclusive_group()
    outliers.add_argument(
        "--show-boxplot-outliers", dest="show_boxplot_outliers", action="store_true", default=True,
        help="Show boxplot outlier points beyond the 1.5×IQR whiskers (default: shown).",
    )
    outliers.add_argument(
        "--hide-boxplot-outliers", dest="show_boxplot_outliers", action="store_false",
        help="Hide boxplot outlier points while retaining the standard 1.5×IQR whiskers.",
    )
    parser.add_argument(
        "--score-agreement-distance-max", type=float, default=100.0,
        help="Upper colour-scale limit for absolute summit distance in score-agreement plots; 0 uses the data maximum (default: 100).",
    )
    parser.add_argument(
        "--score-distance-y-max", type=float, default=0.0,
        help="Optional displayed upper y-axis limit for absolute main-score-versus-distance plots; 0 uses the data range (default: 0).",
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
        "--write-detail-tables", action="store_true",
        help=(
            "Write large record-level supporting tables: one matched-pair TSV per "
            "comparison and the combined percentile-distance TSV. These files are "
            "omitted by default."
        ),
    )
    parser.add_argument(
        "--skip-pairs-tsv", action="store_true", help=argparse.SUPPRESS,
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

    main_text = str(main_value)
    main_spec = _parse_compare_spec(main_text)
    # Main and comparison callsets use the same optional LABEL=BED syntax.
    # Preserve an explicitly resolved main label when partition/serial routing
    # normalizes --main-bed from LABEL=BED to the underlying path.
    if "=" in main_text:
        args._main_label = main_spec.label
    elif not getattr(args, "_main_label", None):
        args._main_label = main_spec.label
    args.main_bed = str(main_spec.path)

    specs = [_parse_compare_spec(str(value)) for value in compare_values]
    labels = [spec.label for spec in specs]
    if len(labels) != len(set(labels)):
        duplicates = sorted({label for label in labels if labels.count(label) > 1})
        raise ValueError("Comparison labels must be unique: " + ", ".join(duplicates))
    return main_spec.path, specs


def _validate_args(args: argparse.Namespace) -> None:
    if not 1 <= int(args.percentile_interval) <= 100:
        raise ValueError("--percentile-interval must be between 1 and 100.")
    if args.max_distance is not None and (not math.isfinite(float(args.max_distance)) or float(args.max_distance) < 0):
        raise ValueError("--max-distance must be finite and zero or greater.")
    if not math.isfinite(float(args.histogram_bin_width)) or float(args.histogram_bin_width) <= 0:
        raise ValueError("--histogram-bin-width must be greater than zero.")
    if not math.isfinite(float(args.histogram_x_min)) or not math.isfinite(float(args.histogram_x_max)) or float(args.histogram_x_max) <= float(args.histogram_x_min):
        raise ValueError("Histogram limits require finite x-min < x-max.")
    if int(args.plot_max_points) < 0:
        raise ValueError("--plot-max-points must be zero or greater.")
    if int(args.dpi) < 1:
        raise ValueError("--dpi must be at least 1.")
    if not math.isfinite(float(args.score_z_limit)) or float(args.score_z_limit) < 0:
        raise ValueError("--score-z-limit must be finite and zero or greater.")
    if not math.isfinite(float(args.percentile_boxplot_y_max)) or float(args.percentile_boxplot_y_max) < 0:
        raise ValueError("--percentile-boxplot-y-max must be finite and zero or greater.")
    if not math.isfinite(float(args.score_agreement_distance_max)) or float(args.score_agreement_distance_max) < 0:
        raise ValueError("--score-agreement-distance-max must be finite and zero or greater.")
    if not math.isfinite(float(args.score_distance_y_max)) or float(args.score_distance_y_max) < 0:
        raise ValueError("--score-distance-y-max must be finite and zero or greater.")
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


def _assign_matched_percentiles(main_scores: np.ndarray, interval: int) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    """Rank matched pairs by main score and return compact percentile groups."""

    n = int(main_scores.size)
    group_names = tuple(label for _low, _high, label in _percentile_group_bounds(interval))
    if n == 0:
        return (
            np.asarray([], dtype=np.uint8),
            np.asarray([], dtype=np.uint8),
            group_names,
        )
    order = np.argsort(main_scores, kind="stable")
    percentile = np.empty(n, dtype=np.uint8)
    ranks = np.arange(1, n + 1, dtype=float)
    percentile[order] = np.ceil(ranks * 100.0 / n).astype(np.uint8)
    group_index = np.empty(n, dtype=np.uint8)
    for index, (lower, upper, _label) in enumerate(_percentile_group_bounds(interval)):
        mask = (percentile > lower) & (percentile <= upper)
        group_index[mask] = index
    return percentile, group_index, group_names


def _parse_integer_coordinate(value: str, context: str) -> int:
    try:
        numeric = float(value)
    except ValueError as exc:
        raise ValueError(f"{context}: expected a numeric summit coordinate, found {value!r}.") from exc
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f"{context}: summit coordinates must be finite integers, found {value!r}.")
    return int(numeric)


def read_compact_positions(
    path: str | Path,
    source: str,
    summit_column: int | None,
    score_column: int,
    *,
    blacklist=None,
    progress: ProgressReporter | None = None,
) -> CompactPeakSet:
    """Read a BED-like position file into compact chromosome-keyed arrays.

    Only the resolved summit, numeric score, and source line number are retained
    for each usable peak. BED start/end coordinates are used transiently for
    validation, midpoint calculation, and blacklist filtering and are then
    discarded.
    """

    input_path = Path(path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input interval file not found: {input_path}")
    if summit_column is not None and int(summit_column) < 1:
        raise ValueError("Summit columns must be one-based and at least 1.")
    if int(score_column) < 1:
        raise ValueError("Score columns must be one-based and at least 1.")
    required_column = max(3, int(score_column), int(summit_column or 0))

    # array.array keeps the construction phase compact too; converting a BED
    # with millions of rows therefore does not first create millions of Python
    # record objects.
    builders: dict[str, tuple[str, array, array, array]] = {}
    excluded = 0
    count = 0
    if progress is not None:
        progress.file_start(source, input_path)
    seen_contigs: set[str] = set()
    with open_interval_text(input_path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            text = raw.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < required_column:
                raise ValueError(
                    f"{input_path}:{line_number}: requires at least {required_column} "
                    f"columns for the selected summit and score fields; found {len(fields)}."
                )
            chrom = fields[0]
            key = canonical_contig_key(chrom)
            if progress is not None and key not in seen_contigs:
                seen_contigs.add(key)
                progress.reading_contig(source, chrom)
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{input_path}:{line_number}: BED start and end must be integers.") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{input_path}:{line_number}: require 0 <= start < end.")
            if blacklist is not None and blacklist.overlaps(chrom, start, end):
                excluded += 1
                continue
            summit = (
                (start + end) // 2
                if summit_column is None
                else _parse_integer_coordinate(fields[int(summit_column) - 1], f"{input_path}:{line_number}")
            )
            try:
                score = float(fields[int(score_column) - 1])
            except ValueError as exc:
                raise ValueError(
                    f"{input_path}:{line_number}: score column {score_column} must be numeric."
                ) from exc
            if not math.isfinite(score):
                raise ValueError(f"{input_path}:{line_number}: score column {score_column} must be finite.")
            if key not in builders:
                builders[key] = (chrom, array("q"), array("d"), array("q"))
            _name, summits, scores, indices = builders[key]
            summits.append(int(summit))
            scores.append(float(score))
            indices.append(int(line_number))
            count += 1

    if count == 0:
        raise ValueError(f"No interval records were read from {input_path}.")

    by_chrom: dict[str, PeakChrom] = {}
    for key, (name, summits_buffer, scores_buffer, indices_buffer) in builders.items():
        summits = np.asarray(summits_buffer, dtype=np.int64)
        scores = np.asarray(scores_buffer, dtype=np.float64)
        indices = np.asarray(indices_buffer, dtype=np.int64)
        order = np.lexsort((indices, summits))
        by_chrom[key] = PeakChrom(
            name=name,
            summits=np.ascontiguousarray(summits[order]),
            scores=np.ascontiguousarray(scores[order]),
            indices=np.ascontiguousarray(indices[order]),
        )
    if progress is not None:
        progress.file_complete(source, input_path, count)
    return CompactPeakSet(input_path, source, by_chrom, count, excluded)


def _coordinate_groups(values: np.ndarray) -> list[tuple[int, int, int]]:
    """Return ``(coordinate, start, end)`` groups for one sorted summit array."""

    if values.size == 0:
        return []
    boundaries = np.flatnonzero(np.diff(values) != 0) + 1
    starts = np.concatenate((np.asarray([0], dtype=np.int64), boundaries))
    ends = np.concatenate((boundaries, np.asarray([values.size], dtype=np.int64)))
    return [(int(values[start]), int(start), int(end)) for start, end in zip(starts, ends)]


def _match_chrom_one_to_one(
    query: PeakChrom,
    target: PeakChrom,
    max_distance: float | None,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Distance-prioritized one-to-one matching for a single chromosome.

    Only compact numeric arrays persist between chromosomes. Temporary frontier
    nodes and the heap are released after each chromosome has been matched.
    """

    q_groups = _coordinate_groups(query.summits)
    t_groups = _coordinate_groups(target.summits)
    exact_q: list[np.ndarray] = []
    exact_t: list[np.ndarray] = []
    nodes: list[_CompactNode] = []
    qi = ti = 0
    while qi < len(q_groups) or ti < len(t_groups):
        if ti >= len(t_groups) or (qi < len(q_groups) and q_groups[qi][0] < t_groups[ti][0]):
            coord, low, high = q_groups[qi]
            nodes.append(_CompactNode("Q", coord, low, high))
            qi += 1
            continue
        if qi >= len(q_groups) or t_groups[ti][0] < q_groups[qi][0]:
            coord, low, high = t_groups[ti]
            nodes.append(_CompactNode("T", coord, low, high))
            ti += 1
            continue
        coord, q_low, q_high = q_groups[qi]
        _coord_t, t_low, t_high = t_groups[ti]
        match_count = min(q_high - q_low, t_high - t_low)
        if match_count:
            exact_q.append(np.arange(q_low, q_low + match_count, dtype=np.int64))
            exact_t.append(np.arange(t_low, t_low + match_count, dtype=np.int64))
        q_low += match_count
        t_low += match_count
        if q_low < q_high:
            nodes.append(_CompactNode("Q", coord, q_low, q_high))
        elif t_low < t_high:
            nodes.append(_CompactNode("T", coord, t_low, t_high))
        qi += 1
        ti += 1

    for index, node in enumerate(nodes):
        node.previous = index - 1 if index else None
        node.following = index + 1 if index + 1 < len(nodes) else None

    heap: list[tuple[int, int, int, int, int, int, int]] = []

    def candidate_indices(left: _CompactNode, right: _CompactNode) -> tuple[int, int]:
        if left.kind == "Q":
            return left.low, right.low
        return right.low, left.high - 1

    def push_boundary(left_index: int | None) -> None:
        if left_index is None:
            return
        left = nodes[left_index]
        right_index = left.following
        if not left.alive or right_index is None:
            return
        right = nodes[right_index]
        if not right.alive or left.kind == right.kind:
            return
        query_index, target_index = candidate_indices(left, right)
        heapq.heappush(
            heap,
            (
                abs(right.coordinate - left.coordinate),
                query_index,
                target_index,
                left_index,
                right_index,
                left.version,
                right.version,
            ),
        )

    for index in range(max(0, len(nodes) - 1)):
        push_boundary(index)

    def unlink(index: int) -> None:
        node = nodes[index]
        previous = node.previous
        following = node.following
        if previous is not None:
            nodes[previous].following = following
        if following is not None:
            nodes[following].previous = previous
        node.alive = False
        node.previous = node.following = None
        node.version += 1

    matched_q: list[int] = []
    matched_t: list[int] = []
    while heap:
        distance, q_index, t_index, left_index, right_index, left_version, right_version = heapq.heappop(heap)
        left = nodes[left_index]
        right = nodes[right_index]
        if (
            not left.alive
            or not right.alive
            or left.following != right_index
            or right.previous != left_index
            or left.version != left_version
            or right.version != right_version
        ):
            continue
        if max_distance is not None and distance > max_distance:
            break
        current_q, current_t = candidate_indices(left, right)
        if current_q != q_index or current_t != t_index:
            continue
        matched_q.append(q_index)
        matched_t.append(t_index)
        affected = {left.previous, left_index, right_index, right.following}
        if left.kind == "Q":
            left.low += 1
            right.low += 1
        else:
            left.high -= 1
            right.low += 1
        left.version += 1
        right.version += 1
        if left.empty:
            unlink(left_index)
        if right.empty:
            unlink(right_index)
        for candidate in list(affected):
            if candidate is None or not nodes[candidate].alive:
                continue
            push_boundary(nodes[candidate].previous)
            push_boundary(candidate)

    remaining_query: list[int] = []
    for node in nodes:
        if node.alive and node.kind == "Q" and not node.empty:
            remaining_query.extend(range(node.low, node.high))
    unmatched_distance = 0
    unmatched_unique = 0
    if max_distance is None:
        unmatched_unique = len(remaining_query)
    elif remaining_query:
        first_target = int(target.summits[0])
        last_target = int(target.summits[-1])
        for q_index in remaining_query:
            summit = int(query.summits[q_index])
            farthest = max(abs(first_target - summit), abs(last_target - summit))
            if farthest > max_distance:
                unmatched_distance += 1
            else:
                unmatched_unique += 1

    q_parts = list(exact_q)
    t_parts = list(exact_t)
    if matched_q:
        q_parts.append(np.asarray(matched_q, dtype=np.int64))
        t_parts.append(np.asarray(matched_t, dtype=np.int64))
    q_out = np.concatenate(q_parts) if q_parts else np.asarray([], dtype=np.int64)
    t_out = np.concatenate(t_parts) if t_parts else np.asarray([], dtype=np.int64)
    return q_out, t_out, unmatched_distance, unmatched_unique


def match_compact_positions(
    main: CompactPeakSet,
    comparison: CompactPeakSet,
    max_distance: float | None,
    *,
    progress: ProgressReporter | None = None,
    progress_stage: str = "Matching one-to-one",
) -> CompactMatch:
    """Match the smaller callset against the larger callset one-to-one."""

    if main.count <= comparison.count:
        query_set, target_set, query_source = main, comparison, "main"
    else:
        query_set, target_set, query_source = comparison, main, "comparison"

    main_chrom_keys = list(main.by_chrom)
    chrom_to_code = {key: index for index, key in enumerate(main_chrom_keys)}
    chrom_names = tuple(main.by_chrom[key].name for key in main_chrom_keys)
    chrom_parts: list[np.ndarray] = []
    main_index_parts: list[np.ndarray] = []
    compare_index_parts: list[np.ndarray] = []
    main_summit_parts: list[np.ndarray] = []
    main_score_parts: list[np.ndarray] = []
    compare_score_parts: list[np.ndarray] = []
    signed_parts: list[np.ndarray] = []
    unmatched_no_target_chrom = 0
    unmatched_distance = 0
    unmatched_unique = 0

    query_items = list(query_set.by_chrom.items())
    for chrom_number, (key, qchrom) in enumerate(query_items, start=1):
        if progress is not None:
            progress.contig(progress_stage, qchrom.name, chrom_number, len(query_items), qchrom.count)
        tchrom = target_set.by_chrom.get(key)
        if tchrom is None:
            unmatched_no_target_chrom += qchrom.count
            continue
        q_idx, t_idx, rejected_distance, rejected_unique = _match_chrom_one_to_one(qchrom, tchrom, max_distance)
        unmatched_distance += rejected_distance
        unmatched_unique += rejected_unique
        if q_idx.size == 0:
            continue
        if query_source == "main":
            mch, cch = qchrom, tchrom
            mi, ci = q_idx, t_idx
        else:
            mch, cch = tchrom, qchrom
            mi, ci = t_idx, q_idx
        code = chrom_to_code[key]
        chrom_parts.append(np.full(mi.size, code, dtype=np.uint16 if len(chrom_names) <= 65535 else np.uint32))
        main_index_parts.append(mch.indices[mi])
        compare_index_parts.append(cch.indices[ci])
        main_summit_parts.append(mch.summits[mi])
        main_score_parts.append(mch.scores[mi])
        compare_score_parts.append(cch.scores[ci])
        signed_parts.append(cch.summits[ci] - mch.summits[mi])

    def concat(parts: list[np.ndarray], dtype) -> np.ndarray:
        return np.concatenate(parts).astype(dtype, copy=False) if parts else np.asarray([], dtype=dtype)

    chrom_code = concat(chrom_parts, np.uint16 if len(chrom_names) <= 65535 else np.uint32)
    main_index = concat(main_index_parts, np.int64)
    compare_index = concat(compare_index_parts, np.int64)
    main_summit = concat(main_summit_parts, np.int64)
    main_scores = concat(main_score_parts, np.float64)
    compare_scores = concat(compare_score_parts, np.float64)
    signed_distance = concat(signed_parts, np.int64)
    if main_scores.size:
        order = np.lexsort((main_index, main_summit, chrom_code))
        chrom_code = chrom_code[order]
        main_index = main_index[order]
        compare_index = compare_index[order]
        main_summit = main_summit[order]
        main_scores = main_scores[order]
        compare_scores = compare_scores[order]
        signed_distance = signed_distance[order]

    return CompactMatch(
        query_source=query_source,
        query_count=query_set.count,
        target_count=target_set.count,
        unmatched_no_target_chrom=unmatched_no_target_chrom,
        unmatched_distance=unmatched_distance,
        unmatched_unique=unmatched_unique,
        chrom_names=chrom_names,
        chrom_code=chrom_code,
        main_index=main_index,
        compare_index=compare_index,
        main_summit=main_summit,
        main_scores=main_scores,
        compare_scores=compare_scores,
        signed_distance=signed_distance,
    )


def _arrays_from_match(
    label: str,
    path: Path,
    result: CompactMatch,
    main_count: int,
    compare_count: int,
    interval: int,
) -> ComparisonArrays:
    absolute = np.abs(result.signed_distance).astype(np.float64)
    percentile, group_index, group_names = _assign_matched_percentiles(result.main_scores, interval)
    return ComparisonArrays(
        label=label,
        path=path,
        main_count=main_count,
        compare_count=compare_count,
        query_source=result.query_source,
        query_count=result.query_count,
        target_count=result.target_count,
        unmatched_no_target_chrom=result.unmatched_no_target_chrom,
        unmatched_distance=result.unmatched_distance,
        unmatched_unique=result.unmatched_unique,
        chrom_names=result.chrom_names,
        chrom_code=result.chrom_code,
        main_line=result.main_index,
        compare_line=result.compare_index,
        main_summit=result.main_summit,
        main_scores=result.main_scores,
        compare_scores=result.compare_scores,
        signed_distance=result.signed_distance,
        absolute_distance=absolute,
        percentile=percentile,
        group_index=group_index,
        group_names=group_names,
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
    counts, edges = np.histogram(result.signed_distance.astype(float), bins=edges)
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
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def _write_pairs(path: Path, result: ComparisonArrays, main_label: str, args: argparse.Namespace) -> Path:
    """Write matched pairs from compact arrays without recreating peak objects."""

    selected_main, selected_compare, _ = _selected_scores(result, args.score_normalization)
    selected_indices = _plot_indices(result.pair_count, args.plot_max_points, args.plot_seed)
    plot_selected = np.zeros(result.pair_count, dtype=np.uint8)
    plot_selected[selected_indices] = 1
    fields = [
        "main_label", "comparison", "pair_id", "query_source", "chrom",
        "main_summit", "main_score", "main_score_normalized", "main_line_number",
        "compare_summit", "compare_score", "compare_score_normalized", "compare_line_number",
        "signed_distance_compare_minus_main", "absolute_distance", "main_score_percentile", "percentile_group",
        "plot_selected", "plot_score_normalization", "plot_score_z_limit", "plot_correlation_method",
        "plot_score_agreement_distance_max", "plot_score_distance_y_max", "plot_label_a", "plot_label_b",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index in range(result.pair_count):
            signed = int(result.signed_distance[index])
            main_summit = int(result.main_summit[index])
            writer.writerow({
                "main_label": main_label,
                "comparison": result.label,
                "pair_id": f"pair_{index + 1:09d}",
                "query_source": result.query_source,
                "chrom": result.chrom_names[int(result.chrom_code[index])],
                "main_summit": main_summit,
                "main_score": float(result.main_scores[index]),
                "main_score_normalized": float(selected_main[index]),
                "main_line_number": int(result.main_line[index]),
                "compare_summit": main_summit + signed,
                "compare_score": float(result.compare_scores[index]),
                "compare_score_normalized": float(selected_compare[index]),
                "compare_line_number": int(result.compare_line[index]),
                "signed_distance_compare_minus_main": signed,
                "absolute_distance": float(result.absolute_distance[index]),
                "main_score_percentile": int(result.percentile[index]),
                "percentile_group": result.group_label(index),
                "plot_selected": int(plot_selected[index]),
                "plot_score_normalization": args.score_normalization,
                "plot_score_z_limit": args.score_z_limit,
                "plot_correlation_method": args.score_correlation,
                "plot_score_agreement_distance_max": args.score_agreement_distance_max,
                "plot_score_distance_y_max": args.score_distance_y_max,
                "plot_label_a": main_label,
                "plot_label_b": result.label,
            })
    return path



def _append_percentile_rows(path: Path, result: ComparisonArrays, main_label: str, *, write_header: bool) -> None:
    """Stream matched percentile rows to disk instead of retaining Python dicts."""

    fields = [
        "main_label", "comparison", "main_line_number", "main_score", "main_score_percentile",
        "percentile_group", "signed_distance_compare_minus_main", "absolute_distance",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if write_header else "a"
    with path.open(mode, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        if write_header:
            writer.writeheader()
        for index in range(result.pair_count):
            writer.writerow({
                "main_label": main_label,
                "comparison": result.label,
                "main_line_number": int(result.main_line[index]),
                "main_score": float(result.main_scores[index]),
                "main_score_percentile": int(result.percentile[index]),
                "percentile_group": result.group_label(index),
                "signed_distance_compare_minus_main": int(result.signed_distance[index]),
                "absolute_distance": float(result.absolute_distance[index]),
            })

def _percentile_summary_rows(results: Sequence[ComparisonArrays], interval: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for lower, upper, label in _percentile_group_bounds(interval):
            mask = result.group_mask(label)
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


def _boxplot_stats(values: np.ndarray) -> tuple[dict[str, float], np.ndarray]:
    """Return Matplotlib-compatible 1.5×IQR box statistics and fliers."""

    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        empty = {
            "minimum_absolute_distance": math.nan,
            "q1_absolute_distance": math.nan,
            "median_absolute_distance": math.nan,
            "mean_absolute_distance": math.nan,
            "q3_absolute_distance": math.nan,
            "maximum_absolute_distance": math.nan,
            "whisker_low_absolute_distance": math.nan,
            "whisker_high_absolute_distance": math.nan,
        }
        return empty, np.asarray([], dtype=float)

    finite.sort()
    q1, median, q3 = np.percentile(finite, [25, 50, 75])
    iqr = float(q3 - q1)
    lower_fence = float(q1 - 1.5 * iqr)
    upper_fence = float(q3 + 1.5 * iqr)
    within = finite[(finite >= lower_fence) & (finite <= upper_fence)]
    whisker_low = float(within[0]) if within.size else float(q1)
    whisker_high = float(within[-1]) if within.size else float(q3)
    fliers = finite[(finite < whisker_low) | (finite > whisker_high)]
    stats = {
        "minimum_absolute_distance": float(finite[0]),
        "q1_absolute_distance": float(q1),
        "median_absolute_distance": float(median),
        "mean_absolute_distance": float(np.mean(finite)),
        "q3_absolute_distance": float(q3),
        "maximum_absolute_distance": float(finite[-1]),
        "whisker_low_absolute_distance": whisker_low,
        "whisker_high_absolute_distance": whisker_high,
    }
    return stats, fliers


def _percentile_boxplot_source_rows(
    results: Sequence[ComparisonArrays],
    main_label: str,
    interval: int,
    statistics_rows: Sequence[dict[str, object]],
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    """Build a compact, self-contained source table for the percentile boxplot.

    The table contains one ``box`` row per comparison/percentile group, only
    the outlier observations required to redraw fliers, and optional ``stat``
    rows for pairwise annotations.  It therefore supports faithful replotting
    without retaining every matched-pair distance.
    """

    fields = (
        "row_type", "main_label", "comparison", "percentile_group",
        "group_lower_percentile", "group_upper_percentile", "matched_pair_count",
        "minimum_absolute_distance", "q1_absolute_distance",
        "median_absolute_distance", "mean_absolute_distance",
        "q3_absolute_distance", "maximum_absolute_distance",
        "whisker_low_absolute_distance", "whisker_high_absolute_distance",
        "outlier_absolute_distance", "comparison_1", "comparison_2", "test",
        "paired", "n_1", "n_2", "n_paired", "statistic", "p_value",
        "p_adjusted", "p_adjustment", "significance", "p_display",
        "show_boxplot_outliers", "percentile_boxplot_y_max",
    )

    rows: list[dict[str, object]] = []
    for result in results:
        for lower, upper, label in _percentile_group_bounds(interval):
            values = result.absolute_distance[result.group_mask(label)]
            box_stats, fliers = _boxplot_stats(values)
            box_row = {key: "" for key in fields}
            box_row.update({
                "row_type": "box",
                "main_label": main_label,
                "comparison": result.label,
                "percentile_group": label,
                "group_lower_percentile": lower,
                "group_upper_percentile": upper,
                "matched_pair_count": int(values.size),
                "show_boxplot_outliers": int(bool(args.show_boxplot_outliers)),
                "percentile_boxplot_y_max": float(args.percentile_boxplot_y_max),
            })
            box_row.update(box_stats)
            rows.append(box_row)
            for value in fliers:
                flier_row = {key: "" for key in fields}
                flier_row.update({
                    "row_type": "flier",
                    "main_label": main_label,
                    "comparison": result.label,
                    "percentile_group": label,
                    "outlier_absolute_distance": float(value),
                    "show_boxplot_outliers": int(bool(args.show_boxplot_outliers)),
                    "percentile_boxplot_y_max": float(args.percentile_boxplot_y_max),
                })
                rows.append(flier_row)

    for source in statistics_rows:
        stat_row = {key: "" for key in fields}
        stat_row.update({
            "row_type": "stat",
            "main_label": main_label,
            "percentile_group": source.get("percentile_group", ""),
            "comparison_1": source.get("comparison_1", ""),
            "comparison_2": source.get("comparison_2", ""),
            "test": source.get("test", ""),
            "paired": source.get("paired", ""),
            "n_1": source.get("n_1", ""),
            "n_2": source.get("n_2", ""),
            "n_paired": source.get("n_paired", ""),
            "statistic": source.get("statistic", ""),
            "p_value": source.get("p_value", ""),
            "p_adjusted": source.get("p_adjusted", ""),
            "p_adjustment": source.get("p_adjustment", ""),
            "significance": source.get("significance", ""),
            "p_display": args.p_display,
            "show_boxplot_outliers": int(bool(args.show_boxplot_outliers)),
            "percentile_boxplot_y_max": float(args.percentile_boxplot_y_max),
        })
        rows.append(stat_row)
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
    """Return paired distances without constructing per-peak Python dictionaries."""

    mask_a = a.group_mask(group)
    mask_b = b.group_mask(group)
    lines_a = a.main_line[mask_a]
    lines_b = b.main_line[mask_b]
    if lines_a.size == 0 or lines_b.size == 0:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    _shared, index_a, index_b = np.intersect1d(
        lines_a, lines_b, assume_unique=True, return_indices=True
    )
    values_a = a.absolute_distance[mask_a]
    values_b = b.absolute_distance[mask_b]
    return values_a[index_a].astype(float, copy=False), values_b[index_b].astype(float, copy=False)


def _pairwise_statistics(results: Sequence[ComparisonArrays], interval: int, family: str, adjustment: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for _lower, _upper, group in _percentile_group_bounds(interval):
        group_rows: list[dict[str, object]] = []
        for first, second in itertools.combinations(results, 2):
            values_first = first.absolute_distance[first.group_mask(group)]
            values_second = second.absolute_distance[second.group_mask(group)]
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


def _write_score_agreement_plot_source(
    path: Path,
    result: ComparisonArrays,
    main_label: str,
    args: argparse.Namespace,
) -> Path:
    """Write only the sampled points and full-data statistics needed for replotting."""

    x_all, y_all, axis_label = _selected_scores(result, args.score_normalization)
    indices = _plot_indices(result.pair_count, args.plot_max_points, args.plot_seed)
    spearman, spearman_p = _safe_corr(x_all, y_all, "spearman")
    pearson, pearson_p = _safe_corr(x_all, y_all, "pearson")
    _slope, _intercept, r2, _regression_p = _linear_stats(x_all, y_all)
    rows: list[dict[str, object]] = []
    for index in indices:
        rows.append({
            "main_label": main_label,
            "comparison": result.label,
            "main_score_plot": float(x_all[index]),
            "compare_score_plot": float(y_all[index]),
            "absolute_distance": float(result.absolute_distance[index]),
            "score_axis_label": axis_label,
            "plot_score_normalization": args.score_normalization,
            "plot_score_z_limit": args.score_z_limit,
            "plot_correlation_method": args.score_correlation,
            "plot_score_agreement_distance_max": args.score_agreement_distance_max,
            "full_matched_pair_count": result.pair_count,
            "full_spearman": spearman,
            "full_spearman_p_value": spearman_p,
            "full_pearson": pearson,
            "full_pearson_p_value": pearson_p,
            "full_linear_r_squared": r2,
        })
    return _write_tsv(path, rows)


def _write_score_distance_plot_source(
    path: Path,
    result: ComparisonArrays,
    main_label: str,
    args: argparse.Namespace,
    stats_row: Mapping[str, object],
) -> Path:
    """Write the sampled score-distance points plus full-data fit statistics."""

    y_all = (
        result.absolute_distance
        if args.score_distance_type == "absolute"
        else result.signed_distance.astype(float)
    )
    indices = _plot_indices(result.pair_count, args.plot_max_points, args.plot_seed)
    rows: list[dict[str, object]] = []
    for index in indices:
        rows.append({
            "main_label": main_label,
            "comparison": result.label,
            "main_score": float(result.main_scores[index]),
            "matched_distance": float(y_all[index]),
            "distance_type": args.score_distance_type,
            "plot_type": args.score_distance_plot,
            "plot_score_distance_y_max": args.score_distance_y_max,
            "full_matched_pair_count": result.pair_count,
            "full_spearman_rho": stats_row.get("spearman_rho", math.nan),
            "full_spearman_p_value": stats_row.get("spearman_p_value", math.nan),
            "full_pearson_r": stats_row.get("pearson_r", math.nan),
            "full_pearson_p_value": stats_row.get("pearson_p_value", math.nan),
            "full_linear_slope": stats_row.get("linear_slope_bp_per_score_unit", math.nan),
            "full_linear_intercept": stats_row.get("linear_intercept", math.nan),
            "full_linear_r_squared": stats_row.get("linear_r_squared", math.nan),
            "full_linear_slope_p_value": stats_row.get("linear_slope_p_value", math.nan),
        })
    return _write_tsv(path, rows)


def _attach_plot_source(plot_path_value: Path, source_table: Path, plot_type: str) -> None:
    """Record the compact table that can recreate an original command plot."""

    from nucleosuite.plotting import write_plot_metadata
    write_plot_metadata(
        plot_path_value,
        extra={"source_table": str(source_table), "detected_plot_type": plot_type},
    )


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
    color_kwargs: dict[str, object] = {}
    if float(args.score_agreement_distance_max) > 0:
        from matplotlib.colors import Normalize
        color_kwargs["norm"] = Normalize(
            vmin=0.0, vmax=float(args.score_agreement_distance_max), clip=True
        )
    scatter = axis.scatter(
        x[indices], y[indices], c=result.absolute_distance[indices], s=9,
        alpha=0.55, linewidths=0, cmap="viridis", rasterized=True, **color_kwargs,
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


def _plot_score_distance(prefix: Path, result: ComparisonArrays, main_label: str, args: argparse.Namespace, stats_row: dict[str, object]) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import plot_path, save_figure

    y_all = result.absolute_distance if args.score_distance_type == "absolute" else result.signed_distance.astype(float)
    indices = _plot_indices(result.pair_count, args.plot_max_points, args.plot_seed)
    x = result.main_scores[indices]
    y = y_all[indices]

    # If a command-specific or shared y-range is explicitly requested, restrict
    # plotting data before hexbin construction. This prevents bins spanning a
    # much larger hidden distance range from appearing as full-height bands.
    from nucleosuite.plotting import get_plot_options
    plot_options = get_plot_options()
    lower = plot_options.y_min
    upper = plot_options.y_max
    if args.score_distance_type == "absolute":
        if lower is None:
            lower = 0.0
        if float(args.score_distance_y_max) > 0 and upper is None:
            upper = float(args.score_distance_y_max)
    mask = np.isfinite(x) & np.isfinite(y)
    if lower is not None:
        mask &= y >= float(lower)
    if upper is not None:
        mask &= y <= float(upper)
    x = x[mask]
    y = y[mask]

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
    axis.set_xlabel(f"{main_label} peak score")
    axis.set_ylabel("Absolute matched distance (bp)" if args.score_distance_type == "absolute" else f"Signed {result.label} − {main_label} distance (bp)")
    if lower is not None or upper is not None:
        bottom, top = axis.get_ylim()
        axis.set_ylim(bottom if lower is None else float(lower), top if upper is None else float(upper))
    axis.set_title(f"{main_label} peak score versus distance: {result.label}")
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

    edges = np.arange(
        float(args.histogram_x_min),
        float(args.histogram_x_max) + float(args.histogram_bin_width),
        float(args.histogram_bin_width),
    )
    colors = _comparison_colors(len(results))
    figure, axis = plt.subplots(figsize=(9.0, 5.8))
    for color, result in zip(colors, results):
        counts, resolved_edges = np.histogram(result.signed_distance.astype(float), bins=edges)
        centres = (resolved_edges[:-1] + resolved_edges[1:]) / 2.0
        axis.plot(centres, counts, label=result.label, color=color, linewidth=1.5)
    axis.set_xlim(float(args.histogram_x_min), float(args.histogram_x_max))
    axis.set_xlabel("Signed summit distance, comparison − main (bp)")
    axis.set_ylabel("Matched pairs")
    axis.set_title("Matched-position distance distributions")
    apply_distance_x_axis(
        axis,
        major_interval=args.distance_x_major_tick,
        minor_interval=args.distance_x_minor_tick,
    )
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
            axis.plot(
                x,
                [float(row["spearman_score_correlation"]) for row in subset],
                marker="o",
                color=color,
                label=(result.label if args.score_correlation == "spearman" else f"{result.label} Spearman"),
            )
        if args.score_correlation in {"pearson", "both"}:
            axis.plot(
                x,
                [float(row["pearson_score_correlation"]) for row in subset],
                marker="s",
                linestyle="--",
                color=color,
                label=(result.label if args.score_correlation == "pearson" else f"{result.label} Pearson"),
            )
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
    main_label: str,
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
            values = result.absolute_distance[result.group_mask(group)]
            pos = float(group_centres[group_index] + offsets[comp_index])
            positions[(group, result.label)] = pos
            if not values.size:
                continue
            box = axis.boxplot(
                [values], positions=[pos], widths=width, patch_artist=True,
                showfliers=bool(args.show_boxplot_outliers), whis=1.5,
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
    axis.set_xlabel(f"{main_label} peak score percentile group")
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
    axis.set_title(f"Matched distance by {main_label} score percentile", y=title_y)
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
    main_records = read_compact_positions(
        main_path, "main", args.main_summit_column, args.main_score_column,
        blacklist=blacklist, progress=reporter,
    )
    main_label = getattr(args, "_main_label", None) or _default_label(main_path)
    prefix = Path(args.output_prefix or f"{_default_label(main_path)}_compare_positions")
    prefix.parent.mkdir(parents=True, exist_ok=True)

    results: list[ComparisonArrays] = []
    summary_rows: list[dict[str, object]] = []
    histogram_rows: list[dict[str, object]] = []
    correlation_rows: list[dict[str, object]] = []
    score_distance_rows: list[dict[str, object]] = []
    percentile_path = Path(f"{prefix}_percentile_distances.tsv")
    outputs: dict[str, Path] = {}
    bounds = _parse_distance_bins(args.distance_bins)

    for index, spec in enumerate(specs, start=1):
        reporter.stage(f"Comparison {index}/{len(specs)}: {spec.label}")
        compare_records = read_compact_positions(
            spec.path, "comparison", args.compare_summit_column, args.compare_score_column,
            blacklist=blacklist, progress=reporter,
        )
        reporter.stage(
            f"One-to-one matching {spec.label}: main={main_records.count:,}; comparison={compare_records.count:,}"
        )
        match = match_compact_positions(
            main_records, compare_records, args.max_distance, progress=reporter,
            progress_stage=f"Matching {spec.label}",
        )
        if not match.pair_count:
            raise ValueError(f"No matched pairs were found for comparison {spec.label!r}.")
        arrays = _arrays_from_match(
            spec.label, spec.path, match, main_records.count, compare_records.count, args.percentile_interval
        )
        results.append(arrays)
        summary = _summary_row(arrays, main_path, main_label)
        summary["blacklist_overlapping_main_records_excluded"] = main_records.excluded_blacklist
        summary["blacklist_overlapping_compare_records_excluded"] = compare_records.excluded_blacklist
        summary_rows.append(summary)
        histogram_rows.extend(
            _histogram_rows(arrays, args.histogram_x_min, args.histogram_x_max, args.histogram_bin_width)
        )
        correlation_rows.extend(_distance_bin_rows(arrays, args.score_normalization, bounds))
        if args.write_detail_tables:
            _append_percentile_rows(percentile_path, arrays, main_label, write_header=(index == 1))
        sd_stats = _score_distance_stats(arrays, args.score_distance_type)
        sd_stats["main_label"] = main_label
        score_distance_rows.append(sd_stats)

        token = _safe_token(spec.label)
        if args.write_detail_tables and not args.skip_pairs_tsv:
            pair_path = Path(f"{prefix}_{token}_pairs.tsv")
            _write_pairs(pair_path, arrays, main_label, args)
            outputs[f"pairs_{token}"] = pair_path

        # Retain only the plotted sample plus full-data statistics so both
        # per-comparison figures remain reproducible without the huge pair TSV.
        score_agreement_source = Path(f"{prefix}_{token}_score_agreement.tsv.gz")
        score_distance_source = Path(f"{prefix}_{token}_main_score_vs_distance.tsv.gz")
        _write_score_agreement_plot_source(score_agreement_source, arrays, main_label, args)
        _write_score_distance_plot_source(score_distance_source, arrays, main_label, args, sd_stats)
        outputs[f"score_agreement_source_{token}"] = score_agreement_source
        outputs[f"score_distance_source_{token}"] = score_distance_source

        score_agreement_plot = _plot_score_agreement(prefix, arrays, main_label, args)
        score_distance_plot = _plot_score_distance(prefix, arrays, main_label, args, sd_stats)
        _attach_plot_source(score_agreement_plot, score_agreement_source, "compare-positions-score")
        _attach_plot_source(score_distance_plot, score_distance_source, "compare-positions-score-distance")
        outputs[f"score_agreement_plot_{token}"] = score_agreement_plot
        outputs[f"score_distance_plot_{token}"] = score_distance_plot

    summary_path = Path(f"{prefix}_summary.tsv")
    histogram_path = Path(f"{prefix}_distance_histogram.tsv")
    correlation_path = Path(f"{prefix}_correlation_by_distance.tsv")
    percentile_summary_path = Path(f"{prefix}_percentile_summary.tsv")
    score_distance_path = Path(f"{prefix}_main_score_vs_distance_statistics.tsv")
    _write_tsv(summary_path, summary_rows)
    _write_tsv(histogram_path, histogram_rows)
    _write_tsv(correlation_path, correlation_rows)
    _write_tsv(percentile_summary_path, _percentile_summary_rows(results, args.percentile_interval))
    _write_tsv(score_distance_path, score_distance_rows)
    outputs.update({
        "summary": summary_path,
        "distance_histogram": histogram_path,
        "correlation_by_distance": correlation_path,
        "percentile_summary": percentile_summary_path,
        "score_distance_statistics": score_distance_path,
    })
    if args.write_detail_tables:
        outputs["percentile_distances"] = percentile_path

    statistics_rows: list[dict[str, object]] = []
    if args.stats and len(results) >= 2:
        statistics_rows = _pairwise_statistics(
            results, args.percentile_interval, args.stats_test, args.p_adjust
        )
        statistics_path = Path(f"{prefix}_percentile_statistics.tsv")
        _write_tsv(statistics_path, statistics_rows)
        outputs["percentile_statistics"] = statistics_path

    # Always retain the compact plot source needed to faithfully recreate the
    # percentile boxplot.  It contains box/whisker statistics, only the actual
    # outliers, and optional statistical annotations rather than every pair.
    percentile_boxplot_source = Path(f"{prefix}_percentile_boxplot.tsv")
    _write_tsv(
        percentile_boxplot_source,
        _percentile_boxplot_source_rows(
            results, main_label, args.percentile_interval, statistics_rows, args
        ),
    )
    outputs["percentile_boxplot_source"] = percentile_boxplot_source

    distance_histogram_plot = _plot_distance_histograms(prefix, results, args)
    correlation_by_distance_plot = _plot_correlation_by_distance(
        prefix, correlation_rows, results, args
    )
    percentile_boxplot = _plot_percentile_boxplot(
        prefix, results, main_label, args.percentile_interval, statistics_rows, args
    )
    _attach_plot_source(distance_histogram_plot, histogram_path, "compare-positions-histogram")
    _attach_plot_source(correlation_by_distance_plot, correlation_path, "compare-positions-correlation")
    _attach_plot_source(percentile_boxplot, percentile_boxplot_source, "compare-positions-percentile-boxplot")
    outputs["distance_histogram_plot"] = distance_histogram_plot
    outputs["correlation_by_distance_plot"] = correlation_by_distance_plot
    outputs["percentile_boxplot"] = percentile_boxplot

    if not args.quiet:
        reporter.emit(f"Main positions: {main_records.count:,}")
        for result in results:
            reporter.emit(
                f"{result.label}: {result.compare_count:,} positions; {result.pair_count:,} one-to-one pairs"
            )
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
