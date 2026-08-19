#!/usr/bin/env python3
"""Measure peak-to-peak distances for one or more neighbour orders.

This module implements the ``nucleosuite distances`` command. It combines the
adjacent-peak distance workflow with higher-order peak comparisons (+1 through
+N), optional chromatin-state stratification, score filtering, percentile
sweeps, and filtered BED output.

Peak input
----------
The peak file must be BED with at least five columns::

    chromosome  start  end  name  score  [strand ...]

By default, the peak position is the integer midpoint of ``start`` and ``end``.
A one-based ``--position-column`` can instead supply an explicit genomic
position. Scores are read from the one-based ``--score-column`` (default: 5).

Chromatin-state input
---------------------
An optional BED state file can be provided with ``--state-bed``. The state label is read from
``--state-label-column`` (default: 4). A peak is assigned according to the
state interval(s) overlapping its point position.

Distance orders
---------------
For position-sorted retained peaks, order +1 compares adjacent peaks, +2
compares peak i with peak i+2, and so on through ``--max-order``. Pooled
"All" distributions include every valid pair. State-specific distributions
include pairs whose two endpoint peaks have the same state label.

Output
------
Each threshold produces:

* ``<prefix>_scorepctX_metadata.tsv``: analysis parameters.
* ``<prefix>_scorepctX_distances.tsv``: tidy distance distributions.
* ``<prefix>_scorepctX_summary.tsv``: mean, median, raw mode, and smoothed mode.
* ``<prefix>_scorepctX_nrl_regression_<scope>.tsv``: order-specific peak distances, fitted values, and residuals when multiple orders are available.
* ``<prefix>_scorepctX_nrl_regression_<scope>.png``: peak-distance regression annotated with NRL and R-squared.
* ``<prefix>_scorepctX_nrl_regression_summary.tsv``: regression statistics across requested scopes.
* ``<prefix>_scorepctX_duplicates.tsv``: duplicate retained peak positions,
  when duplicates are present.
* ``<prefix>_scorepctX_filtered.bed``: optional retained peak records.
* ``<prefix>_percentile_sweep_*_count.png``: raw, unsmoothed count overlays for
  ``--pct-range`` or ``--pct-values``.
* ``<prefix>_percentile_sweep_*_percentage.png``: independently normalized
  raw percentage overlays for ``--pct-range`` or ``--pct-values``.
* ``<prefix>_percentile_sweep_*_peak_counts.png``: standalone retained-peak
  count tables for ``--pct-range`` or ``--pct-values``.

"""

from __future__ import annotations

import argparse
import bisect
import gzip
import math
import re
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, NamedTuple, Sequence, TextIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

from nucleosuite.core.regions import canonical_contig_key, resolve_contig_name
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.dac import categorize_state_name, parse_state_category_rules
from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.progress import ProgressReporter

from nucleosuite.io import open_text as open_interval_text
from nucleosuite.io.intervals import INTERVAL_FORMATS, finalise_interval_files, normalise_chrom_sizes
from nucleosuite.percentiles import (
    equal_rank_bins,
    randomized_score_order,
    rank_bins_from_boundaries,
)
HEADER_PREFIXES = ("#", "track", "browser")

class PeakRecord(NamedTuple):
    """Parsed BED peak retained for distance analysis and BED8 export."""

    position: int
    score: float
    state: str
    start: int
    end: int
    name: str
    strand: str


@dataclass
class RankMembership:
    """Peak records partitioned once into exact equal-rank groups."""

    groups: tuple[dict[str, list[PeakRecord]], ...]


@dataclass(frozen=True)
class StateInterval:
    """One half-open labelled interval from a chromatin-state BED file."""

    start: int
    end: int
    label: str
    rgb: str = "128,128,128"


@dataclass
class IntervalIndex:
    """Bisectable interval index for one chromosome or contig."""

    intervals: list[StateInterval]
    starts: list[int]
    prefix_max_ends: list[int]


@dataclass
class ParseSummary:
    """Parsing counts for one BED input."""

    data_lines: int = 0
    used_lines: int = 0
    skipped_lines: int = 0
    blacklisted_lines: int = 0


@dataclass
class ThresholdSelection:
    """One score threshold to analyse."""

    requested_percentile: float
    effective_percentile: float
    threshold: float
    target_peaks: int | None = None
    percentile_lower: float | None = None
    percentile_upper: float | None = None
    score_upper_bound: float | None = None
    rank_start: int | None = None
    rank_stop: int | None = None
    rank_seed: int | None = None
    rank_score_max: float | None = None
    rank_bin_index: int | None = None
    rank_membership: RankMembership | None = None

    @property
    def is_bin(self) -> bool:
        return self.percentile_lower is not None and self.percentile_upper is not None

    @property
    def is_rank_bin(self) -> bool:
        return self.rank_start is not None and self.rank_stop is not None

    @property
    def label(self) -> str:
        if self.is_bin:
            return (
                f"{percentile_label(self.percentile_lower)}-"
                f"{percentile_label(self.percentile_upper)}"
            )
        return percentile_label(self.effective_percentile)


@dataclass
class DistanceResults:
    """Distance counters and retained-peak metadata for one threshold."""

    chrom_state: dict[int, dict[str, dict[str, Counter[int]]]]
    chrom_all: dict[int, dict[str, Counter[int]]]
    genome_state: dict[int, dict[str, Counter[int]]]
    genome_all: dict[int, Counter[int]]
    duplicates: dict[str, dict[int, list[PeakRecord]]]
    retained_by_chrom: dict[str, list[PeakRecord]]
    threshold_pass_count: int
    retained_count: int


@dataclass
class DistributionStatistics:
    """Summary statistics and smoothed values for one distance distribution."""

    total_pairs: int
    raw_mode: int
    smoothed_mode: int
    median: float
    mean: float
    smoothed_median: float
    smoothed_mean: float
    observed_min: int
    observed_max: int
    distances: np.ndarray
    raw_counts: np.ndarray
    smoothed_counts: np.ndarray
    raw_percent: np.ndarray
    smoothed_percent: np.ndarray


@dataclass(frozen=True)
class OrderPeak:
    """Highest-count distance selected for one neighbour order."""

    order: int
    peak_distance: int
    peak_count: int
    total_pairs: int


@dataclass(frozen=True)
class NRLRegression:
    """Linear regression of order-specific peak distance against neighbour order."""

    scope: str
    chromosome: str
    peaks: tuple[OrderPeak, ...]
    slope: float
    intercept: float
    r_squared: float


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    """Open BED text, BED.gz, or bigBed input."""
    if mode not in {"rt", "r"}:
        raise ValueError("NucleoSuite interval inputs are read-only")
    return open_interval_text(path)


def strip_known_suffix(path: str | Path) -> str:
    """Return a basename with common BED/text compression suffixes removed."""
    name = Path(path).name
    lower = name.lower()
    for suffix in (".bed.gz", ".bigbed", ".bb", ".bedgraph.gz", ".tsv.gz", ".txt.gz", ".gz"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def natural_sort_key(value: str) -> tuple[object, ...]:
    """Sort mixed text/numeric contig names naturally without genome assumptions."""
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", value)
    )


def split_fields(line: str) -> list[str]:
    """Split a BED record, preferring tabs while accepting whitespace."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) == 1:
        fields = line.split()
    return fields


def validate_one_based_column(value: int, option_name: str, *, allow_zero: bool = False) -> int:
    """Validate and convert a user-facing one-based column to a zero-based index."""
    minimum = 0 if allow_zero else 1
    if value < minimum:
        relation = ">= 0" if allow_zero else ">= 1"
        raise ValueError(f"{option_name} must be {relation}")
    return value - 1 if value > 0 else -1


def build_state_indexes(
    state_path: str | Path,
    *,
    label_column: int = 4,
    color_column: int = 9,
    strict: bool = False,
    progress: ProgressReporter | None = None,
) -> tuple[dict[str, IntervalIndex], ParseSummary]:
    """Read a chromatin-state BED/BED.gz and build one interval index per contig."""
    label_index = validate_one_based_column(label_column, "--state-label-column")
    color_index = validate_one_based_column(color_column, "--state-color-column")
    by_chrom: dict[str, list[StateInterval]] = defaultdict(list)
    summary = ParseSummary()
    seen_contigs: set[str] = set()

    with open_text(state_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(HEADER_PREFIXES):
                continue

            summary.data_lines += 1
            fields = split_fields(raw_line)
            error: str | None = None

            if len(fields) <= max(2, label_index):
                error = f"expected state label column {label_column}"
            else:
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    error = "state start/end must be integers"
                else:
                    if start < 0 or end <= start:
                        error = f"invalid half-open interval {start}-{end}"

            if error is not None:
                summary.skipped_lines += 1
                if strict:
                    raise ValueError(f"{state_path}:{line_number}: {error}")
                continue

            label = fields[label_index]
            rgb = fields[color_index] if color_index < len(fields) else "128,128,128"
            chrom = fields[0]
            if progress is not None and chrom not in seen_contigs:
                seen_contigs.add(chrom)
                progress.reading_contig("chromatin states", chrom)
            by_chrom[chrom].append(StateInterval(start, end, label, rgb))
            summary.used_lines += 1

    indexes: dict[str, IntervalIndex] = {}
    for chrom, intervals in by_chrom.items():
        intervals.sort(key=lambda interval: (interval.start, interval.end, interval.label))
        starts = [interval.start for interval in intervals]
        prefix_max_ends: list[int] = []
        running_max = -1
        for interval in intervals:
            running_max = max(running_max, interval.end)
            prefix_max_ends.append(running_max)
        indexes[chrom] = IntervalIndex(intervals, starts, prefix_max_ends)

    if not indexes:
        raise ValueError(f"No valid chromatin-state intervals were found in: {state_path}")
    return indexes, summary


def labels_at_position(index: IntervalIndex, position: int) -> list[str]:
    """Return all interval labels overlapping one 0-based point position."""
    right = bisect.bisect_right(index.starts, position) - 1
    if right < 0:
        return []

    labels: list[str] = []
    i = right
    while i >= 0 and index.prefix_max_ends[i] > position:
        interval = index.intervals[i]
        if interval.start <= position < interval.end:
            labels.append(interval.label)
        i -= 1

    # Return unique labels in genomic interval order.
    labels.reverse()
    return list(dict.fromkeys(labels))


def assign_state_label(
    state_indexes: Mapping[str, IntervalIndex] | None,
    chrom: str,
    position: int,
    *,
    category_rules,
) -> str:
    """Assign a point to a state label and apply optional category rules."""
    if state_indexes is None:
        return "All"

    index = state_indexes.get(chrom)
    if index is None:
        try:
            resolved = resolve_contig_name(
                chrom,
                state_indexes.keys(),
                source_label="chromatin-state BED",
            )
        except (KeyError, ValueError):
            resolved = None
        if resolved is not None:
            index = state_indexes.get(resolved)
    labels = labels_at_position(index, position) if index is not None else []
    if not labels:
        return "NA"

    categorized = [categorize_state_name(label, category_rules) for label in labels]
    unique = sorted(dict.fromkeys(categorized))
    return unique[0] if len(unique) == 1 else "|".join(unique)


def load_peaks(
    peak_path: str | Path,
    *,
    state_indexes: Mapping[str, IntervalIndex] | None,
    category_rules=(),
    position_column: int | None = None,
    score_column: int = 5,
    strict: bool = False,
    blacklist: BlacklistIndex | None = None,
    progress: ProgressReporter | None = None,
) -> tuple[dict[str, list[PeakRecord]], np.ndarray, ParseSummary]:
    """Load, score, annotate, and position-sort peaks from a BED file."""
    score_index = validate_one_based_column(score_column, "--score-column")
    position_index = (
        validate_one_based_column(position_column, "--position-column")
        if position_column is not None
        else None
    )

    peaks_by_chrom: dict[str, list[PeakRecord]] = defaultdict(list)
    summary = ParseSummary()
    seen_contigs: set[str] = set()

    with open_text(peak_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.strip()
            if not stripped or stripped.startswith(HEADER_PREFIXES):
                continue

            summary.data_lines += 1
            fields = split_fields(raw_line)
            required_index = max(4, score_index, position_index or 0)
            error: str | None = None

            if len(fields) <= required_index:
                error = f"expected at least {required_index + 1} BED columns"
            else:
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                    if start < 0 or end <= start:
                        raise ValueError
                except ValueError:
                    error = "BED start/end must define a valid non-negative interval"

            if error is None:
                try:
                    score = float(fields[score_index])
                except ValueError:
                    error = f"score column {score_column} is not numeric"
                else:
                    if not math.isfinite(score):
                        error = f"score column {score_column} is not finite"

            if error is None:
                if position_index is None:
                    position = (start + end) // 2
                else:
                    try:
                        position = int(fields[position_index])
                        if position < start or position >= end:
                            raise ValueError
                    except ValueError:
                        error = (
                            f"position column {position_column} must contain a genomic "
                            "position within the BED interval"
                        )

            if error is not None:
                summary.skipped_lines += 1
                if strict:
                    raise ValueError(f"{peak_path}:{line_number}: {error}")
                continue

            chrom = fields[0]
            if progress is not None and chrom not in seen_contigs:
                seen_contigs.add(chrom)
                progress.reading_contig("peaks", chrom)
            if blacklist is not None and blacklist.overlaps(chrom, start, end):
                summary.skipped_lines += 1
                summary.blacklisted_lines += 1
                continue
            name = fields[3] or f"{chrom}:{position}"
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-", "."} else "."
            state = assign_state_label(
                state_indexes,
                chrom,
                position,
                category_rules=category_rules,
            )
            peaks_by_chrom[chrom].append(
                PeakRecord(position, score, state, start, end, name, strand)
            )
            summary.used_lines += 1

    for chrom in peaks_by_chrom:
        peaks_by_chrom[chrom].sort(key=lambda record: record.position)

    if summary.used_lines == 0:
        raise ValueError(f"No valid scored BED peaks were found in: {peak_path}")

    ordered_peaks = dict(peaks_by_chrom)
    scores = np.fromiter(
        (
            record.score
            for peaks in ordered_peaks.values()
            for record in peaks
        ),
        dtype=float,
        count=summary.used_lines,
    )
    return ordered_peaks, scores, summary

def percentile_values(start: float, stop: float, step: float) -> list[float]:
    """Return an inclusive, floating-point-safe percentile sequence."""
    if step <= 0:
        raise ValueError("--pct-step must be greater than 0")

    lower = min(100.0, max(0.0, start))
    upper = min(100.0, max(0.0, stop))
    if upper < lower:
        lower, upper = upper, lower

    values: list[float] = []
    current = lower
    tolerance = max(1e-12, abs(step) * 1e-10)
    while current <= upper + tolerance:
        values.append(min(upper, current))
        current += step

    if not values or abs(values[-1] - upper) > tolerance:
        values.append(upper)

    deduplicated: list[float] = []
    for value in values:
        if not deduplicated or abs(value - deduplicated[-1]) > tolerance:
            deduplicated.append(value)
    return deduplicated


def explicit_percentile_values(tokens: Sequence[str] | None) -> list[float] | None:
    """Parse ordered, unique percentile values from comma/space-separated tokens."""
    if not tokens:
        return None
    values: list[float] = []
    seen: set[float] = set()
    for token in tokens:
        for part in token.split(","):
            text = part.strip()
            if not text:
                raise ValueError("--pct-values contains an empty percentile")
            try:
                value = float(text)
            except ValueError as exc:
                raise ValueError(
                    f"--pct-values entry is not numeric: {text!r}"
                ) from exc
            if not math.isfinite(value) or not 0.0 <= value <= 100.0:
                raise ValueError("--pct-values entries must be between 0 and 100")
            canonical = round(value, 12)
            if canonical not in seen:
                seen.add(canonical)
                values.append(value)
    if not values:
        raise ValueError("--pct-values requires at least one percentile")
    return values


def percentile_bin_bounds(percentiles: Sequence[float]) -> list[tuple[float, float]]:
    """Convert ordered percentile boundaries into contiguous bins."""
    boundaries = [float(value) for value in percentiles]
    if len(boundaries) < 2:
        raise ValueError("--pct-bins requires at least two percentile boundaries")

    bounds: list[tuple[float, float]] = []
    lower = boundaries[0]
    for upper in boundaries[1:]:
        if upper <= lower:
            raise ValueError(
                "--pct-bins percentile boundaries must be in strictly increasing order"
            )
        bounds.append((lower, upper))
        lower = upper
    return bounds


def choose_thresholds(
    scores: np.ndarray,
    *,
    score_percentile: float,
    target_peaks: int | None,
    pct_range: bool,
    pct_lower: float,
    pct_upper: float,
    pct_step: float,
    pct_values: Sequence[str] | None = None,
    pct_bins: bool = False,
    pct_bin_size: float | None = None,
    pct_bin_seed: int = 1,
    bin_tie_mode: str = "split",
) -> list[ThresholdSelection]:
    """Create one or more percentile-derived score thresholds or score bins."""
    if scores.size == 0:
        raise ValueError("Cannot choose thresholds from an empty score array")
    if bin_tie_mode not in {"split", "keep"}:
        raise ValueError("--bin-tie-mode must be 'split' or 'keep'")

    if target_peaks is not None:
        if pct_bins or pct_bin_size is not None:
            raise ValueError(
                "--target-peaks cannot be combined with --pct-bins or --pct-bin-size"
            )
        if target_peaks < 1:
            raise ValueError("--target-peaks must be at least 1")
        retained_target = min(int(target_peaks), int(scores.size))
        effective = 100.0 * (1.0 - retained_target / float(scores.size))
        threshold = float(np.percentile(scores, effective))
        return [
            ThresholdSelection(
                requested_percentile=score_percentile,
                effective_percentile=effective,
                threshold=threshold,
                target_peaks=target_peaks,
            )
        ]

    fixed_rank_bins = None
    if pct_bin_size is not None:
        if pct_bins or pct_range or pct_values:
            raise ValueError(
                "--pct-bin-size cannot be combined with --pct-bins, --pct-range, "
                "or --pct-values"
            )
        if not math.isclose(score_percentile, 0.0, abs_tol=1e-12):
            raise ValueError("--pct-bin-size cannot be combined with --score-percentile")
        if pct_bin_seed < 0:
            raise ValueError("--pct-bin-seed must be non-negative")
        if not math.isfinite(pct_bin_size) or pct_bin_size <= 0.0 or pct_bin_size > 100.0:
            raise ValueError("--pct-bin-size must be greater than 0 and at most 100")
        bin_count_float = 100.0 / float(pct_bin_size)
        bin_count = int(round(bin_count_float))
        if not math.isclose(
            bin_count * float(pct_bin_size),
            100.0,
            rel_tol=0.0,
            abs_tol=1e-9,
        ):
            raise ValueError("--pct-bin-size must divide 100 into equal percentile ranges")
        boundaries = [round(index * float(pct_bin_size), 12) for index in range(bin_count)]
        boundaries.append(100.0)
        bounds = list(zip(boundaries[:-1], boundaries[1:]))
        if bin_tie_mode == "split":
            fixed_rank_bins = equal_rank_bins(int(scores.size), pct_bin_size)
    else:
        explicit = explicit_percentile_values(pct_values)
        if explicit is not None:
            requested = explicit
        elif pct_range:
            requested = percentile_values(pct_lower, pct_upper, pct_step)
        else:
            if pct_bins:
                raise ValueError("--pct-bins requires --pct-values or --pct-range")
            if not 0.0 <= score_percentile <= 100.0:
                raise ValueError("--score-percentile must be between 0 and 100")
            requested = [float(score_percentile)]

        if not pct_bins:
            thresholds = np.percentile(scores, requested)
            return [
                ThresholdSelection(
                    requested_percentile=float(percentile),
                    effective_percentile=float(percentile),
                    threshold=float(threshold),
                )
                for percentile, threshold in zip(requested, np.atleast_1d(thresholds))
            ]
        bounds = percentile_bin_bounds(requested)

    if bin_tie_mode == "split":
        if pct_bin_seed < 0:
            raise ValueError("--pct-bin-seed must be non-negative")
        rank_bins = fixed_rank_bins or rank_bins_from_boundaries(
            int(scores.size),
            tuple([bounds[0][0], *(upper for _lower, upper in bounds)]),
        )
        ordered_scores = np.sort(scores)
        return [
            ThresholdSelection(
                requested_percentile=rank_bin.percentile_upper,
                effective_percentile=rank_bin.percentile_upper,
                threshold=float(ordered_scores[rank_bin.rank_start]),
                percentile_lower=rank_bin.percentile_lower,
                percentile_upper=rank_bin.percentile_upper,
                score_upper_bound=None,
                rank_start=rank_bin.rank_start,
                rank_stop=rank_bin.rank_stop,
                rank_seed=pct_bin_seed,
                rank_score_max=float(ordered_scores[rank_bin.rank_stop - 1]),
            )
            for rank_bin in rank_bins
        ]

    boundary_values = sorted({value for pair in bounds for value in pair})
    boundary_scores = {
        percentile: float(threshold)
        for percentile, threshold in zip(
            boundary_values,
            np.atleast_1d(np.percentile(scores, boundary_values)),
        )
    }
    return [
        ThresholdSelection(
            requested_percentile=upper,
            effective_percentile=upper,
            threshold=boundary_scores[lower],
            percentile_lower=lower,
            percentile_upper=upper,
            score_upper_bound=(
                None if math.isclose(upper, 100.0) else boundary_scores[upper]
            ),
        )
        for lower, upper in bounds
    ]


def resolve_duplicate_group(
    group: list[PeakRecord],
    *,
    policy: str,
    chrom: str,
) -> list[PeakRecord]:
    """Resolve records sharing one retained peak position."""
    if len(group) == 1:
        return group

    position = group[0].position
    if policy == "keep":
        return group
    if policy == "first":
        return [group[0]]
    if policy == "highest-score":
        return [max(group, key=lambda record: record.score)]
    if policy == "error":
        raise ValueError(
            f"Duplicate retained peak position encountered at {chrom}:{position}"
        )
    raise ValueError(f"Unknown duplicate policy: {policy}")


def attach_rank_bin_records(
    peaks_by_chrom: Mapping[str, Sequence[PeakRecord]],
    selections: Sequence[ThresholdSelection],
    *,
    scores: np.ndarray | None = None,
) -> RankMembership | None:
    """Partition records once and attach shared equal-rank membership."""
    if not selections or not all(selection.is_rank_bin for selection in selections):
        return None
    seed = selections[0].rank_seed
    assert seed is not None
    flat_scores = (
        np.asarray(scores, dtype=float)
        if scores is not None
        else np.fromiter(
            (
                record.score
                for peaks in peaks_by_chrom.values()
                for record in peaks
            ),
            dtype=float,
        )
    )
    expected_peak_count = sum(len(peaks) for peaks in peaks_by_chrom.values())
    if flat_scores.size != expected_peak_count:
        raise ValueError("Rank-bin score array does not match the loaded peak records")
    order = randomized_score_order(flat_scores, seed)
    bin_by_flat_record = np.full(flat_scores.size, -1, dtype=np.int32)
    for bin_index, selection in enumerate(selections):
        assert selection.rank_start is not None and selection.rank_stop is not None
        bin_by_flat_record[order[selection.rank_start : selection.rank_stop]] = bin_index

    groups: list[dict[str, list[PeakRecord]]] = [
        defaultdict(list) for _selection in selections
    ]
    offset = 0
    for chrom, peaks in peaks_by_chrom.items():
        assignments = bin_by_flat_record[offset : offset + len(peaks)]
        for record, bin_index in zip(peaks, assignments):
            if int(bin_index) >= 0:
                groups[int(bin_index)][chrom].append(record)
        offset += len(peaks)

    membership = RankMembership(tuple(dict(group) for group in groups))
    for bin_index, selection in enumerate(selections):
        selection.rank_bin_index = bin_index
        selection.rank_membership = membership
    return membership


def retain_peaks_for_threshold(
    peaks_by_chrom: Mapping[str, Sequence[PeakRecord]],
    *,
    threshold: float,
    score_upper_bound: float | None = None,
    rank_membership: RankMembership | None = None,
    rank_bin_index: int | None = None,
    duplicate_policy: str,
) -> tuple[
    dict[str, list[PeakRecord]],
    dict[str, dict[int, list[PeakRecord]]],
    int,
    int,
]:
    """Apply lower-inclusive/upper-exclusive score bounds and duplicate policy."""
    retained_by_chrom: dict[str, list[PeakRecord]] = {}
    duplicates: dict[str, dict[int, list[PeakRecord]]] = defaultdict(dict)
    threshold_pass_count = 0
    retained_count = 0

    source_peaks = (
        rank_membership.groups[rank_bin_index]
        if rank_membership is not None and rank_bin_index is not None
        else peaks_by_chrom
    )
    for chrom, peaks in source_peaks.items():
        if rank_membership is not None:
            passed = list(peaks)
        else:
            passed = [
                record
                for record in peaks
                if record.score >= threshold
                and (score_upper_bound is None or record.score < score_upper_bound)
            ]
        threshold_pass_count += len(passed)
        if not passed:
            continue

        selected: list[PeakRecord] = []
        group: list[PeakRecord] = [passed[0]]

        for record in passed[1:]:
            if record.position == group[0].position:
                group.append(record)
                continue

            if len(group) > 1:
                duplicates[chrom][group[0].position] = list(group)
            selected.extend(resolve_duplicate_group(group, policy=duplicate_policy, chrom=chrom))
            group = [record]

        if len(group) > 1:
            duplicates[chrom][group[0].position] = list(group)
        selected.extend(resolve_duplicate_group(group, policy=duplicate_policy, chrom=chrom))

        if selected:
            retained_by_chrom[chrom] = selected
            retained_count += len(selected)

    return retained_by_chrom, dict(duplicates), threshold_pass_count, retained_count


def compute_distance_counts(
    peaks_by_chrom: Mapping[str, Sequence[PeakRecord]],
    *,
    threshold: float,
    score_upper_bound: float | None = None,
    rank_membership: RankMembership | None = None,
    rank_bin_index: int | None = None,
    min_distance: int,
    max_distance: int,
    max_order: int,
    duplicate_policy: str,
) -> DistanceResults:
    """Count peak distances for every order from +1 through ``max_order``."""
    if min_distance < 0:
        raise ValueError("--min-distance must be at least 0")
    if max_distance < min_distance:
        raise ValueError("--max-distance must be greater than or equal to --min-distance")
    if max_order < 1:
        raise ValueError("--max-order must be at least 1")

    retained_by_chrom, duplicates, pass_count, retained_count = retain_peaks_for_threshold(
        peaks_by_chrom,
        threshold=threshold,
        score_upper_bound=score_upper_bound,
        rank_membership=rank_membership,
        rank_bin_index=rank_bin_index,
        duplicate_policy=duplicate_policy,
    )

    chrom_state: dict[int, dict[str, dict[str, Counter[int]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(Counter))
    )
    chrom_all: dict[int, dict[str, Counter[int]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    genome_state: dict[int, dict[str, Counter[int]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    genome_all: dict[int, Counter[int]] = defaultdict(Counter)

    for chrom, peaks in retained_by_chrom.items():
        n_peaks = len(peaks)
        for i, first in enumerate(peaks):
            position_1 = first.position
            state_1 = first.state
            highest_order = min(max_order, n_peaks - i - 1)

            for order in range(1, highest_order + 1):
                second = peaks[i + order]
                position_2 = second.position
                state_2 = second.state
                distance = position_2 - position_1

                # Sorted positions guarantee later orders cannot return below max_distance.
                if distance > max_distance:
                    break
                if distance < min_distance or distance == 0:
                    continue

                # Pooled distributions contain every pair, including state transitions.
                chrom_all[order][chrom][distance] += 1
                genome_all[order][distance] += 1

                # State-specific distributions follow the source scripts: endpoint
                # states must match; intervening peaks may have other state labels.
                if state_1 == state_2:
                    chrom_state[order][chrom][state_1][distance] += 1
                    genome_state[order][state_1][distance] += 1

    return DistanceResults(
        chrom_state={
            order: {
                chrom: dict(state_counts)
                for chrom, state_counts in chrom_counts.items()
            }
            for order, chrom_counts in chrom_state.items()
        },
        chrom_all={order: dict(chrom_counts) for order, chrom_counts in chrom_all.items()},
        genome_state={order: dict(state_counts) for order, state_counts in genome_state.items()},
        genome_all=dict(genome_all),
        duplicates=duplicates,
        retained_by_chrom=retained_by_chrom,
        threshold_pass_count=pass_count,
        retained_count=retained_count,
    )


def weighted_mean(counter: Mapping[int, int]) -> float:
    """Return a count-weighted mean distance."""
    total = sum(counter.values())
    if total <= 0:
        raise ValueError("Cannot calculate a weighted mean from zero counts")
    return sum(distance * count for distance, count in counter.items()) / total


def weighted_median(counter: Mapping[int, int]) -> float:
    """Return the conventional midpoint weighted median distance."""
    total = sum(counter.values())
    if total <= 0:
        raise ValueError("Cannot calculate a weighted median from zero counts")

    lower_rank = (total - 1) // 2
    upper_rank = total // 2
    cumulative = 0
    lower_value: int | None = None
    upper_value: int | None = None

    for distance in sorted(counter):
        cumulative += counter[distance]
        if lower_value is None and cumulative > lower_rank:
            lower_value = distance
        if cumulative > upper_rank:
            upper_value = distance
            break

    assert lower_value is not None and upper_value is not None
    return (lower_value + upper_value) / 2.0


def adaptive_savgol(
    values: np.ndarray,
    *,
    window_length: int,
    polyorder: int,
) -> np.ndarray:
    """Smooth values when the requested Savitzky-Golay settings are valid."""
    if values.size == 0 or window_length <= 0:
        return values.astype(float, copy=True)
    if window_length % 2 == 0:
        raise ValueError("Savitzky-Golay window lengths must be odd")
    if polyorder < 0:
        raise ValueError("Savitzky-Golay polynomial orders must be non-negative")

    # Keep short distributions unsmoothed rather than silently changing the
    # requested window, matching the conservative behaviour of the source tools.
    if values.size < window_length or window_length <= polyorder or window_length < 3:
        return values.astype(float, copy=True)

    return savgol_filter(values.astype(float), window_length=window_length, polyorder=polyorder)


def summarize_distribution(
    counter: Mapping[int, int],
    *,
    count_smooth_window: int,
    count_smooth_polyorder: int,
    percent_smooth_window: int,
    percent_smooth_polyorder: int,
) -> DistributionStatistics:
    """Calculate summary statistics and gap-aware smoothed distance curves."""
    positive = {distance: int(count) for distance, count in counter.items() if count > 0}
    if not positive:
        raise ValueError("Cannot summarize an empty distance distribution")

    observed_min = min(positive)
    observed_max = max(positive)
    distances = np.arange(observed_min, observed_max + 1, dtype=int)
    raw_counts = np.asarray([positive.get(int(distance), 0) for distance in distances], dtype=float)
    total_pairs = int(raw_counts.sum())
    raw_percent = raw_counts * (100.0 / total_pairs)

    smoothed_counts = adaptive_savgol(
        raw_counts,
        window_length=count_smooth_window,
        polyorder=count_smooth_polyorder,
    )
    smoothed_counts = np.clip(smoothed_counts, 0.0, None)

    smoothed_percent = adaptive_savgol(
        raw_percent,
        window_length=percent_smooth_window,
        polyorder=percent_smooth_polyorder,
    )
    smoothed_percent = np.clip(smoothed_percent, 0.0, None)

    raw_mode = int(distances[int(np.argmax(raw_counts))])
    smoothed_mode = int(distances[int(np.argmax(smoothed_counts))])

    smoothed_total = float(smoothed_counts.sum())
    if smoothed_total > 0:
        smoothed_mean = float(np.dot(distances, smoothed_counts) / smoothed_total)
        cumulative = np.cumsum(smoothed_counts)
        smoothed_median = float(distances[int(np.searchsorted(cumulative, smoothed_total / 2.0, side="left"))])
    else:
        smoothed_mean = float("nan")
        smoothed_median = float("nan")

    percent_total = float(smoothed_percent.sum())
    if percent_total > 0:
        smoothed_percent = smoothed_percent * (100.0 / percent_total)

    return DistributionStatistics(
        total_pairs=total_pairs,
        raw_mode=raw_mode,
        smoothed_mode=smoothed_mode,
        median=weighted_median(positive),
        mean=weighted_mean(positive),
        smoothed_median=smoothed_median,
        smoothed_mean=smoothed_mean,
        observed_min=observed_min,
        observed_max=observed_max,
        distances=distances,
        raw_counts=raw_counts,
        smoothed_counts=smoothed_counts,
        raw_percent=raw_percent,
        smoothed_percent=smoothed_percent,
    )



def highest_count_distance(counter: Mapping[int, int]) -> OrderPeak:
    """Return the raw mode and its count, preferring the smaller distance on ties."""
    positive = [(int(distance), int(count)) for distance, count in counter.items() if int(count) > 0]
    if not positive:
        raise ValueError("Cannot select a peak distance from an empty distribution")
    peak_distance, peak_count = min(positive, key=lambda item: (-item[1], item[0]))
    return OrderPeak(
        order=0,
        peak_distance=peak_distance,
        peak_count=peak_count,
        total_pairs=sum(count for _, count in positive),
    )


def regress_order_peaks(
    peaks: Sequence[OrderPeak],
    *,
    scope: str,
    chromosome: str,
) -> NRLRegression:
    """Fit peak distance = intercept + NRL × order by ordinary least squares."""
    ordered = tuple(sorted(peaks, key=lambda peak: peak.order))
    if len(ordered) < 2:
        raise ValueError("At least two populated neighbour orders are required for NRL regression")

    orders = np.asarray([peak.order for peak in ordered], dtype=np.float64)
    distances = np.asarray([peak.peak_distance for peak in ordered], dtype=np.float64)
    slope, intercept = np.polyfit(orders, distances, 1)
    fitted = intercept + slope * orders
    residual_sum_squares = float(np.sum((distances - fitted) ** 2))
    total_sum_squares = float(np.sum((distances - np.mean(distances)) ** 2))
    r_squared = (
        1.0 - residual_sum_squares / total_sum_squares
        if total_sum_squares > 0.0
        else float("nan")
    )
    return NRLRegression(
        scope=scope,
        chromosome=chromosome,
        peaks=ordered,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r_squared),
    )


def collect_nrl_regressions(
    results: DistanceResults,
    *,
    max_order: int,
    include_chromosomes: bool,
    include_genome: bool,
    nrl_mode: str = "smoothed",
    count_smooth_window: int = 21,
    count_smooth_polyorder: int = 2,
) -> list[NRLRegression]:
    """Build pooled order-mode regressions from raw or smoothed count modes."""
    if nrl_mode not in {"raw", "smoothed"}:
        raise ValueError("--nrl-mode must be raw or smoothed")

    def select(counter: Mapping[int, int]) -> OrderPeak:
        if nrl_mode == "raw":
            return highest_count_distance(counter)
        stats = summarize_distribution(
            counter,
            count_smooth_window=count_smooth_window,
            count_smooth_polyorder=count_smooth_polyorder,
            percent_smooth_window=0,
            percent_smooth_polyorder=2,
        )
        peak_distance = int(stats.smoothed_mode)
        peak_count = int(counter.get(peak_distance, 0))
        return OrderPeak(0, peak_distance, peak_count, int(stats.total_pairs))
    grouped: dict[tuple[str, str], list[OrderPeak]] = defaultdict(list)

    for order in range(1, max_order + 1):
        if include_genome:
            counter = results.genome_all.get(order, Counter())
            if counter:
                selected = select(counter)
                grouped[("combined_chromosomes", ".")].append(
                    OrderPeak(order, selected.peak_distance, selected.peak_count, selected.total_pairs)
                )

        if include_chromosomes:
            for chromosome in sorted(results.chrom_all.get(order, {}), key=natural_sort_key):
                counter = results.chrom_all[order][chromosome]
                if not counter:
                    continue
                selected = select(counter)
                grouped[("chromosome", chromosome)].append(
                    OrderPeak(order, selected.peak_distance, selected.peak_count, selected.total_pairs)
                )

    regressions: list[NRLRegression] = []
    for (scope, chromosome), peaks in sorted(
        grouped.items(),
        key=lambda item: (0 if item[0][0] == "combined_chromosomes" else 1, natural_sort_key(item[0][1])),
    ):
        if len(peaks) < 2:
            continue
        regressions.append(
            regress_order_peaks(peaks, scope=scope, chromosome=chromosome)
        )
    return regressions


def _safe_output_token(value: str) -> str:
    """Return a filesystem-safe token while retaining recognizable contig names."""
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return token or "unnamed"


def _format_regression_value(value: float, decimals: int = 6) -> str:
    """Format finite regression values for TSV output."""
    return f"{value:.{decimals}f}" if math.isfinite(value) else "NA"


def regression_output_stem(prefix: str | Path, regression: NRLRegression) -> Path:
    """Return the output stem for one combined selected-chromosome or chromosome regression."""
    if regression.scope == "combined_chromosomes":
        token = "combined_chromosomes"
    else:
        token = f"chromosome_{_safe_output_token(regression.chromosome)}"
    return Path(f"{prefix}_nrl_regression_{token}")


def write_nrl_regression_tsv(path: str | Path, regression: NRLRegression) -> None:
    """Write selected order peaks, fitted distances, and residuals."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8") as handle:
        handle.write(
            "scope\tchromosome\torder\tpeak_distance_bp\tpeak_count\t"
            "total_pairs\tfitted_distance_bp\tresidual_bp\n"
        )
        for peak in regression.peaks:
            fitted = regression.intercept + regression.slope * peak.order
            residual = peak.peak_distance - fitted
            handle.write(
                f"{regression.scope}\t{regression.chromosome}\t{peak.order}\t"
                f"{peak.peak_distance}\t{peak.peak_count}\t{peak.total_pairs}\t"
                f"{fitted:.6f}\t{residual:.6f}\n"
            )


def plot_nrl_regression(path: str | Path, regression: NRLRegression) -> Path:
    """Plot order-specific peak distances and the fitted NRL regression line."""
    orders = np.asarray([peak.order for peak in regression.peaks], dtype=np.float64)
    distances = np.asarray([peak.peak_distance for peak in regression.peaks], dtype=np.float64)
    fitted = regression.intercept + regression.slope * orders

    figure, axis = plt.subplots(figsize=(6.5, 6.5))
    axis.scatter(
        orders,
        distances,
        s=34,
        facecolors="none",
        edgecolors="black",
        linewidths=1.0,
        label="Highest-count distance",
    )
    axis.plot(
        orders,
        fitted,
        color="black",
        linewidth=1.5,
        linestyle=":",
        label="Linear regression",
    )
    axis.set_xlabel("Neighbour order")
    axis.set_ylabel("Peak distance (bp)")
    scope_label = "Combined chromosomes" if regression.scope == "combined_chromosomes" else regression.chromosome
    axis.set_title(f"{scope_label} nucleosome repeat length regression")
    annotation = (
        f"NRL (slope) = {regression.slope:.3f} bp\n"
        f"R² = {_format_regression_value(regression.r_squared, 4)}"
    )
    axis.text(
        0.04,
        0.96,
        annotation,
        transform=axis.transAxes,
        va="top",
        ha="left",
    )
    from nucleosuite.plotting import apply_integer_x_axis, apply_integer_y_axis
    apply_integer_x_axis(axis, orders)
    apply_integer_y_axis(axis)
    axis.grid(False)
    axis.legend(frameon=False)
    from nucleosuite.plotting import annotate_points, get_plot_options, save_figure
    options = get_plot_options()
    annotate_points(axis, orders, distances, points_are_peaks=True, options=options)
    figure.tight_layout()
    saved = save_figure(figure, path, default_dpi=220, bbox_inches="tight")
    plt.close(figure)
    return saved


def write_nrl_regression_outputs(
    regressions: Sequence[NRLRegression],
    *,
    prefix: str | Path,
) -> list[Path]:
    """Write one TSV/PNG pair per scope and a combined regression summary TSV."""
    if not regressions:
        return []

    prefix = Path(prefix)
    summary_path = Path(f"{prefix}_nrl_regression_summary.tsv")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = [summary_path]

    with summary_path.open("wt", encoding="utf-8") as summary_handle:
        summary_handle.write(
            "scope\tchromosome\torders_used\tmin_order\tmax_order\t"
            "nrl_bp\tintercept_bp\tr_squared\tpoints_tsv\tplot_png\n"
        )
        for regression in regressions:
            stem = regression_output_stem(prefix, regression)
            tsv_path = Path(f"{stem}.tsv")
            from nucleosuite.plotting import plot_path
            png_path = plot_path(Path(f"{stem}.png"))
            write_nrl_regression_tsv(tsv_path, regression)
            png_path = plot_nrl_regression(png_path, regression)
            outputs.extend([tsv_path, png_path])
            orders = [peak.order for peak in regression.peaks]
            summary_handle.write(
                f"{regression.scope}\t{regression.chromosome}\t{len(orders)}\t"
                f"{min(orders)}\t{max(orders)}\t{regression.slope:.6f}\t"
                f"{regression.intercept:.6f}\t"
                f"{_format_regression_value(regression.r_squared)}\t"
                f"{tsv_path}\t{png_path}\n"
            )
    return outputs


def _parse_rgb(value: str) -> tuple[float, float, float]:
    """Convert a UCSC ``R,G,B`` string to a Matplotlib colour tuple."""
    try:
        channels = [int(token.strip()) for token in value.split(",")]
    except ValueError:
        channels = []
    if len(channels) != 3 or any(channel < 0 or channel > 255 for channel in channels):
        return (0.5, 0.5, 0.5)
    return tuple(channel / 255.0 for channel in channels)  # type: ignore[return-value]


def state_colour_map(
    state_indexes: Mapping[str, IntervalIndex],
    *,
    category_rules=(),
) -> dict[str, str]:
    """Return one RGB value per (optionally categorized) state label."""
    colours: dict[str, str] = {}
    for index in state_indexes.values():
        for interval in index.intervals:
            label = categorize_state_name(interval.label, category_rules)
            colours.setdefault(label, interval.rgb)
    return colours


def compute_same_interval_state_counts(
    retained_by_chrom: Mapping[str, Sequence[PeakRecord]],
    state_indexes: Mapping[str, IntervalIndex],
    *,
    category_rules=(),
    min_distance: int,
    max_distance: int,
) -> dict[str, Counter[int]]:
    """Count adjacent peak distances independently inside every state interval.

    Resetting adjacency at each interval boundary prevents a pair from spanning an
    intervening state, even when the endpoint labels happen to be identical.
    """
    counters: dict[str, Counter[int]] = defaultdict(Counter)

    # Canonical aliases are used only when a key is unique. Exact matches always
    # take precedence, and ambiguous files containing both ``1`` and ``chr1``
    # are not silently collapsed.
    canonical_peak_keys: dict[str, str | None] = {}
    for peak_chrom in retained_by_chrom:
        key = canonical_contig_key(peak_chrom)
        if key in canonical_peak_keys and canonical_peak_keys[key] != peak_chrom:
            canonical_peak_keys[key] = None
        else:
            canonical_peak_keys[key] = peak_chrom

    for chrom, index in state_indexes.items():
        peaks = retained_by_chrom.get(chrom)
        if not peaks:
            alias = canonical_peak_keys.get(canonical_contig_key(chrom))
            peaks = retained_by_chrom.get(alias) if alias is not None else None
        if not peaks:
            continue
        positions = np.asarray([record.position for record in peaks], dtype=np.int64)
        for interval in index.intervals:
            left = int(np.searchsorted(positions, interval.start, side="left"))
            right = int(np.searchsorted(positions, interval.end, side="left"))
            if right - left < 2:
                continue
            distances = np.diff(positions[left:right])
            distances = distances[(distances >= min_distance) & (distances <= max_distance)]
            if distances.size == 0:
                continue
            label = categorize_state_name(interval.label, category_rules)
            counters[label].update(int(value) for value in distances)
    return dict(counters)


def _fixed_distribution(
    counter: Mapping[int, int],
    *,
    minimum: int,
    maximum: int,
    window_length: int,
    polyorder: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    distances = np.arange(minimum, maximum + 1, dtype=int)
    counts = np.asarray([counter.get(int(distance), 0) for distance in distances], dtype=float)
    total = float(counts.sum())
    if total <= 0:
        return distances, counts, counts.copy(), math.nan, math.nan, math.nan
    raw_percent = counts * (100.0 / total)
    smoothed = adaptive_savgol(raw_percent, window_length=window_length, polyorder=polyorder)
    smoothed = np.clip(smoothed, 0.0, None)
    smoothed_total = float(smoothed.sum())
    if smoothed_total > 0:
        smoothed *= 100.0 / smoothed_total
    mode = float(distances[int(np.argmax(smoothed))])
    mean = float(np.dot(distances, smoothed) / smoothed.sum())
    cumulative = np.cumsum(smoothed)
    median = float(distances[int(np.searchsorted(cumulative, 50.0, side="left"))])
    return distances, raw_percent, smoothed, mode, mean, median


def write_state_overlay_outputs(
    counters: Mapping[str, Counter[int]],
    colours: Mapping[str, str],
    *,
    prefix: str | Path,
    minimum: int,
    maximum: int,
    smooth_window: int,
    smooth_polyorder: int,
    title: str | None,
    plot_format: str | None,
    x_major_tick: float | None = None,
    x_minor_tick: float | None = None,
) -> list[Path]:
    """Write state-wise relative distributions, smoothed summaries, and overlay plot."""
    prefix = Path(prefix)
    distribution_path = Path(f"{prefix}_state_relative_percent.tsv")
    summary_path = Path(f"{prefix}_state_relative_percent_summary.tsv")
    distribution_path.parent.mkdir(parents=True, exist_ok=True)
    states = sorted((state for state, counter in counters.items() if sum(counter.values()) > 0), key=natural_sort_key)
    if not states:
        raise ValueError("No within-state adjacent peak distances were available for the overlay")

    curves: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    with distribution_path.open("wt") as dist_handle, summary_path.open("wt") as summary_handle:
        dist_handle.write("state\trgb\tdistance_bp\tcount\traw_percent\tsmoothed_percent\n")
        summary_handle.write(
            "state\trgb\tdistance_count\tsmoothing_window\tsmoothing_polyorder\t"
            "smoothed_mode_bp\tsmoothed_mean_bp\tsmoothed_median_bp\n"
        )
        for state in states:
            distances, raw_percent, smoothed, mode, mean, median = _fixed_distribution(
                counters[state],
                minimum=minimum,
                maximum=maximum,
                window_length=smooth_window,
                polyorder=smooth_polyorder,
            )
            curves[state] = (distances, raw_percent, smoothed)
            rgb = colours.get(state, "128,128,128")
            summary_handle.write(
                f"{state}\t{rgb}\t{sum(counters[state].values())}\t{smooth_window}\t"
                f"{smooth_polyorder}\t{mode:.6f}\t{mean:.6f}\t{median:.6f}\n"
            )
            for distance, raw_value, smooth_value in zip(distances, raw_percent, smoothed):
                dist_handle.write(
                    f"{state}\t{rgb}\t{int(distance)}\t{counters[state].get(int(distance), 0)}\t"
                    f"{raw_value:.9f}\t{smooth_value:.9f}\n"
                )

    fig, ax = plt.subplots(figsize=(13, 8))
    for state in states:
        distances, raw_percent, _smoothed = curves[state]
        ax.plot(
            distances,
            raw_percent,
            label=state,
            color=_parse_rgb(colours.get(state, "128,128,128")),
            linewidth=1.8,
            marker="o",
            markersize=1.8,
            markeredgewidth=0,
        )
    ax.set_xlim(minimum, maximum)
    ax.set_xlabel("Adjacent peak distance (bp)")
    ax.set_ylabel("Raw distance frequency within state (%)")
    ax.set_title(title or "Adjacent peak distances by chromatin state")
    from nucleosuite.plotting import apply_distance_x_axis
    apply_distance_x_axis(
        ax,
        major_interval=x_major_tick,
        minor_interval=x_minor_tick,
    )
    ax.legend(frameon=False, ncol=2, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    from dataclasses import replace
    from nucleosuite.plotting import get_plot_options, plot_path, save_figure
    options = get_plot_options()
    if plot_format is not None:
        options = replace(options, format=plot_format)
    fig.tight_layout()
    plot_output = plot_path(Path(f"{prefix}_state_relative_percent.png"), options)
    plot_output = save_figure(fig, plot_output, default_dpi=220, bbox_inches="tight", options=options)
    plt.close(fig)
    return [distribution_path, summary_path, plot_output]


def iter_distributions(
    results: DistanceResults,
    *,
    max_order: int,
    include_chromosomes: bool,
    include_genome: bool,
    include_state_strata: bool,
) -> Iterable[tuple[int, str, str, str, Counter[int]]]:
    """Yield ``order, scope, chromosome, state, counter`` in stable order."""
    for order in range(1, max_order + 1):
        if include_chromosomes:
            for chrom in sorted(results.chrom_all.get(order, {}), key=natural_sort_key):
                pooled = results.chrom_all[order][chrom]
                if pooled:
                    yield order, "chromosome", chrom, "All", pooled

                if include_state_strata:
                    state_counters = results.chrom_state.get(order, {}).get(chrom, {})
                    for state in sorted(state_counters, key=natural_sort_key):
                        counter = state_counters[state]
                        if counter:
                            yield order, "chromosome", chrom, state, counter

        if include_genome:
            pooled = results.genome_all.get(order, Counter())
            if pooled:
                yield order, "combined_chromosomes", ".", "All", pooled

            if include_state_strata:
                for state in sorted(results.genome_state.get(order, {}), key=natural_sort_key):
                    counter = results.genome_state[order][state]
                    if counter:
                        yield order, "combined_chromosomes", ".", state, counter


def write_threshold_metadata(
    output_path: str | Path,
    *,
    input_path: str | Path,
    state_path: str | Path | None,
    blacklist_path: str | Path | None = None,
    blacklisted_peaks: int = 0,
    selection: ThresholdSelection,
    results: DistanceResults,
    score_column: int,
    min_distance: int,
    max_distance: int,
    max_order: int,
    duplicate_policy: str,
    bin_tie_mode: str = "split",
    nrl_mode: str = "smoothed",
) -> None:
    """Write analysis parameters as a two-column TSV file."""
    percentile_mode = (
        "rank_bin"
        if selection.is_rank_bin
        else "bin"
        if selection.is_bin
        else "threshold"
    )
    rows = [
        ("input", input_path),
        ("state_bed", state_path or ""),
        ("blacklist_bed", blacklist_path or ""),
        ("blacklist_overlapping_peaks_excluded", blacklisted_peaks),
        ("score_column", score_column),
        ("requested_percentile", f"{selection.requested_percentile:.12g}"),
        ("effective_percentile", f"{selection.effective_percentile:.12g}"),
        ("percentile_mode", percentile_mode),
        (
            "bin_tie_mode",
            bin_tie_mode if selection.is_bin else "",
        ),
        ("percentile_range", selection.label if selection.is_bin else ""),
        (
            "percentile_lower",
            "" if selection.percentile_lower is None else f"{selection.percentile_lower:.12g}",
        ),
        (
            "percentile_upper",
            "" if selection.percentile_upper is None else f"{selection.percentile_upper:.12g}",
        ),
        (
            "pct_bin_size",
            (
                ""
                if not selection.is_bin
                else f"{selection.percentile_upper - selection.percentile_lower:.12g}"
            ),
        ),
        ("target_peaks", "" if selection.target_peaks is None else selection.target_peaks),
        ("score_threshold", f"{selection.threshold:.12g}"),
        (
            "score_upper_bound",
            "" if selection.score_upper_bound is None else f"{selection.score_upper_bound:.12g}",
        ),
        (
            "rank_start",
            "" if selection.rank_start is None else selection.rank_start,
        ),
        (
            "rank_stop",
            "" if selection.rank_stop is None else selection.rank_stop,
        ),
        (
            "rank_seed",
            "" if selection.rank_seed is None else selection.rank_seed,
        ),
        (
            "rank_score_min",
            f"{selection.threshold:.12g}" if selection.is_rank_bin else "",
        ),
        (
            "rank_score_max",
            (
                ""
                if selection.rank_score_max is None
                else f"{selection.rank_score_max:.12g}"
            ),
        ),
        (
            "rank_tie_order",
            (
                "global_random_shuffle_then_stable_ascending_score_sort"
                if selection.is_rank_bin
                else ""
            ),
        ),
        (
            "rank_membership_rule",
            "exact_zero_based_start_stop_slice" if selection.is_rank_bin else "",
        ),
        ("threshold_pass_count", results.threshold_pass_count),
        ("retained_count", results.retained_count),
        ("duplicate_policy", duplicate_policy),
        ("min_distance", min_distance),
        ("max_distance", max_distance),
        ("max_order", max_order),
        ("nrl_regression", f"{nrl_mode}_count_mode_by_order" if max_order > 1 else "disabled"),
    ]
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt") as handle:
        handle.write("parameter\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")


def write_distribution_outputs(
    results: DistanceResults,
    *,
    distance_path: str | Path,
    summary_path: str | Path,
    max_order: int,
    include_chromosomes: bool,
    include_genome: bool,
    include_state_strata: bool,
    include_zero_distances: bool,
    count_smooth_window: int,
    count_smooth_polyorder: int,
    percent_smooth_window: int,
    percent_smooth_polyorder: int,
) -> tuple[int, int]:
    """Write tidy distance and summary TSV files; return distribution/row counts."""
    distance_path = Path(distance_path)
    summary_path = Path(summary_path)
    distance_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    distributions_written = 0
    rows_written = 0


    with distance_path.open("wt") as distance_handle, summary_path.open("wt") as summary_handle:
        distance_handle.write(
            "order\tscope\tchromosome\tstate\tdistance_bp\tcount\t"
            "smoothed_count\tpercent\tsmoothed_percent\n"
        )

        summary_handle.write(
            "order\tscope\tchromosome\tstate\ttotal_pairs\traw_mode_bp\t"
            "smoothed_mode_bp\tmedian_bp\tmean_bp\tsmoothed_median_bp\t"
            "smoothed_mean_bp\tobserved_min_bp\tobserved_max_bp\n"
        )

        for order, scope, chrom, state, counter in iter_distributions(
            results,
            max_order=max_order,
            include_chromosomes=include_chromosomes,
            include_genome=include_genome,
            include_state_strata=include_state_strata,
        ):
            stats = summarize_distribution(
                counter,
                count_smooth_window=count_smooth_window,
                count_smooth_polyorder=count_smooth_polyorder,
                percent_smooth_window=percent_smooth_window,
                percent_smooth_polyorder=percent_smooth_polyorder,
            )
            distributions_written += 1

            summary_handle.write(
                f"{order}\t{scope}\t{chrom}\t{state}\t{stats.total_pairs}\t"
                f"{stats.raw_mode}\t{stats.smoothed_mode}\t{stats.median:.6f}\t"
                f"{stats.mean:.6f}\t{stats.smoothed_median:.6f}\t"
                f"{stats.smoothed_mean:.6f}\t{stats.observed_min}\t{stats.observed_max}\n"
            )

            for index, distance in enumerate(stats.distances):
                raw_count = int(stats.raw_counts[index])
                if raw_count == 0 and not include_zero_distances:
                    continue
                distance_handle.write(
                    f"{order}\t{scope}\t{chrom}\t{state}\t{int(distance)}\t"
                    f"{raw_count}\t{stats.smoothed_counts[index]:.6f}\t"
                    f"{stats.raw_percent[index]:.9f}\t{stats.smoothed_percent[index]:.9f}\n"
                )
                rows_written += 1

    return distributions_written, rows_written


def plot_distance_distributions(
    distance_path: str | Path,
    output_path: str | Path,
    *,
    nrl_mode: str = "smoothed",
    label_peaks: bool = True,
) -> Path | None:
    """Plot pooled neighbour-order distance distributions and their selected modes."""
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle, get_plot_options
    configure_unique_category_cycle()

    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    with Path(distance_path).open("rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["scope"] == "combined_chromosomes" and row["state"] == "All":
                grouped[int(row["order"])].append(row)
    if not grouped:
        return None

    fig, ax = plt.subplots(figsize=(10, 5.5))
    options = get_plot_options()
    do_labels = label_peaks and options.label_points != "none"
    for order in sorted(grouped):
        rows = sorted(grouped[order], key=lambda row: int(row["distance_bp"]))
        x = np.asarray([int(row["distance_bp"]) for row in rows], dtype=float)
        raw = np.asarray([float(row["count"]) for row in rows], dtype=float)
        smooth = np.asarray([float(row["smoothed_count"]) for row in rows], dtype=float)
        if nrl_mode == "smoothed":
            ax.plot(x, raw, color="0.78", linewidth=0.8, alpha=0.65, zorder=1)
            line, = ax.plot(x, smooth, linewidth=1.5, label=f"+{order}", zorder=2)
            selected = smooth
        else:
            line, = ax.plot(x, raw, linewidth=1.35, label=f"+{order}", zorder=2)
            selected = raw
        if selected.size:
            idx = int(np.argmax(selected))
            ax.scatter([x[idx]], [selected[idx]], s=22, facecolors="white", edgecolors=line.get_color(), zorder=4)
            if do_labels:
                ax.annotate(f"{x[idx]:g}", (x[idx], selected[idx]), xytext=(0, 5), textcoords="offset points", ha="center", va="bottom", fontsize=8)
    from nucleosuite.plotting import apply_base_pair_x_axis, save_figure
    all_x = [float(row["distance_bp"]) for rows in grouped.values() for row in rows]
    apply_base_pair_x_axis(ax, all_x)
    ax.set_xlabel("Distance (bp)")
    ax.set_ylabel("Count")
    ax.legend(frameon=False, title="Neighbour order")
    fig.tight_layout()
    saved = save_figure(fig, output_path, default_dpi=220, bbox_inches="tight")
    plt.close(fig)
    return saved


class PercentilePlotStore:
    """Disk-backed percentile curves used to keep sweep plotting memory bounded."""

    def __init__(self, base_prefix: str | Path) -> None:
        self.base_prefix = Path(base_prefix)
        self.base_prefix.parent.mkdir(parents=True, exist_ok=True)
        temporary = tempfile.NamedTemporaryFile(
            prefix=f".{self.base_prefix.name}.",
            suffix=".percentile_plots.sqlite",
            dir=self.base_prefix.parent,
            delete=False,
        )
        self.path = Path(temporary.name)
        temporary.close()
        self.connection = sqlite3.connect(self.path)
        self.connection.executescript(
            """
            CREATE TABLE peak_counts (
                neighbour_order INTEGER NOT NULL,
                scope TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                percentile REAL NOT NULL,
                percentile_label TEXT NOT NULL,
                percentile_lower REAL,
                percentile_upper REAL,
                score_threshold REAL NOT NULL,
                score_upper_bound REAL,
                percentile_mode TEXT NOT NULL,
                bin_tie_mode TEXT NOT NULL,
                rank_start INTEGER,
                rank_stop INTEGER,
                rank_seed INTEGER,
                rank_score_max REAL,
                retained_peak_count INTEGER NOT NULL,
                PRIMARY KEY (neighbour_order, scope, chromosome, percentile)
            );
            CREATE TABLE curves (
                neighbour_order INTEGER NOT NULL,
                scope TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                percentile REAL NOT NULL,
                distance_bp INTEGER NOT NULL,
                raw_count REAL NOT NULL,
                raw_percent REAL NOT NULL
            );
            CREATE INDEX curves_identity
                ON curves (neighbour_order, scope, chromosome, percentile, distance_bp);
            """
        )

    def add_threshold(
        self,
        *,
        selection: ThresholdSelection,
        results: DistanceResults,
        chromosomes: Sequence[str],
        max_order: int,
        include_chromosomes: bool,
        include_genome: bool,
    ) -> None:
        """Persist one threshold's raw pooled curves and retained-peak counts."""
        percentile = float(selection.effective_percentile)
        count_rows: list[
            tuple[
                int,
                str,
                str,
                float,
                str,
                float | None,
                float | None,
                float,
                float | None,
                str,
                str,
                int | None,
                int | None,
                int | None,
                float | None,
                int,
            ]
        ] = []
        for order in range(1, max_order + 1):
            if include_genome:
                count_rows.append(
                    (
                        order,
                        "combined_chromosomes",
                        ".",
                        percentile,
                        selection.label,
                        selection.percentile_lower,
                        selection.percentile_upper,
                        float(selection.threshold),
                        selection.score_upper_bound,
                        "rank_bin" if selection.is_rank_bin else "bin" if selection.is_bin else "threshold",
                        "split" if selection.is_rank_bin else "keep" if selection.is_bin else "",
                        selection.rank_start,
                        selection.rank_stop,
                        selection.rank_seed,
                        selection.rank_score_max,
                        int(results.retained_count),
                    )
                )
            if include_chromosomes:
                for chromosome in chromosomes:
                    count_rows.append(
                        (
                            order,
                            "chromosome",
                            chromosome,
                            percentile,
                            selection.label,
                            selection.percentile_lower,
                            selection.percentile_upper,
                            float(selection.threshold),
                            selection.score_upper_bound,
                            "rank_bin" if selection.is_rank_bin else "bin" if selection.is_bin else "threshold",
                            "split" if selection.is_rank_bin else "keep" if selection.is_bin else "",
                            selection.rank_start,
                            selection.rank_stop,
                            selection.rank_seed,
                            selection.rank_score_max,
                            len(results.retained_by_chrom.get(chromosome, ())),
                        )
                    )
        self.connection.executemany(
            "INSERT INTO peak_counts VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            count_rows,
        )

        for order, scope, chromosome, _state, counter in iter_distributions(
            results,
            max_order=max_order,
            include_chromosomes=include_chromosomes,
            include_genome=include_genome,
            include_state_strata=False,
        ):
            minimum = min(counter)
            maximum = max(counter)
            total = float(sum(counter.values()))
            self.connection.executemany(
                "INSERT INTO curves VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    (
                        order,
                        scope,
                        chromosome,
                        percentile,
                        distance,
                        float(counter.get(distance, 0)),
                        float(counter.get(distance, 0)) * 100.0 / total,
                    )
                    for distance in range(minimum, maximum + 1)
                ),
            )
        self.connection.commit()

    def write_outputs(self) -> list[Path]:
        """Render count, percentage, and peak-count figures for every scope/order."""
        outputs: list[Path] = []
        counts_tsv = Path(f"{self.base_prefix}_percentile_sweep_peak_counts.tsv")
        with counts_tsv.open("wt", encoding="utf-8") as handle:
            handle.write(
                "order\tscope\tchromosome\tpercentile_threshold\t"
                "score_threshold\tretained_peak_count\tpercentile_mode\t"
                "bin_tie_mode\tpercentile_range\tpercentile_lower\tpercentile_upper\t"
                "score_upper_bound\trank_start\trank_stop\trank_seed\t"
                "rank_score_max\n"
            )
            for row in self.connection.execute(
                """
                SELECT neighbour_order, scope, chromosome, percentile,
                       score_threshold, retained_peak_count, percentile_label,
                       percentile_lower, percentile_upper, score_upper_bound,
                       percentile_mode, bin_tie_mode, rank_start, rank_stop, rank_seed,
                       rank_score_max
                FROM peak_counts
                ORDER BY neighbour_order,
                         CASE scope WHEN 'combined_chromosomes' THEN 0 ELSE 1 END,
                         chromosome, percentile
                """
            ):
                (
                    order,
                    scope,
                    chromosome,
                    percentile,
                    threshold,
                    peak_count,
                    label,
                    lower,
                    upper,
                    upper_score,
                    mode,
                    bin_tie_mode,
                    rank_start,
                    rank_stop,
                    rank_seed,
                    rank_score_max,
                ) = row
                range_label = label if mode in {"bin", "rank_bin"} else ""
                handle.write(
                    f"{order}\t{scope}\t{chromosome}\t{percentile:.12g}\t"
                    f"{threshold:.12g}\t{peak_count}\t{mode}\t"
                    f"{bin_tie_mode}\t"
                    f"{range_label}\t"
                    f"{'' if lower is None else f'{lower:.12g}'}\t"
                    f"{'' if upper is None else f'{upper:.12g}'}\t"
                    f"{'' if upper_score is None else f'{upper_score:.12g}'}\t"
                    f"{'' if rank_start is None else rank_start}\t"
                    f"{'' if rank_stop is None else rank_stop}\t"
                    f"{'' if rank_seed is None else rank_seed}\t"
                    f"{'' if rank_score_max is None else f'{rank_score_max:.12g}'}\n"
                )
        outputs.append(counts_tsv)

        identities = list(
            self.connection.execute(
                """
                SELECT DISTINCT neighbour_order, scope, chromosome
                FROM peak_counts
                ORDER BY neighbour_order,
                         CASE scope WHEN 'combined_chromosomes' THEN 0 ELSE 1 END,
                         chromosome
                """
            )
        )
        for order, scope, chromosome in identities:
            count_rows = list(
                self.connection.execute(
                    """
                    SELECT percentile, retained_peak_count, percentile_label,
                           percentile_lower
                    FROM peak_counts
                    WHERE neighbour_order = ? AND scope = ? AND chromosome = ?
                    ORDER BY percentile
                    """,
                    (order, scope, chromosome),
                )
            )
            curve_rows = list(
                self.connection.execute(
                    """
                    SELECT percentile, distance_bp, raw_count, raw_percent
                    FROM curves
                    WHERE neighbour_order = ? AND scope = ? AND chromosome = ?
                    ORDER BY percentile, distance_bp
                    """,
                    (order, scope, chromosome),
                )
            )
            stem = percentile_plot_stem(
                self.base_prefix,
                order=int(order),
                scope=str(scope),
                chromosome=str(chromosome),
            )
            identity_counts_path = Path(f"{stem}_peak_counts.tsv")
            with identity_counts_path.open("wt", encoding="utf-8") as handle:
                handle.write(
                    "order\tscope\tchromosome\tpercentile_threshold\t"
                    "retained_peak_count\tpercentile_range\tpercentile_lower\n"
                )
                for percentile, peak_count, label, lower in count_rows:
                    handle.write(
                        f"{order}\t{scope}\t{chromosome}\t{float(percentile):.12g}\t"
                        f"{int(peak_count)}\t{str(label)}\t"
                        f"{'' if lower is None else f'{float(lower):.12g}'}\n"
                    )
            outputs.append(identity_counts_path)

            identity_curves_path = Path(f"{stem}_curves.tsv")
            with identity_curves_path.open("wt", encoding="utf-8") as handle:
                handle.write(
                    "order\tscope\tchromosome\tpercentile_threshold\t"
                    "percentile_range\tpercentile_lower\tdistance_bp\t"
                    "raw_count\traw_percent\n"
                )
                count_metadata = {
                    float(percentile): (str(label), lower)
                    for percentile, _peak_count, label, lower in count_rows
                }
                for percentile, distance, raw_count, raw_percent in curve_rows:
                    label, lower = count_metadata[float(percentile)]
                    handle.write(
                        f"{order}\t{scope}\t{chromosome}\t{float(percentile):.12g}\t"
                        f"{label}\t{'' if lower is None else f'{float(lower):.12g}'}\t"
                        f"{int(distance)}\t{float(raw_count):.12g}\t"
                        f"{float(raw_percent):.12g}\n"
                    )
            outputs.append(identity_curves_path)

            table_path = Path(f"{stem}_peak_counts.png")
            table_path = plot_percentile_peak_count_table(count_rows, table_path)
            outputs.append(table_path)
            if curve_rows:
                count_path = Path(f"{stem}_count.png")
                percentage_path = Path(f"{stem}_percentage.png")
                count_path, percentage_path = plot_percentile_distance_curves(
                    curve_rows,
                    all_percentiles=[float(row[0]) for row in count_rows],
                    all_labels=[str(row[2]) for row in count_rows],
                    percentile_bins=any(row[3] is not None for row in count_rows),
                    count_path=count_path,
                    percentage_path=percentage_path,
                    order=int(order),
                    scope=str(scope),
                    chromosome=str(chromosome),
                )
                outputs.extend([count_path, percentage_path])
        return outputs

    def close(self) -> None:
        self.connection.close()
        self.path.unlink(missing_ok=True)


def percentile_plot_stem(
    base_prefix: str | Path,
    *,
    order: int,
    scope: str,
    chromosome: str,
) -> Path:
    """Return a stable output stem for one percentile-sweep plot set."""
    if scope == "combined_chromosomes":
        scope_token = "combined_chromosomes"
    else:
        scope_token = f"chromosome_{_safe_output_token(chromosome)}"
    return Path(f"{base_prefix}_percentile_sweep_order{order}_{scope_token}")


def _percentile_colour_map(percentiles: Sequence[float]):
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    colour_map = LinearSegmentedColormap.from_list(
        "nucleosuite_percentile",
        ("#ff1f1f", "#a00070", "#2020ff"),
    )
    minimum = min(percentiles)
    maximum = max(percentiles)
    if maximum <= minimum:
        maximum = minimum + 1.0
    return colour_map, Normalize(vmin=minimum, vmax=maximum)


def _scope_plot_label(scope: str, chromosome: str, order: int) -> str:
    scope_label = "Combined chromosomes" if scope == "combined_chromosomes" else chromosome
    return f"{scope_label}; neighbour order +{order}"


def plot_percentile_distance_curves(
    rows: Sequence[tuple[float, int, float, float]],
    *,
    all_percentiles: Sequence[float] | None = None,
    all_labels: Sequence[str] | None = None,
    percentile_bins: bool = False,
    count_path: str | Path,
    percentage_path: str | Path,
    order: int,
    scope: str,
    chromosome: str,
) -> tuple[Path, Path]:
    """Render threshold-coloured raw count and raw percentage overlays."""
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import FuncFormatter
    from nucleosuite.plotting import apply_base_pair_x_axis

    curves: dict[float, tuple[list[int], list[float], list[float]]] = {}
    for percentile, distance, raw_count, raw_percent in rows:
        distances, counts, percentages = curves.setdefault(
            float(percentile), ([], [], [])
        )
        distances.append(int(distance))
        counts.append(float(raw_count))
        percentages.append(float(raw_percent))

    percentiles = sorted(curves)
    colour_scale_percentiles = sorted(all_percentiles or percentiles)
    colour_map, normalizer = _percentile_colour_map(colour_scale_percentiles)
    scalar = ScalarMappable(norm=normalizer, cmap=colour_map)
    scalar.set_array([])
    scope_label = _scope_plot_label(scope, chromosome, order)

    for output_path, value_index, ylabel, title in (
        (
            Path(count_path),
            1,
            "Count",
            f"Peak-distance counts by score threshold\n{scope_label}",
        ),
        (
            Path(percentage_path),
            2,
            "Percentage of distances (%)",
            f"Peak-distance percentages by score threshold\n{scope_label}",
        ),
    ):
        figure, axis = plt.subplots(figsize=(10.5, 6.2))
        all_distances: list[int] = []
        for percentile in percentiles:
            distances, counts, percentages = curves[percentile]
            values = counts if value_index == 1 else percentages
            all_distances.extend(distances)
            axis.plot(
                distances,
                values,
                color=colour_map(normalizer(percentile)),
                linewidth=1.05,
            )
        apply_base_pair_x_axis(axis, all_distances)
        axis.set_xlabel("Distance (bp) between flanking peaks")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(axis="x", alpha=0.35)
        if value_index == 1:
            axis.yaxis.set_major_formatter(
                FuncFormatter(
                    lambda value, _position: (
                        f"{value / 1_000_000:g}M"
                        if abs(value) >= 1_000_000
                        else f"{value / 1_000:g}K"
                        if abs(value) >= 1_000
                        else f"{value:g}"
                    )
                )
            )
        colour_bar = figure.colorbar(scalar, ax=axis, pad=0.025)
        colour_bar.set_label(
            "Peak score percentile range" if percentile_bins
            else "Peak score threshold percentile"
        )
        if len(colour_scale_percentiles) <= 12:
            colour_bar.set_ticks(colour_scale_percentiles)
            if all_labels is not None:
                colour_bar.set_ticklabels(list(all_labels))
        from nucleosuite.plotting import save_figure
        figure.tight_layout()
        saved = save_figure(figure, output_path, default_dpi=220, bbox_inches="tight")
        if value_index == 1:
            count_path = saved
        else:
            percentage_path = saved
        plt.close(figure)
    return Path(count_path), Path(percentage_path)


def plot_percentile_peak_count_table(
    rows: Sequence[tuple[float, int, str, float | None]], output_path: str | Path
) -> Path:
    """Write percentile thresholds and retained peak counts as a standalone image."""
    percentiles = [float(row[0]) for row in rows]
    percentile_bins = any(row[3] is not None for row in rows)
    colour_map, normalizer = _percentile_colour_map(percentiles)
    figure_height = max(2.4, 0.34 * len(rows) + 1.25)
    figure, axis = plt.subplots(figsize=(6.2, figure_height))
    axis.axis("off")
    cell_text = [
        ["━━", str(label), f"{int(count):,}"]
        for percentile, count, label, _lower in rows
    ]
    table = axis.table(
        cellText=cell_text,
        colLabels=(
            "",
            "Percentile range" if percentile_bins else "Percentile threshold",
            "Retained peak count",
        ),
        cellLoc="right",
        colLoc="center",
        loc="center",
        colWidths=(0.12, 0.42, 0.46),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.22)
    for row_index, (percentile, _count, _label, _lower) in enumerate(rows, start=1):
        swatch = table[(row_index, 0)]
        swatch.get_text().set_color(colour_map(normalizer(float(percentile))))
        swatch.get_text().set_ha("center")
    for column in range(3):
        table[(0, column)].get_text().set_weight("bold")
    figure.tight_layout()
    from nucleosuite.plotting import save_figure
    saved = save_figure(figure, output_path, default_dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return saved


def write_duplicate_report(
    duplicates: Mapping[str, Mapping[int, Sequence[PeakRecord]]],
    output_path: str | Path,
) -> int:
    """Write retained duplicate positions and return the number of positions."""
    duplicate_count = sum(len(positions) for positions in duplicates.values())
    if duplicate_count == 0:
        return 0

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt") as handle:
        handle.write("chromosome\tposition\tn_records\tscores\tstates\n")
        for chrom in sorted(duplicates, key=natural_sort_key):
            for position in sorted(duplicates[chrom]):
                records = duplicates[chrom][position]
                scores = ",".join(f"{record.score:.12g}" for record in records)
                states = ",".join(record.state for record in records)
                handle.write(
                    f"{chrom}\t{position}\t{len(records)}\t{scores}\t{states}\n"
                )
    return duplicate_count


def write_filtered_bed(
    retained_by_chrom: Mapping[str, Sequence[PeakRecord]],
    output_path: str | Path,
    chrom_sizes: Mapping[str, int] | Sequence[tuple[str, int]] | str | Path | None = None,
) -> int:
    """Write retained peaks as UCSC BED8 records in the analysis namespace."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    canonical_names = [name for name, _length in normalise_chrom_sizes(chrom_sizes)] if chrom_sizes else []

    with output_path.open("wt") as handle:
        for source_chrom in sorted(retained_by_chrom, key=natural_sort_key):
            chrom = source_chrom
            if canonical_names:
                try:
                    chrom = resolve_contig_name(
                        source_chrom, canonical_names, source_label="chromosome sizes"
                    )
                except KeyError as exc:
                    raise ValueError(str(exc)) from exc
            for record in retained_by_chrom[source_chrom]:
                bed_score = max(0, min(1000, int(round(record.score))))
                thick_start = min(max(record.position, record.start), record.end - 1)
                handle.write(
                    f"{chrom}\t{record.start}\t{record.end}\t{record.name}\t"
                    f"{bed_score}\t{record.strand}\t{thick_start}\t{thick_start + 1}\n"
                )
                written += 1
    return written

def percentile_label(percentile: float) -> str:
    """Format a percentile compactly for filenames."""
    return f"{percentile:.6f}".rstrip("0").rstrip(".")


def default_output_prefix(
    input_path: str | Path,
    state_path: str | Path | None,
) -> Path:
    """Build an automatic output prefix from input and optional state BED names."""
    input_stem = strip_known_suffix(input_path)
    if state_path is None:
        return Path(input_stem)
    return Path(f"{input_stem}_{strip_known_suffix(state_path)}")


def output_grouping_token(args: argparse.Namespace) -> str:
    """Describe the complete score-group definition for output filenames."""

    if args.pct_bin_size is not None:
        return f"binsize-{args.pct_bin_size:g}"
    if args.pct_values:
        raw = "-".join(str(value).replace(",", "-") for value in args.pct_values)
        return f"{'bins' if args.pct_bins else 'values'}-{raw}"
    if args.pct_range:
        return f"range-{args.pct_lower:g}-{args.pct_upper:g}-{args.pct_step:g}"
    if args.target_peaks is not None:
        return f"target-{args.target_peaks}"
    return "threshold"


def validate_smoothing_arguments(args: argparse.Namespace) -> None:
    """Validate Savitzky-Golay CLI settings before processing input."""
    for option, value in (
        ("--count-smooth-window", args.count_smooth_window),
        ("--percent-smooth-window", args.percent_smooth_window),
    ):
        if value < 0 or (value > 0 and value % 2 == 0):
            raise ValueError(f"{option} must be 0 or a positive odd integer")

    for option, value in (
        ("--count-smooth-polyorder", args.count_smooth_polyorder),
        ("--percent-smooth-polyorder", args.percent_smooth_polyorder),
    ):
        if value < 0:
            raise ValueError(f"{option} must be non-negative")

    if args.count_smooth_window and args.count_smooth_polyorder >= args.count_smooth_window:
        raise ValueError("--count-smooth-polyorder must be smaller than --count-smooth-window")
    if args.percent_smooth_window and args.percent_smooth_polyorder >= args.percent_smooth_window:
        raise ValueError("--percent-smooth-polyorder must be smaller than --percent-smooth-window")


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add ``nucleosuite distances`` arguments to an argparse parser."""
    parser.add_argument(
        "input",
        help="Input peak BED/BED.gz file with chromosome, start, end, and score.",
    )
    parser.add_argument(
        "--blacklist-bed",
        help="BED blacklist; complete overlapping peak records are excluded.",
    )
    parser.add_argument(
        "--state-bed",
        help="Optional chromatin-state BED or BED.gz file.",
    )
    parser.add_argument(
        "--state-label-column",
        type=int,
        default=4,
        help="One-based state-label column (default: 4).",
    )
    parser.add_argument(
        "--state-color-column",
        type=int,
        default=9,
        help="One-based RGB colour column used by the state overlay plot (default: 9).",
    )
    parser.add_argument(
        "--state-overlay-plot",
        action="store_true",
        help=(
            "Plot adjacent distances within each individual state interval as relative "
            "percentages on one colour-matched overlay. Requires --state-bed."
        ),
    )
    parser.add_argument(
        "--state-overlay-title",
        help="Optional title for the ChromHMM state distance overlay plot.",
    )
    parser.add_argument(
        "--state-overlay-format",
        choices=("png", "svg"),
        default=None,
        help=("Figure format for the state-overlay plot. When omitted, the shared "
              "plot format is used."),
    )
    parser.add_argument(
        "--state-overlay-smooth-window",
        type=int,
        default=21,
        help=(
            "Odd Savitzky-Golay window for diagnostic state TSV/summary values; "
            "the plotted percentages remain raw (default: 21)."
        ),
    )
    parser.add_argument(
        "--state-overlay-smooth-polyorder",
        type=int,
        default=2,
        help=(
            "Polynomial order for diagnostic state TSV/summary values; the plot "
            "remains raw (default: 2)."
        ),
    )
    parser.add_argument(
        "--state-overlay-x-major-tick",
        type=float,
        default=None,
        help="Major x-axis tick interval in bp for the state overlay; default is automatic.",
    )
    parser.add_argument(
        "--state-overlay-x-minor-tick",
        type=float,
        default=None,
        help=(
            "Minor x-axis tick interval in bp. Automatic defaults use 5 bp for a "
            "10 bp major interval, 10 bp when the major interval is greater than 50 bp, "
            "and half the major interval otherwise."
        ),
    )
    parser.add_argument(
        "--category",
        action="append",
        default=[],
        metavar="NAME=MATCHER[,MATCHER...]",
        help=(
            "Group state labels using exact:, prefix:, or regex: matchers. "
            "May be repeated; unmatched labels retain their original names."
        ),
    )
    parser.add_argument(
        "--position-column",
        type=int,
        help="One-based explicit peak-position column; default uses the BED midpoint.",
    )
    parser.add_argument(
        "--score-column",
        type=int,
        default=5,
        help="One-based numeric peak-score column (default: 5, the BED score field).",
    )
    parser.add_argument(
        "--score-percentile",
        type=float,
        default=0.0,
        help="Retain scores at or above this percentile (default: 0).",
    )
    parser.add_argument(
        "--target-peaks",
        type=int,
        help=(
            "Choose a threshold intended to retain approximately this many peaks. "
            "Score ties can cause the threshold-pass count to be larger."
        ),
    )
    parser.add_argument(
        "--pct-range",
        action="store_true",
        help=(
            "Analyse a percentile sweep using --pct-values or the inclusive "
            "--pct-lower/--pct-upper/--pct-step sequence."
        ),
    )
    parser.add_argument(
        "--pct-values",
        nargs="+",
        metavar="PCT[,PCT...]",
        help=(
            "Explicit score percentiles, for example 10,20,50,90,99. This "
            "activates the sweep and overrides lower/upper/step."
        ),
    )
    parser.add_argument(
        "--pct-bins",
        action="store_true",
        help=(
            "Interpret selected percentiles as consecutive non-cumulative "
            "boundaries; for example 10,20,30 gives 10-20 and 20-30. "
            "Bin membership follows --bin-tie-mode. Without this option, "
            "thresholds remain cumulative."
        ),
    )
    parser.add_argument(
        "--bin-tie-mode",
        choices=("split", "keep"),
        default="split",
        help=(
            "How equal scores are handled at bin boundaries: 'split' makes "
            "rank-based groups containing the requested proportion of peaks "
            "and splits tied scores when needed after reproducible random tie "
            "ordering; 'keep' keeps identical scores together, so bin sizes "
            "may differ from the requested proportion (default: split)."
        ),
    )
    parser.add_argument(
        "--pct-bin-size",
        type=float,
        metavar="PERCENT",
        help=(
            "Create consecutive bins spanning this many percentile points; for "
            "example 1 creates 100 groups and 5 creates 20. Membership follows "
            "--bin-tie-mode."
        ),
    )
    parser.add_argument(
        "--pct-bin-seed",
        type=int,
        default=1,
        help=(
            "Random seed used to order tied scores before split-mode rank "
            "groups are cut (default: 1)."
        ),
    )
    parser.add_argument("--pct-lower", type=float, default=0.0, help="Sweep lower percentile (default: 0).")
    parser.add_argument("--pct-upper", type=float, default=99.0, help="Sweep upper percentile (default: 99).")
    parser.add_argument("--pct-step", type=float, default=1.0, help="Sweep step (default: 1).")
    parser.add_argument(
        "--min-distance",
        type=int,
        default=1,
        help="Minimum non-zero peak distance to count in bp (default: 1).",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=10000,
        help="Maximum peak distance to count in bp (default: 10000).",
    )
    parser.add_argument(
        "--max-order",
        type=int,
        default=10,
        help=("Maximum neighbour order; 1 gives adjacent peaks only. Values above 1 "
              "also produce neighbour-order peak-distance regressions for NRL estimation (default: 10)."),
    )
    parser.add_argument(
        "--nrl-mode",
        choices=("raw", "smoothed"),
        default="smoothed",
        help=("Distance-distribution mode used for neighbour-order NRL regression. "
              "'smoothed' uses the Savitzky-Golay count curve controlled by "
              "--count-smooth-window/--count-smooth-polyorder; 'raw' uses the raw count mode "
              "(default: smoothed)."),
    )
    parser.add_argument(
        "--duplicate-policy",
        choices=("highest-score", "first", "keep", "error"),
        default="highest-score",
        help=(
            "How to handle retained peaks sharing a position. The default prevents "
            "duplicates from consuming neighbour orders (default: highest-score)."
        ),
    )
    parser.add_argument(
        "--scope",
        choices=("all", "combined_chromosomes", "genome", "chromosome"),
        default="all",
        help="Write pooled combined-chromosome, per-chromosome, or both scopes (default: all).",
    )
    parser.add_argument(
        "--include-zero-distances",
        action="store_true",
        help="Emit zero-count rows between the minimum and maximum observed distances.",
    )
    parser.add_argument(
        "--count-smooth-window",
        type=int,
        default=21,
        help=(
            "Odd Savitzky-Golay window for diagnostic count columns/statistics; "
            "percentile overlays remain raw; 0 disables (default: 21)."
        ),
    )
    parser.add_argument(
        "--count-smooth-polyorder",
        type=int,
        default=2,
        help=(
            "Polynomial order for diagnostic count columns/statistics; percentile "
            "overlays remain raw (default: 2)."
        ),
    )
    parser.add_argument(
        "--percent-smooth-window",
        type=int,
        default=21,
        help=(
            "Odd Savitzky-Golay window for diagnostic percentage columns; all "
            "percentage plots remain raw; 0 disables (default: 21)."
        ),
    )
    parser.add_argument(
        "--percent-smooth-polyorder",
        type=int,
        default=2,
        help=(
            "Polynomial order for diagnostic percentage columns; plots remain "
            "raw (default: 2)."
        ),
    )
    parser.add_argument(
        "--write-filtered-bed",
        action="store_true",
        help="Write retained peaks as standard BED8 and/or bigBed according to --interval-format.",
    )
    parser.add_argument(
        "--interval-format", choices=INTERVAL_FORMATS, default="bed",
        help="Filtered interval output: BED, bigBed, or both.",
    )
    parser.add_argument(
        "--interval-chrom-sizes",
        help="Chromosome sizes required when filtered output includes bigBed.",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        help="Output path prefix; default combines the peak and optional state BED basenames.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Stop at the first malformed data line instead of reporting and skipping it.",
    )
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive peak-spacing "
            "analysis (default: 1; independent of suite --cores)."
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    """Create the parser for ``nucleosuite distances``."""
    parser = argparse.ArgumentParser(
        prog="nucleosuite distances",
        description=(
            "Measure +1 through +N peak distances, optionally stratified by "
            "chromatin state and filtered by score percentile."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    add_arguments(parser)
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def register_subcommand(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Register this command with the central NucleoSuite argparse CLI."""
    parser = subparsers.add_parser(
        "distances",
        help="Measure adjacent and higher-order peak distances.",
        description=(
            "Measure +1 through +N peak distances, optionally stratified by "
            "chromatin state and filtered by score percentile."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    add_arguments(parser)
    parser.set_defaults(func=run)
    return parser


def _run_serial(args: argparse.Namespace) -> int:
    """Execute one distances analysis."""
    reporter = ProgressReporter("distances")
    percentile_store: PercentilePlotStore | None = None
    try:
        state_path = args.state_bed

        if args.min_distance < 0:
            raise ValueError("--min-distance must be at least 0")
        if args.max_distance < args.min_distance:
            raise ValueError("--max-distance must be greater than or equal to --min-distance")
        if args.max_order < 1:
            raise ValueError("--max-order must be at least 1")
        percentile_sweep = bool(
            args.pct_range or args.pct_values or args.pct_bin_size is not None
        )
        if args.pct_bins and not percentile_sweep:
            raise ValueError("--pct-bins requires --pct-values or --pct-range")
        if args.target_peaks is not None and percentile_sweep:
            raise ValueError(
                "--target-peaks and percentile sweeps are mutually exclusive"
            )
        if args.write_filtered_bed and args.interval_format != "bed" and not args.interval_chrom_sizes:
            raise ValueError(
                "--interval-chrom-sizes is required when filtered output includes bigBed"
            )

        validate_one_based_column(args.score_column, "--score-column")
        if args.position_column is not None:
            validate_one_based_column(args.position_column, "--position-column")
        validate_one_based_column(args.state_label_column, "--state-label-column")
        validate_one_based_column(args.state_color_column, "--state-color-column")
        if args.state_overlay_plot and state_path is None:
            raise ValueError("--state-overlay-plot requires --state-bed")
        validate_smoothing_arguments(args)
        if args.state_overlay_smooth_window < 3 or args.state_overlay_smooth_window % 2 == 0:
            raise ValueError("--state-overlay-smooth-window must be an odd integer of at least 3")
        if args.state_overlay_smooth_polyorder < 0 or args.state_overlay_smooth_polyorder >= args.state_overlay_smooth_window:
            raise ValueError("--state-overlay-smooth-polyorder must be non-negative and smaller than the window")
        from nucleosuite.plotting import validate_tick_interval
        validate_tick_interval(args.state_overlay_x_major_tick, "--state-overlay-x-major-tick")
        validate_tick_interval(args.state_overlay_x_minor_tick, "--state-overlay-x-minor-tick")

        category_rules = parse_state_category_rules(args.category)
        state_indexes: dict[str, IntervalIndex] | None = None
        state_summary: ParseSummary | None = None
        if state_path is not None:
            reporter.file_start("chromatin states", state_path)
            state_indexes, state_summary = build_state_indexes(
                state_path,
                label_column=args.state_label_column,
                color_column=args.state_color_column,
                strict=args.strict,
                progress=reporter,
            )
            reporter.stage(
                f"Loaded chromatin states: {state_summary.used_lines:,} intervals "
                f"across {len(state_indexes):,} contigs; "
                f"skipped {state_summary.skipped_lines:,} lines"
            )

        if args.blacklist_bed:
            reporter.stage(f"Loading blacklist: {args.blacklist_bed}")
        blacklist = load_blacklist_unbounded(args.blacklist_bed)
        reporter.file_start("peaks", args.input)
        peaks_by_chrom, scores, peak_summary = load_peaks(
            args.input,
            state_indexes=state_indexes,
            category_rules=category_rules,
            position_column=args.position_column,
            score_column=args.score_column,
            strict=args.strict,
            blacklist=blacklist,
            progress=reporter,
        )
        reporter.stage(
            f"Loaded {peak_summary.used_lines:,} peaks across "
            f"{len(peaks_by_chrom):,} contigs; skipped "
            f"{peak_summary.skipped_lines:,} "
            f"({peak_summary.blacklisted_lines:,} blacklisted)"
        )
        input_chromosomes = sorted(peaks_by_chrom, key=natural_sort_key)

        reporter.stage("Selecting score thresholds")
        percentage_bins = (
            args.bin_tie_mode == "split"
            and (args.pct_bin_size is not None or args.pct_bins)
        )
        if percentage_bins:
            reporter.stage(
                f"Randomizing global peak order with seed {args.pct_bin_seed}, "
                "stably sorting by score, and creating percentage rank groups"
            )
        selections = choose_thresholds(
            scores,
            score_percentile=args.score_percentile,
            target_peaks=args.target_peaks,
            pct_range=args.pct_range,
            pct_lower=args.pct_lower,
            pct_upper=args.pct_upper,
            pct_step=args.pct_step,
            pct_values=args.pct_values,
            pct_bins=args.pct_bins,
            pct_bin_size=args.pct_bin_size,
            pct_bin_seed=args.pct_bin_seed,
            bin_tie_mode=args.bin_tie_mode,
        )
        rank_membership = attach_rank_bin_records(
            peaks_by_chrom,
            selections,
            scores=scores,
        )
        del scores
        if rank_membership is not None:
            # Membership groups now own the same record references in their exact
            # rank partitions, so the original list containers are redundant.
            peaks_by_chrom.clear()

        base_prefix = (
            Path(args.output_prefix)
            if args.output_prefix
            else default_output_prefix(args.input, state_path)
        )
        from nucleosuite.output_naming import parameterized_prefix

        grouping = output_grouping_token(args)
        base_prefix = parameterized_prefix(
            base_prefix,
            (
                ("distmin", args.min_distance),
                ("distmax", args.max_distance),
                ("orders", args.max_order),
                ("nrlmode", args.nrl_mode),
                ("countsg", f"{args.count_smooth_window}x{args.count_smooth_polyorder}"),
                ("pctsg", f"{args.percent_smooth_window}x{args.percent_smooth_polyorder}"),
                ("groups", grouping),
                ("ties", args.bin_tie_mode),
            ),
        )
        base_prefix.parent.mkdir(parents=True, exist_ok=True)

        include_chromosomes = args.scope in {"all", "chromosome"}
        include_genome = args.scope in {"all", "combined_chromosomes", "genome"}
        include_state_strata = state_path is not None
        if percentile_sweep:
            percentile_store = PercentilePlotStore(base_prefix)

        for index, selection in enumerate(selections, start=1):
            label = selection.label
            threshold_prefix = Path(f"{base_prefix}_scorepct{label}")

            if selection.is_rank_bin:
                selection_description = (
                    f"percentage bin={selection.label}, ranks "
                    f"{selection.rank_start:,}:{selection.rank_stop:,}, "
                    f"observed score range {selection.threshold:.12g} to "
                    f"{selection.rank_score_max:.12g}"
                )
                selection_kind = "percentage bin"
            elif selection.is_bin:
                score_description = (
                    f"{selection.threshold:.12g}<=score"
                    if selection.score_upper_bound is None
                    else (
                        f"{selection.threshold:.12g}<=score<"
                        f"{selection.score_upper_bound:.12g}"
                    )
                )
                selection_description = f"bin={selection.label}, {score_description}"
                selection_kind = "percentile bin"
            else:
                selection_description = (
                    f"percentile={selection.effective_percentile:.6g}, "
                    f"score>={selection.threshold:.12g}"
                )
                selection_kind = "threshold"

            reporter.stage(
                f"{selection_kind.title()} {index}/{len(selections)}: "
                f"{selection_description}; calculating distances"
            )

            results = compute_distance_counts(
                peaks_by_chrom,
                threshold=selection.threshold,
                score_upper_bound=selection.score_upper_bound,
                rank_membership=selection.rank_membership,
                rank_bin_index=selection.rank_bin_index,
                min_distance=args.min_distance,
                max_distance=args.max_distance,
                max_order=args.max_order,
                duplicate_policy=args.duplicate_policy,
            )

            metadata_path = Path(f"{threshold_prefix}_metadata.tsv")
            distance_path = Path(f"{threshold_prefix}_distances.tsv")
            summary_path = Path(f"{threshold_prefix}_summary.tsv")
            write_threshold_metadata(
                metadata_path,
                input_path=args.input,
                state_path=state_path,
                blacklist_path=args.blacklist_bed,
                blacklisted_peaks=peak_summary.blacklisted_lines,
                selection=selection,
                results=results,
                score_column=args.score_column,
                min_distance=args.min_distance,
                max_distance=args.max_distance,
                max_order=args.max_order,
                duplicate_policy=args.duplicate_policy,
                bin_tie_mode=args.bin_tie_mode,
                nrl_mode=args.nrl_mode,
            )
            distributions, rows = write_distribution_outputs(
                results,
                distance_path=distance_path,
                summary_path=summary_path,
                max_order=args.max_order,
                include_chromosomes=include_chromosomes,
                include_genome=include_genome,
                include_state_strata=include_state_strata,
                include_zero_distances=args.include_zero_distances,
                count_smooth_window=args.count_smooth_window,
                count_smooth_polyorder=args.count_smooth_polyorder,
                percent_smooth_window=args.percent_smooth_window,
                percent_smooth_polyorder=args.percent_smooth_polyorder,
            )

            regression_outputs: list[Path] = []
            if args.max_order > 1:
                regressions = collect_nrl_regressions(
                    results,
                    max_order=args.max_order,
                    include_chromosomes=include_chromosomes,
                    include_genome=include_genome,
                    nrl_mode=args.nrl_mode,
                    count_smooth_window=args.count_smooth_window,
                    count_smooth_polyorder=args.count_smooth_polyorder,
                )
                regression_outputs = write_nrl_regression_outputs(
                    regressions,
                    prefix=threshold_prefix,
                )

            from nucleosuite.plotting import plot_path
            distribution_plot = plot_path(Path(f"{threshold_prefix}_distance_distribution.png"))
            plotted_distribution = plot_distance_distributions(
                distance_path, distribution_plot,
                nrl_mode=args.nrl_mode,
                label_peaks=True,
            )
            if plotted_distribution is not None:
                distribution_plot = plotted_distribution
            output_count = 3 + len(regression_outputs)
            if distribution_plot.exists():
                output_count += 1

            if percentile_store is not None:
                percentile_store.add_threshold(
                    selection=selection,
                    results=results,
                    chromosomes=input_chromosomes,
                    max_order=args.max_order,
                    include_chromosomes=include_chromosomes,
                    include_genome=include_genome,
                )

            if args.state_overlay_plot:
                assert state_indexes is not None
                state_counters = compute_same_interval_state_counts(
                    results.retained_by_chrom,
                    state_indexes,
                    category_rules=category_rules,
                    min_distance=args.min_distance,
                    max_distance=args.max_distance,
                )
                overlay_outputs = write_state_overlay_outputs(
                    state_counters,
                    state_colour_map(state_indexes, category_rules=category_rules),
                    prefix=threshold_prefix,
                    minimum=args.min_distance,
                    maximum=args.max_distance,
                    smooth_window=args.state_overlay_smooth_window,
                    smooth_polyorder=args.state_overlay_smooth_polyorder,
                    title=args.state_overlay_title,
                    plot_format=args.state_overlay_format,
                    x_major_tick=args.state_overlay_x_major_tick,
                    x_minor_tick=args.state_overlay_x_minor_tick,
                )
                output_count += len(overlay_outputs)

            duplicate_path = Path(f"{threshold_prefix}_duplicates.tsv")
            duplicate_positions = write_duplicate_report(results.duplicates, duplicate_path)
            if duplicate_positions:
                output_count += 1

            if args.write_filtered_bed:
                filtered_path = Path(f"{threshold_prefix}_filtered.bed")
                filtered_count = write_filtered_bed(
                    results.retained_by_chrom,
                    filtered_path,
                    args.interval_chrom_sizes,
                )
                interval_outputs = finalise_interval_files(
                    [filtered_path],
                    args.interval_format,
                    args.interval_chrom_sizes or {},
                )
                output_count += len(interval_outputs)

            reporter.stage(
                f"Completed {selection_kind} {index}/{len(selections)}: retained "
                f"{results.retained_count:,} peaks "
                f"({results.threshold_pass_count:,} before duplicate handling); "
                f"{distributions:,} distributions, {rows:,} distance rows, "
                f"{output_count:,} output files"
            )

        if percentile_store is not None:
            reporter.stage("Creating percentile-sweep count, percentage, and peak-count figures")
            percentile_outputs = percentile_store.write_outputs()
            reporter.stage(
                f"Wrote {len(percentile_outputs):,} percentile-sweep plot/table outputs"
            )
            percentile_store.close()
            percentile_store = None

        reporter.stage("Finished distance analysis")

    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    finally:
        if percentile_store is not None:
            percentile_store.close()

    return 0



def run(args: argparse.Namespace) -> int:
    from nucleosuite.partitioned import run_partitioned_command
    from nucleosuite.output_naming import parameterized_prefix

    requested = args.output_prefix or default_output_prefix(args.input, args.state_bed)
    grouping = output_grouping_token(args)
    args.output_prefix = str(
        parameterized_prefix(
            requested,
            (
                ("distmin", args.min_distance),
                ("distmax", args.max_distance),
                ("orders", args.max_order),
                ("nrlmode", args.nrl_mode),
                ("countsg", f"{args.count_smooth_window}x{args.count_smooth_polyorder}"),
                ("pctsg", f"{args.percent_smooth_window}x{args.percent_smooth_polyorder}"),
                ("groups", grouping),
                ("ties", args.bin_tie_mode),
            ),
        )
    )
    percentage_bins = (
        args.bin_tie_mode == "split"
        and (args.pct_bin_size is not None or args.pct_bins)
    )
    if percentage_bins and int(getattr(args, "cores", 1) or 1) > 1:
        print(
            "[distances] percentage binning requires one global ranking across all "
            "contigs; running the analysis serially."
        )
        return _run_serial(args)
    base = Path(args.output_prefix).name if args.output_prefix else default_output_prefix(args.input, args.state_bed).name
    return run_partitioned_command(
        "distances", args, _run_serial,
        runner_module="nucleosuite.distances", runner_function="_run_serial",
        primary_attr="input", output_prefix_attr="output_prefix",
        path_attrs=("state_bed",), base_name=base,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Standalone entry point, also useful for testing before package installation."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
