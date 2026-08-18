#!/usr/bin/env python3
"""Measure peak abundance and score-dependent enrichment by chromatin state."""

from __future__ import annotations

import argparse
import bisect
import math
import re
import sys
from array import array
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.core.regions import expand_contig_tokens, resolve_contig_name
from nucleosuite.io import open_text
from nucleosuite.percentiles import equal_rank_bins, randomized_score_order, rank_bins_from_boundaries
from nucleosuite.progress import ProgressReporter


HEADER_PREFIXES = ("#", "track", "browser")


@dataclass(frozen=True)
class StateInterval:
    start: int
    end: int
    label: str
    source_order: int


@dataclass
class StateIndex:
    intervals: list[StateInterval]
    starts: list[int]
    prefix_max_ends: list[int]


@dataclass
class PeakLoadSummary:
    data_lines: int = 0
    retained_peaks: int = 0
    assigned_peaks: int = 0
    unassigned_peaks: int = 0
    multi_state_peaks: int = 0
    blacklisted_peaks: int = 0
    filtered_contig_peaks: int = 0
    invalid_lines: int = 0


@dataclass(frozen=True)
class ThresholdSelection:
    percentile: float
    score_threshold: float
    percentile_lower: float | None = None
    percentile_upper: float | None = None
    score_upper_bound: float | None = None
    rank_start: int | None = None
    rank_stop: int | None = None
    rank_seed: int | None = None
    rank_score_max: float | None = None

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
                f"{_compact_number(self.percentile_lower)}-"
                f"{_compact_number(self.percentile_upper)}"
            )
        return _compact_number(self.percentile)


@dataclass(frozen=True)
class ThresholdSnapshot:
    counts: np.ndarray
    retained_peaks: int
    assigned_peaks: int
    unassigned_peaks: int


def _split_fields(raw: str) -> list[str]:
    fields = raw.rstrip("\n").split("\t")
    return fields if len(fields) > 1 else raw.split()


def _validate_column(value: int, option: str) -> int:
    if value < 1:
        raise ValueError(f"{option} must be a one-based column number of at least 1")
    return value - 1


def _natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(token) if token.isdigit() else token.casefold()
        for token in re.split(r"(\d+)", value)
    )


def _safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return token or "unnamed"


def _strip_interval_suffix(path: str | Path) -> str:
    name = Path(path).name
    for suffix in (".bed.gz", ".bigbed", ".bb", ".bed", ".tsv.gz", ".tsv", ".gz"):
        if name.casefold().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _normalise_rgb(value: str) -> str:
    try:
        channels = [int(token.strip()) for token in value.split(",")]
    except ValueError:
        channels = []
    if len(channels) != 3 or any(channel < 0 or channel > 255 for channel in channels):
        return "128,128,128"
    return ",".join(str(channel) for channel in channels)


def _rgb_tuple(value: str) -> tuple[float, float, float]:
    return tuple(int(token) / 255.0 for token in _normalise_rgb(value).split(","))  # type: ignore[return-value]


def read_state_intervals(
    path: str | Path,
    *,
    label_column: int,
    color_column: int,
    strict: bool,
    progress: ProgressReporter | None = None,
) -> tuple[
    dict[str, StateIndex],
    list[str],
    dict[str, str],
    dict[str, dict[str, list[tuple[int, int]]]],
    int,
    int,
]:
    """Read state intervals and retain first-seen state and colour ordering."""
    label_index = _validate_column(label_column, "--state-label-column")
    color_index = _validate_column(color_column, "--state-color-column")
    by_contig: dict[str, list[StateInterval]] = defaultdict(list)
    coverage: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    state_order: list[str] = []
    colors: dict[str, str] = {}
    color_conflicts = 0
    invalid = 0
    seen_contigs: set[str] = set()
    source_order = 0

    if progress is not None:
        progress.file_start("chromatin states", path)
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            text = raw.strip()
            if not text or text.startswith(HEADER_PREFIXES):
                continue
            fields = _split_fields(raw)
            error: str | None = None
            if len(fields) <= max(2, label_index):
                error = f"expected state label column {label_column}"
            else:
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    error = "state start/end must be integers"
                else:
                    if start < 0 or end <= start:
                        error = f"invalid half-open interval {start}-{end}"
            if error is not None:
                invalid += 1
                if strict:
                    raise ValueError(f"{path}:{line_number}: {error}")
                continue

            chromosome = fields[0]
            label = fields[label_index]
            if not label:
                invalid += 1
                if strict:
                    raise ValueError(f"{path}:{line_number}: state label is empty")
                continue
            rgb = _normalise_rgb(
                fields[color_index] if color_index < len(fields) else "128,128,128"
            )
            if label not in colors:
                state_order.append(label)
                colors[label] = rgb
            elif colors[label] != rgb:
                color_conflicts += 1
            source_order += 1
            interval = StateInterval(start, end, label, source_order)
            by_contig[chromosome].append(interval)
            coverage[label][chromosome].append((start, end))
            if progress is not None and chromosome not in seen_contigs:
                seen_contigs.add(chromosome)
                progress.reading_contig("chromatin states", chromosome)

    indexes: dict[str, StateIndex] = {}
    for chromosome, intervals in by_contig.items():
        intervals.sort(
            key=lambda interval: (interval.start, interval.end, interval.source_order)
        )
        starts = [interval.start for interval in intervals]
        prefix_max_ends: list[int] = []
        running_end = -1
        for interval in intervals:
            running_end = max(running_end, interval.end)
            prefix_max_ends.append(running_end)
        indexes[chromosome] = StateIndex(intervals, starts, prefix_max_ends)
    if not indexes:
        raise ValueError(f"No valid chromatin-state intervals were found in: {path}")
    return indexes, state_order, colors, coverage, invalid, color_conflicts


def states_at_position(index: StateIndex, position: int) -> list[str]:
    """Return unique state labels overlapping one zero-based genomic position."""
    right = bisect.bisect_right(index.starts, position) - 1
    if right < 0:
        return []
    matched: list[StateInterval] = []
    cursor = right
    while cursor >= 0 and index.prefix_max_ends[cursor] > position:
        interval = index.intervals[cursor]
        if interval.start <= position < interval.end:
            matched.append(interval)
        cursor -= 1
    matched.sort(key=lambda interval: interval.source_order)
    return list(dict.fromkeys(interval.label for interval in matched))


def _merge_intervals(rows: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(rows):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def _subtract_blacklist(
    chromosome: str,
    start: int,
    end: int,
    blacklist: BlacklistIndex | None,
) -> list[tuple[int, int]]:
    if blacklist is None:
        return [(start, end)]
    output: list[tuple[int, int]] = []
    cursor = start
    for left, right in blacklist.overlapping_intervals(chromosome, start, end):
        left = max(start, left)
        right = min(end, right)
        if left > cursor:
            output.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < end:
        output.append((cursor, end))
    return output


def calculate_state_coverage(
    coverage_intervals: Mapping[str, Mapping[str, Sequence[tuple[int, int]]]],
    *,
    selected_contigs: set[str] | None,
    blacklist: BlacklistIndex | None,
) -> tuple[dict[str, int], int, int]:
    """Return per-state, summed-state, and unique annotated coverage in bases."""
    state_coverage: dict[str, int] = {}
    all_by_contig: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for state, by_contig in coverage_intervals.items():
        state_total = 0
        for chromosome, rows in by_contig.items():
            if selected_contigs is not None and chromosome not in selected_contigs:
                continue
            allowed: list[tuple[int, int]] = []
            for start, end in rows:
                allowed.extend(
                    _subtract_blacklist(chromosome, start, end, blacklist)
                )
            merged = _merge_intervals(allowed)
            state_total += sum(end - start for start, end in merged)
            all_by_contig[chromosome].extend(merged)
        state_coverage[state] = state_total
    unique_annotated = sum(
        end - start
        for rows in all_by_contig.values()
        for start, end in _merge_intervals(rows)
    )
    summed_state_coverage = sum(state_coverage.values())
    return state_coverage, summed_state_coverage, unique_annotated


def load_peak_assignments(
    path: str | Path,
    *,
    state_indexes: Mapping[str, StateIndex],
    state_codes: Mapping[str, int],
    position_column: int | None,
    score_column: int,
    selected_contigs: set[str] | None,
    overlap_policy: str,
    blacklist: BlacklistIndex | None,
    strict: bool,
    progress: ProgressReporter | None = None,
) -> tuple[array, array, PeakLoadSummary]:
    """Stream peaks into compact score and state-code buffers."""
    score_index = _validate_column(score_column, "--score-column")
    position_index = (
        _validate_column(position_column, "--position-column")
        if position_column is not None
        else None
    )
    required_index = max(2, score_index, position_index or 0)
    scores = array("d")
    codes = array("i")
    summary = PeakLoadSummary()
    seen_contigs: set[str] = set()
    contig_cache: dict[str, str | None] = {}

    if progress is not None:
        progress.file_start("peaks", path)
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            text = raw.strip()
            if not text or text.startswith(HEADER_PREFIXES):
                continue
            summary.data_lines += 1
            fields = _split_fields(raw)
            error: str | None = None
            if len(fields) <= required_index:
                error = f"expected at least {required_index + 1} columns"
            else:
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    error = "peak start/end must be integers"
                else:
                    if start < 0 or end <= start:
                        error = f"invalid half-open interval {start}-{end}"
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
                        numeric_position = float(fields[position_index])
                    except ValueError:
                        error = f"position column {position_column} is not numeric"
                    else:
                        if (
                            not math.isfinite(numeric_position)
                            or not numeric_position.is_integer()
                        ):
                            error = (
                                f"position column {position_column} must be a finite integer"
                            )
                        else:
                            position = int(numeric_position)
                        if error is None and (position < start or position >= end):
                            error = (
                                f"position column {position_column} must fall inside "
                                "the peak interval"
                            )
            if error is not None:
                summary.invalid_lines += 1
                if strict:
                    raise ValueError(f"{path}:{line_number}: {error}")
                continue

            chromosome = fields[0]
            if progress is not None and chromosome not in seen_contigs:
                seen_contigs.add(chromosome)
                progress.reading_contig("peaks", chromosome)
            if blacklist is not None and blacklist.overlaps(chromosome, start, end):
                summary.blacklisted_peaks += 1
                continue
            if chromosome not in contig_cache:
                if chromosome in state_indexes:
                    resolved: str | None = chromosome
                else:
                    try:
                        resolved = resolve_contig_name(
                            chromosome,
                            list(state_indexes),
                            source_label="chromatin-state BED",
                        )
                    except (KeyError, ValueError):
                        resolved = None
                contig_cache[chromosome] = resolved
            state_chromosome = contig_cache[chromosome]
            if selected_contigs is not None and state_chromosome not in selected_contigs:
                summary.filtered_contig_peaks += 1
                continue

            labels = (
                states_at_position(state_indexes[state_chromosome], position)
                if state_chromosome is not None
                else []
            )
            if len(labels) > 1:
                summary.multi_state_peaks += 1
                if overlap_policy == "error":
                    raise ValueError(
                        f"{path}:{line_number}: peak position overlaps multiple states: "
                        + ", ".join(labels)
                    )
            code = state_codes[labels[0]] if labels else -1
            scores.append(score)
            codes.append(code)
            summary.retained_peaks += 1
            if code < 0:
                summary.unassigned_peaks += 1
            else:
                summary.assigned_peaks += 1

    if not scores:
        raise ValueError(f"No valid peaks remained for peak-state analysis: {path}")
    return scores, codes, summary


def percentile_values(start: float, stop: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("--pct-step must be greater than 0")
    lower = min(100.0, max(0.0, start))
    upper = min(100.0, max(0.0, stop))
    if upper < lower:
        lower, upper = upper, lower
    tolerance = max(1e-12, abs(step) * 1e-10)
    values: list[float] = []
    value = lower
    while value <= upper + tolerance:
        values.append(min(value, upper))
        value += step
    if not values or abs(values[-1] - upper) > tolerance:
        values.append(upper)
    return list(dict.fromkeys(round(value, 12) for value in values))


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


def _compact_number(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6f}".rstrip("0").rstrip(".")


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
    """Create score thresholds or non-overlapping score bins."""
    if scores.size == 0:
        raise ValueError("Cannot choose thresholds from an empty score array")
    if bin_tie_mode not in {"split", "keep"}:
        raise ValueError("--bin-tie-mode must be 'split' or 'keep'")

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
            percentiles = explicit
        elif pct_range:
            percentiles = percentile_values(pct_lower, pct_upper, pct_step)
        else:
            if pct_bins:
                raise ValueError("--pct-bins requires --pct-values or --pct-range")
            if not 0 <= score_percentile <= 100:
                raise ValueError("--score-percentile must be between 0 and 100")
            percentiles = [float(score_percentile)]

        if not pct_bins:
            thresholds = np.percentile(scores, percentiles)
            return [
                ThresholdSelection(float(percentile), float(threshold))
                for percentile, threshold in zip(percentiles, np.atleast_1d(thresholds))
            ]
        bounds = percentile_bin_bounds(percentiles)

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
                percentile=rank_bin.percentile_upper,
                score_threshold=float(ordered_scores[rank_bin.rank_start]),
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
            percentile=upper,
            score_threshold=boundary_scores[lower],
            percentile_lower=lower,
            percentile_upper=upper,
            score_upper_bound=(
                None if math.isclose(upper, 100.0) else boundary_scores[upper]
            ),
        )
        for lower, upper in bounds
    ]


def count_threshold_states(
    scores: np.ndarray,
    state_codes: np.ndarray,
    selections: Sequence[ThresholdSelection],
    *,
    state_count: int,
    progress: ProgressReporter | None = None,
) -> list[ThresholdSnapshot]:
    """Count thresholds in one descending score pass, including complete ties."""
    if selections and any(selection.is_rank_bin for selection in selections):
        if not all(selection.is_rank_bin for selection in selections):
            raise ValueError("Cannot mix rank bins with other percentile selections")
        return count_rank_bin_states(
            scores,
            state_codes,
            selections,
            state_count=state_count,
            progress=progress,
        )
    if selections and any(selection.is_bin for selection in selections):
        if not all(selection.is_bin for selection in selections):
            raise ValueError("Cannot mix percentile thresholds and percentile bins")
        return count_bin_states(
            scores,
            state_codes,
            selections,
            state_count=state_count,
            progress=progress,
        )

    order = np.argsort(scores, kind="quicksort")[::-1]
    ranked = sorted(
        enumerate(selections),
        key=lambda item: (item[1].score_threshold, item[1].percentile),
        reverse=True,
    )
    counts = np.zeros(state_count, dtype=np.int64)
    unassigned = 0
    cursor = 0
    snapshots: list[ThresholdSnapshot | None] = [None] * len(selections)
    for completed, (selection_index, selection) in enumerate(ranked, start=1):
        while cursor < order.size and scores[order[cursor]] >= selection.score_threshold:
            code = int(state_codes[order[cursor]])
            if code < 0:
                unassigned += 1
            else:
                counts[code] += 1
            cursor += 1
        retained = int(cursor)
        assigned = retained - unassigned
        snapshots[selection_index] = ThresholdSnapshot(
            counts=counts.copy(),
            retained_peaks=retained,
            assigned_peaks=assigned,
            unassigned_peaks=unassigned,
        )
        if progress is not None:
            progress.stage(
                f"Counted threshold {completed}/{len(selections)}: "
                f"percentile {selection.percentile:g}, retained {retained:,} peaks"
            )
    return [snapshot for snapshot in snapshots if snapshot is not None]


def count_bin_states(
    scores: np.ndarray,
    state_codes: np.ndarray,
    selections: Sequence[ThresholdSelection],
    *,
    state_count: int,
    progress: ProgressReporter | None = None,
) -> list[ThresholdSnapshot]:
    """Count non-overlapping percentile bins in one ascending score ordering."""
    order = np.argsort(scores, kind="quicksort")
    ordered_scores = scores[order]
    snapshots: list[ThresholdSnapshot] = []
    for completed, selection in enumerate(selections, start=1):
        start = int(np.searchsorted(ordered_scores, selection.score_threshold, side="left"))
        end = (
            int(ordered_scores.size)
            if selection.score_upper_bound is None
            else int(
                np.searchsorted(
                    ordered_scores,
                    selection.score_upper_bound,
                    side="left",
                )
            )
        )
        selected_codes = state_codes[order[start:end]]
        assigned_codes = selected_codes[selected_codes >= 0]
        counts = np.bincount(
            assigned_codes.astype(np.int64, copy=False),
            minlength=state_count,
        )[:state_count]
        unassigned = int(np.count_nonzero(selected_codes < 0))
        retained = int(selected_codes.size)
        snapshots.append(
            ThresholdSnapshot(
                counts=counts.astype(np.int64, copy=False),
                retained_peaks=retained,
                assigned_peaks=retained - unassigned,
                unassigned_peaks=unassigned,
            )
        )
        if progress is not None:
            progress.stage(
                f"Counted bin {completed}/{len(selections)}: percentile "
                f"{selection.label}, retained {retained:,} peaks"
            )
    return snapshots


def count_rank_bin_states(
    scores: np.ndarray,
    state_codes: np.ndarray,
    selections: Sequence[ThresholdSelection],
    *,
    state_count: int,
    progress: ProgressReporter | None = None,
) -> list[ThresholdSnapshot]:
    """Count exact global rank slices after randomized tie ordering."""
    seed = selections[0].rank_seed
    assert seed is not None
    order = randomized_score_order(scores, seed)
    snapshots: list[ThresholdSnapshot] = []
    for completed, selection in enumerate(selections, start=1):
        assert selection.rank_start is not None and selection.rank_stop is not None
        selected_codes = state_codes[order[selection.rank_start : selection.rank_stop]]
        assigned_codes = selected_codes[selected_codes >= 0]
        counts = np.bincount(
            assigned_codes.astype(np.int64, copy=False),
            minlength=state_count,
        )[:state_count]
        unassigned = int(np.count_nonzero(selected_codes < 0))
        retained = int(selected_codes.size)
        snapshots.append(
            ThresholdSnapshot(
                counts=counts.astype(np.int64, copy=False),
                retained_peaks=retained,
                assigned_peaks=retained - unassigned,
                unassigned_peaks=unassigned,
            )
        )
        if progress is not None:
            progress.stage(
                f"Counted equal-rank bin {completed}/{len(selections)}: "
                f"percentile {selection.label}, retained {retained:,} peaks"
            )
    return snapshots


def _format_float(value: float) -> str:
    return f"{value:.9f}" if math.isfinite(value) else "NA"


def _selection_output_fields(selection: ThresholdSelection) -> tuple[str, ...]:
    if selection.is_bin:
        return (
            "rank_bin" if selection.is_rank_bin else "bin",
            selection.label,
            _compact_number(selection.percentile_lower),
            _compact_number(selection.percentile_upper),
            (
                ""
                if selection.score_upper_bound is None
                else f"{selection.score_upper_bound:.12g}"
            ),
            "" if selection.rank_start is None else str(selection.rank_start),
            "" if selection.rank_stop is None else str(selection.rank_stop),
            "" if selection.rank_seed is None else str(selection.rank_seed),
            (
                ""
                if selection.rank_score_max is None
                else f"{selection.rank_score_max:.12g}"
            ),
        )
    return ("threshold", "", "", "", "", "", "", "", "")


def write_outputs(
    prefix: str | Path,
    *,
    peak_path: str | Path,
    state_path: str | Path,
    blacklist_path: str | Path | None,
    selections: Sequence[ThresholdSelection],
    snapshots: Sequence[ThresholdSnapshot],
    state_order: Sequence[str],
    colors: Mapping[str, str],
    state_coverage: Mapping[str, int],
    summed_state_coverage: int,
    unique_annotated_coverage: int,
    peak_summary: PeakLoadSummary,
    invalid_state_lines: int,
    color_conflicts: int,
    position_column: int | None,
    score_column: int,
    bin_tie_mode: str = "split",
) -> dict[str, Path]:
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    coverage_path = Path(f"{prefix}_state_coverage.tsv")
    enrichment_path = Path(f"{prefix}_peak_state_enrichment.tsv")
    threshold_path = Path(f"{prefix}_peak_state_threshold_summary.tsv")
    metadata_path = Path(f"{prefix}_peak_states_metadata.tsv")

    with coverage_path.open("wt", encoding="utf-8") as handle:
        handle.write(
            "state\trgb\tcoverage_bp\tcoverage_mb\t"
            "coverage_pct_of_state_annotations\tcoverage_pct_of_unique_annotated_bases\n"
        )
        for state in state_order:
            bases = int(state_coverage.get(state, 0))
            annotation_pct = (
                100.0 * bases / summed_state_coverage
                if summed_state_coverage > 0
                else math.nan
            )
            unique_pct = (
                100.0 * bases / unique_annotated_coverage
                if unique_annotated_coverage > 0
                else math.nan
            )
            handle.write(
                f"{state}\t{colors[state]}\t{bases}\t{bases / 1_000_000:.9f}\t"
                f"{_format_float(annotation_pct)}\t{_format_float(unique_pct)}\n"
            )

    with enrichment_path.open("wt", encoding="utf-8") as handle:
        handle.write(
            "percentile_threshold\tscore_threshold\tstate\trgb\tpeak_count\t"
            "percentage_of_assigned_peaks\tpercentage_of_all_retained_peaks\t"
            "state_coverage_bp\tstate_coverage_pct\tpeaks_per_mb\t"
            "enrichment_vs_state_coverage\tpercentile_mode\tpercentile_range\t"
            "percentile_lower\tpercentile_upper\tscore_upper_bound\t"
            "rank_start\trank_stop\trank_seed\trank_score_max\n"
        )
        for selection, snapshot in zip(selections, snapshots):
            selection_fields = "\t".join(_selection_output_fields(selection))
            for code, state in enumerate(state_order):
                count = int(snapshot.counts[code])
                bases = int(state_coverage.get(state, 0))
                peak_fraction = (
                    count / snapshot.assigned_peaks
                    if snapshot.assigned_peaks > 0
                    else math.nan
                )
                all_fraction = (
                    count / snapshot.retained_peaks
                    if snapshot.retained_peaks > 0
                    else math.nan
                )
                coverage_fraction = (
                    bases / summed_state_coverage
                    if summed_state_coverage > 0
                    else math.nan
                )
                peaks_per_mb = count * 1_000_000 / bases if bases > 0 else math.nan
                enrichment = (
                    peak_fraction / coverage_fraction
                    if math.isfinite(peak_fraction)
                    and math.isfinite(coverage_fraction)
                    and coverage_fraction > 0
                    else math.nan
                )
                handle.write(
                    f"{selection.percentile:.12g}\t{selection.score_threshold:.12g}\t"
                    f"{state}\t{colors[state]}\t{count}\t"
                    f"{_format_float(100.0 * peak_fraction)}\t"
                    f"{_format_float(100.0 * all_fraction)}\t{bases}\t"
                    f"{_format_float(100.0 * coverage_fraction)}\t"
                    f"{_format_float(peaks_per_mb)}\t{_format_float(enrichment)}\t"
                    f"{selection_fields}\n"
                )

    with threshold_path.open("wt", encoding="utf-8") as handle:
        handle.write(
            "percentile_threshold\tscore_threshold\tretained_peak_count\t"
            "assigned_peak_count\tunassigned_peak_count\tassigned_percent\t"
            "unassigned_percent\tpercentile_mode\tpercentile_range\t"
            "percentile_lower\tpercentile_upper\tscore_upper_bound\t"
            "rank_start\trank_stop\trank_seed\trank_score_max\n"
        )
        for selection, snapshot in zip(selections, snapshots):
            selection_fields = "\t".join(_selection_output_fields(selection))
            assigned_pct = (
                100.0 * snapshot.assigned_peaks / snapshot.retained_peaks
                if snapshot.retained_peaks > 0
                else math.nan
            )
            unassigned_pct = (
                100.0 * snapshot.unassigned_peaks / snapshot.retained_peaks
                if snapshot.retained_peaks > 0
                else math.nan
            )
            handle.write(
                f"{selection.percentile:.12g}\t{selection.score_threshold:.12g}\t"
                f"{snapshot.retained_peaks}\t{snapshot.assigned_peaks}\t"
                f"{snapshot.unassigned_peaks}\t{_format_float(assigned_pct)}\t"
                f"{_format_float(unassigned_pct)}\t"
                f"{selection_fields}\n"
            )

    metadata = (
        ("peak_bed", peak_path),
        ("state_bed", state_path),
        ("blacklist_bed", blacklist_path or ""),
        ("position_column", position_column or "midpoint"),
        ("score_column", score_column),
        ("bin_tie_mode", bin_tie_mode if selections and any(selection.is_bin for selection in selections) else ""),
        (
            "percentile_mode",
            (
                "equal_rank_bins"
                if selections and all(selection.is_rank_bin for selection in selections)
                else "score_bins"
                if selections and all(selection.is_bin for selection in selections)
                else "thresholds"
            ),
        ),
        (
            "pct_bin_size",
            (
                ""
                if not selections or not selections[0].is_rank_bin
                else _compact_number(
                    selections[0].percentile_upper
                    - selections[0].percentile_lower
                )
            ),
        ),
        (
            "pct_bin_seed",
            (
                ""
                if not selections or selections[0].rank_seed is None
                else selections[0].rank_seed
            ),
        ),
        (
            "rank_tie_order",
            (
                "global_random_shuffle_then_stable_ascending_score_sort"
                if selections and selections[0].is_rank_bin
                else ""
            ),
        ),
        (
            "rank_membership_rule",
            (
                "exact_zero_based_start_stop_slice"
                if selections and selections[0].is_rank_bin
                else ""
            ),
        ),
        ("input_data_lines", peak_summary.data_lines),
        ("valid_peaks", peak_summary.retained_peaks),
        ("assigned_peaks_before_score_filter", peak_summary.assigned_peaks),
        ("unassigned_peaks_before_score_filter", peak_summary.unassigned_peaks),
        ("multi_state_peaks", peak_summary.multi_state_peaks),
        ("blacklisted_peaks_excluded", peak_summary.blacklisted_peaks),
        ("contig_filtered_peaks", peak_summary.filtered_contig_peaks),
        ("invalid_peak_lines", peak_summary.invalid_lines),
        ("invalid_state_lines", invalid_state_lines),
        ("state_color_conflicts", color_conflicts),
        ("summed_state_coverage_bp", summed_state_coverage),
        ("unique_annotated_coverage_bp", unique_annotated_coverage),
        (
            "overlapping_state_coverage_bp",
            max(0, summed_state_coverage - unique_annotated_coverage),
        ),
        ("enrichment_denominator", "state_coverage_bp/summed_state_coverage_bp"),
    )
    with metadata_path.open("wt", encoding="utf-8") as handle:
        handle.write("parameter\tvalue\n")
        for name, value in metadata:
            handle.write(f"{name}\t{value}\n")

    return {
        "state_coverage": coverage_path,
        "state_enrichment": enrichment_path,
        "threshold_summary": threshold_path,
        "metadata": metadata_path,
    }


def plot_state_percentages(
    path: str | Path,
    *,
    selections: Sequence[ThresholdSelection],
    snapshots: Sequence[ThresholdSnapshot],
    state_order: Sequence[str],
    colors: Mapping[str, str],
    title: str | None,
    dpi: int,
    x_axis_mode: str = "categorical",
    bar_gap: float = 0.18,
) -> None:
    """Plot assigned-peak percentages as score-percentile stacked bars."""
    if not math.isfinite(bar_gap) or not 0.0 <= bar_gap < 1.0:
        raise ValueError("bar_gap must be at least 0 and less than 1")
    bar_fraction = 1.0 - float(bar_gap)
    labels = [selection.label for selection in selections]
    bins = bool(selections) and all(selection.is_bin for selection in selections)
    if x_axis_mode == "categorical":
        x = np.arange(len(selections), dtype=float)
        width: float | np.ndarray = bar_fraction
        if len(x) <= 15:
            tick_indices = np.arange(len(x), dtype=int)
        else:
            tick_indices = np.unique(np.linspace(0, len(x) - 1, 11, dtype=int))
        tick_positions = x[tick_indices]
        tick_labels = [labels[index] for index in tick_indices]
    elif x_axis_mode == "continuous":
        if bins:
            lower = np.asarray(
                [float(selection.percentile_lower) for selection in selections]
            )
            upper = np.asarray(
                [float(selection.percentile_upper) for selection in selections]
            )
            x = (lower + upper) / 2.0
            width = bar_fraction * (upper - lower)
            boundaries = np.concatenate((lower[:1], upper))
            if len(boundaries) <= 15:
                tick_positions = boundaries
                tick_labels = [_compact_number(value) for value in boundaries]
            else:
                tick_positions = np.asarray([], dtype=float)
                tick_labels = []
        else:
            x = np.asarray(
                [selection.percentile for selection in selections],
                dtype=float,
            )
            if len(x) > 1:
                positive_steps = np.diff(np.unique(x))
                positive_steps = positive_steps[positive_steps > 0]
                width = (
                    bar_fraction * float(positive_steps.min())
                    if positive_steps.size
                    else bar_fraction
                )
            else:
                width = bar_fraction
            if len(x) <= 15:
                tick_positions = x
                tick_labels = labels
            else:
                tick_positions = np.asarray([], dtype=float)
                tick_labels = []
    else:
        raise ValueError("x_axis_mode must be 'categorical' or 'continuous'")

    figure, axis = plt.subplots(figsize=(max(9.5, min(16.0, 7.5 + 0.05 * len(x))), 6.8))
    bottom = np.zeros(len(selections), dtype=float)
    for code, state in enumerate(state_order):
        values = np.asarray(
            [
                100.0 * int(snapshot.counts[code]) / snapshot.assigned_peaks
                if snapshot.assigned_peaks > 0
                else 0.0
                for snapshot in snapshots
            ],
            dtype=float,
        )
        axis.bar(
            x,
            values,
            width=width,
            bottom=bottom,
            color=_rgb_tuple(colors[state]),
            edgecolor="none",
            label=state,
        )
        bottom += values
    axis.set_ylim(0, 100)
    axis.set_xlabel(
        "Peak score percentile range" if bins else "Peak score threshold percentile"
    )
    axis.set_ylabel("Assigned peaks in each chromatin state (%)")
    axis.set_title(
        title
        or (
            "Chromatin-state composition across peak-score percentile ranges"
            if bins
            else "Chromatin-state composition across peak-score thresholds"
        )
    )
    if tick_positions.size:
        axis.set_xticks(tick_positions)
        axis.set_xticklabels(tick_labels)
    axis.grid(axis="y", alpha=0.3)
    axis.legend(
        title="Chromatin state",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        frameon=False,
        fontsize=8,
    )
    figure.tight_layout()
    from nucleosuite.plotting import save_figure
    save_figure(figure, path, default_dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def default_output_prefix(peak_path: str | Path, state_path: str | Path) -> Path:
    return Path(
        f"{_strip_interval_suffix(peak_path)}_{_strip_interval_suffix(state_path)}_peak_states"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite peak-states",
        description=(
            "Count peaks by chromatin state and compare peak abundance with each "
            "state's annotated genomic coverage."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("peaks", help="Peak BED, BED.gz, bigBed, or .bb input.")
    parser.add_argument("--state-bed", required=True, help="Chromatin-state BED, BED.gz, bigBed, or .bb input.")
    parser.add_argument("--blacklist-bed", help="BED blacklist; overlapping complete peaks and state-coverage bases are excluded.")
    parser.add_argument("--position-column", type=int, help="One-based exact peak-position column; default uses the BED interval midpoint.")
    parser.add_argument("--score-column", type=int, default=5, help="One-based numeric peak-score column (default: 5).")
    parser.add_argument("--state-label-column", type=int, default=4, help="One-based chromatin-state label column (default: 4).")
    parser.add_argument("--state-color-column", type=int, default=9, help="One-based R,G,B state colour column (default: 9).")
    parser.add_argument("--contigs", nargs="+", help="Optional state contigs to analyse; supports comma lists and numeric ranges.")
    parser.add_argument("--overlap-policy", choices=("first", "error"), default="first", help="Resolve a peak point overlapping multiple state labels using the first BED state, or stop (default: first).")
    parser.add_argument("--score-percentile", type=float, default=0.0, help="Retain scores at or above this global percentile (default: 0).")
    parser.add_argument("--pct-range", action="store_true", help="Analyse an inclusive global score-percentile sweep.")
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
    parser.add_argument("-o", "--output-prefix", help="Output path prefix; default combines the peak and state BED basenames.")
    parser.add_argument("--plot-title", help="Optional stacked-bar plot title.")
    parser.add_argument(
        "--plot-x-axis",
        choices=("categorical", "continuous"),
        default="categorical",
        help=(
            "Stacked-plot x-axis spacing: equally spaced categories or numeric "
            "percentile positions (default: categorical)."
        ),
    )
    parser.add_argument(
        "--plot-bar-gap",
        type=float,
        default=0.18,
        metavar="FRACTION",
        help=(
            "Fraction of each stacked-bar slot left as white space between bars; "
            "0 removes the gaps entirely (default: 0.18)."
        ),
    )
    parser.add_argument("--dpi", type=int, default=220, help="Stacked-bar plot resolution (default: 220).")
    parser.add_argument("--strict", action="store_true", help="Stop at the first malformed peak or state record.")
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def run(args: argparse.Namespace) -> int:
    if args.dpi < 50:
        raise ValueError("--dpi must be at least 50")
    if not math.isfinite(args.plot_bar_gap) or not 0.0 <= args.plot_bar_gap < 1.0:
        raise ValueError("--plot-bar-gap must be at least 0 and less than 1")
    reporter = ProgressReporter("peak-states")
    blacklist = load_blacklist_unbounded(args.blacklist_bed)
    (
        state_indexes,
        state_order,
        colors,
        coverage_intervals,
        invalid_state_lines,
        color_conflicts,
    ) = read_state_intervals(
        args.state_bed,
        label_column=args.state_label_column,
        color_column=args.state_color_column,
        strict=args.strict,
        progress=reporter,
    )
    selected_contigs = (
        set(expand_contig_tokens(args.contigs, list(state_indexes)))
        if args.contigs
        else None
    )
    if selected_contigs is not None:
        state_indexes = {
            chromosome: index
            for chromosome, index in state_indexes.items()
            if chromosome in selected_contigs
        }
    reporter.stage("Calculating non-blacklisted chromatin-state coverage")
    state_coverage, summed_state_coverage, unique_annotated_coverage = (
        calculate_state_coverage(
            coverage_intervals,
            selected_contigs=selected_contigs,
            blacklist=blacklist,
        )
    )
    if summed_state_coverage <= 0:
        raise ValueError("Chromatin states have zero usable coverage")
    state_order = sorted(
        (state for state in state_order if state_coverage.get(state, 0) > 0),
        key=_natural_key,
    )
    reporter.stage(
        f"Loaded {len(state_order):,} states across {len(state_indexes):,} contigs"
    )
    state_codes = {state: code for code, state in enumerate(state_order)}
    scores_buffer, codes_buffer, peak_summary = load_peak_assignments(
        args.peaks,
        state_indexes=state_indexes,
        state_codes=state_codes,
        position_column=args.position_column,
        score_column=args.score_column,
        selected_contigs=selected_contigs,
        overlap_policy=args.overlap_policy,
        blacklist=blacklist,
        strict=args.strict,
        progress=reporter,
    )
    reporter.stage(
        f"Loaded {peak_summary.retained_peaks:,} peaks: "
        f"{peak_summary.assigned_peaks:,} assigned, "
        f"{peak_summary.unassigned_peaks:,} unassigned"
    )

    scores = np.frombuffer(scores_buffer, dtype=np.float64)
    state_codes_array = np.frombuffer(codes_buffer, dtype=np.int32)
    if (args.pct_bins or args.pct_bin_size is not None) and args.bin_tie_mode == "split":
        reporter.stage(
            f"Randomizing global peak order with seed {args.pct_bin_seed}, "
            "stably sorting by score, and creating percentage rank groups"
        )
    selections = choose_thresholds(
        scores,
        score_percentile=args.score_percentile,
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
    selection_kind = (
        f"{args.bin_tie_mode} bin(s)"
        if args.pct_bins or args.pct_bin_size is not None
        else "threshold(s)"
    )
    reporter.stage(
        f"Ranking scores and counting {len(selections):,} percentile {selection_kind}"
    )
    snapshots = count_threshold_states(
        scores,
        state_codes_array,
        selections,
        state_count=len(state_order),
        progress=reporter,
    )

    prefix = Path(args.output_prefix) if args.output_prefix else default_output_prefix(args.peaks, args.state_bed)
    reporter.stage("Writing state coverage, enrichment, and threshold tables")
    outputs = write_outputs(
        prefix,
        peak_path=args.peaks,
        state_path=args.state_bed,
        blacklist_path=args.blacklist_bed,
        selections=selections,
        snapshots=snapshots,
        state_order=state_order,
        colors=colors,
        state_coverage=state_coverage,
        summed_state_coverage=summed_state_coverage,
        unique_annotated_coverage=unique_annotated_coverage,
        peak_summary=peak_summary,
        invalid_state_lines=invalid_state_lines,
        color_conflicts=color_conflicts,
        position_column=args.position_column,
        score_column=args.score_column,
        bin_tie_mode=args.bin_tie_mode,
    )
    from nucleosuite.plotting import plot_path as resolve_plot_path
    plot_path = resolve_plot_path(Path(f"{prefix}_peak_state_percentages_stacked.png"))
    reporter.stage("Creating chromatin-state percentage stacked-bar plot")
    plot_state_percentages(
        plot_path,
        selections=selections,
        snapshots=snapshots,
        state_order=state_order,
        colors=colors,
        title=args.plot_title,
        dpi=args.dpi,
        x_axis_mode=args.plot_x_axis,
        bar_gap=args.plot_bar_gap,
    )
    outputs["stacked_percentage_plot"] = plot_path
    reporter.stage(
        f"Wrote {len(outputs):,} peak-state outputs; enrichment uses each "
        "state's fraction of summed annotated state coverage"
    )
    for name, path in outputs.items():
        print(f"{name}\t{path}")
    return 0


def validate_argv(argv: Sequence[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    try:
        return run(parser.parse_args(argv))
    except (FileNotFoundError, OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
