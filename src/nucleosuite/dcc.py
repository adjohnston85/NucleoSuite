#!/usr/bin/env python3
"""
Calculate Distance Cross-Correlation (DCC) from either bigWig signals or
BAM- or fragment-interval-derived position signals.

DCC is calculated within each selected region as:

    DCC[lag] = sum_x A[x] * B[x + lag]

where lag = position_B - position_A. Positive signed lags place signal B to the
right/downstream of signal A. For minus-strand BED regions, signals are reversed
before calculation so positive lag remains downstream in feature-oriented
coordinates.

By default, signed lags are collapsed into absolute distances:

    Distance 0 = lag 0
    Distance d = lag -d + lag +d

The script has two input modes:

    dcc.py bigwig ...
    dcc.py bam ...

Region sources:
  * --regions-bed: BED intervals or strand-aware upstream/downstream windows.
  * --chrom-sizes: generated whole-genome or chromosome windows.

BED state names are preserved by default. Repeatable --category rules can group
states using exact, prefix or regular-expression matching. No organism-specific
chromosome or chromatin-state scheme is built in.

Dependencies:
  * numpy
  * pyBigWig for bigWig mode
  * pysam for BAM mode
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

from nucleosuite.core.chrom_sizes import read_chrom_sizes_source
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.io import open_text
from nucleosuite.core.fragment_inputs import IntervalFragmentSource
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.progress import ProgressReporter

import numpy as np


# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Region:
    """A genomic interval assigned to one DCC state/group."""

    chrom: str
    start: int
    end: int
    state: str
    strand: str = "+"
    anchor_start: int | None = None
    anchor_end: int | None = None

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class BigWigTrack:
    """An open bigWig and its chromosome metadata."""

    path: str
    handle: object
    chrom_sizes: Mapping[str, int]
    chrom_cache: MutableMapping[str, Optional[str]]


@dataclass
class BamTrack:
    """An open BAM and its chromosome metadata."""

    path: str
    handle: object
    references: set[str]
    chrom_cache: MutableMapping[str, Optional[str]]


@dataclass(frozen=True)
class FragmentLengthSpec:
    """Inclusive fragment-length range for a BAM-derived signal."""

    min_length: int
    max_length: int

    def matches(self, length: int) -> bool:
        return self.min_length <= length <= self.max_length

    @property
    def label(self) -> str:
        if self.min_length == self.max_length:
            return str(self.min_length)
        return f"{self.min_length}-{self.max_length}"


@dataclass
class DccStats:
    """Processing statistics for one DCC state output."""

    regions_seen: int = 0
    regions_with_a: int = 0
    regions_with_b: int = 0
    regions_with_both: int = 0
    missing_chromosome_a: int = 0
    missing_chromosome_b: int = 0
    short_or_empty: int = 0
    clipped_regions_a: int = 0
    clipped_regions_b: int = 0
    signal_positions_a: int = 0
    signal_positions_b: int = 0
    nonzero_positions_a: int = 0
    nonzero_positions_b: int = 0
    total_signal_a: float = 0.0
    total_signal_b: float = 0.0
    blacklisted_positions: int = 0
    blacklisted_fragments_a: int = 0
    blacklisted_fragments_b: int = 0


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------


def sanitize_filename(value: str) -> str:
    """Return a filesystem-friendly name."""

    value = str(value).strip()
    value = re.sub(r"[\\/\s]+", "_", value)
    value = re.sub(r"[^A-Za-z0-9._+-]+", "_", value)
    value = value.strip("._")
    return value or "unnamed"


def split_comma_values(values: Optional[Sequence[str]]) -> Optional[set[str]]:
    """Expand repeated and comma-separated values."""

    if not values:
        return None

    expanded: set[str] = set()
    for value in values:
        expanded.update(item.strip() for item in value.split(",") if item.strip())
    return expanded or None


def expand_file_inputs(inputs: Sequence[str], kind: str) -> List[str]:
    """Expand glob patterns and return unique existing files."""

    files: List[str] = []
    seen: set[str] = set()

    for item in inputs:
        matches = sorted(glob.glob(item))
        candidates = matches if matches else [item]

        for candidate in candidates:
            path = os.path.abspath(os.path.expanduser(candidate))
            if os.path.isfile(path) and path not in seen:
                seen.add(path)
                files.append(path)

    if not files:
        raise FileNotFoundError(
            f"No {kind} files were found for: " + " ".join(map(str, inputs))
        )

    return files


def compact_input_label(paths: Sequence[str], default_label: str) -> str:
    """Construct a concise label from input paths."""

    if not paths:
        return default_label

    if len(paths) == 1:
        return sanitize_filename(Path(paths[0]).stem)

    stems = [Path(path).stem for path in paths]
    common = os.path.commonprefix(stems).rstrip("._-")
    return sanitize_filename(common or f"combined_{len(paths)}")


def candidate_chromosome_aliases(chrom: str) -> List[str]:
    """Return conservative aliases without assuming a particular genome."""

    aliases = [chrom]

    if chrom.startswith("chr") and len(chrom) > 3:
        aliases.append(chrom[3:])
    else:
        aliases.append(f"chr{chrom}")

    mitochondrial_aliases = {
        "M": ["MT", "chrM", "chrMT"],
        "MT": ["M", "chrM", "chrMT"],
        "chrM": ["M", "MT", "chrMT"],
        "chrMT": ["M", "MT", "chrM"],
    }
    aliases.extend(mitochondrial_aliases.get(chrom, []))

    unique: List[str] = []
    seen: set[str] = set()
    for alias in aliases:
        if alias not in seen:
            seen.add(alias)
            unique.append(alias)
    return unique


def read_chrom_sizes(path: str) -> List[Tuple[str, int]]:
    """Read chromosome sizes from a table, BAM header or CRAM header."""

    return list(read_chrom_sizes_source(path))


def make_windows_from_chrom_sizes(
    chrom_sizes_path: str,
    scope: str,
    selected_chromosomes: Optional[set[str]],
    window_size: int,
    state_name: str,
    min_region_length: int,
) -> List[Region]:
    """Generate non-overlapping analysis windows."""

    if window_size <= 0:
        raise ValueError("--window-size must be greater than zero.")
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    if selected_chromosomes is not None:
        available_names = [chrom for chrom, _size in chrom_sizes]
        resolved_selected: set[str] = set()
        for requested in selected_chromosomes:
            try:
                resolved_selected.add(
                    resolve_contig_name(
                        requested, available_names, source_label="chromosome sizes"
                    )
                )
            except KeyError as exc:
                raise ValueError(str(exc)) from exc
        selected_chromosomes = resolved_selected
    if scope == "chromosome" and not selected_chromosomes:
        raise ValueError(
            "At least one --chromosome is required with --scope chromosome."
        )

    regions: List[Region] = []
    found: set[str] = set()

    for chrom, size in chrom_sizes:
        if selected_chromosomes is not None and chrom not in selected_chromosomes:
            continue

        found.add(chrom)
        for start in range(0, size, window_size):
            end = min(start + window_size, size)
            if end - start >= min_region_length:
                regions.append(Region(chrom, start, end, state_name, "+"))

    if selected_chromosomes:
        missing = selected_chromosomes - found
        if missing:
            raise ValueError(
                "Chromosome(s) not found in chromosome sizes file: "
                + ", ".join(sorted(missing))
            )

    if not regions:
        raise ValueError("No generated windows remained after filtering.")

    return regions


StateMatcher = Tuple[str, object]
StateCategoryRule = Tuple[str, Tuple[StateMatcher, ...]]


def parse_state_category_rules(
    specifications: Optional[Sequence[str]],
) -> List[StateCategoryRule]:
    """Parse CATEGORY=MATCHER[,MATCHER...] state grouping rules."""

    rules: List[StateCategoryRule] = []

    for specification in specifications or []:
        if "=" not in specification:
            raise ValueError(
                f"Invalid --category rule {specification!r}; expected "
                "CATEGORY=MATCHER[,MATCHER...]."
            )

        category, matcher_text = specification.split("=", 1)
        category = category.strip()
        raw_matchers = [item.strip() for item in matcher_text.split(",") if item.strip()]

        if not category:
            raise ValueError(f"Invalid --category rule {specification!r}: empty category.")
        if not raw_matchers:
            raise ValueError(f"Invalid --category rule {specification!r}: no matchers.")

        parsed: List[StateMatcher] = []
        for matcher in raw_matchers:
            if matcher.startswith("exact:"):
                value = matcher[len("exact:") :]
                if not value:
                    raise ValueError(f"Empty exact matcher in {specification!r}.")
                parsed.append(("exact", value))
            elif matcher.startswith("prefix:"):
                value = matcher[len("prefix:") :]
                if not value:
                    raise ValueError(f"Empty prefix matcher in {specification!r}.")
                parsed.append(("prefix", value))
            elif matcher.startswith("regex:"):
                pattern = matcher[len("regex:") :]
                if not pattern:
                    raise ValueError(f"Empty regex matcher in {specification!r}.")
                try:
                    parsed.append(("regex", re.compile(pattern)))
                except re.error as exc:
                    raise ValueError(
                        f"Invalid regex in --category {specification!r}: {exc}"
                    ) from exc
            else:
                parsed.append(("exact", matcher))

        rules.append((category, tuple(parsed)))

    return rules


def categorize_state_name(state: str, rules: Sequence[StateCategoryRule]) -> str:
    """Return the first category matching a state name."""

    for category, matchers in rules:
        for match_type, matcher in matchers:
            if match_type == "exact" and state == matcher:
                return category
            if match_type == "prefix" and state.startswith(str(matcher)):
                return category
            if match_type == "regex" and matcher.search(state):
                return category
    return state


def read_regions_bed(
    path: str,
    state_column: Optional[int],
    strand_column: Optional[int],
    state_name: str,
    region_mode: str,
    extend: int,
    strands: str,
    selected_chromosomes: Optional[set[str]],
    min_region_length: int,
    category_rules: Sequence[StateCategoryRule],
) -> List[Region]:
    """Read BED regions and optionally construct strand-aware windows."""

    if state_column is not None and state_column < 1:
        raise ValueError("--state-column must be a 1-based column number.")
    if strand_column is not None and strand_column < 1:
        raise ValueError("--strand-column must be a 1-based column number.")
    if region_mode != "interval" and extend <= 0:
        raise ValueError("--extend must be greater than zero for anchored windows.")

    state_index = state_column - 1 if state_column is not None else None
    strand_index = strand_column - 1 if strand_column is not None else None
    regions: List[Region] = []

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}: line {line_number} has fewer than three BED columns."
                )

            chrom = fields[0]
            if selected_chromosomes is not None:
                try:
                    resolve_contig_name(
                        chrom,
                        list(selected_chromosomes),
                        source_label="selected chromosomes",
                    )
                except KeyError:
                    continue

            try:
                feature_start = int(fields[1])
                feature_end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}: line {line_number} has non-integer coordinates."
                ) from exc

            if feature_start < 0 or feature_end <= feature_start:
                continue

            state = state_name
            if state_index is not None and state_index < len(fields):
                state = fields[state_index]
            state = categorize_state_name(str(state), category_rules)

            strand = "+"
            if strand_index is not None and strand_index < len(fields):
                candidate = fields[strand_index]
                if candidate in {"+", "-"}:
                    strand = candidate
                elif region_mode != "interval":
                    continue
            elif region_mode != "interval":
                raise ValueError(
                    f"{path}: line {line_number} lacks requested strand column "
                    f"{strand_column}."
                )

            if strands == "plus" and strand != "+":
                continue
            if strands == "minus" and strand != "-":
                continue

            if region_mode == "interval":
                start, end = feature_start, feature_end
            elif region_mode == "downstream":
                if strand == "+":
                    start, end = feature_start, feature_start + extend
                else:
                    start, end = feature_end - extend, feature_end
            elif region_mode == "upstream":
                if strand == "+":
                    start, end = feature_start - extend, feature_start
                else:
                    start, end = feature_end, feature_end + extend
            else:
                raise ValueError(f"Unsupported region mode: {region_mode}")

            start = max(0, start)
            if end - start >= min_region_length:
                regions.append(
                    Region(
                        chrom, start, end, state, strand,
                        anchor_start=feature_start,
                        anchor_end=feature_end,
                    )
                )

    if not regions:
        raise ValueError("No BED regions remained after filtering.")

    return regions


def group_regions_by_state(regions: Iterable[Region]) -> Dict[str, List[Region]]:
    """Group regions by state using deterministic state ordering."""

    grouped: Dict[str, List[Region]] = defaultdict(list)
    for region in regions:
        grouped[region.state].append(region)
    return dict(sorted(grouped.items(), key=lambda item: item[0]))


# -----------------------------------------------------------------------------
# DCC and output helpers
# -----------------------------------------------------------------------------


def signed_opportunity_vector(region_length: int, dmax: int) -> np.ndarray:
    """Return possible base-pair opportunities at each signed lag."""

    lags = np.arange(-dmax, dmax + 1, dtype=np.int64)
    opportunities = region_length - np.abs(lags)
    opportunities[opportunities < 0] = 0
    return opportunities.astype(np.float64)


def signed_opportunity_vector_from_masks(
    valid_a: np.ndarray, valid_b: np.ndarray, dmax: int
) -> np.ndarray:
    """Return signed opportunities for independently masked A/B signals."""
    length = min(valid_a.size, valid_b.size)
    output = np.zeros(2 * dmax + 1, dtype=np.float64)
    if length == 0:
        return output
    a = np.asarray(valid_a[:length], dtype=np.float64)
    b = np.asarray(valid_b[:length], dtype=np.float64)
    full_length = 2 * length - 1
    fft_length = 1 << (full_length - 1).bit_length()
    spectrum_b = np.fft.rfft(b, n=fft_length)
    spectrum_a_reversed = np.fft.rfft(a[::-1], n=fft_length)
    correlation = np.fft.irfft(
        spectrum_b * spectrum_a_reversed, n=fft_length
    )[:full_length]
    maximum = min(dmax, length - 1)
    centre = length - 1
    output[dmax - maximum : dmax + maximum + 1] = np.rint(
        np.maximum(
            0.0,
            correlation[centre - maximum : centre + maximum + 1],
        )
    )
    return output


def collapse_signed_to_absolute(values: np.ndarray, dmax: int) -> np.ndarray:
    """Collapse -dmax..+dmax signed values into 0..dmax distances."""

    absolute = np.zeros(dmax + 1, dtype=np.float64)
    absolute[0] = values[dmax]
    for distance in range(1, dmax + 1):
        absolute[distance] = values[dmax - distance] + values[dmax + distance]
    return absolute


def build_reported_dcc(
    raw_dcc: np.ndarray,
    opportunities: np.ndarray,
    normalize_dcc: bool,
    normalize_by_signal_totals: bool,
    total_signal_a: float,
    total_signal_b: float,
) -> np.ndarray:
    """Build the primary reported DCC vector."""

    reported = raw_dcc.astype(np.float64, copy=True)

    if normalize_dcc:
        normalized = np.zeros_like(reported)
        valid = opportunities > 0
        normalized[valid] = reported[valid] / opportunities[valid]
        reported = normalized

    if normalize_by_signal_totals:
        denominator = total_signal_a * total_signal_b
        if math.isclose(denominator, 0.0, abs_tol=0.0):
            reported = np.zeros_like(reported)
        else:
            reported /= denominator

    return reported


def write_dcc_tsv(
    output_path: str,
    raw_dcc: np.ndarray,
    opportunities: np.ndarray,
    dmax: int,
    signed_lags: bool,
    normalize_dcc: bool,
    normalize_by_signal_totals: bool,
    total_signal_a: float,
    total_signal_b: float,
    cpm_scale: float,
) -> None:
    """Write raw, normalized, percentage and signal-depth-scaled DCC."""

    reported = build_reported_dcc(
        raw_dcc,
        opportunities,
        normalize_dcc,
        normalize_by_signal_totals,
        total_signal_a,
        total_signal_b,
    )

    total = float(np.sum(reported))
    if math.isclose(total, 0.0, abs_tol=0.0):
        percent = np.zeros_like(reported)
    else:
        percent = (reported / total) * 100.0

    signal_pair_denominator = total_signal_a * total_signal_b
    if math.isclose(signal_pair_denominator, 0.0, abs_tol=0.0):
        per_million_pairs = np.zeros_like(raw_dcc)
    else:
        per_million_pairs = (raw_dcc / signal_pair_denominator) * cpm_scale

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    x_label = "Lag" if signed_lags else "Distance"
    x_values = (
        np.arange(-dmax, dmax + 1, dtype=int)
        if signed_lags
        else np.arange(0, dmax + 1, dtype=int)
    )

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                x_label,
                "DCC Value",
                "DCC Value Percent",
                "Raw DCC Value",
                "Opportunities",
                "DCC per million signal-pairs",
            ]
        )

        for row in zip(
            x_values,
            reported,
            percent,
            raw_dcc,
            opportunities,
            per_million_pairs,
        ):
            x, value, pct, raw, opp, scaled = row
            writer.writerow(
                [
                    int(x),
                    f"{value:.12g}",
                    f"{pct:.12g}",
                    f"{raw:.12g}",
                    f"{opp:.12g}",
                    f"{scaled:.12g}",
                ]
            )


def write_shift_summary(
    output_path: str,
    values: np.ndarray,
    dmax: int,
    signed_lags: bool,
    summary_lag_window: int,
) -> None:
    """Write maximum-shift and selected-lag summary values."""

    rows: List[Dict[str, object]] = []

    if signed_lags:
        lags = np.arange(-dmax, dmax + 1, dtype=int)
        if values.size:
            max_idx = int(np.argmax(values))
            rows.append(
                {
                    "Metric": "max_lag_all",
                    "Lag_or_Distance": int(lags[max_idx]),
                    "DCC Value": float(values[max_idx]),
                }
            )

        mask = np.abs(lags) <= summary_lag_window
        if np.any(mask):
            local_lags = lags[mask]
            local_values = values[mask]
            max_idx = int(np.argmax(local_values))
            rows.append(
                {
                    "Metric": f"max_lag_within_{summary_lag_window}bp",
                    "Lag_or_Distance": int(local_lags[max_idx]),
                    "DCC Value": float(local_values[max_idx]),
                }
            )

        for lag in (-20, -10, -5, 0, 5, 10, 20):
            if -dmax <= lag <= dmax:
                rows.append(
                    {
                        "Metric": f"lag_{lag}",
                        "Lag_or_Distance": lag,
                        "DCC Value": float(values[lag + dmax]),
                    }
                )

        for lag in (5, -5, 10, -10, 20, -20):
            if -dmax <= lag <= dmax:
                rows.append(
                    {
                        "Metric": f"lag_{lag}_minus_lag_0",
                        "Lag_or_Distance": lag,
                        "DCC Value": float(values[dmax + lag] - values[dmax]),
                    }
                )
    else:
        distances = np.arange(0, dmax + 1, dtype=int)
        if values.size:
            max_idx = int(np.argmax(values))
            rows.append(
                {
                    "Metric": "max_distance_all",
                    "Lag_or_Distance": int(distances[max_idx]),
                    "DCC Value": float(values[max_idx]),
                }
            )

        mask = distances <= summary_lag_window
        if np.any(mask):
            local_distances = distances[mask]
            local_values = values[mask]
            max_idx = int(np.argmax(local_values))
            rows.append(
                {
                    "Metric": f"max_distance_within_{summary_lag_window}bp",
                    "Lag_or_Distance": int(local_distances[max_idx]),
                    "DCC Value": float(local_values[max_idx]),
                }
            )

        for distance in (0, 5, 10, 20):
            if distance <= dmax:
                rows.append(
                    {
                        "Metric": f"distance_{distance}",
                        "Lag_or_Distance": distance,
                        "DCC Value": float(values[distance]),
                    }
                )

        for distance in (5, 10, 20):
            if distance <= dmax:
                rows.append(
                    {
                        "Metric": f"distance_{distance}_minus_distance_0",
                        "Lag_or_Distance": distance,
                        "DCC Value": float(values[distance] - values[0]),
                    }
                )

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Metric", "Lag_or_Distance", "DCC Value"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary_tsv(output_path: str, rows: Sequence[Mapping[str, object]]) -> None:
    """Write one run-level summary table."""

    if not rows:
        return

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def default_output_prefix(
    inputs_a: Sequence[str],
    inputs_b: Sequence[str],
    regions_bed: Optional[str],
    scope: str,
    selected_chromosomes: Optional[set[str]],
    label_a: Optional[str],
    label_b: Optional[str],
) -> str:
    """Construct an informative output prefix."""

    input_a = label_a or compact_input_label(inputs_a, "A")
    input_b = label_b or compact_input_label(inputs_b, "B")

    if regions_bed:
        region_part = Path(regions_bed).stem
    elif scope == "chromosome" and selected_chromosomes:
        region_part = "_".join(sorted(selected_chromosomes))
    else:
        region_part = "combined_chromosomes"

    return sanitize_filename(f"{input_a}_vs_{input_b}_{region_part}")


# -----------------------------------------------------------------------------
# Sparse and FFT cross-correlation
# -----------------------------------------------------------------------------


def update_sparse_signed_dcc(
    dcc: np.ndarray,
    positions_a: np.ndarray,
    values_a: np.ndarray,
    positions_b: np.ndarray,
    values_b: np.ndarray,
    dmax: int,
) -> None:
    """Update signed DCC from sparse sorted position/value arrays."""

    if positions_a.size == 0 or positions_b.size == 0:
        return

    for pos_a, value_a in zip(positions_a, values_a):
        left = np.searchsorted(positions_b, pos_a - dmax, side="left")
        right = np.searchsorted(positions_b, pos_a + dmax, side="right")
        if right <= left:
            continue

        lags = positions_b[left:right] - pos_a
        products = value_a * values_b[left:right]
        np.add.at(dcc, lags + dmax, products)


def update_dense_dcc_sparse(
    dcc: np.ndarray,
    values_a: np.ndarray,
    values_b: np.ndarray,
    dmax: int,
) -> None:
    """Update signed DCC from dense arrays using only non-zero positions."""

    positions_a = np.flatnonzero(values_a).astype(np.int64)
    positions_b = np.flatnonzero(values_b).astype(np.int64)
    if positions_a.size == 0 or positions_b.size == 0:
        return

    update_sparse_signed_dcc(
        dcc,
        positions_a,
        values_a[positions_a],
        positions_b,
        values_b[positions_b],
        dmax,
    )


def update_dense_dcc_fft(
    dcc: np.ndarray,
    values_a: np.ndarray,
    values_b: np.ndarray,
    dmax: int,
) -> None:
    """Update signed DCC with zero-padded FFT cross-correlation."""

    if values_a.size != values_b.size:
        raise ValueError("FFT DCC requires equal-length arrays.")

    n = values_a.size
    if n == 0:
        return

    max_lag = min(dmax, n - 1)
    full_length = 2 * n - 1
    fft_length = 1 << (full_length - 1).bit_length()

    spectrum_b = np.fft.rfft(values_b, n=fft_length)
    spectrum_a_reversed = np.fft.rfft(values_a[::-1], n=fft_length)
    full_correlation = np.fft.irfft(
        spectrum_b * spectrum_a_reversed,
        n=fft_length,
    )[:full_length]

    centre = n - 1
    dcc[dmax - max_lag : dmax + max_lag + 1] += full_correlation[
        centre - max_lag : centre + max_lag + 1
    ]


def update_dense_dcc(
    dcc: np.ndarray,
    values_a: np.ndarray,
    values_b: np.ndarray,
    dmax: int,
    algorithm: str,
    sparse_threshold: float,
) -> str:
    """Select and run sparse or FFT DCC for dense signal arrays."""

    min_length = min(values_a.size, values_b.size)
    if min_length == 0:
        return "none"

    values_a = values_a[:min_length]
    values_b = values_b[:min_length]

    selected = algorithm
    if algorithm == "auto":
        density_a = float(np.count_nonzero(values_a)) / float(min_length)
        density_b = float(np.count_nonzero(values_b)) / float(min_length)
        selected = "sparse" if max(density_a, density_b) <= sparse_threshold else "fft"

    if selected == "sparse":
        update_dense_dcc_sparse(dcc, values_a, values_b, dmax)
    elif selected == "fft":
        update_dense_dcc_fft(dcc, values_a, values_b, dmax)
    else:
        raise ValueError(f"Unsupported DCC algorithm: {selected}")

    return selected


# -----------------------------------------------------------------------------
# bigWig mode
# -----------------------------------------------------------------------------


def open_bigwig_tracks(paths: Sequence[str]) -> List[BigWigTrack]:
    """Open bigWig files using pyBigWig."""

    try:
        import pyBigWig  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "pyBigWig is required. Install it with 'conda install -c bioconda "
            "pybigwig' or 'pip install pyBigWig'."
        ) from exc

    tracks: List[BigWigTrack] = []
    try:
        for path in paths:
            handle = pyBigWig.open(path)
            if handle is None or not handle.isBigWig():
                raise ValueError(f"Not a readable bigWig: {path}")
            tracks.append(
                BigWigTrack(path, handle, handle.chroms(), {})
            )
    except Exception:
        close_bigwig_tracks(tracks)
        raise

    return tracks


def close_bigwig_tracks(tracks: Iterable[BigWigTrack]) -> None:
    """Close bigWig handles."""

    for track in tracks:
        try:
            track.handle.close()
        except Exception:
            pass


def resolve_bigwig_chromosome(track: BigWigTrack, chrom: str) -> Optional[str]:
    """Resolve a region contig against a bigWig chromosome dictionary."""

    if chrom in track.chrom_cache:
        return track.chrom_cache[chrom]

    resolved = None
    for candidate in candidate_chromosome_aliases(chrom):
        if candidate in track.chrom_sizes:
            resolved = candidate
            break

    track.chrom_cache[chrom] = resolved
    return resolved


def read_bigwig_region(
    track: BigWigTrack,
    region: Region,
    value_limit: Optional[float],
) -> Tuple[Optional[np.ndarray], bool]:
    """Read one bigWig region and report whether it was clipped."""

    bigwig_chrom = resolve_bigwig_chromosome(track, region.chrom)
    if bigwig_chrom is None:
        return None, False

    chrom_size = int(track.chrom_sizes[bigwig_chrom])
    start = max(0, region.start)
    end = min(region.end, chrom_size)
    clipped = start != region.start or end != region.end

    if end <= start:
        return np.empty(0, dtype=np.float64), clipped

    values = track.handle.values(bigwig_chrom, start, end, numpy=True)
    if values is None:
        values = np.zeros(end - start, dtype=np.float64)
    else:
        values = np.asarray(values, dtype=np.float64)
        np.nan_to_num(values, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

    if value_limit is not None:
        np.clip(values, -value_limit, value_limit, out=values)

    return values, clipped


def combine_bigwig_region_signals(
    tracks: Sequence[BigWigTrack],
    region: Region,
    value_limit: Optional[float],
) -> Tuple[Optional[np.ndarray], bool, bool]:
    """Sum matching bigWig signals within one region."""

    combined: Optional[np.ndarray] = None
    matched = False
    clipped = False

    for track in tracks:
        values, track_clipped = read_bigwig_region(track, region, value_limit)
        if values is None:
            continue

        matched = True
        clipped = clipped or track_clipped

        if combined is None:
            combined = values.astype(np.float64, copy=True)
        else:
            min_length = min(combined.size, values.size)
            combined = combined[:min_length]
            combined += values[:min_length]

    if not matched:
        return None, True, False

    if region.strand == "-" and combined is not None:
        combined = combined[::-1].copy()

    return combined, False, clipped


def process_bigwig_group(
    files_a: Sequence[str],
    files_b: Sequence[str],
    regions_by_state: Mapping[str, Sequence[Region]],
    output_prefix: str,
    output_dir: str,
    args: argparse.Namespace,
) -> List[Mapping[str, object]]:
    """Calculate all state outputs for a bigWig A-versus-B comparison."""

    tracks_a = open_bigwig_tracks(files_a)
    tracks_b = open_bigwig_tracks(files_b)
    summary_rows: List[Mapping[str, object]] = []
    blacklist = load_blacklist_unbounded(getattr(args, "blacklist_bed", None))

    try:
        for state, regions in regions_by_state.items():
            if not args.quiet:
                print(
                    f"State/group: {state} ({len(regions):,} regions; "
                    f"A={len(tracks_a):,} bigWig(s); B={len(tracks_b):,})",
                    flush=True,
                )

            raw_signed = np.zeros(2 * args.dmax + 1, dtype=np.float64)
            signed_opportunities = np.zeros_like(raw_signed)
            stats = DccStats()
            algorithm_counts = {"sparse": 0, "fft": 0, "none": 0}

            for region_index, region in enumerate(regions, start=1):
                stats.regions_seen += 1

                values_a, missing_a, clipped_a = combine_bigwig_region_signals(
                    tracks_a, region, args.value_limit
                )
                values_b, missing_b, clipped_b = combine_bigwig_region_signals(
                    tracks_b, region, args.value_limit
                )

                if missing_a:
                    stats.missing_chromosome_a += 1
                if missing_b:
                    stats.missing_chromosome_b += 1
                if clipped_a:
                    stats.clipped_regions_a += 1
                if clipped_b:
                    stats.clipped_regions_b += 1

                if values_a is not None:
                    stats.regions_with_a += 1

                if values_b is not None:
                    stats.regions_with_b += 1

                if values_a is None or values_b is None:
                    continue

                region_length = min(values_a.size, values_b.size)
                if region_length < args.min_region_length:
                    stats.short_or_empty += 1
                    continue
                values_a = values_a[:region_length]
                values_b = values_b[:region_length]
                blacklist_valid = np.ones(region_length, dtype=bool)
                if blacklist is not None:
                    blacklist_valid = blacklist.valid_mask(
                        region.chrom,
                        max(0, region.start),
                        max(0, region.start) + region_length,
                    )
                    if region.strand == "-":
                        blacklist_valid = blacklist_valid[::-1]
                    stats.blacklisted_positions += int(
                        np.count_nonzero(~blacklist_valid)
                    )
                valid_a = np.isfinite(values_a) & blacklist_valid
                valid_b = np.isfinite(values_b) & blacklist_valid
                stats.signal_positions_a += int(np.count_nonzero(valid_a))
                stats.signal_positions_b += int(np.count_nonzero(valid_b))
                stats.nonzero_positions_a += int(np.count_nonzero(values_a[valid_a]))
                stats.nonzero_positions_b += int(np.count_nonzero(values_b[valid_b]))
                stats.total_signal_a += float(np.sum(values_a[valid_a]))
                stats.total_signal_b += float(np.sum(values_b[valid_b]))
                if np.count_nonzero(valid_a) < args.min_region_length or np.count_nonzero(valid_b) < args.min_region_length:
                    stats.short_or_empty += 1
                    continue
                working_a = np.where(valid_a, values_a, 0.0)
                working_b = np.where(valid_b, values_b, 0.0)
                signed_opportunities += signed_opportunity_vector_from_masks(
                    valid_a, valid_b, args.dmax
                )
                selected = update_dense_dcc(
                    raw_signed,
                    working_a,
                    working_b,
                    args.dmax,
                    args.algorithm,
                    args.sparse_threshold,
                )
                algorithm_counts[selected] += 1
                stats.regions_with_both += 1

                if (
                    not args.quiet
                    and args.progress_every > 0
                    and region_index % args.progress_every == 0
                ):
                    print(
                        f"  Processed {region_index:,}/{len(regions):,} regions; "
                        f"used {stats.regions_with_both:,}",
                        flush=True,
                    )

            raw_out, opportunities_out = prepare_output_vectors(
                raw_signed,
                signed_opportunities,
                args.dmax,
                args.signed_lags,
            )
            output_path, shift_path = write_state_outputs(
                mode="bigwig",
                state=state,
                signal_label=None,
                raw_dcc=raw_out,
                opportunities=opportunities_out,
                output_prefix=output_prefix,
                output_dir=output_dir,
                stats=stats,
                args=args,
            )

            summary_rows.append(
                {
                    "State": state,
                    "Output": output_path,
                    "Shift summary": shift_path,
                    "Mode": "bigwig",
                    "A files": len(files_a),
                    "B files": len(files_b),
                    "Regions": stats.regions_seen,
                    "Used regions": stats.regions_with_both,
                    "Regions with A": stats.regions_with_a,
                    "Regions with B": stats.regions_with_b,
                    "Missing chromosome A": stats.missing_chromosome_a,
                    "Missing chromosome B": stats.missing_chromosome_b,
                    "Short or empty": stats.short_or_empty,
                    "Clipped regions A": stats.clipped_regions_a,
                    "Clipped regions B": stats.clipped_regions_b,
                    "Signal positions A": stats.signal_positions_a,
                    "Signal positions B": stats.signal_positions_b,
                    "Non-zero positions A": stats.nonzero_positions_a,
                    "Non-zero positions B": stats.nonzero_positions_b,
                    "Total signal A": f"{stats.total_signal_a:.12g}",
                    "Total signal B": f"{stats.total_signal_b:.12g}",
                    "Blacklisted positions": stats.blacklisted_positions,
                    "Blacklisted anchors excluded": getattr(args, "_blacklisted_anchor_exclusions", 0),
                    "Sparse calculations": algorithm_counts["sparse"],
                    "FFT calculations": algorithm_counts["fft"],
                }
            )

            if not args.quiet:
                print(f"  Wrote: {output_path}")
                print(f"  Shift summary: {shift_path}")
                print(
                    f"  Algorithms: sparse={algorithm_counts['sparse']:,}, "
                    f"FFT={algorithm_counts['fft']:,}",
                    flush=True,
                )
    finally:
        close_bigwig_tracks(tracks_a)
        close_bigwig_tracks(tracks_b)

    return summary_rows


# -----------------------------------------------------------------------------
# BAM mode
# -----------------------------------------------------------------------------


def open_bam_tracks(paths: Sequence[str]) -> List[BamTrack]:
    """Open BAM files using pysam."""

    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "pysam is required. Install it with 'conda install -c bioconda "
            "pysam' or 'pip install pysam'."
        ) from exc

    tracks: List[BamTrack] = []
    try:
        for path in paths:
            handle = pysam.AlignmentFile(path, "rb")
            tracks.append(
                BamTrack(path, handle, set(handle.references), {})
            )
    except Exception:
        close_bam_tracks(tracks)
        raise

    return tracks


def close_bam_tracks(tracks: Iterable[BamTrack]) -> None:
    """Close BAM handles."""

    for track in tracks:
        try:
            track.handle.close()
        except Exception:
            pass


def resolve_bam_chromosome(track: BamTrack, chrom: str) -> Optional[str]:
    """Resolve a region contig against a BAM header."""

    if chrom in track.chrom_cache:
        return track.chrom_cache[chrom]

    resolved = None
    for candidate in candidate_chromosome_aliases(chrom):
        if candidate in track.references:
            resolved = candidate
            break

    track.chrom_cache[chrom] = resolved
    return resolved


def read_passes_filters(
    read,
    mapq: int,
    require_proper_pairs: bool,
    include_duplicate_flag: bool,
) -> bool:
    """Return whether a TLEN-positive alignment represents an accepted fragment."""

    if read.is_unmapped or read.mate_is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    if read.is_duplicate and not include_duplicate_flag:
        return False
    if require_proper_pairs and not read.is_proper_pair:
        return False
    if read.mapping_quality < mapq:
        return False
    if read.template_length <= 0:
        return False
    return True


def fragment_positions_and_weights(
    fragment_start: int,
    fragment_length: int,
    position_type: str,
) -> List[Tuple[int, float]]:
    """Return dyad or fragment-end positions and weights."""

    if position_type == "left_end":
        return [(fragment_start, 1.0)]
    if position_type == "right_end":
        return [(fragment_start + fragment_length - 1, 1.0)]
    if position_type != "dyad":
        raise ValueError(f"Unsupported fragment position: {position_type}")

    if fragment_length % 2:
        return [(fragment_start + fragment_length // 2, 1.0)]

    return [
        (fragment_start + fragment_length // 2 - 1, 0.5),
        (fragment_start + fragment_length // 2, 0.5),
    ]


def group_positions_values(
    positions: np.ndarray,
    values: np.ndarray,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Collapse repeated genomic positions by summing their weights."""

    if positions.size == 0:
        return None, None

    order = np.argsort(positions)
    positions = positions[order]
    values = values[order]
    unique_positions, first_indices = np.unique(positions, return_index=True)
    summed_values = np.add.reduceat(values, first_indices)
    return unique_positions.astype(np.int64), summed_values.astype(np.float64)


def reverse_sparse_signal(
    positions: Optional[np.ndarray],
    values: Optional[np.ndarray],
    region: Region,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Reverse sparse positions within a minus-strand region."""

    if positions is None or values is None or region.strand != "-":
        return positions, values

    reversed_positions = region.start + region.end - 1 - positions
    order = np.argsort(reversed_positions)
    return reversed_positions[order], values[order]


def extract_bam_region_signal(
    tracks: Sequence[BamTrack],
    region: Region,
    length_spec: FragmentLengthSpec,
    position_type: str,
    mapq: int,
    require_proper_pairs: bool,
    include_duplicate_flag: bool,
    max_duplicates: int,
    blacklist: BlacklistIndex | None = None,
    exclusion_counter: list[int] | None = None,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], bool]:
    """Extract one grouped fragment-position signal across BAM files."""

    all_positions: List[np.ndarray] = []
    all_values: List[np.ndarray] = []
    matched_chromosome = False
    fetch_margin = length_spec.max_length + 5
    fetch_start = max(0, region.start - fetch_margin)
    fetch_end = region.end + fetch_margin

    for track in tracks:
        bam_chrom = resolve_bam_chromosome(track, region.chrom)
        if bam_chrom is None:
            continue

        matched_chromosome = True
        coordinate_counts: Dict[Tuple[str, int, int], int] = defaultdict(int)
        positions: List[int] = []
        weights: List[float] = []

        try:
            iterator = track.handle.fetch(bam_chrom, fetch_start, fetch_end)
        except ValueError:
            continue

        for read in iterator:
            if not read_passes_filters(
                read,
                mapq,
                require_proper_pairs,
                include_duplicate_flag,
            ):
                continue

            fragment_length = int(read.template_length)
            if not length_spec.matches(fragment_length):
                continue

            fragment_start = int(read.reference_start)
            fragment_end = fragment_start + fragment_length
            if fragment_end <= fragment_start:
                continue
            if blacklist is not None and blacklist.overlaps(
                region.chrom, fragment_start, fragment_end
            ):
                if exclusion_counter is not None:
                    exclusion_counter[0] += 1
                continue

            if max_duplicates > 0:
                key = (bam_chrom, fragment_start, fragment_end)
                coordinate_counts[key] += 1
                if coordinate_counts[key] > max_duplicates:
                    continue

            for position, weight in fragment_positions_and_weights(
                fragment_start,
                fragment_length,
                position_type,
            ):
                if region.start <= position < region.end:
                    positions.append(position)
                    weights.append(weight)

        if positions:
            all_positions.append(np.asarray(positions, dtype=np.int64))
            all_values.append(np.asarray(weights, dtype=np.float64))

    if not matched_chromosome:
        return None, None, True
    if not all_positions:
        return None, None, False

    grouped_positions, grouped_values = group_positions_values(
        np.concatenate(all_positions),
        np.concatenate(all_values),
    )
    grouped_positions, grouped_values = reverse_sparse_signal(
        grouped_positions,
        grouped_values,
        region,
    )
    return grouped_positions, grouped_values, False


def extract_interval_region_signal(
    source: IntervalFragmentSource,
    region: Region,
    length_spec: FragmentLengthSpec,
    position_type: str,
    max_duplicates: int,
    blacklist: BlacklistIndex | None = None,
    exclusion_counter: list[int] | None = None,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], bool]:
    """Extract a sparse signal from BED/BED.gz/bigBed fragment intervals."""
    resolved = None
    for candidate in candidate_chromosome_aliases(region.chrom):
        if candidate in source.references:
            resolved = candidate
            break
    if resolved is None:
        return None, None, True
    margin = length_spec.max_length + 5
    fragments = source.fetch(
        resolved,
        max(0, region.start - margin),
        region.end + margin,
        max_per_coordinate=max_duplicates,
        subsample=None,
        dedup_scope="per_bam",
    )
    positions: List[int] = []
    weights: List[float] = []
    for fragment_start, fragment_end in fragments:
        if blacklist is not None and blacklist.overlaps(
            region.chrom, fragment_start, fragment_end
        ):
            if exclusion_counter is not None:
                exclusion_counter[0] += 1
            continue
        fragment_length = fragment_end - fragment_start
        if not length_spec.matches(fragment_length):
            continue
        for position, weight in fragment_positions_and_weights(
            fragment_start, fragment_length, position_type
        ):
            if region.start <= position < region.end:
                positions.append(position)
                weights.append(weight)
    if not positions:
        return None, None, False
    grouped_positions, grouped_values = group_positions_values(
        np.asarray(positions, dtype=np.int64),
        np.asarray(weights, dtype=np.float64),
    )
    return (*reverse_sparse_signal(grouped_positions, grouped_values, region), False)


def process_bam_group(
    files_a: Sequence[str],
    files_b: Sequence[str],
    regions_by_state: Mapping[str, Sequence[Region]],
    output_prefix: str,
    output_dir: str,
    args: argparse.Namespace,
) -> List[Mapping[str, object]]:
    """Calculate state outputs from BAM or fragment-coordinate signal pairs."""

    use_fragment_a = bool(getattr(args, "fragments_a", None))
    use_fragment_b = bool(getattr(args, "fragments_b", None))
    tracks_a = IntervalFragmentSource(files_a, chrom_sizes=args.chrom_sizes) if use_fragment_a else open_bam_tracks(files_a)
    tracks_b = IntervalFragmentSource(files_b, chrom_sizes=args.chrom_sizes) if use_fragment_b else open_bam_tracks(files_b)
    spec_a = FragmentLengthSpec(args.min_length_a, args.max_length_a)
    spec_b = FragmentLengthSpec(args.min_length_b, args.max_length_b)
    signal_label = (
        f"A{spec_a.label}_{args.position_a}_"
        f"B{spec_b.label}_{args.position_b}"
    )
    summary_rows: List[Mapping[str, object]] = []
    blacklist = load_blacklist_unbounded(getattr(args, "blacklist_bed", None))

    try:
        for state, regions in regions_by_state.items():
            if not args.quiet:
                print(
                    f"State/group: {state} ({len(regions):,} regions; "
                    f"A={len(files_a):,} input(s); B={len(files_b):,})",
                    flush=True,
                )

            raw_signed = np.zeros(2 * args.dmax + 1, dtype=np.float64)
            signed_opportunities = np.zeros_like(raw_signed)
            stats = DccStats()
            excluded_a = [0]
            excluded_b = [0]

            for region_index, region in enumerate(regions, start=1):
                stats.regions_seen += 1

                if use_fragment_a:
                    positions_a, values_a, missing_a = extract_interval_region_signal(
                        tracks_a, region, spec_a, args.position_a, args.max_duplicates,
                        blacklist, excluded_a,
                    )
                else:
                    positions_a, values_a, missing_a = extract_bam_region_signal(
                        tracks_a,
                        region,
                        spec_a,
                        args.position_a,
                        args.mapq,
                        not args.no_require_proper_pairs,
                        args.include_duplicate_flag,
                        args.max_duplicates,
                        blacklist,
                        excluded_a,
                    )

                if use_fragment_b:
                    positions_b, values_b, missing_b = extract_interval_region_signal(
                        tracks_b, region, spec_b, args.position_b, args.max_duplicates,
                        blacklist, excluded_b,
                    )
                else:
                    positions_b, values_b, missing_b = extract_bam_region_signal(
                        tracks_b,
                        region,
                        spec_b,
                        args.position_b,
                        args.mapq,
                        not args.no_require_proper_pairs,
                        args.include_duplicate_flag,
                        args.max_duplicates,
                        blacklist,
                        excluded_b,
                    )

                if missing_a:
                    stats.missing_chromosome_a += 1
                if missing_b:
                    stats.missing_chromosome_b += 1

                if positions_a is not None and values_a is not None:
                    stats.regions_with_a += 1
                    stats.signal_positions_a += int(positions_a.size)
                    stats.nonzero_positions_a += int(np.count_nonzero(values_a))
                    stats.total_signal_a += float(np.sum(values_a))

                if positions_b is not None and values_b is not None:
                    stats.regions_with_b += 1
                    stats.signal_positions_b += int(positions_b.size)
                    stats.nonzero_positions_b += int(np.count_nonzero(values_b))
                    stats.total_signal_b += float(np.sum(values_b))

                if (
                    positions_a is None
                    or values_a is None
                    or positions_b is None
                    or values_b is None
                ):
                    continue
                if region.length < args.min_region_length:
                    stats.short_or_empty += 1
                    continue

                valid_mask = (
                    blacklist.valid_mask(region.chrom, region.start, region.end)
                    if blacklist is not None
                    else np.ones(region.length, dtype=bool)
                )
                if region.strand == "-":
                    valid_mask = valid_mask[::-1]
                stats.blacklisted_positions += int(np.count_nonzero(~valid_mask))
                signed_opportunities += signed_opportunity_vector_from_masks(
                    valid_mask, valid_mask, args.dmax
                )
                update_sparse_signed_dcc(
                    raw_signed,
                    positions_a,
                    values_a,
                    positions_b,
                    values_b,
                    args.dmax,
                )
                stats.regions_with_both += 1

                if (
                    not args.quiet
                    and args.progress_every > 0
                    and region_index % args.progress_every == 0
                ):
                    print(
                        f"  Processed {region_index:,}/{len(regions):,} regions; "
                        f"used {stats.regions_with_both:,}",
                        flush=True,
                    )

            raw_out, opportunities_out = prepare_output_vectors(
                raw_signed,
                signed_opportunities,
                args.dmax,
                args.signed_lags,
            )
            output_path, shift_path = write_state_outputs(
                mode=(
                    "fragments"
                    if use_fragment_a and use_fragment_b
                    else "bam"
                    if not use_fragment_a and not use_fragment_b
                    else "mixed"
                ),
                state=state,
                signal_label=signal_label,
                raw_dcc=raw_out,
                opportunities=opportunities_out,
                output_prefix=output_prefix,
                output_dir=output_dir,
                stats=stats,
                args=args,
            )

            summary_rows.append(
                {
                    "State": state,
                    "Output": output_path,
                    "Shift summary": shift_path,
                    "Mode": (
                        "fragments"
                        if use_fragment_a and use_fragment_b
                        else "bam"
                        if not use_fragment_a and not use_fragment_b
                        else "mixed"
                    ),
                    "Blacklisted fragments A": excluded_a[0],
                    "Blacklisted fragments B": excluded_b[0],
                    "Blacklisted positions": stats.blacklisted_positions,
                    "Blacklisted anchors excluded": getattr(args, "_blacklisted_anchor_exclusions", 0),
                    "Input type A": "fragments" if use_fragment_a else "bam",
                    "Input type B": "fragments" if use_fragment_b else "bam",
                    "A files": len(files_a),
                    "B files": len(files_b),
                    "Length A": spec_a.label,
                    "Position A": args.position_a,
                    "Length B": spec_b.label,
                    "Position B": args.position_b,
                    "Regions": stats.regions_seen,
                    "Used regions": stats.regions_with_both,
                    "Regions with A": stats.regions_with_a,
                    "Regions with B": stats.regions_with_b,
                    "Missing chromosome A": stats.missing_chromosome_a,
                    "Missing chromosome B": stats.missing_chromosome_b,
                    "Short or empty": stats.short_or_empty,
                    "Signal positions A": stats.signal_positions_a,
                    "Signal positions B": stats.signal_positions_b,
                    "Non-zero positions A": stats.nonzero_positions_a,
                    "Non-zero positions B": stats.nonzero_positions_b,
                    "Total signal A": f"{stats.total_signal_a:.12g}",
                    "Total signal B": f"{stats.total_signal_b:.12g}",
                }
            )

            if not args.quiet:
                print(f"  Wrote: {output_path}")
                print(f"  Shift summary: {shift_path}", flush=True)
    finally:
        if use_fragment_a:
            tracks_a.close()
        else:
            close_bam_tracks(tracks_a)
        if use_fragment_b:
            tracks_b.close()
        else:
            close_bam_tracks(tracks_b)

    return summary_rows


# -----------------------------------------------------------------------------
# Shared state-output helpers
# -----------------------------------------------------------------------------


def prepare_output_vectors(
    raw_signed: np.ndarray,
    signed_opportunities: np.ndarray,
    dmax: int,
    signed_lags: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    """Keep signed values or collapse them to absolute distances."""

    if signed_lags:
        return raw_signed, signed_opportunities
    return (
        collapse_signed_to_absolute(raw_signed, dmax),
        collapse_signed_to_absolute(signed_opportunities, dmax),
    )


def plot_dcc_tsv(tsv_path: str, output_path: str, title: str | None = None) -> Path | None:
    """Plot the primary DCC values written to a TSV."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    data = np.genfromtxt(tsv_path, names=True, delimiter="\t", dtype=None, encoding="utf-8")
    if getattr(data, "size", 0) == 0:
        return
    names = list(data.dtype.names or ())
    x_name = names[0]
    y_name = next((name for name in names if name.replace("_", " ") == "DCC Value"), names[1])
    x = np.atleast_1d(data[x_name]).astype(float)
    y = np.atleast_1d(data[y_name]).astype(float)
    signed_lags = bool(np.any(x < 0))
    plot_mask = np.abs(x) != 1 if signed_lags else x != 1
    x_plot = x[plot_mask]
    y_plot = y[plot_mask]
    if x_plot.size == 0:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x_plot, y_plot, linewidth=1.2, marker="o", markersize=2.0, markeredgewidth=0)
    from nucleosuite.plotting import apply_base_pair_x_axis
    apply_base_pair_x_axis(ax, x_plot)
    ax.set_xlabel("Lag (bp)" if signed_lags else "Distance (bp)")
    ax.set_ylabel("DCC value")
    if title:
        ax.set_title(title)
    ax.grid(axis="x", alpha=0.5)
    from nucleosuite.plotting import annotate_points, get_plot_options, save_figure
    options = get_plot_options()
    if options.label_points == "peaks" and x_plot.size >= 3:
        from nucleosuite.nrl import local_maximum_indices, retain_separated_peaks
        candidates = local_maximum_indices(x_plot, y_plot, min_rise_bp=2.0, min_fall_bp=2.0)
        called = retain_separated_peaks(candidates, x_plot, y_plot, y_plot, min_separation=100.0)
        annotate_points(
            ax, [peak.distance for peak in called], [peak.smoothed_value for peak in called],
            points_are_peaks=True, options=options,
        )
    fig.tight_layout()
    saved = save_figure(fig, output_path, default_dpi=220, bbox_inches="tight")
    plt.close(fig)
    return saved


def write_state_outputs(
    mode: str,
    state: str,
    signal_label: Optional[str],
    raw_dcc: np.ndarray,
    opportunities: np.ndarray,
    output_prefix: str,
    output_dir: str,
    stats: DccStats,
    args: argparse.Namespace,
) -> Tuple[str, str]:
    """Write DCC and shift-summary outputs for one state."""

    normalization = "raw" if args.no_normalize_dcc else "opportunity_normalized"
    lag_mode = "signed_lags" if args.signed_lags else "absolute_distances"
    signal_normalization = (
        "_signal_total_normalized" if args.normalize_by_signal_totals else ""
    )
    signal_component = f"_{sanitize_filename(signal_label)}" if signal_label else ""
    base = (
        f"{sanitize_filename(output_prefix)}_{sanitize_filename(state)}_"
        f"{mode}_DCC{signal_component}_{lag_mode}_{normalization}"
        f"{signal_normalization}"
    )

    output_path = os.path.join(output_dir, f"{base}.tsv")
    shift_path = os.path.join(output_dir, f"{base}_shift_summary.tsv")

    write_dcc_tsv(
        output_path,
        raw_dcc,
        opportunities,
        args.dmax,
        args.signed_lags,
        not args.no_normalize_dcc,
        args.normalize_by_signal_totals,
        stats.total_signal_a,
        stats.total_signal_b,
        args.cpm_scale,
    )
    plot_dcc_tsv(output_path, os.path.splitext(output_path)[0] + ".png", title=base)

    reported = build_reported_dcc(
        raw_dcc,
        opportunities,
        not args.no_normalize_dcc,
        args.normalize_by_signal_totals,
        stats.total_signal_a,
        stats.total_signal_b,
    )
    write_shift_summary(
        shift_path,
        reported,
        args.dmax,
        args.signed_lags,
        args.summary_lag_window,
    )
    return output_path, shift_path


# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------


def add_region_arguments(parser: argparse.ArgumentParser) -> None:
    """Add region-source and BED grouping arguments."""

    region_source = parser.add_mutually_exclusive_group(required=True)
    region_source.add_argument(
        "--regions-bed",
        help=(
            "BED analysis regions. Columns 1-3 are required; column 4 is used "
            "as the state by default when present."
        ),
    )
    region_source.add_argument(
        "--chrom-sizes",
        help="Two-column chromosome-size table, BAM, or CRAM used to generate windows.",
    )

    parser.add_argument(
        "--scope",
        choices=["combined_chromosomes", "genome", "chromosome"],
        default="combined_chromosomes",
        help="Scope for --chrom-sizes mode. Default: combined_chromosomes.",
    )
    parser.add_argument(
        "--chromosome",
        action="append",
        default=None,
        help=(
            "Chromosome/contig to include. May be repeated or comma-separated. "
            "Required with --scope chromosome; also filters BED input."
        ),
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=100_000,
        help="Generated window size. Default: 100000.",
    )
    parser.add_argument(
        "--state-name",
        default="Combined chromosomes",
        help="State for generated windows or BED rows without a state. Default: Combined chromosomes.",
    )
    parser.add_argument(
        "--state-column",
        type=int,
        default=4,
        help="1-based BED state/group column. Default: 4.",
    )
    parser.add_argument(
        "--strand-column",
        type=int,
        default=6,
        help="1-based BED strand column. Default: 6.",
    )
    parser.add_argument(
        "--region-mode",
        choices=["interval", "downstream", "upstream"],
        default="interval",
        help="Use intervals directly or construct strand-aware windows. Default: interval.",
    )
    parser.add_argument(
        "--extend",
        type=int,
        default=2000,
        help="Window length for upstream/downstream modes. Default: 2000.",
    )
    parser.add_argument(
        "--strands",
        choices=["plus", "minus", "both"],
        default="both",
        help="BED strands to retain. Default: both.",
    )
    parser.add_argument(
        "--category",
        action="append",
        default=None,
        metavar="CATEGORY=MATCHER[,MATCHER...]",
        help=(
            "Group BED states using exact:STATE, prefix:PREFIX, regex:PATTERN, "
            "or an unqualified exact state. May be repeated."
        ),
    )


def add_shared_dcc_arguments(parser: argparse.ArgumentParser) -> None:
    """Add calculation, normalization and output arguments."""

    parser.add_argument(
        "--blacklist-bed",
        help=(
            "BED blacklist. Overlapping BED anchors and complete fragments are "
            "excluded; BigWig bases are removed from signal opportunities."
        ),
    )
    parser.add_argument(
        "--dmax",
        type=int,
        default=500,
        help="Maximum signed lag or absolute distance, inclusive. Default: 500.",
    )
    parser.add_argument(
        "--min-region-length",
        type=int,
        default=None,
        help="Minimum region length. Default: max(2, dmax + 1).",
    )
    parser.add_argument(
        "--summary-lag-window",
        type=int,
        default=50,
        help="Window for maximum-lag/distance summary. Default: 50.",
    )
    parser.add_argument(
        "--signed-lags",
        action="store_true",
        help="Retain directional signed lags instead of absolute distances.",
    )
    parser.add_argument(
        "--no-normalize-dcc",
        action="store_true",
        help="Use raw DCC as the main DCC Value instead of opportunity normalization.",
    )
    parser.add_argument(
        "--normalize-by-signal-totals",
        action="store_true",
        help="Also divide the main DCC Value by total signal A times total signal B.",
    )
    parser.add_argument(
        "--cpm-scale",
        type=float,
        default=1_000_000.0,
        help="Scale for DCC per million signal-pairs. Default: 1000000.",
    )
    parser.add_argument("--label-a", default=None, help="Optional label for signal A.")
    parser.add_argument("--label-b", default=None, help="Optional label for signal B.")
    parser.add_argument(
        "--out-prefix",
        default=None,
        help="Output prefix. Default is derived from inputs and regions.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Output directory. Default: current directory.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="Print progress every N regions; 0 disables it. Default: 1000.",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output.")


def build_parser() -> argparse.ArgumentParser:
    """Construct the DCC command parser."""

    parser = argparse.ArgumentParser(
        prog="nucleosuite dcc",
        description="Calculate Distance Cross-Correlation from bigWig or BAM signals.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # bigWig DCC across genome windows
  nucleosuite dcc bigwig \
      --bigwig-a 'sample147_chr*.bw' \
      --bigwig-b 'sample167_chr*.bw' \
      --chrom-sizes genome.chrom.sizes \
      --scope combined_chromosomes --dmax 50

  # BAM dyad DCC for one chromosome
  nucleosuite dcc bam \
      --bam-a '/data/A/*.bam' --length-a 147 --position-a dyad \
      --bam-b '/data/B/*.bam' --length-b 167 --position-b dyad \
      --chrom-sizes genome.chrom.sizes \
      --scope chromosome --chromosome chr20 --dmax 50

  # User-defined state categories
  nucleosuite dcc bigwig \
      --bigwig-a A.bw --bigwig-b B.bw --regions-bed states.bed \
      --category 'Open=exact:Promoter,exact:Enhancer' \
      --category 'Repressed=prefix:13_,prefix:14_'
""",
    )

    subparsers = parser.add_subparsers(dest="mode", required=True)

    bigwig = subparsers.add_parser("bigwig", help="DCC from bigWig signals.")
    bigwig.add_argument(
        "--bigwig-a",
        nargs="+",
        required=True,
        help="bigWig paths or glob patterns for signal A.",
    )
    bigwig.add_argument(
        "--bigwig-b",
        nargs="+",
        required=True,
        help="bigWig paths or glob patterns for signal B.",
    )
    bigwig.add_argument(
        "--value-limit",
        type=float,
        default=None,
        help="Optional absolute cap for bigWig values.",
    )
    bigwig.add_argument(
        "--algorithm",
        choices=["auto", "sparse", "fft"],
        default="auto",
        help="bigWig DCC algorithm. Default: auto.",
    )
    bigwig.add_argument(
        "--sparse-threshold",
        type=float,
        default=0.10,
        help="Maximum non-zero fraction for auto sparse calculation. Default: 0.10.",
    )
    add_region_arguments(bigwig)
    add_shared_dcc_arguments(bigwig)

    bam = subparsers.add_parser("bam", help="DCC from BAM-derived fragment signals.")
    input_a = bam.add_mutually_exclusive_group(required=True)
    input_a.add_argument("--bam-a", nargs="+", help="BAM paths or glob patterns for signal A.")
    input_a.add_argument(
        "--fragments-a", nargs="+",
        help="BED, BED.gz or bigBed fragment files for signal A; only columns 1-3 are required.",
    )
    input_b = bam.add_mutually_exclusive_group(required=True)
    input_b.add_argument("--bam-b", nargs="+", help="BAM paths or glob patterns for signal B.")
    input_b.add_argument(
        "--fragments-b", nargs="+",
        help="BED, BED.gz or bigBed fragment files for signal B; only columns 1-3 are required.",
    )
    bam.add_argument("--length-a", type=int, default=None, help="Exact A fragment length.")
    bam.add_argument("--min-length-a", type=int, default=None, help="Minimum A length.")
    bam.add_argument("--max-length-a", type=int, default=None, help="Maximum A length.")
    bam.add_argument("--length-b", type=int, default=None, help="Exact B fragment length.")
    bam.add_argument("--min-length-b", type=int, default=None, help="Minimum B length.")
    bam.add_argument("--max-length-b", type=int, default=None, help="Maximum B length.")
    bam.add_argument(
        "--position-a",
        choices=["dyad", "left_end", "right_end"],
        default="dyad",
        help="Fragment-derived position for A. Default: dyad.",
    )
    bam.add_argument(
        "--position-b",
        choices=["dyad", "left_end", "right_end"],
        default="dyad",
        help="Fragment-derived position for B. Default: dyad.",
    )
    bam.add_argument("--mapq", type=int, default=0, help="Minimum MAPQ. Default: 0.")
    bam.add_argument(
        "--max-duplicates",
        type=int,
        default=0,
        help=(
            "Maximum identical fragment coordinates retained per BAM. "
            "0 applies no coordinate-count cap. Default: 0."
        ),
    )
    bam.add_argument(
        "--include-duplicate-flag",
        action="store_true",
        help="Include alignments marked duplicate in the BAM flags.",
    )
    bam.add_argument(
        "--no-require-proper-pairs",
        action="store_true",
        help="Do not require the proper-pair BAM flag.",
    )
    add_region_arguments(bam)
    add_shared_dcc_arguments(bam)

    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive signal analysis. "
            "This budget defaults independently to 1."
        ),
    )

    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def resolve_length_spec(
    exact: Optional[int],
    minimum: Optional[int],
    maximum: Optional[int],
    suffix: str,
    parser: argparse.ArgumentParser,
) -> Tuple[int, int]:
    """Resolve exact or ranged BAM fragment-length arguments."""

    if exact is not None:
        if minimum is not None or maximum is not None:
            parser.error(
                f"Use either --length-{suffix} or --min-length-{suffix}/"
                f"--max-length-{suffix}, not both."
            )
        if exact <= 0:
            parser.error(f"--length-{suffix} must be greater than zero.")
        return exact, exact

    if minimum is None or maximum is None:
        parser.error(
            f"Provide --length-{suffix}, or both --min-length-{suffix} and "
            f"--max-length-{suffix}."
        )
    if minimum <= 0 or maximum <= 0:
        parser.error(f"Fragment lengths for signal {suffix.upper()} must be positive.")
    if minimum > maximum:
        parser.error(f"--min-length-{suffix} cannot exceed --max-length-{suffix}.")
    return minimum, maximum


def validate_common_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Validate shared argument combinations."""

    if args.dmax < 0:
        parser.error("--dmax cannot be negative.")
    if args.summary_lag_window < 0:
        parser.error("--summary-lag-window cannot be negative.")
    if args.cpm_scale <= 0:
        parser.error("--cpm-scale must be greater than zero.")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative.")

    if args.min_region_length is None:
        args.min_region_length = max(2, args.dmax + 1)
    elif args.min_region_length < 1:
        parser.error("--min-region-length must be at least 1.")

    if args.regions_bed:
        if args.scope not in {"genome", "combined_chromosomes"}:
            parser.error("--scope is only used with --chrom-sizes.")
    else:
        if args.region_mode != "interval":
            parser.error("--region-mode is only used with --regions-bed.")
        if args.category:
            parser.error("--category is only used with --regions-bed.")


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Validate mode-specific arguments and resolve defaults."""

    validate_common_args(args, parser)

    if args.mode == "bigwig":
        if args.value_limit is not None and args.value_limit <= 0:
            parser.error("--value-limit must be greater than zero.")
        if not 0.0 <= args.sparse_threshold <= 1.0:
            parser.error("--sparse-threshold must be between 0 and 1.")
    elif args.mode == "bam":
        if args.mapq < 0:
            parser.error("--mapq cannot be negative.")
        if args.max_duplicates < 0:
            parser.error("--max-duplicates cannot be negative.")

        args.min_length_a, args.max_length_a = resolve_length_spec(
            args.length_a,
            args.min_length_a,
            args.max_length_a,
            "a",
            parser,
        )
        args.min_length_b, args.max_length_b = resolve_length_spec(
            args.length_b,
            args.min_length_b,
            args.max_length_b,
            "b",
            parser,
        )


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def run(args: argparse.Namespace) -> int:
    """Run the DCC command."""

    parser = build_parser()
    validate_args(args, parser)
    selected_chromosomes = split_comma_values(args.chromosome)
    if args.mode == "bigwig":
        naming_files_a = expand_file_inputs(args.bigwig_a, "bigWig")
        naming_files_b = expand_file_inputs(args.bigwig_b, "bigWig")
    elif args.fragments_a:
        naming_files_a = expand_file_inputs(args.fragments_a, "fragment interval")
        naming_files_b = expand_file_inputs(args.fragments_b, "fragment interval")
    else:
        naming_files_a = expand_file_inputs(args.bam_a, "BAM")
        naming_files_b = expand_file_inputs(args.bam_b, "BAM")
    requested_prefix = args.out_prefix or default_output_prefix(
        naming_files_a,
        naming_files_b,
        args.regions_bed,
        args.scope,
        selected_chromosomes,
        args.label_a,
        args.label_b,
    )
    from nucleosuite.output_naming import parameterized_prefix

    normalization = "raw" if args.no_normalize_dcc else "opportunity"
    if args.normalize_by_signal_totals:
        normalization += "-signal"
    args.out_prefix = str(
        parameterized_prefix(
            requested_prefix,
            (
                ("dmax", args.dmax),
                ("lags", "signed" if args.signed_lags else "absolute"),
                ("norm", normalization),
            ),
        )
    )
    from nucleosuite.parallel import run_region_per_contig
    if not getattr(args, "_per_contig_worker", False) and int(getattr(args, "cores", 1) or 1) > 1:
        return run_region_per_contig("dcc", args, run)

    reporter = ProgressReporter("dcc")
    try:
        category_rules = parse_state_category_rules(args.category)
    except ValueError as exc:
        parser.error(str(exc))

    reporter.stage("Loading analysis regions")
    if args.regions_bed:
        regions = read_regions_bed(
            args.regions_bed,
            args.state_column,
            args.strand_column,
            args.state_name,
            args.region_mode,
            args.extend,
            args.strands,
            selected_chromosomes,
            args.min_region_length,
            category_rules,
        )
    else:
        regions = make_windows_from_chrom_sizes(
            args.chrom_sizes,
            args.scope,
            selected_chromosomes,
            args.window_size,
            args.state_name,
            args.min_region_length,
        )

    for chromosome in dict.fromkeys(region.chrom for region in regions):
        reporter.reading_contig("regions", chromosome)

    blacklist = load_blacklist_unbounded(args.blacklist_bed)
    args._blacklisted_anchor_exclusions = 0
    if blacklist is not None and args.regions_bed:
        retained_regions = []
        for region in regions:
            anchor_start = region.anchor_start if region.anchor_start is not None else region.start
            anchor_end = region.anchor_end if region.anchor_end is not None else region.end
            if blacklist.overlaps(region.chrom, anchor_start, anchor_end):
                args._blacklisted_anchor_exclusions += 1
            else:
                retained_regions.append(region)
        regions = retained_regions
        if not regions:
            parser.error("No regions remained after blacklist filtering")

    regions_by_state = group_regions_by_state(regions)
    os.makedirs(args.output_dir, exist_ok=True)

    if args.mode == "bigwig":
        files_a = expand_file_inputs(args.bigwig_a, "bigWig")
        files_b = expand_file_inputs(args.bigwig_b, "bigWig")
    else:
        if args.fragments_a:
            files_a = expand_file_inputs(args.fragments_a, "fragment interval")
            files_b = expand_file_inputs(args.fragments_b, "fragment interval")
        else:
            files_a = expand_file_inputs(args.bam_a, "BAM")
            files_b = expand_file_inputs(args.bam_b, "BAM")

    reporter.stage(
        f"Prepared {len(regions):,} regions in {len(regions_by_state):,} groups; "
        f"processing {len(files_a):,} A and {len(files_b):,} B input(s)"
    )

    output_prefix = args.out_prefix or default_output_prefix(
        files_a,
        files_b,
        args.regions_bed,
        args.scope,
        selected_chromosomes,
        args.label_a,
        args.label_b,
    )
    from nucleosuite.output_naming import parameterized_prefix

    normalization = "raw" if args.no_normalize_dcc else "opportunity"
    if args.normalize_by_signal_totals:
        normalization += "-signal"
    output_prefix = str(
        parameterized_prefix(
            output_prefix,
            (
                ("dmax", args.dmax),
                ("lags", "signed" if args.signed_lags else "absolute"),
                ("norm", normalization),
            ),
        )
    )

    if not args.quiet:
        print(f"Mode: {args.mode}")
        print(f"Input A files: {len(files_a):,}")
        print(f"Input B files: {len(files_b):,}")
        print(f"Regions: {len(regions):,}")
        print(f"States/groups: {len(regions_by_state):,}")
        print(f"Maximum lag/distance: {args.dmax:,} bp")
        print(
            "Opportunity normalization: "
            + ("off" if args.no_normalize_dcc else "on")
        )

    if args.mode == "bigwig":
        summary_rows = process_bigwig_group(
            files_a,
            files_b,
            regions_by_state,
            output_prefix,
            args.output_dir,
            args,
        )
    else:
        summary_rows = process_bam_group(
            files_a,
            files_b,
            regions_by_state,
            output_prefix,
            args.output_dir,
            args,
        )

    reporter.stage("Writing DCC summary")
    summary_path = os.path.join(
        args.output_dir,
        f"{sanitize_filename(output_prefix)}_DCC_summary.tsv",
    )
    write_summary_tsv(summary_path, summary_rows)

    if not args.quiet:
        print(f"Summary: {summary_path}")

    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
