"""Shared helpers for reproducible equal-count percentile rank bins."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RankBin:
    """One percentile-labelled slice of a globally score-ranked peak list."""

    percentile_lower: float
    percentile_upper: float
    rank_start: int
    rank_stop: int


def equal_rank_bins(peak_count: int, bin_size: float) -> list[RankBin]:
    """Return equal-count rank slices whose percentile widths divide 100."""
    if not math.isfinite(bin_size) or bin_size <= 0.0 or bin_size > 100.0:
        raise ValueError("--pct-bin-size must be greater than 0 and at most 100")

    bin_count_float = 100.0 / float(bin_size)
    bin_count = int(round(bin_count_float))
    if not math.isclose(
        bin_count * float(bin_size),
        100.0,
        rel_tol=0.0,
        abs_tol=1e-9,
    ):
        raise ValueError("--pct-bin-size must divide 100 into equal percentile ranges")
    if peak_count < bin_count:
        raise ValueError(
            f"--pct-bin-size {bin_size:g} creates {bin_count:,} groups, but only "
            f"{peak_count:,} peaks are available"
        )

    quotient, remainder = divmod(int(peak_count), bin_count)
    bins: list[RankBin] = []
    rank_start = 0
    for index in range(bin_count):
        group_size = quotient + (1 if index < remainder else 0)
        rank_stop = rank_start + group_size
        lower = round(index * float(bin_size), 12)
        upper = 100.0 if index == bin_count - 1 else round(
            (index + 1) * float(bin_size), 12
        )
        bins.append(RankBin(lower, upper, rank_start, rank_stop))
        rank_start = rank_stop
    return bins


def rank_bins_from_boundaries(
    peak_count: int,
    boundaries: list[float] | tuple[float, ...],
) -> list[RankBin]:
    """Return exact rank slices for arbitrary global percentage boundaries.

    Counts are allocated across all segments implied by the requested
    boundaries plus the implicit 0 and 100 edges. Integer remainders are
    assigned by the largest-remainder method, with earlier segments winning
    exact ties. This preserves the requested global percentages as closely as
    possible while guaranteeing non-overlapping integer rank slices.
    """
    if peak_count < 1:
        raise ValueError("Rank bins require at least one peak")

    requested = [float(value) for value in boundaries]
    if len(requested) < 2:
        raise ValueError("Rank bins require at least two boundaries")
    if any(not math.isfinite(value) or value < 0.0 or value > 100.0 for value in requested):
        raise ValueError("Percentage-bin boundaries must be between 0 and 100")
    if any(upper <= lower for lower, upper in zip(requested, requested[1:])):
        raise ValueError("Percentage-bin boundaries must be strictly increasing")

    all_boundaries = sorted({0.0, 100.0, *requested})
    widths = [upper - lower for lower, upper in zip(all_boundaries, all_boundaries[1:])]
    exact_counts = [peak_count * width / 100.0 for width in widths]
    counts = [int(math.floor(value)) for value in exact_counts]
    remainder = peak_count - sum(counts)
    if remainder:
        ranked_remainders = sorted(
            range(len(counts)),
            key=lambda index: (-(exact_counts[index] - counts[index]), index),
        )
        for index in ranked_remainders[:remainder]:
            counts[index] += 1

    rank_by_boundary: dict[float, int] = {all_boundaries[0]: 0}
    cursor = 0
    for upper, count in zip(all_boundaries[1:], counts):
        cursor += count
        rank_by_boundary[upper] = cursor

    bins: list[RankBin] = []
    for lower, upper in zip(requested, requested[1:]):
        rank_start = rank_by_boundary[lower]
        rank_stop = rank_by_boundary[upper]
        if rank_stop <= rank_start:
            raise ValueError(
                f"Percentage range {lower:g}-{upper:g} contains no peaks at "
                f"the current total of {peak_count:,}; use wider bins"
            )
        bins.append(RankBin(lower, upper, rank_start, rank_stop))
    return bins


def randomized_score_order(scores: np.ndarray, seed: int) -> np.ndarray:
    """Shuffle globally, then stably sort by score to randomize only score ties."""
    values = np.asarray(scores)
    random_order = np.random.default_rng(seed).permutation(values.size)
    stable_score_order = np.argsort(values[random_order], kind="stable")
    return random_order[stable_score_order]
