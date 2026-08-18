"""WPS peak calling reproduced from the 2016 ``callPeaks.py`` source."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np


def python2_round(value: float) -> int:
    """Reproduce Python 2 round-to-nearest with halves away from zero."""
    value = float(value)
    if value >= 0:
        return int(math.floor(value + 0.5))
    return int(math.ceil(value - 0.5))


def kircher_median(values: Iterable[float]) -> float:
    """Reproduce the custom ``median``/``quantiles`` function in callPeaks.py.

    This is intentionally not NumPy's conventional median. The source uses its
    own quantile rule together with Python 2 rounding, so even-length regions use
    the upper middle value and odd-length regions average the two values around
    the half-index.
    """
    helper = sorted(float(value) for value in values)
    length = len(helper)
    if not length:
        raise ValueError("Cannot calculate a median from no values")
    point = length * 0.5
    rounded = python2_round(point)
    if rounded == int(point):
        return float(helper[int(point)])
    return float((helper[rounded] + helper[int(point)]) * 0.5)


def kircher_continuous_windows(region_pairs):
    """Reproduce ``continousWindows`` including its gap-reset behaviour.

    When a gap is encountered, the source appends the completed window and
    resets without starting a new window from the current pair. The pair at the
    gap is therefore intentionally discarded.
    """
    result = []
    current_start = None
    current_end = None
    current_max = None
    current_sum = 0.0
    for position, value in region_pairs:
        if current_max is None:
            current_start = current_end = position
            current_max = float(value)
            current_sum = float(value)
        elif current_end + 1 == position:
            current_end = position
            current_sum += float(value)
            if current_max < value:
                current_max = float(value)
        else:
            result.append((current_sum, current_start, current_end, current_max))
            current_start = current_end = None
            current_max = None
            current_sum = 0.0
    if current_max is not None:
        result.append((current_sum, current_start, current_end, current_max))
    return result


def evaluate_region(
    chrom: str,
    start_1based: int,
    end_1based: int,
    values,
    min_length: int = 50,
    max_length: int = 150,
    max_region_length: int | None = None,
    score_cutoff: float = 5.0,
):
    """Evaluate one merged above-zero region exactly as callPeaks.py does."""
    if start_1based is None or end_1based is None or not values:
        return []

    if max_region_length is None:
        max_region_length = 3 * max_length

    region_length = len(values)
    if not (min_length <= region_length <= max_region_length):
        return []

    median = kircher_median(values)
    above_median = [
        (position, float(value))
        for position, value in zip(range(start_1based, end_1based + 1), values)
        if value >= median
    ]
    windows = kircher_continuous_windows(above_median)
    if not windows:
        return []

    # For complete regions of 50-150 bp, the source chooses one maximum-sum
    # above-median window. It does NOT impose the 50-bp minimum on that selected
    # subwindow. This distinction is essential for reproducing the call count.
    if region_length <= max_length:
        windows.sort()
        candidate_windows = [windows[-1]]
        require_subwindow_length = False
    else:
        candidate_windows = windows
        require_subwindow_length = True

    records = []
    for score_sum, start, end, maximum in candidate_windows:
        segment_length = end - start + 1
        if require_subwindow_length and not (
            min_length <= segment_length <= max_length
        ):
            continue
        if maximum <= score_cutoff:
            continue

        centre_1based = start + python2_round((end - start) * 0.5)
        centre_0based = centre_1based - 1
        records.append(
            {
                "chrom": chrom,
                "region_start": start - 1,
                "region_end": end,
                "region_centre": centre_0based,
                "raw_peak": centre_0based,
                "peak_score": float(maximum),
                "window_sum": float(score_sum),
                "source_start_1based": int(start),
                "source_end_1based": int(end),
            }
        )
    return records


def call_records(
    scores: np.ndarray,
    chrom: str,
    adjusted_start: int,
    core_start: int,
    core_end: int,
    merge_gap_bp: int = 5,
    min_length: int = 50,
    max_length: int = 150,
    max_region_length: int = 450,
    score_cutoff: float = 5.0,
    flip_scores: bool = False,
):
    """Call WPS peaks and assign each to one core chunk.

    The source joins positive positions separated by no more than five bases and
    inserts zero values for the intervening coordinates. NucleoSuite's chunk
    ownership rule retains a call only where its centre belongs to the unpadded
    core, preventing duplicates without changing the peak evaluator.
    """
    records = []
    current_start = None
    current_end = None
    current_values = []

    def close_region():
        nonlocal current_start, current_end, current_values
        if current_start is None or not current_values:
            return
        called = evaluate_region(
            chrom=chrom,
            start_1based=current_start,
            end_1based=current_end,
            values=current_values,
            min_length=min_length,
            max_length=max_length,
            max_region_length=max_region_length,
            score_cutoff=score_cutoff,
        )
        for record in called:
            if core_start <= record["region_centre"] < core_end:
                if flip_scores:
                    record["peak_score"] *= -1.0
                records.append(record)
        current_start = None
        current_end = None
        current_values = []

    for index, value in enumerate(np.asarray(scores, dtype=float)):
        position_1based = adjusted_start + index + 1
        if value > 0:
            if current_end is not None and position_1based <= current_end + merge_gap_bp:
                while current_end + 1 < position_1based:
                    current_end += 1
                    current_values.append(0.0)
                current_values.append(float(value))
                current_end = position_1based
            else:
                close_region()
                current_start = current_end = position_1based
                current_values = [float(value)]
    close_region()
    return records
