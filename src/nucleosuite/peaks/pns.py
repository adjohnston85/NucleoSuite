"""PNS positive-region and breakpoint peak calling."""

from __future__ import annotations

from typing import Iterable, Iterator, Optional

import numpy as np

from nucleosuite.peaks.common import write_bed8_records


def find_peaks_and_regions(
    scores: np.ndarray,
    genomic_start: int,
    min_length: int = 50,
    max_nonpositive_run: int = 0,
):
    """Identify positive regions, maxima and inter-region minima.

    A positive region may bridge at most ``max_nonpositive_run`` consecutive
    values that are zero or negative. Coordinates returned are absolute,
    zero-based genomic coordinates; regions are half-open.
    """
    if min_length < 1:
        raise ValueError("min_length must be positive")
    if max_nonpositive_run < 0:
        raise ValueError("max_nonpositive_run must be non-negative")

    values = np.asarray(scores, dtype=float)
    regions: list[tuple[int, int]] = []
    region_start: Optional[int] = None
    last_positive: Optional[int] = None
    nonpositive_run = 0

    for index, score in enumerate(values):
        if score > 0:
            if region_start is None:
                region_start = index
            last_positive = index
            nonpositive_run = 0
            continue

        if region_start is None:
            continue
        nonpositive_run += 1
        if nonpositive_run > max_nonpositive_run:
            assert last_positive is not None
            if last_positive - region_start + 1 >= min_length:
                regions.append((region_start, last_positive + 1))
            region_start = None
            last_positive = None
            nonpositive_run = 0

    if region_start is not None and last_positive is not None:
        if last_positive - region_start + 1 >= min_length:
            regions.append((region_start, last_positive + 1))

    positive_peaks: list[int] = []
    positive_scores: list[float] = []
    for start, end in regions:
        local = values[start:end]
        peak_index = start + int(np.argmax(local))
        positive_peaks.append(peak_index + genomic_start)
        positive_scores.append(float(values[peak_index]))

    negative_peaks: list[int] = []
    negative_scores: list[float] = []
    for previous, following in zip(regions, regions[1:]):
        inter_start = previous[1]
        inter_end = following[0]
        if inter_end <= inter_start:
            continue
        local = values[inter_start:inter_end]
        minimum_index = inter_start + int(np.argmin(local))
        negative_peaks.append(minimum_index + genomic_start)
        negative_scores.append(float(values[minimum_index]))

    absolute_regions = [
        (start + genomic_start, end + genomic_start) for start, end in regions
    ]
    region_centres = [(start + end) // 2 for start, end in absolute_regions]
    return {
        "positive_peaks": positive_peaks,
        "positive_peak_scores": positive_scores,
        "negative_peaks": negative_peaks,
        "negative_peak_scores": negative_scores,
        "region_centres": region_centres,
        "regions": absolute_regions,
    }


def iter_peak_records(
    call,
    chrom: str,
    core_start: int,
    core_end: int,
    flip_scores: bool = False,
    coverage_scores: Optional[np.ndarray] = None,
    coverage_start: Optional[int] = None,
) -> Iterator[dict]:
    """Yield records whose region centres belong to the requested core interval."""
    negative_positions = call["negative_peaks"]
    negative_scores = call["negative_peak_scores"]

    for index, centre in enumerate(call["region_centres"]):
        if not core_start <= centre < core_end:
            continue

        upstream_index = None
        downstream_index = None
        for candidate, position in enumerate(negative_positions):
            if position < centre:
                upstream_index = candidate
            elif position > centre:
                downstream_index = candidate
                break

        peak_score = float(call["positive_peak_scores"][index])
        upstream_position = centre
        downstream_position = centre
        upstream_score = peak_score
        downstream_score = peak_score
        if upstream_index is not None:
            upstream_position = int(negative_positions[upstream_index])
            upstream_score = float(negative_scores[upstream_index])
        if downstream_index is not None:
            downstream_position = int(negative_positions[downstream_index])
            downstream_score = float(negative_scores[downstream_index])

        if flip_scores:
            peak_score *= -1.0
            upstream_score *= -1.0
            downstream_score *= -1.0

        maximum_coverage = 0
        maximum_position = 0
        region_start, region_end = call["regions"][index]
        if coverage_scores is not None and coverage_start is not None:
            local_start = max(0, region_start - coverage_start)
            local_end = min(len(coverage_scores), region_end - coverage_start)
            if local_end > local_start:
                local = coverage_scores[local_start:local_end]
                maximum_offset = int(np.argmax(local))
                maximum_coverage = int(local[maximum_offset])
                maximum_position = coverage_start + local_start + maximum_offset

        yield {
            "chrom": chrom,
            "region_start": int(region_start),
            "region_end": int(region_end),
            "region_centre": int(centre),
            "raw_peak": int(call["positive_peaks"][index]),
            "upstream_negative_peak": upstream_position,
            "downstream_negative_peak": downstream_position,
            "upstream_score": upstream_score,
            "downstream_score": downstream_score,
            "peak_score": peak_score,
            "prominence": peak_score,
            "max_coverage": maximum_coverage,
            "max_position": maximum_position,
        }


def call_records(
    scores: np.ndarray,
    chrom: str,
    adjusted_start: int,
    core_start: int,
    core_end: int,
    min_length: int,
    max_nonpositive_run: int,
    flip_scores: bool = False,
    coverage_scores: Optional[np.ndarray] = None,
) -> list[dict]:
    call = find_peaks_and_regions(
        scores=scores,
        genomic_start=adjusted_start,
        min_length=min_length,
        max_nonpositive_run=max_nonpositive_run,
    )
    return list(
        iter_peak_records(
            call=call,
            chrom=chrom,
            core_start=core_start,
            core_end=core_end,
            flip_scores=flip_scores,
            coverage_scores=coverage_scores,
            coverage_start=adjusted_start,
        )
    )



def filter_records_by_coverage(
    records: Iterable[dict],
    coverage_scores: np.ndarray,
    coverage_start: int,
    threshold: float,
) -> tuple[list[dict], int]:
    """Retain nucleosome records with coverage >= threshold at BED column 7.

    PNS BED8 output writes ``region_centre`` to column 7.  The direct
    peak-coverage filter evaluates the internally calculated fragment coverage
    at that same genomic position.
    """
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("coverage threshold must be a finite value of 0 or greater")

    values = np.asarray(coverage_scores)
    retained: list[dict] = []
    filtered = 0
    for record in records:
        position = int(record["region_centre"])
        index = position - int(coverage_start)
        coverage = 0.0
        if 0 <= index < len(values):
            value = float(values[index])
            if np.isfinite(value):
                coverage = value
        if coverage >= threshold:
            retained.append(record)
        else:
            filtered += 1
    return retained, filtered

def write_records(
    path: str,
    records: Iterable[dict],
    label: str,
    score_scale: float,
    mode: str,
) -> None:
    """Write PNS calls as BED8 with six-decimal floating-point scores."""
    write_bed8_records(
        path,
        records,
        label,
        score_scale,
        mode,
        integer_scores=False,
    )
