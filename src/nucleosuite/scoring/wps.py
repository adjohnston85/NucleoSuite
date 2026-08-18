"""Kircher-style windowed protection scoring.

The default smoothed, median-centred track follows the preprocessing implemented
by ``callPeaks.py`` from the 2016 cfDNA study: a 21-base, second-order
Savitzky-Golay smoother followed by subtraction of the running median of the raw
WPS in 1,000-base windows.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Iterable

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

WPS_TRACKS = ("wps", "wps_smoothed", "mWPS", "sm_mWPS")


def wps_kernel_kircher_exact(fragment_length: int, protection: int = 120) -> np.ndarray:
    """Construct the effective Kircher WPS kernel for one true fragment length."""
    half = protection // 2
    if fragment_length <= 0:
        return np.array([], dtype=np.int8)

    flank = protection - 1
    middle = fragment_length - protection + 1
    total_length = 2 * flank + max(middle, 0)
    if total_length <= 0:
        return np.array([], dtype=np.int8)

    if middle <= 0:
        return np.full(total_length, -1, dtype=np.int8)

    kernel = np.empty(total_length, dtype=np.int8)
    kernel[:flank] = -1
    kernel[flank : flank + middle] = 1
    kernel[flank + middle :] = -1
    return kernel


def precompute_distributions(
    fragment_lengths: Iterable[int], protection: int = 120
) -> Dict[int, np.ndarray]:
    return {
        int(length): wps_kernel_kircher_exact(int(length), protection)
        for length in fragment_lengths
    }


def rolling_median(values: np.ndarray, window: int = 1000) -> np.ndarray:
    """Return the WPS running median with fixed edge windows.

    For an even 1,000-base window, the median is attached to the right-middle
    position. Positions before the first centred window use the first complete
    window median; positions after the final centred window use the final
    complete window median. When a signal block is shorter than ``window``, the
    original caller does not emit adjusted values, so this function returns NaN.
    """
    values = np.asarray(values, dtype=float)
    length = values.size
    if window < 1:
        raise ValueError("Rolling-median window must be positive")
    if length < window:
        return np.full(length, np.nan, dtype=float)

    windows = sliding_window_view(values, window_shape=window)
    medians = np.median(windows, axis=1)
    output = np.empty(length, dtype=float)
    half = window // 2

    # Window starting at zero is used for all positions through the first
    # right-middle position, matching the streaming implementation.
    output[: half + 1] = medians[0]
    if medians.size > 1:
        output[half : length - half + 1] = medians
    # The final complete window is retained while the stream is flushed.
    output[length - half + 1 :] = medians[-1]
    return output


@lru_cache(maxsize=None)
def _kircher_savgol_coefficients(window: int, order: int) -> np.ndarray:
    """Return the convolution coefficients used by the 2016 source script."""
    if window < 3 or window % 2 == 0:
        raise ValueError("Savitzky-Golay window must be an odd integer >= 3")
    if order < 0 or order >= window:
        raise ValueError("Savitzky-Golay order must satisfy 0 <= order < window")
    half = (window - 1) // 2
    design = np.asarray(
        [[float(k**i) for i in range(order + 1)] for k in range(-half, half + 1)],
        dtype=float,
    )
    # callPeaks.py obtains the first pseudoinverse row and reverses it before
    # passing it to numpy.convolve.
    return np.asarray(np.linalg.pinv(design)[0], dtype=float)[::-1]


def kircher_savgol(values: np.ndarray, window: int = 21, order: int = 2) -> np.ndarray:
    """Apply the exact Savitzky-Golay edge treatment from ``callPeaks.py``.

    Interior positions are evaluated by convolution. At chromosome/block ends,
    the source script constructs asymmetric reflected values relative to the
    current endpoint; those first and last positions are reproduced explicitly.
    """
    values = np.asarray(values, dtype=float)
    length = values.size
    if length == 0 or window <= 0:
        return values.copy()
    if window < 3 or window % 2 == 0:
        raise ValueError("Savitzky-Golay window must be odd, or 0 to disable")
    if order < 0 or order >= window:
        raise ValueError("Savitzky-Golay order must be smaller than its window")
    if length < window:
        # The original peak caller never processes a block shorter than its
        # 1,000-base median window. Returning the raw values keeps track output
        # well-defined while adjusted peak input is suppressed separately.
        return values.copy()

    half = (window - 1) // 2
    coefficients = _kircher_savgol_coefficients(window, order)
    output = np.empty(length, dtype=float)

    # Full windows. This is the same np.convolve call used by the source script,
    # applied to every complete window at once.
    output[half : length - half] = np.convolve(
        coefficients, values, mode="valid"
    )

    # Left edge: source smoother(values, ctype=True).
    for index in range(half):
        helper = values[index : index + half + 1]
        stop = min(window - helper.size + 1, half + 1)
        reflected = helper[0] - np.abs(helper[1:stop][::-1] - helper[0])
        padded = np.concatenate((reflected, helper))
        output[index] = np.convolve(coefficients, padded, mode="valid")[0]

    # Right edge: source smoother(values, ctype=False).
    for index in range(length - half, length):
        helper = values[max(0, index - half) : index + 1]
        count = min(window - helper.size, half)
        reflected = helper[-1] + np.abs(
            helper[-count - 1 : -1][::-1] - helper[-1]
        )
        padded = np.concatenate((helper, reflected))
        output[index] = np.convolve(coefficients, padded, mode="valid")[0]

    return output


def kircher_adjusted_signal(
    wps_array: np.ndarray,
    baseline_window: int = 1000,
    smooth_window: int = 21,
    smooth_order: int = 2,
) -> np.ndarray:
    """Return the adjusted L-WPS signal used for peak calling.

    A block shorter than the running-median window yields an all-zero adjusted
    signal because the original streaming caller produces no adjusted positions
    for such a block.
    """
    raw = np.asarray(wps_array, dtype=float)
    if raw.size < baseline_window:
        return np.zeros(raw.size, dtype=float)
    smoothed = kircher_savgol(raw, smooth_window, smooth_order)
    baseline = rolling_median(raw, baseline_window)
    adjusted = smoothed - baseline
    if raw.size == baseline_window:
        # With exactly one full buffer, callPeaks.py starts its EOF flush at
        # buffer index 1 and never emits the first coordinate.
        adjusted[0] = 0.0
    return adjusted


def new_array(reference_length: int) -> np.ndarray:
    return np.zeros(reference_length, dtype=np.float64)


def add_fragment(
    wps_array: np.ndarray,
    fragment_start: int,
    fragment_end: int,
    window_start: int,
    protection: int,
    distributions: Dict[int, np.ndarray],
) -> None:
    fragment_length = fragment_end - fragment_start
    kernel = distributions.get(fragment_length)
    if kernel is None or kernel.size == 0:
        return

    half = protection // 2
    kernel_genome_start = fragment_start - half + 1
    kernel_start = kernel_genome_start - window_start
    kernel_end = kernel_start + kernel.size

    array_start = max(kernel_start, 0)
    array_end = min(kernel_end, wps_array.size)
    if array_end <= array_start:
        return
    source_start = array_start - kernel_start
    source_end = source_start + (array_end - array_start)
    wps_array[array_start:array_end] += kernel[source_start:source_end]


def to_scores(
    wps_array: np.ndarray,
    contig: str,
    start: int,
    baseline_window: int = 1000,
    smooth_window: int = 21,
    smooth_order: int = 2,
):
    raw = np.asarray(wps_array, dtype=float)
    if smooth_window > 0:
        smoothed = kircher_savgol(raw, smooth_window, smooth_order)
    else:
        smoothed = raw.copy()

    baseline = rolling_median(raw, baseline_window)
    if np.isnan(baseline).all():
        # Preserve useful raw/smoothed outputs for very short contigs while
        # preventing calls that the original 1-kb caller could not make.
        baseline_for_tracks = np.zeros(raw.size, dtype=float)
        median_centred = np.zeros(raw.size, dtype=float)
        smoothed_median_centred = np.zeros(raw.size, dtype=float)
    else:
        baseline_for_tracks = baseline
        median_centred = raw - baseline_for_tracks
        smoothed_median_centred = smoothed - baseline_for_tracks
        if raw.size == baseline_window:
            median_centred[0] = 0.0
            smoothed_median_centred[0] = 0.0

    return {
        "wps": [(contig, start, raw)],
        "wps_smoothed": [(contig, start, smoothed)],
        "mWPS": [(contig, start, median_centred)],
        "sm_mWPS": [(contig, start, smoothed_median_centred)],
    }
