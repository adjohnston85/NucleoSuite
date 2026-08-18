"""Parity tests derived from the 2016 cfDNA callPeaks.py implementation."""

from __future__ import annotations

import math

import numpy as np

from nucleosuite.peaks.wps import (
    call_records,
    evaluate_region,
    kircher_continuous_windows,
    kircher_median,
)
from nucleosuite.scoring.wps import kircher_adjusted_signal, kircher_savgol


def _python2_round(value: float) -> int:
    return int(math.floor(value + 0.5)) if value >= 0 else int(math.ceil(value - 0.5))


def _reference_region_median(values) -> float:
    helper = sorted(float(value) for value in values)
    point = len(helper) * 0.5
    rounded = _python2_round(point)
    if rounded == int(point):
        return helper[int(point)]
    return 0.5 * (helper[rounded] + helper[int(point)])


def _reference_continuous_windows(region_pairs):
    result = []
    start = end = maximum = None
    total = 0.0
    for position, value in region_pairs:
        if maximum is None:
            start = end = position
            maximum = float(value)
            total = float(value)
        elif end + 1 == position:
            end = position
            total += float(value)
            maximum = max(maximum, float(value))
        else:
            result.append((total, start, end, maximum))
            start = end = maximum = None
            total = 0.0
    if maximum is not None:
        result.append((total, start, end, maximum))
    return result


def _reference_evaluate(values, start=100, minimum=50, maximum=150, cutoff=5.0):
    end = start + len(values) - 1
    if maximum >= len(values) >= minimum:
        median = _reference_region_median(values)
        pairs = [
            (position, value)
            for position, value in zip(range(start, end + 1), values)
            if value >= median
        ]
        windows = _reference_continuous_windows(pairs)
        windows.sort()
        selected = [windows[-1]]
    elif 3 * maximum >= len(values) >= maximum:
        median = _reference_region_median(values)
        pairs = [
            (position, value)
            for position, value in zip(range(start, end + 1), values)
            if value >= median
        ]
        selected = [
            row
            for row in _reference_continuous_windows(pairs)
            if minimum <= row[2] - row[1] + 1 <= maximum
        ]
    else:
        selected = []
    return [row for row in selected if row[3] > cutoff]


def _reference_savgol(values, window=21, order=2):
    values = np.asarray(values, dtype=float)
    half = (window - 1) // 2
    design = np.asarray(
        [[float(k**i) for i in range(order + 1)] for k in range(-half, half + 1)]
    )
    coefficients = np.linalg.pinv(design)[0][::-1]
    output = np.empty(values.size, dtype=float)
    for index in range(values.size):
        if index < half:
            helper = values[index : index + half + 1]
            stop = min(window - helper.size + 1, half + 1)
            pad = helper[0] - np.abs(helper[1:stop][::-1] - helper[0])
            window_values = np.concatenate((pad, helper))
        elif index >= values.size - half:
            helper = values[max(0, index - half) : index + 1]
            count = min(window - helper.size, half)
            pad = helper[-1] + np.abs(
                helper[-count - 1 : -1][::-1] - helper[-1]
            )
            window_values = np.concatenate((helper, pad))
        else:
            window_values = values[index - half : index + half + 1]
        output[index] = np.convolve(coefficients, window_values, mode="valid")[0]
    return output


def _reference_adjusted(values, baseline_window=1000):
    values = np.asarray(values, dtype=float)
    if values.size < baseline_window:
        return np.zeros(values.size, dtype=float)
    smoothed = _reference_savgol(values)
    half = baseline_window // 2
    medians = np.asarray(
        [
            np.median(values[start : start + baseline_window])
            for start in range(values.size - baseline_window + 1)
        ]
    )
    baseline = np.empty(values.size, dtype=float)
    baseline[: half + 1] = medians[0]
    baseline[half : values.size - half + 1] = medians
    baseline[values.size - half + 1 :] = medians[-1]
    return smoothed - baseline


def test_custom_region_median_matches_python2_source_rule():
    for values in ([1, 2, 3, 4], [1, 2, 3, 4, 5], list(range(1, 60))):
        assert kircher_median(values) == _reference_region_median(values)


def test_continuous_window_gap_reset_matches_source():
    pairs = [(1, 10.0), (2, 11.0), (4, 20.0), (5, 21.0), (6, 22.0)]
    assert kircher_continuous_windows(pairs) == _reference_continuous_windows(pairs)
    # Position 4 is intentionally discarded by the source reset behaviour.
    assert kircher_continuous_windows(pairs) == [
        (21.0, 1, 2, 11.0),
        (43.0, 5, 6, 22.0),
    ]


def test_short_region_selected_window_has_no_secondary_minimum_length():
    # The complete positive region is 60 bp, but its selected above-median
    # maximum-sum window is much shorter than 50 bp. callPeaks.py reports it.
    values = list(np.linspace(1, 30, 30)) + list(np.linspace(30, 1, 30))
    expected = _reference_evaluate(values)
    observed = evaluate_region("chr1", 100, 159, values)
    assert len(expected) == 1
    assert len(observed) == 1
    assert observed[0]["region_start"] == expected[0][1] - 1
    assert observed[0]["region_end"] == expected[0][2]
    assert observed[0]["peak_score"] == expected[0][3]
    assert observed[0]["region_end"] - observed[0]["region_start"] < 50


def test_150_bp_boundary_uses_short_region_branch():
    values = list(np.linspace(1, 75, 75)) + list(np.linspace(75, 1, 75))
    expected = _reference_evaluate(values, start=1)
    observed = evaluate_region("chr1", 1, 150, values)
    assert len(expected) == len(observed) == 1
    assert observed[0]["region_start"] == expected[0][1] - 1
    assert observed[0]["region_end"] == expected[0][2]


def test_extended_region_requires_50_to_150_bp_subwindows():
    values = np.r_[
        np.full(60, 10.0),
        np.full(10, 1.0),
        np.full(70, 12.0),
        np.full(60, 1.0),
    ].tolist()
    expected = _reference_evaluate(values, start=1000)
    observed = evaluate_region("chr1", 1000, 1199, values)
    assert [
        (row["region_start"] + 1, row["region_end"], row["peak_score"])
        for row in observed
    ] == [(row[1], row[2], row[3]) for row in expected]


def test_source_savgol_matches_independent_translation():
    rng = np.random.default_rng(17)
    values = rng.integers(-100, 101, size=1200).astype(float)
    np.testing.assert_allclose(
        kircher_savgol(values),
        _reference_savgol(values),
        rtol=0.0,
        atol=1e-12,
    )


def test_adjusted_wps_matches_independent_source_translation():
    x = np.arange(1400, dtype=float)
    values = 15 * np.sin(x / 17.0) + 7 * np.cos(x / 43.0)
    np.testing.assert_allclose(
        kircher_adjusted_signal(values),
        _reference_adjusted(values),
        rtol=0.0,
        atol=1e-10,
    )


def test_full_region_caller_matches_reference_coordinates():
    signal = np.zeros(300, dtype=float)
    signal[40:100] = np.r_[np.linspace(1, 30, 30), np.linspace(30, 1, 30)]
    records = call_records(
        signal,
        "chr1",
        adjusted_start=1000,
        core_start=1000,
        core_end=1300,
    )
    expected = _reference_evaluate(signal[40:100].tolist(), start=1041)
    assert len(records) == len(expected) == 1
    assert records[0]["region_start"] == expected[0][1] - 1
    assert records[0]["region_end"] == expected[0][2]
