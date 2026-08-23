"""Tests for NRL and periodicity estimation."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from nucleosuite.nrl import (
    Peak,
    Regression,
    call_resolution_peaks,
    create_regression_plot,
    local_maximum_indices,
    main,
    moving_average_by_distance,
    regress_peak_distances,
    resolution_smoothing_windows,
    retain_separated_peaks,
    snap_smoothing_window,
)


def test_regression_plot_defaults_to_open_circles_dotted_fit_and_square(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda figure=None: None)
    peaks = [
        Peak(index=i, distance=float(185 * (i + 1)), raw_value=1.0, smoothed_value=1.0)
        for i in range(3)
    ]
    regression = Regression(3, 185.0, 0.0, 1.0, 0.0, 185.0)
    create_regression_plot(tmp_path / "nrl.png", peaks, regression, "NRL", 100)
    figure = plt.gcf()
    axis = figure.axes[0]
    assert tuple(figure.get_size_inches()) == (6.5, 6.5)
    assert axis.lines[0].get_linestyle() == ":"
    assert axis.collections[0].get_facecolors().size == 0
    original_close(figure)


def test_local_maxima_require_persistent_rise_and_fall_and_support_plateaus():
    values = np.asarray([0, 1, 3, 3, 3, 2, 0, 2, 1], dtype=float)
    # The broad plateau has two-base rise and fall spans. The final one-base blip
    # is rejected by the default 2/2 shape requirement.
    assert local_maximum_indices(values) == [3]
    assert local_maximum_indices(values, min_rise_bp=1, min_fall_bp=1) == [3, 7]


def test_profile_peak_persistence_examples():
    assert local_maximum_indices(np.asarray([4, 4, 5, 4, 4], dtype=float)) == []
    assert local_maximum_indices(np.asarray([4, 5, 6, 4, 4], dtype=float)) == []
    assert local_maximum_indices(np.asarray([4, 5, 6, 5, 4], dtype=float)) == [2]
    assert local_maximum_indices(np.asarray([4, 5, 6, 6, 5, 4], dtype=float)) == [2]
    # A fall followed immediately by another rise does not qualify either summit.
    assert local_maximum_indices(np.asarray([4, 5, 6, 5, 6, 5, 4], dtype=float)) == []


def test_profile_peak_persistence_uses_distance_span():
    distances = np.asarray([0, 5, 10], dtype=float)
    values = np.asarray([0, 1, 0], dtype=float)
    assert local_maximum_indices(distances, values, 2, 2) == [1]


def test_peak_separation_keeps_highest_competing_peak():
    distances = np.asarray([0, 10, 20, 30, 40, 50, 60], dtype=float)
    values = np.asarray([0, 5, 0, 9, 0, 4, 0], dtype=float)
    peaks = retain_separated_peaks(
        [1, 3, 5], distances, values, values, min_separation=25
    )
    assert [peak.distance for peak in peaks] == [30.0]


def test_regression_slope_is_peak_spacing():
    distances = np.arange(0, 501, dtype=float)
    values = np.zeros_like(distances)
    indices = [100, 150, 200, 250, 300]
    for index in indices:
        values[index] = 1.0
    peaks = retain_separated_peaks(indices, distances, values, values, 0)
    regression = regress_peak_distances(peaks)
    assert regression.slope == pytest.approx(50.0)
    assert regression.r_squared == 1.0


def test_resolution_windows_snap_down_to_10n_plus_1():
    assert snap_smoothing_window(0) == 1
    assert snap_smoothing_window(10.99) == 1
    assert snap_smoothing_window(11) == 11
    assert snap_smoothing_window(20) == 11
    assert snap_smoothing_window(21) == 21
    assert snap_smoothing_window(53.333) == 51
    assert snap_smoothing_window(60) == 51
    assert snap_smoothing_window(61) == 61
    assert resolution_smoothing_windows(160) == (51, 21)
    assert resolution_smoothing_windows(180) == (51, 21)
    assert resolution_smoothing_windows(200) == (61, 31)


def test_distance_based_smoothing_uses_bp_not_row_count():
    distances = np.asarray([0, 10, 20, 30, 40], dtype=float)
    values = np.asarray([0, 0, 10, 0, 0], dtype=float)
    smoothed = moving_average_by_distance(distances, values, 21)
    assert smoothed[2] == pytest.approx(10 / 3)
    assert smoothed[0] == pytest.approx(0)


def test_long_range_resolution_detects_decaying_wiggly_peaks_and_refines_them():
    distances = np.arange(200, 1201, dtype=float)
    values = np.zeros_like(distances)
    centres = [370, 555, 740, 925, 1110]
    amplitudes = [10.0, 8.0, 6.0, 4.0, 2.5]
    for centre, amplitude in zip(centres, amplitudes):
        values += amplitude * np.exp(-0.5 * ((distances - centre) / 24.0) ** 2)
        values += 0.30 * amplitude * np.exp(-0.5 * ((distances - (centre - 8)) / 3.0) ** 2)
        values += 0.25 * amplitude * np.exp(-0.5 * ((distances - (centre + 9)) / 3.0) ** 2)

    detection_window, local_window = resolution_smoothing_windows(160)
    detection = moving_average_by_distance(distances, values, detection_window)
    local = moving_average_by_distance(distances, values, local_window)
    peaks = call_resolution_peaks(distances, values, local, detection, 160)
    assert detection_window == 51
    assert local_window == 21
    assert [peak.distance for peak in peaks] == [float(x) for x in centres]
    assert all(b.distance - a.distance >= 160 for a, b in zip(peaks, peaks[1:]))


def write_test_profile(path: Path, value_header: str) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Distance", value_header])
        for distance in range(1, 301):
            # Clear 50-bp periodicity with rise-then-fall triangular peaks.
            phase = distance % 50
            value = 25 - abs(phase - 25)
            writer.writerow([distance, value])


def test_nrl_command_writes_tables_and_plots_for_dac(tmp_path):
    input_path = tmp_path / "sample_DAC.tsv"
    write_test_profile(input_path, "DAC Value")
    prefix = tmp_path / "sample_nrl"

    status = main(
        [
            str(input_path),
            "--min-distance", "1",
            "--max-distance", "300",
            "--peak-resolution", "0",
            "--output-prefix", str(prefix),
            "--quiet",
        ]
    )
    assert status == 0
    prefix = Path(f"{prefix}_peakres0_min1_max300_skipfirst0")
    for suffix in (
        "_profile.tsv",
        "_peaks.tsv",
        "_regression.tsv",
        "_profile.png",
        "_regression.png",
    ):
        assert Path(f"{prefix}{suffix}").is_file()

    with Path(f"{prefix}_regression.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["peak_resolution_bp"] == "0"
    assert rows[0]["detection_smoothing_window"] == "1"
    assert rows[0]["local_max_smoothing_window"] == "1"
    assert abs(float(rows[0]["slope_bp_per_peak"]) - 50.0) < 1e-9

    with Path(f"{prefix}_peaks.tsv").open() as handle:
        peak_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert peak_rows
    assert all("detection_peak_distance" in row for row in peak_rows)


def test_nrl_command_auto_detects_dcc_value(tmp_path):
    input_path = tmp_path / "sample_DCC.tsv"
    write_test_profile(input_path, "DCC Value")
    prefix = tmp_path / "sample_dcc_periodicity"
    status = main(
        [
            str(input_path),
            "--min-distance", "1",
            "--max-distance", "300",
            "--peak-resolution", "0",
            "--output-prefix", str(prefix),
            "--quiet",
        ]
    )
    assert status == 0
    prefix = Path(f"{prefix}_peakres0_min1_max300_skipfirst0")
    assert Path(f"{prefix}_regression.tsv").is_file()


def test_nrl_skip_first_peak_keeps_called_peak_numbering(tmp_path):
    input_path = tmp_path / "skip_DAC.tsv"
    write_test_profile(input_path, "DAC Value")
    prefix = tmp_path / "skip_nrl"
    status = main([
        str(input_path), "--min-distance", "1", "--max-distance", "300",
        "--peak-resolution", "0", "--skip-first-peaks", "1",
        "--output-prefix", str(prefix), "--quiet",
    ])
    assert status == 0
    out = Path(f"{prefix}_peakres0_min1_max300_skipfirst1")
    import csv
    with open(f"{out}_peaks.tsv", encoding="utf-8") as handle:
        peaks = list(csv.DictReader(handle, delimiter="\t"))
    with open(f"{out}_regression.tsv", encoding="utf-8") as handle:
        regression = next(csv.DictReader(handle, delimiter="\t"))
    assert len(peaks) >= 2
    assert [int(row["peak_number"]) for row in peaks][:2] == [1, 2]
    assert int(regression["skip_first_peaks"]) == 1
    assert int(regression["called_peak_count"]) == len(peaks)
    assert int(regression["regression_peak_count"]) == len(peaks) - 1
