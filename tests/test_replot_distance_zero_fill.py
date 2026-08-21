"""Regression tests for restoring omitted zero-count distance rows during replotting."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from nucleosuite import replot


def _write_sparse_distance_table(path: Path) -> None:
    fields = [
        "order", "scope", "chromosome", "state", "distance_bp", "count",
        "smoothed_count", "percent", "smoothed_percent", "full_raw_mode_bp",
        "full_smoothed_mode_bp",
    ]
    rows = [
        [1, "combined_chromosomes", ".", "All", 185, 10, 9.0, 66.6666667, 60.0, 185, 185],
        [1, "combined_chromosomes", ".", "All", 190, 5, 6.0, 33.3333333, 40.0, 185, 185],
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        writer.writerows(rows)


def test_replot_distances_can_restore_zero_count_support(tmp_path):
    table = tmp_path / "sample_distances.tsv"
    _write_sparse_distance_table(table)
    headers, rows = replot._read_table(table)
    args = replot.parse_cli_args([
        str(table),
        "--plot-type", "distances",
        "--nrl-mode", "raw",
        "--include-zero-distances",
        "--x-min", "185",
        "--x-max", "190",
    ])
    output = tmp_path / "distance.png"
    _saved, fig = replot._plot_distances(table, headers, rows, args, output, {})
    line = fig.axes[0].lines[0]
    assert np.array_equal(line.get_xdata(), np.arange(185, 191, dtype=float))
    assert np.array_equal(line.get_ydata(), np.asarray([10, 0, 0, 0, 0, 5], dtype=float))


def test_replot_distances_sparse_opt_out_preserves_observed_rows(tmp_path):
    table = tmp_path / "sample_distances.tsv"
    _write_sparse_distance_table(table)
    headers, rows = replot._read_table(table)
    args = replot.parse_cli_args([
        str(table),
        "--plot-type", "distances",
        "--nrl-mode", "raw",
        "--no-include-zero-distances",
    ])
    output = tmp_path / "distance.png"
    _saved, fig = replot._plot_distances(table, headers, rows, args, output, {})
    line = fig.axes[0].lines[0]
    assert np.array_equal(line.get_xdata(), np.asarray([185, 190], dtype=float))
    assert np.array_equal(line.get_ydata(), np.asarray([10, 5], dtype=float))
