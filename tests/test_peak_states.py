"""Tests for chromatin-state peak abundance and percentile enrichment."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from nucleosuite import peak_states


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    states = tmp_path / "states.bed"
    states.write_text(
        "chr1\t0\t100\tA\t0\t.\t0\t100\t255,0,0\n"
        "chr1\t100\t300\tB\t0\t.\t100\t300\t0,0,255\n",
        encoding="utf-8",
    )
    peaks = tmp_path / "peaks.bed"
    rows = [
        (10, 1),
        (110, 2),
        (120, 3),
        (130, 4),
        (140, 5),
        (150, 6),
        (160, 7),
        (350, 8),
        (20, 9),
        (30, 10),
    ]
    peaks.write_text(
        "".join(
            f"chr1\t{position}\t{position + 1}\tp{index}\t{score}\n"
            for index, (position, score) in enumerate(rows, start=1)
        ),
        encoding="utf-8",
    )
    return peaks, states


def test_peak_states_reports_coverage_density_enrichment_and_unassigned(tmp_path):
    peaks, states = _write_inputs(tmp_path)
    prefix = tmp_path / "analysis"

    assert peak_states.main(
        [
            str(peaks),
            "--state-bed", str(states),
            "--pct-values", "0,50,90",
            "--output-prefix", str(prefix),
        ]
    ) == 0
    prefix = Path(f"{prefix}_groupsvalues-0-50-90_tiessplit_overlapfirst")

    coverage = list(
        csv.DictReader(
            Path(f"{prefix}_state_coverage.tsv").open(), delimiter="\t"
        )
    )
    assert [(row["state"], int(row["coverage_bp"]), row["rgb"]) for row in coverage] == [
        ("A", 100, "255,0,0"),
        ("B", 200, "0,0,255"),
    ]

    summary = list(
        csv.DictReader(
            Path(f"{prefix}_peak_state_threshold_summary.tsv").open(),
            delimiter="\t",
        )
    )
    assert [float(row["percentile_threshold"]) for row in summary] == [0, 50, 90]
    assert [int(row["unassigned_peak_count"]) for row in summary] == [1, 1, 0]
    assert [int(row["assigned_peak_count"]) for row in summary] == [9, 4, 1]

    enrichment = list(
        csv.DictReader(
            Path(f"{prefix}_peak_state_enrichment.tsv").open(), delimiter="\t"
        )
    )
    at_zero = {row["state"]: row for row in enrichment if row["percentile_threshold"] == "0"}
    assert int(at_zero["A"]["peak_count"]) == 3
    assert int(at_zero["B"]["peak_count"]) == 6
    assert math.isclose(float(at_zero["A"]["peaks_per_mb"]), 30_000.0)
    assert math.isclose(float(at_zero["B"]["peaks_per_mb"]), 30_000.0)
    assert math.isclose(float(at_zero["A"]["enrichment_vs_state_coverage"]), 1.0)
    assert math.isclose(float(at_zero["B"]["enrichment_vs_state_coverage"]), 1.0)

    plot = Path(f"{prefix}_peak_state_percentages_stacked.png")
    assert plot.is_file() and plot.stat().st_size > 0


def test_explicit_peak_state_percentiles_accept_comma_and_space_tokens():
    values = peak_states.explicit_percentile_values(["10,20", "50", "90,99", "20"])
    assert values == [10, 20, 50, 90, 99]


def test_peak_state_percentile_bins_are_non_cumulative_and_keep_ties_together():
    scores = np.asarray([0, 1, 1, 2, 3, 4, 5, 6, 7, 8], dtype=float)
    state_codes = np.zeros(scores.size, dtype=np.int32)
    selections = peak_states.choose_thresholds(
        scores,
        score_percentile=0,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_values=["10,20,30"],
        pct_bins=True,
        bin_tie_mode="keep",
    )

    assert [selection.label for selection in selections] == ["10-20", "20-30"]
    snapshots = peak_states.count_threshold_states(
        scores,
        state_codes,
        selections,
        state_count=1,
    )
    assert [snapshot.retained_peaks for snapshot in snapshots] == [0, 2]
    assert [int(snapshot.counts[0]) for snapshot in snapshots] == [0, 2]


def test_peak_state_bin_outputs_report_ranges(tmp_path):
    peaks, states = _write_inputs(tmp_path)
    prefix = tmp_path / "bins"

    assert peak_states.main(
        [
            str(peaks),
            "--state-bed", str(states),
            "--pct-values", "10,20,30",
            "--pct-bins",
            "--bin-tie-mode", "keep",
            "--output-prefix", str(prefix),
        ]
    ) == 0
    prefix = Path(f"{prefix}_groupsbins-10-20-30_tieskeep_overlapfirst")

    rows = list(
        csv.DictReader(
            Path(f"{prefix}_peak_state_threshold_summary.tsv").open(),
            delimiter="\t",
        )
    )
    assert [row["percentile_mode"] for row in rows] == ["bin", "bin"]
    assert [row["percentile_range"] for row in rows] == ["10-20", "20-30"]
    assert [float(row["percentile_lower"]) for row in rows] == [10, 20]
    assert [float(row["percentile_upper"]) for row in rows] == [20, 30]


def test_peak_state_equal_rank_bins_split_ties_into_complete_groups():
    scores = np.ones(20, dtype=float)
    state_codes = np.asarray([0] * 10 + [1] * 10, dtype=np.int32)
    selections = peak_states.choose_thresholds(
        scores,
        score_percentile=0,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_bin_size=25,
        pct_bin_seed=7,
    )

    snapshots = peak_states.count_threshold_states(
        scores,
        state_codes,
        selections,
        state_count=2,
    )

    assert [selection.label for selection in selections] == [
        "0-25", "25-50", "50-75", "75-100"
    ]
    assert all(selection.is_rank_bin for selection in selections)
    assert [snapshot.retained_peaks for snapshot in snapshots] == [5, 5, 5, 5]
    assert sum(int(snapshot.counts.sum()) for snapshot in snapshots) == 20
    assert any(tuple(snapshot.counts) not in {(5, 0), (0, 5)} for snapshot in snapshots)


def test_one_percent_equal_rank_bins_keep_all_tied_bins_populated():
    scores = np.ones(100, dtype=float)
    state_codes = np.zeros(100, dtype=np.int32)
    selections = peak_states.choose_thresholds(
        scores,
        score_percentile=0,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_bin_size=1,
        pct_bin_seed=1,
    )
    snapshots = peak_states.count_threshold_states(
        scores,
        state_codes,
        selections,
        state_count=1,
    )

    assert len(selections) == 100
    assert selections[0].label == "0-1"
    assert selections[-1].label == "99-100"
    assert all(snapshot.retained_peaks == 1 for snapshot in snapshots)


def test_peak_state_equal_rank_cli_reports_rank_bins(tmp_path):
    peaks, states = _write_inputs(tmp_path)
    prefix = tmp_path / "rank_bins"

    assert peak_states.main(
        [
            str(peaks),
            "--state-bed", str(states),
            "--pct-bin-size", "50",
            "--pct-bin-seed", "17",
            "--output-prefix", str(prefix),
        ]
    ) == 0
    prefix = Path(f"{prefix}_groupsbinsize-50_tiessplit_overlapfirst")

    rows = list(
        csv.DictReader(
            Path(f"{prefix}_peak_state_threshold_summary.tsv").open(),
            delimiter="\t",
        )
    )
    assert [row["percentile_range"] for row in rows] == ["0-50", "50-100"]
    assert [row["percentile_mode"] for row in rows] == ["rank_bin", "rank_bin"]
    assert [int(row["retained_peak_count"]) for row in rows] == [5, 5]
    assert [int(row["rank_seed"]) for row in rows] == [17, 17]
    assert [row["score_upper_bound"] for row in rows] == ["", ""]
    assert all(row["rank_score_max"] for row in rows)


def test_peak_state_variable_percentage_bins_use_exact_rank_slices():
    scores = np.arange(10, dtype=float)
    state_codes = np.zeros(10, dtype=np.int32)
    selections = peak_states.choose_thresholds(
        scores,
        score_percentile=0,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_values=["0,10,30,100"],
        pct_bins=True,
        bin_tie_mode="split",
        pct_bin_seed=3,
    )
    snapshots = peak_states.count_threshold_states(
        scores, state_codes, selections, state_count=1
    )

    assert [selection.label for selection in selections] == ["0-10", "10-30", "30-100"]
    assert all(selection.is_rank_bin for selection in selections)
    assert [snapshot.retained_peaks for snapshot in snapshots] == [1, 2, 7]
    assert sum(snapshot.retained_peaks for snapshot in snapshots) == 10


def test_peak_state_fixed_width_percentile_mode_keeps_score_ties_together():
    scores = np.asarray([0, 1, 1, 2, 3, 4, 5, 6, 7, 8], dtype=float)
    state_codes = np.zeros(scores.size, dtype=np.int32)
    selections = peak_states.choose_thresholds(
        scores,
        score_percentile=0,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_bin_size=10,
        bin_tie_mode="keep",
    )
    snapshots = peak_states.count_threshold_states(
        scores, state_codes, selections, state_count=1
    )

    assert not any(selection.is_rank_bin for selection in selections)
    assert len(selections) == 10
    assert sum(snapshot.retained_peaks for snapshot in snapshots) == 10


def test_peak_state_plot_uses_equally_spaced_categorical_x_positions(tmp_path, monkeypatch):
    from matplotlib.axes import Axes

    captured_x: list[list[float]] = []
    original_bar = Axes.bar

    def capture_bar(self, x, *args, **kwargs):
        captured_x.append([float(value) for value in x])
        return original_bar(self, x, *args, **kwargs)

    monkeypatch.setattr(Axes, "bar", capture_bar)
    selections = [
        peak_states.ThresholdSelection(percentile, float(percentile))
        for percentile in (0, 90, 99, 99.9)
    ]
    snapshots = [
        peak_states.ThresholdSnapshot(
            counts=np.asarray([1], dtype=np.int64),
            retained_peaks=1,
            assigned_peaks=1,
            unassigned_peaks=0,
        )
        for _selection in selections
    ]
    peak_states.plot_state_percentages(
        tmp_path / "categorical.png",
        selections=selections,
        snapshots=snapshots,
        state_order=["1_State"],
        colors={"1_State": "255,0,0"},
        title=None,
        dpi=100,
    )
    assert captured_x == [[0.0, 1.0, 2.0, 3.0]]


def test_peak_state_plot_can_use_continuous_threshold_positions(tmp_path, monkeypatch):
    from matplotlib.axes import Axes

    captured_x: list[list[float]] = []
    original_bar = Axes.bar

    def capture_bar(self, x, *args, **kwargs):
        captured_x.append([float(value) for value in x])
        return original_bar(self, x, *args, **kwargs)

    monkeypatch.setattr(Axes, "bar", capture_bar)
    selections = [
        peak_states.ThresholdSelection(percentile, float(percentile))
        for percentile in (0, 90, 99, 99.9)
    ]
    snapshots = [
        peak_states.ThresholdSnapshot(
            counts=np.asarray([1], dtype=np.int64),
            retained_peaks=1,
            assigned_peaks=1,
            unassigned_peaks=0,
        )
        for _selection in selections
    ]
    peak_states.plot_state_percentages(
        tmp_path / "continuous.png",
        selections=selections,
        snapshots=snapshots,
        state_order=["1_State"],
        colors={"1_State": "255,0,0"},
        title=None,
        dpi=100,
        x_axis_mode="continuous",
    )
    assert captured_x == [[0.0, 90.0, 99.0, 99.9]]


def test_peak_state_plot_x_axis_defaults_to_categorical():
    args = peak_states.build_parser().parse_args(
        ["peaks.bed", "--state-bed", "states.bed"]
    )
    assert args.plot_x_axis == "categorical"
    assert args.plot_bar_gap == 0.18
    assert args.bin_tie_mode == "split"


def test_peak_state_plot_zero_bar_gap_removes_categorical_spacing(tmp_path, monkeypatch):
    from matplotlib.axes import Axes

    captured_widths: list[float] = []
    original_bar = Axes.bar

    def capture_bar(self, x, *args, **kwargs):
        captured_widths.extend(float(value) for value in np.atleast_1d(kwargs["width"]))
        return original_bar(self, x, *args, **kwargs)

    monkeypatch.setattr(Axes, "bar", capture_bar)
    selections = [
        peak_states.ThresholdSelection(percentile, float(percentile))
        for percentile in (0, 50, 90)
    ]
    snapshots = [
        peak_states.ThresholdSnapshot(np.asarray([1]), 1, 1, 0)
        for _selection in selections
    ]
    peak_states.plot_state_percentages(
        tmp_path / "no_gap.png",
        selections=selections,
        snapshots=snapshots,
        state_order=["1_State"],
        colors={"1_State": "255,0,0"},
        title=None,
        dpi=100,
        bar_gap=0,
    )
    assert captured_widths == [1.0]


def test_continuous_peak_state_bins_use_midpoints_and_proportional_widths(
    tmp_path, monkeypatch
):
    from matplotlib.axes import Axes

    captured: list[tuple[list[float], list[float]]] = []
    original_bar = Axes.bar

    def capture_bar(self, x, *args, **kwargs):
        widths = np.atleast_1d(kwargs["width"])
        captured.append(
            ([float(value) for value in x], [float(value) for value in widths])
        )
        return original_bar(self, x, *args, **kwargs)

    monkeypatch.setattr(Axes, "bar", capture_bar)
    selections = [
        peak_states.ThresholdSelection(10, 1, 0, 10, 2),
        peak_states.ThresholdSelection(50, 2, 10, 50, 6),
    ]
    snapshots = [
        peak_states.ThresholdSnapshot(np.asarray([1]), 1, 1, 0),
        peak_states.ThresholdSnapshot(np.asarray([1]), 1, 1, 0),
    ]
    peak_states.plot_state_percentages(
        tmp_path / "continuous_bins.png",
        selections=selections,
        snapshots=snapshots,
        state_order=["1_State"],
        colors={"1_State": "255,0,0"},
        title=None,
        dpi=100,
        x_axis_mode="continuous",
    )
    assert captured[0][0] == [5.0, 30.0]
    assert np.allclose(captured[0][1], [8.2, 32.8])


def test_peak_states_are_written_in_numerical_label_order(tmp_path):
    states = tmp_path / "unordered_states.bed"
    states.write_text(
        "chr1\t20\t30\t10_State\t0\t.\t20\t30\t0,0,255\n"
        "chr1\t10\t20\t2_State\t0\t.\t10\t20\t0,255,0\n"
        "chr1\t0\t10\t1_State\t0\t.\t0\t10\t255,0,0\n",
        encoding="utf-8",
    )
    peaks = tmp_path / "peaks.bed"
    peaks.write_text(
        "chr1\t1\t2\ta\t1\nchr1\t11\t12\tb\t2\nchr1\t21\t22\tc\t3\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "ordered"
    assert peak_states.main(
        [str(peaks), "--state-bed", str(states), "--output-prefix", str(prefix)]
    ) == 0
    prefix = Path(f"{prefix}_groupsthreshold-0_tiessplit_overlapfirst")
    rows = list(
        csv.DictReader(Path(f"{prefix}_state_coverage.tsv").open(), delimiter="\t")
    )
    assert [row["state"] for row in rows] == ["1_State", "2_State", "10_State"]
