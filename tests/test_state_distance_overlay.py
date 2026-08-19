"""Tests for ChromHMM-coloured relative peak-distance overlays."""

from __future__ import annotations

import csv
from pathlib import Path

from nucleosuite import distances


def test_state_overlay_is_relative_and_resets_at_interval_boundaries(tmp_path):
    states = tmp_path / "states.bed"
    states.write_text(
        "chr1\t0\t500\t1_Active\t0\t.\t0\t500\t255,0,0\n"
        "chr1\t500\t1000\t2_Repressed\t0\t.\t500\t1000\t0,0,255\n"
    )
    peaks = tmp_path / "peaks.bed"
    centres = [100, 200, 350, 450, 550, 750, 900]
    peaks.write_text(
        "".join(
            f"chr1\t{centre-10}\t{centre+11}\tp{index}\t10\t.\t{centre}\t10\n"
            for index, centre in enumerate(centres)
        )
    )
    prefix = tmp_path / "distances"
    assert distances.main(
        [
            str(peaks),
            "--position-column", "7",
            "--score-column", "5",
            "--score-percentile", "0",
            "--min-distance", "1",
            "--max-distance", "500",
            "--max-order", "1",
            "--scope", "genome",
            "--state-bed", str(states),
            "--state-label-column", "4",
            "--state-color-column", "9",
            "--state-overlay-plot",
            "--state-overlay-smooth-window", "21",
            "--state-overlay-smooth-polyorder", "2",
            "--output-prefix", str(prefix),
        ]
    ) == 0
    prefix = Path(
        f"{prefix}_distmin1_distmax500_orders1"
    )

    stem = Path(f"{prefix}_scorepct0")
    distribution = Path(f"{stem}_state_relative_percent.tsv")
    summary = Path(f"{stem}_state_relative_percent_summary.tsv")
    plot = Path(f"{stem}_state_relative_percent.png")
    assert distribution.is_file() and summary.is_file() and plot.is_file()

    rows = list(csv.DictReader(distribution.open(), delimiter="\t"))
    by_state: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_state.setdefault(row["state"], []).append(row)
    assert set(by_state) == {"1_Active", "2_Repressed"}
    for state_rows in by_state.values():
        assert abs(sum(float(row["raw_percent"]) for row in state_rows) - 100.0) < 1e-6
        assert abs(sum(float(row["smoothed_percent"]) for row in state_rows) - 100.0) < 1e-6

    active_counts = {int(row["distance_bp"]): int(row["count"]) for row in by_state["1_Active"]}
    repressed_counts = {int(row["distance_bp"]): int(row["count"]) for row in by_state["2_Repressed"]}
    assert active_counts[100] == 2
    assert active_counts[150] == 1
    assert repressed_counts[150] == 1
    assert repressed_counts[200] == 1
    # The 100 bp pair from 450 to 550 crosses the state boundary and is excluded.
    assert sum(active_counts.values()) + sum(repressed_counts.values()) == 5

    summary_rows = list(csv.DictReader(summary.open(), delimiter="\t"))
    assert {row["rgb"] for row in summary_rows} == {"255,0,0", "0,0,255"}
    assert all(row["smoothing_window"] == "21" for row in summary_rows)
    assert all(row["smoothing_polyorder"] == "2" for row in summary_rows)
    assert all(row["smoothed_mode_bp"] for row in summary_rows)
    assert all(row["smoothed_mean_bp"] for row in summary_rows)
    assert all(row["smoothed_median_bp"] for row in summary_rows)


def test_state_overlay_resolves_chr_prefix_aliases(tmp_path):
    states = tmp_path / "states.bed"
    states.write_text(
        "chr20\t0\t500\t1_Active\t0\t.\t0\t500\t255,0,0\n"
    )
    peaks = tmp_path / "peaks.bed"
    peaks.write_text(
        "20\t90\t111\tp1\t10\t.\t100\t10\n"
        "20\t190\t211\tp2\t10\t.\t200\t10\n"
        "20\t290\t311\tp3\t10\t.\t300\t10\n"
    )
    prefix = tmp_path / "alias_distances"
    assert distances.main(
        [
            str(peaks),
            "--position-column", "7",
            "--score-column", "5",
            "--score-percentile", "0",
            "--min-distance", "1",
            "--max-distance", "500",
            "--max-order", "1",
            "--scope", "genome",
            "--state-bed", str(states),
            "--state-label-column", "4",
            "--state-color-column", "9",
            "--state-overlay-plot",
            "--output-prefix", str(prefix),
        ]
    ) == 0
    prefix = Path(
        f"{prefix}_distmin1_distmax500_orders1"
    )

    distribution = Path(f"{prefix}_scorepct0_state_relative_percent.tsv")
    rows = list(csv.DictReader(distribution.open(), delimiter="\t"))
    active_counts = {
        int(row["distance_bp"]): int(row["count"])
        for row in rows
        if row["state"] == "1_Active"
    }
    assert active_counts[100] == 2
