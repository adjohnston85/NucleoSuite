from __future__ import annotations

import csv
from pathlib import Path

from nucleosuite.profile_plots import plot_profile_overlay


def _write_profile(path: Path, values: list[float]) -> Path:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("relative_position\tscore\n")
        for position, value in zip((-10, 0, 10), values):
            handle.write(f"{position}\t{value}\n")
    return path


def test_profile_overlay_writes_wide_tsv_and_png(tmp_path: Path) -> None:
    active = _write_profile(tmp_path / "active.tsv", [1.0, 2.0, 3.0])
    repressed = _write_profile(tmp_path / "repressed.tsv", [3.0, 2.0, 1.0])
    output_tsv = tmp_path / "combined.tsv"
    output_png = tmp_path / "combined.png"

    plot_profile_overlay(
        [("active_genes", active), ("repressed_genes", repressed)],
        output_tsv,
        output_png,
        xlabel="Position relative to TSS (bp)",
        ylabel="Mean PNS",
        title="TSS profiles",
    )

    with output_tsv.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert list(rows[0]) == ["relative_position", "active_genes", "repressed_genes"]
    assert [row["relative_position"] for row in rows] == ["-10", "0", "10"]
    assert output_png.is_file() and output_png.stat().st_size > 0
