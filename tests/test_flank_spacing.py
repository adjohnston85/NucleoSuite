from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from nucleosuite.flank_spacing import (
    compute_flanking_spacings,
    distribution_curve,
    rank_categories,
)


def test_flanking_spacing_uses_strict_upstream_and_downstream(tmp_path: Path) -> None:
    regions = tmp_path / "regions.bed"
    regions.write_text(
        "chr1\t199\t201\tA\n"
        "chr1\t259\t261\tB\n",
        encoding="utf-8",
    )
    centres = {"chr1": [10, 100, 200, 300, 400]}
    detail, by_category, totals = compute_flanking_spacings(
        regions, centres, category_col=4
    )
    # A is centred exactly on a nucleosome at 200; strict flanks are 100 and 300.
    assert by_category["A"] == [200]
    assert by_category["B"] == [100]
    assert totals == {"A": 1, "B": 1}
    assert detail[0]["upstream_nucleosome_center"] == 100
    assert detail[0]["downstream_nucleosome_center"] == 300


def test_category_ranking_is_lowest_x1_over_x2_first() -> None:
    grid = np.arange(0, 501, dtype=float)
    by_category = {
        "wide": [190] * 2 + [260] * 8,
        "canonical": [190] * 8 + [260] * 2,
    }
    rankings, curves = rank_categories(
        by_category, mode="count", ratio_x1=190, ratio_x2=260, x_grid=grid
    )
    assert [row["category"] for row in rankings] == ["wide", "canonical"]
    assert rankings[0]["rank"] == 1
    assert rankings[0]["ratio_190_to_260"] == 0.25
    assert curves["wide"][260] == 8


def test_flank_spacing_cli_writes_ranked_outputs_and_plot_metadata(tmp_path: Path) -> None:
    from nucleosuite.flank_spacing import main
    nuc = tmp_path / "nucleosomes.bed"
    nuc.write_text(
        "chr1\t99\t101\tn1\n"
        "chr1\t189\t191\tn2\n"
        "chr1\t259\t261\tn3\n"
        "chr1\t449\t451\tn4\n",
        encoding="utf-8",
    )
    regions = tmp_path / "regions.bed"
    regions.write_text(
        "chr1\t199\t201\twide\n"
        "chr1\t299\t301\twide\n"
        "chr1\t149\t151\tcanonical\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    code = main([
        "--nucleosome-bed", str(nuc),
        "--region-bed", str(regions),
        "--distribution", "count",
        "--ratio-x1", "90",
        "--ratio-x2", "190",
        "--x-max", "500",
        "--output-dir", str(out),
        "--output-prefix", "sample",
    ])
    assert code == 0
    ranking = next(out.glob("*_ranking.tsv"))
    plot = next(out.glob("*.png"))
    metadata = plot.with_name(plot.stem + "_metadata.tsv")
    assert metadata.is_file()
    rows = list(csv.DictReader(ranking.open(), delimiter="\t"))
    assert rows
    # Automatic suffixes contain exactly three central parameter tokens.
    assert "_distcount_ratio90to190_xmax500" in plot.stem
