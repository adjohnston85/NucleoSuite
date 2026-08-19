"""Tests for order-mode NRL regression in ``nucleosuite distances``."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pytest

from nucleosuite import distances
from nucleosuite.distances import DistanceResults


def synthetic_results() -> DistanceResults:
    return DistanceResults(
        chrom_state={},
        chrom_all={
            1: {
                "chr1": Counter({182: 15, 181: 2}),
                "chr2": Counter({190: 10, 191: 1}),
            },
            2: {
                "chr1": Counter({367: 14, 366: 2}),
                "chr2": Counter({380: 9, 381: 1}),
            },
            3: {
                "chr1": Counter({552: 13, 551: 2}),
                "chr2": Counter({570: 8, 571: 1}),
            },
        },
        genome_state={},
        genome_all={
            1: Counter({185: 25, 184: 3}),
            2: Counter({370: 23, 369: 3}),
            3: Counter({555: 21, 554: 3}),
        },
        duplicates={},
        retained_by_chrom={},
        threshold_pass_count=100,
        retained_count=100,
    )


def test_highest_count_distance_prefers_smaller_distance_on_ties():
    peak = distances.highest_count_distance(Counter({186: 7, 185: 7, 184: 2}))
    assert peak.peak_distance == 185
    assert peak.peak_count == 7
    assert peak.total_pairs == 16


def test_collect_nrl_regressions_separates_genome_and_chromosomes():
    regressions = distances.collect_nrl_regressions(
        synthetic_results(),
        max_order=3,
        include_chromosomes=True,
        include_genome=True,
    )
    keyed = {(regression.scope, regression.chromosome): regression for regression in regressions}

    genome = keyed[("combined_chromosomes", ".")]
    assert genome.slope == pytest.approx(185.0)
    assert genome.intercept == pytest.approx(0.0)
    assert genome.r_squared == pytest.approx(1.0)

    chr1 = keyed[("chromosome", "chr1")]
    assert chr1.slope == pytest.approx(185.0)
    assert chr1.intercept == pytest.approx(-3.0)
    assert chr1.r_squared == pytest.approx(1.0)

    chr2 = keyed[("chromosome", "chr2")]
    assert chr2.slope == pytest.approx(190.0)
    assert chr2.intercept == pytest.approx(0.0)
    assert chr2.r_squared == pytest.approx(1.0)


def test_write_nrl_regression_outputs_writes_tsv_and_png_per_scope(tmp_path):
    regressions = distances.collect_nrl_regressions(
        synthetic_results(),
        max_order=3,
        include_chromosomes=True,
        include_genome=True,
    )
    prefix = tmp_path / "sample_scorepct0"
    outputs = distances.write_nrl_regression_outputs(regressions, prefix=prefix)

    expected = {
        Path(f"{prefix}_nrl_regression_summary.tsv"),
        Path(f"{prefix}_nrl_regression_combined_chromosomes.tsv"),
        Path(f"{prefix}_nrl_regression_combined_chromosomes.png"),
        Path(f"{prefix}_nrl_regression_chromosome_chr1.tsv"),
        Path(f"{prefix}_nrl_regression_chromosome_chr1.png"),
        Path(f"{prefix}_nrl_regression_chromosome_chr2.tsv"),
        Path(f"{prefix}_nrl_regression_chromosome_chr2.png"),
    }
    assert set(outputs) == expected
    assert all(path.is_file() and path.stat().st_size > 0 for path in expected)

    summary_rows = list(
        csv.DictReader(Path(f"{prefix}_nrl_regression_summary.tsv").open(), delimiter="\t")
    )
    assert len(summary_rows) == 3
    genome = next(row for row in summary_rows if row["scope"] == "combined_chromosomes")
    assert float(genome["nrl_bp"]) == pytest.approx(185.0)
    assert float(genome["r_squared"]) == pytest.approx(1.0)

    points = list(
        csv.DictReader(Path(f"{prefix}_nrl_regression_combined_chromosomes.tsv").open(), delimiter="\t")
    )
    assert [int(row["order"]) for row in points] == [1, 2, 3]
    assert [int(row["peak_distance_bp"]) for row in points] == [185, 370, 555]
    assert all(float(row["residual_bp"]) == pytest.approx(0.0) for row in points)


def test_distance_regression_plot_defaults_to_open_circles_dotted_fit_and_square(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    regression = next(
        item for item in distances.collect_nrl_regressions(
            synthetic_results(), max_order=3, include_chromosomes=False, include_genome=True
        ) if item.scope == "combined_chromosomes"
    )
    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda figure=None: None)
    distances.plot_nrl_regression(tmp_path / "distance_regression.png", regression)
    figure = plt.gcf()
    axis = figure.axes[0]
    assert tuple(figure.get_size_inches()) == (6.5, 6.5)
    assert axis.lines[0].get_linestyle() == ":"
    assert axis.collections[0].get_facecolors().size == 0
    original_close(figure)


def test_distances_main_automatically_writes_regression_for_multiple_orders(tmp_path):
    peak_bed = tmp_path / "peaks.bed"
    # Repeated 185 bp spacing gives exact modes at 185, 370 and 555 bp.
    peak_bed.write_text(
        "".join(
            f"chr1\t{position}\t{position + 1}\tp{index}\t10\t.\t{position}\t{position + 1}\n"
            for index, position in enumerate(range(100, 1210, 185), start=1)
        )
    )
    prefix = tmp_path / "distances"
    assert distances.main(
        [
            str(peak_bed),
            "--position-column", "7",
            "--score-column", "5",
            "--max-order", "3",
            "--max-distance", "1000",
            "--scope", "all",
            "--count-smooth-window", "0",
            "--percent-smooth-window", "0",
            "--output-prefix", str(prefix),
        ]
    ) == 0

    prefix = Path(
        f"{prefix}_distmin1_distmax1000_orders3"
    )
    stem = Path(f"{prefix}_scorepct0")
    summary = Path(f"{stem}_nrl_regression_summary.tsv")
    assert summary.is_file()
    rows = list(csv.DictReader(summary.open(), delimiter="\t"))
    assert {(row["scope"], row["chromosome"]) for row in rows} == {
        ("combined_chromosomes", "."),
        ("chromosome", "chr1"),
    }
    assert all(float(row["nrl_bp"]) == pytest.approx(185.0) for row in rows)
    assert Path(f"{stem}_nrl_regression_combined_chromosomes.png").is_file()
    assert Path(f"{stem}_nrl_regression_chromosome_chr1.png").is_file()


def test_single_order_does_not_write_nrl_regression_outputs(tmp_path):
    peak_bed = tmp_path / "peaks.bed"
    peak_bed.write_text(
        "chr1\t100\t101\tp1\t10\t.\t100\t101\n"
        "chr1\t285\t286\tp2\t10\t.\t285\t286\n"
    )
    prefix = tmp_path / "single"
    assert distances.main(
        [
            str(peak_bed),
            "--position-column", "7",
            "--max-order", "1",
            "--output-prefix", str(prefix),
        ]
    ) == 0
    prefix = Path(
        f"{prefix}_distmin1_distmax10000_orders1_nrlmodesmoothed_countsg21x2_pctsg21x2_"
        "groupsthreshold_tiessplit"
    )
    assert not Path(f"{prefix}_scorepct0_nrl_regression_summary.tsv").exists()


def test_distances_parser_defaults_nrl_to_smoothed_mode():
    args = distances.build_parser().parse_args(["peaks.bed"])
    assert args.nrl_mode == "smoothed"
    assert args.count_smooth_window == 21
    assert args.count_smooth_polyorder == 2


def test_smoothed_nrl_mode_uses_smoothed_distribution_mode(monkeypatch):
    results = synthetic_results()
    original = distances.summarize_distribution

    def fake_summary(counter, **kwargs):
        stats = original(counter, **kwargs)
        # Shift only the pooled order modes so this test distinguishes the path.
        if 185 in counter:
            return stats.__class__(**{**stats.__dict__, "smoothed_mode": 186})
        if 370 in counter:
            return stats.__class__(**{**stats.__dict__, "smoothed_mode": 372})
        if 555 in counter:
            return stats.__class__(**{**stats.__dict__, "smoothed_mode": 558})
        return stats

    monkeypatch.setattr(distances, "summarize_distribution", fake_summary)
    regressions = distances.collect_nrl_regressions(
        results, max_order=3, include_chromosomes=False, include_genome=True,
        nrl_mode="smoothed", count_smooth_window=21, count_smooth_polyorder=2,
    )
    assert [p.peak_distance for p in regressions[0].peaks] == [186, 372, 558]
    assert regressions[0].slope == pytest.approx(186.0)
