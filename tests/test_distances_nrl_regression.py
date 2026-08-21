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
    }
    assert all(float(row["nrl_bp"]) == pytest.approx(185.0) for row in rows)
    assert Path(f"{stem}_nrl_regression_combined_chromosomes.png").is_file()
    assert not Path(f"{stem}_nrl_regression_chromosome_chr1.png").exists()


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


def test_smoothed_peak_calling_uses_full_profile_before_range_filtering():
    # The profile is still rising at the requested 1500-bp boundary, then peaks
    # outside the regression range. 1500 must therefore not be called as a peak.
    counter = Counter({1300: 50, 1299: 20, 1301: 20, 1490: 30, 1500: 40, 1510: 60, 1520: 20})
    peak = distances.select_order_peak_in_range(
        counter, min_distance=1, max_distance=1500, nrl_mode="smoothed",
        count_smooth_window=21, count_smooth_polyorder=2,
    )
    assert peak is not None
    assert peak.peak_distance == 1300


def test_peak_exactly_at_regression_max_is_allowed_when_full_profile_declines_after_it():
    counter = Counter({1490: 20, 1499: 40, 1500: 100, 1501: 40, 1510: 10})
    peak = distances.select_order_peak_in_range(
        counter, min_distance=1, max_distance=1500, nrl_mode="raw",
        count_smooth_window=21, count_smooth_polyorder=2,
    )
    assert peak is not None
    assert peak.peak_distance == 1500


def test_distances_parser_defaults_to_combined_regression_and_1500_bp_max():
    args = distances.build_parser().parse_args(["peaks.bed"])
    assert args.scope == "combined_chromosomes"
    assert args.regression_scope == "combined"
    assert args.max_distance == 1500
    assert args.label_peaks is False


def test_full_distribution_mode_is_not_replaced_by_an_in_range_secondary_peak():
    counter = Counter({
        1299: 30, 1300: 80, 1301: 30,
        1599: 50, 1600: 160, 1601: 50,
    })
    full = distances.select_order_peak_full_profile(
        counter, nrl_mode="smoothed", count_smooth_window=21, count_smooth_polyorder=2
    )
    assert full is not None
    assert full.peak_distance > 1500
    assert distances.select_order_peak_in_range(
        counter, min_distance=1, max_distance=1500, nrl_mode="smoothed",
        count_smooth_window=21, count_smooth_polyorder=2,
    ) is None
    retained = distances.select_order_peak_in_range(
        counter, min_distance=1, max_distance=2000, nrl_mode="smoothed",
        count_smooth_window=21, count_smooth_polyorder=2,
    )
    assert retained is not None
    assert retained.peak_distance == full.peak_distance


def test_distribution_outputs_are_not_truncated_by_plot_maximum(tmp_path):
    results = DistanceResults(
        chrom_state={}, chrom_all={}, genome_state={},
        genome_all={1: Counter({185: 10, 1600: 20})},
        duplicates={}, retained_by_chrom={}, threshold_pass_count=30, retained_count=30,
    )
    distance_path = tmp_path / "distances.tsv"
    summary_path = tmp_path / "summary.tsv"
    distances.write_distribution_outputs(
        results, distance_path=distance_path, summary_path=summary_path, max_order=1,
        include_chromosomes=False, include_genome=True, include_state_strata=False,
        include_zero_distances=False, count_smooth_window=0, count_smooth_polyorder=2,
        percent_smooth_window=0, percent_smooth_polyorder=2, min_distance=1, max_distance=1500,
    )
    rows = list(csv.DictReader(distance_path.open(), delimiter="\t"))
    assert {int(row["distance_bp"]) for row in rows} == {185, 1600}
    summary = next(csv.DictReader(summary_path.open(), delimiter="\t"))
    assert int(summary["raw_mode_bp"]) == 1600
    assert int(summary["smoothed_mode_bp"]) == 1600


def test_distance_peak_label_option_is_explicit():
    args = distances.build_parser().parse_args([
        "peaks.bed", "--label-peaks", "--peak-label-value", "both", "--peak-label-offset", "7"
    ])
    assert args.label_peaks is True
    assert args.peak_label_value == "both"
    assert args.peak_label_offset == pytest.approx(7.0)


def test_distances_zero_filling_is_default_and_can_be_disabled():
    default_args = distances.build_parser().parse_args(["peaks.bed"])
    assert default_args.include_zero_distances is True

    sparse_args = distances.build_parser().parse_args([
        "peaks.bed", "--no-include-zero-distances"
    ])
    assert sparse_args.include_zero_distances is False


def test_zero_filled_distribution_is_bounded_by_requested_reporting_range(tmp_path):
    results = DistanceResults(
        chrom_state={}, chrom_all={}, genome_state={},
        genome_all={1: Counter({185: 10, 8_000_000: 1})},
        duplicates={}, retained_by_chrom={}, threshold_pass_count=11, retained_count=11,
    )
    distance_path = tmp_path / "distances.tsv"
    summary_path = tmp_path / "summary.tsv"
    distributions, rows_written = distances.write_distribution_outputs(
        results,
        distance_path=distance_path,
        summary_path=summary_path,
        max_order=1,
        include_chromosomes=False,
        include_genome=True,
        include_state_strata=False,
        include_zero_distances=True,
        count_smooth_window=0,
        count_smooth_polyorder=2,
        percent_smooth_window=0,
        percent_smooth_polyorder=2,
        min_distance=1,
        max_distance=1500,
    )
    assert distributions == 1
    assert rows_written == 1500

    rows = list(csv.DictReader(distance_path.open(), delimiter="\t"))
    assert len(rows) == 1500
    assert int(rows[0]["distance_bp"]) == 1
    assert int(rows[-1]["distance_bp"]) == 1500
    by_distance = {int(row["distance_bp"]): int(row["count"]) for row in rows}
    assert by_distance[184] == 0
    assert by_distance[185] == 10
    assert by_distance[186] == 0
    assert 8_000_000 not in by_distance

    summary = next(csv.DictReader(summary_path.open(), delimiter="\t"))
    assert int(summary["observed_max_bp"]) == 8_000_000
    assert int(summary["raw_mode_bp"]) == 185


def test_distances_length_filter_is_applied_before_score_percentile(tmp_path):
    peak_bed = tmp_path / "length_peaks.bed"
    peak_bed.write_text(
        "chr1\t0\t50\tshort_high\t1000\n"
        "chr1\t100\t200\ta\t1\n"
        "chr1\t300\t400\tb\t2\n"
        "chr1\t500\t600\tc\t3\n"
    )
    prefix = tmp_path / "length_filter"
    assert distances.main(
        [
            str(peak_bed),
            "--min-length", "100",
            "--max-length", "100",
            "--score-percentile", "50",
            "--max-order", "1",
            "--max-distance", "500",
            "--count-smooth-window", "0",
            "--percent-smooth-window", "0",
            "--output-prefix", str(prefix),
        ]
    ) == 0
    metadata = Path(
        f"{prefix}_distmin1_distmax500_orders1_lenmin100_lenmax100_scorepct50_metadata.tsv"
    )
    assert metadata.is_file()
    rows = dict(
        line.rstrip("\n").split("\t", 1)
        for line in metadata.read_text().splitlines()[1:]
    )
    assert float(rows["score_threshold"]) == pytest.approx(2.0)
    assert rows["min_peak_length"] == "100"
    assert rows["max_peak_length"] == "100"
    assert rows["retained_count"] == "2"
