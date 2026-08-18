"""Regression tests for percentile-sweep distance figures and progress."""

from __future__ import annotations

import csv
import io
from collections import Counter
from pathlib import Path

import numpy as np

from nucleosuite import distances
from nucleosuite.progress import ProgressReporter


def _write_peak_bed(path: Path) -> None:
    rows = []
    index = 0
    for chromosome, spacing in (("chr1", 185), ("chr2", 190)):
        for rank in range(8):
            index += 1
            position = 100 + rank * spacing
            rows.append(
                f"{chromosome}\t{position}\t{position + 1}\tp{index}\t"
                f"{rank + 1}\t.\t{position}\t{position + 1}\n"
            )
    path.write_text("".join(rows), encoding="utf-8")


def test_pct_range_writes_separate_count_percentage_and_peak_count_figures(tmp_path):
    peak_bed = tmp_path / "peaks.bed"
    _write_peak_bed(peak_bed)
    prefix = tmp_path / "distance_sweep"

    assert distances.main(
        [
            str(peak_bed),
            "--position-column", "7",
            "--score-column", "5",
            "--pct-range",
            "--pct-lower", "0",
            "--pct-upper", "75",
            "--pct-step", "25",
            "--min-distance", "1",
            "--max-distance", "500",
            "--max-order", "1",
            "--scope", "all",
            "--count-smooth-window", "0",
            "--percent-smooth-window", "0",
            "--output-prefix", str(prefix),
        ]
    ) == 0

    stems = [
        f"{prefix}_percentile_sweep_order1_combined_chromosomes",
        f"{prefix}_percentile_sweep_order1_chromosome_chr1",
        f"{prefix}_percentile_sweep_order1_chromosome_chr2",
    ]
    for stem in stems:
        for suffix in ("count.png", "percentage.png", "peak_counts.png"):
            path = Path(f"{stem}_{suffix}")
            assert path.is_file() and path.stat().st_size > 0

    count_table = Path(f"{prefix}_percentile_sweep_peak_counts.tsv")
    rows = list(csv.DictReader(count_table.open(), delimiter="\t"))
    assert len(rows) == 12
    combined = [row for row in rows if row["scope"] == "combined_chromosomes"]
    assert [float(row["percentile_threshold"]) for row in combined] == [0, 25, 50, 75]
    assert [int(row["retained_peak_count"]) for row in combined] == [16, 12, 8, 4]
    assert not list(tmp_path.glob("*.percentile_plots.sqlite"))


def test_explicit_percentile_values_activate_sweep_and_override_range(tmp_path):
    peak_bed = tmp_path / "peaks.bed"
    _write_peak_bed(peak_bed)
    prefix = tmp_path / "explicit_sweep"

    assert distances.main(
        [
            str(peak_bed),
            "--position-column", "7",
            "--pct-values", "10,20,50,90,99",
            "--pct-lower", "0",
            "--pct-upper", "5",
            "--max-distance", "500",
            "--max-order", "1",
            "--scope", "combined_chromosomes",
            "--output-prefix", str(prefix),
        ]
    ) == 0

    table = Path(f"{prefix}_percentile_sweep_peak_counts.tsv")
    rows = list(csv.DictReader(table.open(), delimiter="\t"))
    assert [float(row["percentile_threshold"]) for row in rows] == [
        10, 20, 50, 90, 99
    ]


def test_explicit_percentile_bins_write_non_cumulative_ranges(tmp_path):
    peak_bed = tmp_path / "bin_peaks.bed"
    peak_bed.write_text(
        "".join(
            f"chr1\t{100 + index * 185}\t{101 + index * 185}\tp{index}\t{index}\t.\t"
            f"{100 + index * 185}\t{101 + index * 185}\n"
            for index in range(1, 41)
        ),
        encoding="utf-8",
    )
    prefix = tmp_path / "bin_sweep"

    assert distances.main(
        [
            str(peak_bed),
            "--position-column", "7",
            "--pct-values", "10,20,30",
            "--pct-bins",
            "--min-distance", "1",
            "--max-distance", "500",
            "--max-order", "1",
            "--scope", "combined_chromosomes",
            "--output-prefix", str(prefix),
        ]
    ) == 0

    table = Path(f"{prefix}_percentile_sweep_peak_counts.tsv")
    rows = list(csv.DictReader(table.open(), delimiter="\t"))
    assert [row["percentile_mode"] for row in rows] == ["rank_bin", "rank_bin"]
    assert [row["bin_tie_mode"] for row in rows] == ["split", "split"]
    assert [row["percentile_range"] for row in rows] == ["10-20", "20-30"]
    assert [int(row["retained_peak_count"]) for row in rows] == [4, 4]
    for label in ("10-20", "20-30"):
        assert Path(f"{prefix}_scorepct{label}_metadata.tsv").is_file()


def test_equal_rank_distance_bins_are_disjoint_and_cover_every_peak():
    records = {
        "chr1": [
            distances.PeakRecord(position, 1.0, "All", position, position + 1, str(position), ".")
            for position in range(100, 1100, 100)
        ],
        "chr2": [
            distances.PeakRecord(position, 1.0, "All", position, position + 1, str(position), ".")
            for position in range(100, 1100, 100)
        ],
    }
    selections = distances.choose_thresholds(
        np.ones(20, dtype=float),
        score_percentile=0,
        target_peaks=None,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_bin_size=25,
        pct_bin_seed=5,
    )
    distances.attach_rank_bin_records(records, selections)

    membership = selections[0].rank_membership
    assert membership is not None
    selected_sets = [
        {
            id(record)
            for peaks in membership.groups[index].values()
            for record in peaks
        }
        for index in range(len(selections))
    ]
    assert all(len(selected) == 5 for selected in selected_sets)
    assert all(set(group) == {"chr1", "chr2"} for group in membership.groups)
    assert set().union(*selected_sets) == {
        id(record) for peaks in records.values() for record in peaks
    }
    assert sum(len(selected_sets[i] & selected_sets[j]) for i in range(4) for j in range(i)) == 0


def test_distance_percentage_bins_support_varying_widths_and_split_ties_exactly():
    scores = np.ones(20, dtype=float)
    selections = distances.choose_thresholds(
        scores,
        score_percentile=0,
        target_peaks=None,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_values=["0,10,30,100"],
        pct_bins=True,
        pct_bin_seed=11,
    )
    assert [selection.label for selection in selections] == ["0-10", "10-30", "30-100"]
    assert [selection.rank_stop - selection.rank_start for selection in selections] == [2, 4, 14]
    assert all(selection.is_rank_bin for selection in selections)


def test_distance_percentile_bin_mode_keeps_score_ties_together():
    scores = np.ones(10, dtype=float)
    selections = distances.choose_thresholds(
        scores,
        score_percentile=0,
        target_peaks=None,
        pct_range=False,
        pct_lower=0,
        pct_upper=99,
        pct_step=1,
        pct_values=["0,50,100"],
        pct_bins=True,
        bin_tie_mode="keep",
    )
    assert not any(selection.is_rank_bin for selection in selections)
    assert selections[0].threshold == selections[0].score_upper_bound == 1.0
    assert selections[1].threshold == 1.0
    assert selections[1].score_upper_bound is None


def test_distance_bin_tie_mode_defaults_to_split():
    args = distances.build_parser().parse_args(["peaks.bed"])
    assert args.bin_tie_mode == "split"


def test_equal_rank_distance_cli_uses_requested_group_size(tmp_path):
    peak_bed = tmp_path / "rank_peaks.bed"
    peak_bed.write_text(
        "".join(
            f"chr1\t{index * 100}\t{index * 100 + 1}\tp{index}\t1\t.\t"
            f"{index * 100}\t{index * 100 + 1}\n"
            for index in range(40)
        ),
        encoding="utf-8",
    )
    prefix = tmp_path / "rank_sweep"

    assert distances.main(
        [
            str(peak_bed),
            "--position-column", "7",
            "--pct-bin-size", "5",
            "--pct-bin-seed", "23",
            "--min-distance", "1",
            "--max-distance", "5000",
            "--max-order", "1",
            "--scope", "combined_chromosomes",
            "--output-prefix", str(prefix),
        ]
    ) == 0

    rows = list(
        csv.DictReader(
            Path(f"{prefix}_percentile_sweep_peak_counts.tsv").open(),
            delimiter="\t",
        )
    )
    assert len(rows) == 20
    assert all(row["percentile_mode"] == "rank_bin" for row in rows)
    assert all(int(row["retained_peak_count"]) == 2 for row in rows)
    assert rows[0]["percentile_range"] == "0-5"
    assert rows[-1]["percentile_range"] == "95-100"
    assert all(int(row["rank_seed"]) == 23 for row in rows)
    assert all(row["score_upper_bound"] == "" for row in rows)
    assert all(float(row["rank_score_max"]) == 1 for row in rows)
    metadata = dict(
        csv.reader(
            Path(f"{prefix}_scorepct0-5_metadata.tsv").open(),
            delimiter="\t",
        )
    )
    assert metadata["percentile_mode"] == "rank_bin"
    assert metadata["pct_bin_size"] == "5"
    assert metadata["rank_seed"] == "23"
    assert metadata["rank_tie_order"] == (
        "global_random_shuffle_then_stable_ascending_score_sort"
    )


def test_equal_rank_distances_force_one_global_serial_analysis(monkeypatch):
    args = distances.build_parser().parse_args(
        [
            "peaks.bed",
            "--pct-bin-size", "5",
            "--memory-intensive-analysis-cores", "4",
        ]
    )
    observed = []

    def fake_serial(namespace):
        observed.append(namespace)
        return 17

    monkeypatch.setattr(distances, "_run_serial", fake_serial)
    assert distances.run(args) == 17
    assert observed == [args]


def test_percentile_bin_score_bounds_are_lower_inclusive_upper_exclusive():
    records = {
        "chr1": [
            distances.PeakRecord(position, score, "All", position, position + 1, str(position), ".")
            for position, score in ((0, 1.0), (10, 2.0), (20, 2.0), (30, 3.0))
        ]
    }
    results = distances.compute_distance_counts(
        records,
        threshold=1.0,
        score_upper_bound=2.0,
        min_distance=1,
        max_distance=100,
        max_order=1,
        duplicate_policy="highest-score",
    )
    assert results.threshold_pass_count == 1
    assert [record.score for record in results.retained_by_chrom["chr1"]] == [1.0]


def test_distance_percentage_smoothing_defaults_match_count_defaults():
    args = distances.build_parser().parse_args(["peaks.bed"])
    assert (args.count_smooth_window, args.count_smooth_polyorder) == (21, 2)
    assert (args.percent_smooth_window, args.percent_smooth_polyorder) == (21, 2)


def test_percentile_plot_store_persists_raw_unsmoothed_curves(tmp_path):
    prefix = tmp_path / "raw_curves"
    records = [
        distances.PeakRecord(position, 1.0, "All", position, position + 1, str(position), ".")
        for position in (0, 1, 4, 7, 10)
    ]
    results = distances.DistanceResults(
        chrom_state={},
        chrom_all={1: {"chr1": Counter({1: 1, 3: 3})}},
        genome_state={},
        genome_all={1: Counter({1: 1, 3: 3})},
        duplicates={},
        retained_by_chrom={"chr1": records},
        threshold_pass_count=5,
        retained_count=5,
    )
    store = distances.PercentilePlotStore(prefix)
    try:
        store.add_threshold(
            selection=distances.ThresholdSelection(0, 0, 1.0),
            results=results,
            chromosomes=["chr1"],
            max_order=1,
            include_chromosomes=False,
            include_genome=True,
        )
        rows = list(
            store.connection.execute(
                "SELECT distance_bp, raw_count, raw_percent FROM curves "
                "ORDER BY distance_bp"
            )
        )
        assert rows == [(1, 1.0, 25.0), (2, 0.0, 0.0), (3, 3.0, 75.0)]
    finally:
        store.close()


def test_distance_peak_loader_reports_each_contig_once(tmp_path):
    peak_bed = tmp_path / "interleaved.bed"
    peak_bed.write_text(
        "chr1\t10\t11\ta\t1\n"
        "chr2\t20\t21\tb\t2\n"
        "chr1\t30\t31\tc\t3\n",
        encoding="utf-8",
    )
    stream = io.StringIO()
    reporter = ProgressReporter("distances", stream=stream)
    distances.load_peaks(
        peak_bed,
        state_indexes=None,
        progress=reporter,
    )
    output = stream.getvalue()
    assert output.count("Reading peaks contig: chr1") == 1
    assert output.count("Reading peaks contig: chr2") == 1
    assert "records read" not in output


def test_packaged_suites_use_many_to_one_position_matching():
    package_root = Path(__file__).resolve().parents[1] / "src" / "nucleosuite" / "resources"
    for script_name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        script = (package_root / script_name).read_text(encoding="utf-8")
        assert "--matching many-to-one" in script
