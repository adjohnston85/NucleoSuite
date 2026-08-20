"""Tests for multi-callset nearest-position comparison."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from nucleosuite.compare_positions import (
    CompactPeakSet,
    ComparisonArrays,
    PeakChrom,
    _histogram_rows,
    _pairwise_statistics,
    _resolve_inputs,
    _percentile_group_bounds,
    build_parser,
    match_compact_positions,
    read_compact_positions,
    run_comparison,
)


def write_bed(path: Path, rows: list[tuple[object, ...]]) -> None:
    path.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))


def test_default_midpoint_and_score_column_uses_compact_chromosome_arrays(tmp_path: Path):
    bed = tmp_path / "positions.bed"
    write_bed(bed, [("chr1", 10, 20, "p1", 3.5), ("chr1", 30, 41, "p2", 7)])
    records = read_compact_positions(bed, "main", summit_column=None, score_column=5)
    chrom = records.by_chrom["1"]
    assert chrom.summits.tolist() == [15, 35]
    assert chrom.scores.tolist() == [3.5, 7.0]
    assert chrom.indices.tolist() == [1, 2]
    assert records.count == 2
    # Peak storage deliberately contains no per-record BED start/end/name objects.
    assert set(vars(chrom)) == {"name", "summits", "scores", "indices"}
    assert chrom.summits.dtype == np.int64
    assert chrom.scores.dtype == np.float64


def test_smaller_comparison_is_query_but_pairs_keep_main_fields(tmp_path: Path):
    main = tmp_path / "main.bed"
    compare = tmp_path / "compare.bed"
    write_bed(main, [("chr1", 95, 106, "m1", 10, ".", 100), ("chr1", 195, 206, "m2", 20, ".", 200), ("chr1", 295, 306, "m3", 30, ".", 300)])
    write_bed(compare, [("chr1", 201, 212, "c1", 22, ".", 207), ("chr1", 306, 317, "c2", 31, ".", 312)])
    records_main = read_compact_positions(main, "main", summit_column=7, score_column=5)
    records_compare = read_compact_positions(compare, "comparison", summit_column=7, score_column=5)
    result = match_compact_positions(records_main, records_compare, None)
    assert result.query_source == "comparison"
    assert result.main_summit.tolist() == [200, 300]
    assert (result.main_summit + result.signed_distance).tolist() == [207, 312]
    assert result.main_scores.tolist() == [20.0, 30.0]


def _compact_from_summits(summits: list[int], source: str) -> CompactPeakSet:
    values = np.asarray(summits, dtype=np.int64)
    indices = np.arange(1, len(values) + 1, dtype=np.int64)
    scores = indices.astype(np.float64)
    order = np.lexsort((indices, values))
    chrom = PeakChrom("chr1", values[order], scores[order], indices[order])
    return CompactPeakSet(Path(f"{source}.bed"), source, {"1": chrom}, len(values), 0)


def test_one_to_one_matching_is_unique_and_distance_prioritized():
    main = _compact_from_summits([0, 4, 4, 8, 12, 12, 16], "main")
    comparison = _compact_from_summits([2, 2, 6, 10, 14, 14, 18, 100], "comparison")
    result = match_compact_positions(main, comparison, None)
    assert result.query_source == "main"
    assert result.pair_count == main.count
    # No comparison source line can be reused by one-to-one matching.
    assert len(set(result.compare_index.tolist())) == result.pair_count
    # Exact-coordinate duplicates are paired before more distant neighbours.
    exact_main = _compact_from_summits([10, 10, 20], "main")
    exact_compare = _compact_from_summits([10, 10, 11, 20], "comparison")
    exact = match_compact_positions(exact_main, exact_compare, None)
    assert exact.signed_distance.tolist() == [0, 0, 0]


def test_new_parser_defaults_to_main_plus_repeated_compare_and_quartiles():
    args = build_parser().parse_args(["--main-bed", "main.bed", "--compare-bed", "B=b.bed", "--compare-bed", "C=c.bed"])
    assert args.main_bed == "main.bed"
    assert args.compare_beds == ["B=b.bed", "C=c.bed"]
    assert args.percentile_interval == 25
    assert args.stats_test == "nonparametric"
    assert args.p_adjust == "holm"
    assert args.percentile_boxplot_y_max == 200.0
    assert args.score_agreement_distance_max == 100.0
    assert args.show_boxplot_outliers is False
    assert args.stats is False
    assert args.write_detail_tables is False
    assert args.histogram_x_min == -250.0
    assert args.histogram_x_max == 250.0
    assert not any("--main-label" in action.option_strings for action in build_parser()._actions)
    assert _percentile_group_bounds(25) == [(0, 25, "0-25"), (25, 50, "25-50"), (50, 75, "50-75"), (75, 100, "75-100")]


def test_multi_comparison_outputs_use_main_score_percentiles_and_combined_tables(tmp_path: Path):
    main = tmp_path / "main.bed"; b = tmp_path / "b.bed"; c = tmp_path / "c.bed"
    main_rows = []; b_rows = []; c_rows = []
    for index in range(1, 13):
        summit = index * 100
        main_rows.append(("chr1", summit - 5, summit + 5, f"m{index}", index))
        b_rows.append(("chr1", summit - 3, summit + 7, f"b{index}", index * 2))  # +2 midpoint shift
        c_rows.append(("chr1", summit + 3, summit + 13, f"c{index}", index * 3))  # +8 midpoint shift
    # make B larger so main is the query; make C smaller so C is the query
    b_rows.extend([("chr1", 2000, 2010, "b_extra", 1), ("chr1", 2200, 2210, "b_extra2", 1)])
    c_rows = c_rows[:-2]
    write_bed(main, main_rows); write_bed(b, b_rows); write_bed(c, c_rows)
    prefix = tmp_path / "multi"
    args = build_parser().parse_args([
        "--main-bed", str(main), "--compare-bed", f"B={b}", "--compare-bed", f"C={c}",
        "--output-prefix", str(prefix), "--plot-max-points", "100", "--dpi", "40", "--quiet",
    ])
    outputs = run_comparison(args)
    assert outputs["summary"].exists()
    assert outputs["percentile_boxplot_source"].exists()
    assert outputs["percentile_boxplot"].exists()
    assert "percentile_distances" not in outputs
    assert "pairs_B" not in outputs and "pairs_C" not in outputs
    assert outputs["score_agreement_source_B"].exists()
    assert outputs["score_agreement_source_C"].exists()
    with outputs["summary"].open() as handle:
        summary = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["comparison"] for row in summary] == ["B", "C"]
    assert summary[0]["query_source"] == "main"
    assert summary[1]["query_source"] == "comparison"
    with outputs["percentile_boxplot_source"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    box_rows = [row for row in rows if row["row_type"] == "box"]
    assert {row["comparison"] for row in box_rows} == {"B", "C"}
    assert {row["percentile_group"] for row in box_rows} == {"0-25", "25-50", "50-75", "75-100"}
    assert all(row["show_boxplot_outliers"] == "0" for row in box_rows)
    assert outputs["percentile_distance_trend"].exists()
    assert outputs["percentile_distance_trend_plot"].exists()
    for plot_key in ("distance_histogram_plot", "correlation_by_distance_plot", "percentile_boxplot", "percentile_distance_trend_plot", "score_agreement_plot_B", "score_agreement_plot_C"):
        plot = outputs[plot_key]
        assert plot.exists() and plot.stat().st_size > 0
        assert plot.with_name(plot.stem + "_metadata.tsv").exists()


def test_compare_position_large_detail_tables_are_opt_in(tmp_path: Path):
    main = tmp_path / "main.bed"
    compare = tmp_path / "compare.bed"
    write_bed(main, [("chr1", i * 100, i * 100 + 10, f"m{i}", i) for i in range(1, 9)])
    write_bed(compare, [("chr1", i * 100 + 2, i * 100 + 12, f"c{i}", i) for i in range(1, 9)])
    args = build_parser().parse_args([
        "--main-bed", str(main), "--compare-bed", f"C={compare}",
        "--output-prefix", str(tmp_path / "detail"), "--write-detail-tables",
        "--dpi", "40", "--quiet",
    ])
    outputs = run_comparison(args)
    assert outputs["percentile_distances"].exists()
    assert outputs["pairs_C"].exists()


def test_within_percentile_nonparametric_stats_test_all_compare_pairs_and_holm(tmp_path: Path):
    main = tmp_path / "main.bed"; b = tmp_path / "b.bed"; c = tmp_path / "c.bed"; d = tmp_path / "d.bed"
    main_rows=[]; b_rows=[]; c_rows=[]; d_rows=[]
    for index in range(1, 17):
        summit = index * 100
        main_rows.append(("chr1", summit-5, summit+5, f"m{index}", index))
        b_rows.append(("chr1", summit-4, summit+6, f"b{index}", index))   # +1
        c_rows.append(("chr1", summit+5, summit+15, f"c{index}", index))  # +10
        d_rows.append(("chr1", summit+15, summit+25, f"d{index}", index)) # +20
    write_bed(main, main_rows); write_bed(b, b_rows); write_bed(c, c_rows); write_bed(d, d_rows)
    args = build_parser().parse_args([
        "--main-bed", str(main), "--compare-bed", f"B={b}", "--compare-bed", f"C={c}", "--compare-bed", f"D={d}",
        "--output-prefix", str(tmp_path/'stats'), "--stats", "--p-display", "stars", "--dpi", "40", "--quiet",
    ])
    outputs = run_comparison(args)
    with outputs["percentile_statistics"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 12  # B-vs-C, B-vs-D and C-vs-D independently in four quartiles
    assert {row["percentile_group"] for row in rows} == {"0-25", "25-50", "50-75", "75-100"}
    for group in {"0-25", "25-50", "50-75", "75-100"}:
        group_rows = [row for row in rows if row["percentile_group"] == group]
        assert {(row["comparison_1"], row["comparison_2"]) for row in group_rows} == {("B", "C"), ("B", "D"), ("C", "D")}
        assert all(row["p_adjustment"] == "holm" for row in group_rows)
    assert all(row["test"] == "Wilcoxon signed-rank" for row in rows)
    assert all(row["paired"] == "True" for row in rows)
    assert all(int(row["n_paired"]) == 4 for row in rows)


def test_main_bed_accepts_label_equals_path_and_normalizes_for_processing(tmp_path: Path):
    main = tmp_path / "main.bed"
    compare = tmp_path / "compare.bed"
    write_bed(main, [("chr1", 0, 10, "m", 1)])
    write_bed(compare, [("chr1", 1, 11, "c", 1)])
    args = build_parser().parse_args([
        "--main-bed", f"PNS={main}", "--compare-bed", f"DANPOS={compare}"
    ])
    main_path, specs = _resolve_inputs(args)
    assert main_path == main
    assert args.main_bed == str(main)
    assert args._main_label == "PNS"
    assert specs[0].label == "DANPOS"


def test_distance_histogram_uses_signed_distance():
    result = ComparisonArrays(
        label="B", path=Path("b.bed"), main_count=3, compare_count=3,
        query_source="main", query_count=3, target_count=3,
        unmatched_no_target_chrom=0, unmatched_distance=0, unmatched_unique=0,
        chrom_names=("chr1",), chrom_code=np.zeros(3, dtype=np.uint16),
        main_line=np.arange(1, 4, dtype=np.int64),
        compare_line=np.arange(1, 4, dtype=np.int64),
        main_summit=np.asarray([100, 200, 300], dtype=np.int64),
        main_scores=np.asarray([1.0, 2.0, 3.0]),
        compare_scores=np.asarray([1.0, 2.0, 3.0]),
        signed_distance=np.asarray([-10, 0, 10], dtype=np.int64),
        absolute_distance=np.asarray([10.0, 0.0, 10.0]),
        percentile=np.asarray([25, 50, 100], dtype=np.uint8),
        group_index=np.asarray([0, 1, 3], dtype=np.uint8),
        group_names=("0-25", "25-50", "50-75", "75-100"),
    )
    rows = _histogram_rows(result, "PNS", -20, 20, 10)
    counts = {(float(row["bin_start_inclusive"]), float(row["bin_end_exclusive"])): int(row["pair_count"]) for row in rows}
    assert counts[(-10.0, 0.0)] == 1
    assert counts[(0.0, 10.0)] == 1
    assert counts[(10.0, 20.0)] == 1


def test_legacy_two_file_alias_is_accepted(tmp_path: Path):
    main = tmp_path / "main.bed"; compare = tmp_path / "compare.bed"
    write_bed(main, [("chr1", 0, 10, "m", 1)])
    write_bed(compare, [("chr1", 1, 11, "c", 1)])
    args = build_parser().parse_args(["--bed-a", str(main), "--bed-b", str(compare), "--output-prefix", str(tmp_path/'legacy'), "--dpi", "40", "--quiet"])
    outputs = run_comparison(args)
    assert outputs["summary"].exists()


def test_score_bigwig_parser_is_repeatable_and_labelled():
    args = build_parser().parse_args([
        "--main-bed", "main.bed", "--compare-bed", "B=b.bed",
        "--score-bigwig", "Coverage=coverage.bw", "--score-bigwig", "Signal=signal.bigWig",
    ])
    assert args.score_bigwigs == ["Coverage=coverage.bw", "Signal=signal.bigWig"]


def test_score_bigwig_is_score_only_and_samples_main_summits(tmp_path: Path, monkeypatch):
    import sys
    import types

    main = tmp_path / "main.bed"
    compare = tmp_path / "compare.bed"
    bigwig = tmp_path / "coverage.bw"
    write_bed(main, [
        ("chr1", 95, 105, "m1", 1.0),
        ("chr1", 195, 205, "m2", 2.0),
        ("chr1", 295, 305, "m3", 3.0),
        ("chr1", 395, 405, "m4", 4.0),
    ])
    write_bed(compare, [
        ("chr1", 97, 107, "c1", 1.5),
        ("chr1", 197, 207, "c2", 2.5),
        ("chr1", 297, 307, "c3", 3.5),
        ("chr1", 397, 407, "c4", 4.5),
    ])
    bigwig.write_bytes(b"fake")

    class FakeBigWig:
        def isBigWig(self):
            return True
        def chroms(self):
            return {"chr1": 1000}
        def values(self, chrom, start, end, numpy=False):
            values = np.arange(start, end, dtype=float) / 100.0
            # Simulate one missing BigWig base; suite convention converts it to zero.
            if start <= 300 < end:
                values[300 - start] = np.nan
            return values if numpy else values.tolist()
        def close(self):
            pass

    monkeypatch.setitem(sys.modules, "pyBigWig", types.SimpleNamespace(open=lambda path: FakeBigWig()))
    args = build_parser().parse_args([
        "--main-bed", f"PNS={main}", "--compare-bed", f"DANPOS={compare}",
        "--score-bigwig", f"Coverage={bigwig}",
        "--output-prefix", str(tmp_path / "with_signal"), "--plot-max-points", "100",
        "--dpi", "40", "--quiet",
    ])
    outputs = run_comparison(args)

    assert outputs["score_agreement_plot_Coverage"].exists()
    assert outputs["score_agreement_source_Coverage"].exists()
    assert "score_bigwig_values_Coverage" not in outputs

    with outputs["summary"].open() as handle:
        summary = list(csv.DictReader(handle, delimiter="\t"))
    coverage_row = next(row for row in summary if row["comparison"] == "Coverage")
    assert coverage_row["comparison_type"] == "BigWig score comparator"
    assert coverage_row["score_sampling_position"] == "main summit"
    assert coverage_row["matched_pair_count"] == ""
    assert coverage_row["score_observation_count"] == "4"
    assert coverage_row["bigwig_missing_or_nonfinite_values_as_zero"] == "1"

    # BigWig comparators must never appear in distance-dependent summaries.
    with outputs["correlation_by_distance"].open() as handle:
        distance_corr = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["comparison"] for row in distance_corr} == {"DANPOS"}
    with outputs["distance_histogram"].open() as handle:
        histogram = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["comparison"] for row in histogram} == {"DANPOS"}

    import gzip
    with gzip.open(outputs["score_agreement_source_Coverage"], "rt", encoding="utf-8") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert source_rows[0]["score_comparator_type"] == "bigwig"
    assert source_rows[0]["score_sampling_position"] == "main summit"
    assert "absolute_distance" not in source_rows[0]
    metadata = outputs["score_agreement_source_Coverage"].with_name("with_signal_Coverage_score_agreement_metadata.tsv")
    assert metadata.exists()
    assert "detected_plot_type\tcompare-positions-score-signal" in metadata.read_text()


def test_score_bigwig_detail_table_is_opt_in(tmp_path: Path, monkeypatch):
    import sys
    import types

    main = tmp_path / "main.bed"
    compare = tmp_path / "compare.bed"
    bigwig = tmp_path / "coverage.bw"
    write_bed(main, [("chr1", 95, 105, "m1", 1.0), ("chr1", 195, 205, "m2", 2.0)])
    write_bed(compare, [("chr1", 96, 106, "c1", 1.0), ("chr1", 196, 206, "c2", 2.0)])
    bigwig.write_bytes(b"fake")

    class FakeBigWig:
        def isBigWig(self): return True
        def chroms(self): return {"chr1": 1000}
        def values(self, chrom, start, end, numpy=False):
            values = np.ones(end - start, dtype=float) * 5
            return values if numpy else values.tolist()
        def close(self): pass

    monkeypatch.setitem(sys.modules, "pyBigWig", types.SimpleNamespace(open=lambda path: FakeBigWig()))
    args = build_parser().parse_args([
        "--main-bed", str(main), "--compare-bed", str(compare),
        "--score-bigwig", f"Coverage={bigwig}", "--write-detail-tables",
        "--output-prefix", str(tmp_path / "detail_signal"), "--dpi", "40", "--quiet",
    ])
    outputs = run_comparison(args)
    detail = outputs["score_bigwig_values_Coverage"]
    with detail.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["main_summit"]) for row in rows] == [100, 200]
    assert [float(row["bigwig_value"]) for row in rows] == [5.0, 5.0]
