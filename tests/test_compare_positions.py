"""Tests for multi-callset nearest-position comparison."""

from __future__ import annotations

import csv
import heapq
import random
from pathlib import Path

import numpy as np

from nucleosuite.compare_positions import (
    PositionRecord,
    _CandidateCursor,
    _pairwise_statistics,
    _percentile_group_bounds,
    _records_by_chrom,
    build_parser,
    match_positions,
    match_unique,
    read_positions,
    run_comparison,
)


def write_bed(path: Path, rows: list[tuple[object, ...]]) -> None:
    path.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))


def test_default_midpoint_and_score_column(tmp_path: Path):
    bed = tmp_path / "positions.bed"
    write_bed(bed, [("chr1", 10, 20, "p1", 3.5), ("chr1", 30, 41, "p2", 7)])
    records = read_positions(bed, "A", summit_column=None, score_column=5)
    assert [record.summit for record in records] == [15, 35]
    assert [record.score for record in records] == [3.5, 7.0]


def test_smaller_comparison_is_query_but_pairs_keep_main_as_a(tmp_path: Path):
    main = tmp_path / "main.bed"
    compare = tmp_path / "compare.bed"
    write_bed(main, [("chr1", 95, 106, "m1", 10, ".", 100), ("chr1", 195, 206, "m2", 20, ".", 200), ("chr1", 295, 306, "m3", 30, ".", 300)])
    write_bed(compare, [("chr1", 201, 212, "c1", 22, ".", 207), ("chr1", 306, 317, "c2", 31, ".", 312)])
    records_main = read_positions(main, "A", summit_column=7, score_column=5)
    records_compare = read_positions(compare, "B", summit_column=7, score_column=5)
    result = match_positions(records_main, records_compare, "one-to-one", None)
    assert result.query_source == "B"
    assert [(pair.a.summit, pair.b.summit) for pair in result.pairs] == [(200, 207), (300, 312)]


def _reference_one_to_one(queries, targets, max_distance):
    pairs = []
    no_target = beyond_distance = unique_conflict = 0
    target_by_chrom = _records_by_chrom(targets)
    for chrom, chrom_queries in _records_by_chrom(queries).items():
        chrom_targets = target_by_chrom.get(chrom)
        if not chrom_targets:
            no_target += len(chrom_queries)
            continue
        positions = [record.summit for record in chrom_targets]
        cursors = [_CandidateCursor.create(query, chrom_targets, positions) for query in chrom_queries]
        candidates = []
        for qi, cursor in enumerate(cursors):
            candidate = cursor.next()
            if candidate is not None:
                distance, ti = candidate
                heapq.heappush(candidates, (distance, qi, ti))
        assigned_q = set(); assigned_t = set(); rejected = set()
        while candidates:
            distance, qi, ti = heapq.heappop(candidates)
            if qi in assigned_q:
                continue
            if max_distance is not None and distance > max_distance:
                rejected.add(qi); continue
            if ti not in assigned_t:
                assigned_q.add(qi); assigned_t.add(ti)
                pairs.append((chrom_queries[qi].line_number, chrom_targets[ti].line_number))
                continue
            candidate = cursors[qi].next()
            if candidate is not None:
                next_distance, next_target = candidate
                heapq.heappush(candidates, (next_distance, qi, next_target))
        beyond_distance += len(rejected)
        unique_conflict += len(chrom_queries) - len(assigned_q) - len(rejected)
    return sorted(pairs), (no_target, beyond_distance, unique_conflict)


def _records(summits: list[int], source: str) -> list[PositionRecord]:
    return [PositionRecord(source, "chr1", summit, summit + 1, summit, float(index), f"{source}{index}", index) for index, summit in enumerate(summits, start=1)]


def test_one_to_one_matching_preserves_dense_conflict_semantics():
    generator = random.Random(20260819)
    cases = [([0, 4, 4, 8, 12, 12, 16], [2, 2, 6, 10, 14, 14, 18, 100])]
    for _ in range(80):
        cases.append((sorted(generator.randrange(0, 30) for _ in range(generator.randrange(1, 10))), sorted(generator.randrange(0, 30) for _ in range(generator.randrange(1, 10)))))
    for query_summits, target_summits in cases:
        queries = _records(query_summits, "A"); targets = _records(target_summits, "B")
        for max_distance in (None, 0, 3, 10):
            expected_pairs, expected_unmatched = _reference_one_to_one(queries, targets, max_distance)
            observed = match_unique(queries, targets, "A", max_distance)
            assert sorted((pair.a.line_number, pair.b.line_number) for pair in observed.pairs) == expected_pairs
            assert (observed.unmatched_no_target_chrom, observed.unmatched_distance, observed.unmatched_unique) == expected_unmatched


def test_new_parser_defaults_to_main_plus_repeated_compare_and_quartiles():
    args = build_parser().parse_args(["--main-bed", "main.bed", "--compare-bed", "B=b.bed", "--compare-bed", "C=c.bed"])
    assert args.main_bed == "main.bed"
    assert args.compare_beds == ["B=b.bed", "C=c.bed"]
    assert args.percentile_interval == 25
    assert args.stats_test == "nonparametric"
    assert args.p_adjust == "holm"
    assert args.score_distance_type == "absolute"
    assert args.score_distance_correlation == "spearman"
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
    assert outputs["percentile_distances"].exists()
    assert outputs["percentile_boxplot"].exists()
    assert outputs["score_distance_statistics"].exists()
    assert outputs["pairs_B"].exists() and outputs["pairs_C"].exists()
    with outputs["summary"].open() as handle:
        summary = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["comparison"] for row in summary] == ["B", "C"]
    assert summary[0]["query_source"] == "main"
    assert summary[1]["query_source"] == "comparison"
    with outputs["percentile_distances"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["comparison"] for row in rows} == {"B", "C"}
    assert {row["percentile_group"] for row in rows} == {"0-25", "25-50", "50-75", "75-100"}
    # Main scores rise with genomic order, so the first quartile contains the lowest main scores.
    b_first = [row for row in rows if row["comparison"] == "B" and row["percentile_group"] == "0-25"]
    b_last = [row for row in rows if row["comparison"] == "B" and row["percentile_group"] == "75-100"]
    assert max(float(row["main_score"]) for row in b_first) < min(float(row["main_score"]) for row in b_last)
    for plot_key in ("distance_histogram_plot", "correlation_by_distance_plot", "percentile_boxplot", "score_agreement_plot_B", "score_agreement_plot_C", "score_distance_plot_B", "score_distance_plot_C"):
        plot = outputs[plot_key]
        assert plot.exists() and plot.stat().st_size > 0
        assert plot.with_name(plot.stem + "_metadata.tsv").exists()


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


def test_score_distance_statistics_report_spearman_pearson_and_linear_fit(tmp_path: Path):
    main = tmp_path / "main.bed"; compare = tmp_path / "compare.bed"
    main_rows=[]; compare_rows=[]
    # Distance grows monotonically with main score, giving Spearman rho = 1.
    for index in range(1, 11):
        summit = index * 1000
        main_rows.append(("chr1", summit-5, summit+5, f"m{index}", index))
        shift = index
        compare_rows.append(("chr1", summit-5+shift, summit+5+shift, f"c{index}", index*2))
    write_bed(main, main_rows); write_bed(compare, compare_rows)
    args = build_parser().parse_args([
        "--main-bed", str(main), "--compare-bed", f"C={compare}", "--output-prefix", str(tmp_path/'corr'),
        "--score-distance-type", "absolute", "--dpi", "40", "--quiet",
    ])
    outputs = run_comparison(args)
    with outputs["score_distance_statistics"].open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert np.isclose(float(row["spearman_rho"]), 1.0)
    assert float(row["spearman_p_value"]) < 1e-6
    assert np.isclose(float(row["pearson_r"]), 1.0)
    assert np.isclose(float(row["linear_r_squared"]), 1.0)
    assert np.isclose(float(row["linear_slope_bp_per_score_unit"]), 1.0)


def test_legacy_two_file_alias_is_accepted(tmp_path: Path):
    main = tmp_path / "main.bed"; compare = tmp_path / "compare.bed"
    write_bed(main, [("chr1", 0, 10, "m", 1)])
    write_bed(compare, [("chr1", 1, 11, "c", 1)])
    args = build_parser().parse_args(["--bed-a", str(main), "--bed-b", str(compare), "--output-prefix", str(tmp_path/'legacy'), "--dpi", "40", "--quiet"])
    outputs = run_comparison(args)
    assert outputs["summary"].exists()
