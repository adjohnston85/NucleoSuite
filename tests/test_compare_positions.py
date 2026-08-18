"""Tests for nearest-position and score-agreement analysis."""

from __future__ import annotations

import csv
import heapq
import random
from pathlib import Path

import numpy as np

from nucleosuite.compare_positions import (
    PositionRecord,
    _CandidateCursor,
    _records_by_chrom,
    match_positions,
    match_unique,
    read_positions,
    run_comparison,
)


def write_bed(path: Path, rows: list[tuple[object, ...]]) -> None:
    path.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))


def test_default_midpoint_and_score_column(tmp_path: Path):
    bed = tmp_path / "positions.bed"
    write_bed(
        bed,
        [
            ("chr1", 10, 20, "p1", 3.5),
            ("chr1", 30, 41, "p2", 7),
        ],
    )
    records = read_positions(bed, "A", summit_column=None, score_column=5)
    assert [record.summit for record in records] == [15, 35]
    assert [record.score for record in records] == [3.5, 7.0]


def test_explicit_summit_columns_and_smaller_b_is_query(tmp_path: Path):
    bed_a = tmp_path / "a.bed"
    bed_b = tmp_path / "b.bed"
    write_bed(
        bed_a,
        [
            ("chr1", 90, 111, "a1", 10, ".", 100),
            ("chr1", 190, 211, "a2", 20, ".", 200),
            ("chr1", 290, 311, "a3", 30, ".", 300),
        ],
    )
    write_bed(
        bed_b,
        [
            ("chr1", 196, 217, "b1", 22, ".", 207),
            ("chr1", 302, 323, "b2", 31, ".", 312),
        ],
    )
    records_a = read_positions(bed_a, "A", summit_column=7, score_column=5)
    records_b = read_positions(bed_b, "B", summit_column=7, score_column=5)
    result = match_positions(records_a, records_b, "many-to-one", None)
    assert result.query_source == "B"
    assert [(pair.a.summit, pair.b.summit) for pair in result.pairs] == [
        (200, 207),
        (300, 312),
    ]


def test_many_to_one_allows_target_reuse(tmp_path: Path):
    bed_a = tmp_path / "a.bed"
    bed_b = tmp_path / "b.bed"
    write_bed(
        bed_a,
        [
            ("chr1", 98, 103, "a1", 1),
            ("chr1", 103, 108, "a2", 2),
        ],
    )
    write_bed(
        bed_b,
        [
            ("chr1", 100, 105, "b1", 1),
            ("chr1", 500, 505, "b2", 2),
            ("chr1", 800, 805, "b3", 3),
        ],
    )
    result = match_positions(
        read_positions(bed_a, "A", None, 5),
        read_positions(bed_b, "B", None, 5),
        "many-to-one",
        None,
    )
    assert len(result.pairs) == 2
    assert result.pairs[0].b.summit == result.pairs[1].b.summit


def test_unique_matching_prevents_target_reuse(tmp_path: Path):
    bed_a = tmp_path / "a.bed"
    bed_b = tmp_path / "b.bed"
    write_bed(
        bed_a,
        [
            ("chr1", 98, 103, "a1", 1),
            ("chr1", 103, 108, "a2", 2),
        ],
    )
    write_bed(
        bed_b,
        [
            ("chr1", 100, 105, "b1", 1),
            ("chr1", 108, 113, "b2", 2),
            ("chr1", 500, 505, "b3", 3),
        ],
    )
    result = match_positions(
        read_positions(bed_a, "A", None, 5),
        read_positions(bed_b, "B", None, 5),
        "unique",
        None,
    )
    assert len(result.pairs) == 2
    assert len({pair.b.line_number for pair in result.pairs}) == 2


def _reference_one_to_one(
    queries: list[PositionRecord],
    targets: list[PositionRecord],
    max_distance: float | None,
) -> tuple[list[tuple[int, int]], tuple[int, int, int]]:
    """Small regression oracle for dense one-to-one conflict handling."""
    pairs: list[tuple[int, int]] = []
    no_target = beyond_distance = unique_conflict = 0
    target_by_chrom = _records_by_chrom(targets)
    for chrom, chrom_queries in _records_by_chrom(queries).items():
        chrom_targets = target_by_chrom.get(chrom)
        if not chrom_targets:
            no_target += len(chrom_queries)
            continue
        positions = [record.summit for record in chrom_targets]
        cursors = [
            _CandidateCursor.create(query, chrom_targets, positions)
            for query in chrom_queries
        ]
        candidates: list[tuple[int, int, int]] = []
        for query_index, cursor in enumerate(cursors):
            candidate = cursor.next()
            if candidate is not None:
                distance, target_index = candidate
                heapq.heappush(candidates, (distance, query_index, target_index))
        assigned_queries: set[int] = set()
        assigned_targets: set[int] = set()
        rejected: set[int] = set()
        while candidates:
            distance, query_index, target_index = heapq.heappop(candidates)
            if query_index in assigned_queries:
                continue
            if max_distance is not None and distance > max_distance:
                rejected.add(query_index)
                continue
            if target_index not in assigned_targets:
                assigned_queries.add(query_index)
                assigned_targets.add(target_index)
                pairs.append(
                    (
                        chrom_queries[query_index].line_number,
                        chrom_targets[target_index].line_number,
                    )
                )
                continue
            candidate = cursors[query_index].next()
            if candidate is not None:
                next_distance, next_target = candidate
                heapq.heappush(
                    candidates, (next_distance, query_index, next_target)
                )
        beyond_distance += len(rejected)
        unique_conflict += len(chrom_queries) - len(assigned_queries) - len(rejected)
    return sorted(pairs), (no_target, beyond_distance, unique_conflict)


def _records(summits: list[int], source: str) -> list[PositionRecord]:
    return [
        PositionRecord(source, "chr1", summit, summit + 1, summit, float(index),
                       f"{source}{index}", index)
        for index, summit in enumerate(summits, start=1)
    ]


def test_adjacent_boundary_matcher_preserves_dense_conflict_semantics():
    generator = random.Random(20260812)
    cases = [
        ([0, 4, 4, 8, 12, 12, 16], [2, 2, 6, 10, 14, 14, 18, 100]),
    ]
    for _ in range(250):
        cases.append(
            (
                sorted(generator.randrange(0, 30) for _ in range(generator.randrange(1, 12))),
                sorted(generator.randrange(0, 30) for _ in range(generator.randrange(1, 12))),
            )
        )
    for query_summits, target_summits in cases:
        queries = _records(query_summits, "A")
        targets = _records(target_summits, "B")
        for max_distance in (None, 0, 3, 10):
            expected_pairs, expected_unmatched = _reference_one_to_one(
                queries, targets, max_distance
            )
            observed = match_unique(queries, targets, "A", max_distance)
            observed_pairs = sorted(
                (pair.a.line_number, pair.b.line_number) for pair in observed.pairs
            )
            assert observed_pairs == expected_pairs
            assert (
                observed.unmatched_no_target_chrom,
                observed.unmatched_distance,
                observed.unmatched_unique,
            ) == expected_unmatched


def test_one_to_one_default_and_unique_alias_are_identical():
    from nucleosuite.compare_positions import build_parser

    args = build_parser().parse_args(["--bed-a", "a.bed", "--bed-b", "b.bed"])
    assert args.matching == "one-to-one"
    records_a = _records([1, 3, 7], "A")
    records_b = _records([2, 2, 8, 20], "B")
    one_to_one = match_positions(records_a, records_b, "one-to-one", None)
    unique = match_positions(records_a, records_b, "unique", None)
    assert [(pair.a, pair.b) for pair in one_to_one.pairs] == [
        (pair.a, pair.b) for pair in unique.pairs
    ]


def test_maximum_distance_filters_pairs(tmp_path: Path):
    bed_a = tmp_path / "a.bed"
    bed_b = tmp_path / "b.bed"
    write_bed(bed_a, [("chr1", 0, 10, "a", 1)])
    write_bed(bed_b, [("chr1", 100, 110, "b", 1), ("chr1", 200, 210, "c", 2)])
    result = match_positions(
        read_positions(bed_a, "A", None, 5),
        read_positions(bed_b, "B", None, 5),
        "many-to-one",
        20,
    )
    assert not result.pairs
    assert result.unmatched_distance == 1


def test_full_comparison_writes_tables_and_plots(tmp_path: Path):
    bed_a = tmp_path / "method_a.bed"
    bed_b = tmp_path / "method_b.bed"
    rows_a = []
    rows_b = []
    for index, summit in enumerate(range(100, 1100, 100), start=1):
        rows_a.append(("chr1", summit - 5, summit + 5, f"a{index}", index))
        rows_b.append(("chr1", summit - 3, summit + 7, f"b{index}", index * 2))
    rows_b.extend(
        [
            ("chr1", 1500, 1510, "extra1", 1),
            ("chr1", 1700, 1710, "extra2", 1),
        ]
    )
    write_bed(bed_a, rows_a)
    write_bed(bed_b, rows_b)

    prefix = tmp_path / "comparison"
    args = type(
        "Args",
        (),
        {
            "bed_a": str(bed_a),
            "bed_b": str(bed_b),
            "summit_column_a": None,
            "summit_column_b": None,
            "score_column_a": 5,
            "score_column_b": 5,
            "label_a": "A",
            "label_b": "B",
            "matching": "many-to-one",
            "max_distance": None,
            "distance_bins": "5,10,20",
            "plot_max_points": 100,
            "plot_seed": 1,
            "score_normalization": "zscore",
            "correlation_method": "spearman",
            "histogram_bin_width": 1.0,
            "histogram_x_min": 0.0,
            "histogram_x_max": 300.0,
            "score_z_limit": 10.0,
            "percentile_interval": 10,
            "percentile_boxplot_y_max": 500.0,
            "skip_percentile_distance_analysis": False,
            "dpi": 72,
            "output_prefix": str(prefix),
            "quiet": True,
        },
    )()
    outputs = run_comparison(args)
    assert all(path.exists() and path.stat().st_size > 0 for path in outputs.values())

    with outputs["pairs"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 10
    assert all(float(row["absolute_distance"]) == 2 for row in rows)

    summary = {}
    with outputs["summary"].open() as handle:
        next(handle)
        for line in handle:
            key, value = line.rstrip("\n").split("\t")
            summary[key] = value
    assert summary["score_normalization"] == "zscore"
    assert summary["correlation_method"] == "spearman"
    assert float(summary["spearman_score_correlation"]) == 1.0
    assert "pearson_score_correlation" not in summary
    assert "distance_vs_score_difference_plot" not in outputs
    assert {
        "a_percentiles_vs_all_b_distances",
        "a_percentiles_vs_all_b_summary",
        "a_percentiles_vs_all_b_boxplot",
        "b_percentiles_vs_all_a_distances",
        "b_percentiles_vs_all_a_summary",
        "b_percentiles_vs_all_a_boxplot",
    } <= outputs.keys()

    expected_groups = [
        "0-10", "10-20", "20-30", "30-40", "40-50",
        "50-60", "60-70", "70-80", "80-90", "90-100",
    ]
    with outputs["a_percentiles_vs_all_b_summary"].open() as handle:
        a_percentile_summary = list(csv.DictReader(handle, delimiter="\t"))
    with outputs["b_percentiles_vs_all_a_summary"].open() as handle:
        b_percentile_summary = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["percentile_group"] for row in a_percentile_summary] == expected_groups
    assert [row["percentile_group"] for row in b_percentile_summary] == expected_groups
    assert all(row["percentile_source"] == "A" for row in a_percentile_summary)
    assert all(row["target_source"] == "B" for row in a_percentile_summary)
    assert all(row["percentile_source"] == "B" for row in b_percentile_summary)
    assert all(row["target_source"] == "A" for row in b_percentile_summary)

    with outputs["a_percentiles_vs_all_b_distances"].open() as handle:
        a_direction_rows = list(csv.DictReader(handle, delimiter="\t"))
    with outputs["b_percentiles_vs_all_a_distances"].open() as handle:
        b_direction_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(a_direction_rows) == 10
    assert len(b_direction_rows) == 12

    with outputs["distance_histogram"].open() as handle:
        histogram_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(histogram_rows) == 300
    assert histogram_rows[0]["distance_bin"] == "0-1"
    assert histogram_rows[2]["pair_count"] == "10"
    assert sum(int(row["pair_count"]) for row in histogram_rows) == 10


def test_distance_bin_labels_are_integer_ranges():
    from nucleosuite.compare_positions import _bin_label
    assert _bin_label(None, 5.0) == "0-5"
    assert _bin_label(5.0, 10.0) == "6-10"
    assert _bin_label(10.0, 20.0) == "11-20"
    assert _bin_label(100.0, None) == ">100"


def test_percentile_analysis_uses_source_groups_against_all_opposite_positions(tmp_path: Path):
    from nucleosuite.compare_positions import directional_percentile_distance_rows

    bed_a = tmp_path / "a_percentiles.bed"
    bed_b = tmp_path / "b_percentiles.bed"
    rows_a = []
    rows_b = []
    for index in range(20):
        summit = 1000 + index * 100
        rows_a.append(("chr1", summit - 5, summit + 5, f"a{index + 1}", index + 1))
        # Reverse B's score order. Nearest genomic partners therefore usually
        # belong to a different B percentile group than their A query partner.
        rows_b.append(("chr1", summit - 2, summit + 8, f"b{index + 1}", 1000 + (20 - index) * 100))
    write_bed(bed_a, rows_a)
    write_bed(bed_b, rows_b)

    records_a = read_positions(bed_a, "A", None, 5)
    records_b = read_positions(bed_b, "B", None, 5)
    a_details, a_summaries, a_distances = directional_percentile_distance_rows(
        records_a, records_b, "A", "many-to-one", None, 10
    )
    b_details, b_summaries, b_distances = directional_percentile_distance_rows(
        records_a, records_b, "B", "many-to-one", None, 10
    )

    assert len(a_details) == 20
    assert len(b_details) == 20
    assert {row["percentile_group"] for row in a_details} == {
        "0-10", "10-20", "20-30", "30-40", "40-50",
        "50-60", "60-70", "70-80", "80-90", "90-100",
    }
    assert all(row["absolute_distance"] == 3 for row in a_details + b_details)
    assert all(
        int(row["group_lower_percentile"]) < int(row["a_score_percentile"]) <= int(row["group_upper_percentile"])
        for row in a_details
    )
    assert all(
        int(row["group_lower_percentile"]) < int(row["b_score_percentile"]) <= int(row["group_upper_percentile"])
        for row in b_details
    )
    assert any(
        not (
            int(row["group_lower_percentile"])
            < int(row["b_score_percentile"])
            <= int(row["group_upper_percentile"])
        )
        for row in a_details
    )
    assert any(
        not (
            int(row["group_lower_percentile"])
            < int(row["a_score_percentile"])
            <= int(row["group_upper_percentile"])
        )
        for row in b_details
    )
    assert len(a_summaries) == len(b_summaries) == 10
    assert all(row["matched_pair_count"] == 2 for row in a_summaries + b_summaries)
    assert all(row["all_target_position_count"] == 20 for row in a_summaries + b_summaries)
    assert all(values.size == 2 for values in a_distances.values())
    assert all(values.size == 2 for values in b_distances.values())


def test_default_percentile_groups_are_quartiles():
    from nucleosuite.compare_positions import _percentile_group_bounds, build_parser

    args = build_parser().parse_args(["--bed-a", "a.bed", "--bed-b", "b.bed"])
    assert args.percentile_interval == 25
    assert _percentile_group_bounds(args.percentile_interval) == [
        (0, 25, "0-25"),
        (25, 50, "25-50"),
        (50, 75, "50-75"),
        (75, 100, "75-100"),
    ]


def test_percentile_interval_can_be_changed(tmp_path: Path):
    from nucleosuite.compare_positions import assign_score_percentiles

    bed = tmp_path / "scores.bed"
    write_bed(
        bed,
        [("chr1", index * 10, index * 10 + 5, f"p{index}", index) for index in range(1, 21)],
    )
    assignments = assign_score_percentiles(read_positions(bed, "A", None, 5), 25)
    assert {assignment.label for assignment in assignments.values()} == {
        "0-25", "25-50", "50-75", "75-100"
    }


def test_streaming_percentiles_resume_and_optional_pair_tables(
    tmp_path: Path, monkeypatch
):
    import nucleosuite.compare_positions as comparison

    bed_a = tmp_path / "a.bed"
    bed_b = tmp_path / "b.bed"
    write_bed(
        bed_a,
        [("chr1", index * 10, index * 10 + 3, f"a{index}", index)
         for index in range(1, 11)],
    )
    write_bed(
        bed_b,
        [("chr1", index * 10 + 1, index * 10 + 4, f"b{index}", 20 - index)
         for index in range(1, 13)],
    )
    prefix = tmp_path / "streamed"
    args = comparison.build_parser().parse_args(
        [
            "--bed-a", str(bed_a),
            "--bed-b", str(bed_b),
            "--output-prefix", str(prefix),
            "--percentile-interval", "30",
            "--skip-pairs-tsv",
            "--skip-percentile-pairs-tsv",
            "--dpi", "40",
            "--quiet",
        ]
    )
    outputs = comparison.run_comparison(args)
    assert "pairs" not in outputs
    assert not Path(f"{prefix}_pairs.tsv").exists()
    assert "a_percentiles_vs_all_b_distances" not in outputs
    assert "b_percentiles_vs_all_a_distances" not in outputs
    assert not Path(f"{prefix}_A_percentiles_vs_all_B_distances.tsv").exists()
    with outputs["a_percentiles_vs_all_b_summary"].open() as handle:
        summaries = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["percentile_group"] for row in summaries] == [
        "0-30", "30-60", "60-90", "90-100"
    ]
    markers = sorted(Path(f"{prefix}_checkpoints").rglob("*.complete.json"))
    assert len(markers) == 9  # main plus four groups in each direction

    def fail_if_recomputed(**_kwargs):
        raise AssertionError("a completed percentile group was recomputed")

    monkeypatch.setattr(comparison, "_stream_percentile_group", fail_if_recomputed)
    resumed = comparison.run_comparison(args)
    assert resumed["summary"].exists()


def test_plot_defaults_limit_zscores_histogram_and_percentile_axes(tmp_path: Path, monkeypatch):
    import matplotlib.figure

    from nucleosuite.compare_positions import create_plots, plot_percentile_distances

    captured: list[dict[str, object]] = []

    def capture_savefig(figure, path, *args, **kwargs):
        axis = figure.axes[0]
        captured.append(
            {
                "path": str(path),
                "xlim": tuple(float(value) for value in axis.get_xlim()),
                "ylim": tuple(float(value) for value in axis.get_ylim()),
                "line_count": len(axis.lines),
            }
        )

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", capture_savefig)

    arrays = {
        "scores_a": np.asarray([1.0, 2.0, 3.0]),
        "scores_b": np.asarray([2.0, 4.0, 6.0]),
        "z_a": np.asarray([-0.5, 0.0, 700.0]),
        "z_b": np.asarray([-0.4, 0.1, 750.0]),
        "percentile_a": np.asarray([0.0, 0.5, 1.0]),
        "percentile_b": np.asarray([0.0, 0.5, 1.0]),
        "distances": np.asarray([2.0, 20.0, 450.0]),
    }
    bin_rows = [
        {
            "distance_bin": "0-5",
            "spearman_score_correlation": 1.0,
            "pearson_score_correlation": 1.0,
            "pair_count": 1,
        }
    ]
    create_plots(
        tmp_path / "limits",
        arrays,
        bin_rows,
        "A",
        "B",
        0,
        1,
        72,
        "zscore",
        "spearman",
        1.0,
        0.0,
        300.0,
        10.0,
    )
    score_plot = captured[0]
    histogram_plot = captured[1]
    assert score_plot["xlim"] == (-10.0, 10.0)
    assert score_plot["ylim"] == (-10.0, 10.0)
    assert score_plot["line_count"] == 0
    assert histogram_plot["xlim"] == (0.0, 300.0)

    captured.clear()
    plot_percentile_distances(
        tmp_path / "percentiles.png",
        {"0-10": np.asarray([10.0, 600.0])},
        10,
        "A",
        "B",
        72,
        500.0,
    )
    assert captured[0]["ylim"] == (0.0, 500.0)
