"""Tests for WW-type output defaults and streamed coordinate ordering."""

from __future__ import annotations

import argparse

import pytest

from nucleosuite.cli.main import build_parser
from nucleosuite.workflows.ww_types import _ensure_coordinate_order


def _ww_parser():
    parser = build_parser()
    action = next(item for item in parser._actions if isinstance(item, argparse._SubParsersAction))
    return action.choices["ww-types"]


def test_ww_types_dinucleotide_profiles_are_enabled_by_default():
    parser = _ww_parser()
    args = parser.parse_args(["--bam", "sample.bam", "--fasta", "genome.fa"])
    assert args.dinuc_profile is True
    args = parser.parse_args(["--bam", "sample.bam", "--fasta", "genome.fa", "--no-dinuc-profile"])
    assert args.dinuc_profile is False


def test_ww_type_chunk_is_sorted_only_when_an_inversion_is_present():
    ordered = [(10, 20, "type1"), (20, 30, "type2")]
    assert _ensure_coordinate_order(ordered) is ordered
    inverted = [(20, 30, "type2"), (10, 20, "type1"), (10, 18, "type4")]
    assert _ensure_coordinate_order(inverted) == [
        (10, 18, "type4"),
        (10, 20, "type1"),
        (20, 30, "type2"),
    ]


def test_ww_type_length_summary_and_stacked_plot(tmp_path):
    from collections import Counter
    import csv

    from nucleosuite.profile_plots import plot_ww_type_length_stacked
    from nucleosuite.sequence.ww_types import write_length_summary_table

    table = tmp_path / "types_by_length.tsv"
    plot = tmp_path / "types_by_length_stacked.png"
    write_length_summary_table(
        table,
        {
            145: Counter(type1=40, type2=30, type3=20, type4=10, unclassified=5),
            161: Counter(type1=25, type2=25, type3=25, type4=25),
        },
    )
    with table.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["fragment_length"]) for row in rows] == [145, 161]
    assert int(rows[0]["total_fragments"]) == 105
    assert int(rows[0]["classified_fragments"]) == 100
    assert float(rows[0]["type1_percent_of_classified"]) == 40.0
    assert float(rows[0]["unclassified_percent_of_all"]) == pytest.approx(100 * 5 / 105)

    plot_ww_type_length_stacked(table, plot, title="test")
    assert plot.is_file()
    assert plot.stat().st_size > 0


def test_ww_type_length_combination_recalculates_frequencies(tmp_path):
    from collections import Counter
    import csv

    from nucleosuite.combine import _combine_ww_type_by_length
    from nucleosuite.sequence.ww_types import write_length_summary_table

    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    output = tmp_path / "combined_ww_type_by_length.tsv"
    write_length_summary_table(first, {145: Counter(type1=8, type2=2)})
    write_length_summary_table(second, {145: Counter(type1=2, type2=8)})
    _combine_ww_type_by_length([first, second], output)

    with output.open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert int(row["type1_count"]) == 10
    assert int(row["type2_count"]) == 10
    assert float(row["type1_percent_of_classified"]) == 50.0
    assert output.with_name(output.stem + "_stacked.png").is_file()


def test_selected_ww_type_lengths_are_extracted_without_recounting(tmp_path):
    from collections import Counter
    import csv

    from nucleosuite.sequence.ww_types import (
        write_length_summary_table,
        write_selected_length_summary,
    )

    first = tmp_path / "144_146.tsv"
    second = tmp_path / "160_162.tsv"
    third = tmp_path / "166_168.tsv"
    output = tmp_path / "exact.tsv"
    write_length_summary_table(
        first,
        {
            144: Counter(type1=1),
            145: Counter(type1=4, type2=6),
            146: Counter(type2=1),
        },
    )
    write_length_summary_table(second, {161: Counter(type3=7, type4=3)})
    write_length_summary_table(third, {167: Counter(type1=2, type4=8)})

    write_selected_length_summary([first, second, third], [145, 161, 167], output)
    with output.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["fragment_length"]) for row in rows] == [145, 161, 167]
    assert [int(row["total_fragments"]) for row in rows] == [10, 10, 10]
    assert float(rows[0]["type2_percent_of_classified"]) == 60.0


def test_classify_core_sequence_rejects_missing_or_invalid_core():
    from nucleosuite.sequence.ww_types import classify_core_sequence

    assert classify_core_sequence(None) is None
    assert classify_core_sequence("A" * 146) is None
    assert classify_core_sequence("A" * 146 + "N") is None


def test_classify_core_sequence_executes_all_four_classification_branches(monkeypatch):
    from nucleosuite.sequence import ww_types

    cases = {
        "type1": [(36, 0), (32, 0)],
        "type2": [(36, 37), (32, 32)],
        "type3": [(0, 0), (32, 32)],
        "type4": [(0, 37), (32, 32)],
    }
    for expected, counts in cases.items():
        values = iter(counts)
        monkeypatch.setattr(
            ww_types,
            "count_ww_ss_at_sites",
            lambda _core, _sites, values=values: next(values),
        )
        assert ww_types.classify_core_sequence("A" * 147) == expected


def test_classify_fragment_uses_shared_core_classifier(monkeypatch):
    from nucleosuite.sequence import ww_types

    monkeypatch.setattr(ww_types, "extract_reference_sequence", lambda **_kwargs: "A" * 147)
    calls = []

    def classify(core):
        calls.append(core)
        return "type1"

    monkeypatch.setattr(ww_types, "classify_core_sequence", classify)
    assert ww_types.classify_fragment(object(), object(), 10, 155) == "type1"
    assert calls == ["A" * 147]
