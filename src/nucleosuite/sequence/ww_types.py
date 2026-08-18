"""WW/SS groove-site classification using a centred 147-bp core."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Optional, Sequence, Tuple

try:
    import pysam
except ImportError:  # sequence classification requires pysam only when executed
    pysam = None

from nucleosuite.core.reference import (
    ReferenceContext,
    extract_reference_sequence,
    sequence_is_acgt,
)
from nucleosuite.sequence.dinucleotide import SS_DINUCS, WW_DINUCS, fragment_dyad

WW_TYPE_GROUPS = ("type1", "type2", "type3", "type4")
ALL_OUTPUT_GROUPS = ("all",) + WW_TYPE_GROUPS
WW_MAJOR_SITE_COEFFICIENT = 36.0 / 32.0

MINOR_GBS_RELATIVE = (
    (-69, -66), (-59, -56), (-48, -45), (-37, -34),
    (-27, -24), (-17, -14), (14, 17), (24, 27),
    (34, 37), (45, 48), (56, 59), (66, 69),
)
MAJOR_GBS_RELATIVE = (
    (-64, -61), (-53, -51), (-42, -40), (-32, -29),
    (-22, -19), (-12, -9), (9, 12), (19, 22),
    (29, 32), (40, 42), (51, 53), (61, 64),
)


def _dinuc_starts_from_base_ranges(base_ranges):
    return tuple(
        position
        for site_start, site_end in base_ranges
        for position in range(site_start, site_end)
    )


MINOR_GBS_DINUC_STARTS = _dinuc_starts_from_base_ranges(MINOR_GBS_RELATIVE)
MAJOR_GBS_DINUC_STARTS = _dinuc_starts_from_base_ranges(MAJOR_GBS_RELATIVE)

if len(MINOR_GBS_DINUC_STARTS) != 36:
    raise RuntimeError("Expected 36 minor-groove dinucleotide positions")
if len(MAJOR_GBS_DINUC_STARTS) != 32:
    raise RuntimeError("Expected 32 major-groove dinucleotide positions")


def count_ww_ss_at_sites(core_sequence: str, relative_starts) -> Tuple[int, int]:
    ww = 0
    ss = 0
    for relative_start in relative_starts:
        index = relative_start + 73
        dinuc = core_sequence[index : index + 2]
        if dinuc in WW_DINUCS:
            ww += 1
        elif dinuc in SS_DINUCS:
            ss += 1
    return ww, ss




def classify_core_sequence(core: str | None) -> Optional[str]:
    """Classify a pre-fetched 147-bp dyad-centred core sequence."""
    if core is None or len(core) != 147 or not sequence_is_acgt(core):
        return None

    minor_ww, minor_ss = count_ww_ss_at_sites(core, MINOR_GBS_DINUC_STARTS)
    major_ww, major_ss = count_ww_ss_at_sites(core, MAJOR_GBS_DINUC_STARTS)

    ww_minor_enriched = minor_ww >= major_ww * WW_MAJOR_SITE_COEFFICIENT
    ss_minor_enriched = minor_ss > major_ss * WW_MAJOR_SITE_COEFFICIENT

    if ww_minor_enriched and not ss_minor_enriched:
        return "type1"
    if ww_minor_enriched and ss_minor_enriched:
        return "type2"
    if not ww_minor_enriched and not ss_minor_enriched:
        return "type3"
    return "type4"

def classify_fragment(
    fasta: pysam.FastaFile,
    reference_context: ReferenceContext,
    fragment_start: int,
    fragment_end: int,
) -> Optional[str]:
    """Classify a fragment using the reference 147-bp core centred on its dyad."""
    dyad = fragment_dyad(fragment_start, fragment_end)
    core = extract_reference_sequence(
        fasta=fasta,
        context=reference_context,
        seq_start=dyad - 73,
        seq_end=dyad + 74,
    )
    return classify_core_sequence(core)


def group_output_prefix(output_prefix: str, group: str) -> str:
    return output_prefix if group == "all" else f"{output_prefix}_{group}"


def write_summary(
    output_prefix: str,
    type_counts: Counter,
    total_in_range: int,
) -> str:
    path = f"{output_prefix}_ww_type_summary.tsv"
    classified_total = sum(type_counts[group] for group in WW_TYPE_GROUPS)
    with open(path, "w", encoding="utf-8") as output:
        output.write(
            "type\tfragment_count\tpercent_of_all_in_range"
            "\tpercent_of_classified\n"
        )
        for group in WW_TYPE_GROUPS:
            count = int(type_counts[group])
            percent_all = 100.0 * count / total_in_range if total_in_range else float("nan")
            percent_classified = (
                100.0 * count / classified_total
                if classified_total
                else float("nan")
            )
            output.write(
                f"{group}\t{count}\t{percent_all:.8g}\t{percent_classified:.8g}\n"
            )
        unclassified = int(type_counts["unclassified"])
        percent_all = (
            100.0 * unclassified / total_in_range
            if total_in_range
            else float("nan")
        )
        output.write(f"unclassified\t{unclassified}\t{percent_all:.8g}\tNaN\n")
        output.write(f"all\t{int(total_in_range)}\t100\tNaN\n")
    return path


def _percentage(numerator: int, denominator: int) -> float:
    return 100.0 * numerator / denominator if denominator else float("nan")


def write_length_summary_table(
    output_path: str | Path,
    type_counts_by_length: Mapping[int, Counter],
) -> str:
    """Write WW/SS type counts and relative frequencies by fragment length."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "fragment_length",
        "total_fragments",
        "classified_fragments",
        *[f"{group}_count" for group in WW_TYPE_GROUPS],
        "unclassified_count",
        *[f"{group}_percent_of_all" for group in WW_TYPE_GROUPS],
        "unclassified_percent_of_all",
        *[f"{group}_percent_of_classified" for group in WW_TYPE_GROUPS],
    ]
    with path.open("w", encoding="utf-8") as output:
        output.write("\t".join(fieldnames) + "\n")
        for fragment_length in sorted(type_counts_by_length):
            counts = type_counts_by_length[fragment_length]
            classified = sum(int(counts[group]) for group in WW_TYPE_GROUPS)
            unclassified = int(counts["unclassified"])
            total = classified + unclassified
            values: list[str] = [
                str(int(fragment_length)),
                str(total),
                str(classified),
            ]
            values.extend(str(int(counts[group])) for group in WW_TYPE_GROUPS)
            values.append(str(unclassified))
            values.extend(
                f"{_percentage(int(counts[group]), total):.8g}"
                for group in WW_TYPE_GROUPS
            )
            values.append(f"{_percentage(unclassified, total):.8g}")
            values.extend(
                f"{_percentage(int(counts[group]), classified):.8g}"
                for group in WW_TYPE_GROUPS
            )
            output.write("\t".join(values) + "\n")
    return str(path)


def write_length_summary(
    output_prefix: str,
    type_counts_by_length: Mapping[int, Counter],
) -> str:
    return write_length_summary_table(
        f"{output_prefix}_ww_type_by_length.tsv",
        type_counts_by_length,
    )


def write_selected_length_summary(
    input_tables: Sequence[str | Path],
    selected_lengths: Sequence[int],
    output_path: str | Path,
) -> str:
    """Select non-overlapping fragment lengths from by-length summary tables."""

    selected = [int(value) for value in selected_lengths]
    selected_set = set(selected)
    counts_by_length: dict[int, Counter] = defaultdict(Counter)
    source_by_length: dict[int, str] = {}
    for table_value in input_tables:
        table = Path(table_value)
        with table.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                fragment_length = int(float(row["fragment_length"]))
                if fragment_length not in selected_set:
                    continue
                if fragment_length in source_by_length:
                    raise ValueError(
                        f"Fragment length {fragment_length} occurs in both "
                        f"{source_by_length[fragment_length]} and {table}; "
                        "source range classes must not overlap"
                    )
                source_by_length[fragment_length] = str(table)
                for group in WW_TYPE_GROUPS:
                    counts_by_length[fragment_length][group] = int(
                        float(row[f"{group}_count"])
                    )
                counts_by_length[fragment_length]["unclassified"] = int(
                    float(row["unclassified_count"])
                )
    missing = [value for value in selected if value not in counts_by_length]
    if missing:
        raise ValueError(
            "No WW/SS type-frequency row was found for fragment length(s): "
            + ", ".join(map(str, missing))
        )
    ordered = {value: counts_by_length[value] for value in selected}
    return write_length_summary_table(output_path, ordered)
