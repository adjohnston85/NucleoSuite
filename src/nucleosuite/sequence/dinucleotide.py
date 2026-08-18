"""Observed dinucleotide profiles aligned to fragment dyads."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import List, Optional, Tuple

from nucleosuite.core.randomization import DINUCS
from nucleosuite.core.reference import sequence_is_acgt

WW_DINUCS = frozenset({"AA", "AT", "TA", "TT"})
SS_DINUCS = frozenset({"CC", "CG", "GC", "GG"})


def fragment_dyad(fragment_start: int, fragment_end: int) -> int:
    """Return the right-hand central base for even-length fragments."""
    return fragment_start + ((fragment_end - fragment_start) // 2)


def expected_profile_positions(fragment_lower: int, fragment_upper: int) -> List[int]:
    positions: set[int] = set()
    for fragment_length in range(fragment_lower, fragment_upper + 1):
        dyad_offset = fragment_length // 2
        for index in range(fragment_length - 1):
            positions.add(index - dyad_offset)
    return sorted(positions)


def new_accumulator():
    return {
        "counts": defaultdict(Counter),
        "n_valid": Counter(),
        "fragments_used": 0,
        "fragments_skipped": 0,
    }


def add_fragment(
    accumulator,
    fragment_sequence: Optional[str],
    fragment_start: int,
    fragment_end: int,
) -> bool:
    """Add one observed fragment sequence to a dyad-aligned profile."""
    expected_length = fragment_end - fragment_start
    if fragment_sequence is None or len(fragment_sequence) != expected_length:
        accumulator["fragments_skipped"] += 1
        return False
    if not sequence_is_acgt(fragment_sequence):
        accumulator["fragments_skipped"] += 1
        return False

    dyad = fragment_dyad(fragment_start, fragment_end)
    for index in range(len(fragment_sequence) - 1):
        dinuc = fragment_sequence[index : index + 2]
        relative_position = (fragment_start + index) - dyad
        accumulator["counts"][relative_position][dinuc] += 1
        accumulator["n_valid"][relative_position] += 1
    accumulator["fragments_used"] += 1
    return True


def write_counts(
    output_path: str,
    accumulator,
    positions: List[int],
) -> None:
    """Write exact dinucleotide counts used to combine per-contig profiles."""
    header = ["position", "n_valid"]
    header.extend(f"{dinuc}_count" for dinuc in DINUCS)
    header.extend(("fragments_used", "fragments_skipped"))
    with open(output_path, "w", encoding="utf-8") as output:
        output.write("\t".join(header) + "\n")
        for row_index, relative_position in enumerate(positions):
            counts = accumulator["counts"].get(relative_position, Counter())
            row = [
                str(relative_position),
                str(int(accumulator["n_valid"].get(relative_position, 0))),
            ]
            row.extend(str(int(counts[dinuc])) for dinuc in DINUCS)
            row.extend((
                str(int(accumulator["fragments_used"])) if row_index == 0 else "0",
                str(int(accumulator["fragments_skipped"])) if row_index == 0 else "0",
            ))
            output.write("\t".join(row) + "\n")


def write_profile(
    output_path: str,
    accumulator,
    positions: List[int],
    fraction: bool = False,
    write_count_table: bool = True,
) -> None:
    multiplier = 1.0 if fraction else 100.0
    suffix = "frac" if fraction else "pct"
    header = ["position", "n_valid"]
    header.extend(f"{dinuc}_{suffix}" for dinuc in DINUCS)
    header.extend((f"WW_{suffix}", f"SS_{suffix}"))

    with open(output_path, "w", encoding="utf-8") as output:
        output.write("\t".join(header) + "\n")
        for relative_position in positions:
            n_valid = int(accumulator["n_valid"].get(relative_position, 0))
            counts = accumulator["counts"].get(relative_position, Counter())
            row = [str(relative_position), str(n_valid)]
            if n_valid == 0:
                row.extend(["NaN"] * (len(DINUCS) + 2))
            else:
                row.extend(
                    f"{(counts[dinuc] / n_valid) * multiplier:.8g}"
                    for dinuc in DINUCS
                )
                ww_count = sum(counts[dinuc] for dinuc in WW_DINUCS)
                ss_count = sum(counts[dinuc] for dinuc in SS_DINUCS)
                row.extend(
                    (
                        f"{(ww_count / n_valid) * multiplier:.8g}",
                        f"{(ss_count / n_valid) * multiplier:.8g}",
                    )
                )
            output.write("\t".join(row) + "\n")

    if write_count_table:
        counts_path = (
            output_path[:-4] + "_counts.tsv"
            if output_path.endswith(".tsv")
            else output_path + "_counts.tsv"
        )
        write_counts(counts_path, accumulator, positions)
