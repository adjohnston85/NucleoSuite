"""Common peak-record writers."""

from __future__ import annotations

import os
from typing import Iterable, Mapping


def bed_score(value: float, scale: float = 1.0) -> int:
    """Convert a signal score into a valid non-negative BED score."""
    return max(0, min(1000, int(round(abs(float(value)) * scale))))


def write_bed8_records(
    path: str,
    records: Iterable[Mapping],
    label: str,
    score_scale: float = 1.0,
    mode: str = "a",
    *,
    integer_scores: bool = True,
) -> None:
    with open(path, mode, encoding="utf-8") as output:
        for record in records:
            centre = int(record["region_centre"])
            if integer_scores:
                score = str(bed_score(float(record["peak_score"]), score_scale))
            else:
                score = f"{abs(float(record['peak_score'])) * score_scale:.6f}"
            name = f'{record["chrom"]}:{centre}_{label}'
            output.write(
                f'{record["chrom"]}\t{int(record["region_start"])}\t'
                f'{int(record["region_end"])}\t{name}\t{score}\t.\t'
                f'{centre}\t{centre + 1}\n'
            )


def prepare_output(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    for candidate in (path, os.path.splitext(path)[0] + ".bb"):
        if os.path.exists(candidate):
            os.remove(candidate)
