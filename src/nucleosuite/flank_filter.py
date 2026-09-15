#!/usr/bin/env python3
"""Filter BED regions to records flanked by features on both sides."""

from __future__ import annotations

import argparse
import bisect
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, TextIO

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.io import open_text


HEADER_PREFIXES = ("#", "track", "browser")


@dataclass(frozen=True)
class BedRecord:
    """A BED row with its representative genomic position."""

    chrom: str
    start: int
    end: int
    position: int
    raw_line: str


def _split_fields(raw: str) -> list[str]:
    fields = raw.rstrip("\n").split("\t")
    return fields if len(fields) > 1 else raw.split()


def _parse_position(
    fields: Sequence[str],
    *,
    start: int,
    end: int,
    position_column: int | None,
    path: str | Path,
    line_number: int,
) -> int:
    """Return the requested position or the NucleoSuite BED8 centre when available."""

    if position_column is not None:
        if position_column < 1:
            raise ValueError("Position columns must be one-based integers of at least 1")
        index = position_column - 1
        if index >= len(fields):
            raise ValueError(
                f"{path}:{line_number}: requested position column {position_column} is missing"
            )
        try:
            return int(fields[index])
        except ValueError as error:
            raise ValueError(
                f"{path}:{line_number}: position column {position_column} must contain an integer"
            ) from error

    # NucleoSuite peak BED8 uses thickStart (BED column 7) as the representative
    # call centre. Prefer it automatically when it is present and numeric.
    if len(fields) >= 7:
        try:
            return int(fields[6])
        except ValueError:
            pass

    return (start + end) // 2


def read_bed_records(
    path: str | Path,
    *,
    position_column: int | None = None,
) -> list[BedRecord]:
    """Read BED3+ records and assign one representative position to each row."""

    records: list[BedRecord] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith(HEADER_PREFIXES):
                continue
            fields = _split_fields(raw)
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least BED3")
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: BED start/end must be integers"
                ) from error
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: invalid half-open BED interval {start}-{end}"
                )
            position = _parse_position(
                fields,
                start=start,
                end=end,
                position_column=position_column,
                path=path,
                line_number=line_number,
            )
            if position < 0:
                raise ValueError(f"{path}:{line_number}: representative position must be non-negative")
            records.append(
                BedRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    position=position,
                    raw_line=raw.rstrip("\n"),
                )
            )
    return records


def build_flank_index(records: Sequence[BedRecord]) -> dict[str, list[int]]:
    """Return sorted representative flank positions for each chromosome."""

    by_chrom: dict[str, list[int]] = {}
    for record in records:
        by_chrom.setdefault(record.chrom, []).append(record.position)
    for positions in by_chrom.values():
        positions.sort()
    return by_chrom


def filter_flanked_regions(
    regions: Sequence[BedRecord],
    flank_positions: dict[str, list[int]],
    *,
    max_flank_distance: int | None = None,
) -> list[BedRecord]:
    """Keep only clean FLANK--REGION--FLANK arrangements.

    The nearest strictly upstream and downstream flank positions define the
    bracketing interval for each candidate region. Exactly one ``--regions``
    record must occur inside that interval: the candidate itself. This rejects
    consecutive region calls such as FLANK--REGION--REGION--FLANK while still
    allowing consecutive flank calls.
    """

    if max_flank_distance is not None and max_flank_distance < 0:
        raise ValueError("--max-flank-distance must be at least 0")

    # Index all candidate-region positions as well as the flank positions.
    # This lets us enforce the adjacency rule without merging or rewriting the
    # original BED records. Duplicate region positions count as multiple calls.
    region_positions_by_chrom = build_flank_index(regions)

    retained: list[BedRecord] = []
    for region in regions:
        positions = flank_positions.get(region.chrom)
        if not positions:
            continue

        # Flanks must be strictly on opposite sides of the candidate position.
        # A flank at exactly the region position does not count for either side.
        left_index = bisect.bisect_left(positions, region.position) - 1
        right_index = bisect.bisect_right(positions, region.position)
        if left_index < 0 or right_index >= len(positions):
            continue

        left = positions[left_index]
        right = positions[right_index]

        # A clean flanked call contains exactly one region between the nearest
        # bracketing flanks. Thus BRK--NUC--BRK is retained, whereas
        # BRK--NUC--NUC--BRK rejects both nucleosome calls.
        region_positions = region_positions_by_chrom[region.chrom]
        first_inside = bisect.bisect_right(region_positions, left)
        first_at_or_after_right = bisect.bisect_left(region_positions, right)
        if first_at_or_after_right - first_inside != 1:
            continue

        if max_flank_distance is not None:
            if region.position - left > max_flank_distance:
                continue
            if right - region.position > max_flank_distance:
                continue
        retained.append(region)
    return retained


def _open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("w", encoding="utf-8")


def write_bed(records: Sequence[BedRecord], output: str | Path) -> Path:
    """Write retained records without altering their BED columns."""

    path = Path(output)
    with _open_output(path) as handle:
        for record in records:
            handle.write(record.raw_line + "\n")
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite flank-filter",
        description=(
            "Retain BED regions only when the nearest surrounding events form a clean "
            "flank--region--flank arrangement with no second region between the two flanks."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument(
        "--regions",
        required=True,
        help=(
            "BED3+, BED.gz, or bigBed records to filter. Retained rows are written "
            "unchanged to the output."
        ),
    )
    parser.add_argument(
        "--flanks",
        required=True,
        help=(
            "BED3+, BED.gz, or bigBed features used as the required immediate "
            "upstream/downstream flank type."
        ),
    )
    parser.add_argument(
        "--out",
        "--output",
        dest="output",
        required=True,
        help="Filtered BED output path; .gz writes gzip-compressed BED text.",
    )
    parser.add_argument(
        "--region-position-column",
        type=int,
        default=None,
        help=(
            "One-based absolute position column for --regions. By default, BED column 7 "
            "is used when present and numeric; otherwise the interval midpoint is used."
        ),
    )
    parser.add_argument(
        "--flank-position-column",
        type=int,
        default=None,
        help=(
            "One-based absolute position column for --flanks. By default, BED column 7 "
            "is used when present and numeric; otherwise the interval midpoint is used."
        ),
    )
    parser.add_argument(
        "--max-flank-distance",
        type=int,
        default=None,
        help=(
            "Optional maximum distance in bp from the region position to each nearest "
            "flanking feature; omit for no distance limit."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.region_position_column is not None and args.region_position_column < 1:
        raise ValueError("--region-position-column must be at least 1")
    if args.flank_position_column is not None and args.flank_position_column < 1:
        raise ValueError("--flank-position-column must be at least 1")
    if args.max_flank_distance is not None and args.max_flank_distance < 0:
        raise ValueError("--max-flank-distance must be at least 0")

    regions = read_bed_records(
        args.regions,
        position_column=args.region_position_column,
    )
    flanks = read_bed_records(
        args.flanks,
        position_column=args.flank_position_column,
    )
    retained = filter_flanked_regions(
        regions,
        build_flank_index(flanks),
        max_flank_distance=args.max_flank_distance,
    )
    output = write_bed(retained, args.output)

    print(f"Regions read: {len(regions):,}")
    print(f"Flank features read: {len(flanks):,}")
    print(f"Regions retained: {len(retained):,}")
    print(f"Wrote: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
