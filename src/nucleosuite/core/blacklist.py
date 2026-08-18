"""BED blacklist parsing and interval-overlap helpers."""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.io import open_text


@dataclass(frozen=True)
class BlacklistSummary:
    path: str
    interval_count: int
    blacklisted_bases: int
    ignored_interval_count: int = 0


class BlacklistIndex:
    """Merged half-open blacklist intervals indexed by canonical contig."""

    def __init__(self, intervals: Mapping[str, Sequence[tuple[int, int]]], *, path: str, ignored_interval_count: int = 0):
        self.path = str(path)
        self.intervals: dict[str, tuple[tuple[int, int], ...]] = {
            chrom: tuple(rows) for chrom, rows in intervals.items() if rows
        }
        self._starts: dict[str, tuple[int, ...]] = {
            chrom: tuple(start for start, _end in rows)
            for chrom, rows in self.intervals.items()
        }
        self._chrom_cache: dict[str, str | None] = {}
        self.summary = BlacklistSummary(
            path=self.path,
            interval_count=sum(len(rows) for rows in self.intervals.values()),
            blacklisted_bases=sum(
                end - start for rows in self.intervals.values() for start, end in rows
            ),
            ignored_interval_count=int(ignored_interval_count),
        )

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        if end <= start:
            return False
        chrom = self._canonical_chrom(chrom)
        if chrom is None:
            return False
        rows = self.intervals.get(chrom)
        if not rows:
            return False
        starts = self._starts[chrom]
        index = bisect_left(starts, end)
        if index > 0 and rows[index - 1][1] > start:
            return True
        return index < len(rows) and rows[index][0] < end and rows[index][1] > start

    def overlap_bases(self, chrom: str, start: int, end: int) -> int:
        if end <= start:
            return 0
        chrom = self._canonical_chrom(chrom)
        if chrom is None:
            return 0
        rows = self.intervals.get(chrom)
        if not rows:
            return 0
        starts = self._starts[chrom]
        index = max(0, bisect_left(starts, start) - 1)
        total = 0
        while index < len(rows):
            left, right = rows[index]
            if left >= end:
                break
            total += max(0, min(end, right) - max(start, left))
            index += 1
        return total

    def overlapping_intervals(
        self, chrom: str, start: int, end: int
    ) -> Iterator[tuple[int, int]]:
        """Yield merged blacklist intervals intersecting a half-open region."""
        if end <= start:
            return
        chrom = self._canonical_chrom(chrom)
        if chrom is None:
            return
        rows = self.intervals.get(chrom)
        if not rows:
            return
        starts = self._starts[chrom]
        index = max(0, bisect_left(starts, start) - 1)
        while index < len(rows):
            left, right = rows[index]
            if left >= end:
                break
            if right > start:
                yield left, right
            index += 1

    def _canonical_chrom(self, chrom: str) -> str | None:
        """Resolve conservative ``chr``/mitochondrial aliases once."""
        if chrom in self._chrom_cache:
            return self._chrom_cache[chrom]
        if chrom in self.intervals:
            resolved: str | None = chrom
        else:
            try:
                resolved = resolve_contig_name(
                    chrom, list(self.intervals), source_label="blacklist"
                )
            except KeyError:
                resolved = None
        self._chrom_cache[chrom] = resolved
        return resolved

    def mask_values(
        self,
        chrom: str,
        region_start: int,
        values,
        *,
        masked_value: float = float("nan"),
    ) -> int:
        """Mask overlapping positions in a one-dimensional NumPy-like array.

        Returns the number of positions masked. The array is modified in place.
        """
        region_end = region_start + len(values)
        masked = 0
        for left, right in self.overlapping_intervals(chrom, region_start, region_end):
            local_start = max(left, region_start) - region_start
            local_end = min(right, region_end) - region_start
            if local_end > local_start:
                values[local_start:local_end] = masked_value
                masked += local_end - local_start
        return masked

    def valid_mask(self, chrom: str, start: int, end: int):
        """Return a Boolean NumPy mask that is false at blacklisted bases."""
        import numpy as np

        mask = np.ones(max(0, end - start), dtype=bool)
        for left, right in self.overlapping_intervals(chrom, start, end):
            mask[max(left, start) - start : min(right, end) - start] = False
        return mask

    def allowed_start_mask(self, chrom: str, region_start: int, region_end: int, length: int) -> list[bool]:
        """Return validity for each fragment start fitting in the region."""
        maximum = region_end - length
        if maximum < region_start:
            return []
        return [
            not self.overlaps(chrom, start, start + length)
            for start in range(region_start, maximum + 1)
        ]


def _merge(rows: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    output: list[tuple[int, int]] = []
    for start, end in sorted(rows):
        if not output or start > output[-1][1]:
            output.append((start, end))
        else:
            output[-1] = (output[-1][0], max(output[-1][1], end))
    return output


def load_blacklist(
    path: str | Path | None,
    references: Sequence[str],
    lengths: Sequence[int],
) -> BlacklistIndex | None:
    """Load, canonicalise, clip and merge a BED blacklist."""
    if path is None:
        return None
    blacklist_path = Path(path)
    if not blacklist_path.exists():
        raise FileNotFoundError(blacklist_path)
    length_by_chrom = dict(zip(references, map(int, lengths)))
    collected: dict[str, list[tuple[int, int]]] = {chrom: [] for chrom in references}
    ignored_interval_count = 0
    with open_text(blacklist_path) as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{blacklist_path}:{line_number}: expected at least three BED columns"
                )
            try:
                chrom = resolve_contig_name(
                    fields[0], references, source_label="analysis contigs"
                )
            except KeyError:
                # Whole-assembly blacklist files commonly contain contigs that are
                # absent from a chromosome-split BAM or selected-contig run. Such
                # intervals are irrelevant to this source and are ignored.
                ignored_interval_count += 1
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{blacklist_path}:{line_number}: start and end must be integers"
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"{blacklist_path}:{line_number}: require 0 <= start < end"
                )
            chrom_length = length_by_chrom[chrom]
            if start >= chrom_length:
                continue
            collected[chrom].append((start, min(end, chrom_length)))
    merged = {chrom: _merge(rows) for chrom, rows in collected.items() if rows}
    return BlacklistIndex(merged, path=str(blacklist_path.resolve()), ignored_interval_count=ignored_interval_count)


def load_blacklist_unbounded(path: str | Path | None) -> BlacklistIndex | None:
    """Load and merge a blacklist without requiring reference lengths.

    This variant is for interval-only analyses whose input does not carry a
    chromosome-size dictionary. Coordinates remain exactly as supplied.
    """
    if path is None:
        return None
    blacklist_path = Path(path)
    if not blacklist_path.exists():
        raise FileNotFoundError(blacklist_path)
    collected: dict[str, list[tuple[int, int]]] = {}
    with open_text(blacklist_path) as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{blacklist_path}:{line_number}: expected at least three BED columns"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{blacklist_path}:{line_number}: start and end must be integers"
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"{blacklist_path}:{line_number}: require 0 <= start < end"
                )
            collected.setdefault(fields[0], []).append((start, end))
    return BlacklistIndex(
        {chrom: _merge(rows) for chrom, rows in collected.items()},
        path=str(blacklist_path.resolve()),
    )


def write_canonical_blacklist(index: BlacklistIndex, output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wt", encoding="utf-8") as handle:
        for chrom, rows in index.intervals.items():
            for start, end in rows:
                handle.write(f"{chrom}\t{start}\t{end}\n")
    return output


# UCSC hg19/GRCh37 primary-sequence lengths. A selected-contig table is
# accepted when every recognised primary contig has its exact hg19 length and
# at least one such contig is present. Unknown alternate contigs are ignored.
HG19_PRIMARY_LENGTHS: Mapping[str, int] = {
    "1": 249250621, "2": 243199373, "3": 198022430,
    "4": 191154276, "5": 180915260, "6": 171115067,
    "7": 159138663, "8": 146364022, "9": 141213431,
    "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248,
    "19": 59128983, "20": 63025520, "21": 48129895,
    "22": 51304566, "X": 155270560, "Y": 59373566,
    "M": 16571,
}


def _primary_contig_token(name: str) -> str:
    token = str(name).strip()
    if token.lower().startswith("chr"):
        token = token[3:]
    token = token.upper()
    return "M" if token == "MT" else token


def is_hg19_reference(
    references: Sequence[str] | Mapping[str, int],
    lengths: Sequence[int] | None = None,
) -> bool:
    """Return true only when selected primary-contig lengths match hg19."""
    rows = (
        references.items()
        if isinstance(references, Mapping)
        else zip(references, lengths or ())
    )
    recognised = 0
    for chrom, raw_length in rows:
        token = _primary_contig_token(str(chrom))
        expected = HG19_PRIMARY_LENGTHS.get(token)
        if expected is None:
            continue
        recognised += 1
        if int(raw_length) != expected:
            return False
    return recognised > 0
