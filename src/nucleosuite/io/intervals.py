"""Standard BED and bigBed interval helpers.

NucleoSuite writes ordinary BED text while analyses are running, then optionally
converts the completed file to UCSC bigBed with ``bedToBigBed``.  BigBed inputs
are materialised temporarily with ``bigBedToBed`` so commands that already read
BED text can accept either representation without duplicating parsers.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TextIO
from nucleosuite.core.chrom_sizes import read_chrom_sizes_source
from nucleosuite.core.regions import resolve_contig_name


INTERVAL_FORMATS = ("bed", "bigbed", "both")
BIGBED_SUFFIXES = (".bb", ".bigbed")


def is_bigbed_path(path: str | Path) -> bool:
    return str(path).lower().endswith(BIGBED_SUFFIXES)


def _require_executable(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        package = "ucsc-bedtobigbed" if name == "bedToBigBed" else "ucsc-bigbedtobed"
        raise RuntimeError(
            f"{name} is required for bigBed support. Install it with "
            f"'mamba install -c bioconda {package}'."
        )
    return executable


class ManagedTextHandle:
    """Text handle that also owns any temporary bigBed conversion directory."""

    def __init__(self, path: str | Path):
        self.input_path = Path(path)
        self._temp_dir: tempfile.TemporaryDirectory | None = None
        if is_bigbed_path(self.input_path):
            executable = _require_executable("bigBedToBed")
            self._temp_dir = tempfile.TemporaryDirectory(prefix="nucleosuite_bigbed_")
            bed_path = Path(self._temp_dir.name) / (self.input_path.stem + ".bed")
            completed = subprocess.run(
                [executable, str(self.input_path), str(bed_path)],
                capture_output=True,
                text=True,
                check=False,
            )
            if completed.returncode != 0:
                self._temp_dir.cleanup()
                detail = (completed.stderr or completed.stdout).strip()
                raise OSError(f"bigBedToBed failed for {self.input_path}: {detail}")
            self.handle = bed_path.open("rt", encoding="utf-8")
        elif self.input_path.suffix.lower() == ".gz":
            self.handle = gzip.open(self.input_path, "rt", encoding="utf-8")
        else:
            self.handle = self.input_path.open("rt", encoding="utf-8")

    def __enter__(self) -> TextIO:
        return self.handle

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def __iter__(self):
        return iter(self.handle)

    def __getattr__(self, name):
        return getattr(self.handle, name)

    def close(self) -> None:
        try:
            self.handle.close()
        finally:
            if self._temp_dir is not None:
                self._temp_dir.cleanup()
                self._temp_dir = None


def open_interval_text(path: str | Path) -> ManagedTextHandle:
    """Open BED text, BED.gz, or bigBed for iteration and context management."""
    return ManagedTextHandle(path)


def normalise_chrom_sizes(
    chrom_sizes: Mapping[str, int] | Sequence[tuple[str, int]] | str | Path,
) -> list[tuple[str, int]]:
    """Return ordered ``(chromosome, length)`` records."""
    if isinstance(chrom_sizes, Mapping):
        rows = [(str(name), int(length)) for name, length in chrom_sizes.items()]
    elif isinstance(chrom_sizes, (str, Path)):
        rows = read_chrom_sizes_source(chrom_sizes)
    else:
        rows = [(str(name), int(length)) for name, length in chrom_sizes]

    if not rows:
        raise ValueError("Chromosome sizes are required for bigBed output")
    seen: set[str] = set()
    cleaned: list[tuple[str, int]] = []
    for name, length in rows:
        if name in seen:
            raise ValueError(f"Duplicate chromosome in sizes: {name}")
        if length < 1:
            raise ValueError(f"Chromosome length must be positive: {name}={length}")
        seen.add(name)
        cleaned.append((name, length))
    return cleaned


def write_chrom_sizes(
    chrom_sizes: Mapping[str, int] | Sequence[tuple[str, int]] | str | Path,
    output_path: str | Path,
) -> Path:
    rows = normalise_chrom_sizes(chrom_sizes)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wt", encoding="utf-8") as handle:
        for chrom, length in rows:
            handle.write(f"{chrom}\t{length}\n")
    return output


def _read_bed_records(
    path: str | Path, *, bigbed_score_multiplier: float = 1.0
) -> tuple[list[list[str]], int]:
    records: list[list[str]] = []
    column_count: int | None = None
    with open_interval_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            text = raw.rstrip("\n")
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < 3 or len(fields) > 12:
                raise ValueError(
                    f"{path}:{line_no}: bigBed conversion supports standard BED3-BED12; "
                    f"found {len(fields)} columns"
                )
            if column_count is None:
                column_count = len(fields)
            elif len(fields) != column_count:
                raise ValueError(
                    f"{path}:{line_no}: inconsistent BED column count "
                    f"({len(fields)} versus {column_count})"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid BED coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: require 0 <= start < end")
            if len(fields) >= 5:
                try:
                    score = int(round(float(fields[4]) * bigbed_score_multiplier))
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_no}: invalid BED score") from exc
                fields[4] = str(max(0, min(1000, score)))
            records.append(fields)
    return records, int(column_count or 3)


def convert_bed_to_bigbed(
    bed_path: str | Path,
    chrom_sizes: Mapping[str, int] | Sequence[tuple[str, int]] | str | Path,
    output_path: str | Path | None = None,
    *,
    remove_bed: bool = False,
    bigbed_score_multiplier: float = 1.0,
) -> Path | None:
    """Sort and convert a standard BED3-BED12 file to bigBed.

    Empty BED files are retained and return ``None`` because UCSC bigBed cannot
    represent an interval collection with no records.
    """
    bed = Path(bed_path)
    if output_path is None:
        output = bed.with_suffix(".bb")
    else:
        output = Path(output_path)
    rows = normalise_chrom_sizes(chrom_sizes)
    rank = {name: index for index, (name, _length) in enumerate(rows)}
    lengths = dict(rows)
    records, column_count = _read_bed_records(
        bed, bigbed_score_multiplier=bigbed_score_multiplier
    )
    empty_marker = Path(str(output) + ".empty")
    if not records:
        output.unlink(missing_ok=True)
        empty_marker.parent.mkdir(parents=True, exist_ok=True)
        empty_marker.write_text(
            "status\tempty\nrecord_count\t0\n", encoding="utf-8"
        )
        print(
            f"[WARNING] {bed} contains no intervals; wrote empty-output marker "
            f"{empty_marker}."
        )
        return None
    empty_marker.unlink(missing_ok=True)

    for fields in records:
        source_chrom = fields[0]
        try:
            chrom = resolve_contig_name(
                source_chrom, list(rank), source_label="chromosome sizes"
            )
        except KeyError as exc:
            raise ValueError(
                f"BED chromosome {source_chrom!r} is absent from the chromosome sizes"
            ) from exc
        fields[0] = chrom
        if int(fields[2]) > lengths[chrom]:
            raise ValueError(
                f"BED interval {chrom}:{fields[1]}-{fields[2]} exceeds chromosome "
                f"length {lengths[chrom]}"
            )
    records.sort(key=lambda fields: (rank[fields[0]], int(fields[1]), int(fields[2])))

    executable = _require_executable("bedToBigBed")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="nucleosuite_bed2bb_") as temp_dir:
        temp = Path(temp_dir)
        sorted_bed = temp / "sorted.bed"
        sizes_path = temp / "chrom.sizes"
        with sorted_bed.open("wt", encoding="utf-8") as handle:
            for fields in records:
                handle.write("\t".join(fields) + "\n")
        write_chrom_sizes(rows, sizes_path)
        completed = subprocess.run(
            [
                executable,
                f"-type=bed{column_count}",
                str(sorted_bed),
                str(sizes_path),
                str(output),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout).strip()
            raise OSError(f"bedToBigBed failed for {bed}: {detail}")
    if remove_bed:
        bed.unlink(missing_ok=True)
    return output


def finalise_interval_files(
    bed_paths: Iterable[str | Path],
    interval_format: str,
    chrom_sizes: Mapping[str, int] | Sequence[tuple[str, int]] | str | Path,
    *,
    bigbed_score_multiplier: float = 1.0,
) -> list[Path]:
    """Apply ``bed``, ``bigbed`` or ``both`` output policy to completed BEDs."""
    if interval_format not in INTERVAL_FORMATS:
        raise ValueError(f"Unknown interval format: {interval_format}")
    paths = [Path(path) for path in bed_paths]
    if interval_format == "bed":
        return paths
    outputs: list[Path] = []
    for bed_path in paths:
        if not bed_path.exists():
            continue
        bigbed = convert_bed_to_bigbed(
            bed_path,
            chrom_sizes,
            remove_bed=interval_format == "bigbed",
            bigbed_score_multiplier=bigbed_score_multiplier,
        )
        if interval_format == "both":
            outputs.append(bed_path)
        if bigbed is not None:
            outputs.append(bigbed)
        else:
            marker = Path(str(bed_path.with_suffix(".bb")) + ".empty")
            if marker.exists():
                outputs.append(marker)
            if interval_format == "bigbed":
                outputs.append(bed_path)
    return outputs
