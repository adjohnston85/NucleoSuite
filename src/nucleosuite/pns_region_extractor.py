#!/usr/bin/env python3
"""
Extract PNS-related signal and flanking peaks around BED regions.

Signal inputs
-------------
BigWig tracks are queried directly with pyBigWig. The script can extract the
standard PNS coverage and PNS score tracks, as well as additional named signal
tracks.

Peak inputs
-----------
Peak tracks may be bigBed, BED, or BED.GZ files. For each input BED region, the
nearest upstream and downstream peak centres within ``--peak-flank-bp`` are
reported. The default peak centre is BED column 7 when present and numeric;
otherwise the interval midpoint is used. The default score is BED column 5.

Examples
--------
Extract the standard four PNS outputs::

    nucleosuite region-extract \
        --bed regions.bed \
        --coverage-bw sample_coverage.bw \
        --score-bw sample_pns.bw \
        --nucleosome-peaks sample_nucleosome_regions.bb \
        --breakpoint-peaks sample_breakpoint_peaks.bb \
        --peak-flank-bp 2000 \
        --out-prefix sample_regions

Extract only PNS signal::

    nucleosuite region-extract \
        --bed regions.bed \
        --score-bw sample_pns.bw

Add arbitrary named tracks::

    nucleosuite region-extract \
        --bed regions.bed \
        --signal-track dyad=sample_dyad.bw \
        --peak-track protected=sample_protected_regions.bb

Outputs
-------
Depending on the supplied tracks, the script writes:

* ``<prefix>_flanking_peaks.tsv``
* ``<prefix>_<signal-name>_signal.tsv``
* ``<prefix>_skipped_lines.tsv`` when invalid BED records were skipped

Required Python package
-----------------------
pyBigWig
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, TextIO

from nucleosuite.io import open_text as open_interval_text
from nucleosuite.parallel import add_parallel_arguments
from nucleosuite.partitioned import run_partitioned_command
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
import numpy as np

try:
    import pyBigWig
except ImportError:  # pragma: no cover - depends on runtime environment
    pyBigWig = None


def require_pybigwig() -> None:
    if pyBigWig is None:
        raise ExtractionError(
            "pyBigWig is required for bigWig and bigBed inputs. Install it with "
            "'conda install -c bioconda pybigwig' or 'python -m pip install pyBigWig'."
        )


STANDARD_BED_COLUMNS = ("chrom", "start", "end", "name", "score", "strand")
TRACK_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


@dataclass(frozen=True)
class BedRegion:
    """One valid input BED interval, retaining every original column."""

    index: int
    line_no: int
    fields: tuple[str, ...]

    @property
    def chrom(self) -> str:
        return self.fields[0]

    @property
    def start(self) -> int:
        return int(self.fields[1])

    @property
    def end(self) -> int:
        return int(self.fields[2])

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2

    @property
    def region_id(self) -> str:
        if len(self.fields) >= 4 and self.fields[3] not in {"", "."}:
            return self.fields[3]
        return f"region_{self.index}"


@dataclass(frozen=True)
class SkippedBedLine:
    line_no: int
    reason: str
    raw_line: str


@dataclass(frozen=True)
class PeakRecord:
    chrom: str
    start: int
    end: int
    center: int
    score: str
    fields: tuple[str, ...]


@dataclass(frozen=True)
class NamedPath:
    name: str
    path: str


@dataclass(frozen=True)
class ExtractionOutputs:
    flanking_peaks: Optional[str]
    signal_matrices: dict[str, str]
    skipped_lines: Optional[str]


class ExtractionError(RuntimeError):
    """Raised when an input track cannot be queried safely."""


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def open_text(path: str) -> TextIO:
    """Open BED text, BED.gz, bigBed, or standard input."""
    if path == "-":
        return sys.stdin
    return open_interval_text(path)


def parse_bed(
    bed_path: str,
    *,
    strict: bool = False,
    blacklist: BlacklistIndex | None = None,
) -> tuple[list[BedRegion], list[SkippedBedLine], int]:
    """Read a BED file and retain all columns from valid records."""
    regions: list[BedRegion] = []
    skipped: list[SkippedBedLine] = []
    max_columns = 3

    handle = open_text(bed_path)
    should_close = handle is not sys.stdin
    try:
        for line_no, raw in enumerate(handle, start=1):
            raw_line = raw.rstrip("\r\n")
            stripped = raw_line.strip()
            if not stripped:
                continue
            lowered = stripped.lower()
            if stripped.startswith("#") or lowered.startswith("track") or lowered.startswith("browser"):
                continue

            fields = tuple(stripped.split())
            reason: Optional[str] = None
            if len(fields) < 3:
                reason = "fewer than 3 columns"
            else:
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    reason = "non-integer start or end"
                else:
                    if start < 0:
                        reason = "start is negative"
                    elif end <= start:
                        reason = "end is not greater than start"

            if reason is not None:
                if strict:
                    raise ValueError(f"Invalid BED record at line {line_no}: {reason}: {raw_line}")
                skipped.append(SkippedBedLine(line_no, reason, raw_line))
                continue

            if blacklist is not None and blacklist.overlaps(fields[0], start, end):
                skipped.append(SkippedBedLine(line_no, "overlaps blacklist", raw_line))
                continue

            region = BedRegion(
                index=len(regions) + 1,
                line_no=line_no,
                fields=fields,
            )
            regions.append(region)
            max_columns = max(max_columns, len(fields))
    finally:
        if should_close:
            handle.close()

    if not regions:
        raise ValueError("No valid BED regions were found")

    return regions, skipped, max_columns


def bed_headers(max_columns: int) -> list[str]:
    headers: list[str] = []
    for index in range(1, max_columns + 1):
        if index <= len(STANDARD_BED_COLUMNS):
            headers.append(STANDARD_BED_COLUMNS[index - 1])
        else:
            headers.append(f"bed_col_{index}")
    return headers


def padded_bed_fields(region: BedRegion, max_columns: int) -> list[str]:
    return list(region.fields) + [""] * (max_columns - len(region.fields))


def parse_named_path(value: str, option_name: str) -> NamedPath:
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"{option_name} must use NAME=PATH")
    name, path = value.split("=", 1)
    name = name.strip()
    path = path.strip()
    if not name or not path:
        raise argparse.ArgumentTypeError(f"{option_name} must use non-empty NAME=PATH")
    if not TRACK_NAME_RE.fullmatch(name):
        raise argparse.ArgumentTypeError(
            f"Invalid track name {name!r}; use letters, numbers, '.', '_' or '-', beginning with a letter"
        )
    return NamedPath(name=name, path=path)


def output_safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._-") or "track"


def unique_named_paths(items: Iterable[NamedPath], kind: str) -> list[NamedPath]:
    seen: set[str] = set()
    result: list[NamedPath] = []
    for item in items:
        key = item.name.lower()
        if key in seen:
            raise ValueError(f"Duplicate {kind} track name: {item.name}")
        seen.add(key)
        result.append(item)
    return result


def validate_regular_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} does not exist or is not a regular file: {path}")


def auto_output_prefix(bed_path: str) -> str:
    if bed_path == "-":
        return "pns_region_extractor"
    name = Path(bed_path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".bed", ".tsv", ".txt"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return str(Path(bed_path).with_name(name + "_pns_region_extractor"))


def ensure_output_paths_available(paths: Sequence[str], overwrite: bool) -> None:
    for path in paths:
        parent = os.path.dirname(os.path.abspath(path))
        os.makedirs(parent, exist_ok=True)
        if os.path.exists(path) and not overwrite:
            raise FileExistsError(f"Output already exists: {path}. Use --overwrite to replace it.")


def resolve_chrom_name(
    requested: str,
    available: Iterable[str],
    *,
    mode: str,
) -> Optional[str]:
    """Resolve exact names, with optional conservative chr-prefix aliases."""
    available_set = set(available)
    if requested in available_set:
        return requested
    if mode == "strict":
        return None

    candidates: list[str] = []
    if requested.startswith("chr"):
        candidates.append(requested[3:])
    else:
        candidates.append("chr" + requested)

    mitochondrial_aliases = {
        "M": ("MT", "chrM", "chrMT"),
        "MT": ("M", "chrM", "chrMT"),
        "chrM": ("M", "MT", "chrMT"),
        "chrMT": ("M", "MT", "chrM"),
    }
    candidates.extend(mitochondrial_aliases.get(requested, ()))

    matches = [candidate for candidate in candidates if candidate in available_set]
    if len(matches) == 1:
        return matches[0]
    return None


def parse_peak_fields(
    fields: Sequence[str],
    *,
    center_column: int,
    score_column: int,
) -> Optional[PeakRecord]:
    if len(fields) < 3:
        return None
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError:
        return None
    if start < 0 or end <= start:
        return None

    center = (start + end) // 2
    if center_column > 0 and len(fields) >= center_column:
        try:
            center = int(float(fields[center_column - 1]))
        except ValueError:
            pass

    score = "."
    if score_column > 0 and len(fields) >= score_column and fields[score_column - 1] != "":
        score = fields[score_column - 1]

    return PeakRecord(
        chrom=fields[0],
        start=start,
        end=end,
        center=center,
        score=score,
        fields=tuple(fields),
    )


class SignalTrack:
    """An open bigWig track."""

    def __init__(self, named_path: NamedPath, chrom_mode: str):
        self.name = named_path.name
        self.path = named_path.path
        self.chrom_mode = chrom_mode
        validate_regular_file(self.path, f"Signal track {self.name!r}")
        require_pybigwig()
        try:
            self.handle = pyBigWig.open(self.path)
        except Exception as exc:
            raise ExtractionError(f"Could not open signal track {self.path}: {exc}") from exc
        if self.handle is None or not self.handle.isBigWig():
            if self.handle is not None:
                self.handle.close()
            raise ExtractionError(f"Signal track is not a bigWig file: {self.path}")
        self.chroms: dict[str, int] = dict(self.handle.chroms())

    def close(self) -> None:
        self.handle.close()

    def values(
        self,
        region: BedRegion,
        *,
        missing_chrom: str,
    ) -> tuple[Optional[str], list[float]]:
        resolved = resolve_chrom_name(region.chrom, self.chroms, mode=self.chrom_mode)
        width = region.end - region.start
        if resolved is None:
            if missing_chrom == "error":
                raise ExtractionError(
                    f"Chromosome {region.chrom!r} from BED line {region.line_no} is absent from signal track "
                    f"{self.name!r} ({self.path})"
                )
            return None, [math.nan] * width

        chrom_size = self.chroms[resolved]
        query_start = min(max(region.start, 0), chrom_size)
        query_end = min(max(region.end, 0), chrom_size)

        left_pad = max(0, query_start - region.start)
        right_pad = max(0, region.end - query_end)
        values: list[float] = [math.nan] * left_pad
        if query_end > query_start:
            try:
                queried = self.handle.values(resolved, query_start, query_end, numpy=False)
            except RuntimeError as exc:
                raise ExtractionError(
                    f"Failed to query signal track {self.name!r} at "
                    f"{resolved}:{query_start}-{query_end}: {exc}"
                ) from exc
            values.extend(float(value) if value is not None else math.nan for value in queried)
        values.extend([math.nan] * right_pad)

        if len(values) < width:
            values.extend([math.nan] * (width - len(values)))
        elif len(values) > width:
            values = values[:width]
        return resolved, values


class PeakTrack:
    """Interface shared by bigBed and text BED peak tracks."""

    name: str
    path: str

    def query(
        self,
        region: BedRegion,
        flank_bp: int,
        *,
        missing_chrom: str,
    ) -> tuple[Optional[str], list[PeakRecord]]:
        raise NotImplementedError

    def close(self) -> None:
        return None


class BigBedPeakTrack(PeakTrack):
    def __init__(
        self,
        named_path: NamedPath,
        *,
        chrom_mode: str,
        center_column: int,
        score_column: int,
    ):
        self.name = named_path.name
        self.path = named_path.path
        self.chrom_mode = chrom_mode
        self.center_column = center_column
        self.score_column = score_column
        validate_regular_file(self.path, f"Peak track {self.name!r}")
        require_pybigwig()
        try:
            self.handle = pyBigWig.open(self.path)
        except Exception as exc:
            raise ExtractionError(f"Could not open peak track {self.path}: {exc}") from exc
        if self.handle is None or not self.handle.isBigBed():
            if self.handle is not None:
                self.handle.close()
            raise ExtractionError(f"Peak track is not a bigBed file: {self.path}")
        self.chroms: dict[str, int] = dict(self.handle.chroms())

    def close(self) -> None:
        self.handle.close()

    def query(
        self,
        region: BedRegion,
        flank_bp: int,
        *,
        missing_chrom: str,
    ) -> tuple[Optional[str], list[PeakRecord]]:
        resolved = resolve_chrom_name(region.chrom, self.chroms, mode=self.chrom_mode)
        if resolved is None:
            if missing_chrom == "error":
                raise ExtractionError(
                    f"Chromosome {region.chrom!r} from BED line {region.line_no} is absent from peak track "
                    f"{self.name!r} ({self.path})"
                )
            return None, []

        chrom_size = self.chroms[resolved]
        query_start = max(0, region.center - flank_bp)
        query_end = min(chrom_size, region.center + flank_bp + 1)
        if query_end <= query_start:
            return resolved, []

        try:
            entries = self.handle.entries(resolved, query_start, query_end) or []
        except RuntimeError as exc:
            raise ExtractionError(
                f"Failed to query peak track {self.name!r} at "
                f"{resolved}:{query_start}-{query_end}: {exc}"
            ) from exc

        records: list[PeakRecord] = []
        for start, end, rest in entries:
            extra_fields: list[str]
            if rest is None or rest == "":
                extra_fields = []
            elif isinstance(rest, bytes):
                extra_fields = rest.decode("utf-8").split("\t")
            else:
                extra_fields = str(rest).split("\t")
            fields = [resolved, str(start), str(end), *extra_fields]
            record = parse_peak_fields(
                fields,
                center_column=self.center_column,
                score_column=self.score_column,
            )
            if record is not None:
                records.append(record)

        records.sort(key=lambda record: (record.center, record.start, record.end))
        return resolved, records


class BedPeakTrack(PeakTrack):
    def __init__(
        self,
        named_path: NamedPath,
        *,
        chrom_mode: str,
        center_column: int,
        score_column: int,
        strict: bool,
    ):
        self.name = named_path.name
        self.path = named_path.path
        self.chrom_mode = chrom_mode
        validate_regular_file(self.path, f"Peak track {self.name!r}")

        by_chrom: dict[str, list[PeakRecord]] = {}
        skipped = 0
        handle = open_text(self.path)
        try:
            for line_no, raw in enumerate(handle, start=1):
                stripped = raw.strip()
                if not stripped:
                    continue
                lowered = stripped.lower()
                if stripped.startswith("#") or lowered.startswith("track") or lowered.startswith("browser"):
                    continue
                fields = stripped.split()
                record = parse_peak_fields(
                    fields,
                    center_column=center_column,
                    score_column=score_column,
                )
                if record is None:
                    if strict:
                        raise ValueError(
                            f"Invalid peak BED record in {self.path} at line {line_no}: {stripped}"
                        )
                    skipped += 1
                    continue
                by_chrom.setdefault(record.chrom, []).append(record)
        finally:
            handle.close()

        self.records: dict[str, list[PeakRecord]] = {}
        self.starts: dict[str, list[int]] = {}
        self.prefix_max_ends: dict[str, list[int]] = {}
        for chrom, records in by_chrom.items():
            records.sort(key=lambda record: (record.start, record.end, record.center))
            self.records[chrom] = records
            self.starts[chrom] = [record.start for record in records]
            prefix: list[int] = []
            max_end = -1
            for record in records:
                max_end = max(max_end, record.end)
                prefix.append(max_end)
            self.prefix_max_ends[chrom] = prefix

        self.chrom_names = set(self.records)
        if skipped:
            log(f"[{self.name}] skipped {skipped:,} invalid peak BED record(s)")

    def query(
        self,
        region: BedRegion,
        flank_bp: int,
        *,
        missing_chrom: str,
    ) -> tuple[Optional[str], list[PeakRecord]]:
        resolved = resolve_chrom_name(region.chrom, self.chrom_names, mode=self.chrom_mode)
        if resolved is None:
            if missing_chrom == "error":
                raise ExtractionError(
                    f"Chromosome {region.chrom!r} from BED line {region.line_no} is absent from peak track "
                    f"{self.name!r} ({self.path})"
                )
            return None, []

        query_start = max(0, region.center - flank_bp)
        query_end = region.center + flank_bp + 1
        records = self.records[resolved]
        starts = self.starts[resolved]
        prefix_max_ends = self.prefix_max_ends[resolved]

        first = bisect.bisect_right(prefix_max_ends, query_start)
        last = bisect.bisect_left(starts, query_end)
        matches = [
            record
            for record in records[first:last]
            if record.end > query_start and record.start < query_end
        ]
        matches.sort(key=lambda record: (record.center, record.start, record.end))
        return resolved, matches


def open_peak_track(
    named_path: NamedPath,
    *,
    chrom_mode: str,
    center_column: int,
    score_column: int,
    strict_peak_bed: bool,
) -> PeakTrack:
    lower = named_path.path.lower()
    if lower.endswith((".bb", ".bigbed")):
        return BigBedPeakTrack(
            named_path,
            chrom_mode=chrom_mode,
            center_column=center_column,
            score_column=score_column,
        )
    return BedPeakTrack(
        named_path,
        chrom_mode=chrom_mode,
        center_column=center_column,
        score_column=score_column,
        strict=strict_peak_bed,
    )


def nearest_flanking_peaks(
    peaks: Sequence[PeakRecord],
    center: int,
) -> tuple[Optional[PeakRecord], Optional[PeakRecord]]:
    """Return nearest centre strictly upstream and nearest centre at/downstream."""
    if not peaks:
        return None, None
    centers = [peak.center for peak in peaks]
    insertion = bisect.bisect_left(centers, center)
    upstream = peaks[insertion - 1] if insertion > 0 else None
    downstream = peaks[insertion] if insertion < len(peaks) else None
    return upstream, downstream


def peak_output_headers(track_name: str) -> list[str]:
    headers = [f"{track_name}_query_chrom", f"{track_name}_peaks_found_in_query"]
    for direction in ("upstream", "downstream"):
        prefix = f"{track_name}_{direction}"
        headers.extend(
            [
                f"{prefix}_chrom",
                f"{prefix}_start",
                f"{prefix}_end",
                f"{prefix}_center",
                f"{prefix}_score",
                f"{prefix}_distance_from_region_center",
            ]
        )
    return headers


def peak_output_values(
    resolved_chrom: Optional[str],
    peaks: Sequence[PeakRecord],
    region_center: int,
    missing_value: str,
) -> list[str | int]:
    values: list[str | int] = [resolved_chrom or missing_value, len(peaks)]
    upstream, downstream = nearest_flanking_peaks(peaks, region_center)
    for peak in (upstream, downstream):
        if peak is None:
            values.extend([missing_value] * 6)
        else:
            values.extend(
                [
                    peak.chrom,
                    peak.start,
                    peak.end,
                    peak.center,
                    peak.score,
                    peak.center - region_center,
                ]
            )
    return values


def format_signal_value(value: float, missing_value: str, precision: int) -> str:
    if not math.isfinite(value):
        return missing_value
    return format(value, f".{precision}g")


def write_skipped_lines(path: str, skipped: Sequence[SkippedBedLine]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["line_number", "reason", "raw_line"])
        for item in skipped:
            writer.writerow([item.line_no, item.reason, item.raw_line])


def write_flanking_peak_table(
    path: str,
    regions: Sequence[BedRegion],
    max_bed_columns: int,
    peak_tracks: Sequence[PeakTrack],
    *,
    flank_bp: int,
    missing_value: str,
    missing_chrom: str,
    progress_every: int,
    blacklist: BlacklistIndex | None = None,
) -> None:
    headers = ["region_index", "region_id", *bed_headers(max_bed_columns), "region_center", "peak_query_flank_bp"]
    for track in peak_tracks:
        headers.extend(peak_output_headers(track.name))

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        total = len(regions)
        for completed, region in enumerate(regions, start=1):
            row: list[str | int] = [
                region.index,
                region.region_id,
                *padded_bed_fields(region, max_bed_columns),
                region.center,
                flank_bp,
            ]
            for track in peak_tracks:
                resolved, peaks = track.query(
                    region,
                    flank_bp,
                    missing_chrom=missing_chrom,
                )
                if blacklist is not None:
                    peaks = [
                        peak for peak in peaks
                        if not blacklist.overlaps(peak.chrom, peak.start, peak.end)
                    ]
                row.extend(peak_output_values(resolved, peaks, region.center, missing_value))
            writer.writerow(row)
            if progress_every > 0 and (completed % progress_every == 0 or completed == total):
                log(f"[peaks] processed {completed:,}/{total:,} regions")


def write_signal_matrix(
    path: str,
    regions: Sequence[BedRegion],
    max_bed_columns: int,
    signal_track: SignalTrack,
    *,
    missing_value: str,
    missing_chrom: str,
    precision: int,
    progress_every: int,
    blacklist: BlacklistIndex | None = None,
) -> None:
    max_width = max(region.end - region.start for region in regions)
    headers = [
        "region_index",
        "region_id",
        *bed_headers(max_bed_columns),
        "region_center",
        f"{signal_track.name}_query_chrom",
        *[f"offset_{offset}" for offset in range(max_width)],
    ]

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        total = len(regions)
        for completed, region in enumerate(regions, start=1):
            resolved, values = signal_track.values(region, missing_chrom=missing_chrom)
            if blacklist is not None:
                array = np.asarray(values, dtype=float)
                blacklist.mask_values(region.chrom, region.start, array)
                values = array.tolist()
            formatted = [
                format_signal_value(value, missing_value, precision)
                for value in values
            ]
            if len(formatted) < max_width:
                formatted.extend([missing_value] * (max_width - len(formatted)))
            row: list[str | int] = [
                region.index,
                region.region_id,
                *padded_bed_fields(region, max_bed_columns),
                region.center,
                resolved or missing_value,
                *formatted,
            ]
            writer.writerow(row)
            if progress_every > 0 and (completed % progress_every == 0 or completed == total):
                log(f"[{signal_track.name}] processed {completed:,}/{total:,} regions")


def run_extraction(
    *,
    bed_path: str,
    signal_specs: Sequence[NamedPath],
    peak_specs: Sequence[NamedPath],
    out_prefix: str,
    peak_flank_bp: int = 2000,
    peak_center_column: int = 7,
    peak_score_column: int = 5,
    missing_value: str = "nan",
    chrom_mode: str = "auto",
    missing_chrom: str = "error",
    strict_bed: bool = False,
    strict_peak_bed: bool = False,
    precision: int = 8,
    progress_every: int = 1000,
    overwrite: bool = False,
    blacklist_bed: str | None = None,
) -> ExtractionOutputs:
    """Run extraction and return the paths that were written."""
    if peak_flank_bp < 0:
        raise ValueError("peak_flank_bp must be at least 0")
    if peak_center_column < 0 or peak_score_column < 0:
        raise ValueError("Peak centre and score columns must be at least 0")
    if precision < 1:
        raise ValueError("precision must be at least 1")
    if progress_every < 0:
        raise ValueError("progress_every must be at least 0")
    if not out_prefix:
        raise ValueError("out_prefix must not be empty")
    if not signal_specs and not peak_specs:
        raise ValueError("At least one signal or peak track must be supplied")

    if bed_path != "-":
        validate_regular_file(bed_path, "Input BED")

    signal_specs = unique_named_paths(signal_specs, "signal")
    peak_specs = unique_named_paths(peak_specs, "peak")
    blacklist = load_blacklist_unbounded(blacklist_bed)
    regions, skipped, max_bed_columns = parse_bed(
        bed_path, strict=strict_bed, blacklist=blacklist
    )
    log(f"Loaded {len(regions):,} valid BED region(s)")
    if skipped:
        log(f"Skipped {len(skipped):,} invalid BED line(s)")

    flanking_path = f"{out_prefix}_flanking_peaks.tsv" if peak_specs else None
    signal_paths = {
        spec.name: f"{out_prefix}_{output_safe_name(spec.name)}_signal.tsv"
        for spec in signal_specs
    }
    skipped_path = f"{out_prefix}_skipped_lines.tsv" if skipped else None
    all_outputs = [path for path in [flanking_path, skipped_path, *signal_paths.values()] if path]
    ensure_output_paths_available(all_outputs, overwrite)

    signal_tracks: list[SignalTrack] = []
    peak_tracks: list[PeakTrack] = []
    try:
        for spec in signal_specs:
            signal_tracks.append(SignalTrack(spec, chrom_mode=chrom_mode))
        for spec in peak_specs:
            peak_tracks.append(
                open_peak_track(
                    spec,
                    chrom_mode=chrom_mode,
                    center_column=peak_center_column,
                    score_column=peak_score_column,
                    strict_peak_bed=strict_peak_bed,
                )
            )

        if flanking_path is not None:
            log(f"Writing {flanking_path}")
            write_flanking_peak_table(
                flanking_path,
                regions,
                max_bed_columns,
                peak_tracks,
                flank_bp=peak_flank_bp,
                missing_value=missing_value,
                missing_chrom=missing_chrom,
                progress_every=progress_every,
                blacklist=blacklist,
            )

        for track in signal_tracks:
            path = signal_paths[track.name]
            log(f"Writing {path}")
            write_signal_matrix(
                path,
                regions,
                max_bed_columns,
                track,
                missing_value=missing_value,
                missing_chrom=missing_chrom,
                precision=precision,
                progress_every=progress_every,
                blacklist=blacklist,
            )

        if skipped_path is not None:
            write_skipped_lines(skipped_path, skipped)
            log(f"Wrote {skipped_path}")
    except Exception:
        # Remove partial outputs so a failed run is not mistaken for completion.
        for path in all_outputs:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except OSError:
                pass
        raise
    finally:
        for track in signal_tracks:
            track.close()
        for track in peak_tracks:
            track.close()

    return ExtractionOutputs(
        flanking_peaks=flanking_path,
        signal_matrices=signal_paths,
        skipped_lines=skipped_path,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite region-extract",
        description=(
            "Extract per-base bigWig signal and nearest upstream/downstream peaks "
            "for regions in a BED file."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--bed", required=True, help="Input BED/BED.GZ file, or '-' for stdin")
    parser.add_argument(
        "--blacklist-bed",
        help=(
            "BED blacklist. Overlapping anchors and peaks are excluded; "
            "overlapping signal positions are written as missing."
        ),
    )

    signal_group = parser.add_argument_group("signal tracks")
    signal_group.add_argument("--coverage-bw", help="Coverage bigWig")
    signal_group.add_argument("--score-bw", "--pns-bw", dest="pns_bw", help="Nucleosome-score bigWig; --pns-bw is retained as a compatibility alias.")
    signal_group.add_argument(
        "--signal-track",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Additional named bigWig; may be supplied multiple times",
    )

    peak_group = parser.add_argument_group("peak tracks")
    peak_group.add_argument("--nucleosome-peaks", help="Nucleosome peak bigBed, BED, or BED.GZ")
    peak_group.add_argument("--breakpoint-peaks", help="Breakpoint peak bigBed, BED, or BED.GZ")
    peak_group.add_argument(
        "--peak-track",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Additional named bigBed/BED/BED.GZ; may be supplied multiple times",
    )
    peak_group.add_argument(
        "--peak-flank-bp",
        type=int,
        default=2000,
        help="Distance queried on each side of the input region centre (default: 2000)",
    )
    peak_group.add_argument(
        "--peak-center-column",
        type=int,
        default=7,
        help="1-based peak BED column containing an absolute centre; 0 always uses interval midpoint",
    )
    peak_group.add_argument(
        "--peak-score-column",
        type=int,
        default=5,
        help="1-based peak BED score column; 0 reports '.'",
    )
    peak_group.add_argument(
        "--strict-peak-bed",
        action="store_true",
        help="Stop at the first malformed record in text BED peak files",
    )

    output_group = parser.add_argument_group("output")
    output_group.add_argument(
        "--out-prefix",
        help="Output path prefix; derived from --bed when omitted",
    )
    output_group.add_argument("--missing-value", default="nan", help="Text used for unavailable values")
    output_group.add_argument(
        "--precision",
        type=int,
        default=8,
        help="Significant digits written for signal values",
    )
    output_group.add_argument("--overwrite", action="store_true", help="Replace existing output files")

    validation_group = parser.add_argument_group("validation and progress")
    validation_group.add_argument(
        "--chrom-mode",
        choices=("auto", "strict"),
        default="auto",
        help="Allow conservative chr-prefix/mitochondrial chromosome aliases, or require exact names",
    )
    validation_group.add_argument(
        "--missing-chrom",
        choices=("error", "fill"),
        default="error",
        help="Stop when a BED chromosome is absent from a track, or fill that track with missing values",
    )
    validation_group.add_argument(
        "--strict-bed",
        action="store_true",
        help="Stop at the first malformed input BED record rather than reporting skipped lines",
    )
    validation_group.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="Print progress after this many regions; 0 disables progress messages",
    )
    add_parallel_arguments(parser)
    return parser


def collect_track_specs(args: argparse.Namespace) -> tuple[list[NamedPath], list[NamedPath]]:
    signal_specs: list[NamedPath] = []
    peak_specs: list[NamedPath] = []

    if args.coverage_bw:
        signal_specs.append(NamedPath("coverage", args.coverage_bw))
    if args.pns_bw:
        signal_specs.append(NamedPath("pns", args.pns_bw))
    for value in args.signal_track:
        signal_specs.append(parse_named_path(value, "--signal-track"))

    if args.nucleosome_peaks:
        peak_specs.append(NamedPath("nucleosome", args.nucleosome_peaks))
    if args.breakpoint_peaks:
        peak_specs.append(NamedPath("breakpoint", args.breakpoint_peaks))
    for value in args.peak_track:
        peak_specs.append(parse_named_path(value, "--peak-track"))

    return signal_specs, peak_specs


def _run_serial(args: argparse.Namespace) -> int:
    signal_specs, peak_specs = collect_track_specs(args)
    out_prefix = args.out_prefix or auto_output_prefix(args.bed)
    outputs = run_extraction(
        bed_path=args.bed,
        signal_specs=signal_specs,
        peak_specs=peak_specs,
        out_prefix=out_prefix,
        peak_flank_bp=args.peak_flank_bp,
        peak_center_column=args.peak_center_column,
        peak_score_column=args.peak_score_column,
        missing_value=args.missing_value,
        chrom_mode=args.chrom_mode,
        missing_chrom=args.missing_chrom,
        strict_bed=args.strict_bed,
        strict_peak_bed=args.strict_peak_bed,
        precision=args.precision,
        progress_every=args.progress_every,
        overwrite=args.overwrite,
        blacklist_bed=args.blacklist_bed,
    )
    log("Completed successfully")
    if outputs.flanking_peaks:
        print(outputs.flanking_peaks)
    for path in outputs.signal_matrices.values():
        print(path)
    if outputs.skipped_lines:
        print(outputs.skipped_lines)
    return 0


def run(args: argparse.Namespace) -> int:
    if not args.out_prefix:
        args.out_prefix = auto_output_prefix(args.bed)
    return run_partitioned_command(
        "region-extract",
        args,
        _run_serial,
        runner_module="nucleosuite.pns_region_extractor",
        runner_function="_run_serial",
        primary_attr="bed",
        output_prefix_attr="out_prefix",
        path_attrs=("blacklist_bed",),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (
        argparse.ArgumentTypeError,
        ExtractionError,
        FileExistsError,
        FileNotFoundError,
        OSError,
        RuntimeError,
        ValueError,
    ) as exc:
        log(f"ERROR: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
