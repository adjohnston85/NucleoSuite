#!/usr/bin/env python3
"""Filter BED peaks by BigWig coverage at a representative peak position."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.io.intervals import open_interval_text
from nucleosuite.progress import ProgressReporter


@dataclass(frozen=True)
class CoverageFilterSummary:
    total_peaks: int
    retained_peaks: int
    filtered_peaks: int
    missing_values: int

    @property
    def retained_percent(self) -> float:
        if self.total_peaks == 0:
            return 0.0
        return self.retained_peaks / self.total_peaks * 100.0


def _threshold_token(value: float) -> str:
    """Return a compact deterministic token for an output filename."""
    return format(float(value), ".12g")


def default_output_path(input_bed: str | Path, threshold: float) -> Path:
    """Derive ``<input>_coverage_ge<threshold>.bed`` from a BED-like input."""
    path = Path(input_bed)
    name = path.name
    lower = name.lower()
    for suffix in (".bed.gz", ".bigbed", ".bed", ".bb", ".gz"):
        if lower.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return path.with_name(f"{name}_coverage_ge{_threshold_token(threshold)}.bed")


def default_summary_path(output_bed: str | Path) -> Path:
    path = Path(output_bed)
    return path.with_name(f"{path.stem}_summary.tsv")


def bed_position(fields: Sequence[str], position_column: int | None) -> int:
    """Return the zero-based point used to query coverage for one BED record."""
    try:
        start = int(fields[1])
        end = int(fields[2])
    except (IndexError, ValueError) as exc:
        raise ValueError("invalid BED coordinates") from exc
    if start < 0 or end <= start:
        raise ValueError("require 0 <= start < end")

    if position_column is None:
        return (start + end) // 2

    index = position_column - 1
    if index >= len(fields):
        raise ValueError(
            f"position column {position_column} is absent from a {len(fields)}-column record"
        )
    try:
        position = int(fields[index])
    except ValueError as exc:
        raise ValueError(f"position column {position_column} is not an integer") from exc
    if position < 0:
        raise ValueError(f"position column {position_column} must be zero or greater")
    return position


class BigWigCoverageReader:
    """Read point values with a small fixed genomic cache for efficient BED scans."""

    def __init__(self, handle, chunk_size: int = 1_000_000) -> None:
        if chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")
        self.handle = handle
        self.chunk_size = int(chunk_size)
        self.chrom_sizes = {str(k): int(v) for k, v in handle.chroms().items()}
        if not self.chrom_sizes:
            raise ValueError("The input BigWig contains no chromosomes")
        self._resolved: dict[str, str] = {}
        self._cache_chrom: str | None = None
        self._cache_start = 0
        self._cache_values = np.empty(0, dtype=float)

    def _resolve(self, chrom: str) -> str:
        resolved = self._resolved.get(chrom)
        if resolved is not None:
            return resolved
        try:
            resolved = resolve_contig_name(
                chrom, list(self.chrom_sizes), source_label="BigWig"
            )
        except KeyError as exc:
            raise ValueError(str(exc)) from exc
        self._resolved[chrom] = resolved
        return resolved

    def value(self, chrom: str, position: int) -> tuple[float, bool]:
        resolved = self._resolve(chrom)
        chrom_length = self.chrom_sizes[resolved]
        if position < 0 or position >= chrom_length:
            raise ValueError(
                f"Peak position {chrom}:{position} lies outside BigWig chromosome "
                f"length {chrom_length}"
            )

        cache_end = self._cache_start + len(self._cache_values)
        if not (
            self._cache_chrom == resolved
            and self._cache_start <= position < cache_end
        ):
            chunk_start = (position // self.chunk_size) * self.chunk_size
            chunk_end = min(chrom_length, chunk_start + self.chunk_size)
            try:
                values = self.handle.values(
                    resolved, chunk_start, chunk_end, numpy=True
                )
            except TypeError:  # older pyBigWig builds
                values = self.handle.values(resolved, chunk_start, chunk_end)
            self._cache_chrom = resolved
            self._cache_start = chunk_start
            self._cache_values = np.asarray(values, dtype=float)

        value = float(self._cache_values[position - self._cache_start])
        if not math.isfinite(value):
            return 0.0, True
        return value, False


def write_summary(
    output_path: str | Path,
    *,
    input_bed: str | Path,
    bigwig: str | Path,
    filtered_bed: str | Path,
    threshold: float,
    position_column: int | None,
    summary: CoverageFilterSummary,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["input_bed", str(input_bed)])
        writer.writerow(["coverage_bigwig", str(bigwig)])
        writer.writerow(["output_bed", str(filtered_bed)])
        writer.writerow(["coverage_threshold", _threshold_token(threshold)])
        writer.writerow(
            [
                "position_source",
                "interval_midpoint" if position_column is None else f"bed_column_{position_column}",
            ]
        )
        writer.writerow(["total_peaks", summary.total_peaks])
        writer.writerow(["retained_peaks", summary.retained_peaks])
        writer.writerow(["filtered_peaks", summary.filtered_peaks])
        writer.writerow(["retained_percent", f"{summary.retained_percent:.6f}"])
        writer.writerow(["missing_bigwig_values_treated_as_zero", summary.missing_values])
    return output


def filter_bed_by_coverage(
    input_bed: str | Path,
    bigwig_path: str | Path,
    *,
    coverage_threshold: float,
    position_column: int | None = None,
    output_bed: str | Path | None = None,
    summary_output: str | Path | None = None,
    chunk_size: int = 1_000_000,
    progress: ProgressReporter | None = None,
) -> tuple[Path, Path, CoverageFilterSummary]:
    if not math.isfinite(coverage_threshold) or coverage_threshold < 0:
        raise ValueError("--coverage-threshold must be a finite value of 0 or greater")
    if position_column is not None and position_column < 1:
        raise ValueError("--position-column must be a one-based column number of 1 or greater")
    if chunk_size < 1:
        raise ValueError("--chunk-size must be at least 1")

    input_path = Path(input_bed)
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    bigwig = Path(bigwig_path)
    if not bigwig.exists():
        raise FileNotFoundError(bigwig)

    output = Path(output_bed) if output_bed is not None else default_output_path(
        input_path, coverage_threshold
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    summary_path = (
        Path(summary_output)
        if summary_output is not None
        else default_summary_path(output)
    )

    try:
        import pyBigWig  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "filter-coverage requires pyBigWig. Install it with "
            "'conda install -c bioconda pybigwig' or 'pip install pyBigWig'."
        ) from exc

    if progress is not None:
        progress.file_start("peaks", input_path)
        progress.stage(f"Opening coverage track: {bigwig}")

    total = retained = missing = 0
    seen_contigs: set[str] = set()
    handle = pyBigWig.open(str(bigwig))
    if handle is None:
        raise OSError(f"Could not open BigWig: {bigwig}")
    try:
        reader = BigWigCoverageReader(handle, chunk_size=chunk_size)
        with open_interval_text(input_path) as source, output.open(
            "wt", encoding="utf-8"
        ) as destination:
            for line_number, raw in enumerate(source, 1):
                text = raw.rstrip("\n")
                stripped = text.strip()
                if not stripped or stripped.startswith(("#", "track", "browser")):
                    destination.write(text + "\n")
                    continue
                fields = text.split("\t") if "\t" in text else text.split()
                if len(fields) < 3:
                    raise ValueError(
                        f"{input_path}:{line_number}: expected at least BED3"
                    )
                chrom = fields[0]
                if progress is not None and chrom not in seen_contigs:
                    seen_contigs.add(chrom)
                    progress.reading_contig("peaks", chrom)
                try:
                    position = bed_position(fields, position_column)
                except ValueError as exc:
                    raise ValueError(f"{input_path}:{line_number}: {exc}") from exc
                try:
                    coverage, was_missing = reader.value(chrom, position)
                except ValueError as exc:
                    raise ValueError(f"{input_path}:{line_number}: {exc}") from exc
                total += 1
                missing += int(was_missing)
                if coverage >= coverage_threshold:
                    destination.write(text + "\n")
                    retained += 1
    finally:
        handle.close()

    result = CoverageFilterSummary(
        total_peaks=total,
        retained_peaks=retained,
        filtered_peaks=total - retained,
        missing_values=missing,
    )
    write_summary(
        summary_path,
        input_bed=input_path,
        bigwig=bigwig,
        filtered_bed=output,
        threshold=coverage_threshold,
        position_column=position_column,
        summary=result,
    )
    return output, summary_path, result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite filter-coverage",
        description=(
            "Filter BED peaks by the BigWig coverage value at each peak position. "
            "The BED interval midpoint is used by default; --position-column can "
            "supply an explicit zero-based genomic summit position from the BED."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("bed", help="Input BED, BED.gz, or bigBed peak file.")
    parser.add_argument("--bigwig", required=True, help="Coverage BigWig track.")
    parser.add_argument(
        "--coverage-threshold",
        required=True,
        type=float,
        help=(
            "Minimum coverage required at the selected peak position. Peaks are "
            "retained when coverage is greater than or equal to this value."
        ),
    )
    parser.add_argument(
        "--position-column",
        type=int,
        help=(
            "One-based BED column containing the zero-based genomic summit position. "
            "Default: use the midpoint between BED start and end."
        ),
    )
    parser.add_argument(
        "--output",
        help=(
            "Filtered BED output. Default: <input>_coverage_ge<threshold>.bed beside "
            "the input BED."
        ),
    )
    parser.add_argument(
        "--summary-output",
        help="Summary TSV path. Default: <filtered BED stem>_summary.tsv.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="BigWig cache chunk size in bp (default: 1000000).",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    reporter = ProgressReporter("filter-coverage")
    output, summary_path, result = filter_bed_by_coverage(
        args.bed,
        args.bigwig,
        coverage_threshold=args.coverage_threshold,
        position_column=args.position_column,
        output_bed=args.output,
        summary_output=args.summary_output,
        chunk_size=args.chunk_size,
        progress=reporter,
    )
    print(f"filtered_bed\t{output}")
    print(f"summary\t{summary_path}")
    print(
        f"[INFO] Coverage filter >= {_threshold_token(args.coverage_threshold)}: "
        f"retained {result.retained_peaks:,}/{result.total_peaks:,} peaks "
        f"({result.retained_percent:.2f}%); filtered {result.filtered_peaks:,}."
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
