#!/usr/bin/env python3
"""Filter nucleosome/peak intervals by score and interval length."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence, TextIO

import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.io.intervals import (
    convert_bed_to_bigbed,
    is_bigbed_path,
    open_interval_text,
)
from nucleosuite.progress import ProgressReporter

HEADER_PREFIXES = ("#", "track", "browser")
OUTPUT_FORMATS = ("same", "bed", "bed.gz", "bigbed")


@dataclass(frozen=True)
class PeakFilterSummary:
    total_records: int
    valid_records: int
    length_eligible_records: int
    coverage_eligible_records: int
    retained_records: int
    malformed_records: int
    length_filtered_records: int
    coverage_filtered_records: int
    score_filtered_records: int
    missing_coverage_values: int
    percentile_threshold: float | None

    @property
    def retained_percent(self) -> float:
        if self.valid_records == 0:
            return 0.0
        return self.retained_records / self.valid_records * 100.0


class BigWigCoverageReader:
    """Read point values with a fixed genomic cache for efficient peak scans."""

    def __init__(self, handle, chunk_size: int = 1_000_000) -> None:
        if chunk_size < 1:
            raise ValueError("--coverage-chunk-size must be at least 1")
        self.handle = handle
        self.chunk_size = int(chunk_size)
        self.chrom_sizes = {str(k): int(v) for k, v in handle.chroms().items()}
        if not self.chrom_sizes:
            raise ValueError("The coverage BigWig contains no chromosomes")
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
                chrom, list(self.chrom_sizes), source_label="coverage BigWig"
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
                f"Peak position {chrom}:{position} lies outside coverage BigWig "
                f"chromosome length {chrom_length}"
            )
        cache_end = self._cache_start + len(self._cache_values)
        if not (
            self._cache_chrom == resolved
            and self._cache_start <= position < cache_end
        ):
            chunk_start = (position // self.chunk_size) * self.chunk_size
            chunk_end = min(chrom_length, chunk_start + self.chunk_size)
            try:
                values = self.handle.values(resolved, chunk_start, chunk_end, numpy=True)
            except TypeError:  # older pyBigWig builds
                values = self.handle.values(resolved, chunk_start, chunk_end)
            self._cache_chrom = resolved
            self._cache_start = chunk_start
            self._cache_values = np.asarray(values, dtype=float)
        value = float(self._cache_values[position - self._cache_start])
        if not math.isfinite(value):
            return 0.0, True
        return value, False


def bed_position(fields: Sequence[str], position_column: int | None) -> int:
    """Return the zero-based genomic position used to sample coverage."""
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
            f"coverage position column {position_column} is absent from a "
            f"{len(fields)}-column record"
        )
    try:
        position = int(fields[index])
    except ValueError as exc:
        raise ValueError(
            f"coverage position column {position_column} is not an integer"
        ) from exc
    if position < 0:
        raise ValueError(
            f"coverage position column {position_column} must be zero or greater"
        )
    return position


def _split_fields(text: str) -> list[str]:
    fields = text.split("\t")
    if len(fields) == 1:
        fields = text.split()
    return fields


def _format_number(value: float) -> str:
    return f"{float(value):.12g}"


def _filter_token(value: float | int | None) -> str:
    if value is None:
        return "none"
    text = _format_number(float(value))
    if text.startswith("-"):
        text = "neg" + text[1:]
    return text.replace(".", "p").replace("+", "")


def input_interval_format(path: str | Path) -> str:
    text = str(path).lower()
    if text.endswith((".bb", ".bigbed")):
        return "bigbed"
    if text.endswith(".bed.gz") or text.endswith(".gz"):
        return "bed.gz"
    return "bed"


def _strip_interval_suffix(path: str | Path) -> tuple[str, str]:
    name = Path(path).name
    lower = name.lower()
    for suffix in (".bed.gz", ".bigbed", ".bed", ".bb", ".gz"):
        if lower.endswith(suffix):
            return name[: -len(suffix)], suffix
    return Path(name).stem, Path(name).suffix


def _safe_path_token(path: str | Path) -> str:
    name = Path(path).name
    lower = name.lower()
    for suffix in (".bigwig", ".bw", ".bedgraph", ".bed.gz", ".bed", ".gz"):
        if lower.endswith(suffix):
            name = name[: -len(suffix)]
            break
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_.-")
    return token or "track"


def resolved_output_format(input_path: str | Path, requested: str) -> str:
    if requested not in OUTPUT_FORMATS:
        raise ValueError(f"Unknown output format: {requested}")
    return input_interval_format(input_path) if requested == "same" else requested


def automatic_output_path(
    input_path: str | Path,
    *,
    output_format: str,
    min_score: float | None,
    max_score: float | None,
    score_percentile: float | None,
    min_length: int | None,
    max_length: int | None,
    abs_score: bool,
    score_scale: float,
    coverage_bigwig: str | Path | None,
    min_coverage: float | None,
    coverage_position_column: int | None,
) -> Path:
    """Build a collision-resistant automatic output filename."""
    path = Path(input_path)
    stem, _suffix = _strip_interval_suffix(path)
    tokens: list[str] = []
    if score_percentile is not None:
        tokens.append(f"scorepct{_filter_token(score_percentile)}")
    else:
        if min_score is not None:
            tokens.append(f"scoremin{_filter_token(min_score)}")
        if max_score is not None:
            tokens.append(f"scoremax{_filter_token(max_score)}")
    if min_length is not None:
        tokens.append(f"lenmin{min_length}")
    if max_length is not None:
        tokens.append(f"lenmax{max_length}")
    if coverage_bigwig is not None:
        tokens.append(f"cov{_safe_path_token(coverage_bigwig)}")
        tokens.append(f"covmin{_filter_token(min_coverage)}")
        if coverage_position_column is not None:
            tokens.append(f"covposcol{coverage_position_column}")
    if abs_score:
        tokens.append("absscore")
    if not math.isclose(score_scale, 1.0, rel_tol=0.0, abs_tol=1e-15):
        tokens.append(f"scorescale{_filter_token(score_scale)}")
    detail = "_".join(tokens) if tokens else "all"
    extension = {"bed": ".bed", "bed.gz": ".bed.gz", "bigbed": ".bb"}[output_format]
    return path.with_name(f"{stem}_filtered_{detail}{extension}")


def default_summary_path(output_path: str | Path) -> Path:
    path = Path(output_path)
    stem, _suffix = _strip_interval_suffix(path)
    return path.with_name(f"{stem}_summary.tsv")


def transformed_score(raw_score: float, *, abs_score: bool, score_scale: float) -> float:
    value = abs(raw_score) if abs_score else raw_score
    return value * score_scale


def _score_required(
    *,
    min_score: float | None,
    max_score: float | None,
    score_percentile: float | None,
    abs_score: bool,
    score_scale: float,
) -> bool:
    return (
        min_score is not None
        or max_score is not None
        or score_percentile is not None
        or abs_score
        or not math.isclose(score_scale, 1.0, rel_tol=0.0, abs_tol=1e-15)
    )


def _validate_filters(
    *,
    score_column: int,
    min_score: float | None,
    max_score: float | None,
    score_percentile: float | None,
    min_length: int | None,
    max_length: int | None,
    score_scale: float,
    coverage_bigwig: str | Path | None,
    min_coverage: float | None,
    coverage_position_column: int | None,
    coverage_chunk_size: int,
) -> None:
    if score_column < 1:
        raise ValueError("--score-column must be a one-based column number of 1 or greater")
    for option, value in (("--min-score", min_score), ("--max-score", max_score)):
        if value is not None and not math.isfinite(value):
            raise ValueError(f"{option} must be finite")
    if min_score is not None and max_score is not None and max_score < min_score:
        raise ValueError("--max-score must be greater than or equal to --min-score")
    if score_percentile is not None:
        if min_score is not None or max_score is not None:
            raise ValueError("--score-percentile cannot be combined with --min-score or --max-score")
        if not math.isfinite(score_percentile) or not 0.0 <= score_percentile <= 100.0:
            raise ValueError("--score-percentile must be between 0 and 100")
    for option, value in (("--min-length", min_length), ("--max-length", max_length)):
        if value is not None and value < 1:
            raise ValueError(f"{option} must be at least 1 bp")
    if min_length is not None and max_length is not None and max_length < min_length:
        raise ValueError("--max-length must be greater than or equal to --min-length")
    if not math.isfinite(score_scale) or score_scale <= 0:
        raise ValueError("--score-scale must be a finite value greater than 0")
    if (coverage_bigwig is None) != (min_coverage is None):
        raise ValueError("--coverage-bigwig and --min-coverage must be supplied together")
    if min_coverage is not None and (
        not math.isfinite(min_coverage) or min_coverage < 0
    ):
        raise ValueError("--min-coverage must be a finite value of 0 or greater")
    if coverage_position_column is not None:
        if coverage_bigwig is None:
            raise ValueError("--coverage-position-column requires --coverage-bigwig")
        if coverage_position_column < 1:
            raise ValueError("--coverage-position-column must be a one-based column number")
    if coverage_chunk_size < 1:
        raise ValueError("--coverage-chunk-size must be at least 1")


def _parse_interval_record(text: str) -> tuple[list[str], int]:
    fields = _split_fields(text)
    if len(fields) < 3:
        raise ValueError("expected at least BED3")
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as exc:
        raise ValueError("BED start/end must be integers") from exc
    if start < 0 or end <= start:
        raise ValueError("BED start/end must satisfy 0 <= start < end")
    return fields, end - start


def _parse_score(fields: Sequence[str], score_column: int) -> float:
    score_index = score_column - 1
    if score_index >= len(fields):
        raise ValueError(
            f"score column {score_column} is absent from a {len(fields)}-column record"
        )
    try:
        score = float(fields[score_index])
    except ValueError as exc:
        raise ValueError(f"score column {score_column} is not numeric") from exc
    if not math.isfinite(score):
        raise ValueError(f"score column {score_column} is not finite")
    return score


def _open_coverage_reader(path: str | Path, chunk_size: int):
    try:
        import pyBigWig  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("Coverage filtering requires pyBigWig") from exc
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise OSError(f"Could not open coverage BigWig: {path}")
    try:
        reader = BigWigCoverageReader(handle, chunk_size=chunk_size)
    except Exception:
        handle.close()
        raise
    return handle, reader


def collect_filter_eligible_scores(
    input_path: str | Path,
    *,
    score_column: int,
    min_length: int | None,
    max_length: int | None,
    abs_score: bool,
    require_score: bool,
    coverage_bigwig: str | Path | None,
    min_coverage: float | None,
    coverage_position_column: int | None,
    coverage_chunk_size: int,
    strict: bool,
    progress: ProgressReporter | None = None,
) -> tuple[np.ndarray, int, int, int, int, int, int, int, int]:
    """Collect score values after length and optional coverage eligibility."""
    scores: list[float] = []
    total = valid = malformed = length_filtered = coverage_filtered = 0
    length_eligible = coverage_eligible = missing_coverage = 0
    seen_contigs: set[str] = set()
    coverage_handle = coverage_reader = None
    if coverage_bigwig is not None:
        coverage_handle, coverage_reader = _open_coverage_reader(
            coverage_bigwig, coverage_chunk_size
        )
    try:
        with open_interval_text(input_path) as source:
            for line_number, raw in enumerate(source, 1):
                text = raw.rstrip("\n")
                stripped = text.strip()
                if not stripped or stripped.startswith(HEADER_PREFIXES):
                    continue
                total += 1
                try:
                    fields, length = _parse_interval_record(text)
                    raw_score = _parse_score(fields, score_column) if require_score else None
                except ValueError as exc:
                    malformed += 1
                    if strict:
                        raise ValueError(f"{input_path}:{line_number}: {exc}") from exc
                    continue
                valid += 1
                chrom = fields[0]
                if progress is not None and chrom not in seen_contigs:
                    seen_contigs.add(chrom)
                    progress.reading_contig("peaks", chrom)
                if min_length is not None and length < min_length:
                    length_filtered += 1
                    continue
                if max_length is not None and length > max_length:
                    length_filtered += 1
                    continue
                length_eligible += 1
                if coverage_reader is not None:
                    try:
                        position = bed_position(fields, coverage_position_column)
                        coverage, was_missing = coverage_reader.value(chrom, position)
                    except ValueError as exc:
                        raise ValueError(f"{input_path}:{line_number}: {exc}") from exc
                    missing_coverage += int(was_missing)
                    if coverage < float(min_coverage):
                        coverage_filtered += 1
                        continue
                coverage_eligible += 1
                if raw_score is not None:
                    scores.append(abs(raw_score) if abs_score else raw_score)
    finally:
        if coverage_handle is not None:
            coverage_handle.close()
    return (
        np.asarray(scores, dtype=float),
        total,
        valid,
        malformed,
        length_filtered,
        length_eligible,
        coverage_filtered,
        coverage_eligible,
        missing_coverage,
    )


def infer_bigbed_chrom_sizes(path: str | Path) -> dict[str, int]:
    """Read chromosome sizes embedded in a bigBed through pyBigWig."""
    try:
        import pyBigWig  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "Automatic chromosome-size discovery for bigBed input requires pyBigWig"
        ) from exc
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise OSError(f"Could not open bigBed to obtain chromosome sizes: {path}")
    try:
        chroms = handle.chroms()
    finally:
        handle.close()
    if not chroms:
        raise ValueError(f"No chromosome sizes were found in bigBed: {path}")
    return {str(chrom): int(length) for chrom, length in chroms.items()}


def _open_output_text(path: Path, *, compressed: bool) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compressed:
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("wt", encoding="utf-8", newline="")


def write_summary(
    output_path: str | Path,
    *,
    input_path: str | Path,
    filtered_path: str | Path,
    score_column: int,
    min_score: float | None,
    max_score: float | None,
    score_percentile: float | None,
    percentile_threshold: float | None,
    min_length: int | None,
    max_length: int | None,
    abs_score: bool,
    score_scale: float,
    coverage_bigwig: str | Path | None,
    min_coverage: float | None,
    coverage_position_column: int | None,
    coverage_chunk_size: int,
    summary: PeakFilterSummary,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        ("input", input_path),
        ("output", filtered_path),
        ("score_column", score_column),
        ("min_score", "" if min_score is None else _format_number(min_score)),
        ("max_score", "" if max_score is None else _format_number(max_score)),
        ("score_percentile", "" if score_percentile is None else _format_number(score_percentile)),
        ("percentile_score_threshold", "" if percentile_threshold is None else _format_number(percentile_threshold)),
        ("min_length", "" if min_length is None else min_length),
        ("max_length", "" if max_length is None else max_length),
        ("coverage_bigwig", "" if coverage_bigwig is None else coverage_bigwig),
        ("min_coverage", "" if min_coverage is None else _format_number(min_coverage)),
        (
            "coverage_position_source",
            "" if coverage_bigwig is None else (
                "interval_midpoint" if coverage_position_column is None
                else f"bed_column_{coverage_position_column}"
            ),
        ),
        ("coverage_chunk_size", "" if coverage_bigwig is None else coverage_chunk_size),
        ("absolute_score", str(bool(abs_score)).lower()),
        ("score_scale", _format_number(score_scale)),
        ("total_records", summary.total_records),
        ("valid_records", summary.valid_records),
        ("length_eligible_records", summary.length_eligible_records),
        ("coverage_eligible_records", summary.coverage_eligible_records),
        ("retained_records", summary.retained_records),
        ("malformed_records", summary.malformed_records),
        ("length_filtered_records", summary.length_filtered_records),
        ("coverage_filtered_records", summary.coverage_filtered_records),
        ("score_filtered_records", summary.score_filtered_records),
        ("missing_coverage_values_treated_as_zero", summary.missing_coverage_values),
        ("retained_percent_of_valid", f"{summary.retained_percent:.6f}"),
    ]
    with output.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)
    return output


def filter_peaks(
    input_path: str | Path,
    *,
    output: str | Path | None = None,
    output_format: str = "same",
    score_column: int = 5,
    min_score: float | None = None,
    max_score: float | None = None,
    score_percentile: float | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    abs_score: bool = False,
    score_scale: float = 1.0,
    coverage_bigwig: str | Path | None = None,
    min_coverage: float | None = None,
    coverage_position_column: int | None = None,
    coverage_chunk_size: int = 1_000_000,
    chrom_sizes: str | Path | None = None,
    summary_output: str | Path | None = None,
    strict: bool = False,
    progress: ProgressReporter | None = None,
) -> tuple[Path, Path, PeakFilterSummary]:
    """Filter BED/BED.gz/bigBed peaks by score, length, and/or coverage."""
    _validate_filters(
        score_column=score_column,
        min_score=min_score,
        max_score=max_score,
        score_percentile=score_percentile,
        min_length=min_length,
        max_length=max_length,
        score_scale=score_scale,
        coverage_bigwig=coverage_bigwig,
        min_coverage=min_coverage,
        coverage_position_column=coverage_position_column,
        coverage_chunk_size=coverage_chunk_size,
    )
    source_path = Path(input_path)
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    if coverage_bigwig is not None and not Path(coverage_bigwig).exists():
        raise FileNotFoundError(coverage_bigwig)
    fmt = resolved_output_format(source_path, output_format)
    require_score = _score_required(
        min_score=min_score,
        max_score=max_score,
        score_percentile=score_percentile,
        abs_score=abs_score,
        score_scale=score_scale,
    )
    if fmt == "bigbed" and score_column != 5 and (
        abs_score or not math.isclose(score_scale, 1.0, rel_tol=0.0, abs_tol=1e-15)
    ):
        raise ValueError(
            "bigBed score transformation requires --score-column 5 because BED column 5 "
            "is the standard bigBed score field"
        )

    if progress is not None:
        progress.file_start("peaks", source_path)
        progress.stage("Reading peaks and applying length/coverage eligibility")

    (
        scores,
        total,
        valid,
        malformed,
        length_filtered,
        length_eligible,
        coverage_filtered,
        coverage_eligible,
        missing_coverage,
    ) = collect_filter_eligible_scores(
        source_path,
        score_column=score_column,
        min_length=min_length,
        max_length=max_length,
        abs_score=abs_score,
        require_score=require_score,
        coverage_bigwig=coverage_bigwig,
        min_coverage=min_coverage,
        coverage_position_column=coverage_position_column,
        coverage_chunk_size=coverage_chunk_size,
        strict=strict,
        progress=progress,
    )
    if valid == 0:
        raise ValueError(f"No valid BED records were found in: {source_path}")
    percentile_threshold: float | None = None
    if score_percentile is not None:
        if coverage_eligible == 0 or scores.size == 0:
            raise ValueError(
                "No scored peaks remain after the requested length/coverage filters "
                "for percentile calculation"
            )
        percentile_threshold = float(np.percentile(scores, score_percentile))

    output_path = Path(output) if output is not None else automatic_output_path(
        source_path,
        output_format=fmt,
        min_score=min_score,
        max_score=max_score,
        score_percentile=score_percentile,
        min_length=min_length,
        max_length=max_length,
        abs_score=abs_score,
        score_scale=score_scale,
        coverage_bigwig=coverage_bigwig,
        min_coverage=min_coverage,
        coverage_position_column=coverage_position_column,
    )
    summary_path = Path(summary_output) if summary_output is not None else default_summary_path(output_path)

    if progress is not None:
        filters: list[str] = []
        if min_length is not None or max_length is not None:
            filters.append("length")
        if coverage_bigwig is not None:
            filters.append("coverage")
        if percentile_threshold is not None:
            filters.append(
                f"score percentile {score_percentile:g} -> {percentile_threshold:.12g}"
            )
        elif min_score is not None or max_score is not None:
            filters.append("absolute score")
        progress.stage(
            "Writing filtered intervals" + (f" ({', '.join(filters)})" if filters else "")
        )

    temporary_dir: tempfile.TemporaryDirectory[str] | None = None
    if fmt == "bigbed":
        temporary_dir = tempfile.TemporaryDirectory(prefix="nucleosuite_filter_peaks_")
        text_output = Path(temporary_dir.name) / "filtered.bed"
        compressed = False
    else:
        text_output = output_path
        compressed = fmt == "bed.gz"

    retained = 0
    score_filtered = 0
    score_index = score_column - 1
    coverage_handle = coverage_reader = None
    if coverage_bigwig is not None:
        coverage_handle, coverage_reader = _open_coverage_reader(
            coverage_bigwig, coverage_chunk_size
        )
    try:
        with open_interval_text(source_path) as source, _open_output_text(text_output, compressed=compressed) as destination:
            for line_number, raw in enumerate(source, 1):
                text = raw.rstrip("\n")
                stripped = text.strip()
                if not stripped or stripped.startswith(HEADER_PREFIXES):
                    if fmt != "bigbed":
                        destination.write(text + "\n")
                    continue
                try:
                    fields, length = _parse_interval_record(text)
                    raw_score = _parse_score(fields, score_column) if require_score else None
                except ValueError as exc:
                    if strict:
                        raise ValueError(f"{source_path}:{line_number}: {exc}") from exc
                    continue
                if min_length is not None and length < min_length:
                    continue
                if max_length is not None and length > max_length:
                    continue
                if coverage_reader is not None:
                    try:
                        position = bed_position(fields, coverage_position_column)
                        coverage, _was_missing = coverage_reader.value(fields[0], position)
                    except ValueError as exc:
                        raise ValueError(f"{source_path}:{line_number}: {exc}") from exc
                    if coverage < float(min_coverage):
                        continue

                keep = True
                if raw_score is not None:
                    filter_score = abs(raw_score) if abs_score else raw_score
                    if percentile_threshold is not None:
                        keep = filter_score >= percentile_threshold
                    else:
                        if min_score is not None and filter_score < min_score:
                            keep = False
                        if max_score is not None and filter_score > max_score:
                            keep = False
                if not keep:
                    score_filtered += 1
                    continue

                if raw_score is not None and (
                    abs_score
                    or not math.isclose(score_scale, 1.0, rel_tol=0.0, abs_tol=1e-15)
                ):
                    output_score = transformed_score(
                        raw_score, abs_score=abs_score, score_scale=score_scale
                    )
                    fields[score_index] = _format_number(output_score)
                    destination.write("\t".join(fields) + "\n")
                else:
                    destination.write(text + "\n")
                retained += 1

        if fmt == "bigbed":
            sizes = chrom_sizes
            if sizes is None and is_bigbed_path(source_path):
                sizes = infer_bigbed_chrom_sizes(source_path)
            if sizes is None:
                raise ValueError(
                    "--chrom-sizes is required for bigBed output when chromosome sizes "
                    "cannot be inherited from a bigBed input"
                )
            converted = convert_bed_to_bigbed(
                text_output,
                sizes,
                output_path,
                bigbed_score_multiplier=1.0,
            )
            if converted is None:
                output_path = Path(str(output_path) + ".empty")
    finally:
        if coverage_handle is not None:
            coverage_handle.close()
        if temporary_dir is not None:
            temporary_dir.cleanup()

    summary = PeakFilterSummary(
        total_records=total,
        valid_records=valid,
        length_eligible_records=length_eligible,
        coverage_eligible_records=coverage_eligible,
        retained_records=retained,
        malformed_records=malformed,
        length_filtered_records=length_filtered,
        coverage_filtered_records=coverage_filtered,
        score_filtered_records=score_filtered,
        missing_coverage_values=missing_coverage,
        percentile_threshold=percentile_threshold,
    )
    write_summary(
        summary_path,
        input_path=source_path,
        filtered_path=output_path,
        score_column=score_column,
        min_score=min_score,
        max_score=max_score,
        score_percentile=score_percentile,
        percentile_threshold=percentile_threshold,
        min_length=min_length,
        max_length=max_length,
        abs_score=abs_score,
        score_scale=score_scale,
        coverage_bigwig=coverage_bigwig,
        min_coverage=min_coverage,
        coverage_position_column=coverage_position_column,
        coverage_chunk_size=coverage_chunk_size,
        summary=summary,
    )
    return output_path, summary_path, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite filter-peaks",
        description=(
            "Filter nucleosome/peak BED, BED.gz, or bigBed intervals by score, "
            "score percentile, region length, and/or BigWig coverage. The filtered "
            "output uses the same interval format as the input by default."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("input", help="Input peak BED, BED.gz, or bigBed file.")
    parser.add_argument(
        "--score-column", type=int, default=5,
        help=(
            "One-based numeric score column (default: 5). Required only when a "
            "score filter or score transformation is requested."
        ),
    )
    score_group = parser.add_argument_group("score filtering")
    score_group.add_argument(
        "--min-score", type=float,
        help="Minimum score to retain (inclusive). Cannot be combined with --score-percentile.",
    )
    score_group.add_argument(
        "--max-score", type=float,
        help="Maximum score to retain (inclusive). Cannot be combined with --score-percentile.",
    )
    score_group.add_argument(
        "--score-percentile", type=float,
        help=(
            "Retain peaks at or above this score percentile (0-100), calculated "
            "after region-length and coverage filtering."
        ),
    )
    score_group.add_argument(
        "--abs-score", action="store_true",
        help=(
            "Use absolute score values for filtering and write absolute scores to "
            "the output. Off by default."
        ),
    )
    score_group.add_argument(
        "--score-scale", type=float, default=1.0,
        help=(
            "Multiply retained output scores by this factor (default: 1). For "
            "bigBed output, BED column-5 scores are rounded and clamped to 0-1000."
        ),
    )
    length_group = parser.add_argument_group("region-length filtering")
    length_group.add_argument(
        "--min-length", type=int,
        help="Minimum BED interval length in bp to retain (inclusive).",
    )
    length_group.add_argument(
        "--max-length", type=int,
        help="Maximum BED interval length in bp to retain (inclusive).",
    )
    coverage_group = parser.add_argument_group("coverage filtering")
    coverage_group.add_argument(
        "--coverage-bigwig",
        help="BigWig track sampled to filter peaks by coverage.",
    )
    coverage_group.add_argument(
        "--min-coverage", type=float,
        help=(
            "Minimum coverage required at the selected peak position (inclusive). "
            "Requires --coverage-bigwig."
        ),
    )
    coverage_group.add_argument(
        "--coverage-position-column", type=int,
        help=(
            "One-based BED column containing the zero-based genomic position to "
            "sample from --coverage-bigwig. Default: interval midpoint."
        ),
    )
    coverage_group.add_argument(
        "--coverage-chunk-size", type=int, default=1_000_000,
        help="BigWig cache chunk size in bp (default: 1000000).",
    )
    parser.add_argument(
        "--output",
        help="Filtered interval output path. Default: an automatic parameterized name beside the input.",
    )
    parser.add_argument(
        "--output-format", choices=OUTPUT_FORMATS, default="same",
        help="Output interval format (default: same as input).",
    )
    parser.add_argument(
        "--chrom-sizes",
        help=(
            "Chromosome sizes for bigBed output. A bigBed input supplies its own "
            "sizes automatically when possible."
        ),
    )
    parser.add_argument(
        "--summary-output",
        help="Summary TSV path. Default: <filtered-output-stem>_summary.tsv.",
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="Stop at the first malformed data record instead of skipping it.",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    reporter = ProgressReporter("filter-peaks")
    try:
        output, summary_path, result = filter_peaks(
            args.input,
            output=args.output,
            output_format=args.output_format,
            score_column=args.score_column,
            min_score=args.min_score,
            max_score=args.max_score,
            score_percentile=args.score_percentile,
            min_length=args.min_length,
            max_length=args.max_length,
            abs_score=args.abs_score,
            score_scale=args.score_scale,
            coverage_bigwig=args.coverage_bigwig,
            min_coverage=args.min_coverage,
            coverage_position_column=args.coverage_position_column,
            coverage_chunk_size=args.coverage_chunk_size,
            chrom_sizes=args.chrom_sizes,
            summary_output=args.summary_output,
            strict=args.strict,
            progress=reporter,
        )
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
        print(f"[ERROR] {exc}")
        return 2
    print(f"filtered_peaks\t{output}")
    print(f"summary\t{summary_path}")
    print(
        f"[INFO] Retained {result.retained_records:,}/{result.valid_records:,} valid peaks "
        f"({result.retained_percent:.2f}%); length-filtered "
        f"{result.length_filtered_records:,}, coverage-filtered "
        f"{result.coverage_filtered_records:,}, score-filtered "
        f"{result.score_filtered_records:,}."
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
