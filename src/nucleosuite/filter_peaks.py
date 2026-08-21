#!/usr/bin/env python3
"""Filter nucleosome/peak intervals by score and interval length."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence, TextIO

import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
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
    retained_records: int
    malformed_records: int
    length_filtered_records: int
    score_filtered_records: int
    percentile_threshold: float | None

    @property
    def retained_percent(self) -> float:
        if self.valid_records == 0:
            return 0.0
        return self.retained_records / self.valid_records * 100.0


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


def _validate_filters(
    *,
    score_column: int,
    min_score: float | None,
    max_score: float | None,
    score_percentile: float | None,
    min_length: int | None,
    max_length: int | None,
    score_scale: float,
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


def _parse_record(
    path: str | Path,
    line_number: int,
    text: str,
    *,
    score_column: int,
) -> tuple[list[str], int, float]:
    fields = _split_fields(text)
    score_index = score_column - 1
    if len(fields) < 3:
        raise ValueError("expected at least BED3")
    if score_index >= len(fields):
        raise ValueError(
            f"score column {score_column} is absent from a {len(fields)}-column record"
        )
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as exc:
        raise ValueError("BED start/end must be integers") from exc
    if start < 0 or end <= start:
        raise ValueError("BED start/end must satisfy 0 <= start < end")
    try:
        score = float(fields[score_index])
    except ValueError as exc:
        raise ValueError(f"score column {score_column} is not numeric") from exc
    if not math.isfinite(score):
        raise ValueError(f"score column {score_column} is not finite")
    return fields, end - start, score


def collect_length_eligible_scores(
    input_path: str | Path,
    *,
    score_column: int,
    min_length: int | None,
    max_length: int | None,
    abs_score: bool,
    strict: bool,
    progress: ProgressReporter | None = None,
) -> tuple[np.ndarray, int, int, int, int]:
    """First pass: collect scores after the region-length filter."""
    scores: list[float] = []
    total = valid = malformed = length_filtered = 0
    seen_contigs: set[str] = set()
    with open_interval_text(input_path) as source:
        for line_number, raw in enumerate(source, 1):
            text = raw.rstrip("\n")
            stripped = text.strip()
            if not stripped or stripped.startswith(HEADER_PREFIXES):
                continue
            total += 1
            try:
                fields, length, score = _parse_record(
                    input_path, line_number, text, score_column=score_column
                )
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
            scores.append(abs(score) if abs_score else score)
    return np.asarray(scores, dtype=float), total, valid, malformed, length_filtered


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
        ("absolute_score", str(bool(abs_score)).lower()),
        ("score_scale", _format_number(score_scale)),
        ("total_records", summary.total_records),
        ("valid_records", summary.valid_records),
        ("length_eligible_records", summary.length_eligible_records),
        ("retained_records", summary.retained_records),
        ("malformed_records", summary.malformed_records),
        ("length_filtered_records", summary.length_filtered_records),
        ("score_filtered_records", summary.score_filtered_records),
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
    chrom_sizes: str | Path | None = None,
    summary_output: str | Path | None = None,
    strict: bool = False,
    progress: ProgressReporter | None = None,
) -> tuple[Path, Path, PeakFilterSummary]:
    """Filter a BED/BED.gz/bigBed while preserving its interval representation."""
    _validate_filters(
        score_column=score_column,
        min_score=min_score,
        max_score=max_score,
        score_percentile=score_percentile,
        min_length=min_length,
        max_length=max_length,
        score_scale=score_scale,
    )
    source_path = Path(input_path)
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    fmt = resolved_output_format(source_path, output_format)
    if fmt == "bigbed" and score_column != 5 and (
        abs_score or not math.isclose(score_scale, 1.0, rel_tol=0.0, abs_tol=1e-15)
    ):
        raise ValueError(
            "bigBed score transformation requires --score-column 5 because BED column 5 "
            "is the standard bigBed score field"
        )

    if progress is not None:
        progress.file_start("peaks", source_path)
        progress.stage("Reading peak scores and applying region-length eligibility")

    scores, total, valid, malformed, length_filtered = collect_length_eligible_scores(
        source_path,
        score_column=score_column,
        min_length=min_length,
        max_length=max_length,
        abs_score=abs_score,
        strict=strict,
        progress=progress,
    )
    if valid == 0:
        raise ValueError(f"No valid scored BED records were found in: {source_path}")
    if scores.size == 0:
        raise ValueError("No peaks remain after the requested region-length filter")

    percentile_threshold: float | None = None
    if score_percentile is not None:
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
    )
    summary_path = Path(summary_output) if summary_output is not None else default_summary_path(output_path)

    if progress is not None:
        detail = (
            f"score percentile {score_percentile:g} -> threshold {percentile_threshold:.12g}"
            if percentile_threshold is not None
            else "absolute score bounds"
        )
        progress.stage(f"Writing filtered intervals using {detail}")

    temporary_dir: tempfile.TemporaryDirectory[str] | None = None
    if fmt == "bigbed":
        temporary_dir = tempfile.TemporaryDirectory(prefix="nucleosuite_filter_peaks_")
        text_output = Path(temporary_dir.name) / "filtered.bed"
        compressed = False
    else:
        text_output = output_path
        compressed = fmt == "bed.gz"

    retained = 0
    length_eligible = 0
    score_filtered = 0
    score_index = score_column - 1
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
                    fields, length, raw_score = _parse_record(
                        source_path, line_number, text, score_column=score_column
                    )
                except ValueError as exc:
                    if strict:
                        raise ValueError(f"{source_path}:{line_number}: {exc}") from exc
                    continue
                if min_length is not None and length < min_length:
                    continue
                if max_length is not None and length > max_length:
                    continue
                length_eligible += 1
                filter_score = abs(raw_score) if abs_score else raw_score
                keep = True
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

                output_score = transformed_score(
                    raw_score, abs_score=abs_score, score_scale=score_scale
                )
                if abs_score or not math.isclose(score_scale, 1.0, rel_tol=0.0, abs_tol=1e-15):
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
                # Keep an explicit BED marker path from the converter semantics but
                # report the requested bigBed path for a predictable CLI contract.
                output_path = Path(str(output_path) + ".empty")
    finally:
        if temporary_dir is not None:
            temporary_dir.cleanup()

    summary = PeakFilterSummary(
        total_records=total,
        valid_records=valid,
        length_eligible_records=length_eligible,
        retained_records=retained,
        malformed_records=malformed,
        length_filtered_records=length_filtered,
        score_filtered_records=score_filtered,
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
        summary=summary,
    )
    return output_path, summary_path, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite filter-peaks",
        description=(
            "Filter nucleosome/peak BED, BED.gz, or bigBed intervals by absolute "
            "score, score percentile, and region length. The filtered output uses "
            "the same interval format as the input by default."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("input", help="Input peak BED, BED.gz, or bigBed file.")
    parser.add_argument(
        "--score-column", type=int, default=5,
        help="One-based numeric score column (default: 5).",
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
            "after region-length filtering."
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
            chrom_sizes=args.chrom_sizes,
            summary_output=args.summary_output,
            strict=args.strict,
            progress=reporter,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"[ERROR] {exc}")
        return 2
    print(f"filtered_peaks\t{output}")
    print(f"summary\t{summary_path}")
    print(
        f"[INFO] Retained {result.retained_records:,}/{result.valid_records:,} valid peaks "
        f"({result.retained_percent:.2f}%); length-filtered "
        f"{result.length_filtered_records:,}, score-filtered {result.score_filtered_records:,}."
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
