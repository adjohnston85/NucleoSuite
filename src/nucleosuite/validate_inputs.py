"""Preflight validation for suite inputs and reference resources."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.chrom_sizes import read_chrom_sizes_source
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.io import open_text
from nucleosuite.progress import ProgressReporter


@dataclass(frozen=True)
class ValidationRow:
    kind: str
    path: str
    status: str
    records: int | str
    detail: str


def _validate_interval_file(
    path: str | Path,
    *,
    kind: str,
    max_records: int | None,
    require_sorted: bool,
) -> ValidationRow:
    source = Path(path)
    if not source.is_file():
        return ValidationRow(kind, str(source), "FAIL", 0, "file not found")
    records = 0
    previous: tuple[str, int, int] | None = None
    try:
        with open_text(source) as handle:
            for line_number, raw in enumerate(handle, 1):
                text = raw.strip()
                if not text or text.startswith(("#", "track", "browser")):
                    continue
                fields = text.split("\t") if "\t" in text else text.split()
                if len(fields) < 3:
                    raise ValueError(f"line {line_number}: expected at least BED3")
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as exc:
                    raise ValueError(
                        f"line {line_number}: non-integer start/end"
                    ) from exc
                if start < 0 or end <= start:
                    raise ValueError(
                        f"line {line_number}: require 0 <= start < end"
                    )
                key = (fields[0], start, end)
                if require_sorted and previous is not None and key < previous:
                    raise ValueError(f"line {line_number}: records are not sorted")
                previous = key
                records += 1
                if max_records is not None and records >= max_records:
                    break
        if records == 0:
            raise ValueError("no interval records")
    except (OSError, EOFError, gzip.BadGzipFile, ValueError) as exc:
        return ValidationRow(kind, str(source), "FAIL", records, str(exc))
    detail = "sampled" if max_records is not None else "validated to EOF"
    return ValidationRow(kind, str(source.resolve()), "PASS", records, detail)


def _validate_chrom_sizes(path: str | Path) -> tuple[ValidationRow, dict[str, int]]:
    source = Path(path)
    try:
        rows = list(read_chrom_sizes_source(str(source)))
        if not rows:
            raise ValueError("no chromosome-size records")
        names = [name for name, _length in rows]
        if len(names) != len(set(names)):
            raise ValueError("duplicate chromosome names")
        if any(int(length) <= 0 for _name, length in rows):
            raise ValueError("chromosome lengths must be positive")
    except (OSError, ValueError) as exc:
        return ValidationRow("chrom_sizes", str(source), "FAIL", 0, str(exc)), {}
    return (
        ValidationRow(
            "chrom_sizes", str(source.resolve()), "PASS", len(rows), "valid names and lengths"
        ),
        dict(rows),
    )


def _validate_bam(path: str | Path, *, require_index: bool) -> tuple[ValidationRow, dict[str, int]]:
    source = Path(path)
    if not source.is_file():
        return ValidationRow("bam", str(source), "FAIL", 0, "file not found"), {}
    try:
        import pysam
    except ImportError:
        return ValidationRow("bam", str(source), "FAIL", 0, "pysam is not installed"), {}
    try:
        with pysam.AlignmentFile(str(source), "rb") as handle:
            references = list(handle.references)
            lengths = list(map(int, handle.lengths))
            if not references:
                raise ValueError("BAM header has no references")
            if require_index and not handle.has_index():
                raise ValueError("BAM index is required but unavailable")
            detail = "header and index valid" if handle.has_index() else "header valid; index absent"
    except (OSError, ValueError) as exc:
        return ValidationRow("bam", str(source), "FAIL", 0, str(exc)), {}
    return ValidationRow("bam", str(source.resolve()), "PASS", len(references), detail), dict(zip(references, lengths))


def _validate_fasta(path: str | Path) -> tuple[ValidationRow, dict[str, int]]:
    source = Path(path)
    if not source.is_file():
        return ValidationRow("fasta", str(source), "FAIL", 0, "file not found"), {}
    try:
        import pysam
    except ImportError:
        return ValidationRow("fasta", str(source), "FAIL", 0, "pysam is not installed"), {}
    try:
        with pysam.FastaFile(str(source)) as handle:
            references = list(handle.references)
            lengths = [int(handle.get_reference_length(name)) for name in references]
            if not references:
                raise ValueError("FASTA index has no references")
    except (OSError, ValueError) as exc:
        return ValidationRow("fasta", str(source), "FAIL", 0, str(exc)), {}
    return ValidationRow("fasta", str(source.resolve()), "PASS", len(references), "FASTA index readable"), dict(zip(references, lengths))


def _reference_compatibility(reference_sets: Sequence[tuple[str, dict[str, int]]]) -> ValidationRow | None:
    populated = [(label, rows) for label, rows in reference_sets if rows]
    if len(populated) < 2:
        return None
    baseline_label, baseline = populated[0]
    for label, rows in populated[1:]:
        shared: list[tuple[str, str]] = []
        for chrom in baseline:
            try:
                resolved = resolve_contig_name(
                    chrom, list(rows), source_label=label
                )
            except KeyError:
                continue
            shared.append((chrom, resolved))
        for baseline_chrom, other_chrom in shared:
            if int(baseline[baseline_chrom]) != int(rows[other_chrom]):
                return ValidationRow(
                    "reference_compatibility",
                    f"{baseline_label} vs {label}",
                    "FAIL",
                    len(shared),
                    f"length mismatch for {baseline_chrom}/{other_chrom}: "
                    f"{baseline[baseline_chrom]} vs {rows[other_chrom]}",
                )
        if not shared:
            return ValidationRow(
                "reference_compatibility",
                f"{baseline_label} vs {label}",
                "FAIL",
                0,
                "no compatible chromosome names in common",
            )
    return ValidationRow(
        "reference_compatibility",
        "all reference-bearing inputs",
        "PASS",
        len(populated),
        "shared chromosome lengths agree after conservative alias resolution",
    )


def _write_report(path: Path, rows: Iterable[ValidationRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["kind", "path", "status", "records", "detail"])
            for row in rows:
                writer.writerow([row.kind, row.path, row.status, row.records, row.detail])
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite validate-inputs",
        description="Validate suite inputs, compressed interval integrity and reference compatibility.",
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--bam", action="append", default=[], help="BAM input to validate; repeat for multiple files.")
    parser.add_argument("--fragments", action="append", default=[], help="Fragment BED/BED.gz input to validate; repeat for multiple files.")
    parser.add_argument("--bed", action="append", default=[], help="Annotation BED/BED.gz input to validate; repeat for multiple files.")
    parser.add_argument("--blacklist-bed", help="Blacklist BED/BED.gz input to validate.")
    parser.add_argument("--fasta", help="Indexed FASTA to validate and compare with other reference-bearing inputs.")
    parser.add_argument("--chrom-sizes", help="Two-column chromosome-size table to validate and compare with other inputs.")
    parser.add_argument("--require-bam-index", action="store_true", help="Fail BAM validation when an index is unavailable.")
    parser.add_argument("--require-sorted-fragments", action="store_true", help="Require fragment interval records to be coordinate sorted.")
    parser.add_argument(
        "--max-records",
        type=int,
        default=0,
        help="Maximum interval records sampled per file; 0 validates to EOF.",
    )
    parser.add_argument("--report", type=Path, help="Optional atomically written TSV containing every validation result.")
    return parser


def run(args: argparse.Namespace) -> int:
    if args.max_records < 0:
        raise ValueError("--max-records must be zero or greater")
    if not any((args.bam, args.fragments, args.bed, args.blacklist_bed, args.fasta, args.chrom_sizes)):
        raise ValueError("Provide at least one input to validate")
    reporter = ProgressReporter("validate-inputs")
    maximum = args.max_records or None
    rows: list[ValidationRow] = []
    reference_sets: list[tuple[str, dict[str, int]]] = []
    for path in args.bam:
        reporter.stage(f"Validating BAM: {path}")
        row, references = _validate_bam(path, require_index=args.require_bam_index)
        rows.append(row)
        reference_sets.append((f"BAM:{path}", references))
    for path in args.fragments:
        reporter.stage(f"Validating fragment intervals: {path}")
        rows.append(
            _validate_interval_file(
                path,
                kind="fragments",
                max_records=maximum,
                require_sorted=args.require_sorted_fragments,
            )
        )
    for path in args.bed:
        reporter.stage(f"Validating BED: {path}")
        rows.append(
            _validate_interval_file(
                path, kind="bed", max_records=maximum, require_sorted=False
            )
        )
    if args.blacklist_bed:
        reporter.stage(f"Validating blacklist: {args.blacklist_bed}")
        rows.append(
            _validate_interval_file(
                args.blacklist_bed,
                kind="blacklist",
                max_records=maximum,
                require_sorted=False,
            )
        )
    if args.fasta:
        reporter.stage(f"Validating FASTA: {args.fasta}")
        row, references = _validate_fasta(args.fasta)
        rows.append(row)
        reference_sets.append((f"FASTA:{args.fasta}", references))
    if args.chrom_sizes:
        reporter.stage(f"Validating chromosome sizes: {args.chrom_sizes}")
        row, references = _validate_chrom_sizes(args.chrom_sizes)
        rows.append(row)
        reference_sets.append((f"chrom-sizes:{args.chrom_sizes}", references))
    compatibility = _reference_compatibility(reference_sets)
    if compatibility is not None:
        rows.append(compatibility)
    if args.report:
        reporter.stage(f"Writing validation report: {args.report}")
        _write_report(args.report, rows)
    for row in rows:
        print(f"{row.status}\t{row.kind}\t{row.path}\t{row.records}\t{row.detail}")
    return 0 if all(row.status == "PASS" for row in rows) else 2


def validate_argv(argv: Sequence[str] | None = None) -> None:
    build_parser().parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    try:
        return run(parser.parse_args(argv))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
