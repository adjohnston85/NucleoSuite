#!/usr/bin/env python3
"""Mean-normalize BigWig signal or BED-family interval scores."""

from __future__ import annotations

import argparse
import gzip
import math
import tempfile
from pathlib import Path
from typing import Sequence, TextIO

import numpy as np

from nucleosuite.io.intervals import (
    convert_bed_to_bigbed,
    is_bigbed_path,
    open_interval_text,
)
from nucleosuite.output_naming import compact_parameter
from nucleosuite.progress import ProgressReporter

try:
    import pyBigWig
except ImportError:  # pragma: no cover - exercised through runtime validation
    pyBigWig = None


_CHUNK_BP = 1_000_000
_INTERVAL_FORMATS = ("bed", "bed.gz", "bigbed")
_OUTPUT_FORMATS = ("same", "bigwig", *_INTERVAL_FORMATS)
_HEADER_PREFIXES = ("#", "track", "browser")


def _require_pybigwig() -> None:
    if pyBigWig is None:
        raise RuntimeError(
            "mean-scale requires pyBigWig. Install it with "
            "'conda install -c bioconda pybigwig' or 'pip install pyBigWig'."
        )


def _input_kind(path: str | Path) -> str:
    lower = Path(path).name.lower()
    if lower.endswith((".bw", ".bigwig")):
        return "bigwig"
    if lower.endswith((".bed", ".bed.gz", ".bb", ".bigbed")):
        return "interval"
    raise ValueError(
        "mean-scale input must be BigWig (.bw/.bigWig), BED, BED.gz, or bigBed (.bb/.bigBed)"
    )


def _interval_format(path: str | Path) -> str:
    lower = Path(path).name.lower()
    if lower.endswith((".bb", ".bigbed")):
        return "bigbed"
    if lower.endswith(".bed.gz"):
        return "bed.gz"
    return "bed"


def _strip_input_suffix(path: str | Path) -> str:
    name = Path(path).name
    lower = name.lower()
    for suffix in (".bed.gz", ".bigwig", ".bigbed", ".bed", ".bw", ".bb"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _region_score_mean(path: str | Path, score_column: int) -> tuple[float, int, int]:
    """Return the unweighted mean of finite region scores."""
    if score_column < 1:
        raise ValueError("--score-column must be at least 1")
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    index = score_column - 1
    total = 0.0
    count = 0
    nonfinite = 0
    with open_interval_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(_HEADER_PREFIXES):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if index >= len(fields):
                raise ValueError(
                    f"{path}:{line_number}: score column {score_column} is absent "
                    f"from a {len(fields)}-column record"
                )
            try:
                value = float(fields[index])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid score in column {score_column}"
                ) from exc
            if not math.isfinite(value):
                nonfinite += 1
                continue
            total += value
            count += 1
    if count == 0:
        raise ValueError(f"No finite region scores were found in {path}")
    mean = total / count
    if not math.isfinite(mean) or mean == 0:
        raise ValueError("The region-score reference mean must be finite and non-zero")
    return float(mean), count, nonfinite


def _bigwig_nonzero_mean(handle, reporter: ProgressReporter | None = None) -> tuple[float, int]:
    """Return the exact mean across finite, non-zero BigWig bases."""
    chroms = handle.chroms()
    if not chroms:
        raise ValueError("Input BigWig contains no chromosomes")
    total = 0.0
    count = 0
    for chrom, length in chroms.items():
        if reporter is not None:
            reporter.reading_contig("BigWig mean", chrom)
        length = int(length)
        for start in range(0, length, _CHUNK_BP):
            end = min(length, start + _CHUNK_BP)
            values = np.asarray(handle.values(chrom, start, end, numpy=True), dtype=np.float64)
            mask = np.isfinite(values) & (values != 0)
            if np.any(mask):
                selected = values[mask]
                total += float(np.sum(selected, dtype=np.float64))
                count += int(selected.size)
    if count == 0:
        raise ValueError("Input BigWig contains no finite non-zero values")
    mean = total / count
    if not math.isfinite(mean) or mean == 0:
        raise ValueError("The BigWig-derived reference mean must be finite and non-zero")
    return float(mean), count


def bigwig_nonzero_mean(
    path: str | Path, reporter: ProgressReporter | None = None
) -> tuple[float, int]:
    """Return the finite non-zero mean for a BigWig path."""
    _require_pybigwig()
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    handle = pyBigWig.open(str(source))
    if handle is None:
        raise OSError(f"Could not open BigWig: {source}")
    try:
        return _bigwig_nonzero_mean(handle, reporter)
    finally:
        handle.close()


def _iter_scaled_intervals(handle, chrom: str, length: int, factor: float):
    """Yield clipped scaled BigWig intervals without materialising a chromosome."""
    for chunk_start in range(0, int(length), _CHUNK_BP):
        chunk_end = min(int(length), chunk_start + _CHUNK_BP)
        intervals = handle.intervals(chrom, chunk_start, chunk_end)
        if not intervals:
            continue
        for raw_start, raw_end, raw_value in intervals:
            start = max(int(raw_start), chunk_start)
            end = min(int(raw_end), chunk_end)
            if end <= start:
                continue
            value = float(raw_value)
            if not math.isfinite(value):
                continue
            yield start, end, value * factor


def _write_scaled_bigwig(
    input_handle,
    output_path: str | Path,
    *,
    reference_mean: float,
    scale: float,
    reporter: ProgressReporter | None = None,
) -> Path:
    _require_pybigwig()
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    partial.unlink(missing_ok=True)
    factor = float(scale) / float(reference_mean)
    chroms = input_handle.chroms()
    writer = pyBigWig.open(str(partial), "w")
    committed = False
    try:
        writer.addHeader([(str(chrom), int(length)) for chrom, length in chroms.items()])
        for chrom, length in chroms.items():
            if reporter is not None:
                reporter.reading_contig("mean-scaled BigWig", chrom)
            starts: list[int] = []
            ends: list[int] = []
            values: list[float] = []
            for start, end, value in _iter_scaled_intervals(
                input_handle, str(chrom), int(length), factor
            ):
                starts.append(start)
                ends.append(end)
                values.append(float(value))
                if len(starts) >= 100_000:
                    writer.addEntries(
                        [str(chrom)] * len(starts), starts, ends=ends, values=values
                    )
                    starts.clear(); ends.clear(); values.clear()
            if starts:
                writer.addEntries(
                    [str(chrom)] * len(starts), starts, ends=ends, values=values
                )
        writer.close()
        writer = None
        partial.replace(output)
        committed = True
    finally:
        if writer is not None:
            try:
                writer.close()
            except Exception:
                pass
        if not committed:
            partial.unlink(missing_ok=True)
    return output


def scale_bigwig_by_reference(
    input_path: str | Path,
    reference_path: str | Path,
    output_path: str | Path,
    *,
    scale: float = 1.0,
    reporter: ProgressReporter | None = None,
) -> tuple[Path, float, int]:
    """Scale one BigWig by the finite non-zero mean of another BigWig."""
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be a finite value greater than zero")
    reference_mean, count = bigwig_nonzero_mean(reference_path, reporter)
    _require_pybigwig()
    source = pyBigWig.open(str(input_path))
    if source is None:
        raise OSError(f"Could not open BigWig: {input_path}")
    try:
        output = _write_scaled_bigwig(
            source,
            output_path,
            reference_mean=reference_mean,
            scale=scale,
            reporter=reporter,
        )
    finally:
        source.close()
    return output, reference_mean, count


def _infer_bigbed_chrom_sizes(path: str | Path) -> dict[str, int]:
    """Read chromosome sizes embedded in a bigBed."""
    _require_pybigwig()
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


def _open_text_output(path: Path, *, compressed: bool) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compressed:
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("wt", encoding="utf-8", newline="")


def _format_score(value: float, integer_scores: bool) -> str:
    if integer_scores:
        return str(int(round(value)))
    return f"{value:.12g}"


def _effective_interval_controls(
    *,
    output_format: str,
    integer_scores: bool,
    clamp_min: float | None,
    clamp_max: float | None,
) -> tuple[bool, float | None, float | None]:
    """Return score controls, enforcing the standard bigBed score domain."""
    if output_format == "bigbed":
        integer_scores = True
        clamp_min = 0.0 if clamp_min is None else max(0.0, float(clamp_min))
        clamp_max = 1000.0 if clamp_max is None else min(1000.0, float(clamp_max))
    if clamp_min is not None and (not math.isfinite(clamp_min)):
        raise ValueError("--clamp-min must be finite")
    if clamp_max is not None and (not math.isfinite(clamp_max)):
        raise ValueError("--clamp-max must be finite")
    if clamp_min is not None and clamp_max is not None and clamp_max < clamp_min:
        raise ValueError("--clamp-max must be greater than or equal to --clamp-min")
    return integer_scores, clamp_min, clamp_max


def _write_scaled_intervals(
    input_path: str | Path,
    output_path: str | Path,
    *,
    reference_mean: float,
    scale: float,
    score_column: int,
    output_format: str,
    integer_scores: bool,
    clamp_min: float | None,
    clamp_max: float | None,
    chrom_sizes: str | Path | dict[str, int] | None,
    reporter: ProgressReporter | None = None,
) -> Path:
    """Scale one BED score column while preserving every other BED field."""
    source = Path(input_path)
    destination = Path(output_path)
    integer_scores, clamp_min, clamp_max = _effective_interval_controls(
        output_format=output_format,
        integer_scores=integer_scores,
        clamp_min=clamp_min,
        clamp_max=clamp_max,
    )
    factor = float(scale) / float(reference_mean)
    score_index = score_column - 1

    temp_dir: tempfile.TemporaryDirectory | None = None
    if output_format == "bigbed":
        temp_dir = tempfile.TemporaryDirectory(prefix="nucleosuite_meanscale_")
        text_output = Path(temp_dir.name) / "scaled.bed"
        compressed = False
    else:
        text_output = destination
        compressed = output_format == "bed.gz"

    seen_contigs: set[str] = set()
    try:
        with open_interval_text(source) as src, _open_text_output(
            text_output, compressed=compressed
        ) as dst:
            for line_number, raw in enumerate(src, 1):
                text = raw.rstrip("\n")
                stripped = text.strip()
                if not stripped:
                    continue
                if stripped.startswith(_HEADER_PREFIXES):
                    if output_format != "bigbed":
                        dst.write(text + "\n")
                    continue
                fields = text.split("\t") if "\t" in text else text.split()
                if len(fields) < 3:
                    raise ValueError(f"{source}:{line_number}: expected at least BED3")
                if score_index >= len(fields):
                    raise ValueError(
                        f"{source}:{line_number}: score column {score_column} is absent "
                        f"from a {len(fields)}-column record"
                    )
                try:
                    value = float(fields[score_index])
                except ValueError as exc:
                    raise ValueError(
                        f"{source}:{line_number}: invalid score in column {score_column}"
                    ) from exc
                if not math.isfinite(value):
                    raise ValueError(
                        f"{source}:{line_number}: score in column {score_column} is not finite"
                    )
                scaled = value * factor
                if clamp_min is not None:
                    scaled = max(float(clamp_min), scaled)
                if clamp_max is not None:
                    scaled = min(float(clamp_max), scaled)
                fields[score_index] = _format_score(scaled, integer_scores)
                dst.write("\t".join(fields) + "\n")
                chrom = fields[0]
                if reporter is not None and chrom not in seen_contigs:
                    seen_contigs.add(chrom)
                    reporter.reading_contig("mean-scaled intervals", chrom)

        if output_format == "bigbed":
            sizes = chrom_sizes
            if sizes is None and is_bigbed_path(source):
                sizes = _infer_bigbed_chrom_sizes(source)
            if sizes is None:
                raise ValueError(
                    "--chrom-sizes is required for bigBed output when chromosome sizes "
                    "cannot be inherited from a bigBed input"
                )
            converted = convert_bed_to_bigbed(
                text_output,
                sizes,
                destination,
                bigbed_score_multiplier=1.0,
            )
            if converted is None:
                return Path(str(destination) + ".empty")
        return destination
    finally:
        if temp_dir is not None:
            temp_dir.cleanup()


def _reference_token(
    mode: str,
    *,
    reference_mean: float,
    regions: Path | None,
    score_column: int,
) -> str:
    if mode == "bigwig-nonzero-mean":
        return "bwnonzero"
    if mode == "reference-bigwig-nonzero-mean":
        return "reference-bwnonzero"
    if mode == "input-score-mean":
        return f"scores-col{score_column}"
    if mode == "region-score-mean":
        region_token = compact_parameter(_strip_input_suffix(regions) if regions else "regions")
        return f"regions-{region_token}-col{score_column}"
    return f"mean-{compact_parameter(reference_mean)}"


def _default_output(
    input_path: Path,
    *,
    input_kind: str = "bigwig",
    output_format: str = "bigwig",
    mode: str,
    reference_mean: float,
    scale: float,
    regions: Path | None,
    score_column: int,
    integer_scores: bool = False,
    clamp_min: float | None = None,
    clamp_max: float | None = None,
) -> Path:
    stem = _strip_input_suffix(input_path)
    source_token = _reference_token(
        mode,
        reference_mean=reference_mean,
        regions=regions,
        score_column=score_column,
    )
    tokens = [f"{stem}_meanscale_{source_token}_x{compact_parameter(scale)}"]
    if input_kind == "interval":
        if integer_scores or output_format == "bigbed":
            tokens.append("int")
        if clamp_min is not None or output_format == "bigbed":
            effective = 0.0 if clamp_min is None and output_format == "bigbed" else clamp_min
            tokens.append(f"min{compact_parameter(effective)}")
        if clamp_max is not None or output_format == "bigbed":
            effective = 1000.0 if clamp_max is None and output_format == "bigbed" else clamp_max
            tokens.append(f"max{compact_parameter(effective)}")
    extension = {
        "bigwig": ".bw",
        "bed": ".bed",
        "bed.gz": ".bed.gz",
        "bigbed": ".bb",
    }[output_format]
    return input_path.with_name("_".join(tokens) + extension)


def _write_summary(
    output_path: Path,
    *,
    input_path: Path,
    input_kind: str,
    output_format: str,
    mode: str,
    reference_mean: float,
    scale: float,
    region_path: Path | None,
    score_column: int,
    contributing_values: int | None,
    nonfinite_region_scores: int | None,
    integer_scores: bool,
    clamp_min: float | None,
    clamp_max: float | None,
) -> Path:
    name = output_path.name
    lower = name.lower()
    for suffix in (".bed.gz", ".bigwig", ".bigbed", ".bed", ".bw", ".bb", ".empty"):
        if lower.endswith(suffix):
            name = name[: -len(suffix)]
            break
    summary = output_path.with_name(name + "_mean_scale_summary.tsv")
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"input\t{input_path}\n")
        handle.write(f"output\t{output_path}\n")
        handle.write(f"input_kind\t{input_kind}\n")
        handle.write(f"output_format\t{output_format}\n")
        handle.write(f"reference_mode\t{mode}\n")
        handle.write(f"reference_mean\t{reference_mean:.12g}\n")
        handle.write(f"scale\t{scale:.12g}\n")
        handle.write(f"multiplier\t{(scale / reference_mean):.12g}\n")
        if input_kind == "interval" or region_path is not None:
            handle.write(f"score_column\t{score_column}\n")
        if region_path is not None:
            handle.write(f"regions\t{region_path}\n")
        if contributing_values is not None:
            key = (
                "finite_nonzero_bigwig_bases"
                if mode in {"bigwig-nonzero-mean", "reference-bigwig-nonzero-mean"}
                else "finite_region_scores"
            )
            handle.write(f"{key}\t{contributing_values}\n")
        if nonfinite_region_scores is not None:
            handle.write(f"nonfinite_region_scores_excluded\t{nonfinite_region_scores}\n")
        if input_kind == "interval":
            handle.write(f"integer_scores\t{str(integer_scores).lower()}\n")
            handle.write(f"clamp_min\t{'' if clamp_min is None else f'{clamp_min:.12g}'}\n")
            handle.write(f"clamp_max\t{'' if clamp_max is None else f'{clamp_max:.12g}'}\n")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite mean-scale",
        description=(
            "Scale BigWig signal or BED-family scores as value / reference_mean × scale. "
            "BigWig input defaults to the finite non-zero signal mean; BED, BED.gz and "
            "bigBed input default to the mean of the selected score column."
        ),
    )
    parser.add_argument(
        "input",
        help="Input BigWig, BED, BED.gz, or bigBed to mean-scale.",
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--reference-mean", "--normalization-mean", dest="reference_mean", type=float,
        help="Use this finite non-zero reference mean directly.",
    )
    source.add_argument(
        "--regions",
        help="BED, BED.gz or bigBed whose score-column mean defines the reference mean.",
    )
    source.add_argument(
        "--reference-bigwig",
        help=(
            "BigWig whose finite non-zero mean defines the reference mean. Useful "
            "for scaling a signal by a chosen non-negative reference track."
        ),
    )
    parser.add_argument(
        "--score-column", type=int, default=5,
        help="1-based score column for BED-family input/reference files (default: 5).",
    )
    parser.add_argument(
        "--scale", type=float, default=100.0,
        help="Target mean after scaling (default: 100).",
    )
    parser.add_argument(
        "--integer-scores", action="store_true",
        help="Round scaled BED-family scores to integers. Automatic for bigBed output.",
    )
    parser.add_argument(
        "--clamp-min", type=float,
        help="Clamp scaled BED-family scores to at least this value.",
    )
    parser.add_argument(
        "--clamp-max", type=float,
        help="Clamp scaled BED-family scores to at most this value.",
    )
    parser.add_argument(
        "--output-format", choices=_OUTPUT_FORMATS, default="same",
        help=(
            "Output format (default: same as input). BED-family inputs may be written "
            "as bed, bed.gz, or bigbed; BigWig input writes bigwig."
        ),
    )
    parser.add_argument(
        "--chrom-sizes",
        help=(
            "Chromosome sizes for bigBed output. Existing bigBed input supplies its "
            "own chromosome sizes when possible."
        ),
    )
    parser.add_argument(
        "--output", "-o",
        help="Output path. Default name records the input basename, reference mode, and scale.",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    input_kind = _input_kind(input_path)
    if args.score_column < 1:
        raise ValueError("--score-column must be at least 1")
    if not math.isfinite(args.scale) or args.scale <= 0:
        raise ValueError("--scale must be a finite value greater than zero")
    for option, value in (("--clamp-min", args.clamp_min), ("--clamp-max", args.clamp_max)):
        if value is not None and not math.isfinite(value):
            raise ValueError(f"{option} must be finite")
    if args.clamp_min is not None and args.clamp_max is not None and args.clamp_max < args.clamp_min:
        raise ValueError("--clamp-max must be greater than or equal to --clamp-min")

    if input_kind == "bigwig":
        if args.output_format not in ("same", "bigwig"):
            raise ValueError("BigWig input can only be written as BigWig")
        if args.integer_scores or args.clamp_min is not None or args.clamp_max is not None:
            raise ValueError("--integer-scores and --clamp-min/--clamp-max apply only to BED-family input")
        output_format = "bigwig"
        _require_pybigwig()
    else:
        if args.output_format == "bigwig":
            raise ValueError("BED-family input cannot be written as BigWig")
        output_format = _interval_format(input_path) if args.output_format == "same" else args.output_format

    reporter = ProgressReporter("mean-scale")
    region_path = Path(args.regions) if args.regions else None
    reference_bigwig_path = (
        Path(args.reference_bigwig) if args.reference_bigwig else None
    )
    contributing_values: int | None = None
    nonfinite_region_scores: int | None = None

    if args.reference_mean is not None:
        reference_mean = float(args.reference_mean)
        if not math.isfinite(reference_mean) or reference_mean == 0:
            raise ValueError("--reference-mean must be finite and non-zero")
        mode = "supplied-reference-mean"
    elif reference_bigwig_path is not None:
        if input_kind != "bigwig":
            raise ValueError("--reference-bigwig requires BigWig input")
        reporter.stage("Calculating reference mean from positive-score BigWig")
        reference_mean, contributing_values = bigwig_nonzero_mean(
            reference_bigwig_path, reporter
        )
        mode = "reference-bigwig-nonzero-mean"
    elif region_path is not None:
        reporter.stage("Calculating reference mean from region scores")
        reference_mean, contributing_values, nonfinite_region_scores = _region_score_mean(
            region_path, args.score_column
        )
        mode = "region-score-mean"
    elif input_kind == "interval":
        reporter.stage("Calculating reference mean from input interval scores")
        reference_mean, contributing_values, nonfinite_region_scores = _region_score_mean(
            input_path, args.score_column
        )
        mode = "input-score-mean"
    else:
        reporter.stage("Calculating reference mean from finite non-zero BigWig values")
        reader = pyBigWig.open(str(input_path))
        try:
            reference_mean, contributing_values = _bigwig_nonzero_mean(reader, reporter)
        finally:
            reader.close()
        mode = "bigwig-nonzero-mean"

    effective_integer = bool(args.integer_scores)
    effective_clamp_min = args.clamp_min
    effective_clamp_max = args.clamp_max
    if input_kind == "interval":
        effective_integer, effective_clamp_min, effective_clamp_max = _effective_interval_controls(
            output_format=output_format,
            integer_scores=effective_integer,
            clamp_min=effective_clamp_min,
            clamp_max=effective_clamp_max,
        )

    output_path = Path(args.output) if args.output else _default_output(
        input_path,
        input_kind=input_kind,
        output_format=output_format,
        mode=mode,
        reference_mean=reference_mean,
        scale=args.scale,
        regions=region_path,
        score_column=args.score_column,
        integer_scores=effective_integer,
        clamp_min=effective_clamp_min,
        clamp_max=effective_clamp_max,
    )
    if output_path.resolve() == input_path.resolve():
        raise ValueError("Output must differ from the input")

    reporter.stage(f"Scaling by {args.scale:.12g} / {reference_mean:.12g}")
    if input_kind == "bigwig":
        reader = pyBigWig.open(str(input_path))
        try:
            _write_scaled_bigwig(
                reader,
                output_path,
                reference_mean=reference_mean,
                scale=args.scale,
                reporter=reporter,
            )
        finally:
            reader.close()
    else:
        output_path = _write_scaled_intervals(
            input_path,
            output_path,
            reference_mean=reference_mean,
            scale=args.scale,
            score_column=args.score_column,
            output_format=output_format,
            integer_scores=effective_integer,
            clamp_min=effective_clamp_min,
            clamp_max=effective_clamp_max,
            chrom_sizes=args.chrom_sizes,
            reporter=reporter,
        )

    summary = _write_summary(
        output_path,
        input_path=input_path,
        input_kind=input_kind,
        output_format=output_format,
        mode=mode,
        reference_mean=reference_mean,
        scale=args.scale,
        region_path=region_path,
        score_column=args.score_column,
        contributing_values=contributing_values,
        nonfinite_region_scores=nonfinite_region_scores,
        integer_scores=effective_integer,
        clamp_min=effective_clamp_min,
        clamp_max=effective_clamp_max,
    )
    if reference_bigwig_path is not None:
        with summary.open("at", encoding="utf-8") as handle:
            handle.write(f"reference_bigwig\t{reference_bigwig_path}\n")
    print(f"Reference mean: {reference_mean:.12g}")
    print(f"Scale: {args.scale:.12g}")
    print(f"Wrote: {output_path}")
    print(f"Wrote: {summary}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
