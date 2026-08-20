#!/usr/bin/env python3
"""Scale a BigWig relative to a reference mean."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from nucleosuite.io.intervals import open_interval_text
from nucleosuite.output_naming import compact_parameter
from nucleosuite.progress import ProgressReporter

try:
    import pyBigWig
except ImportError:  # pragma: no cover - exercised through runtime validation
    pyBigWig = None


_CHUNK_BP = 1_000_000


def _require_pybigwig() -> None:
    if pyBigWig is None:
        raise RuntimeError(
            "mean-scale requires pyBigWig. Install it with "
            "'conda install -c bioconda pybigwig' or 'pip install pyBigWig'."
        )


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
            if not text or text.startswith(("#", "track", "browser")):
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


def _iter_scaled_intervals(handle, chrom: str, length: int, factor: float):
    """Yield clipped scaled intervals without materialising a whole chromosome."""
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


def _default_output(
    input_path: Path,
    *,
    mode: str,
    reference_mean: float,
    scale: float,
    regions: Path | None,
    score_column: int,
) -> Path:
    base = input_path.name
    lower = base.lower()
    if lower.endswith(".bigwig"):
        stem = base[:-7]
    elif lower.endswith(".bw"):
        stem = base[:-3]
    else:
        stem = input_path.stem
    if mode == "bigwig-nonzero-mean":
        source_token = "bwnonzero"
    elif mode == "region-score-mean":
        region_token = compact_parameter(regions.stem if regions is not None else "regions")
        source_token = f"regions-{region_token}-col{score_column}"
    else:
        source_token = f"mean-{compact_parameter(reference_mean)}"
    filename = (
        f"{stem}_meanscale_{source_token}_x{compact_parameter(scale)}.bw"
    )
    return input_path.with_name(filename)


def _write_summary(
    output_bigwig: Path,
    *,
    input_bigwig: Path,
    mode: str,
    reference_mean: float,
    scale: float,
    region_path: Path | None,
    score_column: int,
    contributing_values: int | None,
    nonfinite_region_scores: int | None,
) -> Path:
    summary = output_bigwig.with_name(output_bigwig.stem + "_mean_scale_summary.tsv")
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"input_bigwig\t{input_bigwig}\n")
        handle.write(f"output_bigwig\t{output_bigwig}\n")
        handle.write(f"reference_mode\t{mode}\n")
        handle.write(f"reference_mean\t{reference_mean:.12g}\n")
        handle.write(f"scale\t{scale:.12g}\n")
        handle.write(f"multiplier\t{(scale / reference_mean):.12g}\n")
        if region_path is not None:
            handle.write(f"regions\t{region_path}\n")
            handle.write(f"score_column\t{score_column}\n")
        if contributing_values is not None:
            key = "finite_region_scores" if mode == "region-score-mean" else "finite_nonzero_bigwig_bases"
            handle.write(f"{key}\t{contributing_values}\n")
        if nonfinite_region_scores is not None:
            handle.write(f"nonfinite_region_scores_excluded\t{nonfinite_region_scores}\n")
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite mean-scale",
        description=(
            "Scale a BigWig as value / reference_mean × scale. The reference mean "
            "can be supplied directly, calculated from region scores, or calculated "
            "from finite non-zero BigWig values."
        ),
    )
    parser.add_argument("bigwig", help="Input BigWig to scale.")
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--reference-mean", "--normalization-mean", dest="reference_mean", type=float,
        help="Use this finite non-zero reference mean directly.",
    )
    source.add_argument(
        "--regions",
        help="BED, BED.gz or bigBed whose region scores define the reference mean.",
    )
    parser.add_argument(
        "--score-column", type=int, default=5,
        help="1-based score column used with --regions (default: 5).",
    )
    parser.add_argument(
        "--scale", type=float, default=100.0,
        help="Scale applied after division by the reference mean (default: 100).",
    )
    parser.add_argument(
        "--output", "-o",
        help="Output BigWig. Default name records the reference mode and scale.",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    _require_pybigwig()
    input_path = Path(args.bigwig)
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    if args.score_column < 1:
        raise ValueError("--score-column must be at least 1")
    if not math.isfinite(args.scale) or args.scale <= 0:
        raise ValueError("--scale must be a finite value greater than zero")

    reporter = ProgressReporter("mean-scale")
    region_path = Path(args.regions) if args.regions else None
    contributing_values: int | None = None
    nonfinite_region_scores: int | None = None

    if args.reference_mean is not None:
        reference_mean = float(args.reference_mean)
        if not math.isfinite(reference_mean) or reference_mean == 0:
            raise ValueError("--reference-mean must be finite and non-zero")
        mode = "supplied-reference-mean"
    elif region_path is not None:
        reporter.stage("Calculating reference mean from region scores")
        reference_mean, contributing_values, nonfinite_region_scores = _region_score_mean(
            region_path, args.score_column
        )
        mode = "region-score-mean"
    else:
        reporter.stage("Calculating reference mean from finite non-zero BigWig values")
        reader = pyBigWig.open(str(input_path))
        try:
            reference_mean, contributing_values = _bigwig_nonzero_mean(reader, reporter)
        finally:
            reader.close()
        mode = "bigwig-nonzero-mean"

    output_path = Path(args.output) if args.output else _default_output(
        input_path,
        mode=mode,
        reference_mean=reference_mean,
        scale=args.scale,
        regions=region_path,
        score_column=args.score_column,
    )
    if output_path.resolve() == input_path.resolve():
        raise ValueError("Output BigWig must differ from the input BigWig")

    reporter.stage(
        f"Scaling BigWig by {args.scale:.12g} / {reference_mean:.12g}"
    )
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

    summary = _write_summary(
        output_path,
        input_bigwig=input_path,
        mode=mode,
        reference_mean=reference_mean,
        scale=args.scale,
        region_path=region_path,
        score_column=args.score_column,
        contributing_values=contributing_values,
        nonfinite_region_scores=nonfinite_region_scores,
    )
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
