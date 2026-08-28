#!/usr/bin/env python3
"""Empirical peak FDR from observed and fragment-randomized peak callsets."""

from __future__ import annotations

import argparse
import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, TextIO

import numpy as np

from nucleosuite.io.intervals import open_interval_text
from nucleosuite.output_naming import compact_parameter
from nucleosuite.progress import ProgressReporter


_HEADER_PREFIXES = ("#", "track", "browser")


@dataclass(frozen=True)
class PeakRow:
    """One BED-family record retained verbatim apart from appended annotations."""

    fields: tuple[str, ...]
    score: float
    source_line: int


@dataclass(frozen=True)
class PeakFdrResult:
    """Paths and counts produced by one empirical-FDR comparison."""

    annotated_path: Path
    significant_path: Path | None
    summary_path: Path
    sample_peaks: int
    random_peak_sets: int
    random_peaks: int
    significant_peaks: int | None


def _read_peak_rows(
    path: str | Path,
    score_column: int,
    *,
    allow_empty: bool = False,
) -> list[PeakRow]:
    if score_column < 1:
        raise ValueError("--score-column must be at least 1")
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    index = score_column - 1
    rows: list[PeakRow] = []
    with open_interval_text(source) as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.rstrip("\n")
            stripped = text.strip()
            if not stripped or stripped.startswith(_HEADER_PREFIXES):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < 3:
                raise ValueError(f"{source}:{line_number}: expected at least BED3")
            if index >= len(fields):
                raise ValueError(
                    f"{source}:{line_number}: score column {score_column} is absent "
                    f"from a {len(fields)}-column record"
                )
            try:
                score = float(fields[index])
            except ValueError as exc:
                raise ValueError(
                    f"{source}:{line_number}: invalid score in column {score_column}"
                ) from exc
            if not math.isfinite(score):
                raise ValueError(
                    f"{source}:{line_number}: score in column {score_column} is not finite"
                )
            if score < 0:
                raise ValueError(
                    f"{source}:{line_number}: empirical peak FDR requires non-negative scores"
                )
            rows.append(PeakRow(tuple(fields), score, line_number))
    if not rows and not allow_empty:
        raise ValueError(f"No scored peak records were found in {source}")
    return rows


def empirical_peak_pvalues(
    sample_scores: Sequence[float] | np.ndarray,
    randomized_score_sets: Sequence[Sequence[float] | np.ndarray],
) -> np.ndarray:
    """Return pooled empirical upper-tail p-values for observed peak scores.

    The randomized callsets are pooled as positional-null peak scores. For an
    observed score ``s``, the p-value is ``(1 + R(score >= s)) / (1 + N_R)``.
    The +1 pseudocount prevents a zero p-value when an observed score exceeds
    every randomized peak.
    """

    observed = np.asarray(sample_scores, dtype=np.float64)
    if observed.ndim != 1 or observed.size == 0:
        raise ValueError("sample_scores must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(observed)) or np.any(observed < 0):
        raise ValueError("sample_scores must contain finite non-negative values")
    if not randomized_score_sets:
        raise ValueError("At least one randomized peak score set is required")
    random_arrays = [np.asarray(values, dtype=np.float64) for values in randomized_score_sets]
    for values in random_arrays:
        if values.ndim != 1:
            raise ValueError("Randomized peak scores must be one-dimensional")
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("Randomized peak scores must be finite and non-negative")
    pooled = np.concatenate(random_arrays) if random_arrays else np.empty(0, dtype=np.float64)
    if pooled.size == 0:
        return np.ones(observed.size, dtype=np.float64)
    pooled.sort()
    counts = pooled.size - np.searchsorted(pooled, observed, side="left")
    return (1.0 + counts.astype(np.float64)) / (1.0 + float(pooled.size))


def empirical_peak_qvalues(
    sample_scores: Sequence[float] | np.ndarray,
    randomized_score_sets: Sequence[Sequence[float] | np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return monotonic q-values and threshold counts for sample peak scores.

    At threshold ``s``, the estimated false-discovery rate is::

        (1 + sum_b R_b(score >= s)) / (B * S(score >= s))

    where ``B`` is the number of independently randomized peak callsets.  The
    returned q-value for each observed peak is the minimum estimated FDR among
    all score thresholds that retain that peak.
    """

    observed = np.asarray(sample_scores, dtype=np.float64)
    if observed.ndim != 1 or observed.size == 0:
        raise ValueError("sample_scores must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(observed)) or np.any(observed < 0):
        raise ValueError("sample_scores must contain finite non-negative values")
    if not randomized_score_sets:
        raise ValueError("At least one randomized peak score set is required")
    random_arrays = [np.asarray(values, dtype=np.float64) for values in randomized_score_sets]
    for values in random_arrays:
        if values.ndim != 1:
            raise ValueError("Randomized peak scores must be one-dimensional")
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("Randomized peak scores must be finite and non-negative")

    thresholds = np.unique(observed)[::-1]
    sample_sorted = np.sort(observed)
    random_sorted = [np.sort(values) for values in random_arrays]
    sample_counts = observed.size - np.searchsorted(sample_sorted, thresholds, side="left")
    random_counts = np.zeros(thresholds.size, dtype=np.float64)
    for values in random_sorted:
        random_counts += values.size - np.searchsorted(values, thresholds, side="left")

    fdr = (1.0 + random_counts) / (
        float(len(random_arrays)) * np.maximum(sample_counts.astype(float), 1.0)
    )
    fdr = np.clip(fdr, 0.0, 1.0)
    # thresholds are high to low.  A peak can be retained at its own threshold
    # or any less-stringent threshold below it, so take a reverse cumulative min.
    threshold_q = np.minimum.accumulate(fdr[::-1])[::-1]
    threshold_index = {float(score): index for index, score in enumerate(thresholds)}
    qvalues = np.asarray(
        [threshold_q[threshold_index[float(score)]] for score in observed],
        dtype=np.float64,
    )
    return qvalues, sample_counts.astype(np.int64), random_counts / len(random_arrays)


def _strip_bed_suffix(path: str | Path) -> str:
    name = Path(path).name
    lower = name.lower()
    for suffix in (".bed.gz", ".bigbed", ".bed", ".bb"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("wt", encoding="utf-8", newline="")


def _default_annotated_path(sample_path: Path, output_prefix: str | Path | None) -> Path:
    if output_prefix is not None:
        return Path(f"{output_prefix}_empirical_peak_fdr.bed")
    return sample_path.with_name(f"{_strip_bed_suffix(sample_path)}_empirical_peak_fdr.bed")


def annotate_peak_fdr(
    sample_peaks: str | Path,
    randomized_peaks: Sequence[str | Path],
    *,
    score_column: int = 5,
    fdr_threshold: float | None = None,
    output_prefix: str | Path | None = None,
    output_path: str | Path | None = None,
) -> PeakFdrResult:
    """Append empirical p-value and FDR to every observed BED record and optionally filter."""

    if fdr_threshold is not None and (
        not math.isfinite(fdr_threshold) or not 0.0 <= fdr_threshold <= 1.0
    ):
        raise ValueError("--fdr must be between 0 and 1")
    if not randomized_peaks:
        raise ValueError("At least one randomized peak BED is required")

    reporter = ProgressReporter("empirical-peak-fdr")
    reporter.stage("Reading observed peaks")
    observed_rows = _read_peak_rows(sample_peaks, score_column, allow_empty=True)
    random_rows: list[list[PeakRow]] = []
    for path in randomized_peaks:
        reporter.file_start("randomized peaks", path)
        random_rows.append(_read_peak_rows(path, score_column, allow_empty=True))

    if observed_rows:
        observed_scores = [row.score for row in observed_rows]
        randomized_scores = [[row.score for row in rows] for rows in random_rows]
        pvalues = empirical_peak_pvalues(observed_scores, randomized_scores)
        qvalues, _sample_counts, _mean_random_counts = empirical_peak_qvalues(
            observed_scores, randomized_scores
        )
    else:
        pvalues = np.empty(0, dtype=np.float64)
        qvalues = np.empty(0, dtype=np.float64)
    annotated = (
        Path(output_path)
        if output_path is not None
        else _default_annotated_path(Path(sample_peaks), output_prefix)
    )
    reporter.stage("Writing observed peaks with empirical p-values and FDR")
    with _open_output(annotated) as handle:
        for row, pvalue, qvalue in zip(observed_rows, pvalues, qvalues):
            handle.write(
                "\t".join((*row.fields, f"{pvalue:.12g}", f"{qvalue:.12g}")) + "\n"
            )

    significant_path: Path | None = None
    significant_count: int | None = None
    if fdr_threshold is not None:
        significant_count = int(np.count_nonzero(qvalues <= fdr_threshold))
        base = (
            str(output_prefix)
            if output_prefix is not None
            else str(annotated.with_suffix(""))
        )
        significant_path = Path(
            f"{base}_fdr{compact_parameter(fdr_threshold)}_significant.bed"
        )
        with _open_output(significant_path) as handle:
            for row, pvalue, qvalue in zip(observed_rows, pvalues, qvalues):
                if qvalue <= fdr_threshold:
                    handle.write(
                        "\t".join((*row.fields, f"{pvalue:.12g}", f"{qvalue:.12g}")) + "\n"
                    )

    summary_path = annotated.with_name(
        f"{_strip_bed_suffix(annotated)}_summary.tsv"
    )
    with summary_path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"sample_peaks\t{Path(sample_peaks)}\n")
        for index, path in enumerate(randomized_peaks, 1):
            handle.write(f"randomized_peaks_{index}\t{Path(path)}\n")
        handle.write(f"score_column\t{score_column}\n")
        handle.write(f"sample_peak_count\t{len(observed_rows)}\n")
        handle.write(f"randomized_peak_sets\t{len(random_rows)}\n")
        handle.write(f"randomized_peak_count_total\t{sum(map(len, random_rows))}\n")
        handle.write(f"fdr_pseudocount\t1\n")
        handle.write(f"annotated_output\t{annotated}\n")
        if fdr_threshold is not None:
            handle.write(f"fdr_threshold\t{fdr_threshold:.12g}\n")
            handle.write(f"significant_peak_count\t{significant_count}\n")
            handle.write(f"significant_output\t{significant_path}\n")

    return PeakFdrResult(
        annotated_path=annotated,
        significant_path=significant_path,
        summary_path=summary_path,
        sample_peaks=len(observed_rows),
        random_peak_sets=len(random_rows),
        random_peaks=sum(map(len, random_rows)),
        significant_peaks=significant_count,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite empirical-peak-fdr",
        description=(
            "Assign empirical p-values and FDR values to observed peaks using one or more "
            "fragment-randomized peak callsets as positional null controls."
        ),
    )
    parser.add_argument("sample_peaks", help="Observed sample peak BED/BED.gz/bigBed.")
    parser.add_argument(
        "randomized_peaks", nargs="+",
        help="One or more matched fragment-randomized peak callsets.",
    )
    parser.add_argument(
        "--score-column", type=int, default=5,
        help="1-based peak-score column in every input (default: 5).",
    )
    parser.add_argument(
        "--fdr", type=float,
        help=(
            "Optional FDR cutoff. Every peak is always written with its empirical "
            "p-value and FDR; when set, an additional filtered BED is written."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        help="Output prefix. Default: observed input basename in its input directory.",
    )
    parser.add_argument(
        "--output",
        help="Exact path for the complete annotated BED; overrides its default name.",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    result = annotate_peak_fdr(
        args.sample_peaks,
        args.randomized_peaks,
        score_column=args.score_column,
        fdr_threshold=args.fdr,
        output_prefix=args.output_prefix,
        output_path=args.output,
    )
    print(f"annotated_peaks\t{result.annotated_path}")
    if result.significant_path is not None:
        print(f"significant_peaks\t{result.significant_path}")
    print(f"summary\t{result.summary_path}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
