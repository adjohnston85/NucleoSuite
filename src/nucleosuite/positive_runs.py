#!/usr/bin/env python3
"""Measure uninterrupted positive-score runs in a single BigWig track.

A run begins when the per-base BigWig value is strictly greater than the
selected threshold and ends at a zero, negative, missing, or non-finite value,
a genomic gap, a selected-region boundary, or a contig boundary. Coordinates
are reported as zero-based half-open intervals.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nucleosuite.plotting import configure_unique_category_cycle
configure_unique_category_cycle()
import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.progress import ProgressReporter


@dataclass
class ScanSummary:
    scanned_bases: int = 0
    positive_bases: int = 0
    nonpositive_bases: int = 0
    missing_bases: int = 0
    blacklisted_bases: int = 0
    total_runs_observed: int = 0
    total_runs_retained: int = 0
    retained_positive_bases: int = 0


@dataclass
class ActiveRun:
    chrom: str
    start: int
    end: int
    score_sum: float
    maximum: float

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class SelectedRegion:
    chrom: str
    start: int
    end: int


class RunCollector:
    def __init__(
        self,
        handle,
        *,
        min_run_length: int,
        max_run_length: int | None,
    ) -> None:
        self.handle = handle
        self.min_run_length = int(min_run_length)
        self.max_run_length = max_run_length
        self.counts: Counter[int] = Counter()
        self.summary = ScanSummary()

    def add(self, run: ActiveRun) -> None:
        length = run.length
        self.summary.total_runs_observed += 1
        if length < self.min_run_length:
            return
        if self.max_run_length is not None and length > self.max_run_length:
            return
        self.summary.total_runs_retained += 1
        self.summary.retained_positive_bases += length
        self.counts[length] += 1
        mean_score = run.score_sum / length
        self.handle.write(
            f"{run.chrom}\t{run.start}\t{run.end}\t{length}\t"
            f"{run.maximum:.10g}\t{mean_score:.10g}\t{run.score_sum:.10g}\n"
        )


class PositiveRunScanner:
    """Stateful chunk consumer that preserves runs across chunk boundaries."""

    def __init__(self, collector: RunCollector, threshold: float) -> None:
        self.collector = collector
        self.threshold = float(threshold)
        self.active: ActiveRun | None = None

    def _finish_active(self) -> None:
        if self.active is not None:
            self.collector.add(self.active)
            self.active = None

    def finish_region(self) -> None:
        self._finish_active()

    def consume(self, chrom: str, offset: int, values: np.ndarray) -> None:
        array = np.asarray(values, dtype=float)
        n = int(array.size)
        if n == 0:
            return

        finite = np.isfinite(array)
        positive = finite & (array > self.threshold)
        summary = self.collector.summary
        summary.scanned_bases += n
        positive_count = int(np.count_nonzero(positive))
        missing_count = int(np.count_nonzero(~finite))
        summary.positive_bases += positive_count
        summary.missing_bases += missing_count
        summary.nonpositive_bases += n - positive_count - missing_count

        if not positive.any():
            self._finish_active()
            return

        padded = np.concatenate(([False], positive, [False]))
        changes = np.diff(padded.astype(np.int8))
        starts = np.flatnonzero(changes == 1)
        ends = np.flatnonzero(changes == -1)

        for index, (local_start, local_end) in enumerate(zip(starts, ends)):
            local_start = int(local_start)
            local_end = int(local_end)
            genome_start = int(offset + local_start)
            genome_end = int(offset + local_end)
            segment = array[local_start:local_end]

            continues_previous = (
                index == 0
                and local_start == 0
                and self.active is not None
                and self.active.chrom == chrom
                and self.active.end == genome_start
            )
            if continues_previous:
                self.active.end = genome_end
                self.active.score_sum += float(np.sum(segment, dtype=np.float64))
                self.active.maximum = max(self.active.maximum, float(np.max(segment)))
            else:
                if self.active is not None:
                    self._finish_active()
                self.active = ActiveRun(
                    chrom=chrom,
                    start=genome_start,
                    end=genome_end,
                    score_sum=float(np.sum(segment, dtype=np.float64)),
                    maximum=float(np.max(segment)),
                )

            if local_end < n:
                self._finish_active()


def _natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", value)
    )


def _split_tokens(tokens: Sequence[str] | None) -> list[str]:
    if not tokens:
        return ["all"]
    output: list[str] = []
    for token in tokens:
        output.extend(piece.strip() for piece in token.split(",") if piece.strip())
    return output or ["all"]


def select_regions(
    tokens: Sequence[str] | None,
    chrom_sizes: Mapping[str, int],
) -> list[SelectedRegion]:
    available = list(chrom_sizes)
    available_set = set(available)
    requested: list[SelectedRegion] = []

    for token in _split_tokens(tokens):
        lower = token.lower()
        if lower == "all":
            requested.extend(SelectedRegion(chrom, 0, int(chrom_sizes[chrom])) for chrom in available)
            continue
        if lower == "autosomes":
            autosomes = [
                chrom
                for chrom in available
                if re.fullmatch(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2])", chrom, flags=re.IGNORECASE)
            ]
            requested.extend(SelectedRegion(chrom, 0, int(chrom_sizes[chrom])) for chrom in autosomes)
            continue

        if ":" in token:
            chrom, coordinates = token.split(":", 1)
            try:
                chrom = resolve_contig_name(chrom, available, source_label="BigWig")
            except KeyError as exc:
                raise ValueError(str(exc)) from exc
            match = re.fullmatch(r"(\d+)-(\d+)", coordinates)
            if match is None:
                raise ValueError(
                    f"Invalid genomic interval {token!r}; expected, for example, chr20:1000-2000"
                )
            start, end = int(match.group(1)), int(match.group(2))
            if start < 0 or end <= start or end > int(chrom_sizes[chrom]):
                raise ValueError(
                    f"Invalid interval {token!r} for {chrom} length {chrom_sizes[chrom]}"
                )
            requested.append(SelectedRegion(chrom, start, end))
            continue

        match = re.fullmatch(r"(chr)?(\d+)-(?:chr)?(\d+)", token, flags=re.IGNORECASE)
        if match:
            prefix = "chr" if match.group(1) else ""
            first, last = int(match.group(2)), int(match.group(3))
            step = 1 if last >= first else -1
            for number in range(first, last + step, step):
                requested_name = f"{prefix}{number}"
                try:
                    chrom = resolve_contig_name(requested_name, available, source_label="BigWig")
                except KeyError as exc:
                    raise ValueError(str(exc)) from exc
                requested.append(SelectedRegion(chrom, 0, int(chrom_sizes[chrom])))
            continue

        try:
            chrom = resolve_contig_name(token, available, source_label="BigWig")
        except KeyError as exc:
            raise ValueError(str(exc)) from exc
        requested.append(SelectedRegion(chrom, 0, int(chrom_sizes[chrom])))

    seen: set[tuple[str, int, int]] = set()
    unique: list[SelectedRegion] = []
    rank = {name: index for index, name in enumerate(available)}
    for region in requested:
        key = (region.chrom, region.start, region.end)
        if key not in seen:
            seen.add(key)
            unique.append(region)
    unique.sort(key=lambda region: (rank.get(region.chrom, len(rank)), region.start, region.end))
    merged: list[SelectedRegion] = []
    for region in unique:
        if (
            merged
            and merged[-1].chrom == region.chrom
            and region.start <= merged[-1].end
        ):
            previous = merged[-1]
            merged[-1] = SelectedRegion(
                previous.chrom,
                previous.start,
                max(previous.end, region.end),
            )
        else:
            merged.append(region)
    return merged


def weighted_quantile(counts: Counter[int], fraction: float) -> float:
    if not counts:
        return math.nan
    total = sum(counts.values())
    target = fraction * (total - 1)
    lower_rank = int(math.floor(target))
    upper_rank = int(math.ceil(target))

    def value_at(rank: int) -> int:
        cumulative = 0
        for value in sorted(counts):
            cumulative += counts[value]
            if cumulative > rank:
                return value
        return max(counts)

    lower = value_at(lower_rank)
    upper = value_at(upper_rank)
    if lower_rank == upper_rank:
        return float(lower)
    weight = target - lower_rank
    return float(lower + (upper - lower) * weight)


def write_length_counts(
    path: Path,
    counts: Counter[int],
) -> None:
    total = sum(counts.values())
    with path.open("wt", encoding="utf-8") as handle:
        handle.write("run_length_bp\tcount\tfraction\tpercent\n")
        for length in sorted(counts):
            count = counts[length]
            fraction = count / total if total else math.nan
            percent = fraction * 100.0 if total else math.nan
            handle.write(f"{length}\t{count}\t{fraction:.12g}\t{percent:.12g}\n")


def write_summary(
    path: Path,
    *,
    bigwig: Path,
    regions: Sequence[SelectedRegion],
    threshold: float,
    chunk_size: int,
    min_run_length: int,
    max_run_length: int | None,
    counts: Counter[int],
    scan: ScanSummary,
) -> None:
    total = sum(counts.values())
    mean = (
        sum(length * count for length, count in counts.items()) / total
        if total
        else math.nan
    )
    mode = min(
        (length for length, count in counts.items() if count == max(counts.values())),
        default=math.nan,
    )
    rows = {
        "bigwig": str(bigwig.resolve()),
        "selected_regions": len(regions),
        "threshold": threshold,
        "chunk_size_bp": chunk_size,
        "minimum_retained_run_length_bp": min_run_length,
        "maximum_retained_run_length_bp": max_run_length if max_run_length is not None else "none",
        "scanned_bases": scan.scanned_bases,
        "positive_bases": scan.positive_bases,
        "nonpositive_bases": scan.nonpositive_bases,
        "missing_bases": scan.missing_bases,
        "blacklisted_bases": scan.blacklisted_bases,
        "total_runs_observed": scan.total_runs_observed,
        "total_runs_retained": scan.total_runs_retained,
        "retained_positive_bases": scan.retained_positive_bases,
        "minimum_run_length_bp": min(counts) if counts else math.nan,
        "q1_run_length_bp": weighted_quantile(counts, 0.25),
        "median_run_length_bp": weighted_quantile(counts, 0.5),
        "mean_run_length_bp": mean,
        "q3_run_length_bp": weighted_quantile(counts, 0.75),
        "maximum_run_length_bp": max(counts) if counts else math.nan,
        "mode_run_length_bp": mode,
    }
    with path.open("wt", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in rows.items():
            handle.write(f"{key}\t{value}\n")


def plot_distribution(
    path: Path,
    counts: Counter[int],
    *,
    normalization: str,
    plot_x_max: int | None,
    title: str,
) -> Path:
    lengths = np.asarray(sorted(counts), dtype=int)
    raw_counts = np.asarray([counts[int(length)] for length in lengths], dtype=float)
    if normalization == "count":
        values = raw_counts
        ylabel = "Count"
    elif normalization == "fraction":
        values = raw_counts / raw_counts.sum() if raw_counts.size and raw_counts.sum() else raw_counts
        ylabel = "Fraction of runs"
    else:
        values = (
            raw_counts / raw_counts.sum() * 100.0
            if raw_counts.size and raw_counts.sum()
            else raw_counts
        )
        ylabel = "Runs (%)"

    figure, axis = plt.subplots(figsize=(10, 6))
    if lengths.size:
        axis.plot(lengths, values, linewidth=1.4, marker="o", markersize=2.0, markeredgewidth=0)
    from nucleosuite.plotting import apply_base_pair_x_axis, apply_integer_y_axis
    apply_base_pair_x_axis(axis, lengths)
    if normalization == "count":
        apply_integer_y_axis(axis)
    axis.set_xlabel("Contiguous positive run length (bp)")
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.grid(True, alpha=0.25)
    if lengths.size:
        data_max = max(1, int(lengths.max()))
        axis.set_xlim(0, min(data_max, plot_x_max) if plot_x_max is not None else data_max)
    elif plot_x_max is not None:
        axis.set_xlim(0, plot_x_max)
    figure.tight_layout()
    from nucleosuite.plotting import save_figure
    saved = save_figure(figure, path, default_dpi=300)
    plt.close(figure)
    return saved


def run_analysis(
    bigwig_path: str | Path,
    output_prefix: str | Path,
    *,
    contigs: Sequence[str] | None = None,
    threshold: float = 0.0,
    chunk_size: int = 1_000_000,
    min_run_length: int = 1,
    max_run_length: int | None = None,
    plot_x_max: int | None = 550,
    normalization: str = "count",
    title: str | None = None,
    blacklist_bed: str | Path | None = None,
    progress: ProgressReporter | None = None,
) -> dict[str, Path]:
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("pyBigWig is required for positive-runs") from exc

    bigwig_path = Path(bigwig_path)
    if not bigwig_path.is_file():
        raise FileNotFoundError(bigwig_path)
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = {
        "runs": Path(f"{prefix}_runs.tsv.gz"),
        "counts": Path(f"{prefix}_length_counts.tsv"),
        "summary": Path(f"{prefix}_summary.tsv"),
        "plot": __import__("nucleosuite.plotting", fromlist=["plot_path"]).plot_path(f"{prefix}_run_length_distribution.png"),
    }

    bigwig = pyBigWig.open(str(bigwig_path))
    if bigwig is None:
        raise OSError(f"Unable to open BigWig: {bigwig_path}")
    try:
        chrom_sizes = {name: int(length) for name, length in bigwig.chroms().items()}
        from nucleosuite.core.blacklist import load_blacklist
        blacklist = load_blacklist(
            blacklist_bed, list(chrom_sizes), list(chrom_sizes.values())
        )
        regions = select_regions(contigs, chrom_sizes)
        with gzip.open(outputs["runs"], "wt", encoding="utf-8") as run_handle:
            run_handle.write(
                "chrom\tstart\tend\trun_length_bp\tmaximum_score\tmean_score\tscore_area\n"
            )
            collector = RunCollector(
                run_handle,
                min_run_length=min_run_length,
                max_run_length=max_run_length,
            )
            scanner = PositiveRunScanner(collector, threshold)
            for region_index, region in enumerate(regions, start=1):
                if progress is not None:
                    progress.contig(
                        "Scanning signal",
                        region.chrom,
                        region_index,
                        len(regions),
                    )
                scanner.finish_region()
                for start in range(region.start, region.end, chunk_size):
                    end = min(start + chunk_size, region.end)
                    values = np.asarray(
                        bigwig.values(region.chrom, start, end, numpy=True),
                        dtype=float,
                    )
                    if blacklist is not None:
                        collector.summary.blacklisted_bases += blacklist.mask_values(
                            region.chrom, start, values
                        )
                    scanner.consume(region.chrom, start, values)
                scanner.finish_region()
    finally:
        close = getattr(bigwig, "close", None)
        if callable(close):
            close()

    write_length_counts(outputs["counts"], collector.counts)
    write_summary(
        outputs["summary"],
        bigwig=bigwig_path,
        regions=regions,
        threshold=threshold,
        chunk_size=chunk_size,
        min_run_length=min_run_length,
        max_run_length=max_run_length,
        counts=collector.counts,
        scan=collector.summary,
    )
    saved_plot = plot_distribution(
        outputs["plot"],
        collector.counts,
        normalization=normalization,
        plot_x_max=plot_x_max,
        title=title or f"Positive run lengths: {bigwig_path.stem}",
    )
    if saved_plot is not None:
        outputs["plot"] = Path(saved_plot)
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite positive-runs",
        description=(
            "Count uninterrupted per-base runs above a threshold in one BigWig track. "
            "Zero, negative, missing and non-finite values terminate a run."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--bigwig", required=True, help="Input BigWig signal track.")
    parser.add_argument(
        "--blacklist-bed",
        help="BED blacklist; masked bases terminate positive runs and count as missing.",
    )
    parser.add_argument(
        "--output-prefix",
        help="Output prefix. Default: the input BigWig basename in the current directory.",
    )
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=["all"],
        help=(
            "Contigs or genomic ranges to scan. Supports comma lists, chr1-22, autosomes, "
            "all and chr20:START-END."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="A base belongs to a run only when its value is strictly greater than this threshold.",
    )
    parser.add_argument("--chunk-size", type=int, default=1_000_000, help="BigWig read chunk size in bp.")
    parser.add_argument("--min-run-length", type=int, default=1, help="Minimum retained run length in bp.")
    parser.add_argument(
        "--max-run-length",
        type=int,
        default=0,
        help="Maximum retained run length in bp; 0 disables the maximum.",
    )
    parser.add_argument(
        "--plot-x-max",
        type=int,
        default=550,
        help="Displayed x-axis maximum in bp; 0 uses the full observed range.",
    )
    parser.add_argument(
        "--normalization",
        choices=["count", "fraction", "percent"],
        default="count",
        help="Y-axis measure for the distribution plot.",
    )
    parser.add_argument("--title", help="Optional plot title.")
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser)
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def run(args: argparse.Namespace) -> int:
    parser = build_parser()
    if not math.isfinite(args.threshold):
        parser.error("--threshold must be finite")
    if args.chunk_size < 1:
        parser.error("--chunk-size must be at least 1")
    if args.min_run_length < 1:
        parser.error("--min-run-length must be at least 1")
    if args.max_run_length < 0:
        parser.error("--max-run-length must be at least 0")
    if args.max_run_length and args.max_run_length < args.min_run_length:
        parser.error("--max-run-length cannot be smaller than --min-run-length")
    if args.plot_x_max < 0:
        parser.error("--plot-x-max must be at least 0")

    from nucleosuite.parallel import run_bigwig_per_contig

    def serial(namespace: argparse.Namespace) -> int:
        reporter = ProgressReporter("positive-runs")
        reporter.stage(f"Opening signal track: {namespace.bigwig}")
        output_prefix = namespace.output_prefix or Path(namespace.bigwig).name.removesuffix(".bw").removesuffix(".bigWig")
        outputs = run_analysis(
            namespace.bigwig,
            output_prefix,
            contigs=namespace.contigs,
            threshold=namespace.threshold,
            chunk_size=namespace.chunk_size,
            min_run_length=namespace.min_run_length,
            max_run_length=namespace.max_run_length or None,
            plot_x_max=namespace.plot_x_max or None,
            normalization=namespace.normalization,
            title=namespace.title,
            blacklist_bed=namespace.blacklist_bed,
            progress=reporter,
        )
        reporter.stage("Writing run-length tables and plot")
        for name, path in outputs.items():
            print(f"{name}\t{path}")
        return 0

    return run_bigwig_per_contig(
        "positive-runs", args, serial,
        bigwig_attr="bigwig", selector_attr="contigs", prefix_attr="output_prefix"
    )


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
