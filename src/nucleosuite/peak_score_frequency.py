#!/usr/bin/env python3
"""Plot score-frequency distributions from one or more BED-like peak files."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from nucleosuite.io.intervals import open_interval_text
from nucleosuite.parallel import add_parallel_arguments
from nucleosuite.partitioned import run_partitioned_command
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.progress import ProgressReporter


@dataclass(frozen=True)
class PeakScoreSet:
    label: str
    path: Path
    scores: np.ndarray
    records: list[tuple[str, int, int, str, float]]
    blacklisted_records: int = 0


def _parse_labelled_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        path = Path(value)
        label = path.name
        for suffix in (".bed.gz", ".bigbed", ".bed", ".bb", ".gz"):
            if label.lower().endswith(suffix):
                label = label[: -len(suffix)]
                break
        return label, path
    label, raw_path = value.split("=", 1)
    label = label.strip()
    if not label:
        raise argparse.ArgumentTypeError("Dataset labels cannot be empty")
    return label, Path(raw_path)


def _read_scores(
    label: str,
    path: Path,
    score_column: int,
    blacklist: BlacklistIndex | None = None,
    progress: ProgressReporter | None = None,
) -> PeakScoreSet:
    if not path.exists():
        raise FileNotFoundError(path)
    index = score_column - 1
    records: list[tuple[str, int, int, str, float]] = []
    blacklisted_records = 0
    seen_contigs: set[str] = set()
    if progress is not None:
        progress.file_start(f"peaks ({label})", path)
    with open_interval_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least BED3")
            if index >= len(fields):
                raise ValueError(
                    f"{path}:{line_number}: score column {score_column} is absent "
                    f"from a {len(fields)}-column record"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
                score = float(fields[index])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid coordinate or score") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: require 0 <= start < end")
            if progress is not None and fields[0] not in seen_contigs:
                seen_contigs.add(fields[0])
                progress.reading_contig(f"peaks ({label})", fields[0])
            if blacklist is not None and blacklist.overlaps(fields[0], start, end):
                blacklisted_records += 1
                continue
            if not math.isfinite(score):
                continue
            name = fields[3] if len(fields) >= 4 else f"{fields[0]}:{start}-{end}"
            records.append((fields[0], start, end, name, score))
    if not records:
        raise ValueError(f"No finite scores were found in {path}")
    return PeakScoreSet(
        label=label,
        path=path,
        scores=np.asarray([record[4] for record in records], dtype=float),
        records=records,
        blacklisted_records=blacklisted_records,
    )


def _histogram_edges(
    datasets: Sequence[PeakScoreSet],
    *,
    bins: int | None,
    bin_width: float | None,
    score_min: float | None,
    score_max: float | None,
    integer_bins: bool,
    score_scale: float = 1.0,
) -> tuple[np.ndarray, np.ndarray | None]:
    all_scores = np.concatenate([dataset.scores for dataset in datasets]) * float(score_scale)
    if integer_bins:
        if bins is not None or bin_width is not None:
            raise ValueError("--integer-bins cannot be combined with --bins or --bin-width")
        rounded = np.floor(all_scores + 0.5).astype(np.int64)
        observed_lower = int(np.min(rounded))
        observed_upper = int(np.max(rounded))
        lower = observed_lower if score_min is None else int(math.floor(float(score_min) + 0.5))
        upper = observed_upper if score_max is None else int(math.floor(float(score_max) + 0.5))
        # Peak scores are normally non-negative. Keep zero in the shared range,
        # while still retaining unexpected negative scores rather than dropping them.
        lower = min(0, lower)
        upper = max(0, upper)
        if upper < lower:
            raise ValueError("--score-max must be greater than or equal to --score-min")
        labels = np.arange(lower, upper + 1, dtype=np.int64)
        edges = np.arange(lower - 0.5, upper + 1.5, 1.0, dtype=float)
        return edges, labels

    lower = float(np.min(all_scores)) if score_min is None else float(score_min)
    upper = float(np.max(all_scores)) if score_max is None else float(score_max)
    if upper < lower:
        raise ValueError("--score-max must be greater than or equal to --score-min")
    if upper == lower:
        padding = max(abs(lower) * 0.01, 0.5)
        lower -= padding
        upper += padding
    if bin_width is not None:
        if bin_width <= 0:
            raise ValueError("--bin-width must be positive")
        first = math.floor(lower / bin_width) * bin_width
        last = math.ceil(upper / bin_width) * bin_width
        if last <= first:
            last = first + bin_width
        edge_count = int(round((last - first) / bin_width)) + 1
        if edge_count > 100_001:
            raise ValueError("Requested bin width would create more than 100,000 bins")
        return first + np.arange(edge_count, dtype=float) * bin_width, None
    if bins is None:
        raise ValueError("A continuous histogram requires --bins or --bin-width")
    if bins < 1:
        raise ValueError("--bins must be positive")
    return np.linspace(lower, upper, bins + 1, dtype=float), None


def _normalised_values(counts: np.ndarray, edges: np.ndarray, mode: str) -> np.ndarray:
    total = float(counts.sum())
    if mode == "count":
        return counts.astype(float)
    if total <= 0:
        return np.zeros_like(counts, dtype=float)
    if mode == "fraction":
        return counts / total
    if mode == "percent":
        return counts / total * 100.0
    if mode == "density":
        widths = np.diff(edges)
        return counts / total / widths
    raise ValueError(f"Unknown normalization: {mode}")


def _write_values(datasets: Sequence[PeakScoreSet], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8") as handle:
        handle.write("dataset\tsource\tchrom\tstart\tend\tname\tscore\n")
        for dataset in datasets:
            for chrom, start, end, name, score in dataset.records:
                handle.write(
                    f"{dataset.label}\t{dataset.path}\t{chrom}\t{start}\t{end}\t{name}\t{score:.12g}\n"
                )


def _write_summary(datasets: Sequence[PeakScoreSet], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "dataset", "source", "count", "minimum", "q1", "median", "mean",
                "q3", "maximum", "standard_deviation", "blacklisted_records_excluded",
            ]
        )
        for dataset in datasets:
            scores = dataset.scores
            q1, median, q3 = np.quantile(scores, [0.25, 0.5, 0.75])
            writer.writerow(
                [
                    dataset.label,
                    dataset.path,
                    scores.size,
                    f"{np.min(scores):.12g}",
                    f"{q1:.12g}",
                    f"{median:.12g}",
                    f"{np.mean(scores):.12g}",
                    f"{q3:.12g}",
                    f"{np.max(scores):.12g}",
                    f"{np.std(scores, ddof=0):.12g}",
                    dataset.blacklisted_records,
                ]
            )


def _write_frequency(
    datasets: Sequence[PeakScoreSet],
    edges: np.ndarray,
    output: Path,
    *,
    integer_labels: np.ndarray | None = None,
    score_scale: float = 1.0,
) -> dict[str, np.ndarray]:
    output.parent.mkdir(parents=True, exist_ok=True)
    histograms: dict[str, np.ndarray] = {}
    with output.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if integer_labels is not None:
            writer.writerow(
                [
                    "dataset", "score", "count", "fraction", "percent",
                    "cumulative_count", "cumulative_fraction", "cumulative_percent",
                ]
            )
        else:
            writer.writerow(
                [
                    "dataset", "bin_index", "bin_start", "bin_end", "bin_midpoint",
                    "count", "fraction", "percent", "density",
                ]
            )
        widths = np.diff(edges)
        for dataset in datasets:
            scaled_scores = dataset.scores * float(score_scale)
            counts, _ = np.histogram(scaled_scores, bins=edges)
            histograms[dataset.label] = counts
            total = counts.sum()
            fractions = counts / total if total else np.zeros_like(counts, dtype=float)
            if integer_labels is not None:
                cumulative_counts = np.cumsum(counts)
                cumulative_fractions = np.cumsum(fractions)
                for score, count, fraction, cumulative_count, cumulative_fraction in zip(
                    integer_labels, counts, fractions, cumulative_counts, cumulative_fractions
                ):
                    writer.writerow(
                        [
                            dataset.label,
                            int(score),
                            int(count),
                            f"{fraction:.12g}",
                            f"{(fraction * 100):.12g}",
                            int(cumulative_count),
                            f"{cumulative_fraction:.12g}",
                            f"{(cumulative_fraction * 100):.12g}",
                        ]
                    )
            else:
                densities = fractions / widths
                for index, count in enumerate(counts):
                    writer.writerow(
                        [
                            dataset.label,
                            index + 1,
                            f"{edges[index]:.12g}",
                            f"{edges[index + 1]:.12g}",
                            f"{((edges[index] + edges[index + 1]) / 2):.12g}",
                            int(count),
                            f"{fractions[index]:.12g}",
                            f"{(fractions[index] * 100):.12g}",
                            f"{densities[index]:.12g}",
                        ]
                    )
    return histograms


def _plot_frequency(
    datasets: Sequence[PeakScoreSet],
    edges: np.ndarray,
    histograms: dict[str, np.ndarray],
    output: Path,
    *,
    normalization: str,
    title: str | None,
    log_y: bool,
    plot_x_min: float | None,
    plot_x_max: float | None,
    integer_labels: np.ndarray | None = None,
    score_scale: float = 1.0,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    fig, ax = plt.subplots(figsize=(10, 6))
    midpoints = integer_labels.astype(float) if integer_labels is not None else (edges[:-1] + edges[1:]) / 2
    for dataset in datasets:
        values = _normalised_values(histograms[dataset.label], edges, normalization)
        ax.step(midpoints, values, where="mid", linewidth=1.5, label=dataset.label)
    ax.set_xlabel(
        "Peak score" if math.isclose(score_scale, 1.0)
        else f"Peak score (×{score_scale:g})"
    )
    ylabel = {
        "count": "Frequency (count)",
        "fraction": "Frequency (fraction)",
        "percent": "Frequency (%)",
        "density": "Probability density",
    }[normalization]
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if len(datasets) > 1:
        ax.legend(frameon=False)
    if log_y:
        ax.set_yscale("log")
    if plot_x_min is not None or plot_x_max is not None:
        left, right = ax.get_xlim()
        ax.set_xlim(
            left if plot_x_min is None else plot_x_min,
            right if plot_x_max is None else plot_x_max,
        )
    if integer_labels is not None:
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    from nucleosuite.plotting import save_figure
    saved = save_figure(fig, output, default_dpi=220, bbox_inches="tight")
    plt.close(fig)
    return saved


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite peak-score-frequency",
        description=(
            "Plot binned score frequencies from one or more BED, BED.gz or bigBed "
            "peak files. Repeated labelled inputs share the same histogram bins."
        ),
    )
    parser.add_argument(
        "--peaks",
        action="append",
        required=True,
        metavar="LABEL=FILE",
        help="Peak interval input. Repeat to overlay observed and control callsets.",
    )
    parser.add_argument("--output-prefix", "-o", required=True, help="Path prefix for score tables, summaries, and plot output.")
    parser.add_argument(
        "--blacklist-bed",
        help="BED blacklist; complete overlapping peak records are excluded.",
    )
    parser.add_argument("--score-column", type=int, default=5, help="1-based score column. Default: 5")
    parser.add_argument(
        "--score-scale", type=float, default=1.0,
        help=(
            "Multiply scores by this value before histogram binning and plotting "
            "(default: 1). Raw score summaries/detail tables remain unscaled."
        ),
    )
    parser.add_argument(
        "--integer-bins",
        action="store_true",
        help=(
            "Round scores to the nearest integer and report every integer score from "
            "the shared selected minimum through the shared selected maximum. This is "
            "the default when neither --bins nor --bin-width is supplied."
        ),
    )
    parser.add_argument(
        "--bins",
        type=int,
        help="Use this many shared continuous histogram bins instead of integer-score bins.",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        help="Use fixed-width continuous bins instead of integer-score bins.",
    )
    parser.add_argument("--score-min", type=float, help="Lower histogram boundary after --score-scale. Default: observed scaled minimum.")
    parser.add_argument("--score-max", type=float, help="Upper histogram boundary after --score-scale. Default: observed scaled maximum.")
    parser.add_argument(
        "--normalization",
        choices=("count", "fraction", "percent", "density"),
        default="count",
        help="Y-axis and frequency-table measure (default: count).",
    )
    parser.add_argument("--plot-x-min", type=float, help="Displayed x-axis minimum; histogram tables retain the full selected range.")
    parser.add_argument("--plot-x-max", type=float, help="Displayed x-axis maximum; histogram tables retain the full selected range.")
    parser.add_argument("--log-y", action="store_true", help="Use a logarithmic y-axis for the frequency plot.")
    parser.add_argument("--title", help="Optional plot title.")
    parser.add_argument(
        "--write-detail-tables", action="store_true",
        help="Write the individual finite peak-score table; omitted by default.",
    )
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive score analysis "
            "(default: 1; independent of suite --cores)."
        ),
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def _run_serial(args: argparse.Namespace) -> int:
    if args.score_column < 1:
        raise ValueError("--score-column must be at least 1")
    if not math.isfinite(args.score_scale) or args.score_scale <= 0:
        raise ValueError("--score-scale must be a finite value greater than zero")
    labelled = [_parse_labelled_path(value) for value in args.peaks]
    labels = [label for label, _ in labelled]
    if len(set(labels)) != len(labels):
        raise ValueError("Every --peaks label must be unique")
    reporter = ProgressReporter("peak-score-frequency")
    blacklist = load_blacklist_unbounded(args.blacklist_bed)
    datasets = [
        _read_scores(label, path, args.score_column, blacklist, reporter)
        for label, path in labelled
    ]
    reporter.stage(
        f"Loaded {sum(dataset.scores.size for dataset in datasets):,} scores "
        f"from {len(datasets):,} callset(s); calculating shared bins"
    )
    integer_bins = bool(args.integer_bins or (args.bins is None and args.bin_width is None))
    edges, integer_labels = _histogram_edges(
        datasets,
        bins=args.bins,
        bin_width=args.bin_width,
        score_min=args.score_min,
        score_max=args.score_max,
        integer_bins=integer_bins,
        score_scale=args.score_scale,
    )
    from nucleosuite.output_naming import parameterized_prefix

    bin_parameter = (
        ("bins", args.bins)
        if args.bins is not None
        else ("binwidth", args.bin_width)
        if args.bin_width is not None
        else ("bins", "integer")
    )
    prefix = parameterized_prefix(
        args.output_prefix,
        (
            bin_parameter,
            ("scorescale", args.score_scale),
            ("scoremin", args.score_min),
            ("scoremax", args.score_max),
            ("norm", args.normalization),
        ),
    )
    values_path = Path(f"{prefix}_scores.tsv.gz")
    frequency_path = Path(f"{prefix}_score_frequency.tsv")
    summary_path = Path(f"{prefix}_score_summary.tsv")
    plot_path = Path(f"{prefix}_score_frequency.png")
    reporter.stage("Writing score tables and frequency plot")
    if args.write_detail_tables:
        _write_values(datasets, values_path)
    histograms = _write_frequency(
        datasets, edges, frequency_path, integer_labels=integer_labels,
        score_scale=args.score_scale,
    )
    _write_summary(datasets, summary_path)
    plot_path = _plot_frequency(
        datasets,
        edges,
        histograms,
        plot_path,
        normalization=args.normalization,
        title=args.title,
        log_y=args.log_y,
        plot_x_min=args.plot_x_min,
        plot_x_max=args.plot_x_max,
        integer_labels=integer_labels,
        score_scale=args.score_scale,
    )
    from nucleosuite.plotting import write_plot_metadata
    write_plot_metadata(
        plot_path,
        extra={
            "source_table": str(frequency_path),
            "detected_plot_type": "peak-score-frequency",
            "score_scale": args.score_scale,
            "x_label": (
                "Peak score" if math.isclose(args.score_scale, 1.0)
                else f"Peak score (×{args.score_scale:g})"
            ),
        },
    )
    if args.write_detail_tables:
        print(f"Wrote: {values_path}")
    for path in (frequency_path, summary_path, plot_path):
        print(f"Wrote: {path}")
    return 0


def run(args: argparse.Namespace) -> int:
    from nucleosuite.output_naming import parameterized_prefix

    bin_parameter = (
        ("bins", args.bins)
        if args.bins is not None
        else ("binwidth", args.bin_width)
        if args.bin_width is not None
        else ("bins", "integer")
    )
    args.output_prefix = str(
        parameterized_prefix(
            args.output_prefix,
            (
                bin_parameter,
                ("scorescale", args.score_scale),
                ("scoremin", args.score_min),
                ("scoremax", args.score_max),
                ("norm", args.normalization),
            ),
        )
    )
    return run_partitioned_command(
        "peak-score-frequency",
        args,
        _run_serial,
        runner_module="nucleosuite.peak_score_frequency",
        runner_function="_run_serial",
        primary_attr=None,
        primary_list_attr="peaks",
        output_prefix_attr="output_prefix",
        named_path_list_attrs=("peaks",),
        path_attrs=("blacklist_bed",),
    )


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
