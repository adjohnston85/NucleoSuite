"""Fragment statistics, console histograms and BED fragment writers."""

from __future__ import annotations

from collections import Counter
from typing import Iterable, Mapping, Optional, TextIO, Tuple
import sys

from nucleosuite.profile_plots import plot_count_profile


def print_fragment_length_histogram(
    length_counts: Mapping[int, int],
    *,
    label: str | None = None,
    width: int = 100,
    stream: TextIO | None = None,
) -> None:
    """Print a proportional ASCII fragment-length histogram.

    The most abundant fragment length receives ``width`` dashes. Other bars are
    scaled relative to that maximum and rounded to the nearest dash. A positive
    count always receives at least one dash.
    """
    if width < 1:
        raise ValueError("Histogram width must be at least 1")
    stream = stream or sys.stdout

    title = "Fragment-length distribution"
    if label:
        title += f": {label}"
    print(f"\n{title}", file=stream)
    print(
        f"{'length':>6}  {'distribution':<{width}}  count",
        file=stream,
    )

    positive = {
        int(length): int(count)
        for length, count in length_counts.items()
        if int(count) > 0
    }
    if not positive:
        print("No fragments were counted.", file=stream)
        print(file=stream)
        return

    max_count = max(positive.values())
    for length in sorted(positive):
        count = positive[length]
        bar_length = max(1, int(round((count / max_count) * width)))
        print(
            f"{length:>6}  {'-' * bar_length:<{width}}  {count}",
            file=stream,
        )
    print(file=stream)


def write_fragment_outputs(
    output_prefix: str,
    total_fragments_filtered: int,
    total_fragments_used: int,
    unique_bases_covered: int,
    length_counts: Counter,
    dedup_scope: str,
    max_duplicates: int | None = None,
    max_per_coordinate: int = 0,
    *,
    print_histogram: bool = True,
    histogram_width: int = 100,
) -> tuple[str, str]:
    """Write fragment summary/count TSVs and optionally print the distribution.

    ``max_duplicates`` is the input-fragment limit. When it is omitted,
    ``max_per_coordinate`` supplies that limit and the output-coordinate cap is
    disabled.
    """
    if max_duplicates is None:
        max_duplicates = int(max_per_coordinate)
        max_per_coordinate = 0
    summary_path = f"{output_prefix}_fragment_summary.tsv"
    lengths_path = f"{output_prefix}_fragment_length_counts.tsv"
    with open(summary_path, "w", encoding="utf-8") as output:
        output.write("metric\tvalue\n")
        output.write(f"total_fragments_filtered_all\t{total_fragments_filtered}\n")
        output.write(f"total_fragments_used_in_range\t{total_fragments_used}\n")
        output.write(
            f"unique_bases_covered_by_used_fragments\t{unique_bases_covered}\n"
        )
        output.write(f"dedup_scope\t{dedup_scope}\n")
        output.write(f"max_duplicates\t{max_duplicates}\n")
        output.write(f"max_per_coordinate\t{max_per_coordinate}\n")

    with open(lengths_path, "w", encoding="utf-8") as output:
        output.write("fragment_length\tcount\n")
        for length in sorted(length_counts):
            output.write(f"{int(length)}\t{int(length_counts[length])}\n")

    plot_count_profile(
        lengths_path,
        f"{output_prefix}_fragment_length_distribution.png",
        x_column="fragment_length",
        y_column="count",
        xlabel="Fragment length (bp)",
        ylabel="Fragment count",
        title="Fragment-length distribution",
    )

    if print_histogram:
        print_fragment_length_histogram(
            length_counts,
            label=output_prefix,
            width=histogram_width,
        )
    return summary_path, lengths_path


def write_fragment_bed_rows(
    output_path: str,
    contig: str,
    records: Iterable[Tuple[int, int, Optional[str]]],
    mode: str = "a",
    include_type: bool = False,
) -> None:
    """Write BED3 or BED4 fragment records."""
    with open(output_path, mode, encoding="utf-8") as output:
        for start, end, ww_type in records:
            if include_type:
                output.write(f"{contig}\t{start}\t{end}\t{ww_type or 'unclassified'}\n")
            else:
                output.write(f"{contig}\t{start}\t{end}\n")
