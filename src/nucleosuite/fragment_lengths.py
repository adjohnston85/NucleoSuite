#!/usr/bin/env python3
"""Count paired-end fragment lengths, optionally stratified by BED intervals.

Without ``--bed``, one genome-wide fragment-length distribution is produced.
With ``--bed``, each fragment is assigned by its midpoint to the BED interval
containing that midpoint and distributions are reported by BED label.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence, TextIO

import numpy as np

from nucleosuite.io.summaries import print_fragment_length_histogram
from nucleosuite.io import open_text as open_interval_text
from nucleosuite.core.fragment_inputs import IntervalFragmentSource
from nucleosuite.core.blacklist import load_blacklist
from nucleosuite.progress import ProgressReporter
from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.nrl import Peak, Regression

try:
    import pysam
except ImportError:  # allow command help without runtime dependencies installed
    pysam = None


@dataclass(frozen=True)
class BedInterval:
    start: int
    end: int
    label: str


@dataclass
class IntervalIndex:
    intervals: list[BedInterval]
    starts: list[int]
    prefix_max_ends: list[int]


@dataclass
class CountSummary:
    reads_examined: int = 0
    fragments_counted: int = 0
    fragments_unassigned: int = 0
    fragments_multiassigned: int = 0


@dataclass(frozen=True)
class FragmentSizeNRLResult:
    """Peak calls and regression for one fragment-length distribution."""

    label: str
    fragment_lengths: np.ndarray
    counts: np.ndarray
    density: np.ndarray
    local_values: np.ndarray
    detection_values: np.ndarray
    peaks: tuple[Peak, ...]
    regression: Regression
    minimum: int
    maximum: int
    peak_resolution: float
    detection_window: int
    local_max_window: int


def observed_maximum_length(counts: Mapping[str, Counter[int]]) -> int | None:
    """Return the longest positive-count fragment across all profiles."""

    observed = [
        int(length)
        for profile in counts.values()
        for length, count in profile.items()
        if count > 0
    ]
    return max(observed) if observed else None


def effective_plot_maximum(
    counts: Mapping[str, Counter[int]], requested_maximum: int = 1000
) -> int:
    """Cap a plot at both its requested maximum and longest counted fragment."""

    observed = observed_maximum_length(counts)
    if observed is None:
        raise ValueError("No fragment lengths were counted")
    return min(int(requested_maximum), observed)


def open_text(path: str | Path) -> TextIO:
    """Open BED text, BED.gz, or bigBed input."""
    return open_interval_text(path)


def sanitize_filename(text: str) -> str:
    """Return a filesystem-safe representation of a BED label."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    return cleaned.strip("._") or "unnamed"


def parse_bed(
    bed_path: str | Path,
    label_column: int = 4,
    default_label: str = "regions",
) -> dict[str, IntervalIndex]:
    """Read a BED/ BED.gz file and build a midpoint-query index per contig.

    Parameters
    ----------
    bed_path
        BED file with at least three columns.
    label_column
        One-based label column. Set to 0 to assign all intervals the same
        ``default_label``.
    default_label
        Label used when ``label_column`` is 0.
    """
    if label_column < 0:
        raise ValueError("label_column must be 0 or a positive one-based column")

    by_contig: dict[str, list[BedInterval]] = defaultdict(list)

    with open_text(bed_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{bed_path}:{line_number}: expected at least 3 columns")

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{bed_path}:{line_number}: BED start/end must be integers"
                ) from exc

            if start < 0 or end <= start:
                raise ValueError(
                    f"{bed_path}:{line_number}: invalid half-open interval {start}-{end}"
                )

            if label_column == 0:
                label = default_label
            else:
                idx = label_column - 1
                if idx >= len(fields):
                    raise ValueError(
                        f"{bed_path}:{line_number}: label column {label_column} is absent"
                    )
                label = fields[idx]

            by_contig[fields[0]].append(BedInterval(start, end, label))

    indexes: dict[str, IntervalIndex] = {}
    for contig, intervals in by_contig.items():
        intervals.sort(key=lambda iv: (iv.start, iv.end, iv.label))
        starts = [iv.start for iv in intervals]
        prefix_max_ends: list[int] = []
        running_max = -1
        for interval in intervals:
            running_max = max(running_max, interval.end)
            prefix_max_ends.append(running_max)
        indexes[contig] = IntervalIndex(intervals, starts, prefix_max_ends)

    if not indexes:
        raise ValueError(f"No valid intervals were found in BED file: {bed_path}")
    return indexes


def labels_at_position(index: IntervalIndex, position: int) -> list[str]:
    """Return labels of all half-open intervals containing ``position``."""
    right = bisect.bisect_right(index.starts, position) - 1
    if right < 0:
        return []

    labels: list[str] = []
    i = right
    while i >= 0 and index.prefix_max_ends[i] > position:
        interval = index.intervals[i]
        if interval.start <= position < interval.end:
            labels.append(interval.label)
        i -= 1
    labels.reverse()
    return labels


def parse_contigs(value: str | Sequence[str] | None) -> set[str] | None:
    """Parse one or more exact contig selectors; ``all`` means all input contigs."""
    if value is None:
        return None
    values = [value] if isinstance(value, str) else list(value)
    contigs: set[str] = set()
    for token in values:
        for item in str(token).split(","):
            item = item.strip()
            if not item:
                continue
            if item.casefold() == "all":
                return None
            contigs.add(item)
    if not contigs:
        raise ValueError("--contigs did not contain any contig names")
    return contigs


def fragment_from_read1(read: pysam.AlignedSegment) -> tuple[int, int, int] | None:
    """Return ``(start, end, length)`` for a same-contig paired fragment."""
    if read.template_length == 0 or read.next_reference_id != read.reference_id:
        return None

    length = abs(read.template_length)
    start = min(read.reference_start, read.next_reference_start)
    end = start + length
    if start < 0 or end <= start:
        return None
    return start, end, length


def read_passes_filters(
    read: pysam.AlignedSegment,
    *,
    min_mapq: int,
    include_duplicates: bool,
    require_proper_pair: bool,
) -> bool:
    """Apply common fragment-level read filters."""
    if not read.is_paired or not read.is_read1:
        return False
    if read.is_unmapped or read.mate_is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    if read.mapping_quality < min_mapq:
        return False
    if read.is_duplicate and not include_duplicates:
        return False
    if require_proper_pair and not read.is_proper_pair:
        return False
    return True


def count_fragment_lengths(
    bam_paths: str | Path | Sequence[str | Path],
    *,
    bed_index: Mapping[str, IntervalIndex] | None = None,
    contigs: set[str] | None = None,
    min_mapq: int = 0,
    include_duplicates: bool = False,
    require_proper_pair: bool = True,
    min_length: int = 1,
    max_length: int | None = None,
    overlap_policy: str = "all",
    include_unassigned: bool = False,
    unassigned_label: str = "unassigned",
    blacklist_bed: str | None = None,
    progress: ProgressReporter | None = None,
) -> tuple[dict[str, Counter[int]], CountSummary]:
    """Count fragment lengths genome-wide or by BED label.

    ``overlap_policy`` controls fragments whose midpoint lies in overlapping BED
    intervals: ``all`` assigns to every distinct label, ``first`` assigns to the
    first interval in genomic sort order, and ``error`` raises an exception.
    """
    if min_length < 1:
        raise ValueError("min_length must be at least 1")
    if max_length is not None and max_length < min_length:
        raise ValueError("max_length must be greater than or equal to min_length")
    if overlap_policy not in {"all", "first", "error"}:
        raise ValueError("overlap_policy must be one of: all, first, error")

    if pysam is None:
        raise RuntimeError("pysam is required for fragment-length counting")

    counts: dict[str, Counter[int]] = defaultdict(Counter)
    summary = CountSummary()
    seen_contigs: set[str] = set()

    if isinstance(bam_paths, (str, Path)):
        paths = [Path(bam_paths)]
    else:
        paths = [Path(path) for path in bam_paths]
    if not paths:
        raise ValueError("At least one BAM file is required")

    header_handles = [pysam.AlignmentFile(str(path), "rb") for path in paths]
    try:
        merged_header = merge_bam_reference_headers_with_aliases(header_handles)
    finally:
        for handle in header_handles:
            close = getattr(handle, "close", None)
            if close is not None:
                close()
    blacklist = load_blacklist(
        blacklist_bed, merged_header.references, merged_header.lengths
    )
    if contigs is not None:
        resolved_contigs: set[str] = set()
        for requested in contigs:
            try:
                resolved_contigs.add(
                    resolve_contig_name(
                        requested,
                        merged_header.references,
                        source_label="BAM headers",
                    )
                )
            except KeyError as exc:
                raise ValueError(str(exc)) from exc
        contigs = resolved_contigs

    for path, canonical_to_source in zip(paths, merged_header.source_contigs):
        source_to_canonical = {source: canonical for canonical, source in canonical_to_source.items()}
        with pysam.AlignmentFile(str(path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                summary.reads_examined += 1
                if not read_passes_filters(
                    read,
                    min_mapq=min_mapq,
                    include_duplicates=include_duplicates,
                    require_proper_pair=require_proper_pair,
                ):
                    continue

                contig = source_to_canonical.get(read.reference_name, read.reference_name)
                if contigs is not None and contig not in contigs:
                    continue
                if progress is not None and contig not in seen_contigs:
                    seen_contigs.add(contig)
                    progress.reading_contig("fragments", contig)

                fragment = fragment_from_read1(read)
                if fragment is None:
                    continue
                start, end, length = fragment
                if length < min_length or (max_length is not None and length > max_length):
                    continue
                if blacklist is not None and blacklist.overlaps(contig, start, end):
                    continue

                if bed_index is None:
                    labels = ["all"]
                else:
                    try:
                        bed_contig = resolve_contig_name(contig, list(bed_index), source_label="region BED")
                    except KeyError:
                        bed_contig = None
                    index = bed_index.get(bed_contig) if bed_contig is not None else None
                    midpoint = (start + end) // 2
                    labels = labels_at_position(index, midpoint) if index is not None else []

                    # Avoid double-counting when overlapping intervals share a label.
                    labels = list(dict.fromkeys(labels))
                    if len(labels) > 1:
                        summary.fragments_multiassigned += 1
                        if overlap_policy == "error":
                            raise ValueError(
                                f"Fragment midpoint {contig}:{midpoint} overlaps multiple labels: "
                                + ", ".join(labels)
                            )
                        if overlap_policy == "first":
                            labels = labels[:1]

                    if not labels:
                        summary.fragments_unassigned += 1
                        if include_unassigned:
                            labels = [unassigned_label]
                        else:
                            continue

                for label in labels:
                    counts[label][length] += 1
                summary.fragments_counted += 1

    return dict(counts), summary


def count_interval_fragment_lengths(
    fragment_paths: Sequence[str | Path],
    *,
    bed_index: Mapping[str, IntervalIndex] | None = None,
    contigs: set[str] | None = None,
    min_length: int = 1,
    max_length: int | None = None,
    overlap_policy: str = "all",
    include_unassigned: bool = False,
    unassigned_label: str = "unassigned",
    max_per_coordinate: int = 0,
    dedup_scope: str = "all_bams",
    chrom_sizes: str | None = None,
    blacklist_bed: str | None = None,
    progress: ProgressReporter | None = None,
) -> tuple[dict[str, Counter[int]], CountSummary]:
    """Count lengths from BED/BED.gz/bigBed fragment intervals."""
    source = IntervalFragmentSource([str(path) for path in fragment_paths], chrom_sizes=chrom_sizes)
    blacklist = load_blacklist(blacklist_bed, source.references, source.lengths)
    counts: dict[str, Counter[int]] = defaultdict(Counter)
    summary = CountSummary()
    seen_contigs: set[str] = set()
    try:
        available = set(source.references)
        if contigs is not None:
            resolved_contigs: set[str] = set()
            for requested in contigs:
                try:
                    resolved_contigs.add(resolve_contig_name(requested, source.references, source_label="fragment inputs"))
                except KeyError as exc:
                    raise ValueError(str(exc)) from exc
            contigs = resolved_contigs
        for chrom, start, end in source.iter_all(
            contigs=contigs,
            min_length=min_length,
            max_length=max_length,
            max_per_coordinate=max_per_coordinate,
            dedup_scope=dedup_scope,
        ):
            summary.reads_examined += 1
            if progress is not None and chrom not in seen_contigs:
                seen_contigs.add(chrom)
                progress.reading_contig("fragments", chrom)
            if blacklist is not None and blacklist.overlaps(chrom, start, end):
                continue
            length = end - start
            if bed_index is None:
                labels = ["all"]
            else:
                bed_chrom = chrom
                if chrom not in bed_index:
                    try:
                        bed_chrom = resolve_contig_name(
                            chrom, list(bed_index), source_label="region BED"
                        )
                    except KeyError:
                        bed_chrom = chrom
                index = bed_index.get(bed_chrom)
                midpoint = (start + end) // 2
                labels = labels_at_position(index, midpoint) if index is not None else []
                labels = list(dict.fromkeys(labels))
                if len(labels) > 1:
                    summary.fragments_multiassigned += 1
                    if overlap_policy == "error":
                        raise ValueError(
                            f"Fragment midpoint {chrom}:{midpoint} overlaps multiple labels: "
                            + ", ".join(labels)
                        )
                    if overlap_policy == "first":
                        labels = labels[:1]
                if not labels:
                    summary.fragments_unassigned += 1
                    if include_unassigned:
                        labels = [unassigned_label]
            if labels:
                for label in labels:
                    counts[label][length] += 1
                summary.fragments_counted += 1
    finally:
        source.close()
    return dict(counts), summary


def write_long_tsv(
    counts: Mapping[str, Counter[int]],
    output_path: str | Path,
    *,
    include_label: bool,
) -> None:
    """Write all distributions to a single tidy TSV."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt") as handle:
        if include_label:
            handle.write("label\tfragment_length\tcount\n")
            for label in sorted(counts):
                for length in sorted(counts[label]):
                    handle.write(f"{label}\t{length}\t{counts[label][length]}\n")
        else:
            handle.write("fragment_length\tcount\n")
            for length in sorted(counts.get("all", Counter())):
                handle.write(f"{length}\t{counts['all'][length]}\n")


def write_separate_tsvs(
    counts: Mapping[str, Counter[int]], output_dir: str | Path, prefix: str
) -> list[Path]:
    """Write one two-column TSV per label and return the created paths."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    used_names: set[str] = set()

    for label in sorted(counts):
        safe_label = sanitize_filename(label)
        candidate = safe_label
        suffix = 2
        while candidate in used_names:
            candidate = f"{safe_label}_{suffix}"
            suffix += 1
        used_names.add(candidate)

        path = output_dir / f"{prefix}_{candidate}_fragment_lengths.tsv"
        with path.open("wt") as handle:
            handle.write("fragment_length\tcount\n")
            for length in sorted(counts[label]):
                handle.write(f"{length}\t{counts[label][length]}\n")
        paths.append(path)
    return paths


def plot_distributions(
    counts: Mapping[str, Counter[int]],
    output_path: str | Path,
    *,
    minimum: int,
    maximum: int,
    density: bool,
) -> Path:
    """Plot fragment-length count or within-window density curves."""
    try:
        import matplotlib.pyplot as plt
        from nucleosuite.plotting import configure_unique_category_cycle
        configure_unique_category_cycle()
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib") from exc

    maximum = effective_plot_maximum(counts, maximum)
    if maximum < minimum:
        raise ValueError(
            "No fragments fall inside the requested plot range: the longest "
            f"counted fragment is {maximum} bp, below --plot-min {minimum}."
        )
    xs = list(range(minimum, maximum + 1))
    plotted = 0
    for label in sorted(counts):
        ys = [counts[label].get(x, 0) for x in xs]
        total = sum(ys)
        if total == 0:
            continue
        if density:
            ys = [value / total for value in ys]
        plt.plot(xs, ys, label=label, marker="o", markersize=2.0, markeredgewidth=0)
        plotted += 1

    if plotted == 0:
        raise ValueError("No fragments fall inside the requested plot range")
    axis = plt.gca()
    from nucleosuite.plotting import apply_base_pair_x_axis, apply_integer_y_axis
    apply_base_pair_x_axis(axis, xs)
    if not density:
        apply_integer_y_axis(axis)
    plt.xlabel("Fragment length (bp)")
    plt.ylabel("Density" if density else "Count")
    if plotted > 1:
        plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    plt.tight_layout()
    from nucleosuite.plotting import save_figure
    figure = plt.gcf()
    saved = save_figure(figure, output_path, default_dpi=200)
    plt.close(figure)
    return saved


def analyse_fragment_size_nrl(
    counts: Mapping[str, Counter[int]],
    *,
    minimum: int = 100,
    maximum: int = 1000,
    peak_resolution: float = 160.0,
) -> list[FragmentSizeNRLResult]:
    """Call multinucleosome fragment-size peaks and regress their summits.

    Peak detection uses the same resolution-derived smoothing and refinement as
    the DAC/DCC NRL command.  Each BED label is analysed independently.
    """

    from nucleosuite.nrl import (
        call_resolution_peaks,
        moving_average_by_distance,
        regress_peak_distances,
        resolution_smoothing_windows,
    )

    if minimum < 0:
        raise ValueError("--nrl-min-length must be zero or greater")
    if maximum <= minimum:
        raise ValueError("--nrl-max-length must be greater than --nrl-min-length")
    if peak_resolution < 0:
        raise ValueError("--nrl-peak-resolution must be zero or greater")

    detection_window, local_max_window = resolution_smoothing_windows(peak_resolution)
    results: list[FragmentSizeNRLResult] = []
    for label in sorted(counts):
        profile = counts[label]
        positive_lengths = [int(length) for length, count in profile.items() if count > 0]
        if not positive_lengths:
            continue
        profile_maximum = min(maximum, max(positive_lengths))
        if profile_maximum - minimum < 2:
            continue
        fragment_lengths = np.arange(minimum, profile_maximum + 1, dtype=np.float64)
        raw_counts = np.asarray(
            [profile.get(int(length), 0) for length in fragment_lengths],
            dtype=np.float64,
        )
        total = float(np.sum(raw_counts))
        if total <= 0:
            continue
        density = raw_counts / total
        detection_values = moving_average_by_distance(
            fragment_lengths, density, detection_window
        )
        local_values = moving_average_by_distance(
            fragment_lengths, density, local_max_window
        )
        peaks = tuple(
            call_resolution_peaks(
                fragment_lengths,
                density,
                local_values,
                detection_values,
                peak_resolution,
            )
        )
        regression = regress_peak_distances(peaks)
        results.append(
            FragmentSizeNRLResult(
                label=label,
                fragment_lengths=fragment_lengths,
                counts=raw_counts,
                density=density,
                local_values=local_values,
                detection_values=detection_values,
                peaks=peaks,
                regression=regression,
                minimum=minimum,
                maximum=profile_maximum,
                peak_resolution=peak_resolution,
                detection_window=detection_window,
                local_max_window=local_max_window,
            )
        )
    return results


def fragment_size_nrl_quality(result: FragmentSizeNRLResult) -> str:
    """Return a concise diagnostic status for a fragment-size NRL fit."""

    regression = result.regression
    if regression.n_peaks < 3:
        return "insufficient_peaks"
    if not math.isfinite(regression.slope):
        return "fit_failed"
    if not math.isfinite(regression.r_squared) or regression.r_squared < 0.9:
        return "low_r_squared"
    return "pass"


def _format_nrl_float(value: float) -> str:
    return "NaN" if not math.isfinite(value) else f"{value:.12g}"


def _fragment_size_output_prefixes(
    output_path: Path, results: Sequence[FragmentSizeNRLResult]
) -> list[tuple[FragmentSizeNRLResult, Path]]:
    base = output_path.with_suffix("") if output_path.suffix else output_path
    used: set[str] = set()
    labelled = len(results) > 1 or (results and results[0].label != "all")
    prefixed: list[tuple[FragmentSizeNRLResult, Path]] = []
    for result in results:
        suffix = ""
        if labelled:
            safe = sanitize_filename(result.label)
            candidate = safe
            number = 2
            while candidate in used:
                candidate = f"{safe}_{number}"
                number += 1
            used.add(candidate)
            suffix = f"_{candidate}"
        prefixed.append(
            (result, Path(f"{base}{suffix}_fragment_size_nrl"))
        )
    return prefixed


def write_fragment_size_nrl_outputs(
    results: Sequence[FragmentSizeNRLResult],
    output_path: str | Path,
    *,
    dpi: int = 200,
) -> list[Path]:
    """Write re-plottable fragment-size NRL tables and figures."""

    from nucleosuite.nrl import create_profile_plot, create_regression_plot

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    summary_path = Path(f"{output_path.with_suffix('')}_fragment_size_nrl_summary.tsv")
    with summary_path.open("w", encoding="utf-8", newline="") as summary_handle:
        summary_writer = csv.writer(summary_handle, delimiter="\t")
        summary_writer.writerow(
            [
                "label",
                "nrl_method",
                "min_fragment_length",
                "max_fragment_length",
                "peak_resolution_bp",
                "detection_smoothing_window",
                "local_max_smoothing_window",
                "peak_count",
                "nrl_bp",
                "intercept_bp",
                "r_squared",
                "slope_standard_error",
                "mean_adjacent_peak_spacing_bp",
                "quality_status",
            ]
        )

        for result, prefix in _fragment_size_output_prefixes(output_path, results):
            prefix.parent.mkdir(parents=True, exist_ok=True)
            profile_tsv = Path(f"{prefix}_profile.tsv")
            peaks_tsv = Path(f"{prefix}_peaks.tsv")
            regression_tsv = Path(f"{prefix}_regression.tsv")

            peak_number_by_index = {
                peak.index: number for number, peak in enumerate(result.peaks, start=1)
            }
            detection_indices = {
                peak.detection_index
                for peak in result.peaks
                if peak.detection_index is not None
            }
            with profile_tsv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
                        "label",
                        "fragment_length",
                        "count",
                        "unsmoothed_density",
                        "local_max_smoothed_density",
                        "detection_smoothed_density",
                        "detection_smoothing_window",
                        "local_max_smoothing_window",
                        "is_detection_peak",
                        "is_peak",
                        "peak_number",
                    ]
                )
                for index, (length, count, density, local, detection) in enumerate(
                    zip(
                        result.fragment_lengths,
                        result.counts,
                        result.density,
                        result.local_values,
                        result.detection_values,
                    )
                ):
                    peak_number = peak_number_by_index.get(index)
                    writer.writerow(
                        [
                            result.label,
                            f"{length:.12g}",
                            f"{count:.12g}",
                            f"{density:.12g}",
                            f"{local:.12g}",
                            f"{detection:.12g}",
                            result.detection_window,
                            result.local_max_window,
                            1 if index in detection_indices else 0,
                            1 if peak_number is not None else 0,
                            peak_number if peak_number is not None else "",
                        ]
                    )

            with peaks_tsv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
                        "label",
                        "peak_number",
                        "fragment_length",
                        "fitted_fragment_length",
                        "residual_bp",
                        "unsmoothed_density",
                        "local_max_smoothed_density",
                        "detection_peak_fragment_length",
                        "detection_smoothed_density",
                    ]
                )
                for number, peak in enumerate(result.peaks, start=1):
                    fitted = (
                        result.regression.intercept + result.regression.slope * number
                        if math.isfinite(result.regression.slope)
                        else float("nan")
                    )
                    writer.writerow(
                        [
                            result.label,
                            number,
                            _format_nrl_float(peak.distance),
                            _format_nrl_float(fitted),
                            _format_nrl_float(peak.distance - fitted),
                            _format_nrl_float(peak.raw_value),
                            _format_nrl_float(peak.smoothed_value),
                            _format_nrl_float(peak.detection_distance),
                            _format_nrl_float(peak.detection_value),
                        ]
                    )

            quality = fragment_size_nrl_quality(result)
            with regression_tsv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    [
                        "label",
                        "nrl_method",
                        "min_fragment_length",
                        "max_fragment_length",
                        "peak_resolution_bp",
                        "detection_smoothing_window",
                        "local_max_smoothing_window",
                        "peak_count",
                        "nrl_bp",
                        "intercept_bp",
                        "r_squared",
                        "slope_standard_error",
                        "mean_adjacent_peak_spacing_bp",
                        "quality_status",
                    ]
                )
                summary_row = [
                    result.label,
                    "fragment_size_distribution",
                    result.minimum,
                    result.maximum,
                    _format_nrl_float(result.peak_resolution),
                    result.detection_window,
                    result.local_max_window,
                    result.regression.n_peaks,
                    _format_nrl_float(result.regression.slope),
                    _format_nrl_float(result.regression.intercept),
                    _format_nrl_float(result.regression.r_squared),
                    _format_nrl_float(result.regression.slope_standard_error),
                    _format_nrl_float(result.regression.mean_adjacent_spacing),
                    quality,
                ]
                writer.writerow(summary_row)
                summary_writer.writerow(summary_row)

            title_label = (
                "Fragment-size NRL" if result.label == "all" else f"{result.label}: fragment-size NRL"
            )
            profile_plot = create_profile_plot(
                Path(f"{prefix}_profile.png"),
                result.fragment_lengths,
                result.density,
                result.local_values,
                result.detection_values,
                result.peaks,
                result.local_max_window,
                result.detection_window,
                "Fragment length",
                "Density",
                title_label,
                dpi,
            )
            regression_plot = create_regression_plot(
                Path(f"{prefix}_regression.png"),
                result.peaks,
                result.regression,
                f"{title_label} regression",
                dpi,
                y_label="Peak fragment length (bp)",
                slope_label="Fragment-size NRL",
            )
            outputs.extend(
                [profile_tsv, peaks_tsv, regression_tsv, profile_plot, regression_plot]
            )

    outputs.insert(0, summary_path)
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite fragment-lengths",
        description=(
            "Count BAM- or interval-derived fragment lengths, optionally stratified by "
            "the BED interval containing each fragment midpoint."
        ),
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-b", "--bam", "--bamfiles", dest="bamfiles", nargs="+",
        help="One or more paired-end BAM files. Chromosome-split BAMs are supported.",
    )
    inputs.add_argument(
        "--fragments", "--fragment-bed", dest="fragment_files", nargs="+",
        help="Fragment BED, BED.gz or bigBed file(s); only the first three columns are required.",
    )
    parser.add_argument("--chrom-sizes", help="Optional chromosome-size table, BAM or CRAM for fragment BED input.")
    parser.add_argument(
        "--blacklist-bed",
        help="Optional BED; complete fragments overlapping excluded regions are omitted.",
    )
    parser.add_argument(
        "--max-duplicates", dest="max_duplicates", type=int, default=0,
        help="Maximum identical fragment coordinates retained for interval input; 0 disables.",
    )
    parser.add_argument("--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams", help="Apply interval-coordinate limits across all fragment files or within each file (default: all_bams).")
    parser.add_argument(
        "--bed",
        help=(
            "Optional BED or BED.gz file. When supplied, fragment lengths are "
            "counted separately by BED label."
        ),
    )
    parser.add_argument(
        "--bed-label-column",
        type=int,
        default=4,
        help="One-based BED label column; use 0 for one pooled label (default: 4).",
    )
    parser.add_argument(
        "--bed-default-label",
        default="regions",
        help="Label used when --bed-label-column 0 (default: regions).",
    )
    parser.add_argument(
        "--overlap-policy",
        choices=("all", "first", "error"),
        default="all",
        help="How to handle midpoint overlaps with multiple labels (default: all).",
    )
    parser.add_argument(
        "--include-unassigned",
        action="store_true",
        help="Include fragments whose midpoint is outside all BED intervals.",
    )
    parser.add_argument(
        "--unassigned-label",
        default="unassigned",
        help="Label for fragments outside BED intervals (default: unassigned).",
    )
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=["all"],
        help=(
            "Exact contig names to include. Accepts space-separated values, comma lists, or 'all' "
            "(default: all)."
        ),
    )
    parser.add_argument("--mapq", type=int, default=0, help="Minimum MAPQ (default: 0).")
    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include reads marked as PCR/optical duplicates.",
    )
    parser.add_argument(
        "--allow-improper-pairs",
        action="store_true",
        help="Do not require the proper-pair SAM flag.",
    )
    parser.add_argument("--min-length", type=int, default=1, help="Minimum length (default: 1).")
    parser.add_argument("--max-length", type=int, help="Optional maximum fragment length.")
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Combined output TSV. Default: derived from the fragment input name(s), "
            "with the BED name appended when --bed is supplied."
        ),
    )
    parser.add_argument(
        "--separate-files",
        action="store_true",
        help="Also write one two-column TSV per BED label.",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for --separate-files (default: directory of --output).",
    )
    parser.add_argument("--plot", help="Optional plot output stem/path; final extension follows --plot-format (png or svg).")
    parser.add_argument("--plot-min", type=int, default=0, help="Plot minimum (default: 0).")
    parser.add_argument(
        "--plot-max",
        type=int,
        default=1000,
        help=(
            "Upper plot limit. The displayed axis stops at the longest counted "
            "fragment when that is shorter (default: 1000)."
        ),
    )
    parser.add_argument(
        "--plot-counts",
        action="store_true",
        help="Plot raw counts instead of within-window density.",
    )
    parser.add_argument(
        "--nrl-min-length",
        type=int,
        default=100,
        help="Minimum fragment length used for fragment-size NRL peak calling (default: 100).",
    )
    parser.add_argument(
        "--nrl-max-length",
        type=int,
        default=1000,
        help=(
            "Maximum fragment length used for fragment-size NRL peak calling; "
            "the longest counted fragment is used when shorter (default: 1000)."
        ),
    )
    parser.add_argument(
        "--nrl-peak-resolution",
        type=float,
        default=160.0,
        help=(
            "Minimum spacing of multinucleosome fragment-size peaks. This gives "
            "51 bp detection smoothing and 21 bp summit refinement at the default "
            "160 bp resolution."
        ),
    )
    parser.add_argument(
        "--no-fragment-size-nrl",
        action="store_true",
        help="Do not call multinucleosome fragment-size peaks or fit fragment-size NRL.",
    )
    parser.add_argument(
        "--no-console-histogram",
        action="store_true",
        help="Do not print proportional ASCII fragment-length histograms.",
    )
    parser.add_argument(
        "--histogram-width",
        type=int,
        default=100,
        help="Maximum ASCII histogram bar width (default: 100).",
    )
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser)
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def default_output_path(
    bam_paths: str | Path | Sequence[str | Path],
    bed_path: str | Path | None,
) -> Path:
    if isinstance(bam_paths, (str, Path)):
        paths = [Path(bam_paths)]
    else:
        paths = [Path(path) for path in bam_paths]
    if not paths:
        raise ValueError("At least one fragment input is required")
    if len(paths) == 1:
        input_stem = paths[0].name
        lower = input_stem.lower()
        for suffix in (".bed.gz", ".bigbed", ".bam", ".bed", ".bb", ".gz"):
            if lower.endswith(suffix):
                input_stem = input_stem[: -len(suffix)]
                break
        else:
            input_stem = Path(input_stem).stem
    else:
        input_stem = (
            "combined_bams"
            if all(path.name.lower().endswith(".bam") for path in paths)
            else "combined_fragments"
        )

    if bed_path is None:
        return Path(f"{input_stem}_fragment_lengths.tsv")

    bed_name = Path(bed_path).name
    lower = bed_name.lower()
    if lower.endswith(".bed.gz"):
        bed_stem = bed_name[:-7]
    elif lower.endswith(".bigbed"):
        bed_stem = bed_name[:-7]
    elif lower.endswith(".bb"):
        bed_stem = bed_name[:-3]
    elif lower.endswith(".bed"):
        bed_stem = bed_name[:-4]
    else:
        bed_stem = Path(bed_name).stem
    return Path(f"{input_stem}_{bed_stem}_fragment_lengths.tsv")


def run(args: argparse.Namespace) -> int:
    parser = build_parser()
    reporter = ProgressReporter("fragment-lengths")
    from nucleosuite.parallel import run_native_per_contig
    if not getattr(args, "_per_contig_worker", False) and int(getattr(args, "cores", 1) or 1) > 1:
        return run_native_per_contig("fragment-lengths", args, run)

    try:
        selected_contigs = parse_contigs(args.contigs)
        if args.bed:
            reporter.file_start("region annotations", args.bed)
        bed_index = (
            parse_bed(
                args.bed,
                label_column=args.bed_label_column,
                default_label=args.bed_default_label,
            )
            if args.bed
            else None
        )
        reporter.stage(
            "Counting fragment lengths from BAM inputs"
            if args.bamfiles
            else "Counting fragment lengths from interval inputs"
        )
        if args.bamfiles:
            counts, summary = count_fragment_lengths(
                args.bamfiles,
                bed_index=bed_index,
                contigs=selected_contigs,
                min_mapq=args.mapq,
                include_duplicates=args.include_duplicates,
                require_proper_pair=not args.allow_improper_pairs,
                min_length=args.min_length,
                max_length=args.max_length,
                overlap_policy=args.overlap_policy,
                include_unassigned=args.include_unassigned,
                unassigned_label=args.unassigned_label,
                blacklist_bed=args.blacklist_bed,
                progress=reporter,
            )
        else:
            counts, summary = count_interval_fragment_lengths(
                args.fragment_files,
                bed_index=bed_index,
                contigs=selected_contigs,
                min_length=args.min_length,
                max_length=args.max_length,
                overlap_policy=args.overlap_policy,
                include_unassigned=args.include_unassigned,
                unassigned_label=args.unassigned_label,
                max_per_coordinate=args.max_duplicates,
                dedup_scope=args.dedup_scope,
                chrom_sizes=args.chrom_sizes,
                blacklist_bed=args.blacklist_bed,
                progress=reporter,
            )
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))

    if args.histogram_width < 1:
        parser.error("--histogram-width must be at least 1")
    if args.plot_min < 0:
        parser.error("--plot-min must be zero or greater")
    if args.plot_max <= args.plot_min:
        parser.error("--plot-max must be greater than --plot-min")
    if args.nrl_min_length < 0:
        parser.error("--nrl-min-length must be zero or greater")
    if args.nrl_max_length <= args.nrl_min_length:
        parser.error("--nrl-max-length must be greater than --nrl-min-length")
    if args.nrl_peak_resolution < 0:
        parser.error("--nrl-peak-resolution must be zero or greater")

    output_path = (
        Path(args.output)
        if args.output
        else default_output_path(args.bamfiles or args.fragment_files, args.bed)
    )
    write_long_tsv(counts, output_path, include_label=args.bed is not None)
    reporter.stage(
        f"Writing fragment-length outputs: {summary.fragments_counted:,} fragments, "
        f"{len(counts):,} profiles"
    )
    print(f"Wrote: {output_path}")

    if not args.no_console_histogram:
        for label in sorted(counts):
            histogram_label = (
                str(output_path) if label == "all" and len(counts) == 1
                else f"{output_path} [{label}]"
            )
            print_fragment_length_histogram(
                counts[label],
                label=histogram_label,
                width=args.histogram_width,
            )

    if args.separate_files:
        output_dir = Path(args.output_dir) if args.output_dir else output_path.parent
        prefix = output_path.name
        if prefix.lower().endswith(".tsv"):
            prefix = prefix[:-4]
        paths = write_separate_tsvs(counts, output_dir, prefix)
        print(f"Wrote separate label files: {len(paths):,}")

    if args.plot:
        plot_output = plot_distributions(
            counts,
            args.plot,
            minimum=args.plot_min,
            maximum=args.plot_max,
            density=not args.plot_counts,
        )
        print(f"Wrote plot: {plot_output}")

    if not args.no_fragment_size_nrl:
        reporter.stage(
            "Calling multinucleosome fragment-size peaks "
            f"(resolution {args.nrl_peak_resolution:g} bp)"
        )
        nrl_results = analyse_fragment_size_nrl(
            counts,
            minimum=args.nrl_min_length,
            maximum=args.nrl_max_length,
            peak_resolution=args.nrl_peak_resolution,
        )
        if nrl_results:
            nrl_outputs = write_fragment_size_nrl_outputs(nrl_results, output_path)
            for result in nrl_results:
                quality = fragment_size_nrl_quality(result)
                if math.isfinite(result.regression.slope):
                    print(
                        f"Fragment-size NRL [{result.label}]: "
                        f"{result.regression.slope:.6g} bp "
                        f"({result.regression.n_peaks} peaks; "
                        f"R-squared {result.regression.r_squared:.6g}; {quality})"
                    )
                else:
                    print(
                        f"Fragment-size NRL [{result.label}]: NaN "
                        f"({result.regression.n_peaks} peaks; {quality})"
                    )
            for nrl_output in nrl_outputs:
                print(f"Wrote: {nrl_output}")
        else:
            print(
                "WARNING: fragment-size NRL was not written because no profile "
                "contains data across at least three positions in the requested range."
            )

    print(f"Reads examined: {summary.reads_examined:,}")
    print(f"Fragments counted: {summary.fragments_counted:,}")
    if args.bed:
        print(f"Fragments outside BED intervals: {summary.fragments_unassigned:,}")
        print(f"Fragments overlapping multiple labels: {summary.fragments_multiassigned:,}")
    print(f"Labels with counts: {len(counts):,}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
