"""Per-base BigWig alignment, metaprofile and heatmap generation."""

from __future__ import annotations

import logging
import csv
import math
from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Literal

import numpy as np

from .io import open_text, strip_known_suffix
from .core.regions import resolve_contig_name
from nucleosuite.progress import ProgressReporter
from nucleosuite.nrl import Peak, Regression

LOGGER = logging.getLogger(__name__)
SortMode = Literal["center", "rise_after_min", "mean_absolute", "absmean", "max", "unsorted"]
SubsampleMode = Literal["first", "random"]
MissingStrandMode = Literal["forward", "random", "error"]


@dataclass(frozen=True)
class ChromosomeIntervals:
    starts: list[int]
    intervals: list[tuple[int, int]]


@dataclass(frozen=True)
class AlignmentConfig:
    bigwig: Path
    region_bed: Path
    blacklist_bed: Path | None = None
    nucleosome_bed: Path | None = None
    nucleosome_offset: int = 1
    state_bed: Path | None = None
    output_dir: Path = Path(".")
    output_prefix: str | None = None
    heatmap_output: Path | None = None
    heatmap_matrix_output: Path | None = None
    aggregate_output: Path | None = None
    plotted_mean_output: Path | None = None
    mean_plot_output: Path | None = None
    summary_output: Path | None = None
    window_half: int = 2500
    chrom_col: int = 1
    start_col: int = 2
    end_col: int = 3
    strand_col: int = 6
    point_col: int = 0
    skip_header: bool = False
    missing_strand: MissingStrandMode = "forward"
    zero_thresh: int = 5
    max_score: float | None = 300.0
    nan_to_zero: bool = True
    max_heatmap_rows: int | None = None
    subsample_mode: SubsampleMode = "first"
    stop_after_valid: int | None = None
    seed: int | None = None
    breadth: float = 1.0
    vmin: float | None = None
    vmax: float | None = None
    sort_mode: SortMode = "mean_absolute"
    axis_label: str = "Distance from reference-site centre"
    mean_ylim: float | None = None
    colorbar_label: str = "Score"
    mean_ylabel: str = "Mean score"
    dpi: int = 300
    nrl: bool = True
    nrl_peak_resolution: float = 160.0
    nrl_regression_min: float = 0.0
    nrl_regression_max: float | None = None
    nrl_exclusion: bool = True
    nrl_regression_exclusion_start: float | None = None
    nrl_regression_exclusion_end: float | None = None


@dataclass(frozen=True)
class AggregateNRLResult:
    """Unified aggregate peak calls and two outward-direction regressions."""

    positions: np.ndarray
    raw_values: np.ndarray
    local_values: np.ndarray
    detection_values: np.ndarray
    peaks: tuple[Peak, ...]
    central_peak: Peak | None
    positive_peaks: tuple[Peak, ...]
    negative_peaks: tuple[Peak, ...]
    positive_peak_numbers: tuple[int, ...]
    negative_peak_numbers: tuple[int, ...]
    positive_regression: Regression
    negative_regression: Regression
    peak_resolution: float
    detection_window: int
    local_max_window: int
    regression_min: float
    regression_max: float
    exclusion_start: float | None
    exclusion_end: float | None


@dataclass
class AlignmentStats:
    input_lines: int = 0
    skipped_header: int = 0
    skipped_blank: int = 0
    skipped_comment: int = 0
    skipped_missing_chromosome: int = 0
    skipped_chromosome_not_in_bigwig: int = 0
    skipped_invalid_point: int = 0
    skipped_invalid_coordinates: int = 0
    skipped_invalid_interval: int = 0
    skipped_invalid_interval_for_state_filter: int = 0
    skipped_no_state_overlap: int = 0
    skipped_blacklisted_anchor: int = 0
    skipped_missing_requested_nucleosome: int = 0
    skipped_missing_or_invalid_strand: int = 0
    randomly_oriented: int = 0
    unstranded_kept_forward: int = 0
    skipped_bigwig_error: int = 0
    skipped_nan: int = 0
    skipped_nonfinite: int = 0
    skipped_consecutive_zeros: int = 0
    skipped_above_max_score: int = 0
    nan_converted_vectors: int = 0
    blacklisted_window_bases: int = 0
    aligned_to_region: int = 0
    aligned_to_nucleosome: int = 0
    valid_total: int = 0
    selected_for_plot: int = 0
    stopped_after_valid_limit: int = 0

    def increment(self, name: str, amount: int = 1) -> None:
        if name not in {item.name for item in fields(self)}:
            raise KeyError(f"Unknown statistic: {name}")
        setattr(self, name, getattr(self, name) + amount)


def no_valid_regions_message(stats: AlignmentStats, config: AlignmentConfig) -> str:
    """Describe why aggregate windows were rejected and suggest likely fixes."""

    rejection_fields = (
        ("missing chromosome field", stats.skipped_missing_chromosome),
        ("missing chromosome", stats.skipped_chromosome_not_in_bigwig),
        ("invalid or out-of-range coordinates", stats.skipped_invalid_coordinates),
        ("missing BigWig data", stats.skipped_nan),
        ("non-finite signal", stats.skipped_nonfinite),
        ("consecutive zeros", stats.skipped_consecutive_zeros),
        ("score above maximum", stats.skipped_above_max_score),
        ("BigWig read error", stats.skipped_bigwig_error),
        ("blacklisted anchor", stats.skipped_blacklisted_anchor),
        ("missing requested nucleosome", stats.skipped_missing_requested_nucleosome),
        ("no state overlap", stats.skipped_no_state_overlap),
    )
    observed = [f"{label}={count:,}" for label, count in rejection_fields if count]
    lines = ["No valid regions remained after aggregate filtering."]
    if observed:
        lines.append("Rejection counts: " + "; ".join(observed) + ".")
    else:
        lines.append(f"Region lines examined: {stats.input_lines:,}.")
    suggestions: list[str] = []
    if stats.skipped_consecutive_zeros:
        suggestions.append(
            f"--zero-thresh 0 disables consecutive-zero rejection (current: {config.zero_thresh})"
        )
    if stats.skipped_above_max_score:
        current = "inf" if config.max_score is None or math.isinf(config.max_score) else f"{config.max_score:g}"
        suggestions.append(
            f"--max-score inf disables upper-signal rejection (current: {current})"
        )
    if suggestions:
        lines.append("Possible filter adjustments: " + "; ".join(suggestions) + ".")
    return " ".join(lines)


def validate_config(config: AlignmentConfig) -> None:
    for path, label in (
        (config.bigwig, "BigWig"),
        (config.region_bed, "region BED"),
        (config.blacklist_bed, "blacklist BED"),
        (config.nucleosome_bed, "nucleosome BED"),
        (config.state_bed, "state BED"),
    ):
        if path is not None and not Path(path).is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
    if config.nucleosome_offset == 0:
        raise ValueError("nucleosome_offset cannot be 0")
    for name in ("chrom_col", "start_col", "end_col", "strand_col"):
        if getattr(config, name) < 1:
            raise ValueError(f"{name} must be >= 1")
    if config.point_col < 0 or config.window_half < 0 or config.zero_thresh < 0:
        raise ValueError("point_col, window_half and zero_thresh must be >= 0")
    if config.max_heatmap_rows is not None and config.max_heatmap_rows < 1:
        raise ValueError("max_heatmap_rows must be >= 1")
    if config.stop_after_valid is not None and config.stop_after_valid < 1:
        raise ValueError("stop_after_valid must be >= 1")
    if not 0 < config.breadth <= 1:
        raise ValueError("breadth must be > 0 and <= 1")
    if config.mean_ylim is not None and config.mean_ylim <= 0:
        raise ValueError("mean_ylim must be > 0")
    if config.dpi < 1:
        raise ValueError("dpi must be >= 1")
    if config.nrl_peak_resolution < 0:
        raise ValueError("nrl_peak_resolution must be >= 0")
    if config.nrl_regression_min < 0:
        raise ValueError("nrl_regression_min must be >= 0")
    if (
        config.nrl_regression_max is not None
        and config.nrl_regression_max <= config.nrl_regression_min
    ):
        raise ValueError("nrl_regression_max must be greater than nrl_regression_min")
    if (
        config.nrl_regression_max is not None
        and config.nrl_regression_max > config.window_half
    ):
        raise ValueError("nrl_regression_max cannot exceed window_half")
    exclusion_values = (
        config.nrl_regression_exclusion_start,
        config.nrl_regression_exclusion_end,
    )
    if (exclusion_values[0] is None) != (exclusion_values[1] is None):
        raise ValueError(
            "nrl_regression_exclusion_start and nrl_regression_exclusion_end "
            "must be supplied together"
        )
    if exclusion_values[0] is not None and exclusion_values[1] is not None:
        exclusion_start, exclusion_end = exclusion_values
        if not math.isfinite(exclusion_start) or not math.isfinite(exclusion_end):
            raise ValueError("NRL regression exclusion positions must be finite")
        if exclusion_end < exclusion_start:
            raise ValueError(
                "nrl_regression_exclusion_end must be greater than or equal to "
                "nrl_regression_exclusion_start"
            )
        if exclusion_start < -config.window_half or exclusion_end > config.window_half:
            raise ValueError(
                "NRL regression exclusion positions must lie within the aggregate window"
            )


def infer_aggregate_track_options(path: str | Path) -> dict[str, object]:
    """Infer filtering and labels from NucleoSuite BigWig filename suffixes."""

    name = Path(path).name.lower()
    suffixes: tuple[tuple[str, str, str], ...] = (
        ("_sm_mwps", "Smoothed median-adjusted WPS", "Mean smoothed median-adjusted WPS"),
        ("_mwps", "Median-adjusted WPS", "Mean median-adjusted WPS"),
        ("_pospns", "Positive probabilistic nucleosome score (posPNS)", "Mean positive probabilistic nucleosome score (posPNS)"),
        ("_pns_smoothed", "Probabilistic nucleosome score (PNS)", "Mean probabilistic nucleosome score (PNS)"),
        ("_pns", "Probabilistic nucleosome score (PNS)", "Mean probabilistic nucleosome score (PNS)"),
        ("_posbns", "Positive boxcar nucleosome score (posBNS)", "Mean positive boxcar nucleosome score (posBNS)"),
        ("_bns_smoothed", "Boxcar nucleosome score (BNS)", "Mean boxcar nucleosome score (BNS)"),
        ("_bns", "Boxcar nucleosome score (BNS)", "Mean boxcar nucleosome score (BNS)"),
        ("_postns", "Positive triangular nucleosome score (posTNS)", "Mean positive triangular nucleosome score (posTNS)"),
        ("_tns_smoothed", "Triangular nucleosome score (TNS)", "Mean triangular nucleosome score (TNS)"),
        ("_tns", "Triangular nucleosome score (TNS)", "Mean triangular nucleosome score (TNS)"),
        ("_wps", "Windowed protection score (WPS)", "Mean windowed protection score (WPS)"),
        ("_fragment_left_ends", "Left fragment-end count", "Mean left fragment-end count"),
        ("_fragment_right_ends", "Right fragment-end count", "Mean right fragment-end count"),
        ("_fragment_ends", "Fragment-end count", "Mean fragment-end count"),
        ("_fragment_coverage", "Fragment coverage", "Mean fragment coverage"),
        ("_coverage", "Fragment coverage", "Mean fragment coverage"),
        ("_dyad", "Dyad count", "Mean dyad count"),
    )
    stem = name
    for extension in (".bigwig", ".bw"):
        if stem.endswith(extension):
            stem = stem[: -len(extension)]
            break
    inferred: dict[str, object] = {
        "track_type": "generic",
        "zero_thresh": 5,
        "max_score": 300.0,
        "colorbar_label": "Score",
        "mean_ylabel": "Mean score",
    }
    for suffix, colorbar, mean_ylabel in suffixes:
        if stem.endswith(suffix):
            inferred.update(
                track_type=suffix.removeprefix("_"),
                colorbar_label=colorbar,
                mean_ylabel=mean_ylabel,
            )
            if suffix == "_dyad":
                inferred.update(zero_thresh=0, max_score=float("inf"))
            break
    return inferred


def resolve_nrl_exclusion(config: AlignmentConfig) -> tuple[float | None, float | None]:
    """Return the effective regression-only exclusion interval."""

    start = config.nrl_regression_exclusion_start
    end = config.nrl_regression_exclusion_end
    if start is not None and end is not None:
        return float(start), float(end)
    if not config.nrl_exclusion:
        return None, None
    half_resolution = float(config.nrl_peak_resolution) / 2.0
    return -half_resolution, half_resolution


def make_output_prefix(config: AlignmentConfig) -> str:
    if config.output_prefix:
        base = config.output_prefix
    else:
        parts = [strip_known_suffix(config.region_bed), strip_known_suffix(config.bigwig)]
        if config.nucleosome_bed is not None:
            label = (
                f"plus{config.nucleosome_offset}"
                if config.nucleosome_offset > 0
                else f"minus{abs(config.nucleosome_offset)}"
            )
            parts.extend([strip_known_suffix(config.nucleosome_bed), label])
        if config.state_bed is not None:
            parts.append(strip_known_suffix(config.state_bed))
        base = "_".join(parts)
    from nucleosuite.output_naming import parameter_range, parameterized_prefix

    parameters: list[tuple[str, object]] = [
        ("win", config.window_half),
        ("zero", config.zero_thresh),
        ("maxscore", config.max_score),
        ("missing", "zero" if config.nan_to_zero else "reject"),
        ("sort", config.sort_mode),
    ]
    if config.nrl:
        exclusion_start, exclusion_end = resolve_nrl_exclusion(config)
        regression_max = (
            config.window_half
            if config.nrl_regression_max is None
            else config.nrl_regression_max
        )
        parameters.extend(
            [
                ("nrlres", config.nrl_peak_resolution),
                ("nrlmin", config.nrl_regression_min),
                ("nrlmax", regression_max),
                (
                    "excl",
                    "none"
                    if exclusion_start is None
                    else parameter_range(exclusion_start, exclusion_end),
                ),
            ]
        )
    else:
        parameters.append(("nrl", "off"))
    return str(parameterized_prefix(base, parameters))


def resolve_output_paths(config: AlignmentConfig) -> dict[str, Path]:
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = make_output_prefix(config)
    from nucleosuite.plotting import plot_path
    outputs = {
        "heatmap": plot_path(config.heatmap_output or output_dir / f"{prefix}_heatmap.png"),
        "heatmap_matrix": config.heatmap_matrix_output or output_dir / f"{prefix}_heatmap_matrix.tsv.gz",
        "aggregate": config.aggregate_output or output_dir / f"{prefix}_aggregate_all.tsv",
        "plotted_mean": config.plotted_mean_output or output_dir / f"{prefix}_heatmap_mean.tsv",
        "mean_plot": plot_path(config.mean_plot_output or output_dir / f"{prefix}_heatmap_mean.png"),
        "summary": config.summary_output or output_dir / f"{prefix}_summary.tsv",
        "row_metadata": output_dir / f"{prefix}_heatmap_rows.tsv",
    }
    if config.nrl:
        outputs.update(
            {
                "nrl_profile": output_dir / f"{prefix}_aggregate_nrl_profile.tsv",
                "nrl_peaks": output_dir / f"{prefix}_aggregate_nrl_peaks.tsv",
                "nrl_summary": output_dir / f"{prefix}_aggregate_nrl_summary.tsv",
                "nrl_profile_plot": plot_path(
                    output_dir / f"{prefix}_aggregate_nrl_profile.png"
                ),
                "nrl_positive_regression": output_dir
                / f"{prefix}_aggregate_nrl_positive_regression.tsv",
                "nrl_positive_regression_plot": plot_path(
                    output_dir / f"{prefix}_aggregate_nrl_positive_regression.png"
                ),
                "nrl_negative_regression": output_dir
                / f"{prefix}_aggregate_nrl_negative_regression.tsv",
                "nrl_negative_regression_plot": plot_path(
                    output_dir / f"{prefix}_aggregate_nrl_negative_regression.png"
                ),
            }
        )
    return outputs


def parse_integer(value: str | None) -> int | None:
    if value is None or value.strip() == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def get_field(values: list[str], one_based_column: int) -> str | None:
    index = one_based_column - 1
    return values[index] if 0 <= index < len(values) else None


def load_bed_intervals(path: str | Path) -> dict[str, ChromosomeIntervals]:
    grouped: dict[str, list[tuple[int, int]]] = {}
    with open_text(path) as handle:
        for raw in handle:
            values = raw.strip().split()
            if not values or raw.lstrip().startswith("#") or len(values) < 3:
                continue
            start = parse_integer(values[1])
            end = parse_integer(values[2])
            if start is None or end is None or start < 0 or end <= start:
                continue
            grouped.setdefault(values[0], []).append((start, end))

    indexed: dict[str, ChromosomeIntervals] = {}
    for chromosome, intervals in grouped.items():
        intervals.sort()
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        indexed[chromosome] = ChromosomeIntervals(
            starts=[start for start, _ in merged], intervals=merged
        )
    return indexed


def overlaps_any(
    chromosome_data: ChromosomeIntervals | None,
    region_start: int,
    region_end: int,
) -> bool:
    if chromosome_data is None:
        return False
    index = bisect_left(chromosome_data.starts, region_end)
    if index == 0:
        return False
    start, end = chromosome_data.intervals[index - 1]
    return end > region_start and start < region_end


def load_nucleosome_centers(path: str | Path) -> dict[str, list[int]]:
    centers: dict[str, list[int]] = {}
    with open_text(path) as handle:
        for raw in handle:
            values = raw.strip().split()
            if not values or raw.lstrip().startswith("#") or len(values) < 3:
                continue
            start = parse_integer(values[1])
            end = parse_integer(values[2])
            if start is None or end is None or start < 0 or end < start:
                continue
            centers.setdefault(values[0], []).append((start + end) // 2)
    for chromosome_centers in centers.values():
        chromosome_centers.sort()
    return centers


def find_relative_nucleosome_center(
    centers: list[int] | None,
    region_center: int,
    strand: str,
    signed_offset: int,
) -> int | None:
    """Find a strict +N/-N nucleosome relative to region orientation."""
    if not centers or signed_offset == 0:
        return None
    relative_direction = 1 if signed_offset > 0 else -1
    genomic_direction = relative_direction if strand == "+" else -relative_direction
    rank = abs(signed_offset)
    if genomic_direction > 0:
        target = bisect_right(centers, region_center) + rank - 1
    else:
        target = bisect_left(centers, region_center) - rank
    return centers[target] if 0 <= target < len(centers) else None


def has_consecutive_zeros(values: np.ndarray, threshold: int) -> bool:
    if threshold <= 0 or values.size < threshold:
        return False
    zero_indices = np.flatnonzero(values == 0.0)
    if zero_indices.size < threshold:
        return False
    run = 1
    previous = int(zero_indices[0])
    for current_raw in zero_indices[1:]:
        current = int(current_raw)
        run = run + 1 if current == previous + 1 else 1
        if run >= threshold:
            return True
        previous = current
    return False


def extract_window(
    bigwig: Any,
    chromosome: str,
    center: int,
    window_half: int,
    chromosome_length: int,
) -> np.ndarray:
    length = 2 * window_half + 1
    requested_start = center - window_half
    requested_end = center + window_half + 1
    output = np.full(length, np.nan, dtype=np.float64)
    clipped_start = max(0, requested_start)
    clipped_end = min(chromosome_length, requested_end)
    if clipped_start >= clipped_end:
        return output
    extracted = np.asarray(
        bigwig.values(chromosome, clipped_start, clipped_end, numpy=True),
        dtype=np.float64,
    )
    destination_start = clipped_start - requested_start
    output[destination_start:destination_start + extracted.size] = extracted
    return output


def add_heatmap_row(
    selected_rows: list[np.ndarray],
    values: np.ndarray,
    valid_index: int,
    max_rows: int | None,
    mode: SubsampleMode,
    rng: np.random.Generator,
) -> None:
    if max_rows is None:
        selected_rows.append(values.copy())
    elif mode == "first":
        if len(selected_rows) < max_rows:
            selected_rows.append(values.copy())
    elif len(selected_rows) < max_rows:
        selected_rows.append(values.copy())
    else:
        replacement = int(rng.integers(0, valid_index))
        if replacement < max_rows:
            selected_rows[replacement] = values.copy()


def central_crop(matrix: np.ndarray, breadth: float) -> tuple[np.ndarray, np.ndarray]:
    columns = matrix.shape[1]
    x_values = np.arange(-(columns // 2), -(columns // 2) + columns)
    if breadth >= 1.0:
        return matrix, x_values
    keep = max(1, min(int(columns * breadth), columns))
    if keep % 2 == 0:
        keep = min(columns, keep + 1)
    centre = columns // 2
    half = keep // 2
    start = centre - half
    end = centre + half + 1
    return matrix[:, start:end], x_values[start:end]


def compute_sort_key_rise_after_min(row: np.ndarray, center_index: int) -> float:
    left = row[:center_index + 1]
    minimum_index = int(np.argmin(left))
    if row[minimum_index] >= 0:
        return np.inf
    crossings = np.flatnonzero(row[minimum_index + 1:] >= 0)
    return (
        float(minimum_index + 1 + int(crossings[0]) - minimum_index)
        if crossings.size
        else np.inf
    )


def sort_matrix(matrix: np.ndarray, mode: SortMode) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sort rows and return the matrix, original-row order and sorting scores."""
    center = matrix.shape[1] // 2
    if mode == "center":
        scores = matrix[:, center].astype(float)
        order = np.argsort(-scores, kind="stable")
    elif mode == "rise_after_min":
        scores = np.asarray(
            [compute_sort_key_rise_after_min(row, center) for row in matrix],
            dtype=np.float64,
        )
        order = np.lexsort((-matrix[:, center], scores))
    elif mode in {"mean_absolute", "absmean"}:
        scores = np.nanmean(np.abs(matrix), axis=1)
        order = np.argsort(-scores, kind="stable")
    elif mode == "max":
        scores = np.nanmax(matrix, axis=1)
        order = np.argsort(-scores, kind="stable")
    else:
        scores = np.full(matrix.shape[0], np.nan, dtype=float)
        order = np.arange(matrix.shape[0])
    return matrix[order], order, scores[order]


def write_heatmap_row_metadata(
    path: Path,
    original_order: np.ndarray,
    sort_scores: np.ndarray,
    sort_mode: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8") as output:
        output.write("heatmap_row\toriginal_selected_row\tsort_mode\tsort_score\n")
        for heatmap_row, (original_index, score) in enumerate(
            zip(original_order, sort_scores), start=1
        ):
            score_text = "NaN" if not np.isfinite(score) else f"{float(score):.10g}"
            output.write(
                f"{heatmap_row}\t{int(original_index) + 1}\t{sort_mode}\t{score_text}\n"
            )

def write_heatmap_matrix(matrix: np.ndarray, x_values: np.ndarray, path: Path) -> None:
    """Write the exact sorted/subsampled matrix used for the heatmap plot."""

    import gzip

    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8") as output:
        header = ["row_index", *[str(int(value)) for value in x_values]]
        output.write("\t".join(header) + "\n")
        for row_index, row in enumerate(matrix, start=1):
            output.write(
                str(row_index)
                + "\t"
                + "\t".join(f"{float(value):.10g}" for value in row)
                + "\n"
            )


def write_profile(values: np.ndarray, path: Path, x_values: np.ndarray | None = None) -> None:
    if x_values is None:
        centre = values.size // 2
        x_values = np.arange(-centre, -centre + values.size)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8") as output:
        output.write("relative_position\tscore\n")
        for position, score in zip(x_values, values):
            output.write(f"{int(position)}\t{float(score):.6f}\n")


def _outward_peak(peak: Peak) -> Peak:
    """Return a peak whose distance fields are absolute offsets from position 0."""

    return Peak(
        index=peak.index,
        distance=abs(peak.distance),
        raw_value=peak.raw_value,
        smoothed_value=peak.smoothed_value,
        rise_bp=peak.rise_bp,
        fall_bp=peak.fall_bp,
        detection_index=peak.detection_index,
        detection_distance=abs(peak.detection_distance),
        detection_value=peak.detection_value,
    )


def analyse_aggregate_nrl(
    values: np.ndarray,
    *,
    positions: np.ndarray | None = None,
    peak_resolution: float = 160.0,
    regression_min: float = 0.0,
    regression_max: float | None = None,
    exclusion_start: float | None = None,
    exclusion_end: float | None = None,
) -> AggregateNRLResult:
    """Call peaks once across a complete aggregate and fit each side separately."""

    from nucleosuite.nrl import (
        call_resolution_peaks,
        moving_average_by_distance,
        regress_peak_distances,
        resolution_smoothing_windows,
    )

    values = np.asarray(values, dtype=np.float64)
    if positions is None:
        centre = values.size // 2
        positions = np.arange(-centre, -centre + values.size, dtype=np.float64)
    else:
        positions = np.asarray(positions, dtype=np.float64)
    if values.shape != positions.shape:
        raise ValueError("Aggregate positions and values must have the same shape")
    if peak_resolution < 0:
        raise ValueError("peak_resolution must be zero or greater")
    if regression_min < 0:
        raise ValueError("regression_min must be zero or greater")
    if (exclusion_start is None) != (exclusion_end is None):
        raise ValueError("exclusion_start and exclusion_end must be supplied together")
    if exclusion_start is not None and exclusion_end is not None:
        if not math.isfinite(exclusion_start) or not math.isfinite(exclusion_end):
            raise ValueError("Regression exclusion positions must be finite")
        if exclusion_end < exclusion_start:
            raise ValueError(
                "exclusion_end must be greater than or equal to exclusion_start"
            )

    finite = np.isfinite(positions) & np.isfinite(values)
    x = positions[finite]
    raw = values[finite]
    if x.size < 3:
        raise ValueError("Aggregate NRL analysis requires at least three finite positions")
    order = np.argsort(x, kind="stable")
    x, raw = x[order], raw[order]
    if np.any(np.diff(x) <= 0):
        raise ValueError("Aggregate relative positions must be unique and increasing")

    available_max = float(max(abs(float(x[0])), abs(float(x[-1]))))
    effective_max = available_max if regression_max is None else float(regression_max)
    if effective_max <= regression_min:
        raise ValueError("regression_max must be greater than regression_min")
    if effective_max > available_max:
        raise ValueError("regression_max exceeds the aggregate alignment range")
    if exclusion_start is not None and exclusion_end is not None:
        if exclusion_start < float(x[0]) or exclusion_end > float(x[-1]):
            raise ValueError(
                "Regression exclusion positions must lie within the aggregate alignment range"
            )

    detection_window, local_max_window = resolution_smoothing_windows(peak_resolution)
    # Smoothing and peak calling deliberately operate across the complete signed
    # profile. The profile is divided by direction only after peaks are called.
    detection = moving_average_by_distance(x, raw, detection_window)
    local = moving_average_by_distance(x, raw, local_max_window)
    peaks = tuple(
        call_resolution_peaks(x, raw, local, detection, peak_resolution)
    )

    central_radius = float(peak_resolution) / 2.0 if peak_resolution > 0 else 0.0
    central_candidates = [
        peak for peak in peaks if abs(peak.distance) <= central_radius
    ]
    central_peak = (
        min(
            central_candidates,
            key=lambda peak: (
                abs(peak.distance),
                -peak.smoothed_value,
                -peak.raw_value,
                peak.distance,
            ),
        )
        if central_candidates
        else None
    )

    def excluded(peak: Peak) -> bool:
        return bool(
            exclusion_start is not None
            and exclusion_end is not None
            and exclusion_start <= peak.distance <= exclusion_end
        )

    def selected_outward(direction: str) -> tuple[tuple[Peak, ...], tuple[int, ...]]:
        directional = [
            peak
            for peak in peaks
            if peak is not central_peak
            and (peak.distance > 0 if direction == "positive" else peak.distance < 0)
        ]
        directional.sort(key=lambda peak: abs(peak.distance))
        numbered: list[tuple[int, Peak]] = []
        if central_peak is not None:
            numbered.append((0, central_peak))
        numbered.extend(enumerate(directional, start=1))
        selected = [
            (number, peak)
            for number, peak in numbered
            if regression_min <= abs(peak.distance) <= effective_max
            and not excluded(peak)
        ]
        return (
            tuple(_outward_peak(peak) for _, peak in selected),
            tuple(number for number, _ in selected),
        )

    positive, positive_numbers = selected_outward("positive")
    negative, negative_numbers = selected_outward("negative")
    return AggregateNRLResult(
        positions=x,
        raw_values=raw,
        local_values=local,
        detection_values=detection,
        peaks=peaks,
        central_peak=central_peak,
        positive_peaks=positive,
        negative_peaks=negative,
        positive_peak_numbers=positive_numbers,
        negative_peak_numbers=negative_numbers,
        positive_regression=regress_peak_distances(
            positive, peak_numbers=positive_numbers
        ),
        negative_regression=regress_peak_distances(
            negative, peak_numbers=negative_numbers
        ),
        peak_resolution=float(peak_resolution),
        detection_window=detection_window,
        local_max_window=local_max_window,
        regression_min=float(regression_min),
        regression_max=effective_max,
        exclusion_start=exclusion_start,
        exclusion_end=exclusion_end,
    )


def _format_nrl_value(value: float) -> str:
    return "NaN" if not math.isfinite(value) else f"{value:.12g}"


def _aggregate_nrl_quality(regression: Regression) -> str:
    if regression.n_peaks < 3:
        return "insufficient_peaks"
    if not math.isfinite(regression.slope):
        return "fit_failed"
    if not math.isfinite(regression.r_squared) or regression.r_squared < 0.9:
        return "low_r_squared"
    return "pass"


def _create_aggregate_nrl_profile_plot(
    path: Path,
    result: AggregateNRLResult,
    *,
    title: str,
    x_label: str,
    y_label: str,
    dpi: int,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(10, 5.5))
    if result.exclusion_start is not None and result.exclusion_end is not None:
        axis.axvspan(
            result.exclusion_start,
            result.exclusion_end,
            color="0.9",
            alpha=0.55,
            linewidth=0,
            label="Regression exclusion zone",
            zorder=0,
        )
    axis.plot(
        result.positions,
        result.raw_values,
        color="0.72",
        linewidth=0.9,
        label="Unsmoothed",
    )
    axis.plot(
        result.positions,
        result.local_values,
        color="black",
        linewidth=1.5,
        label=(
            f"Local maxima ({result.local_max_window} bp)"
            if result.local_max_window > 1
            else "Local maxima signal"
        ),
    )
    if result.detection_window != result.local_max_window:
        axis.plot(
            result.positions,
            result.detection_values,
            color="0.4",
            linewidth=1.2,
            linestyle="--",
            label=(
                f"Peak detection ({result.detection_window} bp)"
                if result.detection_window > 1
                else "Peak detection signal"
            ),
        )
    if result.peaks:
        axis.scatter(
            [peak.distance for peak in result.peaks],
            [peak.smoothed_value for peak in result.peaks],
            s=28,
            facecolors="white",
            edgecolors="black",
            linewidths=1.0,
            zorder=5,
            label="Called peaks",
        )
    axis.axvline(0, color="0.55", linewidth=0.8, zorder=0)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_title(title)
    axis.legend(frameon=False)
    from nucleosuite.plotting import apply_distance_x_axis, annotate_points, get_plot_options

    apply_distance_x_axis(axis)
    if result.peaks:
        annotate_points(
            axis,
            [peak.distance for peak in result.peaks],
            [peak.smoothed_value for peak in result.peaks],
            points_are_peaks=True,
            options=get_plot_options(),
        )
    from nucleosuite.plotting import save_figure

    figure.tight_layout()
    saved = save_figure(figure, path, default_dpi=dpi)
    plt.close(figure)
    return saved


def _write_direction_regression(
    path: Path,
    *,
    direction: str,
    peaks: tuple[Peak, ...],
    peak_numbers: tuple[int, ...],
    signed_by_index: dict[int, float],
    regression: Regression,
    plot_title: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "direction",
                "peak_number",
                "signed_position_bp",
                "distance_from_zero_bp",
                "fitted_distance_from_zero_bp",
                "residual_bp",
                "plot_title",
            ]
        )
        if not peaks:
            writer.writerow([direction, "", "NaN", "NaN", "NaN", "NaN", plot_title])
            return
        for number, peak in zip(peak_numbers, peaks):
            fitted = (
                regression.intercept + regression.slope * number
                if math.isfinite(regression.slope)
                else float("nan")
            )
            writer.writerow(
                [
                    direction,
                    number,
                    _format_nrl_value(signed_by_index[peak.index]),
                    _format_nrl_value(peak.distance),
                    _format_nrl_value(fitted),
                    _format_nrl_value(peak.distance - fitted),
                    plot_title,
                ]
            )


def write_aggregate_nrl_outputs(
    result: AggregateNRLResult,
    config: AlignmentConfig,
    outputs: dict[str, Path],
) -> None:
    """Write the unified peak profile and separate directional regressions."""

    if not config.nrl:
        return
    prefix = make_output_prefix(config)
    profile_title = f"{prefix}: aggregate repeat peaks"
    x_label = f"{config.axis_label} (bp)"

    peak_by_index = {peak.index: number for number, peak in enumerate(result.peaks, start=1)}
    detection_indices = {
        peak.detection_index for peak in result.peaks if peak.detection_index is not None
    }
    positive_indices = {peak.index for peak in result.positive_peaks}
    negative_indices = {peak.index for peak in result.negative_peaks}
    central_index = result.central_peak.index if result.central_peak is not None else None
    excluded_indices = {
        peak.index
        for peak in result.peaks
        if result.exclusion_start is not None
        and result.exclusion_end is not None
        and result.exclusion_start <= peak.distance <= result.exclusion_end
    }
    exclusion_start_text = (
        "" if result.exclusion_start is None else _format_nrl_value(result.exclusion_start)
    )
    exclusion_end_text = (
        "" if result.exclusion_end is None else _format_nrl_value(result.exclusion_end)
    )

    profile_path = outputs["nrl_profile"]
    profile_path.parent.mkdir(parents=True, exist_ok=True)
    with profile_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "relative_position",
                "unsmoothed_value",
                "local_max_smoothed_value",
                "detection_smoothed_value",
                "detection_smoothing_window",
                "local_max_smoothing_window",
                "regression_exclusion_start_bp",
                "regression_exclusion_end_bp",
                "is_detection_peak",
                "is_peak",
                "peak_number",
                "is_shared_central_peak",
                "excluded_by_regression_zone",
                "included_in_positive_regression",
                "included_in_negative_regression",
                "plot_title",
                "x_label",
                "y_label",
            ]
        )
        for index, (position, raw, local, detection) in enumerate(
            zip(
                result.positions,
                result.raw_values,
                result.local_values,
                result.detection_values,
            )
        ):
            writer.writerow(
                [
                    _format_nrl_value(position),
                    _format_nrl_value(raw),
                    _format_nrl_value(local),
                    _format_nrl_value(detection),
                    result.detection_window,
                    result.local_max_window,
                    exclusion_start_text,
                    exclusion_end_text,
                    1 if index in detection_indices else 0,
                    1 if index in peak_by_index else 0,
                    peak_by_index.get(index, ""),
                    1 if index == central_index else 0,
                    1 if index in excluded_indices else 0,
                    1 if index in positive_indices else 0,
                    1 if index in negative_indices else 0,
                    profile_title,
                    x_label,
                    config.mean_ylabel,
                ]
            )

    signed_by_index = {peak.index: peak.distance for peak in result.peaks}
    positive_order = dict(
        zip((peak.index for peak in result.positive_peaks), result.positive_peak_numbers)
    )
    negative_order = dict(
        zip((peak.index for peak in result.negative_peaks), result.negative_peak_numbers)
    )
    with outputs["nrl_peaks"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "peak_number_full_profile",
                "relative_position",
                "distance_from_zero_bp",
                "direction",
                "direction_peak_number",
                "included_in_regression",
                "included_in_positive_regression",
                "included_in_negative_regression",
                "is_shared_central_peak",
                "excluded_by_regression_zone",
                "unsmoothed_value",
                "local_max_smoothed_value",
                "detection_peak_position",
                "detection_smoothed_value",
            ]
        )
        for number, peak in enumerate(result.peaks, start=1):
            included_positive = peak.index in positive_order
            included_negative = peak.index in negative_order
            if peak.index == central_index:
                direction = "shared_centre"
                direction_number = (
                    positive_order.get(peak.index, negative_order.get(peak.index, ""))
                )
            elif peak.distance > 0:
                direction = "positive"
                direction_number = positive_order.get(peak.index, "")
            elif peak.distance < 0:
                direction = "negative"
                direction_number = negative_order.get(peak.index, "")
            else:
                direction = "centre"
                direction_number = ""
            writer.writerow(
                [
                    number,
                    _format_nrl_value(peak.distance),
                    _format_nrl_value(abs(peak.distance)),
                    direction,
                    direction_number,
                    1 if included_positive or included_negative else 0,
                    1 if included_positive else 0,
                    1 if included_negative else 0,
                    1 if peak.index == central_index else 0,
                    1 if peak.index in excluded_indices else 0,
                    _format_nrl_value(peak.raw_value),
                    _format_nrl_value(peak.smoothed_value),
                    _format_nrl_value(peak.detection_distance),
                    _format_nrl_value(peak.detection_value),
                ]
            )

    direction_data = (
        (
            "positive",
            result.positive_peaks,
            result.positive_peak_numbers,
            result.positive_regression,
            outputs["nrl_positive_regression"],
            outputs["nrl_positive_regression_plot"],
        ),
        (
            "negative",
            result.negative_peaks,
            result.negative_peak_numbers,
            result.negative_regression,
            outputs["nrl_negative_regression"],
            outputs["nrl_negative_regression_plot"],
        ),
    )
    from nucleosuite.nrl import create_regression_plot

    for (
        direction,
        peaks,
        peak_numbers,
        regression,
        table_path,
        plot_path_value,
    ) in direction_data:
        title = f"{prefix}: {direction}-direction repeat-length regression"
        _write_direction_regression(
            table_path,
            direction=direction,
            peaks=peaks,
            peak_numbers=peak_numbers,
            signed_by_index=signed_by_index,
            regression=regression,
            plot_title=title,
        )
        outputs[f"nrl_{direction}_regression_plot"] = create_regression_plot(
            plot_path_value,
            peaks,
            regression,
            title,
            config.dpi,
            peak_numbers=peak_numbers,
            y_label="Distance from position 0 (bp)",
            slope_label=f"{direction.title()} repeat length",
        )

    with outputs["nrl_summary"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "direction",
                "peak_calling_scope",
                "regression_min_distance_bp",
                "regression_max_distance_bp",
                "peak_resolution_bp",
                "detection_smoothing_window",
                "local_max_smoothing_window",
                "central_peak_position_bp",
                "central_peak_max_abs_position_bp",
                "regression_exclusion_start_bp",
                "regression_exclusion_end_bp",
                "peak_count",
                "repeat_length_bp",
                "intercept_bp",
                "r_squared",
                "slope_standard_error",
                "mean_adjacent_peak_spacing_bp",
                "quality_status",
            ]
        )
        for (
            direction,
            _peaks,
            _peak_numbers,
            regression,
            _table,
            _plot,
        ) in direction_data:
            writer.writerow(
                [
                    direction,
                    "complete_aggregate_alignment",
                    _format_nrl_value(result.regression_min),
                    _format_nrl_value(result.regression_max),
                    _format_nrl_value(result.peak_resolution),
                    result.detection_window,
                    result.local_max_window,
                    _format_nrl_value(
                        result.central_peak.distance
                        if result.central_peak is not None
                        else float("nan")
                    ),
                    _format_nrl_value(
                        result.peak_resolution / 2.0
                        if result.peak_resolution > 0
                        else 0.0
                    ),
                    exclusion_start_text or "NaN",
                    exclusion_end_text or "NaN",
                    regression.n_peaks,
                    _format_nrl_value(regression.slope),
                    _format_nrl_value(regression.intercept),
                    _format_nrl_value(regression.r_squared),
                    _format_nrl_value(regression.slope_standard_error),
                    _format_nrl_value(regression.mean_adjacent_spacing),
                    _aggregate_nrl_quality(regression),
                ]
            )

    outputs["nrl_profile_plot"] = _create_aggregate_nrl_profile_plot(
        outputs["nrl_profile_plot"],
        result,
        title=profile_title,
        x_label=x_label,
        y_label=config.mean_ylabel,
        dpi=config.dpi,
    )


def write_summary(
    path: Path,
    config: AlignmentConfig,
    stats: AlignmentStats,
    outputs: dict[str, Path],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8") as output:
        output.write("section\tkey\tvalue\n")
        for key, value in asdict(stats).items():
            output.write(f"statistics\t{key}\t{value}\n")
        for key, value in asdict(config).items():
            output.write(f"parameters\t{key}\t{'' if value is None else value}\n")
        for key, value in outputs.items():
            output.write(f"outputs\t{key}\t{value}\n")


def plot_outputs(
    matrix: np.ndarray,
    x_values: np.ndarray,
    config: AlignmentConfig,
    paths: dict[str, Path],
) -> np.ndarray:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    automatic = float(np.nanmax(np.abs(matrix))) if np.isfinite(matrix).any() else 1.0
    automatic = automatic or 1.0
    vmin = -automatic if config.vmin is None else config.vmin
    vmax = automatic if config.vmax is None else config.vmax
    if vmin >= vmax:
        raise ValueError("vmin must be less than vmax")

    paths["heatmap"].parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(15, 5))
    image = axis.imshow(
        matrix,
        aspect="auto",
        cmap="seismic",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        origin="upper",
        extent=[float(x_values[0]) - 0.5, float(x_values[-1]) + 0.5, matrix.shape[0], 0],
    )
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label(config.colorbar_label)
    from nucleosuite.plotting import apply_base_pair_x_axis
    apply_base_pair_x_axis(axis, x_values)
    axis.set_xlabel(f"{config.axis_label} (bp)")
    axis.set_ylabel("Region index (sorted)")
    from nucleosuite.plotting import save_figure
    figure.tight_layout()
    paths["heatmap"] = save_figure(figure, paths["heatmap"], default_dpi=config.dpi)
    plt.close(figure)

    finite_counts = np.sum(np.isfinite(matrix), axis=0)
    mean_profile = np.divide(
        np.nansum(matrix, axis=0),
        finite_counts,
        out=np.full(matrix.shape[1], np.nan, dtype=float),
        where=finite_counts > 0,
    )
    figure, axis = plt.subplots(figsize=(15, 5))
    axis.plot(x_values, mean_profile, color="black", linewidth=1.5, marker="o", markersize=1.6, markeredgewidth=0)
    from nucleosuite.plotting import apply_base_pair_x_axis
    apply_base_pair_x_axis(axis, x_values)
    axis.set_xlabel(f"{config.axis_label} (bp)")
    axis.set_ylabel(config.mean_ylabel)
    axis.set_xlim(float(x_values[0]), float(x_values[-1]))
    if config.mean_ylim is not None:
        axis.set_ylim(-config.mean_ylim, config.mean_ylim)
    figure.tight_layout()
    paths["mean_plot"] = save_figure(figure, paths["mean_plot"], default_dpi=config.dpi)
    plt.close(figure)
    return mean_profile


def run_alignment(
    config: AlignmentConfig,
    *,
    progress: ProgressReporter | None = None,
) -> dict[str, Path]:
    """Execute an alignment run and return all generated output paths."""
    validate_config(config)
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("pyBigWig is required to run 'nucleosuite aggregate'") from exc

    if progress is not None:
        progress.stage("Opening signal track and optional annotations")
    outputs = resolve_output_paths(config)
    stats = AlignmentStats()
    sample_rng = np.random.default_rng(config.seed)
    orientation_rng = np.random.default_rng(config.seed)
    nucleosomes = (
        load_nucleosome_centers(config.nucleosome_bed)
        if config.nucleosome_bed is not None else None
    )
    states = load_bed_intervals(config.state_bed) if config.state_bed is not None else None
    running_sum: np.ndarray | None = None
    running_count: np.ndarray | None = None
    selected_rows: list[np.ndarray] = []

    with pyBigWig.open(str(config.bigwig)) as bigwig:
        chromosome_lengths = bigwig.chroms()
        if not chromosome_lengths:
            raise RuntimeError("The BigWig contains no chromosome entries")
        from .core.blacklist import load_blacklist
        blacklist = load_blacklist(
            config.blacklist_bed,
            list(chromosome_lengths),
            [int(chromosome_lengths[name]) for name in chromosome_lengths],
        )
        seen_contigs: set[str] = set()
        if progress is not None:
            progress.file_start("regions", config.region_bed)
        with open_text(config.region_bed) as handle:
            for line_number, raw in enumerate(handle, start=1):
                stats.input_lines += 1
                if config.skip_header and line_number == 1:
                    stats.skipped_header += 1
                    continue
                stripped = raw.strip()
                if not stripped:
                    stats.skipped_blank += 1
                    continue
                if stripped.startswith("#"):
                    stats.skipped_comment += 1
                    continue
                values = stripped.split()
                chromosome = get_field(values, config.chrom_col)
                if not chromosome:
                    stats.skipped_missing_chromosome += 1
                    continue
                if progress is not None and chromosome not in seen_contigs:
                    seen_contigs.add(chromosome)
                    progress.reading_contig("regions", chromosome)
                try:
                    bigwig_chromosome = resolve_contig_name(
                        chromosome, list(chromosome_lengths), source_label="BigWig"
                    )
                except KeyError:
                    stats.skipped_chromosome_not_in_bigwig += 1
                    continue
                chromosome_length = chromosome_lengths[bigwig_chromosome]

                start = parse_integer(get_field(values, config.start_col))
                end = parse_integer(get_field(values, config.end_col))
                if start is None or end is None:
                    stats.skipped_invalid_coordinates += 1
                    continue
                if start < 0 or end < start:
                    stats.skipped_invalid_interval += 1
                    continue
                if config.point_col > 0:
                    region_center = parse_integer(get_field(values, config.point_col))
                    if region_center is None:
                        stats.skipped_invalid_point += 1
                        continue
                else:
                    region_center = (start + end) // 2

                if blacklist is not None and blacklist.overlaps(
                    bigwig_chromosome, start, max(end, start + 1)
                ):
                    stats.skipped_blacklisted_anchor += 1
                    continue

                if states is not None:
                    if end <= start:
                        stats.skipped_invalid_interval_for_state_filter += 1
                        continue
                    try:
                        state_chromosome = resolve_contig_name(
                            chromosome, list(states), source_label="state BED"
                        )
                    except KeyError:
                        stats.skipped_no_state_overlap += 1
                        continue
                    if not overlaps_any(states.get(state_chromosome), start, end):
                        stats.skipped_no_state_overlap += 1
                        continue

                strand = get_field(values, config.strand_col)
                if strand not in {"+", "-"}:
                    if config.missing_strand == "error":
                        stats.skipped_missing_or_invalid_strand += 1
                        continue
                    if config.missing_strand == "random":
                        strand = "+" if int(orientation_rng.integers(0, 2)) == 0 else "-"
                        stats.randomly_oriented += 1
                    else:
                        strand = "+"
                        stats.unstranded_kept_forward += 1

                alignment_center = region_center
                if nucleosomes is not None:
                    try:
                        nucleosome_chromosome = resolve_contig_name(
                            chromosome, list(nucleosomes), source_label="nucleosome BED"
                        )
                    except KeyError:
                        stats.skipped_missing_requested_nucleosome += 1
                        continue
                    found = find_relative_nucleosome_center(
                        nucleosomes.get(nucleosome_chromosome),
                        region_center,
                        strand,
                        config.nucleosome_offset,
                    )
                    if found is None:
                        stats.skipped_missing_requested_nucleosome += 1
                        continue
                    alignment_center = found
                    stats.aligned_to_nucleosome += 1
                else:
                    stats.aligned_to_region += 1

                try:
                    row = extract_window(
                        bigwig, bigwig_chromosome, alignment_center, config.window_half, chromosome_length
                    )
                except RuntimeError:
                    stats.skipped_bigwig_error += 1
                    continue
                blacklist_mask = np.zeros(row.size, dtype=bool)
                if blacklist is not None:
                    window_start = alignment_center - config.window_half
                    before = row.copy()
                    stats.blacklisted_window_bases += blacklist.mask_values(
                        bigwig_chromosome, window_start, row
                    )
                    blacklist_mask = np.isnan(row) & ~np.isnan(before)
                    # Also distinguish blacklist positions that were already NaN
                    # in the BigWig from ordinary missing signal.
                    blacklist_mask |= ~blacklist.valid_mask(
                        bigwig_chromosome, window_start, window_start + row.size
                    )
                if strand == "-":
                    row = row[::-1]
                    blacklist_mask = blacklist_mask[::-1]
                if config.nan_to_zero:
                    ordinary_nan = np.isnan(row) & ~blacklist_mask
                    if ordinary_nan.any():
                        stats.nan_converted_vectors += 1
                    row[ordinary_nan] = 0.0
                elif np.any(np.isnan(row) & ~blacklist_mask):
                    stats.skipped_nan += 1
                    continue
                if np.any(~np.isfinite(row) & ~blacklist_mask):
                    stats.skipped_nonfinite += 1
                    continue
                if has_consecutive_zeros(row, config.zero_thresh):
                    stats.skipped_consecutive_zeros += 1
                    continue
                if config.max_score is not None and np.any(row > config.max_score):
                    stats.skipped_above_max_score += 1
                    continue

                if running_sum is None:
                    running_sum = np.zeros_like(row)
                    running_count = np.zeros_like(row, dtype=np.int64)
                finite = np.isfinite(row)
                running_sum[finite] += row[finite]
                assert running_count is not None
                running_count[finite] += 1
                stats.valid_total += 1
                add_heatmap_row(
                    selected_rows, row, stats.valid_total,
                    config.max_heatmap_rows, config.subsample_mode, sample_rng
                )
                if config.stop_after_valid is not None and stats.valid_total >= config.stop_after_valid:
                    stats.stopped_after_valid_limit = 1
                    break

    if running_sum is None or stats.valid_total == 0:
        write_summary(outputs["summary"], config, stats, outputs)
        raise RuntimeError(no_valid_regions_message(stats, config))
    if not selected_rows:
        raise RuntimeError("No vectors were retained for plotting")

    assert running_count is not None
    if progress is not None:
        progress.stage(
            f"Creating aggregate and heatmap from {stats.valid_total:,} valid regions"
        )
    full_mean = np.divide(
        running_sum,
        running_count,
        out=np.full_like(running_sum, np.nan),
        where=running_count > 0,
    )
    matrix = np.vstack(selected_rows)
    stats.selected_for_plot = matrix.shape[0]
    write_profile(full_mean, outputs["aggregate"])
    if config.nrl:
        if progress is not None:
            progress.stage(
                "Calling aggregate peaks across the complete alignment and fitting "
                "positive/negative repeat lengths"
            )
        exclusion_start, exclusion_end = resolve_nrl_exclusion(config)
        nrl_result = analyse_aggregate_nrl(
            full_mean,
            peak_resolution=config.nrl_peak_resolution,
            regression_min=config.nrl_regression_min,
            regression_max=config.nrl_regression_max,
            exclusion_start=exclusion_start,
            exclusion_end=exclusion_end,
        )
        write_aggregate_nrl_outputs(nrl_result, config, outputs)
    matrix, x_values = central_crop(matrix, config.breadth)
    matrix, original_order, sort_scores = sort_matrix(matrix, config.sort_mode)
    write_heatmap_matrix(matrix, x_values, outputs["heatmap_matrix"])
    write_heatmap_row_metadata(
        outputs["row_metadata"], original_order, sort_scores, config.sort_mode
    )
    plotted_mean = plot_outputs(matrix, x_values, config, outputs)
    write_profile(plotted_mean, outputs["plotted_mean"], x_values)
    write_summary(outputs["summary"], config, stats, outputs)

    LOGGER.info("Valid vectors: %s", f"{stats.valid_total:,}")
    LOGGER.info("Heatmap rows: %s", f"{stats.selected_for_plot:,}")
    for name, path in outputs.items():
        LOGGER.info("%s: %s", name.replace("_", " ").title(), path)
    if progress is not None:
        progress.stage(f"Wrote {len(outputs):,} aggregate output files")
    return outputs
