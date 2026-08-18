#!/usr/bin/env python3
"""Estimate nucleosome repeat length or other periodicities from DAC/DCC tables."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.progress import ProgressReporter


@dataclass(frozen=True)
class Peak:
    """One called local maximum in a derived distance profile."""

    index: int
    distance: float
    raw_value: float
    smoothed_value: float
    rise_bp: float = float("nan")
    fall_bp: float = float("nan")
    detection_index: int | None = None
    detection_distance: float = float("nan")
    detection_value: float = float("nan")


@dataclass(frozen=True)
class Regression:
    """Linear regression of peak number against peak distance."""

    n_peaks: int
    slope: float
    intercept: float
    r_squared: float
    slope_standard_error: float
    mean_adjacent_spacing: float


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite nrl",
        description=(
            "Estimate nucleosome repeat length or another periodicity from the "
            "normalized value column of a DAC or DCC TSV."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("input", help="DAC or DCC TSV file.")
    parser.add_argument(
        "--distance-column",
        default=None,
        help="Distance/lag column. Auto-detected as Distance or Lag when omitted.",
    )
    parser.add_argument(
        "--value-column",
        default=None,
        help="Signal column. Auto-detected as DAC Value or DCC Value when omitted.",
    )
    parser.add_argument(
        "--min-distance",
        type=float,
        default=200.0,
        help="Inclusive lower distance or lag bound.",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=1200.0,
        help="Inclusive upper distance or lag bound.",
    )
    parser.add_argument(
        "--peak-resolution",
        type=float,
        default=160.0,
        help=(
            "Minimum distance between long-range peaks in bp. This also sets the "
            "peak-detection and local-maximum smoothing windows: resolution/3 and "
            "resolution/6, each rounded down to 10n+1 bp windows. Values below "
            "11 bp use no smoothing. Default: 160 bp (51 bp detection, 21 bp refinement)."
        ),
    )
    parser.add_argument(
        "--x-major-tick",
        type=float,
        default=None,
        help="Major x-axis tick interval in bp for the profile plot; default is automatic.",
    )
    parser.add_argument(
        "--x-minor-tick",
        type=float,
        default=None,
        help=(
            "Minor x-axis tick interval in bp. By default, 10 bp major ticks use 5 bp "
            "minor ticks, major intervals over 50 bp use 10 bp minor ticks, and other "
            "major intervals are divided in half."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix. Default: INPUT stem plus analysis range.",
    )
    parser.add_argument("--title", default=None, help="Optional plot title prefix.")
    parser.add_argument("--dpi", type=int, default=300, help="Figure resolution setting; --plot-dpi takes precedence when supplied.")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages.")
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser, label_points_default="peaks")
    return parser


def _pick_column(fieldnames: Sequence[str], requested: str | None, candidates: Sequence[str], kind: str) -> str:
    if requested is not None:
        if requested not in fieldnames:
            raise ValueError(
                f"Requested {kind} column {requested!r} was not found. "
                f"Available columns: {', '.join(fieldnames)}"
            )
        return requested

    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        f"Could not auto-detect the {kind} column. Expected one of: "
        + ", ".join(candidates)
    )


def read_profile(
    path: str | Path,
    distance_column: str | None,
    value_column: str | None,
    min_distance: float,
    max_distance: float,
) -> tuple[np.ndarray, np.ndarray, str, str]:
    """Read and range-filter a DAC or DCC profile."""

    input_path = Path(path)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {input_path}")

    distances: list[float] = []
    values: list[float] = []

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Input TSV has no header: {input_path}")
        x_column = _pick_column(
            reader.fieldnames,
            distance_column,
            ("Distance", "Lag"),
            "distance",
        )
        y_column = _pick_column(
            reader.fieldnames,
            value_column,
            ("DAC Value", "DCC Value"),
            "value",
        )

        for line_number, row in enumerate(reader, start=2):
            try:
                x = float(row[x_column])
                y = float(row[y_column])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{input_path}: line {line_number} contains a non-numeric "
                    f"{x_column!r} or {y_column!r} value."
                ) from exc
            if not (math.isfinite(x) and math.isfinite(y)):
                continue
            if min_distance <= x <= max_distance:
                distances.append(x)
                values.append(y)

    if len(distances) < 3:
        raise ValueError(
            f"Fewer than three finite rows fall within {min_distance:g} to "
            f"{max_distance:g} in {input_path}."
        )

    x_array = np.asarray(distances, dtype=np.float64)
    y_array = np.asarray(values, dtype=np.float64)
    order = np.argsort(x_array, kind="stable")
    x_array = x_array[order]
    y_array = y_array[order]

    if np.any(np.diff(x_array) <= 0):
        raise ValueError(
            "Distance/lag values must be unique and strictly increasing after sorting."
        )

    return x_array, y_array, x_column, y_column


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    """Return a centred moving average with partial edge windows."""

    if window < 1:
        raise ValueError("--smooth-window must be at least 1.")
    if window == 1:
        return values.astype(np.float64, copy=True)
    if window > values.size:
        raise ValueError(
            f"--smooth-window ({window}) exceeds the selected profile length "
            f"({values.size})."
        )

    kernel = np.ones(window, dtype=np.float64)
    numerator = np.convolve(values, kernel, mode="same")
    denominator = np.convolve(np.ones(values.size, dtype=np.float64), kernel, mode="same")
    return numerator / denominator


def moving_average_by_distance(
    distances: np.ndarray, values: np.ndarray, window_bp: int
) -> np.ndarray:
    """Return a centred moving average over a genomic-distance window in bp."""

    distances = np.asarray(distances, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if distances.shape != values.shape:
        raise ValueError("Distance and value arrays must have the same shape.")
    if window_bp < 1:
        raise ValueError("Smoothing window must be at least 1 bp.")
    if window_bp == 1:
        return values.astype(np.float64, copy=True)
    if window_bp % 2 == 0:
        raise ValueError("Smoothing window must be odd.")

    half_width = (float(window_bp) - 1.0) / 2.0
    cumulative = np.concatenate(([0.0], np.cumsum(values, dtype=np.float64)))
    output = np.empty(values.size, dtype=np.float64)
    left = 0
    right = 0
    n = values.size
    for i, centre in enumerate(distances):
        lower = centre - half_width
        upper = centre + half_width
        while left < n and distances[left] < lower:
            left += 1
        if right < left:
            right = left
        while right < n and distances[right] <= upper:
            right += 1
        count = right - left
        output[i] = (cumulative[right] - cumulative[left]) / count if count else values[i]
    return output


def snap_smoothing_window(target_window: float) -> int:
    """Round down to 11, 21, 31, ...; return 1 for no smoothing below 11."""

    if target_window < 0:
        raise ValueError("Smoothing-window target must be 0 or greater.")
    if target_window < 11:
        return 1
    return int(math.floor((float(target_window) - 1.0) / 10.0) * 10 + 1)


def resolution_smoothing_windows(peak_resolution: float) -> tuple[int, int]:
    """Derive detection and local-maximum smoothing windows from peak resolution."""

    if peak_resolution < 0:
        raise ValueError("--peak-resolution must be 0 or greater.")
    return (
        snap_smoothing_window(float(peak_resolution) / 3.0),
        snap_smoothing_window(float(peak_resolution) / 6.0),
    )


def _peak_from_refined_index(
    refined_index: int,
    detection_index: int,
    distances: np.ndarray,
    raw_values: np.ndarray,
    local_values: np.ndarray,
    detection_values: np.ndarray,
) -> Peak:
    return Peak(
        index=int(refined_index),
        distance=float(distances[refined_index]),
        raw_value=float(raw_values[refined_index]),
        smoothed_value=float(local_values[refined_index]),
        detection_index=int(detection_index),
        detection_distance=float(distances[detection_index]),
        detection_value=float(detection_values[detection_index]),
    )


def call_resolution_peaks(
    distances: np.ndarray,
    raw_values: np.ndarray,
    local_values: np.ndarray,
    detection_values: np.ndarray,
    peak_resolution: float,
) -> list[Peak]:
    """Detect broad peaks, then refine each one on the finer-smoothed profile."""

    if peak_resolution < 0:
        raise ValueError("--peak-resolution must be 0 or greater.")

    candidate_indices = local_maximum_indices(
        distances, detection_values, min_rise_bp=0.0, min_fall_bp=0.0
    )
    detected = retain_separated_peaks(
        candidate_indices, distances, raw_values, detection_values, peak_resolution
    )
    detection_indices = [peak.index for peak in detected]
    if not detection_indices:
        return []

    fine_candidates = local_maximum_indices(
        distances, local_values, min_rise_bp=0.0, min_fall_bp=0.0
    )
    half_resolution = float(peak_resolution) / 2.0
    refined: list[Peak] = []

    for detection_index in detection_indices:
        if peak_resolution == 0:
            eligible = [idx for idx in fine_candidates if idx == detection_index]
        else:
            eligible = [
                idx for idx in fine_candidates
                if abs(float(distances[idx] - distances[detection_index])) <= half_resolution
            ]

        if eligible:
            refined_index = min(
                eligible,
                key=lambda idx: (
                    -float(local_values[idx]),
                    abs(float(distances[idx] - distances[detection_index])),
                    -float(raw_values[idx]),
                    float(distances[idx]),
                ),
            )
        else:
            if peak_resolution == 0:
                neighbourhood = [detection_index]
            else:
                neighbourhood = np.flatnonzero(
                    np.abs(distances - distances[detection_index]) <= half_resolution
                ).tolist() or [detection_index]
            refined_index = min(
                neighbourhood,
                key=lambda idx: (
                    -float(local_values[idx]),
                    abs(float(distances[idx] - distances[detection_index])),
                    -float(raw_values[idx]),
                    float(distances[idx]),
                ),
            )

        refined.append(
            _peak_from_refined_index(
                refined_index, detection_index, distances, raw_values, local_values, detection_values
            )
        )

    if peak_resolution == 0:
        return sorted(refined, key=lambda peak: peak.distance)

    ranked = sorted(
        refined,
        key=lambda peak: (
            -peak.detection_value, -peak.smoothed_value, -peak.raw_value, peak.distance
        ),
    )
    kept: list[Peak] = []
    for peak in ranked:
        if all(abs(peak.distance - other.distance) >= peak_resolution for other in kept):
            kept.append(peak)
    return sorted(kept, key=lambda peak: peak.distance)


def _plateau_bounds(values: np.ndarray, index: int) -> tuple[int, int]:
    """Return inclusive bounds of the equal-valued plateau containing ``index``."""

    start = int(index)
    end = int(index)
    while start > 0 and values[start - 1] == values[index]:
        start -= 1
    last = values.size - 1
    while end < last and values[end + 1] == values[index]:
        end += 1
    return start, end


def peak_shape_spans(
    distances: np.ndarray,
    values: np.ndarray,
    index: int,
) -> tuple[float, float]:
    """Return strict contiguous rise and fall spans around a summit plateau.

    Equal values are allowed within the summit plateau but break monotonic rise
    and fall runs elsewhere. Distances, rather than row counts, are used so the
    rule remains meaningful for profiles sampled at intervals other than 1 bp.
    """

    plateau_start, plateau_end = _plateau_bounds(values, index)

    left = plateau_start
    while left > 0 and values[left] > values[left - 1]:
        left -= 1
    rise_bp = float(distances[plateau_start] - distances[left])

    right = plateau_end
    last = values.size - 1
    while right < last and values[right] > values[right + 1]:
        right += 1
    fall_bp = float(distances[right] - distances[plateau_end])
    return rise_bp, fall_bp


def local_maximum_indices(
    distances: np.ndarray,
    values: np.ndarray | None = None,
    min_rise_bp: float = 2.0,
    min_fall_bp: float = 2.0,
) -> list[int]:
    """Call plateau-aware maxima with persistent strict rise and fall.

    ``distances`` and ``values`` normally contain the profile coordinates and
    analysis values. For compatibility with direct Python use, passing only one
    array treats it as ``values`` sampled at 1 bp intervals.

    A summit is retained only when the immediately adjacent monotonic run rises
    for at least ``min_rise_bp`` and falls for at least ``min_fall_bp``. A flat
    summit is represented by its midpoint. Equal values outside that summit
    plateau break the rise or fall, preventing momentary one-position blips from
    being called.
    """

    if values is None:
        values = np.asarray(distances, dtype=np.float64)
        distances = np.arange(values.size, dtype=np.float64)
    else:
        distances = np.asarray(distances, dtype=np.float64)
        values = np.asarray(values, dtype=np.float64)

    if distances.shape != values.shape:
        raise ValueError("Distance and value arrays must have the same shape.")
    if min_rise_bp < 0 or min_fall_bp < 0:
        raise ValueError("Minimum rise and fall spans must be 0 or greater.")

    candidates: list[int] = []
    index = 0
    last_index = values.size - 1

    while index <= last_index:
        plateau_start = index
        plateau_end = index
        while plateau_end < last_index and values[plateau_end + 1] == values[index]:
            plateau_end += 1

        if plateau_start > 0 and plateau_end < last_index:
            is_local_maximum = (
                values[plateau_start - 1] < values[plateau_start]
                and values[plateau_end + 1] < values[plateau_end]
            )
            if is_local_maximum:
                centre = (plateau_start + plateau_end) // 2
                rise_bp, fall_bp = peak_shape_spans(distances, values, centre)
                if rise_bp >= min_rise_bp and fall_bp >= min_fall_bp:
                    candidates.append(centre)

        index = plateau_end + 1

    return candidates


def retain_separated_peaks(
    candidate_indices: Sequence[int],
    distances: np.ndarray,
    raw_values: np.ndarray,
    smoothed_values: np.ndarray,
    min_separation: float,
) -> list[Peak]:
    """Keep the highest local maximum in each competing distance neighbourhood."""

    if min_separation < 0:
        raise ValueError("--min-peak-separation must be 0 or greater.")

    if min_separation == 0:
        retained_indices = list(candidate_indices)
    else:
        ranked = sorted(
            candidate_indices,
            key=lambda idx: (
                -float(smoothed_values[idx]),
                -float(raw_values[idx]),
                float(distances[idx]),
            ),
        )
        retained_indices: list[int] = []
        for idx in ranked:
            if all(
                abs(float(distances[idx] - distances[kept])) >= min_separation
                for kept in retained_indices
            ):
                retained_indices.append(idx)
        retained_indices.sort(key=lambda idx: float(distances[idx]))

    return [
        Peak(
            index=idx,
            distance=float(distances[idx]),
            raw_value=float(raw_values[idx]),
            smoothed_value=float(smoothed_values[idx]),
            rise_bp=peak_shape_spans(distances, smoothed_values, idx)[0],
            fall_bp=peak_shape_spans(distances, smoothed_values, idx)[1],
        )
        for idx in retained_indices
    ]


def regress_peak_distances(
    peaks: Sequence[Peak],
    peak_numbers: Sequence[float] | None = None,
) -> Regression:
    """Regress peak number (x) against peak distance (y)."""

    n_peaks = len(peaks)
    if peak_numbers is None:
        numbers = np.arange(1, n_peaks + 1, dtype=np.float64)
    else:
        numbers = np.asarray(peak_numbers, dtype=np.float64)
        if numbers.shape != (n_peaks,):
            raise ValueError("peak_numbers must contain one value for each peak")
        if not np.isfinite(numbers).all():
            raise ValueError("peak_numbers must be finite")
        if numbers.size > 1 and np.any(np.diff(numbers) <= 0):
            raise ValueError("peak_numbers must be strictly increasing")
    if n_peaks < 2:
        return Regression(
            n_peaks=n_peaks,
            slope=float("nan"),
            intercept=float("nan"),
            r_squared=float("nan"),
            slope_standard_error=float("nan"),
            mean_adjacent_spacing=float("nan"),
        )

    peak_distances = np.asarray([peak.distance for peak in peaks], dtype=np.float64)
    slope, intercept = np.polyfit(numbers, peak_distances, 1)
    fitted = intercept + slope * numbers
    residuals = peak_distances - fitted
    total_sum_squares = float(np.sum((peak_distances - np.mean(peak_distances)) ** 2))
    residual_sum_squares = float(np.sum(residuals**2))
    r_squared = (
        1.0 - residual_sum_squares / total_sum_squares
        if total_sum_squares > 0
        else 1.0
    )

    if n_peaks > 2:
        x_ss = float(np.sum((numbers - np.mean(numbers)) ** 2))
        residual_variance = residual_sum_squares / float(n_peaks - 2)
        slope_standard_error = math.sqrt(residual_variance / x_ss) if x_ss > 0 else float("nan")
    else:
        slope_standard_error = float("nan")

    mean_spacing = float(np.mean(np.diff(peak_distances) / np.diff(numbers)))
    return Regression(
        n_peaks=n_peaks,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r_squared),
        slope_standard_error=float(slope_standard_error),
        mean_adjacent_spacing=mean_spacing,
    )


def _format_float(value: float) -> str:
    return "NaN" if not math.isfinite(value) else f"{value:.12g}"


def write_profile_tsv(
    path: Path,
    distances: np.ndarray,
    raw_values: np.ndarray,
    local_values: np.ndarray,
    detection_values: np.ndarray,
    peaks: Sequence[Peak],
) -> None:
    peak_number_by_index = {peak.index: number for number, peak in enumerate(peaks, start=1)}
    detection_indices = {peak.detection_index for peak in peaks if peak.detection_index is not None}
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "distance_or_lag", "unsmoothed_value", "local_max_smoothed_value",
            "detection_smoothed_value", "is_detection_peak", "is_peak", "peak_number"
        ])
        for idx, (distance, raw, local, detection) in enumerate(
            zip(distances, raw_values, local_values, detection_values)
        ):
            peak_number = peak_number_by_index.get(idx)
            writer.writerow([
                f"{float(distance):.12g}", f"{float(raw):.12g}", f"{float(local):.12g}",
                f"{float(detection):.12g}", 1 if idx in detection_indices else 0,
                1 if peak_number is not None else 0, peak_number if peak_number is not None else ""
            ])

def write_peaks_tsv(path: Path, peaks: Sequence[Peak]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "peak_number", "distance_or_lag", "unsmoothed_value",
            "local_max_smoothed_value", "detection_peak_distance", "detection_smoothed_value"
        ])
        for number, peak in enumerate(peaks, start=1):
            writer.writerow([
                number, f"{peak.distance:.12g}", f"{peak.raw_value:.12g}",
                f"{peak.smoothed_value:.12g}", _format_float(peak.detection_distance),
                _format_float(peak.detection_value)
            ])

def write_regression_tsv(
    path: Path, input_path: Path, distance_column: str, value_column: str,
    min_distance: float, max_distance: float, peak_resolution: float,
    detection_window: int, local_max_window: int, regression: Regression,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "input", "distance_column", "value_column", "min_distance", "max_distance",
            "peak_resolution_bp", "detection_smoothing_window", "local_max_smoothing_window",
            "peak_count", "slope_bp_per_peak", "intercept_bp", "r_squared",
            "slope_standard_error", "mean_adjacent_peak_spacing_bp"
        ])
        writer.writerow([
            str(input_path), distance_column, value_column, f"{min_distance:.12g}", f"{max_distance:.12g}",
            f"{peak_resolution:.12g}", detection_window, local_max_window, regression.n_peaks,
            _format_float(regression.slope), _format_float(regression.intercept),
            _format_float(regression.r_squared), _format_float(regression.slope_standard_error),
            _format_float(regression.mean_adjacent_spacing)
        ])

def create_profile_plot(
    path: Path, distances: np.ndarray, raw_values: np.ndarray, local_values: np.ndarray,
    detection_values: np.ndarray, peaks: Sequence[Peak], local_max_window: int,
    detection_window: int, distance_column: str, value_column: str, title: str, dpi: int,
    x_major_tick: float | None = None, x_minor_tick: float | None = None,
) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()
    figure, axis = plt.subplots(figsize=(10, 5.5))
    axis.plot(distances, raw_values, color="0.72", linewidth=0.9, label="Unsmoothed")
    axis.plot(
        distances, local_values, color="black", linewidth=1.5,
        label=(f"Local maxima ({local_max_window} bp)" if local_max_window > 1 else "Local maxima signal"),
    )
    if detection_window != local_max_window:
        axis.plot(
            distances, detection_values, color="0.4", linewidth=1.2, linestyle="--",
            label=(f"Peak detection ({detection_window} bp)" if detection_window > 1 else "Peak detection signal"),
        )
    if peaks:
        axis.scatter(
            [peak.distance for peak in peaks], [peak.smoothed_value for peak in peaks],
            s=28, facecolors="white", edgecolors="black", linewidths=1.0, zorder=4, label="Called peaks"
        )
    axis.set_xlabel(f"{distance_column} (bp)")
    axis.set_ylabel(value_column)
    axis.set_title(title)
    axis.legend(frameon=False)
    from nucleosuite.plotting import apply_distance_x_axis
    apply_distance_x_axis(axis, major_interval=x_major_tick, minor_interval=x_minor_tick)
    from nucleosuite.plotting import annotate_points, get_plot_options, save_figure
    options = get_plot_options()
    if peaks:
        annotate_points(
            axis, [peak.distance for peak in peaks], [peak.smoothed_value for peak in peaks],
            points_are_peaks=True, options=options
        )
    figure.tight_layout()
    saved = save_figure(figure, path, default_dpi=dpi)
    plt.close(figure)
    return saved

def create_regression_plot(
    path: Path,
    peaks: Sequence[Peak],
    regression: Regression,
    title: str,
    dpi: int,
    *,
    peak_numbers: Sequence[float] | None = None,
    x_label: str = "Peak number",
    y_label: str = "Peak distance (bp)",
    slope_label: str = "Slope",
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    figure, axis = plt.subplots(figsize=(6.5, 6.5))
    if peaks:
        numbers = (
            np.arange(1, len(peaks) + 1, dtype=np.float64)
            if peak_numbers is None
            else np.asarray(peak_numbers, dtype=np.float64)
        )
        if numbers.shape != (len(peaks),):
            raise ValueError("peak_numbers must contain one value for each peak")
        peak_distances = np.asarray([peak.distance for peak in peaks], dtype=np.float64)
        axis.scatter(
            numbers,
            peak_distances,
            s=34,
            facecolors="none",
            edgecolors="black",
            linewidths=1.0,
            label="Called peaks",
        )
        if math.isfinite(regression.slope):
            fitted = regression.intercept + regression.slope * numbers
            axis.plot(
                numbers,
                fitted,
                color="black",
                linewidth=1.5,
                linestyle=":",
                label="Linear regression",
            )
            annotation = (
                f"{slope_label} = {regression.slope:.3f} bp/peak\n"
                f"$R^2$ = {regression.r_squared:.4f}"
            )
            axis.text(
                0.04,
                0.96,
                annotation,
                transform=axis.transAxes,
                ha="left",
                va="top",
            )
    else:
        axis.text(0.5, 0.5, "No peaks called", transform=axis.transAxes, ha="center", va="center")

    from nucleosuite.plotting import apply_integer_x_axis, apply_integer_y_axis
    if peak_numbers is None:
        displayed_numbers = range(1, len(peaks) + 1)
    else:
        displayed_numbers = peak_numbers
    apply_integer_x_axis(axis, displayed_numbers)
    apply_integer_y_axis(axis)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_title(title)
    if peaks:
        axis.legend(frameon=False)
    axis.grid(False)
    from nucleosuite.plotting import save_figure
    figure.tight_layout()
    saved = save_figure(figure, path, default_dpi=dpi)
    plt.close(figure)
    return saved


def default_output_prefix(input_path: Path, min_distance: float, max_distance: float) -> Path:
    def compact(value: float) -> str:
        return str(int(value)) if float(value).is_integer() else f"{value:g}".replace(".", "p")

    return input_path.with_name(
        f"{input_path.stem}_nrl_{compact(min_distance)}_{compact(max_distance)}"
    )


def run(args: argparse.Namespace) -> int:
    if args.max_distance <= args.min_distance:
        raise ValueError("--max-distance must be greater than --min-distance.")
    if args.peak_resolution < 0:
        raise ValueError("--peak-resolution must be 0 or greater.")
    from nucleosuite.plotting import validate_tick_interval
    validate_tick_interval(args.x_major_tick, "--x-major-tick")
    validate_tick_interval(args.x_minor_tick, "--x-minor-tick")
    if args.dpi < 1:
        raise ValueError("--dpi must be at least 1.")

    reporter = ProgressReporter("nrl")
    input_path = Path(args.input).resolve()
    reporter.stage(f"Loading distance profile: {input_path}")
    distances, raw_values, distance_column, value_column = read_profile(
        input_path, args.distance_column, args.value_column, args.min_distance, args.max_distance
    )
    detection_window, local_max_window = resolution_smoothing_windows(args.peak_resolution)
    reporter.stage(
        f"Calling periodic peaks (resolution {args.peak_resolution:g} bp; "
        f"detection {detection_window} bp; local maxima {local_max_window} bp)"
    )
    detection_values = moving_average_by_distance(distances, raw_values, detection_window)
    local_values = moving_average_by_distance(distances, raw_values, local_max_window)
    peaks = call_resolution_peaks(
        distances, raw_values, local_values, detection_values, args.peak_resolution
    )
    regression = regress_peak_distances(peaks)

    prefix = Path(args.output_prefix).resolve() if args.output_prefix else default_output_prefix(
        input_path, args.min_distance, args.max_distance
    )
    prefix.parent.mkdir(parents=True, exist_ok=True)
    profile_tsv = Path(f"{prefix}_profile.tsv")
    peaks_tsv = Path(f"{prefix}_peaks.tsv")
    regression_tsv = Path(f"{prefix}_regression.tsv")
    profile_png = Path(f"{prefix}_profile.png")
    regression_png = Path(f"{prefix}_regression.png")

    reporter.stage("Writing NRL tables and plots")
    write_profile_tsv(profile_tsv, distances, raw_values, local_values, detection_values, peaks)
    write_peaks_tsv(peaks_tsv, peaks)
    write_regression_tsv(
        regression_tsv, input_path, distance_column, value_column, args.min_distance, args.max_distance,
        args.peak_resolution, detection_window, local_max_window, regression
    )
    plot_title = args.title or input_path.stem
    profile_png = create_profile_plot(
        profile_png, distances, raw_values, local_values, detection_values, peaks,
        local_max_window, detection_window, distance_column, value_column,
        f"{plot_title}: {args.min_distance:g}–{args.max_distance:g} bp",
        args.dpi, args.x_major_tick, args.x_minor_tick
    )
    regression_png = create_regression_plot(
        regression_png, peaks, regression, f"{plot_title}: peak-spacing regression", args.dpi
    )

    if not args.quiet:
        print(f"Input: {input_path}")
        print(f"Range: {args.min_distance:g} to {args.max_distance:g} bp")
        print(f"Peak resolution: {args.peak_resolution:g} bp")
        print("Detection smoothing window: " + (f"{detection_window} bp" if detection_window > 1 else "none"))
        print("Local-max smoothing window: " + (f"{local_max_window} bp" if local_max_window > 1 else "none"))
        print(f"Called peaks: {len(peaks)}")
        if math.isfinite(regression.slope):
            print(f"Estimated periodicity/NRL: {regression.slope:.6g} bp")
            print(f"R-squared: {regression.r_squared:.6g}")
        else:
            print("Estimated periodicity/NRL: NaN (fewer than two peaks)")
        for output in (profile_tsv, peaks_tsv, regression_tsv, profile_png, regression_png):
            print(f"Wrote: {output}")
    return 0

def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
