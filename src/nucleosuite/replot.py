"""Recreate and customise plots from NucleoSuite tabular outputs.

The ``nucleosuite plot`` command is deliberately file-driven: it reads TSV/TSV.GZ
outputs already written by NucleoSuite, detects the plot family from the filename
and table columns, and creates a new figure without rerunning the genomic analysis.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


PLOT_TYPES = (
    "auto",
    "dac",
    "dcc",
    "nrl-profile",
    "nrl-regression",
    "fragment-size-nrl-profile",
    "fragment-size-nrl-regression",
    "aggregate-nrl-profile",
    "aggregate-nrl-regression",
    "distances",
    "distance-state-overlay",
    "distance-percentile-curves",
    "distance-percentile-peak-counts",
    "aggregate-profile",
    "heatmap",
    "fragment-lengths",
    "positive-runs",
    "peak-score-frequency",
    "peak-states",
    "flank-spacing",
    "compare-positions",
    "compare-positions-score",
    "compare-positions-histogram",
    "compare-positions-correlation",
    "compare-positions-percentile-boxplot",
    "compare-positions-score-distance",
    "gene-expression",
    "gene-expression-spacing",
    "gene-expression-spacing-scatter",
    "gene-expression-fft-trajectory",
    "gene-expression-ranking",
    "tss-expression",
    "dinucleotide-profile",
    "ww-ss-profile",
    "ww-type-summary",
    "ww-type-by-length",
    "gene-sets-venn",
    "profile-overlay",
    "count-profile",
    "generic-line",
    "generic-scatter",
    "generic-bar",
    "generic-heatmap",
)

COMMAND_TYPES = {
    "dac": "dac",
    "dcc": "dcc",
    "nrl": "nrl-profile",
    "distances": "distances",
    "aggregate": "aggregate-profile",
    "fragment-heatmap": "heatmap",
    "fragment-lengths": "fragment-lengths",
    "fragments": "fragment-lengths",
    "positive-runs": "positive-runs",
    "peak-score-frequency": "peak-score-frequency",
    "peak-states": "peak-states",
    "flank-spacing": "flank-spacing",
    "compare-positions": "compare-positions",
    "gene-expression": "gene-expression",
    "tss-expression-quintiles": "tss-expression",
    "dinuc-profile": "dinucleotide-profile",
    "ww-types": "ww-type-summary",
    "gene-sets": "gene-sets-venn",
    "randomize-fragments": "count-profile",
}


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).lower().endswith(".gz") else path.open("rt", encoding="utf-8", newline="")


def _clean_header(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def _read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No tabular header was found in {path}")
        headers = [str(h) for h in reader.fieldnames]
        rows = [{str(k): ("" if v is None else str(v)) for k, v in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"No data rows were found in {path}")
    return headers, rows


def _header_map(headers: Sequence[str]) -> dict[str, str]:
    return {_clean_header(h): h for h in headers}


def _column(rows: Sequence[Mapping[str, str]], name: str) -> list[str]:
    return [str(row.get(name, "")) for row in rows]


def _numeric(rows: Sequence[Mapping[str, str]], name: str) -> np.ndarray:
    values: list[float] = []
    for text in _column(rows, name):
        try:
            values.append(float(text))
        except (TypeError, ValueError):
            values.append(math.nan)
    return np.asarray(values, dtype=float)


def _first_numeric_columns(headers: Sequence[str], rows: Sequence[Mapping[str, str]]) -> list[str]:
    found: list[str] = []
    for header in headers:
        values = _numeric(rows, header)
        if values.size and np.isfinite(values).sum() >= max(1, int(0.6 * values.size)):
            found.append(header)
    return found


def detect_plot_type(path: Path, headers: Sequence[str]) -> str:
    names = set(_header_map(headers))
    lower_name = path.name.lower()
    if {"distance", "dac_value"}.issubset(names):
        return "dac"
    if "dcc_value" in names and ("lag" in names or "distance" in names):
        return "dcc"
    if {"distance_or_lag", "unsmoothed_value", "local_max_smoothed_value"}.issubset(names):
        return "nrl-profile"
    if {
        "relative_position",
        "unsmoothed_value",
        "local_max_smoothed_value",
        "detection_smoothed_value",
    }.issubset(names):
        return "aggregate-nrl-profile"
    if {
        "direction",
        "peak_number",
        "signed_position_bp",
        "distance_from_zero_bp",
        "fitted_distance_from_zero_bp",
    }.issubset(names):
        return "aggregate-nrl-regression"
    if {"peak_number", "fragment_length", "fitted_fragment_length"}.issubset(names):
        return "fragment-size-nrl-regression"
    if {
        "fragment_length",
        "unsmoothed_density",
        "local_max_smoothed_density",
    }.issubset(names):
        return "fragment-size-nrl-profile"
    if {"order", "peak_distance_bp", "fitted_distance_bp"}.issubset(names):
        return "nrl-regression"
    if {"peak_number", "distance_or_lag"}.issubset(names):
        return "nrl-regression"
    if "position" in names and any(re.fullmatch(r"[acgt]{2}_(?:pct|frac)", name) for name in names):
        if "ww_ss" in lower_name or lower_name.endswith("_ww_ss_profile.tsv"):
            return "ww-ss-profile"
        return "dinucleotide-profile"
    if "fragment_length" in names and any(name.endswith("_percent_of_classified") for name in names):
        return "ww-type-by-length"
    if "candidate_sets" in names and "gene_id" in names:
        return "gene-sets-venn"
    if "ww_type" in names and "count" in names:
        return "ww-type-summary"
    if {"quintile", "relative_position", "mean_signal"}.issubset(names):
        return "tss-expression"
    if {"state", "rgb", "distance_bp", "raw_percent", "smoothed_percent"}.issubset(names):
        return "distance-state-overlay"
    if {"percentile_threshold", "distance_bp", "raw_count", "raw_percent", "order", "scope"}.issubset(names):
        return "distance-percentile-curves"
    if {"retained_peak_count", "percentile_threshold", "order", "scope"}.issubset(names):
        return "distance-percentile-peak-counts"
    if {"relative_position", "score"}.issubset(names):
        return "aggregate-profile"
    if "relative_position" in names and len(headers) >= 3:
        return "profile-overlay"
    if "heatmap_matrix" in lower_name or "normalised_matrix" in lower_name or "normalized_matrix" in lower_name:
        return "heatmap"
    if {"fragment_length", "count"}.issubset(names):
        return "fragment-lengths"
    if {"run_length_bp", "count"}.issubset(names):
        return "positive-runs"
    if "dataset" in names and "count" in names and ("score" in names or "bin_midpoint" in names):
        return "peak-score-frequency"
    if {"percentile_threshold", "state", "percentage_of_assigned_peaks"}.issubset(names):
        return "peak-states"
    if {"category", "spacing_bp", "value", "distribution", "rank", "highlighted"}.issubset(names):
        return "flank-spacing"
    if {"order", "scope", "distance_bp"}.issubset(names) and ({"count", "raw_count"} & names):
        return "distances"
    if {"row_type", "percentile_group", "q1_absolute_distance", "whisker_low_absolute_distance"}.issubset(names):
        return "compare-positions-percentile-boxplot"
    if {"main_score_plot", "compare_score_plot", "absolute_distance"}.issubset(names):
        return "compare-positions-score"
    if {"main_score", "matched_distance", "distance_type"}.issubset(names):
        return "compare-positions-score-distance"
    if {"a_score", "b_score", "absolute_distance"}.issubset(names):
        if "percentile_group" in names:
            return "compare-positions-percentile-boxplot"
        return "compare-positions-score"
    if {"main_score", "compare_score", "absolute_distance"}.issubset(names):
        return "compare-positions-score"
    if {"percentile_group", "absolute_distance"}.issubset(names):
        return "compare-positions-percentile-boxplot"
    if {"bin_start_inclusive", "bin_end_exclusive", "pair_count"}.issubset(names):
        return "compare-positions-histogram"
    if {"distance_bin", "pair_count"}.issubset(names) and ({"spearman_score_correlation", "pearson_score_correlation"} & names):
        return "compare-positions-correlation"
    if {"sample", "profile", "body_median_spacing_bp", "transformed_expression"}.issubset(names):
        return "gene-expression-spacing-scatter"
    if {"sample", "region", "subset", "profile", "correlation"}.issubset(names):
        return "gene-expression-spacing"
    if {"sample", "profile", "period_bp", "correlation"}.issubset(names):
        return "gene-expression-fft-trajectory"
    if {"sample", "rank", "profile", "correlation", "ranking_periods"}.issubset(names):
        return "gene-expression-ranking"
    if {"sample", "profile", "correlation"}.issubset(names):
        return "gene-expression"
    if {"quintile", "relative_position", "mean_signal"}.issubset(names):
        return "tss-expression"
    if "fragment" in lower_name and "heatmap" in lower_name:
        return "heatmap"
    if {"relocation_bp", "count"}.issubset(names):
        return "count-profile"
    if "distance" in names and "count" in names:
        return "distances"
    return "generic-line"


def _parse_value(text: str) -> Any:
    lowered = text.strip().lower()
    if lowered in {"true", "yes", "on"}:
        return True
    if lowered in {"false", "no", "off"}:
        return False
    if lowered in {"none", "null"}:
        return None
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return text


def _parse_key_values(values: Sequence[str] | None, *, option: str) -> dict[str, Any]:
    parsed: dict[str, Any] = {}
    for item in values or []:
        if "=" not in item:
            raise ValueError(f"{option} expects KEY=VALUE, received {item!r}")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"{option} received an empty key")
        parsed[key] = _parse_value(value.strip())
    return parsed


def _parse_artist_values(values: Sequence[str] | None) -> dict[str, dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = defaultdict(dict)
    for item in values or []:
        if "=" not in item or "." not in item.split("=", 1)[0]:
            raise ValueError("--mpl-kw expects TARGET.KEY=VALUE, for example raw.color=0.7")
        lhs, value = item.split("=", 1)
        target, key = lhs.split(".", 1)
        target = target.strip().lower()
        if target not in {"line", "raw", "smooth", "points", "bar", "heatmap", "legend"}:
            raise ValueError(f"Unsupported --mpl-kw target {target!r}")
        groups[target][key.strip()] = _parse_value(value.strip())
    return dict(groups)


def _apply_rc(values: Sequence[str] | None) -> None:
    import matplotlib as mpl
    for key, value in _parse_key_values(values, option="--mpl-rc").items():
        if key not in mpl.rcParams:
            raise ValueError(f"Unknown Matplotlib rcParam: {key}")
        mpl.rcParams[key] = value


def _resolve_output(input_path: Path, output: Path | None, fmt: str) -> Path:
    if output is not None:
        path = Path(output)
        if path.suffix.lower() not in {".png", ".svg", ".pdf"}:
            path = path.with_suffix(f".{fmt}")
        return path
    name = input_path.name
    if name.lower().endswith(".tsv.gz"):
        stem = name[:-7]
    else:
        stem = input_path.stem
    return input_path.parent / f"{stem}_replot.{fmt}"


def _nice_bp_major(span: float) -> float:
    if not math.isfinite(span) or span <= 0:
        return 10.0
    raw = max(10.0, span / 8.0)
    power = 10 ** math.floor(math.log10(raw))
    for multiplier in (1, 2, 5, 10):
        candidate = multiplier * power
        if candidate >= raw:
            return float(max(10.0, candidate))
    return float(max(10.0, 10 * power))


def _configure_ticks_and_grids(ax, args, *, bp_x: bool = False, default_x_major: float | None = None, default_x_minor: float | None = None) -> None:
    from matplotlib.ticker import MultipleLocator

    x0, x1 = ax.get_xlim()
    span = abs(float(x1 - x0))
    x_major = args.x_major_tick
    x_minor = args.x_minor_tick
    if bp_x:
        if x_major is None:
            x_major = default_x_major if default_x_major is not None else _nice_bp_major(span)
        if x_minor is None:
            if default_x_minor is not None:
                x_minor = default_x_minor
            elif x_major > 50:
                x_minor = 10.0
            elif math.isclose(x_major, 10.0):
                x_minor = 5.0
            else:
                x_minor = x_major / 2.0
    if x_major is not None:
        if x_major <= 0:
            raise ValueError("--x-major-tick must be greater than zero")
        ax.xaxis.set_major_locator(MultipleLocator(float(x_major)))
    if x_minor is not None:
        if x_minor <= 0:
            raise ValueError("--x-minor-tick must be greater than zero")
        ax.xaxis.set_minor_locator(MultipleLocator(float(x_minor)))
    if args.y_major_tick is not None:
        if args.y_major_tick <= 0:
            raise ValueError("--y-major-tick must be greater than zero")
        ax.yaxis.set_major_locator(MultipleLocator(float(args.y_major_tick)))
    if args.y_minor_tick is not None:
        if args.y_minor_tick <= 0:
            raise ValueError("--y-minor-tick must be greater than zero")
        ax.yaxis.set_minor_locator(MultipleLocator(float(args.y_minor_tick)))

    # New plot command defaults to restrained vertical guides only. Major and
    # minor grid families can be enabled/disabled independently on each axis.
    defaults = {
        "x_major_grid": False,
        "x_minor_grid": False,
        "y_major_grid": False,
        "y_minor_grid": False,
    }
    for axis_name in ("x", "y"):
        for which in ("major", "minor"):
            dest = f"{axis_name}_{which}_grid"
            enabled = getattr(args, dest)
            if enabled is None:
                enabled = defaults[dest]
            if not enabled:
                continue
            is_major = which == "major"
            ax.grid(
                True,
                axis=axis_name,
                which=which,
                color=args.major_grid_color if is_major else args.minor_grid_color,
                alpha=args.major_grid_alpha if is_major else args.minor_grid_alpha,
                linewidth=args.major_grid_width if is_major else args.minor_grid_width,
                linestyle=args.major_grid_style if is_major else args.minor_grid_style,
            )
    ax.set_axisbelow(True)
    ax.tick_params(which="minor", labelbottom=False, labelleft=False)


def _finish(
    ax,
    fig,
    args,
    output: Path,
    *,
    bp_x: bool = False,
    default_x_major: float | None = None,
    default_x_minor: float | None = None,
    legend: bool = False,
    artist_kw: Mapping[str, Mapping[str, Any]] | None = None,
    default_size: tuple[float, float] = (10.0, 5.5),
    preserve_canvas: bool = False,
) -> Path:
    if args.title is not None:
        ax.set_title(args.title)
    elif args.no_title:
        ax.set_title("")
    if args.x_label is not None:
        ax.set_xlabel(args.x_label)
    if args.y_label is not None:
        ax.set_ylabel(args.y_label)
    if args.x_min is not None or args.x_max is not None:
        left, right = ax.get_xlim()
        ax.set_xlim(left if args.x_min is None else args.x_min, right if args.x_max is None else args.x_max)
    if args.y_min is not None or args.y_max is not None:
        bottom, top = ax.get_ylim()
        ax.set_ylim(bottom if args.y_min is None else args.y_min, top if args.y_max is None else args.y_max)
    _configure_ticks_and_grids(ax, args, bp_x=bp_x, default_x_major=default_x_major, default_x_minor=default_x_minor)
    if args.x_tick_rotation is not None:
        for label in ax.get_xticklabels():
            label.set_rotation(args.x_tick_rotation)
    if args.y_tick_rotation is not None:
        for label in ax.get_yticklabels():
            label.set_rotation(args.y_tick_rotation)
    if legend and not args.no_legend:
        kwargs = {"frameon": False}
        kwargs.update((artist_kw or {}).get("legend", {}))
        ax.legend(**kwargs)
    if args.axes_facecolor is not None:
        ax.set_facecolor(args.axes_facecolor)
    width = default_size[0] if args.width is None else args.width
    height = default_size[1] if args.height is None else args.height
    fig.set_size_inches(width, height, forward=True)
    if preserve_canvas and ax.get_title():
        import textwrap

        ax.set_title(textwrap.fill(ax.get_title(), width=70))
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output,
        dpi=args.dpi,
        bbox_inches=None if preserve_canvas else "tight",
        transparent=args.transparent,
    )
    return output


def _running_average(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values.astype(float, copy=True)
    if window % 2 == 0:
        raise ValueError("--smooth-window must be odd")
    finite = np.isfinite(values)
    clean = np.where(finite, values, 0.0)
    kernel = np.ones(window, dtype=float)
    numerator = np.convolve(clean, kernel, mode="same")
    denominator = np.convolve(finite.astype(float), kernel, mode="same")
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)


def _detect_peaks(x: np.ndarray, y: np.ndarray, min_distance_bp: float) -> np.ndarray:
    from scipy.signal import find_peaks
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 3:
        return np.asarray([], dtype=int)
    indices = np.flatnonzero(finite)
    xf = x[finite]
    yf = y[finite]
    steps = np.diff(np.unique(xf))
    step = float(np.median(steps[steps > 0])) if np.any(steps > 0) else 1.0
    distance_samples = max(1, int(round(float(min_distance_bp) / step)))
    local, _ = find_peaks(yf, distance=distance_samples)
    return indices[local]


def _label_peaks(
    ax, x: np.ndarray, y: np.ndarray, indices: Iterable[int], args, *, default: bool
) -> None:
    if args.label_peaks == "none" or (args.label_peaks == "auto" and not default):
        return
    for idx in indices:
        value = x[idx] if args.peak_label_value == "x" else y[idx]
        if args.peak_label_value == "both":
            text = f"{x[idx]:g}, {y[idx]:g}"
        else:
            text = f"{value:g}"
        ax.annotate(
            text,
            (x[idx], y[idx]),
            xytext=(0, args.peak_label_offset),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=max(6.0, args.font_size * 0.78),
            annotation_clip=False,
        )


def _add_nrl_inset(ax, x: np.ndarray, y: np.ndarray, indices: np.ndarray, args, artist_kw) -> None:
    if args.nrl_inset == "off" or indices.size < 3:
        return
    if args.nrl_inset == "auto" and (np.nanmax(x) - np.nanmin(x)) < 500:
        return
    peak_x = x[indices]
    mask = np.isfinite(peak_x) & (peak_x >= args.nrl_min_distance)
    peak_x = peak_x[mask]
    if peak_x.size < 3:
        return
    order = np.arange(1, peak_x.size + 1, dtype=float)
    slope, intercept = np.polyfit(order, peak_x, 1)
    fitted = intercept + slope * order
    ss_res = float(np.sum((peak_x - fitted) ** 2))
    ss_tot = float(np.sum((peak_x - np.mean(peak_x)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
    inset = ax.inset_axes(args.inset_bounds)
    point_kw = {"s": 18, "facecolors": "none", "edgecolors": "black", "linewidths": 0.8}
    point_kw.update(artist_kw.get("points", {}))
    inset.scatter(order, peak_x, **point_kw)
    line_kw = {"color": "black", "linewidth": 1.0, "linestyle": ":"}
    line_kw.update(artist_kw.get("line", {}))
    inset.plot(order, fitted, **line_kw)
    inset.set_xlabel("Peak #", fontsize=7)
    inset.set_ylabel("Distance (bp)", fontsize=7)
    inset.tick_params(labelsize=7)
    inset.text(0.06, 0.94, f"Slope\n= {slope:.2f}\n\n$R^2$ = {r2:.4f}", transform=inset.transAxes, va="top", fontsize=7)


def _plot_dac(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    x = _numeric(rows, hm["distance"])
    y = _numeric(rows, hm["dac_value"])
    mask = np.isfinite(x) & np.isfinite(y) & (x != 1)
    x, y = x[mask], y[mask]
    fig, ax = plt.subplots()

    if args.detect_peaks and x.size:
        # When peak detection is requested, reproduce the same profile layers
        # and resolution-driven caller used by `nucleosuite nrl`.
        from nucleosuite.nrl import (
            call_resolution_peaks, moving_average_by_distance,
            resolution_smoothing_windows,
        )
        detection_window, local_window = resolution_smoothing_windows(args.peak_resolution)
        local_values = moving_average_by_distance(x, y, local_window)
        detection_values = moving_average_by_distance(x, y, detection_window)
        peaks = call_resolution_peaks(x, y, local_values, detection_values, args.peak_resolution)
        indices = np.asarray([peak.index for peak in peaks], dtype=int)

        raw_kw = {"color": "0.72", "linewidth": 0.9, "label": "Unsmoothed"}
        raw_kw.update(artist_kw.get("raw", {}))
        ax.plot(x, y, **raw_kw)

        smooth_kw = {
            "color": "black", "linewidth": 1.5,
            "label": (f"Local maxima ({local_window} bp)" if local_window > 1 else "Local maxima signal"),
        }
        smooth_kw.update(artist_kw.get("smooth", artist_kw.get("line", {})))
        ax.plot(x, local_values, **smooth_kw)

        if detection_window != local_window:
            detect_kw = {
                "color": "0.4", "linewidth": 1.2, "linestyle": "--",
                "label": (f"Peak detection ({detection_window} bp)" if detection_window > 1 else "Peak detection signal"),
            }
            detect_kw.update(artist_kw.get("line", {}))
            ax.plot(x, detection_values, **detect_kw)

        if indices.size:
            points_kw = {
                "s": 28, "facecolors": "white", "edgecolors": "black",
                "linewidths": 1.0, "zorder": 5, "label": "Called peaks",
            }
            points_kw.update(artist_kw.get("points", {}))
            ax.scatter(x[indices], local_values[indices], **points_kw)
            _label_peaks(ax, x, local_values, indices, args, default=True)
            _add_nrl_inset(ax, x, local_values, indices, args, artist_kw)
        legend = True
    else:
        # Default DAC replot: the measured DAC profile only. No smoothing,
        # peak detection, peak markers, or peak labels are added.
        line_kw = {
            "linewidth": 1.2,
            "marker": "o",
            "markersize": 2.0,
            "markeredgewidth": 0,
        }
        line_kw.update(artist_kw.get("raw", artist_kw.get("line", {})))
        ax.plot(x, y, **line_kw)
        legend = False

    ax.set_xlabel("Distance (bp)")
    ax.set_ylabel("DAC value")
    ax.set_title(path.stem.replace("_", " "))
    artist_kw = {key: dict(value) for key, value in artist_kw.items()}
    if legend:
        artist_kw.setdefault("legend", {}).setdefault("loc", "lower center")
        artist_kw.setdefault("legend", {}).setdefault("ncol", 3)
    if args.x_major_grid is None and args.x_minor_grid is None:
        ax.grid(axis="x", alpha=0.5)
    return _finish(
        ax,
        fig,
        args,
        output,
        bp_x=True,
        legend=legend,
        artist_kw=artist_kw,
        default_size=(10.0, 5.0),
    ), fig

def _plot_dcc(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    x_key = hm.get("lag", hm.get("distance"))
    x = _numeric(rows, x_key)
    y = _numeric(rows, hm["dcc_value"])
    signed_lags = bool(np.any(x < 0))
    exclude_one = np.abs(x) != 1 if signed_lags else x != 1
    mask = np.isfinite(x) & np.isfinite(y) & exclude_one
    x, y = x[mask], y[mask]
    fig, ax = plt.subplots()
    line_kw = {"linewidth": 1.2, "marker": "o", "markersize": 2.0, "markeredgewidth": 0}
    line_kw.update(artist_kw.get("line", {}))
    ax.plot(x, y, **line_kw)
    ax.set_xlabel("Distance (bp)" if _clean_header(x_key) == "distance" else "Lag (bp)")
    ax.set_ylabel("DCC value")
    ax.set_title(path.stem.replace("_", " "))
    if args.x_major_grid is None and args.x_minor_grid is None:
        ax.grid(axis="x", alpha=0.5)
    return _finish(
        ax,
        fig,
        args,
        output,
        bp_x=True,
        artist_kw=artist_kw,
        default_size=(10.0, 5.0),
    ), fig


def _plot_nrl_profile(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    x = _numeric(rows, hm["distance_or_lag"])
    raw = _numeric(rows, hm["unsmoothed_value"])
    local = _numeric(rows, hm["local_max_smoothed_value"])
    detection = _numeric(rows, hm.get("detection_smoothed_value", hm["local_max_smoothed_value"]))
    is_peak = _numeric(rows, hm.get("is_peak", hm["unsmoothed_value"])) if "is_peak" in hm else np.zeros(x.size)
    fig, ax = plt.subplots()
    raw_kw = {"color": "0.72", "linewidth": 0.9, "label": "Unsmoothed"}; raw_kw.update(artist_kw.get("raw", {}))
    smooth_kw = {"color": "black", "linewidth": 1.5, "label": "Local-max smoothed"}; smooth_kw.update(artist_kw.get("smooth", {}))
    ax.plot(x, raw, **raw_kw); ax.plot(x, local, **smooth_kw)
    if not np.allclose(local, detection, equal_nan=True):
        detect_kw = {"linewidth": 1.0, "linestyle": "--", "label": "Detection smoothed"}; detect_kw.update(artist_kw.get("line", {}))
        ax.plot(x, detection, **detect_kw)
    indices = np.flatnonzero(is_peak > 0)
    if indices.size:
        pkw = {"s": 28, "facecolors": "white", "edgecolors": "black", "linewidths": 1.0, "zorder": 5, "label": "Called peaks"}; pkw.update(artist_kw.get("points", {}))
        ax.scatter(x[indices], local[indices], **pkw); _label_peaks(ax, x, local, indices, args, default=True)
        _add_nrl_inset(ax, x, local, indices, args, artist_kw)
    ax.set_xlabel("Distance or lag (bp)"); ax.set_ylabel("Signal"); ax.set_title(path.stem.replace("_", " "))
    return _finish(ax, fig, args, output, bp_x=True, legend=True, artist_kw=artist_kw), fig


def _plot_nrl_regression(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    if "order" in hm and "peak_distance_bp" in hm:
        order = _numeric(rows, hm["order"]); y = _numeric(rows, hm["peak_distance_bp"])
        fitted = _numeric(rows, hm["fitted_distance_bp"]) if "fitted_distance_bp" in hm else np.full_like(y, np.nan)
    else:
        order = _numeric(rows, hm["peak_number"]); y = _numeric(rows, hm["distance_or_lag"]); fitted = np.full_like(y, np.nan)
    mask = np.isfinite(order) & np.isfinite(y); order, y, fitted = order[mask], y[mask], fitted[mask]
    if not np.isfinite(fitted).all() and order.size >= 2:
        fit_slope, fit_intercept = np.polyfit(order, y, 1)
        fitted = fit_intercept + fit_slope * order
    fig, ax = plt.subplots()
    from_distance_command = "scope" in hm or "peak_distance_bp" in hm
    point_label = "Highest-count distance" if from_distance_command else "Called peaks"
    pkw = {
        "s": 34,
        "facecolors": "none",
        "edgecolors": "black",
        "linewidths": 1.0,
        "label": point_label,
    }
    pkw.update(artist_kw.get("points", {}))
    ax.scatter(order, y, **pkw)
    lkw = {
        "color": "black",
        "linewidth": 1.5,
        "linestyle": ":",
        "label": "Linear regression",
    }
    lkw.update(artist_kw.get("line", {}))
    ax.plot(order, fitted, **lkw)
    if order.size >= 2:
        slope, intercept = np.polyfit(order, y, 1); predicted = intercept + slope * order; ssr = np.sum((y-predicted)**2); sst=np.sum((y-np.mean(y))**2); r2=1-ssr/sst if sst else math.nan
        prefix = "NRL (slope)" if from_distance_command else "Slope"
        ax.text(
            0.04,
            0.96,
            f"{prefix} = {slope:.3f} bp\n$R^2$ = {r2:.4f}",
            transform=ax.transAxes,
            va="top",
        )
    ax.set_xlabel("Neighbour order" if from_distance_command else "Peak number")
    ax.set_ylabel("Peak distance (bp)")
    if from_distance_command and "scope" in hm:
        scopes = [str(value) for value in _column(rows, hm["scope"]) if str(value)]
        chromosomes = (
            [str(value) for value in _column(rows, hm["chromosome"]) if str(value)]
            if "chromosome" in hm else []
        )
        if scopes and scopes[0] == "combined_chromosomes":
            scope_label = "Combined chromosomes"
        elif chromosomes:
            scope_label = chromosomes[0]
        else:
            scope_label = scopes[0] if scopes else "Distance"
        ax.set_title(f"{scope_label} nucleosome repeat length regression")
    else:
        ax.set_title(path.stem.replace("_", " "))
    from nucleosuite.plotting import apply_integer_x_axis, apply_integer_y_axis
    apply_integer_x_axis(ax, order)
    apply_integer_y_axis(ax)
    ax.grid(False)
    return _finish(
        ax,
        fig,
        args,
        output,
        legend=True,
        artist_kw=artist_kw,
        default_size=(6.5, 6.5),
        preserve_canvas=True,
    ), fig


def _plot_fragment_size_nrl_profile(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt

    hm = _header_map(headers)
    x = _numeric(rows, hm["fragment_length"])
    raw = _numeric(rows, hm["unsmoothed_density"])
    local = _numeric(rows, hm["local_max_smoothed_density"])
    detection = _numeric(
        rows,
        hm.get("detection_smoothed_density", hm["local_max_smoothed_density"]),
    )
    is_peak = (
        _numeric(rows, hm["is_peak"])
        if "is_peak" in hm
        else np.zeros(x.size, dtype=float)
    )
    mask = np.isfinite(x) & np.isfinite(raw) & np.isfinite(local)
    x, raw, local, detection, is_peak = (
        x[mask],
        raw[mask],
        local[mask],
        detection[mask],
        is_peak[mask],
    )
    fig, ax = plt.subplots()
    raw_kw = {"color": "0.72", "linewidth": 0.9, "label": "Unsmoothed"}
    raw_kw.update(artist_kw.get("raw", {}))
    smooth_kw = {
        "color": "black",
        "linewidth": 1.5,
        "label": "Local maxima signal",
    }
    if "local_max_smoothing_window" in hm:
        local_window = int(float(rows[0][hm["local_max_smoothing_window"]]))
        if local_window > 1:
            smooth_kw["label"] = f"Local maxima ({local_window} bp)"
    smooth_kw.update(artist_kw.get("smooth", {}))
    ax.plot(x, raw, **raw_kw)
    ax.plot(x, local, **smooth_kw)
    if not np.allclose(local, detection, equal_nan=True):
        detect_kw = {
            "color": "0.4",
            "linewidth": 1.2,
            "linestyle": "--",
            "label": "Peak detection signal",
        }
        if "detection_smoothing_window" in hm:
            detection_window = int(float(rows[0][hm["detection_smoothing_window"]]))
            if detection_window > 1:
                detect_kw["label"] = f"Peak detection ({detection_window} bp)"
        detect_kw.update(artist_kw.get("line", {}))
        ax.plot(x, detection, **detect_kw)
    indices = np.flatnonzero(is_peak > 0)
    if indices.size:
        point_kw = {
            "s": 28,
            "facecolors": "white",
            "edgecolors": "black",
            "linewidths": 1.0,
            "zorder": 5,
            "label": "Called peaks",
        }
        point_kw.update(artist_kw.get("points", {}))
        ax.scatter(x[indices], local[indices], **point_kw)
        _label_peaks(ax, x, local, indices, args, default=False)
    labels = [str(value) for value in _column(rows, hm["label"])] if "label" in hm else []
    label = next((value for value in labels if value), "all")
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Density")
    ax.set_title("Fragment-size NRL" if label == "all" else f"{label}: fragment-size NRL")
    return _finish(
        ax,
        fig,
        args,
        output,
        bp_x=True,
        legend=True,
        artist_kw=artist_kw,
        default_size=(10.0, 5.5),
    ), fig


def _plot_fragment_size_nrl_regression(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt

    hm = _header_map(headers)
    order = _numeric(rows, hm["peak_number"])
    y = _numeric(rows, hm["fragment_length"])
    fitted = _numeric(rows, hm["fitted_fragment_length"])
    mask = np.isfinite(order) & np.isfinite(y)
    order, y, fitted = order[mask], y[mask], fitted[mask]
    if not np.isfinite(fitted).all() and order.size >= 2:
        slope, intercept = np.polyfit(order, y, 1)
        fitted = intercept + slope * order

    fig, ax = plt.subplots()
    point_kw = {
        "s": 34,
        "facecolors": "none",
        "edgecolors": "black",
        "linewidths": 1.0,
        "label": "Called peaks",
    }
    point_kw.update(artist_kw.get("points", {}))
    ax.scatter(order, y, **point_kw)
    line_kw = {
        "color": "black",
        "linewidth": 1.5,
        "linestyle": ":",
        "label": "Linear regression",
    }
    line_kw.update(artist_kw.get("line", {}))
    ax.plot(order, fitted, **line_kw)
    if order.size >= 2:
        slope, intercept = np.polyfit(order, y, 1)
        predicted = intercept + slope * order
        residual_sum = float(np.sum((y - predicted) ** 2))
        total_sum = float(np.sum((y - np.mean(y)) ** 2))
        r_squared = 1.0 - residual_sum / total_sum if total_sum else math.nan
        ax.text(
            0.04,
            0.96,
            f"Fragment-size NRL = {slope:.3f} bp/peak\n$R^2$ = {r_squared:.4f}",
            transform=ax.transAxes,
            va="top",
        )
    labels = [str(value) for value in _column(rows, hm["label"])] if "label" in hm else []
    label = next((value for value in labels if value), "all")
    ax.set_xlabel("Peak number")
    ax.set_ylabel("Peak fragment length (bp)")
    title = "Fragment-size NRL regression"
    ax.set_title(title if label == "all" else f"{label}: {title}")
    from nucleosuite.plotting import apply_integer_x_axis, apply_integer_y_axis
    apply_integer_x_axis(ax, order)
    apply_integer_y_axis(ax)
    ax.grid(False)
    return _finish(
        ax,
        fig,
        args,
        output,
        legend=True,
        artist_kw=artist_kw,
        default_size=(6.5, 6.5),
        preserve_canvas=True,
    ), fig


def _plot_aggregate_nrl_profile(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt

    hm = _header_map(headers)
    x = _numeric(rows, hm["relative_position"])
    raw = _numeric(rows, hm["unsmoothed_value"])
    local = _numeric(rows, hm["local_max_smoothed_value"])
    detection = _numeric(rows, hm["detection_smoothed_value"])
    is_peak = _numeric(rows, hm["is_peak"]) if "is_peak" in hm else np.zeros(x.size)
    mask = np.isfinite(x) & np.isfinite(raw) & np.isfinite(local) & np.isfinite(detection)
    x, raw, local, detection, is_peak = (
        x[mask], raw[mask], local[mask], detection[mask], is_peak[mask]
    )
    fig, ax = plt.subplots()
    if (
        "regression_exclusion_start_bp" in hm
        and "regression_exclusion_end_bp" in hm
    ):
        exclusion_starts = _numeric(rows, hm["regression_exclusion_start_bp"])
        exclusion_ends = _numeric(rows, hm["regression_exclusion_end_bp"])
        finite_starts = exclusion_starts[np.isfinite(exclusion_starts)]
        finite_ends = exclusion_ends[np.isfinite(exclusion_ends)]
        if finite_starts.size and finite_ends.size:
            ax.axvspan(
                float(finite_starts[0]),
                float(finite_ends[0]),
                color="0.9",
                alpha=0.55,
                linewidth=0,
                label="Regression exclusion zone",
                zorder=0,
            )
    raw_kw = {"color": "0.72", "linewidth": 0.9, "label": "Unsmoothed"}
    raw_kw.update(artist_kw.get("raw", {}))
    local_window = (
        int(float(rows[0][hm["local_max_smoothing_window"]]))
        if "local_max_smoothing_window" in hm
        else None
    )
    smooth_kw = {
        "color": "black",
        "linewidth": 1.5,
        "label": (
            f"Local maxima ({local_window} bp)"
            if local_window and local_window > 1
            else "Local maxima signal"
        ),
    }
    smooth_kw.update(artist_kw.get("smooth", {}))
    detection_window = (
        int(float(rows[0][hm["detection_smoothing_window"]]))
        if "detection_smoothing_window" in hm
        else None
    )
    detection_kw = {
        "color": "0.4",
        "linewidth": 1.2,
        "linestyle": "--",
        "label": (
            f"Peak detection ({detection_window} bp)"
            if detection_window and detection_window > 1
            else "Peak detection signal"
        ),
    }
    detection_kw.update(artist_kw.get("line", {}))
    ax.plot(x, raw, **raw_kw)
    ax.plot(x, local, **smooth_kw)
    if not np.allclose(local, detection, equal_nan=True):
        ax.plot(x, detection, **detection_kw)
    indices = np.flatnonzero(is_peak > 0)
    if indices.size:
        point_kw = {
            "s": 28,
            "facecolors": "white",
            "edgecolors": "black",
            "linewidths": 1.0,
            "zorder": 5,
            "label": "Called peaks",
        }
        point_kw.update(artist_kw.get("points", {}))
        ax.scatter(x[indices], local[indices], **point_kw)
        _label_peaks(ax, x, local, indices, args, default=False)
    ax.axvline(0, color="0.55", linewidth=0.8, zorder=0)
    title = rows[0].get(hm.get("plot_title", ""), path.stem.replace("_", " "))
    x_label = rows[0].get(hm.get("x_label", ""), "Relative position (bp)")
    y_label = rows[0].get(hm.get("y_label", ""), "Mean signal")
    ax.set_title(str(title))
    ax.set_xlabel(str(x_label))
    ax.set_ylabel(str(y_label))
    return _finish(
        ax,
        fig,
        args,
        output,
        bp_x=True,
        legend=True,
        artist_kw=artist_kw,
        default_size=(10.0, 5.5),
    ), fig


def _plot_aggregate_nrl_regression(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt

    hm = _header_map(headers)
    order = _numeric(rows, hm["peak_number"])
    distance = _numeric(rows, hm["distance_from_zero_bp"])
    fitted = _numeric(rows, hm["fitted_distance_from_zero_bp"])
    mask = np.isfinite(order) & np.isfinite(distance)
    order, distance, fitted = order[mask], distance[mask], fitted[mask]
    direction_values = [
        str(value) for value in _column(rows, hm["direction"]) if str(value)
    ]
    direction = direction_values[0] if direction_values else "directional"
    title = rows[0].get(
        hm.get("plot_title", ""),
        f"{direction}-direction repeat-length regression",
    )

    fig, ax = plt.subplots()
    if order.size:
        point_kw = {
            "s": 34,
            "facecolors": "none",
            "edgecolors": "black",
            "linewidths": 1.0,
            "label": "Called peaks",
        }
        point_kw.update(artist_kw.get("points", {}))
        ax.scatter(order, distance, **point_kw)
        if not np.isfinite(fitted).all() and order.size >= 2:
            slope, intercept = np.polyfit(order, distance, 1)
            fitted = intercept + slope * order
        if np.isfinite(fitted).all():
            line_kw = {
                "color": "black",
                "linewidth": 1.5,
                "linestyle": ":",
                "label": "Linear regression",
            }
            line_kw.update(artist_kw.get("line", {}))
            ax.plot(order, fitted, **line_kw)
        if order.size >= 2:
            slope, intercept = np.polyfit(order, distance, 1)
            predicted = intercept + slope * order
            residual_sum = float(np.sum((distance - predicted) ** 2))
            total_sum = float(np.sum((distance - np.mean(distance)) ** 2))
            r_squared = 1.0 - residual_sum / total_sum if total_sum else math.nan
            ax.text(
                0.04,
                0.96,
                f"{direction.title()} repeat length = {slope:.3f} bp/peak\n"
                f"$R^2$ = {r_squared:.4f}",
                transform=ax.transAxes,
                va="top",
            )
    else:
        ax.text(
            0.5,
            0.5,
            "No peaks called",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )
    from nucleosuite.plotting import apply_integer_x_axis, apply_integer_y_axis

    apply_integer_x_axis(ax, order)
    apply_integer_y_axis(ax)
    ax.set_xlabel("Peak number")
    ax.set_ylabel("Distance from position 0 (bp)")
    ax.set_title(str(title))
    ax.grid(False)
    return _finish(
        ax,
        fig,
        args,
        output,
        legend=bool(order.size),
        artist_kw=artist_kw,
        default_size=(6.5, 6.5),
        preserve_canvas=True,
    ), fig


def _grouped_line(
    path,
    headers,
    rows,
    args,
    output,
    artist_kw,
    *,
    x_candidates,
    y_candidates,
    group_candidates=(),
    xlabel=None,
    ylabel=None,
    bp_x=False,
    line_defaults: Mapping[str, Any] | None = None,
    default_size: tuple[float, float] = (10.0, 5.5),
    default_title: str | None = None,
    legend_defaults: Mapping[str, Any] | None = None,
):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    x_key = args.x_column or next((hm[k] for k in x_candidates if k in hm), None)
    y_key = args.y_column or next((hm[k] for k in y_candidates if k in hm), None)
    if x_key is None or y_key is None:
        raise ValueError(f"Could not determine x/y columns for {path}; use --x-column and --y-column")
    group_key = args.group_column or next((hm[k] for k in group_candidates if k in hm), None)
    fig, ax = plt.subplots()
    groups: dict[str, list[int]] = defaultdict(list)
    if group_key:
        for idx, row in enumerate(rows): groups[str(row.get(group_key, ""))].append(idx)
    else:
        groups[""] = list(range(len(rows)))
    for group, indices in groups.items():
        subrows = [rows[i] for i in indices]
        x = _numeric(subrows, x_key); y = _numeric(subrows, y_key)
        mask=np.isfinite(x)&np.isfinite(y); x,y=x[mask],y[mask]; order=np.argsort(x); x,y=x[order],y[order]
        lkw={"linewidth":1.4, "marker":"o", "markersize":2.2, "markeredgewidth":0}
        lkw.update(line_defaults or {})
        lkw.update(artist_kw.get("line", {}))
        if group: lkw["label"] = group
        ax.plot(x,y,**lkw)
    ax.set_xlabel(xlabel or x_key)
    ax.set_ylabel(ylabel or y_key)
    if default_title is not None:
        ax.set_title(default_title)
    configured_artist_kw = {key: dict(value) for key, value in artist_kw.items()}
    if legend_defaults:
        configured_artist_kw.setdefault("legend", {})
        configured_artist_kw["legend"] = {
            **legend_defaults,
            **configured_artist_kw["legend"],
        }
    return _finish(
        ax,
        fig,
        args,
        output,
        bp_x=bp_x,
        legend=len(groups)>1,
        artist_kw=configured_artist_kw,
        default_size=default_size,
    ), fig


def _plot_distances(path, headers, rows, args, output, artist_kw):
    """Recreate the standard solid-line neighbour-order distribution plot."""
    import matplotlib.pyplot as plt

    hm = _header_map(headers)
    x_key = args.x_column or hm.get("distance_bp") or hm.get("distance")
    if x_key is None:
        raise ValueError("Distance replots require a distance_bp or distance column")

    selected_rows = list(rows)
    if "scope" in hm:
        combined = [row for row in selected_rows if row.get(hm["scope"]) == "combined_chromosomes"]
        if combined:
            selected_rows = combined
    if "state" in hm:
        pooled = [row for row in selected_rows if row.get(hm["state"]) == "All"]
        if pooled:
            selected_rows = pooled

    order_key = args.group_column or hm.get("order")
    groups: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    if order_key is None:
        groups[""] = selected_rows
    else:
        for row in selected_rows:
            groups[str(row.get(order_key, ""))].append(row)

    def group_sort(value: str):
        try:
            return (0, float(value))
        except ValueError:
            return (1, _natural_key(value))

    fig, ax = plt.subplots()
    for group in sorted(groups, key=group_sort):
        group_rows = groups[group]
        x = _numeric(group_rows, x_key)
        raw_key = args.y_column or hm.get("count") or hm.get("raw_count")
        smooth_key = None if args.y_column else hm.get("smoothed_count")
        if raw_key is None:
            raw_key = hm.get("raw_percent") or hm.get("percent") or hm.get("smoothed_percent")
        if raw_key is None:
            raise ValueError("Distance replots require a count or percentage column")
        raw = _numeric(group_rows, raw_key)
        smooth = _numeric(group_rows, smooth_key) if smooth_key is not None else raw
        mask = np.isfinite(x) & np.isfinite(raw) & np.isfinite(smooth)
        x, raw, smooth = x[mask], raw[mask], smooth[mask]
        order = np.argsort(x)
        x, raw, smooth = x[order], raw[order], smooth[order]
        if not x.size:
            continue

        if smooth_key is not None:
            raw_kw = {"color": "0.78", "linewidth": 0.8, "alpha": 0.65, "zorder": 1}
            raw_kw.update(artist_kw.get("raw", {}))
            ax.plot(x, raw, **raw_kw)

        line_kw = {"linewidth": 1.5 if smooth_key is not None else 1.35, "linestyle": "-", "zorder": 2}
        line_kw.update(artist_kw.get("line", {}))
        if group:
            line_kw["label"] = f"+{group}" if re.fullmatch(r"\d+(?:\.0+)?", group) else group
        line, = ax.plot(x, smooth, **line_kw)

        peak_index = int(np.argmax(smooth))
        point_kw = {
            "s": 22,
            "facecolors": "white",
            "edgecolors": line.get_color(),
            "linewidths": 0.8,
            "zorder": 4,
        }
        point_kw.update(artist_kw.get("points", {}))
        ax.scatter([x[peak_index]], [smooth[peak_index]], **point_kw)
        _label_peaks(ax, x, smooth, [peak_index], args, default=False)

    ax.set_xlabel("Distance (bp)")
    if args.y_column:
        ax.set_ylabel(args.y_column.replace("_", " ").title())
    else:
        ax.set_ylabel("Count" if (hm.get("count") or hm.get("raw_count")) else "Percentage")
    return _finish(
        ax,
        fig,
        args,
        output,
        bp_x=True,
        legend=len(groups) > 1,
        artist_kw={
            **artist_kw,
            "legend": {
                "title": "Neighbour order",
                **artist_kw.get("legend", {}),
            },
        },
        default_size=(10.0, 5.5),
    ), fig


def _plot_fragment_lengths(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import apply_integer_y_axis

    hm = _header_map(headers)
    x_key = args.x_column or hm.get("fragment_length")
    y_key = args.y_column or hm.get("count")
    if x_key is None or y_key is None:
        raise ValueError("Fragment-length replots require fragment_length and count columns")
    group_key = args.group_column or hm.get("label")
    groups: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[str(row.get(group_key, "")) if group_key else ""].append(row)

    filename = path.name.lower()
    from_fragments = args.from_command == "fragments" or "fragment_length_counts" in filename
    normalisation = args.normalization
    if normalisation == "auto":
        normalisation = "count" if from_fragments else "density"

    fig, ax = plt.subplots()
    if from_fragments:
        # The `fragments` command plots only the observed support.
        plot_minimum = plot_maximum = None
    else:
        # `fragment-lengths --plot` displays from zero through the longest
        # counted fragment, capped at 1000 bp by default, and includes zero-count
        # positions before calculating within-window density.
        all_lengths = _numeric(rows, x_key)
        finite_lengths = all_lengths[np.isfinite(all_lengths)]
        plot_minimum = 0
        plot_maximum = int(min(float(np.max(finite_lengths)), 1000)) if finite_lengths.size else 1000

    plotted = 0
    for label in sorted(groups):
        x = _numeric(groups[label], x_key)
        y = _numeric(groups[label], y_key)
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if not x.size:
            continue
        if plot_minimum is not None and plot_maximum is not None:
            values = {int(round(xx)): float(yy) for xx, yy in zip(x, y)}
            x = np.arange(plot_minimum, plot_maximum + 1, dtype=float)
            y = np.asarray([values.get(int(xx), 0.0) for xx in x], dtype=float)
        else:
            order = np.argsort(x)
            x, y = x[order], y[order]
        if normalisation == "density":
            total = float(np.sum(y))
            if total > 0:
                y = y / total
        line_kw = {"linewidth": 1.5 if not from_fragments else 1.3, "marker": "o", "markersize": 2.0 if not from_fragments else 2.2, "markeredgewidth": 0}
        line_kw.update(artist_kw.get("line", {}))
        if label:
            line_kw["label"] = label
        ax.plot(x, y, **line_kw)
        plotted += 1

    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Density" if normalisation == "density" else ("Fragment count" if from_fragments else "Count"))
    if normalisation == "count":
        apply_integer_y_axis(ax)
    if from_fragments:
        ax.set_title("Fragment-length distribution")
        ax.grid(axis="y", alpha=0.25)
    legend_kw = {key: dict(value) for key, value in artist_kw.items()}
    legend_kw.setdefault("legend", {})
    legend_kw["legend"] = {
        "loc": "center left", "bbox_to_anchor": (1.02, 0.5),
        **legend_kw["legend"],
    }
    default_size = (10.0, 5.0) if from_fragments else tuple(plt.rcParams["figure.figsize"])
    return _finish(
        ax, fig, args, output, bp_x=True, legend=plotted > 1,
        artist_kw=legend_kw, default_size=default_size,
    ), fig


def _plot_positive_runs(path, headers, rows, args, output, artist_kw):
    saved, fig = _grouped_line(
        path, headers, rows, args, output, artist_kw,
        x_candidates=("run_length_bp",),
        y_candidates=("count", "percent", "fraction"),
        xlabel="Contiguous positive run length (bp)", ylabel=None, bp_x=True,
        line_defaults={"linewidth": 1.4, "markersize": 2.0},
        default_size=(10.0, 6.0),
    )
    fig.set_size_inches(
        10.0 if args.width is None else args.width,
        6.0 if args.height is None else args.height,
        forward=True,
    )
    ax = fig.axes[0]
    ax.grid(True, alpha=0.25)
    if ax.lines and ax.lines[0].get_xdata().size:
        maximum = float(np.nanmax(ax.lines[0].get_xdata()))
        ax.set_xlim(0, min(maximum, 550.0))
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight", transparent=args.transparent)
    return saved, fig


def _plot_peak_score_frequency(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    x_key = args.x_column or hm.get("bin_midpoint") or hm.get("score")
    y_key = args.y_column or hm.get("count") or hm.get("percent") or hm.get("fraction") or hm.get("density")
    if x_key is None or y_key is None:
        raise ValueError("Peak-score frequency replots require score/bin_midpoint and frequency columns")
    dataset_key = args.group_column or hm.get("dataset")
    groups: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[str(row.get(dataset_key, "")) if dataset_key else ""].append(row)
    fig, ax = plt.subplots()
    for label in groups:
        x = _numeric(groups[label], x_key); y = _numeric(groups[label], y_key)
        mask = np.isfinite(x) & np.isfinite(y); x, y = x[mask], y[mask]
        order = np.argsort(x); x, y = x[order], y[order]
        kw = {"where": "mid", "linewidth": 1.5}
        kw.update(artist_kw.get("line", {}))
        if label: kw["label"] = label
        ax.step(x, y, **kw)
    ylabel = {
        "count": "Frequency (count)", "fraction": "Frequency (fraction)",
        "percent": "Frequency (%)", "density": "Probability density",
    }.get(_clean_header(y_key), y_key)
    ax.set_xlabel("Peak score"); ax.set_ylabel(ylabel); ax.set_title(path.stem.replace("_", " "))
    ax.grid(axis="y", alpha=0.25)
    return _finish(
        ax, fig, args, output, legend=len(groups) > 1, artist_kw=artist_kw,
        default_size=(10.0, 6.0),
    ), fig


def _plot_gene_expression(path, headers, rows, args, output, artist_kw):
    hm = _header_map(headers)
    if "period_bp" in hm or "rank" in hm:
        return _grouped_line(path, headers, rows, args, output, artist_kw, x_candidates=("period_bp","rank"), y_candidates=("correlation",), group_candidates=("profile","sample"), xlabel=None, ylabel="Correlation", bp_x="period_bp" in hm)
    import matplotlib.pyplot as plt
    correlation_key = hm.get("correlation")
    sample_key = hm.get("sample")
    if correlation_key is None:
        raise ValueError("Gene-expression replots require a correlation column")
    by_sample: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        try:
            value = float(row.get(correlation_key, "nan"))
        except ValueError:
            continue
        if math.isfinite(value):
            by_sample[str(row.get(sample_key, "sample")) if sample_key else "sample"].append(value)
    fig, ax = plt.subplots()
    for sample, values in sorted(by_sample.items()):
        values = sorted(values)
        line_kw = {"linewidth": 1.4, "marker": "o", "markersize": 2.2, "markeredgewidth": 0, "label": sample}
        line_kw.update(artist_kw.get("line", {}))
        ax.plot(np.arange(1, len(values)+1), values, **line_kw)
    ax.axhline(0, color="0.5", linewidth=0.8)
    ax.set_xlabel("Expression profile rank by correlation"); ax.set_ylabel("Correlation"); ax.set_title(path.stem.replace("_", " "))
    return _finish(ax, fig, args, output, legend=len(by_sample)>1, artist_kw=artist_kw), fig


def _plot_gene_expression_spacing(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers); sample_key = hm["sample"]; correlation_key = hm["correlation"]
    selected = [
        row for row in rows
        if str(row.get(hm["region"], "")) == "body" and str(row.get(hm["subset"], "")) == "all"
    ]
    by_sample: dict[str, list[float]] = defaultdict(list)
    for row in selected:
        try: value = float(row.get(correlation_key, "nan"))
        except ValueError: continue
        if math.isfinite(value): by_sample[str(row.get(sample_key, ""))].append(value)
    fig, ax = plt.subplots()
    for sample in sorted(by_sample):
        values = sorted(by_sample[sample])
        kw = {"marker": "o", "markersize": 2.2, "markeredgewidth": 0, "label": sample}
        kw.update(artist_kw.get("line", {})); ax.plot(np.arange(1, len(values) + 1), values, **kw)
    method = str(selected[0].get(hm.get("correlation_method", ""), "Pearson")) if selected else "Pearson"
    ax.axhline(0, linewidth=0.8); ax.set_xlabel("Expression profile rank by correlation")
    ax.set_ylabel(f"{method.title()} correlation: median spacing vs expression")
    ax.set_title("Gene-body nucleosome spacing and expression")
    return _finish(ax, fig, args, output, legend=bool(by_sample), artist_kw=artist_kw, default_size=(11.0, 6.0)), fig


def _plot_gene_expression_spacing_scatter(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    panels = []
    for row in rows:
        key = (str(row[hm["sample"]]), str(row[hm["profile"]]))
        if key not in panels: panels.append(key)
    panel_count = max(1, len(panels)); columns = min(3, panel_count); row_count = int(math.ceil(panel_count / columns))
    fig, axes = plt.subplots(row_count, columns, squeeze=False)
    for ax, (sample, profile) in zip(axes.flat, panels):
        selected = [row for row in rows if str(row[hm["sample"]]) == sample and str(row[hm["profile"]]) == profile]
        x = _numeric(selected, hm["body_median_spacing_bp"]); y = _numeric(selected, hm["transformed_expression"])
        mask = np.isfinite(x) & np.isfinite(y); x, y = x[mask], y[mask]
        kw = {"s": 8, "alpha": 0.35, "linewidths": 0}; kw.update(artist_kw.get("points", {})); ax.scatter(x, y, **kw)
        transform = str(selected[0].get(hm.get("expression_transform", ""), "expression")) if selected else "expression"
        correlation = float(selected[0].get(hm.get("correlation", ""), "nan")) if selected else math.nan
        matched = int(float(selected[0].get(hm.get("matched_genes", ""), len(x)))) if selected else len(x)
        ax.set_xlabel("Median adjacent peak spacing in gene body (bp)"); ax.set_ylabel(f"{transform} expression")
        ax.set_title(f"{sample} vs {profile}\nr={correlation:.3f}, n={matched:,}")
    for ax in axes.flat[len(panels):]: ax.set_axis_off()
    return _finish(
        axes.flat[0], fig, args, output, artist_kw=artist_kw,
        default_size=(5.2 * columns, 4.5 * row_count),
    ), fig


def _plot_gene_expression_fft_trajectory(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers); samples = sorted({str(row[hm["sample"]]) for row in rows}, key=_natural_key)
    fig, axes = plt.subplots(len(samples), 1, squeeze=False)
    for ax, sample in zip(axes.flat, samples):
        sample_rows = [row for row in rows if str(row[hm["sample"]]) == sample]
        by_profile: dict[str, list[Mapping[str, str]]] = defaultdict(list)
        for row in sample_rows: by_profile[str(row[hm["profile"]])].append(row)
        ranked = []
        for profile, profile_rows in by_profile.items():
            values = _numeric(profile_rows, hm["correlation"]); finite = values[np.isfinite(values)]
            ranked.append((float(np.mean(finite)) if finite.size else math.inf, profile))
        highlighted_values = {
            str(row[hm["profile"]])
            for row in sample_rows
            if "plot_highlight" in hm and str(row.get(hm["plot_highlight"], "")).lower() in {"1", "true", "yes"}
        }
        if not highlighted_values and ranked:
            highlighted_values = {min(ranked)[1]}
        for profile, profile_rows in by_profile.items():
            profile_rows = sorted(profile_rows, key=lambda row: float(row[hm["period_bp"]]))
            x = _numeric(profile_rows, hm["period_bp"]); y = _numeric(profile_rows, hm["correlation"])
            kw = {"linewidth": 0.6, "alpha": 0.25, "color": "0.5", "marker": "o", "markersize": 1.2, "markeredgewidth": 0}
            kw.update(artist_kw.get("raw", {})); ax.plot(x, y, **kw)
        for highlighted in sorted(highlighted_values, key=_natural_key):
            profile_rows = sorted(by_profile[highlighted], key=lambda row: float(row[hm["period_bp"]]))
            x = _numeric(profile_rows, hm["period_bp"]); y = _numeric(profile_rows, hm["correlation"])
            kw = {"linewidth": 2.0, "marker": "o", "markersize": 2.5, "label": highlighted}
            kw.update(artist_kw.get("line", {})); ax.plot(x, y, **kw)
        if highlighted_values: ax.legend(frameon=False)
        method = str(sample_rows[0].get(hm.get("correlation_method", ""), "Pearson")) if sample_rows else "Pearson"
        ax.axhline(0, linewidth=0.8)
        if "ranking_periods" in hm and sample_rows:
            ranking_periods = [float(value) for value in str(sample_rows[0].get(hm["ranking_periods"], "")).split(",") if value]
            if ranking_periods: ax.axvspan(min(ranking_periods), max(ranking_periods), alpha=0.08)
        ax.set_xlabel("Nucleosome period (bp)")
        ax.set_ylabel(f"{method.title()} correlation"); ax.set_title(f"{sample}: FFT intensity versus expression")
        from nucleosuite.plotting import apply_base_pair_x_axis
        all_periods = _numeric(sample_rows, hm["period_bp"]); apply_base_pair_x_axis(ax, all_periods)
    first = axes.flat[0]
    return _finish(first, fig, args, output, artist_kw=artist_kw, default_size=(10.0, 5.0 * len(samples))), fig


def _plot_gene_expression_ranking(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers); samples = sorted({str(row[hm["sample"]]) for row in rows}, key=_natural_key)
    top_profiles = 30
    fig, axes = plt.subplots(len(samples), 1, squeeze=False)
    for ax, sample in zip(axes.flat, samples):
        selected = [row for row in rows if str(row[hm["sample"]]) == sample and str(row.get(hm["rank"], "")).strip()]
        selected.sort(key=lambda row: int(float(row[hm["rank"]])))
        if "plot_selected" in hm:
            selected = [row for row in selected if str(row.get(hm["plot_selected"], "")).lower() in {"1", "true", "yes"}]
        else:
            selected = selected[:top_profiles]
        labels = [str(row[hm["profile"]]) for row in selected][::-1]
        values = [float(row[hm["correlation"]]) for row in selected][::-1]
        kw = {}; kw.update(artist_kw.get("bar", {})); ax.barh(np.arange(len(selected)), values, **kw)
        ax.set_yticks(np.arange(len(selected)), labels=labels); ax.axvline(0, linewidth=0.8)
        periods = str(selected[0].get(hm.get("ranking_periods", ""), "")) if selected else ""
        method = str(selected[0].get(hm.get("correlation_method", ""), "Pearson")) if selected else "Pearson"
        ax.set_xlabel(f"{method.title()} correlation with mean FFT intensity at {periods.replace(',', ', ')} bp")
        ax.set_title(f"{sample}: expression profiles ranked from most negative correlation")
    first = axes.flat[0]
    return _finish(first, fig, args, output, artist_kw=artist_kw, default_size=(10.0, max(5.0, 0.28 * top_profiles) * len(samples))), fig


def _plot_tss_expression(path, headers, rows, args, output, artist_kw):
    saved, fig = _grouped_line(
        path, headers, rows, args, output, artist_kw,
        x_candidates=("relative_position",), y_candidates=("mean_signal",),
        group_candidates=("quintile_label", "quintile"),
        xlabel="Position relative to TSS (bp)", ylabel="Mean signal", bp_x=True,
        line_defaults={"linewidth": 1.4, "marker": "None"},
        default_size=(10.0, 5.5),
        legend_defaults={"title": "Expression quintile"},
    )
    ax = fig.axes[0]
    ax.axvline(0, linewidth=0.8, linestyle="--")
    positions = _numeric(rows, _header_map(headers)["relative_position"])
    finite = positions[np.isfinite(positions)]
    if finite.size:
        ax.set_xlim(float(np.min(finite)), float(np.max(finite)))
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight", transparent=args.transparent)
    return saved, fig


def _plot_aggregate_profile(path, headers, rows, args, output, artist_kw):
    saved, fig = _grouped_line(
        path, headers, rows, args, output, artist_kw,
        x_candidates=("relative_position",), y_candidates=("score",),
        xlabel="Relative position (bp)", ylabel="Mean signal", bp_x=True,
        line_defaults={
            "color": "black", "linewidth": 1.5, "marker": "o",
            "markersize": 1.6, "markeredgewidth": 0,
        },
        default_size=(15.0, 5.0),
    )
    ax = fig.axes[0]
    if ax.lines and ax.lines[0].get_xdata().size:
        values = np.asarray(ax.lines[0].get_xdata(), dtype=float)
        ax.set_xlim(float(np.nanmin(values)), float(np.nanmax(values)))
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight", transparent=args.transparent)
    return saved, fig


def _plot_compare_positions(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    if "bin_start_inclusive" in hm and "bin_end_exclusive" in hm:
        start_key, end_key = hm["bin_start_inclusive"], hm["bin_end_exclusive"]
    elif "lower_exclusive" in hm and "upper_inclusive" in hm:
        start_key, end_key = hm["lower_exclusive"], hm["upper_inclusive"]
    else:
        raise ValueError("compare-positions histogram replots require numeric bin boundaries")
    comparison_key = hm.get("comparison")
    groups = []
    if comparison_key:
        for row in rows:
            label = str(row.get(comparison_key, ""))
            if label not in groups:
                groups.append(label)
    else:
        groups = [""]
    fig, ax = plt.subplots()
    for label in groups:
        subset = rows if comparison_key is None else [row for row in rows if str(row.get(comparison_key, "")) == label]
        start = _numeric(subset, start_key); end = _numeric(subset, end_key); counts = _numeric(subset, hm["pair_count"])
        mid = (start + end) / 2.0
        mask = np.isfinite(mid) & np.isfinite(counts)
        mid, counts = mid[mask], counts[mask]
        order = np.argsort(mid); mid, counts = mid[order], counts[order]
        kw = {"linewidth": 1.5}
        if label:
            kw["label"] = label
        kw.update(artist_kw.get("line", {}))
        ax.plot(mid, counts, **kw)
    all_start = _numeric(rows, start_key); all_end = _numeric(rows, end_key)
    finite_start = all_start[np.isfinite(all_start)]; finite_end = all_end[np.isfinite(all_end)]
    minimum = float(np.nanmin(finite_start)) if finite_start.size else 0.0
    maximum = float(np.nanmax(finite_end)) if finite_end.size else 300.0
    ax.set_xlim(minimum, maximum)
    ax.set_xlabel("Signed summit distance, comparison − main (bp)"); ax.set_ylabel("Matched pairs"); ax.set_title("Matched-position distance distributions")
    from nucleosuite.plotting import apply_distance_x_axis, apply_integer_y_axis
    if all(getattr(args, name) is None for name in ("x_major_grid", "x_minor_grid")):
        apply_distance_x_axis(ax, major_interval=args.x_major_tick, minor_interval=args.x_minor_tick); bp_x = False
    else:
        bp_x = True
    apply_integer_y_axis(ax)
    return _finish(ax, fig, args, output, bp_x=bp_x, legend=bool(comparison_key), artist_kw=artist_kw, default_size=(9.0, 5.8)), fig


def _plot_compare_positions_score(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    metadata_row = rows[0] if rows else {}
    normalization = str(metadata_row.get(hm.get("plot_score_normalization", ""), "zscore") or "zscore")

    compact_source = "main_score_plot" in hm and "compare_score_plot" in hm
    if compact_source:
        x_key, y_key = hm["main_score_plot"], hm["compare_score_plot"]
        score_label = str(metadata_row.get(hm.get("score_axis_label", ""), "score") or "score")
    elif "main_score" in hm and "compare_score" in hm:
        if normalization == "raw":
            x_key, y_key, score_label = hm["main_score"], hm["compare_score"], "raw score"
        else:
            x_key = hm.get("main_score_normalized", hm["main_score"])
            y_key = hm.get("compare_score_normalized", hm["compare_score"])
            score_label = "score percentile rank" if normalization == "percentile" else "score z-score"
    else:
        default_columns = {
            "raw": ("a_score", "b_score", "raw score"),
            "zscore": ("a_score_z", "b_score_z", "score z-score"),
            "percentile": ("a_score_percentile", "b_score_percentile", "score percentile rank"),
        }
        x_name, y_name, score_label = default_columns.get(normalization, default_columns["zscore"])
        x_key = args.x_column or hm.get(x_name) or hm.get("a_score_z") or hm.get("a_score")
        y_key = args.y_column or hm.get(y_name) or hm.get("b_score_z") or hm.get("b_score")
    distance_key = hm.get("absolute_distance")
    if x_key is None or y_key is None or distance_key is None:
        raise ValueError("Score-agreement replots require main/A score, comparison/B score, and absolute_distance columns")

    x_all = _numeric(rows, x_key)
    y_all = _numeric(rows, y_key)
    distance_all = _numeric(rows, distance_key)
    mask = np.isfinite(x_all) & np.isfinite(y_all) & np.isfinite(distance_all)
    x_all, y_all, distance_all = x_all[mask], y_all[mask], distance_all[mask]
    if "plot_selected" in hm:
        selected = _numeric(rows, hm["plot_selected"])[mask] > 0
    else:
        selected = np.ones(x_all.size, dtype=bool)
    x, y, distance = x_all[selected], y_all[selected], distance_all[selected]

    fig, ax = plt.subplots()
    kw = {"c": distance, "s": 9, "alpha": 0.55, "linewidths": 0, "cmap": "viridis", "rasterized": True}
    if "plot_score_agreement_distance_max" in hm:
        try:
            distance_max = float(metadata_row.get(hm["plot_score_agreement_distance_max"], 0.0))
        except (TypeError, ValueError):
            distance_max = 0.0
        if distance_max > 0:
            from matplotlib.colors import Normalize
            kw["norm"] = Normalize(vmin=0.0, vmax=distance_max, clip=True)
    kw.update(artist_kw.get("points", {}))
    scatter = ax.scatter(x, y, **kw)
    colorbar = fig.colorbar(scatter, ax=ax)
    colorbar.set_label("Absolute summit distance (bp)")

    if normalization == "zscore":
        z_limit = 10.0
        if "plot_score_z_limit" in hm:
            try:
                z_limit = float(metadata_row.get(hm["plot_score_z_limit"], 10.0))
            except (TypeError, ValueError):
                pass
        if z_limit > 0:
            ax.set_xlim(-z_limit, z_limit)
            ax.set_ylim(-z_limit, z_limit)

    label_a = str(metadata_row.get(hm.get("main_label", ""), metadata_row.get(hm.get("plot_label_a", ""), "Main")) or "Main")
    label_b = str(metadata_row.get(hm.get("comparison", ""), metadata_row.get(hm.get("plot_label_b", ""), "Comparison")) or "Comparison")
    ax.set_xlabel(f"{label_a} {score_label}")
    ax.set_ylabel(f"{label_b} {score_label}")
    ax.set_title("Score agreement coloured by summit distance")

    annotation: list[str] = []
    method = str(metadata_row.get(hm.get("plot_correlation_method", ""), "spearman") or "spearman")
    if compact_source and "full_spearman" in hm:
        try:
            spearman = float(metadata_row.get(hm["full_spearman"], "nan"))
            pearson = float(metadata_row.get(hm.get("full_pearson", ""), "nan"))
            r_squared = float(metadata_row.get(hm.get("full_linear_r_squared", ""), "nan"))
            n_full = int(float(metadata_row.get(hm.get("full_matched_pair_count", ""), len(x_all))))
        except (TypeError, ValueError):
            spearman = pearson = r_squared = math.nan
            n_full = len(x_all)
    elif x_all.size >= 2:
        try:
            from scipy.stats import pearsonr, spearmanr
            pearson = float(pearsonr(x_all, y_all).statistic)
            spearman = float(spearmanr(x_all, y_all).statistic)
        except Exception:
            pearson = spearman = math.nan
        slope, intercept = np.polyfit(x_all, y_all, 1)
        fitted = intercept + slope * x_all
        ss_res = float(np.sum((y_all - fitted) ** 2))
        ss_tot = float(np.sum((y_all - np.mean(y_all)) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot else math.nan
        n_full = len(x_all)
    else:
        spearman = pearson = r_squared = math.nan
        n_full = len(x_all)

    if method in {"spearman", "both"}:
        annotation.append(f"Spearman ρ = {spearman:.3f}")
    if method in {"pearson", "both"}:
        annotation.append(f"Pearson r = {pearson:.3f}")
    annotation.append(f"R² = {r_squared:.3f}")
    annotation.append(f"n = {n_full:,}")
    ax.text(0.02, 0.98, "\n".join(annotation), transform=ax.transAxes, va="top", ha="left")
    return _finish(ax, fig, args, output, artist_kw=artist_kw, default_size=(7.5, 6.5)), fig

def _plot_compare_positions_score_distance(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    metadata_row = rows[0] if rows else {}
    compact_source = "matched_distance" in hm and "distance_type" in hm

    x_key = args.x_column or hm.get("main_score") or hm.get("a_score")
    if x_key is None:
        raise ValueError("Score-distance replots require main_score (or a_score)")

    if args.y_column:
        y_key = args.y_column
    elif compact_source:
        y_key = hm["matched_distance"]
    else:
        y_key = hm.get("absolute_distance")
    if y_key is None:
        raise ValueError("Score-distance replots require matched_distance, absolute_distance, or --y-column")

    x_all = _numeric(rows, x_key)
    y_all = _numeric(rows, y_key)
    mask = np.isfinite(x_all) & np.isfinite(y_all)
    x_all, y_all = x_all[mask], y_all[mask]
    if "plot_selected" in hm:
        selected = _numeric(rows, hm["plot_selected"])[mask] > 0
        xp, yp = x_all[selected], y_all[selected]
    else:
        xp, yp = x_all, y_all

    distance_type = str(metadata_row.get(hm.get("distance_type", ""), "absolute") or "absolute")
    source_y_max = 0.0
    if "plot_score_distance_y_max" in hm:
        try:
            source_y_max = float(metadata_row.get(hm["plot_score_distance_y_max"], 0.0))
        except (TypeError, ValueError):
            source_y_max = 0.0

    # Filter points before binning whenever a visible y-range is requested.
    # This avoids giant hidden bins being clipped into vertical bands.
    lower = args.y_min
    upper = args.y_max
    if lower is None and distance_type == "absolute":
        lower = 0.0
    if upper is None and source_y_max > 0:
        upper = source_y_max
    visible = np.isfinite(xp) & np.isfinite(yp)
    if lower is not None:
        visible &= yp >= float(lower)
    if upper is not None:
        visible &= yp <= float(upper)
    xp, yp = xp[visible], yp[visible]

    fig, ax = plt.subplots()
    render = str(metadata_row.get(hm.get("plot_type", ""), "hexbin") or "hexbin")
    if render == "scatter":
        kw = {"s": 7, "alpha": 0.35, "linewidths": 0, "rasterized": True}
        kw.update(artist_kw.get("points", {}))
        ax.scatter(xp, yp, **kw)
    else:
        kw = {"gridsize": 60, "mincnt": 1, "bins": "log", "rasterized": True}
        artist = ax.hexbin(xp, yp, **kw)
        fig.colorbar(artist, ax=ax).set_label("log10 plotted pair count")

    if compact_source and "full_linear_slope" in hm:
        try:
            slope = float(metadata_row.get(hm["full_linear_slope"], "nan"))
            intercept = float(metadata_row.get(hm["full_linear_intercept"], "nan"))
            rho = float(metadata_row.get(hm.get("full_spearman_rho", ""), "nan"))
            r2 = float(metadata_row.get(hm.get("full_linear_r_squared", ""), "nan"))
            n_full = int(float(metadata_row.get(hm.get("full_matched_pair_count", ""), len(x_all))))
        except (TypeError, ValueError):
            slope = intercept = rho = r2 = math.nan
            n_full = len(x_all)
    elif x_all.size >= 2 and np.nanmax(x_all) != np.nanmin(x_all):
        try:
            from scipy.stats import linregress, spearmanr
            reg = linregress(x_all, y_all)
            rho_result = spearmanr(x_all, y_all)
            slope, intercept = float(reg.slope), float(reg.intercept)
            rho, r2 = float(rho_result.statistic), float(reg.rvalue**2)
        except Exception:
            slope = intercept = rho = r2 = math.nan
        n_full = len(x_all)
    else:
        slope = intercept = rho = r2 = math.nan
        n_full = len(x_all)

    if x_all.size and math.isfinite(slope) and math.isfinite(intercept):
        endpoints = np.asarray([np.nanmin(x_all), np.nanmax(x_all)], dtype=float)
        line_kw = {"linestyle": ":", "linewidth": 1.2}
        line_kw.update(artist_kw.get("line", {}))
        ax.plot(endpoints, intercept + slope * endpoints, **line_kw)
    ax.text(
        0.02, 0.98,
        f"Spearman ρ = {rho:.3f}\nLinear R² = {r2:.3f}\nn = {n_full:,}",
        transform=ax.transAxes, va="top", ha="left",
    )

    main_label = str(metadata_row.get(hm.get("main_label", ""), metadata_row.get(hm.get("plot_label_a", ""), "Main")) or "Main")
    label = str(metadata_row.get(hm.get("comparison", ""), "comparison")) if rows else "comparison"
    ax.set_xlabel(f"{main_label} peak score")
    if distance_type == "signed":
        ax.set_ylabel(f"Signed {label} − {main_label} distance (bp)")
    else:
        ax.set_ylabel("Absolute matched distance (bp)")
    if lower is not None or upper is not None:
        bottom, top = ax.get_ylim()
        ax.set_ylim(bottom if lower is None else float(lower), top if upper is None else float(upper))
    ax.set_title(f"{main_label} peak score versus distance: {label}")
    return _finish(ax, fig, args, output, artist_kw=artist_kw, default_size=(8.0, 6.0)), fig

def _plot_compare_positions_correlation(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    comparison_key = hm.get("comparison")
    groups = []
    if comparison_key:
        for row in rows:
            label = str(row.get(comparison_key, ""))
            if label not in groups: groups.append(label)
    else:
        groups = [""]
    fig, ax = plt.subplots(); plotted = 0
    method = str(rows[0].get(hm.get("plot_correlation_method", ""), "spearman") or "spearman") if rows else "spearman"
    base_labels = []
    for group in groups:
        subset = rows if comparison_key is None else [row for row in rows if str(row.get(comparison_key, "")) == group]
        labels = [str(row.get(hm["distance_bin"], "")) for row in subset]
        if not base_labels: base_labels = labels
        x = np.arange(len(labels), dtype=float)
        for key, corr_label, marker, linestyle in (
            ("spearman_score_correlation", "Spearman", "o", "-"),
            ("pearson_score_correlation", "Pearson", "s", "--"),
        ):
            if key not in hm or (method != "both" and corr_label.lower() != method):
                continue
            values = _numeric(subset, hm[key])
            legend_label = corr_label if not group else (group if method != "both" else f"{group} {corr_label}")
            kw = {"marker": marker, "linestyle": linestyle, "label": legend_label}; kw.update(artist_kw.get("line", {}))
            ax.plot(x, values, **kw); plotted += 1
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(np.arange(len(base_labels), dtype=float), base_labels, rotation=35, ha="right"); ax.set_ylim(-1.05, 1.05)
    ax.set_xlabel("Absolute summit-distance bin (bp)"); ax.set_ylabel("Main/comparison score correlation")
    ax.set_title("Score correlation by summit-distance bin")
    return _finish(ax, fig, args, output, legend=plotted > 1 or bool(comparison_key), artist_kw=artist_kw, default_size=(9.0, 5.8)), fig


def _plot_compare_positions_percentile_boxplot(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    hm = _header_map(headers)
    group_key = hm.get("percentile_group")
    comparison_key = hm.get("comparison")
    row_type_key = hm.get("row_type")
    compact = (
        row_type_key is not None
        and hm.get("q1_absolute_distance") is not None
        and hm.get("whisker_low_absolute_distance") is not None
    )

    if group_key is None:
        raise ValueError("Percentile boxplots require a percentile_group column")

    # Preserve source order for percentile groups and comparison callsets.
    group_labels: list[str] = []
    comparison_labels: list[str] = []
    if compact:
        box_rows = [row for row in rows if str(row.get(row_type_key, "")).lower() == "box"]
        for row in box_rows:
            group = str(row.get(group_key, ""))
            comparison = str(row.get(comparison_key, "comparison")) if comparison_key else "comparison"
            if group and group not in group_labels:
                group_labels.append(group)
            if comparison and comparison not in comparison_labels:
                comparison_labels.append(comparison)
    else:
        distance_key = hm.get("absolute_distance")
        if distance_key is None:
            raise ValueError(
                "Percentile boxplots require either compact boxplot statistics or an absolute_distance column"
            )
        for row in rows:
            group = str(row.get(group_key, ""))
            comparison = str(row.get(comparison_key, "comparison")) if comparison_key else "comparison"
            if group and group not in group_labels:
                group_labels.append(group)
            if comparison and comparison not in comparison_labels:
                comparison_labels.append(comparison)

    if not group_labels or not comparison_labels:
        raise ValueError("No percentile boxplot groups were found")

    from nucleosuite.plotting import category_colors
    colors = category_colors(len(comparison_labels))
    fig, ax = plt.subplots()
    centres = np.arange(1, len(group_labels) + 1, dtype=float)
    if len(comparison_labels) <= 1:
        offsets = np.asarray([0.0])
        width = 0.48
    else:
        spread = min(0.70, 0.14 * len(comparison_labels))
        offsets = np.linspace(-spread / 2.0, spread / 2.0, len(comparison_labels))
        width = min(0.16, 0.70 / len(comparison_labels))

    # None means follow the compact source setting when available, otherwise
    # use the current NucleoSuite default (outliers shown).
    show_fliers = args.show_boxplot_outliers
    if show_fliers is None:
        show_fliers = True
        if compact and hm.get("show_boxplot_outliers"):
            for row in box_rows:
                raw = str(row.get(hm["show_boxplot_outliers"], "")).strip().lower()
                if raw in {"0", "false", "no", "off"}:
                    show_fliers = False
                    break
                if raw in {"1", "true", "yes", "on"}:
                    show_fliers = True
                    break

    positions: dict[tuple[str, str], float] = {}
    if compact:
        box_index = {
            (str(row.get(group_key, "")), str(row.get(comparison_key, "comparison"))): row
            for row in box_rows
        }
        flier_values: dict[tuple[str, str], list[float]] = defaultdict(list)
        outlier_key = hm.get("outlier_absolute_distance")
        if outlier_key:
            for row in rows:
                if str(row.get(row_type_key, "")).lower() != "flier":
                    continue
                group = str(row.get(group_key, ""))
                comparison = str(row.get(comparison_key, "comparison"))
                try:
                    value = float(row.get(outlier_key, "nan"))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(value):
                    flier_values[(group, comparison)].append(value)

        required = {
            "q1": hm["q1_absolute_distance"],
            "med": hm["median_absolute_distance"],
            "q3": hm["q3_absolute_distance"],
            "whislo": hm["whisker_low_absolute_distance"],
            "whishi": hm["whisker_high_absolute_distance"],
        }
        for ci, (comparison, color) in enumerate(zip(comparison_labels, colors)):
            for gi, group in enumerate(group_labels):
                row = box_index.get((group, comparison))
                if row is None:
                    continue
                try:
                    stats = {key: float(row.get(column, "nan")) for key, column in required.items()}
                except (TypeError, ValueError):
                    continue
                if not all(math.isfinite(value) for value in stats.values()):
                    continue
                stats["fliers"] = flier_values.get((group, comparison), [])
                pos = float(centres[gi] + offsets[ci])
                positions[(group, comparison)] = pos
                box = ax.bxp(
                    [stats], positions=[pos], widths=width, patch_artist=True,
                    showfliers=bool(show_fliers), manage_ticks=False,
                )
                box["boxes"][0].set_facecolor(color)
                box["boxes"][0].set_edgecolor(color)
                for key in ("whiskers", "caps"):
                    for artist in box[key]:
                        artist.set_color(color)
                for artist in box["medians"]:
                    artist.set_color("black")
                for artist in box.get("fliers", []):
                    artist.set_markeredgecolor(color)
                    artist.set_markerfacecolor("none")
                    artist.set_markersize(3.5)
                    artist.set_alpha(0.7)
    else:
        distance_key = hm["absolute_distance"]
        grouped = defaultdict(list)
        for row in rows:
            group = str(row.get(group_key, ""))
            comparison = str(row.get(comparison_key, "comparison")) if comparison_key else "comparison"
            try:
                value = float(row.get(distance_key, "nan"))
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                grouped[(group, comparison)].append(value)
        for ci, (comparison, color) in enumerate(zip(comparison_labels, colors)):
            for gi, group in enumerate(group_labels):
                values = grouped[(group, comparison)]
                if not values:
                    continue
                pos = float(centres[gi] + offsets[ci])
                positions[(group, comparison)] = pos
                box = ax.boxplot(
                    [values], positions=[pos], widths=width, patch_artist=True,
                    showfliers=bool(show_fliers), whis=1.5,
                )
                box["boxes"][0].set_facecolor(color)
                box["boxes"][0].set_edgecolor(color)
                for key in ("whiskers", "caps"):
                    for artist in box[key]:
                        artist.set_color(color)
                for artist in box["medians"]:
                    artist.set_color("black")
                for artist in box.get("fliers", []):
                    artist.set_markeredgecolor(color)
                    artist.set_markerfacecolor("none")
                    artist.set_markersize(3.5)
                    artist.set_alpha(0.7)

    ax.set_xticks(centres, group_labels, rotation=0)
    ax.set_xlim(0.45, len(group_labels) + 0.55)

    source_y_max = 200.0
    if compact and hm.get("percentile_boxplot_y_max"):
        for row in box_rows:
            try:
                candidate = float(row.get(hm["percentile_boxplot_y_max"], "nan"))
            except (TypeError, ValueError):
                continue
            if math.isfinite(candidate):
                source_y_max = candidate
                break
    if source_y_max > 0:
        ax.set_ylim(0.0, source_y_max)

    main_label = "Main"
    if hm.get("main_label"):
        for row in rows:
            value = str(row.get(hm["main_label"], "")).strip()
            if value:
                main_label = value
                break
    ax.set_xlabel(f"{main_label} peak score percentile group")
    ax.set_ylabel("Absolute matched distance (bp)")

    # Recreate within-percentile statistical brackets from compact source rows.
    max_levels = 0
    if compact and row_type_key:
        stat_rows = [row for row in rows if str(row.get(row_type_key, "")).lower() == "stat"]
        by_group: dict[str, list[dict[str, str]]] = {group: [] for group in group_labels}
        for row in stat_rows:
            by_group.setdefault(str(row.get(group_key, "")), []).append(row)
        for group in group_labels:
            group_stats = by_group.get(group, [])
            max_levels = max(max_levels, len(group_stats))
            for level, row in enumerate(group_stats):
                c1 = str(row.get(hm.get("comparison_1", ""), ""))
                c2 = str(row.get(hm.get("comparison_2", ""), ""))
                x1 = positions.get((group, c1))
                x2 = positions.get((group, c2))
                if x1 is None or x2 is None:
                    continue
                y = 1.03 + level * 0.065
                transform = ax.get_xaxis_transform()
                ax.plot(
                    [x1, x1, x2, x2], [y - 0.015, y, y, y - 0.015],
                    transform=transform, clip_on=False, color="black", linewidth=0.8,
                )
                display_mode = str(row.get(hm.get("p_display", ""), "value") or "value")
                significance = str(row.get(hm.get("significance", ""), ""))
                if display_mode == "stars" and significance:
                    text = significance
                else:
                    p_adjustment = str(row.get(hm.get("p_adjustment", ""), "none") or "none")
                    p_col = hm.get("p_adjusted") if p_adjustment == "holm" else hm.get("p_value")
                    try:
                        p_value = float(row.get(p_col, "nan")) if p_col else math.nan
                    except (TypeError, ValueError):
                        p_value = math.nan
                    if not math.isfinite(p_value):
                        text = "NA"
                    elif p_value < 0.0001:
                        text = "p_adj<1e-4" if p_adjustment == "holm" else "p<1e-4"
                    else:
                        text = ("p_adj=" if p_adjustment == "holm" else "p=") + f"{p_value:.3g}"
                ax.text(
                    (x1 + x2) / 2.0, y + 0.006, text, transform=transform,
                    ha="center", va="bottom", fontsize=7, clip_on=False,
                )

    title_y = 1.02 if max_levels == 0 else 1.09 + max_levels * 0.065
    ax.set_title(f"Matched distance by {main_label} score percentile", y=title_y)

    # Create the comparison legend exactly once. The previous replot path first
    # supplied proxy handles and then called ax.legend() a second time, which
    # caused the "No artists with labels" warning seen on percentile replots.
    if not args.no_legend:
        proxies = [Patch(facecolor=color, edgecolor=color) for color in colors]
        legend_kwargs = {"frameon": False}
        legend_kwargs.update(artist_kw.get("legend", {}))
        ax.legend(proxies, comparison_labels, **legend_kwargs)

    return _finish(
        ax, fig, args, output, legend=False, artist_kw=artist_kw,
        default_size=(max(9.0, 1.8 * len(group_labels) + 3.0), 6.5),
    ), fig

def _plot_distance_state_overlay(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers); state_key = hm["state"]; x_key = hm["distance_bp"]
    value_key = args.y_column or hm.get("raw_percent") or hm.get("smoothed_percent")
    states = sorted({str(row.get(state_key, "")) for row in rows}, key=_natural_key)
    fig, ax = plt.subplots()
    for state in states:
        state_rows = [row for row in rows if str(row.get(state_key, "")) == state]
        x = _numeric(state_rows, x_key); y = _numeric(state_rows, value_key)
        mask = np.isfinite(x) & np.isfinite(y); x, y = x[mask], y[mask]; order = np.argsort(x); x, y = x[order], y[order]
        colour = None
        if "rgb" in hm and state_rows:
            colour = _parse_rgb(str(state_rows[0].get(hm["rgb"], "")))
        kw = {"linewidth": 1.8, "marker": "o", "markersize": 1.8, "markeredgewidth": 0, "label": state}
        if colour is not None: kw["color"] = colour
        kw.update(artist_kw.get("line", {})); ax.plot(x, y, **kw)
    all_x = _numeric(rows, x_key); finite_x = all_x[np.isfinite(all_x)]
    if finite_x.size: ax.set_xlim(float(np.min(finite_x)), float(np.max(finite_x)))
    ax.set_xlabel("Adjacent peak distance (bp)"); ax.set_ylabel("Raw distance frequency within state (%)")
    ax.set_title("Adjacent peak distances by chromatin state")
    from nucleosuite.plotting import apply_distance_x_axis
    if args.x_major_grid is None and args.x_minor_grid is None:
        apply_distance_x_axis(ax, major_interval=args.x_major_tick, minor_interval=args.x_minor_tick); bp_x = False
    else: bp_x = True
    return _finish(
        ax, fig, args, output, bp_x=bp_x, legend=True,
        artist_kw={**artist_kw, "legend": {"ncol": 2, "fontsize": 8, "bbox_to_anchor": (1.02, 1), "loc": "upper left", **artist_kw.get("legend", {})}},
        default_size=(13.0, 8.0),
    ), fig


def _percentile_colour_scale(values: Sequence[float]):
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    colour_map = LinearSegmentedColormap.from_list(
        "nucleosuite_percentile", ("#ff1f1f", "#a00070", "#2020ff")
    )
    minimum = min(values); maximum = max(values)
    if maximum <= minimum: maximum = minimum + 1.0
    return colour_map, Normalize(vmin=minimum, vmax=maximum)


def _plot_distance_percentile_curves(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import FuncFormatter
    hm = _header_map(headers)
    percentile_key = hm["percentile_threshold"]; x_key = hm["distance_bp"]
    y_key = args.y_column or hm.get("raw_count") or hm.get("raw_percent")
    if y_key is None:
        raise ValueError("Percentile distance curves require raw_count or raw_percent")
    percentiles = sorted({float(row[percentile_key]) for row in rows})
    colour_map, normalizer = _percentile_colour_scale(percentiles)
    fig, ax = plt.subplots()
    for percentile in percentiles:
        subset = [row for row in rows if float(row[percentile_key]) == percentile]
        x = _numeric(subset, x_key); y = _numeric(subset, y_key)
        mask = np.isfinite(x) & np.isfinite(y); x, y = x[mask], y[mask]; order = np.argsort(x)
        kw = {"color": colour_map(normalizer(percentile)), "linewidth": 1.05}
        kw.update(artist_kw.get("line", {})); ax.plot(x[order], y[order], **kw)
    scalar = ScalarMappable(norm=normalizer, cmap=colour_map); scalar.set_array([])
    colourbar = fig.colorbar(scalar, ax=ax, pad=0.025)
    is_bins = "percentile_lower" in hm and any(str(row.get(hm["percentile_lower"], "")).strip() for row in rows)
    colourbar.set_label("Peak score percentile range" if is_bins else "Peak score threshold percentile")
    if len(percentiles) <= 12:
        colourbar.set_ticks(percentiles)
        if "percentile_range" in hm:
            labels = []
            for percentile in percentiles:
                match = next(row for row in rows if float(row[percentile_key]) == percentile)
                labels.append(str(match.get(hm["percentile_range"], "")) or f"{percentile:g}")
            colourbar.set_ticklabels(labels)
    is_percent = _clean_header(y_key) == "raw_percent"
    ax.set_xlabel("Distance (bp) between flanking peaks")
    ax.set_ylabel("Percentage of distances (%)" if is_percent else "Count")
    scope = str(rows[0].get(hm.get("scope", ""), "")) if rows else ""
    chromosome = str(rows[0].get(hm.get("chromosome", ""), "")) if rows else ""
    order_value = str(rows[0].get(hm.get("order", ""), "")) if rows else ""
    scope_label = "Combined chromosomes" if scope == "combined_chromosomes" else chromosome
    metric = "percentages" if is_percent else "counts"
    ax.set_title(f"Peak-distance {metric} by score threshold\n{scope_label}; neighbour order +{order_value}")
    ax.grid(axis="x", alpha=0.35)
    if not is_percent:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _position: f"{value / 1_000_000:g}M" if abs(value) >= 1_000_000 else f"{value / 1_000:g}K" if abs(value) >= 1_000 else f"{value:g}"))
    return _finish(ax, fig, args, output, bp_x=True, artist_kw=artist_kw, default_size=(10.5, 6.2)), fig


def _plot_distance_percentile_peak_counts(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers); percentile_key = hm["percentile_threshold"]; count_key = hm["retained_peak_count"]
    ordered = sorted(rows, key=lambda row: float(row[percentile_key]))
    percentiles = [float(row[percentile_key]) for row in ordered]
    is_bins = "percentile_lower" in hm and any(str(row.get(hm["percentile_lower"], "")).strip() for row in ordered)
    labels = [
        (str(row.get(hm.get("percentile_range", ""), "")) or f"{float(row[percentile_key]):g}")
        for row in ordered
    ]
    colour_map, normalizer = _percentile_colour_scale(percentiles)
    fig, ax = plt.subplots(figsize=(6.2, max(2.4, 0.34 * len(ordered) + 1.25))); ax.axis("off")
    table = ax.table(
        cellText=[["━━", label, f"{int(float(row[count_key])):,}"] for label, row in zip(labels, ordered)],
        colLabels=("", "Percentile range" if is_bins else "Percentile threshold", "Retained peak count"),
        cellLoc="right", colLoc="center", loc="center", colWidths=(0.12, 0.42, 0.46),
    )
    table.auto_set_font_size(False); table.set_fontsize(10); table.scale(1.0, 1.22)
    for row_index, percentile in enumerate(percentiles, start=1):
        table[(row_index, 0)].get_text().set_color(colour_map(normalizer(percentile)))
        table[(row_index, 0)].get_text().set_ha("center")
    for column in range(3): table[(0, column)].get_text().set_weight("bold")
    return _finish(
        ax, fig, args, output, artist_kw=artist_kw,
        default_size=(6.2, max(2.4, 0.34 * len(ordered) + 1.25)),
    ), fig


def _parse_rgb(value: str):
    try:
        parts=[int(x) for x in value.split(",")];
        if len(parts)==3: return tuple(x/255 for x in parts)
    except Exception: pass
    return None


def _natural_key(text: str):
    return [int(p) if p.isdigit() else p.casefold() for p in re.split(r"(\d+)", text)]


def _plot_flank_spacing(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    category_key = hm["category"]
    x_key = hm["spacing_bp"]
    y_key = args.y_column or hm["value"]
    rank_key = hm.get("rank")
    highlighted_key = hm.get("highlighted")

    categories: list[str] = []
    by_category: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        category = str(row.get(category_key, ""))
        if category not in categories:
            categories.append(category)
        by_category[category].append(row)

    def rank_of(category: str) -> int:
        subset = by_category[category]
        if rank_key and subset:
            try:
                return int(float(subset[0].get(rank_key, "nan")))
            except (TypeError, ValueError):
                pass
        return categories.index(category) + 1

    selected: list[str] = []
    for category in categories:
        subset = by_category[category]
        highlighted = False
        if highlighted_key and subset:
            raw = str(subset[0].get(highlighted_key, "")).strip().lower()
            highlighted = raw in {"1", "true", "yes", "on"}
        if highlighted:
            selected.append(category)
    selected.sort(key=rank_of)

    from nucleosuite.plotting import category_colors
    colors = category_colors(len(selected))
    color_by_category = {category: colors[index] for index, category in enumerate(selected)}

    fig, ax = plt.subplots()
    line_overrides = dict(artist_kw.get("line", {}))
    for category in categories:
        if category in color_by_category:
            continue
        subset = by_category[category]
        x = _numeric(subset, x_key)
        y = _numeric(subset, y_key)
        mask = np.isfinite(x) & np.isfinite(y)
        order = np.argsort(x[mask])
        kw = {"color": "0.72", "linewidth": 1.0, "alpha": 0.75, "zorder": 1}
        kw.update(line_overrides)
        ax.plot(x[mask][order], y[mask][order], **kw)

    # Lower-priority highlighted categories first; rank 1 is drawn last.
    handles: dict[str, object] = {}
    for category in reversed(selected):
        subset = by_category[category]
        x = _numeric(subset, x_key)
        y = _numeric(subset, y_key)
        mask = np.isfinite(x) & np.isfinite(y)
        order = np.argsort(x[mask])
        rank = rank_of(category)
        kw = {
            "color": color_by_category[category], "linewidth": 1.6,
            "label": category, "zorder": 10 + (len(selected) - rank),
        }
        kw.update(line_overrides)
        line, = ax.plot(x[mask][order], y[mask][order], **kw)
        handles[category] = line

    first = rows[0] if rows else {}
    distribution = str(first.get(hm.get("distribution", ""), "density") or "density")
    ratio_x1 = str(first.get(hm.get("ratio_x1", ""), "190") or "190")
    ratio_x2 = str(first.get(hm.get("ratio_x2", ""), "260") or "260")
    try:
        source_x_min = float(first.get(hm.get("x_min", ""), np.nan))
        source_x_max = float(first.get(hm.get("x_max", ""), np.nan))
    except (TypeError, ValueError):
        source_x_min = source_x_max = math.nan
    if math.isfinite(source_x_min) and math.isfinite(source_x_max):
        ax.set_xlim(source_x_min, source_x_max)
    ax.set_xlabel("Distance (bp) between flanking nucleosome centres")
    ax.set_ylabel("Density" if distribution == "density" else "Count")
    ax.set_title(f"Flanking nucleosome spacing ({ratio_x1}/{ratio_x2} ratio ranking)")
    ax.spines[["top", "right"]].set_visible(False)
    if selected and not args.no_legend:
        legend_kwargs = {"frameon": False, "title": "Category"}
        legend_kwargs.update(artist_kw.get("legend", {}))
        ax.legend([handles[c] for c in selected], selected, **legend_kwargs)
    return _finish(ax, fig, args, output, legend=False, artist_kw=artist_kw, default_size=(6.5, 4.6)), fig


def _plot_peak_states(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm=_header_map(headers)
    selection_key = hm.get("percentile_range", hm.get("percentile_threshold")); state_key=hm["state"]; value_key=args.y_column or hm["percentage_of_assigned_peaks"]
    labels=[]; states=sorted({r[state_key] for r in rows},key=_natural_key)
    for row in rows:
        label=row.get(selection_key,"") or row.get(hm.get("percentile_threshold",""),"")
        if label not in labels: labels.append(label)
    x=np.arange(len(labels),dtype=float); fig,ax=plt.subplots(); bottom=np.zeros(len(labels)); gap=min(max(args.bar_gap,0),0.95); width=1-gap
    for state in states:
        vals=[]; colour=None
        for label in labels:
            match=next((r for r in rows if r[state_key]==state and ((r.get(selection_key,"") or r.get(hm.get("percentile_threshold",""),""))==label)),None)
            vals.append(float(match.get(value_key,0)) if match else 0.0)
            if match and "rgb" in hm: colour=_parse_rgb(match.get(hm["rgb"],""))
        bkw={"width":width,"bottom":bottom,"label":state,"edgecolor":"none"}; bkw.update(artist_kw.get("bar", {}));
        if colour is not None and "color" not in bkw: bkw["color"]=colour
        ax.bar(x,vals,**bkw); bottom += np.asarray(vals)
    bins = "percentile_range" in hm and any(str(row.get(hm["percentile_range"], "")).strip() for row in rows)
    ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylim(0,100)
    ax.set_xlabel("Peak score percentile range" if bins else "Peak score threshold percentile")
    ax.set_ylabel("Assigned peaks in each chromatin state (%)")
    ax.set_title(
        "Chromatin-state composition across peak-score percentile ranges"
        if bins else "Chromatin-state composition across peak-score thresholds"
    )
    ax.grid(axis="y", alpha=0.3)
    return _finish(
        ax, fig, args, output, legend=True,
        artist_kw={
            **artist_kw,
            "legend": {
                "title": "Chromatin state", "bbox_to_anchor": (1.02, 1.0),
                "loc": "upper left", "fontsize": 8,
                **artist_kw.get("legend", {}),
            },
        },
        default_size=(max(9.5, min(16.0, 7.5 + 0.05 * len(x))), 6.8),
    ), fig


def _matrix_from_table(headers, rows) -> tuple[list[str], np.ndarray, np.ndarray]:
    if len(headers) < 2: raise ValueError("Heatmap tables require a row-label column plus numeric columns")
    x=[]; columns=[]
    for h in headers[1:]:
        try: x.append(float(h)); columns.append(h)
        except ValueError: continue
    if not columns: raise ValueError("Heatmap column headers must contain numeric x positions")
    names=[str(r.get(headers[0],i+1)) for i,r in enumerate(rows)]
    matrix=np.asarray([[float(r.get(h,"nan")) for h in columns] for r in rows],dtype=float)
    return names,np.asarray(x,dtype=float),matrix


def _fragment_heatmap_prefix(path: Path) -> Path:
    name = path.name
    for suffix in ("_normalised_matrix.tsv", "_normalized_matrix.tsv", "_normalised_matrix.tsv.gz", "_normalized_matrix.tsv.gz"):
        if name.lower().endswith(suffix):
            return path.with_name(name[: -len(suffix)])
    return path.with_name(path.stem)


def _fragment_heatmap_sidecars(path: Path):
    """Load optional fragment-heatmap presentation and clustering sidecars."""
    prefix = _fragment_heatmap_prefix(path)
    settings: dict[str, str] = {}
    setting_path = Path(f"{prefix}_heatmap_plot_metadata.tsv")
    if setting_path.is_file():
        _, setting_rows = _read_table(setting_path)
        settings = {
            str(row.get("setting", "")): str(row.get("value", ""))
            for row in setting_rows if str(row.get("setting", ""))
        }

    profile_rows: list[dict[str, str]] = []
    profile_path = Path(f"{prefix}_clustered_profiles.tsv")
    if profile_path.is_file():
        _, profile_rows = _read_table(profile_path)

    linkage_matrix = None
    linkage_path = Path(f"{prefix}_heatmap_linkage.tsv")
    if linkage_path.is_file():
        try:
            _, linkage_rows = _read_table(linkage_path)
        except ValueError:
            linkage_rows = []
        values = []
        for row in linkage_rows:
            try:
                values.append([
                    float(row["left_child"]), float(row["right_child"]),
                    float(row["distance"]), float(row["member_count"]),
                ])
            except (KeyError, TypeError, ValueError):
                continue
        if values:
            linkage_matrix = np.asarray(values, dtype=float)
    return settings, profile_rows, linkage_matrix


def _optional_float(settings: Mapping[str, str], key: str) -> float | None:
    text = str(settings.get(key, "")).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _plot_fragment_heatmap(path, names, x, matrix, args, output, artist_kw):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize, TwoSlopeNorm
    from matplotlib.patches import Patch
    from nucleosuite.fragment_heatmap import colour_limits, observed_category_colours

    settings, profile_rows, linkage_matrix = _fragment_heatmap_sidecars(path)
    profile_by_name = {str(row.get("profile", "")): row for row in profile_rows}

    # New fragment-heatmap outputs retain original row indices, allowing the
    # exact linkage matrix and its leaves to be reconstructed from the already
    # clustered matrix table.
    if linkage_matrix is not None and profile_rows and all(
        str(profile_by_name.get(name, {}).get("original_index", "")).strip() for name in names
    ):
        unclustered = np.empty_like(matrix)
        valid = True
        for clustered_index, name in enumerate(names):
            try:
                original_index = int(profile_by_name[name]["original_index"])
            except (KeyError, TypeError, ValueError):
                valid = False
                break
            if not 0 <= original_index < len(names):
                valid = False
                break
            unclustered[original_index] = matrix[clustered_index]
        if valid:
            matrix_for_order = unclustered
            from scipy.cluster.hierarchy import dendrogram
            order = np.asarray(dendrogram(linkage_matrix, no_plot=True)["leaves"], dtype=int)
            matrix = matrix_for_order[order]
            names = [next(name for name in names if int(profile_by_name[name]["original_index"]) == index) for index in order]
        else:
            linkage_matrix = None

    categories = [str(profile_by_name.get(name, {}).get("category", "")) for name in names]
    have_categories = any(categories)
    requested_colours = {
        key.split(":", 1)[1]: value
        for key, value in settings.items() if key.startswith("category_colour:")
    }
    category_map = observed_category_colours(categories, requested_colours)

    height = max(4.5, min(0.24 * len(names) + 2.5, 24.0))
    width = max(14.0, min(0.045 * len(x) + 12.0, 40.0))
    width = width if args.width is None else args.width
    height = height if args.height is None else args.height
    fig = plt.figure(figsize=(width, height))

    ratios = []
    if linkage_matrix is not None:
        ratios.append(2.2)
    if have_categories:
        ratios.append(0.18)
    label_gutter = _optional_float(settings, "label_gutter") or 1.2
    ratios.extend([label_gutter, 10.0, 0.38])
    grid = fig.add_gridspec(1, len(ratios), width_ratios=ratios, wspace=0.08)
    column = 0

    if linkage_matrix is not None:
        from scipy.cluster.hierarchy import dendrogram
        dendrogram_axis = fig.add_subplot(grid[0, column]); column += 1
        dendrogram_colour = settings.get("dendrogram_colour", "black") or "black"
        dendrogram(
            linkage_matrix, orientation="left", no_labels=True, color_threshold=0,
            above_threshold_color=dendrogram_colour,
            link_color_func=lambda _: dendrogram_colour, ax=dendrogram_axis,
        )
        dendrogram_axis.invert_yaxis(); dendrogram_axis.axis("off")

    if have_categories:
        category_axis = fig.add_subplot(grid[0, column]); column += 1
        keys = list(category_map)
        codes = {name: i for i, name in enumerate(keys)}
        values = np.asarray([codes.get(category, 0) for category in categories])
        category_axis.imshow(
            values.reshape(-1, 1), aspect="auto", interpolation="nearest", origin="upper",
            cmap=ListedColormap([category_map[key] for key in keys]),
            vmin=-0.5, vmax=max(len(keys) - 0.5, 0.5),
        )
        category_axis.set_xticks([]); category_axis.set_yticks([])
        category_axis.set_title("Category", fontsize=9, pad=6)

    gutter = fig.add_subplot(grid[0, column]); column += 1; gutter.axis("off")
    ax = fig.add_subplot(grid[0, column]); column += 1
    colourbar_axis = fig.add_subplot(grid[0, column])

    finite = matrix[np.isfinite(matrix)]
    centre = _optional_float(settings, "heatmap_centre")
    if "heatmap_centre" not in settings and finite.size and np.nanmin(finite) < 0 < np.nanmax(finite):
        centre = 0.0
    percentile = _optional_float(settings, "colour_percentile") or 99.0
    explicit_min = _optional_float(settings, "colour_min")
    explicit_max = _optional_float(settings, "colour_max")
    vmin, vmax = colour_limits(matrix, centre, percentile, explicit_min, explicit_max)
    if args.vmin is not None:
        vmin = args.vmin
    if args.vmax is not None:
        vmax = args.vmax
    if not vmin < vmax:
        raise ValueError("--vmin must be less than --vmax")
    norm = TwoSlopeNorm(vmin=vmin, vcenter=centre, vmax=vmax) if centre is not None and vmin < centre < vmax else Normalize(vmin=vmin, vmax=vmax)

    cmap_name = settings.get("heatmap_cmap", "").strip()
    cmap = (
        plt.get_cmap(cmap_name) if cmap_name else
        LinearSegmentedColormap.from_list(
            "blue_white_orange",
            [settings.get("low_colour", "#2166AC") or "#2166AC",
             settings.get("mid_colour", "#FFFFFF") or "#FFFFFF",
             settings.get("high_colour", "#F58518") or "#F58518"],
            N=256,
        )
    )
    heatmap_kw = {"aspect": "auto", "interpolation": "nearest", "origin": "upper", "cmap": cmap, "norm": norm}
    heatmap_kw.update(artist_kw.get("heatmap", {}))
    image = ax.imshow(matrix, **heatmap_kw)
    colourbar = fig.colorbar(image, cax=colourbar_axis)
    colourbar.set_label(settings.get("colourbar_label", "Heatmap value") or "Heatmap value")
    colourbar.ax.tick_params(labelsize=8)

    ax.set_xlabel(args.x_label or "Fragment length (bp)")
    max_yticks = int(_optional_float(settings, "max_yticks") or 80)
    step = 1 if len(names) <= max_yticks else int(np.ceil(len(names) / max_yticks))
    y_positions = np.arange(0, len(names), step)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([names[i] for i in y_positions], fontsize=8)
    ax.tick_params(axis="y", length=0, pad=6)
    first_tick = int(math.ceil(x[0] / 50.0) * 50)
    tick_values = np.arange(first_tick, x[-1] + 1, 50, dtype=int)
    ax.set_xticks(tick_values - x[0])
    ax.set_xticklabels([str(value) for value in tick_values], rotation=45, ha="right", fontsize=8)
    title = settings.get("title", "")
    if args.title is not None:
        title = args.title
    if args.no_title:
        title = ""
    if title:
        ax.set_title(title, pad=18)
    if args.x_min is not None or args.x_max is not None:
        left, right = ax.get_xlim(); ax.set_xlim(left if args.x_min is None else args.x_min - x[0], right if args.x_max is None else args.x_max - x[0])
    if args.y_min is not None or args.y_max is not None:
        bottom, top = ax.get_ylim(); ax.set_ylim(bottom if args.y_min is None else args.y_min, top if args.y_max is None else args.y_max)
    ax.grid(False)
    _configure_ticks_and_grids(ax, args, bp_x=False)
    for label in ax.get_xticklabels(): label.set_rotation(45 if args.x_tick_rotation is None else args.x_tick_rotation)
    if args.y_tick_rotation is not None:
        for label in ax.get_yticklabels(): label.set_rotation(args.y_tick_rotation)
    if args.axes_facecolor is not None: ax.set_facecolor(args.axes_facecolor)

    if have_categories and category_map and not args.no_legend:
        handles = [Patch(facecolor=colour, edgecolor="none", label=name) for name, colour in category_map.items()]
        legend_kw = {"loc": "lower left", "bbox_to_anchor": (0.0, 1.005), "ncol": min(4, len(handles)), "frameon": False, "fontsize": 9}
        legend_kw.update(artist_kw.get("legend", {}))
        ax.legend(handles=handles, **legend_kw)

    fig.text(0.01, 0.5, args.y_label or "Profile", rotation=90, va="center", fontsize=10)
    plt.subplots_adjust(left=0.06, right=0.95, top=0.91, bottom=0.12, wspace=0.08)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight", transparent=args.transparent)
    return output, fig


def _plot_heatmap(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    names,x,matrix=_matrix_from_table(headers,rows)
    is_fragment_heatmap = "normalised_matrix" in path.name.lower() or "normalized_matrix" in path.name.lower()
    if is_fragment_heatmap:
        return _plot_fragment_heatmap(path, names, x, matrix, args, output, artist_kw)
    fig,ax=plt.subplots(); finite=matrix[np.isfinite(matrix)]; automatic=float(np.nanmax(np.abs(finite))) if finite.size else 1.0
    # Keep symmetric automatic limits for matrices crossing zero; otherwise use data range.
    vmin = 0.0
    if args.vmin is None:
        if finite.size and np.nanmin(finite)<0<np.nanmax(finite): vmin=-automatic
        else: vmin=float(np.nanmin(finite)) if finite.size else 0.0
    vmax=args.vmax if args.vmax is not None else (automatic if finite.size and np.nanmin(finite)<0<np.nanmax(finite) else (float(np.nanmax(finite)) if finite.size else 1.0))
    if not vmin < vmax: raise ValueError("--vmin must be less than --vmax")
    default_size = (15.0, 5.0)
    hkw={"aspect":"auto","interpolation":"nearest","origin":"upper","vmin":vmin,"vmax":vmax,"extent":[x[0]-0.5,x[-1]+0.5,len(names),0],"cmap":"seismic"}
    hkw.update(artist_kw.get("heatmap", {})); image=ax.imshow(matrix,**hkw); colourbar=fig.colorbar(image,ax=ax)
    ax.set_xlabel("Position (bp)"); ax.set_ylabel("Region index (sorted)"); colourbar.set_label("Signal")
    if len(names)<=40: ax.set_yticks(np.arange(len(names))+0.5); ax.set_yticklabels(names,fontsize=max(5,args.font_size*0.72))
    ax.set_title(path.stem.replace("_"," "))
    return _finish(ax,fig,args,output,bp_x=True,artist_kw=artist_kw,default_size=default_size),fig



def _plot_multi_numeric_profile(
    path, headers, rows, args, output, artist_kw, *,
    x_key: str, columns: Sequence[str], xlabel: str, ylabel: str,
    bp_x: bool = True, linewidth: float = 1.25,
    markersize: float = 1.8, default_size: tuple[float, float] = (10.0, 5.5),
    vertical_zero: bool = False, exact_x_limits: bool = False,
    legend_defaults: Mapping[str, Any] | None = None,
):
    import matplotlib.pyplot as plt
    x = _numeric(rows, x_key)
    fig, ax = plt.subplots()
    plotted = 0
    for column in columns:
        y = _numeric(rows, column)
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            continue
        xx, yy = x[mask], y[mask]
        order = np.argsort(xx)
        kw = {"linewidth": linewidth, "marker": "o", "markersize": markersize, "markeredgewidth": 0, "label": column.rsplit("_", 1)[0]}
        kw.update(artist_kw.get("line", {}))
        ax.plot(xx[order], yy[order], **kw)
        plotted += 1
    if not plotted:
        raise ValueError(f"No numeric profile columns could be plotted from {path}")
    if vertical_zero:
        ax.axvline(0, linewidth=0.8, alpha=0.5)
    if exact_x_limits:
        finite_x = x[np.isfinite(x)]
        if finite_x.size:
            ax.set_xlim(float(np.min(finite_x)), float(np.max(finite_x)))
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(path.stem.replace("_", " "))
    configured_artist_kw = {key: dict(value) for key, value in artist_kw.items()}
    if legend_defaults:
        configured_artist_kw.setdefault("legend", {})
        configured_artist_kw["legend"] = {
            **legend_defaults,
            **configured_artist_kw["legend"],
        }
    return _finish(
        ax, fig, args, output, bp_x=bp_x, legend=plotted > 1,
        artist_kw=configured_artist_kw, default_size=default_size,
    ), fig


def _plot_dinucleotide_profile(path, headers, rows, args, output, artist_kw):
    hm = _header_map(headers)
    x_key = hm.get("position")
    if x_key is None:
        raise ValueError("Dinucleotide profiles require a position column")
    columns = [h for h in headers if re.fullmatch(r"[ACGT]{2}_(?:pct|frac)", h, re.IGNORECASE)]
    suffix = "fraction" if columns and columns[0].lower().endswith("_frac") else "percentage"
    return _plot_multi_numeric_profile(
        path, headers, rows, args, output, artist_kw,
        x_key=x_key, columns=columns, xlabel="Position relative to dyad (bp)",
        ylabel=f"Dinucleotide {suffix}", linewidth=1.1, markersize=1.8,
        default_size=(12.0, 5.0), vertical_zero=True,
        legend_defaults={"ncol": 4, "fontsize": "small"},
    )


def _plot_ww_ss_profile(path, headers, rows, args, output, artist_kw):
    hm = _header_map(headers)
    x_key = hm.get("position")
    if x_key is None:
        raise ValueError("WW/SS profiles require a position column")
    columns = [hm[name] for name in ("ww_pct", "ss_pct", "ww_frac", "ss_frac") if name in hm]
    suffix = "fraction" if columns and columns[0].lower().endswith("_frac") else "percentage"
    return _plot_multi_numeric_profile(
        path, headers, rows, args, output, artist_kw,
        x_key=x_key, columns=columns, xlabel="Position relative to dyad (bp)",
        ylabel=f"Dinucleotide {suffix}", linewidth=1.4, markersize=2.0,
        default_size=(12.0, 5.0), vertical_zero=True,
    )


def _plot_profile_overlay(path, headers, rows, args, output, artist_kw):
    hm = _header_map(headers)
    x_key = hm.get("relative_position")
    if x_key is None:
        raise ValueError("Profile overlays require a relative_position column")
    columns = [h for h in _first_numeric_columns(headers, rows) if h != x_key]
    return _plot_multi_numeric_profile(
        path, headers, rows, args, output, artist_kw,
        x_key=x_key, columns=columns, xlabel="Position relative to feature (bp)",
        ylabel="Mean signal", linewidth=1.4, markersize=1.8,
        default_size=(12.0, 5.0), vertical_zero=True, exact_x_limits=True,
    )


def _plot_ww_type_summary(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    category_key = hm.get("ww_type") or headers[0]
    count_key = hm.get("count") or next((h for h in _first_numeric_columns(headers, rows) if h != category_key), None)
    if count_key is None:
        raise ValueError("WW-type summaries require a count column")
    labels = _column(rows, category_key); values = _numeric(rows, count_key)
    mask = np.isfinite(values); labels = [label for label, keep in zip(labels, mask) if keep]; values = values[mask]
    fig, ax = plt.subplots()
    from nucleosuite.plotting import category_colors, apply_integer_y_axis
    kw = {"color": category_colors(len(labels))}; kw.update(artist_kw.get("bar", {})); ax.bar(np.arange(len(labels)), values, **kw)
    apply_integer_y_axis(ax)
    ax.set_xticks(np.arange(len(labels))); ax.set_xticklabels(labels, rotation=45); ax.set_xlabel(category_key.replace("_", " ").title()); ax.set_ylabel(count_key.replace("_", " ").title()); ax.set_title(path.stem.replace("_", " "))
    return _finish(ax, fig, args, output, artist_kw=artist_kw, default_size=(9.0, 5.0)), fig


def _plot_ww_type_by_length(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    hm = _header_map(headers)
    x_key = hm.get("fragment_length")
    if x_key is None:
        raise ValueError("WW-type length profiles require fragment_length")
    group_columns = [h for h in headers if _clean_header(h).endswith("_percent_of_classified")]
    x = _numeric(rows, x_key); mask = np.isfinite(x); x = x[mask]
    row_subset = [row for row, keep in zip(rows, mask) if keep]
    order = np.argsort(x); x = x[order]; row_subset = [row_subset[i] for i in order]
    positions = np.arange(len(x), dtype=float)
    fig, ax = plt.subplots(); bottom = np.zeros(len(x))
    for column in group_columns:
        vals = _numeric(row_subset, column)
        kw = {"bottom": bottom, "label": re.sub(r"_percent_of_classified$", "", column, flags=re.I)}; kw.update(artist_kw.get("bar", {}))
        ax.bar(positions, vals, **kw); bottom += np.nan_to_num(vals)
    ax.set_xticks(positions, [f"{value:g}" for value in x])
    ax.set_xlabel("Fragment length (bp)"); ax.set_ylabel("Relative frequency among classified fragments (%)"); ax.set_ylim(0, 100); ax.set_title(path.stem.replace("_", " "))
    ax.grid(axis="y", alpha=0.25)
    return _finish(
        ax, fig, args, output, legend=True,
        artist_kw={**artist_kw, "legend": {"ncol": 4, **artist_kw.get("legend", {})}},
        default_size=(max(7.0, 1.1 * len(x) + 4.0), 5.5),
    ), fig


def _plot_gene_sets_venn(path, headers, rows, args, output, artist_kw):
    import matplotlib.pyplot as plt
    try:
        from matplotlib_venn import venn2, venn3
    except ImportError as exc:
        raise RuntimeError("Gene-set Venn replots require matplotlib-venn") from exc
    hm = _header_map(headers); id_key = hm.get("gene_id"); sets_key = hm.get("candidate_sets")
    if id_key is None or sets_key is None:
        raise ValueError("Gene-set Venn replots require gene_id and candidate_sets columns")
    membership: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        gene = str(row.get(id_key, "")).strip()
        for name in filter(None, (part.strip() for part in str(row.get(sets_key, "")).split(","))):
            membership[name].add(gene)
    names = sorted(membership, key=_natural_key)
    if len(names) > 3:
        names = sorted(names, key=lambda name: (-len(membership[name]), _natural_key(name)))[:3]
    if len(names) not in {2, 3}:
        raise ValueError("A gene-set Venn replot requires two or three candidate sets")
    fig = plt.figure(); selected = [membership[name] for name in names]
    if len(names) == 2: venn2(selected, set_labels=names)
    else: venn3(selected, set_labels=names)
    ax = plt.gca(); ax.set_title("Candidate gene-set overlap")
    return _finish(ax, fig, args, output, artist_kw=artist_kw, default_size=(8.0, 8.0)), fig


def _plot_count_profile(path, headers, rows, args, output, artist_kw):
    hm = _header_map(headers)
    x_candidates = ("relocation_bp", "fragment_length", "distance_bp", "distance", "run_length_bp")
    x_key = args.x_column or next((hm[name] for name in x_candidates if name in hm), None)
    y_key = args.y_column or hm.get("count")
    if x_key is None or y_key is None:
        raise ValueError("Count profiles require a numeric coordinate and count column")
    title = None
    xlabel = x_key.replace("_", " ").title()
    if _clean_header(x_key) == "relocation_bp":
        xlabel = "Fragment relocation (bp)"
        title = "Randomized fragment relocation distances"
    saved, fig = _grouped_line(
        path, headers, rows, args, output, artist_kw,
        x_candidates=(_clean_header(x_key),), y_candidates=(_clean_header(y_key),),
        xlabel=xlabel, ylabel="Count", bp_x=True,
        line_defaults={"linewidth": 1.3, "markersize": 2.2},
        default_size=(10.0, 5.0), default_title=title,
    )
    ax = fig.axes[0]
    if _clean_header(x_key) == "relocation_bp":
        ax.axvline(0, linewidth=0.8, alpha=0.5)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(output, dpi=args.dpi, bbox_inches="tight", transparent=args.transparent)
    return saved, fig


def _generic_plot(path, headers, rows, args, output, artist_kw, kind: str):
    import matplotlib.pyplot as plt
    numeric=_first_numeric_columns(headers,rows); x_key=args.x_column or (numeric[0] if numeric else None); y_key=args.y_column or (numeric[1] if len(numeric)>1 else None)
    if x_key is None or y_key is None: raise ValueError("Generic plots require two numeric columns; specify --x-column and --y-column")
    x=_numeric(rows,x_key); y=_numeric(rows,y_key); mask=np.isfinite(x)&np.isfinite(y); x,y=x[mask],y[mask]; fig,ax=plt.subplots()
    if kind=="generic-scatter":
        kw={"s":20}; kw.update(artist_kw.get("points",{})); ax.scatter(x,y,**kw)
    elif kind=="generic-bar":
        kw={}; kw.update(artist_kw.get("bar",{})); ax.bar(x,y,**kw)
    else:
        kw={"linewidth":1.4,"marker":"o","markersize":2.2,"markeredgewidth":0}; kw.update(artist_kw.get("line",{})); ax.plot(x,y,**kw)
    ax.set_xlabel(x_key); ax.set_ylabel(y_key); ax.set_title(path.stem.replace("_"," "))
    return _finish(ax,fig,args,output,bp_x=args.bp_x,artist_kw=artist_kw),fig


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite plot",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", type=Path, help="NucleoSuite TSV/TSV.GZ output used to recreate a figure.")
    parser.add_argument("--plot-type", choices=PLOT_TYPES, default="auto", help="Plot family; default auto-detects from filename and columns.")
    parser.add_argument("--from-command", choices=tuple(sorted(COMMAND_TYPES)), help="Identify the NucleoSuite command that produced the input when auto-detection is ambiguous.")
    parser.add_argument("--output", type=Path, help="Output image; default <input>_replot.<format>.")
    parser.add_argument("--format", choices=("png", "svg", "pdf"), default="png", help="Output image format (default: png).")
    parser.add_argument("--width", type=float, default=None, help="Figure width in inches. Default: 10, or 6.5 for regression plots.")
    parser.add_argument("--height", type=float, default=None, help="Figure height in inches. Default: 5.5, or 6.5 for regression plots.")
    parser.add_argument("--dpi", type=int, default=300, help="Raster resolution in dots per inch (default: 300).")
    parser.add_argument("--title", help="Override the automatically generated plot title.")
    parser.add_argument("--no-title", action="store_true", help="Remove the plot title.")
    parser.add_argument("--x-label", help="Override the x-axis label.")
    parser.add_argument("--y-label", help="Override the y-axis label.")
    parser.add_argument("--font-size", type=float, default=10.0, help="Base font size in points (default: 10).")
    parser.add_argument("--x-min", type=float, help="Displayed x-axis minimum.")
    parser.add_argument("--x-max", type=float, help="Displayed x-axis maximum.")
    parser.add_argument("--y-min", type=float, help="Displayed y-axis minimum.")
    parser.add_argument("--y-max", type=float, help="Displayed y-axis maximum.")
    parser.add_argument("--x-major-tick", type=float, help="Major x tick increment.")
    parser.add_argument("--x-minor-tick", type=float, help="Minor x tick increment; minor tick labels remain hidden.")
    parser.add_argument("--y-major-tick", type=float, help="Major y tick increment.")
    parser.add_argument("--y-minor-tick", type=float, help="Minor y tick increment; minor tick labels remain hidden.")
    for axis in ("x", "y"):
        for which in ("major", "minor"):
            group = parser.add_mutually_exclusive_group()
            group.add_argument(
                f"--{axis}-{which}-grid", dest=f"{axis}_{which}_grid",
                action="store_true", default=None,
                help=f"Show {which} {axis}-axis grid lines.",
            )
            group.add_argument(
                f"--no-{axis}-{which}-grid", dest=f"{axis}_{which}_grid",
                action="store_false",
                help=f"Hide {which} {axis}-axis grid lines.",
            )
    parser.add_argument("--major-grid-color", default="0.65", help="Matplotlib color for major grid lines (default: 0.65).")
    parser.add_argument("--minor-grid-color", default="0.65", help="Matplotlib color for minor grid lines (default: 0.65).")
    parser.add_argument("--major-grid-alpha", type=float, default=0.65, help="Major grid-line opacity (default: 0.65).")
    parser.add_argument("--minor-grid-alpha", type=float, default=0.55, help="Minor grid-line opacity (default: 0.55).")
    parser.add_argument("--major-grid-width", type=float, default=0.75, help="Major grid-line width in points (default: 0.75).")
    parser.add_argument("--minor-grid-width", type=float, default=0.65, help="Minor grid-line width in points (default: 0.65).")
    parser.add_argument("--major-grid-style", default="-", help="Matplotlib linestyle for major grids (default: -).")
    parser.add_argument("--minor-grid-style", default="--", help="Matplotlib linestyle for minor grids (default: --).")
    parser.add_argument("--x-tick-rotation", type=float, default=None, help="Rotate x tick labels by this many degrees; defaults to the source command style.")
    parser.add_argument("--y-tick-rotation", type=float, default=None, help="Rotate y tick labels by this many degrees; defaults to the source command style.")
    parser.add_argument("--axes-facecolor", help="Matplotlib color for the axes background.")
    parser.add_argument("--transparent", action="store_true", help="Save with a transparent figure background.")
    parser.add_argument("--no-legend", action="store_true", help="Hide the legend.")
    boxplot_outliers = parser.add_mutually_exclusive_group()
    boxplot_outliers.add_argument(
        "--show-boxplot-outliers", dest="show_boxplot_outliers", action="store_true", default=None,
        help="Show boxplot outlier points beyond the 1.5×IQR whiskers (default: source setting, otherwise shown).",
    )
    boxplot_outliers.add_argument(
        "--hide-boxplot-outliers", dest="show_boxplot_outliers", action="store_false",
        help="Hide boxplot outlier points while retaining the standard 1.5×IQR whiskers.",
    )
    parser.add_argument("--x-column", help="Explicit x column for generic or ambiguous tables.")
    parser.add_argument("--y-column", help="Explicit y column for generic or ambiguous tables.")
    parser.add_argument("--group-column", help="Column used to split a table into multiple plotted series.")
    parser.add_argument(
        "--normalization", "--normalisation", dest="normalization",
        choices=("auto", "count", "density"), default="auto",
        help=(
            "Fragment-length y normalization. Auto reproduces the source command: "
            "density for fragment-lengths and counts for fragments (default: auto)."
        ),
    )
    parser.add_argument("--bp-x", action="store_true", help="Treat a generic x axis as base pairs and use bp tick defaults.")
    parser.add_argument("--smooth-window", type=int, default=21, help="DAC display-smoothing window. NRL-style DAC peak detection uses --peak-resolution.")
    parser.add_argument("--show-raw", action="store_true", help="Retain the raw DAC profile; raw DAC is shown by default.")
    parser.add_argument("--detect-peaks", action="store_true", help="Detect and display DAC peaks using the same resolution-driven method as `nucleosuite nrl` (default: off).")
    parser.add_argument("--peak-resolution", type=float, default=160.0, help="Peak resolution for opt-in DAC peak detection; also controls NRL-style smoothing windows (default: 160 bp).")
    parser.add_argument("--peak-min-distance", type=float, default=100.0, help="Display peak-spacing setting; DAC detection uses --peak-resolution.")
    parser.add_argument("--label-peaks", choices=("auto", "none", "peaks"), default="auto", help="Peak labels: auto follows the source command (shown for NRL, hidden for distance distributions).")
    parser.add_argument("--peak-label-value", choices=("x", "y", "both"), default="x", help="Value written for peak labels (default: x).")
    parser.add_argument("--peak-label-offset", type=float, default=5.0, help="Vertical peak-label offset in points (default: 5).")
    parser.add_argument("--nrl-inset", choices=("auto", "on", "off"), default="auto", help="Embed peak-number versus distance regression in DAC/NRL plots; auto uses it for long-range profiles.")
    parser.add_argument("--nrl-min-distance", type=float, default=100.0, help="Minimum peak distance included in an auto NRL inset (default: 100 bp).")
    parser.add_argument("--inset-bounds", type=float, nargs=4, metavar=("X", "Y", "W", "H"), default=(0.70, 0.57, 0.27, 0.36), help="Inset position as axes fractions X Y WIDTH HEIGHT.")
    parser.add_argument("--vmin", type=float, help="Heatmap lower saturation value.")
    parser.add_argument("--vmax", type=float, help="Heatmap upper saturation value.")
    parser.add_argument("--bar-gap", type=float, default=0.18, help="Fractional gap between peak-state stacked bars (default: 0.18).")
    parser.add_argument("--mpl-rc", action="append", default=[], metavar="KEY=VALUE", help="Set any Matplotlib rcParam; repeatable.")
    parser.add_argument("--mpl-kw", action="append", default=[], metavar="TARGET.KEY=VALUE", help="Pass Matplotlib artist keyword arguments. TARGET is line, raw, smooth, points, bar, heatmap, or legend; repeatable.")
    return parser

def validate_argv(argv: Sequence[str] | None = None) -> None:
    args=build_parser().parse_args(argv)
    if ((args.width is not None and args.width <= 0)
            or (args.height is not None and args.height <= 0)
            or args.dpi <= 0 or args.font_size <= 0):
        raise ValueError("Figure dimensions, DPI, and font size must be greater than zero")
    if args.smooth_window<1 or args.smooth_window%2==0: raise ValueError("--smooth-window must be a positive odd integer")
    if args.peak_resolution < 0: raise ValueError("--peak-resolution must be zero or greater")
    if args.vmin is not None and args.vmax is not None and args.vmin>=args.vmax: raise ValueError("--vmin must be less than --vmax")


def main(argv: Sequence[str] | None = None) -> int:
    args=build_parser().parse_args(argv); validate_argv(argv)
    import matplotlib
    matplotlib.use("Agg")
    _apply_rc(args.mpl_rc); matplotlib.rcParams["font.size"] = args.font_size
    headers,rows=_read_table(args.input); detected=detect_plot_type(args.input,headers)
    plot_type=args.plot_type
    if plot_type=="auto": plot_type=COMMAND_TYPES.get(args.from_command,detected) if args.from_command else detected
    # Specific table families should win over broad source-command hints.
    if args.plot_type=="auto" and (
        detected in {
            "nrl-regression", "fragment-size-nrl-profile",
            "fragment-size-nrl-regression", "heatmap", "distance-state-overlay",
            "aggregate-nrl-profile", "aggregate-nrl-regression",
            "distance-percentile-curves", "distance-percentile-peak-counts",
            "gene-expression-spacing", "gene-expression-fft-trajectory",
            "gene-expression-spacing-scatter", "gene-expression-ranking",
        }
        or detected.startswith("compare-positions-")
    ):
        plot_type=detected
    output=_resolve_output(args.input,args.output,args.format); artist_kw=_parse_artist_values(args.mpl_kw)
    dispatch={
        "dac":_plot_dac,"dcc":_plot_dcc,"nrl-profile":_plot_nrl_profile,"nrl-regression":_plot_nrl_regression,
        "fragment-size-nrl-profile":_plot_fragment_size_nrl_profile,
        "fragment-size-nrl-regression":_plot_fragment_size_nrl_regression,
        "aggregate-nrl-profile":_plot_aggregate_nrl_profile,
        "aggregate-nrl-regression":_plot_aggregate_nrl_regression,
        "distances":_plot_distances,"distance-state-overlay":_plot_distance_state_overlay,
        "distance-percentile-curves":_plot_distance_percentile_curves,
        "distance-percentile-peak-counts":_plot_distance_percentile_peak_counts,
        "aggregate-profile":_plot_aggregate_profile,"heatmap":_plot_heatmap,"generic-heatmap":_plot_heatmap,
        "fragment-lengths":_plot_fragment_lengths,"positive-runs":_plot_positive_runs,"peak-score-frequency":_plot_peak_score_frequency,
        "peak-states":_plot_peak_states,"flank-spacing":_plot_flank_spacing,"compare-positions":_plot_compare_positions,
        "compare-positions-histogram":_plot_compare_positions,
        "compare-positions-score":_plot_compare_positions_score,
        "compare-positions-correlation":_plot_compare_positions_correlation,
        "compare-positions-percentile-boxplot":_plot_compare_positions_percentile_boxplot,
        "compare-positions-score-distance":_plot_compare_positions_score_distance,
        "gene-expression":_plot_gene_expression,
        "gene-expression-spacing":_plot_gene_expression_spacing,
        "gene-expression-spacing-scatter":_plot_gene_expression_spacing_scatter,
        "gene-expression-fft-trajectory":_plot_gene_expression_fft_trajectory,
        "gene-expression-ranking":_plot_gene_expression_ranking,
        "tss-expression":_plot_tss_expression,
        "dinucleotide-profile":_plot_dinucleotide_profile,"ww-ss-profile":_plot_ww_ss_profile,"ww-type-summary":_plot_ww_type_summary,
        "ww-type-by-length":_plot_ww_type_by_length,"gene-sets-venn":_plot_gene_sets_venn,"profile-overlay":_plot_profile_overlay,"count-profile":_plot_count_profile,
    }
    if plot_type in dispatch: saved,fig=dispatch[plot_type](args.input,headers,rows,args,output,artist_kw)
    elif plot_type in {"generic-line","generic-scatter","generic-bar"}: saved,fig=_generic_plot(args.input,headers,rows,args,output,artist_kw,plot_type)
    else: raise ValueError(f"Unsupported plot type: {plot_type}")
    from nucleosuite.plotting import write_plot_metadata
    write_plot_metadata(saved, extra={"detected_plot_type": plot_type, "source_table": str(args.input)})
    import matplotlib.pyplot as plt; plt.close(fig)
    print(f"Detected plot type: {plot_type}"); print(f"Wrote: {saved}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
