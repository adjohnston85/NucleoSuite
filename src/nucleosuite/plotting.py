"""Shared plotting helpers for numeric genomic-distance axes."""

from __future__ import annotations

import argparse
import math
from typing import Iterable, Optional, Tuple

import numpy as np



def _finite_values(values: Iterable[float] | None) -> np.ndarray:
    """Return finite numeric values as a one-dimensional array."""

    if values is None:
        return np.asarray([], dtype=float)
    array = np.asarray(list(values), dtype=float).reshape(-1)
    return array[np.isfinite(array)]


def apply_integer_x_axis(axis, values: Iterable[float] | None = None) -> None:
    """Use integer major ticks for a discrete genomic x-axis.

    A one-value series receives a small integer-width view so the point and its
    coordinate remain visible instead of being expanded to decimal tick labels.
    """

    from matplotlib.ticker import MaxNLocator, StrMethodFormatter

    finite = _finite_values(values)
    if finite.size and float(np.min(finite)) == float(np.max(finite)):
        centre = int(round(float(finite[0])))
        axis.set_xlim(centre - 1, centre + 1)
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))


def _nice_base_pair_interval(span: float, target_intervals: int = 8) -> float:
    """Return an automatic base-pair tick interval divisible by 10.

    Intervals follow the restrained 1/2/5 progression (10, 20, 50, 100,
    200, 500, ...), preventing automatic labels such as 8, 16, 24 bp.
    """

    span = abs(float(span))
    if not math.isfinite(span) or span <= 0:
        return 10.0
    raw = max(10.0, span / max(1, int(target_intervals)))
    magnitude = 10.0 ** math.floor(math.log10(raw))
    scaled = raw / magnitude
    for step in (1.0, 2.0, 5.0, 10.0):
        if scaled <= step:
            interval = step * magnitude
            break
    else:  # pragma: no cover - the final step always catches finite input
        interval = 10.0 * magnitude
    return float(max(10.0, round(interval / 10.0) * 10.0))


def apply_base_pair_x_axis(axis, values: Iterable[float] | None = None) -> float:
    """Use automatic major ticks at multiples of 10 on a base-pair axis."""

    from matplotlib.ticker import MultipleLocator, StrMethodFormatter

    finite = _finite_values(values)
    if finite.size:
        minimum = float(np.min(finite))
        maximum = float(np.max(finite))
    else:
        minimum, maximum = map(float, axis.get_xlim())

    if math.isclose(minimum, maximum, rel_tol=1e-12, abs_tol=1e-12):
        centre = float(minimum)
        axis.set_xlim(centre - 5.0, centre + 5.0)
        span = 10.0
    else:
        span = maximum - minimum

    interval = _nice_base_pair_interval(span)
    axis.xaxis.set_major_locator(MultipleLocator(interval))
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    return interval


def apply_integer_y_axis(axis) -> None:
    """Use integer major ticks for an axis that represents raw counts."""

    from matplotlib.ticker import MaxNLocator, StrMethodFormatter

    axis.yaxis.set_major_locator(MaxNLocator(integer=True))
    axis.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))


def discrete_line_kwargs(*, markersize: float = 2.2) -> dict[str, object]:
    """Return restrained point-marker settings for discrete line profiles."""

    return {
        "marker": "o",
        "markersize": float(markersize),
        "markeredgewidth": 0.0,
    }



def category_colors(count: int):
    """Return ``count`` distinct category colours for one plot.

    The palette is sampled without cycling so plots with more than Matplotlib's
    default ten categories (for example all 16 dinucleotides) never reuse a
    colour within the same figure.
    """
    if count <= 0:
        return []
    import matplotlib.pyplot as plt
    if count <= 20:
        cmap = plt.get_cmap("tab20")
        return [cmap(index) for index in range(count)]
    if count <= 60:
        colours = []
        for name in ("tab20", "tab20b", "tab20c"):
            cmap = plt.get_cmap(name)
            colours.extend(cmap(index) for index in range(20))
        return colours[:count]
    cmap = plt.get_cmap("turbo")
    denominator = max(1, count - 1)
    return [cmap(index / denominator) for index in range(count)]


def configure_unique_category_cycle(count: int = 60) -> None:
    """Set a non-repeating categorical line/bar cycle for the current process."""
    import matplotlib as mpl
    from matplotlib import cycler

    mpl.rcParams["axes.prop_cycle"] = cycler(color=category_colors(count))


def validate_tick_interval(value: Optional[float], option: str) -> Optional[float]:
    """Validate an optional positive tick interval."""

    if value is None:
        return None
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{option} must be a finite value greater than zero.")
    return value


def infer_major_tick_interval(axis) -> float:
    """Infer the current positive major-tick spacing from a Matplotlib axis."""

    ticks = np.asarray(axis.get_xticks(), dtype=float)
    differences = np.diff(ticks)
    differences = differences[np.isfinite(differences) & (differences > 0)]
    if differences.size == 0:
        return 1.0
    return float(np.median(differences))


def default_minor_tick_interval(major_interval: float) -> float:
    """Return the default minor spacing requested for distance plots.

    A 10 bp major interval receives 5 bp minor ticks. Major intervals greater
    than 50 bp receive 10 bp minor ticks. Other intervals are divided in half.
    """

    major_interval = float(major_interval)
    if not math.isfinite(major_interval) or major_interval <= 0:
        raise ValueError("Major tick interval must be finite and greater than zero.")
    if math.isclose(major_interval, 10.0, rel_tol=1e-9, abs_tol=1e-9):
        return 5.0
    if major_interval > 50.0:
        return 10.0
    return major_interval / 2.0


def apply_distance_x_axis(
    axis,
    *,
    major_interval: Optional[float] = None,
    minor_interval: Optional[float] = None,
    major_grid_alpha: float = 0.5,
    minor_grid_alpha: float = 0.5,
) -> Tuple[float, float]:
    """Configure major/minor ticks and vertical grid lines on a distance axis.

    When ``major_interval`` is omitted, a readable 10/20/50 progression is
    selected so labelled base-pair ticks remain multiples of 10. When
    ``minor_interval`` is omitted, it is derived from the major spacing using
    :func:`default_minor_tick_interval`.
    """

    from matplotlib.ticker import MultipleLocator, StrMethodFormatter

    major_interval = validate_tick_interval(major_interval, "--x-major-tick")
    minor_interval = validate_tick_interval(minor_interval, "--x-minor-tick")

    if major_interval is not None:
        axis.xaxis.set_major_locator(MultipleLocator(major_interval))
        resolved_major = major_interval
    else:
        x0, x1 = axis.dataLim.intervalx
        if not (math.isfinite(x0) and math.isfinite(x1)):
            x0, x1 = axis.get_xlim()
        if math.isclose(x0, x1, rel_tol=1e-12, abs_tol=1e-12):
            centre = float(x0)
            axis.set_xlim(centre - 5.0, centre + 5.0)
            span = 10.0
        else:
            span = float(x1 - x0)
        resolved_major = _nice_base_pair_interval(span)
        axis.xaxis.set_major_locator(MultipleLocator(resolved_major))
    if math.isclose(resolved_major, round(resolved_major), rel_tol=1e-9, abs_tol=1e-9):
        axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))

    resolved_minor = (
        minor_interval
        if minor_interval is not None
        else default_minor_tick_interval(resolved_major)
    )
    axis.xaxis.set_minor_locator(MultipleLocator(resolved_minor))
    axis.tick_params(axis="x", which="minor", length=3)
    axis.set_axisbelow(True)
    axis.grid(
        axis="x",
        which="major",
        color="0.5",
        alpha=float(major_grid_alpha),
        linestyle="-",
        linewidth=0.8,
    )
    axis.grid(
        axis="x",
        which="minor",
        color="0.5",
        alpha=float(minor_grid_alpha),
        linestyle=":",
        linewidth=0.7,
    )
    return float(resolved_major), float(resolved_minor)

# ---------------------------------------------------------------------------
# Shared user-configurable plotting interface
# ---------------------------------------------------------------------------

from dataclasses import dataclass, replace
import os
from pathlib import Path
from typing import Any, Sequence


_POINT_MARKERS = {
    "circle": "o",
    "square": "s",
    "triangle": "^",
    "diamond": "D",
}

_LEGEND_LOCATIONS = {
    "best": "best",
    "upper-right": "upper right",
    "upper-left": "upper left",
    "lower-right": "lower right",
    "lower-left": "lower left",
    "outside-right": "center left",
}


@dataclass(frozen=True)
class PlotOptions:
    """Optional overrides shared by every NucleoSuite figure.

    ``None`` means that the plotting function keeps its command-specific
    default.  The file format is the one deliberate exception: NucleoSuite
    Generated figures use PNG unless ``--plot-format svg``
    is requested.
    """

    format: str = "png"
    width: float | None = None
    height: float | None = None
    dpi: int | None = None
    title: str | None = None
    no_title: bool = False
    x_label: str | None = None
    y_label: str | None = None
    font_size: float | None = None
    grid: str | None = None
    grid_color: str | None = None
    grid_alpha: float | None = None
    grid_width: float | None = None
    x_min: float | None = None
    x_max: float | None = None
    y_min: float | None = None
    y_max: float | None = None
    line_width: float | None = None
    line_color: str | None = None
    fill_color: str | None = None
    points: bool | None = None
    point_size: float | None = None
    point_fill: str | None = None
    point_edge: str | None = None
    point_edge_width: float | None = None
    point_shape: str | None = None
    label_points: str = "none"
    point_label_value: str = "x"
    point_label_offset: float | None = None
    legend: bool | None = None
    legend_position: str | None = None
    x_tick_rotation: float | None = None
    y_tick_rotation: float | None = None
    transparent: bool = False


_CURRENT_PLOT_OPTIONS = PlotOptions()


def _parser_has_option(parser, option: str) -> bool:
    return any(option in action.option_strings for action in parser._actions)


def _add_if_missing(parser, *flags: str, **kwargs):
    if any(_parser_has_option(parser, flag) for flag in flags):
        return None
    return parser.add_argument(*flags, **kwargs)


class _ExpandedPlottingHelpAction(argparse.Action):
    """Print command help with the normally hidden plot controls expanded."""

    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, **kwargs):
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            default=default,
            **kwargs,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        for action, original_help in getattr(parser, "_nucleosuite_plotting_help", []):
            action.help = original_help
        parser.print_help()
        parser.exit()


def _hide_plotting_help_actions(parser) -> None:
    """Hide ``--plot-*`` controls while retaining their parser behaviour."""
    stored = list(getattr(parser, "_nucleosuite_plotting_help", []))
    seen = {id(action) for action, _help in stored}
    for action in parser._actions:
        if id(action) in seen or not action.option_strings:
            continue
        if any(flag.startswith("--plot-") for flag in action.option_strings):
            stored.append((action, action.help))
            action.help = argparse.SUPPRESS
            seen.add(id(action))
    parser._nucleosuite_plotting_help = stored


def add_plotting_arguments(
    parser, *, label_points_default: str = "none", show_in_help: bool = False
) -> None:
    """Add the shared plotting options to an argparse parser.

    Plot customization is hidden from ordinary command help so analysis help
    stays compact.  ``--help-plotting`` expands the complete plot controls.
    The dedicated ``nucleosuite plot`` command can request visible plot help
    by passing ``show_in_help=True``.
    """

    group = parser.add_argument_group("plot customization")
    _add_if_missing(
        group, "--plot-format", choices=("png", "svg"), default="png",
        help="Figure format for generated plots (default: png).",
    )
    _add_if_missing(group, "--plot-width", type=float, help="Figure width in inches.")
    _add_if_missing(group, "--plot-height", type=float, help="Figure height in inches.")
    _add_if_missing(group, "--plot-dpi", type=int, help="Raster resolution in dots per inch; relevant mainly to PNG output.")
    _add_if_missing(group, "--plot-title", help="Override the automatically generated plot title.")
    _add_if_missing(group, "--no-plot-title", action="store_true", help="Remove plot titles.")
    _add_if_missing(group, "--plot-x-label", help="Override the x-axis label.")
    _add_if_missing(group, "--plot-y-label", help="Override the y-axis label.")
    _add_if_missing(group, "--plot-font-size", type=float, help="Base plot font size in points.")
    _add_if_missing(
        group, "--plot-grid", choices=("none", "x", "y", "both"),
        help="Show tick-aligned grid lines on neither axis, x only, y only, or both axes.",
    )
    _add_if_missing(group, "--plot-grid-color", help="Matplotlib color for grid lines.")
    _add_if_missing(group, "--plot-grid-alpha", type=float, help="Grid-line opacity from 0 to 1.")
    _add_if_missing(group, "--plot-grid-width", type=float, help="Grid-line width in points.")
    _add_if_missing(group, "--plot-x-min", type=float, help="Displayed x-axis minimum.")
    _add_if_missing(group, "--plot-x-max", type=float, help="Displayed x-axis maximum.")
    _add_if_missing(group, "--plot-y-min", type=float, help="Displayed y-axis minimum.")
    _add_if_missing(group, "--plot-y-max", type=float, help="Displayed y-axis maximum.")
    _add_if_missing(group, "--plot-line-width", type=float, help="Width of plotted data lines in points.")
    _add_if_missing(group, "--plot-line-color", help="Color override for single-series line plots.")
    _add_if_missing(group, "--plot-fill-color", help="Fill-color override for single-series filled/bar plots.")

    if not _parser_has_option(parser, "--plot-points") and not _parser_has_option(parser, "--no-plot-points"):
        points = group.add_mutually_exclusive_group()
        points.add_argument("--plot-points", dest="plot_points", action="store_true", default=None, help="Show point markers on line plots.")
        points.add_argument("--no-plot-points", dest="plot_points", action="store_false", help="Hide point markers on line plots.")
    _add_if_missing(group, "--plot-point-size", type=float, help="Line-plot point marker size in points.")
    _add_if_missing(group, "--plot-point-fill", help="Fill color of line-plot point markers.")
    _add_if_missing(group, "--plot-point-edge", help="Outline color of line-plot point markers.")
    _add_if_missing(group, "--plot-point-edge-width", type=float, help="Outline width of line-plot point markers.")
    _add_if_missing(group, "--plot-point-shape", choices=tuple(_POINT_MARKERS), help="Point marker shape: circle, square, triangle, or diamond.")
    _add_if_missing(
        group, "--plot-label-points", choices=("none", "peaks", "all"),
        default=label_points_default,
        help=("Label no points, only command-defined peak/local-maximum calls, or all plotted points "
              f"(default: {label_points_default})."),
    )
    _add_if_missing(
        group, "--plot-point-label-value", choices=("x", "y", "both"), default="x",
        help="Value written above labelled points: x coordinate, y coordinate, or both (default: x).",
    )
    _add_if_missing(
        group, "--plot-point-label-offset", type=float,
        help="Vertical label offset above the point in typographic points; default: 4.",
    )

    if not _parser_has_option(parser, "--plot-legend") and not _parser_has_option(parser, "--no-plot-legend"):
        legend = group.add_mutually_exclusive_group()
        legend.add_argument("--plot-legend", dest="plot_legend", action="store_true", default=None, help="Show a legend when labelled plot elements are available.")
        legend.add_argument("--no-plot-legend", dest="plot_legend", action="store_false", help="Hide the plot legend.")
    _add_if_missing(
        group, "--plot-legend-position",
        choices=("best", "upper-right", "upper-left", "lower-right", "lower-left", "outside-right"),
        help="Legend position.",
    )
    _add_if_missing(group, "--plot-x-tick-rotation", type=float, help="Rotate x-axis tick labels by this many degrees.")
    _add_if_missing(group, "--plot-y-tick-rotation", type=float, help="Rotate y-axis tick labels by this many degrees.")
    _add_if_missing(group, "--plot-transparent", action="store_true", help="Save figures with a transparent background.")

    if not show_in_help:
        _hide_plotting_help_actions(parser)
        if not _parser_has_option(parser, "--help-plotting"):
            group.add_argument(
                "--help-plotting",
                action=_ExpandedPlottingHelpAction,
                help="Show this command's full plot customization options and exit.",
            )


def _option_value(args: Any, name: str, fallback: Any = None) -> Any:
    return getattr(args, name, fallback) if args is not None else fallback


def _validate_color(value: str | None, option: str) -> None:
    if value is None:
        return
    from matplotlib.colors import is_color_like
    if not is_color_like(value):
        raise ValueError(f"{option} is not a valid Matplotlib color: {value!r}")


def plot_options_from_namespace(args: Any) -> PlotOptions:
    """Build and validate :class:`PlotOptions` from an argparse namespace."""

    existing = get_plot_options()
    options = PlotOptions(
        format=str(_option_value(args, "plot_format", existing.format) or existing.format),
        width=_option_value(args, "plot_width", None),
        height=_option_value(args, "plot_height", None),
        dpi=_option_value(args, "plot_dpi", None),
        title=_option_value(args, "plot_title", None),
        no_title=bool(_option_value(args, "no_plot_title", False)),
        x_label=_option_value(args, "plot_x_label", None),
        y_label=_option_value(args, "plot_y_label", None),
        font_size=_option_value(args, "plot_font_size", None),
        grid=_option_value(args, "plot_grid", None),
        grid_color=_option_value(args, "plot_grid_color", None),
        grid_alpha=_option_value(args, "plot_grid_alpha", None),
        grid_width=_option_value(args, "plot_grid_width", None),
        x_min=_option_value(args, "plot_x_min", None),
        x_max=_option_value(args, "plot_x_max", None),
        y_min=_option_value(args, "plot_y_min", None),
        y_max=_option_value(args, "plot_y_max", None),
        line_width=_option_value(args, "plot_line_width", None),
        line_color=_option_value(args, "plot_line_color", None),
        fill_color=_option_value(args, "plot_fill_color", None),
        points=_option_value(args, "plot_points", None),
        point_size=_option_value(args, "plot_point_size", None),
        point_fill=_option_value(args, "plot_point_fill", None),
        point_edge=_option_value(args, "plot_point_edge", None),
        point_edge_width=_option_value(args, "plot_point_edge_width", None),
        point_shape=_option_value(args, "plot_point_shape", None),
        label_points=str(_option_value(args, "plot_label_points", existing.label_points) or existing.label_points),
        point_label_value=str(_option_value(args, "plot_point_label_value", existing.point_label_value) or existing.point_label_value),
        point_label_offset=_option_value(args, "plot_point_label_offset", None),
        legend=_option_value(args, "plot_legend", None),
        legend_position=_option_value(args, "plot_legend_position", None),
        x_tick_rotation=_option_value(args, "plot_x_tick_rotation", None),
        y_tick_rotation=_option_value(args, "plot_y_tick_rotation", None),
        transparent=bool(_option_value(args, "plot_transparent", False)),
    )
    validate_plot_options(options)
    return options


def validate_plot_options(options: PlotOptions) -> None:
    if options.format not in {"png", "svg"}:
        raise ValueError("--plot-format must be png or svg.")
    for name, value in (("--plot-width", options.width), ("--plot-height", options.height), ("--plot-font-size", options.font_size), ("--plot-grid-width", options.grid_width), ("--plot-line-width", options.line_width), ("--plot-point-size", options.point_size), ("--plot-point-edge-width", options.point_edge_width)):
        if value is not None and (not math.isfinite(float(value)) or float(value) <= 0):
            raise ValueError(f"{name} must be a finite value greater than zero.")
    if options.dpi is not None and int(options.dpi) < 1:
        raise ValueError("--plot-dpi must be at least 1.")
    if options.grid_alpha is not None and (not math.isfinite(float(options.grid_alpha)) or not 0.0 <= float(options.grid_alpha) <= 1.0):
        raise ValueError("--plot-grid-alpha must be between 0 and 1.")
    for name, low, high in (("x", options.x_min, options.x_max), ("y", options.y_min, options.y_max)):
        if low is not None and not math.isfinite(float(low)):
            raise ValueError(f"--plot-{name}-min must be finite.")
        if high is not None and not math.isfinite(float(high)):
            raise ValueError(f"--plot-{name}-max must be finite.")
        if low is not None and high is not None and float(low) >= float(high):
            raise ValueError(f"--plot-{name}-min must be less than --plot-{name}-max.")
    if options.label_points not in {"none", "peaks", "all"}:
        raise ValueError("--plot-label-points must be none, peaks, or all.")
    if options.point_label_value not in {"x", "y", "both"}:
        raise ValueError("--plot-point-label-value must be x, y, or both.")
    if options.point_label_offset is not None and not math.isfinite(float(options.point_label_offset)):
        raise ValueError("--plot-point-label-offset must be finite.")
    if options.no_title and options.title is not None:
        raise ValueError("--plot-title and --no-plot-title cannot be used together.")
    for option, value in (("--plot-grid-color", options.grid_color), ("--plot-line-color", options.line_color), ("--plot-fill-color", options.fill_color), ("--plot-point-fill", options.point_fill), ("--plot-point-edge", options.point_edge)):
        _validate_color(value, option)


def configure_plot_options(args: Any | None = None, *, options: PlotOptions | None = None) -> PlotOptions:
    """Set process-wide plotting overrides for the current command invocation."""

    global _CURRENT_PLOT_OPTIONS
    if options is None:
        options = plot_options_from_namespace(args)
    validate_plot_options(options)
    _CURRENT_PLOT_OPTIONS = options
    return options


def _plot_options_from_environment() -> PlotOptions | None:
    prefix = "NUCLEOSUITE_PLOT_"
    if not any(key.startswith(prefix) for key in os.environ):
        return None

    def env(name: str) -> str | None:
        value = os.environ.get(prefix + name)
        return value if value not in {None, ""} else None

    def as_float(name: str) -> float | None:
        value = env(name)
        return None if value is None else float(value)

    def as_int(name: str) -> int | None:
        value = env(name)
        return None if value is None else int(value)

    def as_bool(name: str) -> bool | None:
        value = env(name)
        if value is None:
            return None
        return value.lower() in {"1", "true", "yes", "on"}

    options = PlotOptions(
        format=env("FORMAT") or "png",
        width=as_float("WIDTH"), height=as_float("HEIGHT"), dpi=as_int("DPI"),
        title=env("TITLE"), no_title=bool(as_bool("NO_TITLE")),
        x_label=env("X_LABEL"), y_label=env("Y_LABEL"), font_size=as_float("FONT_SIZE"),
        grid=env("GRID"), grid_color=env("GRID_COLOR"), grid_alpha=as_float("GRID_ALPHA"), grid_width=as_float("GRID_WIDTH"),
        x_min=as_float("X_MIN"), x_max=as_float("X_MAX"), y_min=as_float("Y_MIN"), y_max=as_float("Y_MAX"),
        line_width=as_float("LINE_WIDTH"), line_color=env("LINE_COLOR"), fill_color=env("FILL_COLOR"),
        points=as_bool("POINTS"), point_size=as_float("POINT_SIZE"), point_fill=env("POINT_FILL"), point_edge=env("POINT_EDGE"), point_edge_width=as_float("POINT_EDGE_WIDTH"), point_shape=env("POINT_SHAPE"),
        label_points=env("LABEL_POINTS") or "none", point_label_value=env("POINT_LABEL_VALUE") or "x", point_label_offset=as_float("POINT_LABEL_OFFSET"),
        legend=as_bool("LEGEND"), legend_position=env("LEGEND_POSITION"),
        x_tick_rotation=as_float("X_TICK_ROTATION"), y_tick_rotation=as_float("Y_TICK_ROTATION"),
        transparent=bool(as_bool("TRANSPARENT")),
    )
    validate_plot_options(options)
    return options


def get_plot_options() -> PlotOptions:
    environment = _plot_options_from_environment()
    return environment if environment is not None else _CURRENT_PLOT_OPTIONS


def plot_environment(options: PlotOptions | None = None) -> dict[str, str]:
    """Return environment variables that propagate plotting overrides to subprocesses."""

    options = options or get_plot_options()
    values = {
        "FORMAT": options.format, "WIDTH": options.width, "HEIGHT": options.height,
        "DPI": options.dpi, "TITLE": options.title, "NO_TITLE": options.no_title,
        "X_LABEL": options.x_label, "Y_LABEL": options.y_label, "FONT_SIZE": options.font_size,
        "GRID": options.grid, "GRID_COLOR": options.grid_color, "GRID_ALPHA": options.grid_alpha,
        "GRID_WIDTH": options.grid_width, "X_MIN": options.x_min, "X_MAX": options.x_max,
        "Y_MIN": options.y_min, "Y_MAX": options.y_max, "LINE_WIDTH": options.line_width,
        "LINE_COLOR": options.line_color, "FILL_COLOR": options.fill_color, "POINTS": options.points,
        "POINT_SIZE": options.point_size, "POINT_FILL": options.point_fill, "POINT_EDGE": options.point_edge,
        "POINT_EDGE_WIDTH": options.point_edge_width, "POINT_SHAPE": options.point_shape,
        "LABEL_POINTS": options.label_points, "POINT_LABEL_VALUE": options.point_label_value,
        "POINT_LABEL_OFFSET": options.point_label_offset,
        "LEGEND": options.legend, "LEGEND_POSITION": options.legend_position,
        "X_TICK_ROTATION": options.x_tick_rotation, "Y_TICK_ROTATION": options.y_tick_rotation,
        "TRANSPARENT": options.transparent,
    }
    output: dict[str, str] = {}
    for key, value in values.items():
        if value is None:
            continue
        if isinstance(value, bool):
            rendered = "1" if value else "0"
        else:
            rendered = str(value)
        output[f"NUCLEOSUITE_PLOT_{key}"] = rendered
    return output


def plot_path(path: str | Path, options: PlotOptions | None = None) -> Path:
    """Return ``path`` with the configured PNG/SVG extension."""

    options = options or get_plot_options()
    path = Path(path)
    return path.with_suffix(f".{options.format}")


def plot_file_variants(path: str | Path) -> tuple[Path, Path]:
    """Return the PNG and SVG variants of a generated figure path.

    Cleanup/resume code uses this helper so switching ``--plot-format`` never
    leaves a stale figure from an earlier run in the alternate format.
    """

    path = Path(path)
    return path.with_suffix(".png"), path.with_suffix(".svg")


def remove_plot_variants(path: str | Path) -> None:
    """Remove both supported figure-format variants if they exist."""

    for candidate in plot_file_variants(path):
        candidate.unlink(missing_ok=True)


def _primary_axes(figure) -> list[Any]:
    axes = []
    for axis in figure.axes:
        if getattr(axis, "get_label", lambda: "")() == "<colorbar>":
            continue
        axes.append(axis)
    return axes


def _data_lines(axis) -> list[Any]:
    lines = []
    for line in axis.get_lines():
        try:
            x = np.asarray(line.get_xdata(), dtype=float).reshape(-1)
            y = np.asarray(line.get_ydata(), dtype=float).reshape(-1)
        except (TypeError, ValueError):
            continue
        if x.size < 2 or y.size < 2:
            continue
        # Two-point constant lines are normally reference lines from axvline/
        # axhline rather than the scientific data series.
        if x.size <= 2 and (np.allclose(x, x[0]) or np.allclose(y, y[0])):
            continue
        lines.append(line)
    return lines


def _format_point_label_value(value: float) -> str:
    value = float(value)
    if math.isfinite(value) and math.isclose(value, round(value), rel_tol=1e-9, abs_tol=1e-9):
        return str(int(round(value)))
    return f"{value:.6g}"


def point_label_text(x: float, y: float, *, options: PlotOptions | None = None) -> str:
    """Return the configured text for one labelled point."""

    options = options or get_plot_options()
    if options.point_label_value == "y":
        return _format_point_label_value(y)
    if options.point_label_value == "both":
        return f"{_format_point_label_value(x)}, {_format_point_label_value(y)}"
    return _format_point_label_value(x)


def annotate_points(
    axis,
    x_values: Sequence[float],
    y_values: Sequence[float],
    *,
    points_are_peaks: bool = False,
    options: PlotOptions | None = None,
) -> int:
    """Annotate selected scientific points directly above their coordinates.

    ``--plot-label-points peaks`` labels the supplied coordinates only when
    ``points_are_peaks`` is true. ``all`` labels all supplied coordinates.
    Labels are centred above their points and rotated vertically (90 degrees).
    The function returns the number of labels added.
    """

    options = options or get_plot_options()
    mode = options.label_points
    if mode == "none" or (mode == "peaks" and not points_are_peaks):
        return 0
    x = np.asarray(list(x_values), dtype=float).reshape(-1)
    y = np.asarray(list(y_values), dtype=float).reshape(-1)
    if x.size != y.size:
        raise ValueError("Point-label x and y arrays must have the same length.")
    offset = 4.0 if options.point_label_offset is None else float(options.point_label_offset)
    count = 0
    for xv, yv in zip(x, y):
        if not (math.isfinite(float(xv)) and math.isfinite(float(yv))):
            continue
        axis.annotate(
            point_label_text(float(xv), float(yv), options=options),
            xy=(float(xv), float(yv)),
            xytext=(0.0, offset),
            textcoords="offset points",
            ha="center",
            va="bottom",
            rotation=90,
            rotation_mode="default",
            clip_on=False,
            annotation_clip=False,
        )
        count += 1
    return count


def apply_plot_options(figure, *, options: PlotOptions | None = None) -> None:
    """Apply shared overrides after a plot has established its own defaults."""

    options = options or get_plot_options()
    axes = _primary_axes(figure)
    if not axes:
        return

    if options.width is not None or options.height is not None:
        current_width, current_height = figure.get_size_inches()
        figure.set_size_inches(
            float(options.width if options.width is not None else current_width),
            float(options.height if options.height is not None else current_height),
            forward=True,
        )

    if options.no_title:
        if getattr(figure, "_suptitle", None) is not None:
            figure._suptitle.set_text("")
        for axis in axes:
            axis.set_title("")
    elif options.title is not None:
        if len(axes) == 1:
            axes[0].set_title(options.title)
        else:
            figure.suptitle(options.title)

    if options.x_label is not None:
        for axis in axes:
            axis.set_xlabel(options.x_label)
    if options.y_label is not None:
        for axis in axes:
            axis.set_ylabel(options.y_label)

    if options.font_size is not None:
        size = float(options.font_size)
        if getattr(figure, "_suptitle", None) is not None:
            figure._suptitle.set_fontsize(size * 1.15)
        for axis in axes:
            axis.title.set_fontsize(size * 1.15)
            axis.xaxis.label.set_fontsize(size)
            axis.yaxis.label.set_fontsize(size)
            axis.tick_params(axis="both", labelsize=size * 0.9)
            for text in axis.texts:
                text.set_fontsize(size * 0.9)
            legend = axis.get_legend()
            if legend is not None:
                for text in legend.get_texts():
                    text.set_fontsize(size * 0.9)

    if options.grid is not None:
        color = options.grid_color or "0.5"
        alpha = 0.5 if options.grid_alpha is None else float(options.grid_alpha)
        width = 0.8 if options.grid_width is None else float(options.grid_width)
        for axis in axes:
            axis.grid(False, which="both", axis="both")
            if options.grid != "none":
                target = "both" if options.grid == "both" else options.grid
                axis.set_axisbelow(True)
                axis.grid(True, which="major", axis=target, color=color, alpha=alpha, linewidth=width)
    elif any(value is not None for value in (options.grid_color, options.grid_alpha, options.grid_width)):
        # Style existing major grid lines without changing whether the command
        # elected to show them.
        for axis in axes:
            for line in list(axis.get_xgridlines()) + list(axis.get_ygridlines()):
                if options.grid_color is not None:
                    line.set_color(options.grid_color)
                if options.grid_alpha is not None:
                    line.set_alpha(float(options.grid_alpha))
                if options.grid_width is not None:
                    line.set_linewidth(float(options.grid_width))

    for axis in axes:
        if options.x_min is not None or options.x_max is not None:
            left, right = axis.get_xlim()
            axis.set_xlim(
                left if options.x_min is None else float(options.x_min),
                right if options.x_max is None else float(options.x_max),
            )
        if options.y_min is not None or options.y_max is not None:
            bottom, top = axis.get_ylim()
            axis.set_ylim(
                bottom if options.y_min is None else float(options.y_min),
                top if options.y_max is None else float(options.y_max),
            )

        lines = _data_lines(axis)
        if options.line_width is not None:
            for line in lines:
                line.set_linewidth(float(options.line_width))
        if options.line_color is not None and len(lines) == 1:
            lines[0].set_color(options.line_color)

        marker_requested = options.points is True or any(
            value is not None for value in (
                options.point_size, options.point_fill, options.point_edge,
                options.point_edge_width, options.point_shape,
            )
        )
        if options.points is False:
            for line in lines:
                line.set_marker("")
        elif marker_requested:
            marker = _POINT_MARKERS.get(options.point_shape or "circle", "o")
            for line in lines:
                line.set_marker(marker)
                if options.point_size is not None:
                    line.set_markersize(float(options.point_size))
                if options.point_fill is not None:
                    line.set_markerfacecolor(options.point_fill)
                elif options.points is True and line.get_markerfacecolor() in {"none", "None"}:
                    line.set_markerfacecolor(line.get_color())
                if options.point_edge is not None:
                    line.set_markeredgecolor(options.point_edge)
                if options.point_edge_width is not None:
                    line.set_markeredgewidth(float(options.point_edge_width))

        if options.fill_color is not None:
            # Only override an axis whose filled artists already form a single
            # visual series. Multi-state/stacked/category palettes are retained.
            patches = list(axis.patches)
            collections = [collection for collection in axis.collections if collection.__class__.__name__ == "PolyCollection"]
            if patches:
                existing = {tuple(round(float(c), 6) for c in patch.get_facecolor()) for patch in patches}
                if len(existing) <= 1:
                    for patch in patches:
                        patch.set_facecolor(options.fill_color)
            if len(collections) == 1:
                collections[0].set_facecolor(options.fill_color)

        if options.x_tick_rotation is not None:
            for label in axis.get_xticklabels():
                label.set_rotation(float(options.x_tick_rotation))
        if options.y_tick_rotation is not None:
            for label in axis.get_yticklabels():
                label.set_rotation(float(options.y_tick_rotation))

        if options.label_points == "all":
            for line in lines:
                try:
                    lx = np.asarray(line.get_xdata(), dtype=float).reshape(-1)
                    ly = np.asarray(line.get_ydata(), dtype=float).reshape(-1)
                except (TypeError, ValueError):
                    continue
                if lx.size == ly.size and lx.size:
                    annotate_points(axis, lx, ly, points_are_peaks=False, options=options)

        legend = axis.get_legend()
        if options.legend is False:
            if legend is not None:
                legend.remove()
        elif options.legend is True or options.legend_position is not None:
            if legend is not None:
                handles, labels = axis.get_legend_handles_labels()
                legend.remove()
            else:
                handles, labels = axis.get_legend_handles_labels()
            if handles and labels:
                position = options.legend_position or "best"
                if position == "outside-right":
                    axis.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
                else:
                    axis.legend(handles, labels, loc=_LEGEND_LOCATIONS[position], frameon=False)


_CURRENT_PLOT_METADATA: dict[str, object] = {}


def configure_plot_metadata(command: str, argv, parameters=None) -> None:
    """Set command-level provenance written beside every generated plot."""
    import shlex
    from nucleosuite.command_logging import serializable_parameters
    global _CURRENT_PLOT_METADATA
    _CURRENT_PLOT_METADATA = {
        "command": str(command),
        "invocation": shlex.join(["nucleosuite", *list(argv)]),
        "parameters": serializable_parameters(parameters),
    }


def plot_metadata_path(plot_path_value: str | Path) -> Path:
    path = Path(plot_path_value)
    return path.with_name(path.stem + "_metadata.tsv")


def write_plot_metadata(plot_path_value: str | Path, *, extra=None) -> Path:
    """Write a complete parameter sidecar for one plot output."""
    import json
    from nucleosuite import __version__
    path = Path(plot_path_value)
    metadata = plot_metadata_path(path)
    metadata.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        ("nucleosuite_version", __version__),
        ("plot_output", str(path)),
        ("command", _CURRENT_PLOT_METADATA.get("command", "")),
        ("invocation", _CURRENT_PLOT_METADATA.get("invocation", "")),
    ]
    params = _CURRENT_PLOT_METADATA.get("parameters", {})
    if isinstance(params, dict):
        for key, value in sorted(params.items()):
            if isinstance(value, (dict, list, tuple)):
                value = json.dumps(value, sort_keys=True, separators=(",", ":"))
            rows.append((f"parameter.{key}", value))
    for key, value in sorted((extra or {}).items()):
        if isinstance(value, (dict, list, tuple)):
            value = json.dumps(value, sort_keys=True, separators=(",", ":"))
        rows.append((str(key), value))
    with metadata.open("w", encoding="utf-8", newline="") as handle:
        handle.write("field\tvalue\n")
        for key, value in rows:
            text = "" if value is None else str(value).replace("\t", " ").replace("\n", " ")
            handle.write(f"{key}\t{text}\n")
    return metadata


def save_figure(
    figure,
    output_path: str | Path,
    *,
    default_dpi: int = 220,
    bbox_inches: str | None = "tight",
    options: PlotOptions | None = None,
    **savefig_kwargs: Any,
) -> Path:
    """Apply shared options, resolve the extension, save, and return the path."""

    options = options or get_plot_options()
    apply_plot_options(figure, options=options)
    output = plot_path(output_path, options)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Keep one canonical figure per stem. If the user changes format between
    # runs, remove the alternate PNG/SVG rather than leaving a stale plot that
    # could be mistaken for current output by resume/combine workflows.
    for candidate in plot_file_variants(output):
        if candidate != output:
            candidate.unlink(missing_ok=True)
    kwargs: dict[str, Any] = dict(savefig_kwargs)
    if bbox_inches is not None:
        kwargs["bbox_inches"] = bbox_inches
    kwargs["transparent"] = bool(options.transparent)
    # DPI can still matter for rasterized artists embedded inside an SVG.
    kwargs["dpi"] = int(options.dpi if options.dpi is not None else default_dpi)
    figure.savefig(output, **kwargs)
    write_plot_metadata(output)
    return output


def extract_plotting_argv(argv: Sequence[str]) -> tuple[list[str], dict[str, str]]:
    """Remove explicitly supplied shared plot options and return env overrides.

    This is used by suite wrappers whose public argument parsing is implemented
    by packaged shell workflows. Only options the user actually supplied are
    exported, so command-specific defaults (for example DAC peak labels) remain
    active when the suite user does not override them.
    """
    import argparse

    parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
    add_plotting_arguments(parser, show_in_help=True)
    for action in parser._actions:
        if action.dest != "help":
            action.default = argparse.SUPPRESS
    known, remaining = parser.parse_known_args(list(argv))
    values = vars(known)
    mapping = {
        "plot_format": "FORMAT", "plot_width": "WIDTH", "plot_height": "HEIGHT",
        "plot_dpi": "DPI", "plot_title": "TITLE", "no_plot_title": "NO_TITLE",
        "plot_x_label": "X_LABEL", "plot_y_label": "Y_LABEL", "plot_font_size": "FONT_SIZE",
        "plot_grid": "GRID", "plot_grid_color": "GRID_COLOR", "plot_grid_alpha": "GRID_ALPHA",
        "plot_grid_width": "GRID_WIDTH", "plot_x_min": "X_MIN", "plot_x_max": "X_MAX",
        "plot_y_min": "Y_MIN", "plot_y_max": "Y_MAX", "plot_line_width": "LINE_WIDTH",
        "plot_line_color": "LINE_COLOR", "plot_fill_color": "FILL_COLOR", "plot_points": "POINTS",
        "plot_point_size": "POINT_SIZE", "plot_point_fill": "POINT_FILL", "plot_point_edge": "POINT_EDGE",
        "plot_point_edge_width": "POINT_EDGE_WIDTH", "plot_point_shape": "POINT_SHAPE",
        "plot_label_points": "LABEL_POINTS", "plot_point_label_value": "POINT_LABEL_VALUE",
        "plot_point_label_offset": "POINT_LABEL_OFFSET", "plot_legend": "LEGEND",
        "plot_legend_position": "LEGEND_POSITION", "plot_x_tick_rotation": "X_TICK_ROTATION",
        "plot_y_tick_rotation": "Y_TICK_ROTATION", "plot_transparent": "TRANSPARENT",
    }
    env: dict[str, str] = {}
    for dest, value in values.items():
        key = mapping.get(dest)
        if key is None:
            continue
        if isinstance(value, bool):
            rendered = "1" if value else "0"
        else:
            rendered = str(value)
        env[f"NUCLEOSUITE_PLOT_{key}"] = rendered
    return remaining, env
