"""Tests for shared numeric distance-axis settings."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from pathlib import Path

from nucleosuite.plotting import (
    apply_distance_x_axis,
    apply_integer_x_axis,
    default_minor_tick_interval,
)
from nucleosuite.profile_plots import (
    plot_count_profile,
    plot_dinucleotide_profile,
    plot_ww_ss_profile,
)


def test_default_minor_tick_intervals():
    assert default_minor_tick_interval(10) == 5
    assert default_minor_tick_interval(100) == 10
    assert default_minor_tick_interval(50) == 25
    assert default_minor_tick_interval(20) == 10


def test_distance_axis_applies_requested_ticks_and_grids():
    figure, axis = plt.subplots()
    axis.plot([0, 100], [0, 1])
    major, minor = apply_distance_x_axis(
        axis,
        major_interval=20,
        minor_interval=5,
    )
    assert major == 20
    assert minor == 5
    assert axis.xaxis.get_minor_locator() is not None
    major_gridlines = axis.get_xgridlines()
    assert major_gridlines
    assert any(line.get_alpha() == pytest.approx(0.5) for line in major_gridlines)
    plt.close(figure)


def test_count_profile_writes_png(tmp_path: Path):
    table = tmp_path / "counts.tsv"
    table.write_text("fragment_length\tcount\n100\t2\n101\t5\n102\t3\n")
    output = tmp_path / "counts.png"
    plot_count_profile(
        str(table),
        str(output),
        x_column="fragment_length",
        y_column="count",
        xlabel="Fragment length (bp)",
    )
    assert output.exists()
    assert output.stat().st_size > 0


def test_dinucleotide_and_ww_ss_profiles_write_pngs(tmp_path: Path):
    table = tmp_path / "dinuc.tsv"
    columns = ["position", "n_valid"] + [
        f"{a}{b}_pct" for a in "ACGT" for b in "ACGT"
    ] + ["WW_pct", "SS_pct"]
    row = ["0", "10"] + ["6.25"] * 16 + ["25", "25"]
    table.write_text("\t".join(columns) + "\n" + "\t".join(row) + "\n")
    dinuc_png = tmp_path / "dinuc.png"
    ww_png = tmp_path / "ww.png"
    plot_dinucleotide_profile(str(table), str(dinuc_png))
    plot_ww_ss_profile(str(table), str(ww_png))
    assert dinuc_png.exists() and dinuc_png.stat().st_size > 0
    assert ww_png.exists() and ww_png.stat().st_size > 0


def test_single_discrete_value_uses_integer_ticks_and_remains_visible():
    figure, axis = plt.subplots()
    line, = axis.plot([145], [7], marker="o", markersize=2)
    apply_integer_x_axis(axis, [145])
    figure.canvas.draw()
    labels = [label.get_text() for label in axis.get_xticklabels() if label.get_text()]
    assert labels
    assert all("." not in label for label in labels)
    assert axis.get_xlim()[0] < 145 < axis.get_xlim()[1]
    assert line.get_marker() == "o"
    plt.close(figure)


def test_automatic_base_pair_ticks_use_multiples_of_ten():
    from nucleosuite.plotting import apply_base_pair_x_axis

    figure, axis = plt.subplots()
    axis.plot([0, 64], [0, 1])
    interval = apply_base_pair_x_axis(axis, [0, 64])
    figure.canvas.draw()
    assert interval == 10
    visible_ticks = [
        tick for tick in axis.get_xticks()
        if axis.get_xlim()[0] <= tick <= axis.get_xlim()[1]
    ]
    assert visible_ticks
    assert all(tick % 10 == pytest.approx(0) for tick in visible_ticks)
    plt.close(figure)


def test_automatic_distance_ticks_use_nice_multiples_of_ten():
    figure, axis = plt.subplots()
    axis.plot([0, 64], [0, 1])
    major, _minor = apply_distance_x_axis(axis)
    assert major == 10
    figure.canvas.draw()
    visible_ticks = [
        tick for tick in axis.get_xticks()
        if axis.get_xlim()[0] <= tick <= axis.get_xlim()[1]
    ]
    assert visible_ticks
    assert all(tick % 10 == pytest.approx(0) for tick in visible_ticks)
    plt.close(figure)


def test_shared_plot_options_can_save_svg_and_remove_stale_png(tmp_path: Path):
    from nucleosuite.plotting import PlotOptions, save_figure

    png = tmp_path / "figure.png"
    png.write_bytes(b"stale")
    figure, axis = plt.subplots()
    axis.plot([0, 1, 2], [0, 2, 1])
    saved = save_figure(
        figure,
        png,
        options=PlotOptions(format="svg"),
    )
    plt.close(figure)

    assert saved == tmp_path / "figure.svg"
    assert saved.exists() and saved.stat().st_size > 0
    assert not png.exists()


def test_shared_plot_overrides_apply_title_labels_size_grid_and_markers():
    from nucleosuite.plotting import PlotOptions, apply_plot_options

    figure, axis = plt.subplots(figsize=(4, 3))
    line, = axis.plot([0, 1, 2], [1, 3, 2])
    options = PlotOptions(
        width=8,
        height=5,
        title="Custom title",
        x_label="Custom x",
        y_label="Custom y",
        grid="both",
        grid_color="0.25",
        grid_alpha=0.4,
        grid_width=1.2,
        line_width=2.5,
        line_color="red",
        points=True,
        point_size=6,
        point_fill="white",
        point_edge="black",
        point_edge_width=1.5,
        point_shape="diamond",
        x_tick_rotation=30,
    )
    apply_plot_options(figure, options=options)
    figure.canvas.draw()

    assert tuple(figure.get_size_inches()) == pytest.approx((8, 5))
    assert axis.get_title() == "Custom title"
    assert axis.get_xlabel() == "Custom x"
    assert axis.get_ylabel() == "Custom y"
    assert line.get_linewidth() == pytest.approx(2.5)
    assert line.get_color() == "red"
    assert line.get_marker() == "D"
    assert line.get_markersize() == pytest.approx(6)
    assert line.get_markerfacecolor() == "white"
    assert line.get_markeredgecolor() == "black"
    assert line.get_markeredgewidth() == pytest.approx(1.5)
    assert any(item.get_visible() for item in axis.get_xgridlines())
    assert any(item.get_visible() for item in axis.get_ygridlines())
    plt.close(figure)


def test_selective_point_labels_are_vertical_above_points_and_default_to_x():
    from nucleosuite.plotting import PlotOptions, annotate_points

    figure, axis = plt.subplots()
    axis.plot([10, 20, 30], [1, 4, 2])
    options = PlotOptions(label_points="peaks", point_label_value="x")
    count = annotate_points(
        axis,
        [20],
        [4],
        points_are_peaks=True,
        options=options,
    )
    assert count == 1
    annotation = axis.texts[0]
    assert annotation.get_text() == "20"
    assert annotation.get_rotation() == pytest.approx(90)
    assert annotation.get_ha() == "center"
    assert annotation.get_va() == "bottom"
    assert annotation.get_rotation_mode() == "default"
    assert tuple(annotation.xy) == pytest.approx((20, 4))
    assert tuple(annotation.get_position()) == pytest.approx((0, 4))
    plt.close(figure)


def test_point_label_can_show_y_or_both_values():
    from nucleosuite.plotting import PlotOptions, point_label_text

    assert point_label_text(185, 0.1234567, options=PlotOptions(point_label_value="y")) == "0.123457"
    assert point_label_text(185, 2.5, options=PlotOptions(point_label_value="both")) == "185, 2.5"


def test_extract_plotting_argv_exports_only_explicit_overrides():
    from nucleosuite.plotting import extract_plotting_argv

    remaining, environment = extract_plotting_argv([
        "--bam", "sample.bam",
        "--plot-format", "svg",
        "--plot-label-points", "peaks",
        "--plot-point-shape", "triangle",
    ])
    assert remaining == ["--bam", "sample.bam"]
    assert environment == {
        "NUCLEOSUITE_PLOT_FORMAT": "svg",
        "NUCLEOSUITE_PLOT_LABEL_POINTS": "peaks",
        "NUCLEOSUITE_PLOT_POINT_SHAPE": "triangle",
    }


def test_dac_peak_labels_are_opt_in_and_nrl_labels_peaks_by_default():
    from nucleosuite import dac, dcc, nrl

    def default_for(module):
        parser = module.build_parser()
        return next(
            action.default for action in parser._actions
            if action.dest == "plot_label_points"
        )

    assert default_for(dac) == "none"
    assert default_for(dcc) == "none"
    assert default_for(nrl) == "peaks"
