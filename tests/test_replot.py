from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from nucleosuite.replot import build_parser, detect_plot_type, main


def _write(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def test_detects_major_plot_tables(tmp_path: Path) -> None:
    cases = [
        ("dac.tsv", ["Distance", "DAC Value"], "dac"),
        ("dcc.tsv", ["Lag", "DCC Value"], "dcc"),
        ("nrl_profile.tsv", ["distance_or_lag", "unsmoothed_value", "local_max_smoothed_value"], "nrl-profile"),
        ("aggregate_nrl_profile.tsv", ["relative_position", "unsmoothed_value", "local_max_smoothed_value", "detection_smoothed_value"], "aggregate-nrl-profile"),
        ("aggregate_nrl_positive_regression.tsv", ["direction", "peak_number", "signed_position_bp", "distance_from_zero_bp", "fitted_distance_from_zero_bp"], "aggregate-nrl-regression"),
        ("dist.tsv", ["order", "scope", "distance_bp", "count"], "distances"),
        ("aggregate.tsv", ["relative_position", "score"], "aggregate-profile"),
        ("fragment.tsv", ["fragment_length", "count"], "fragment-lengths"),
        ("runs.tsv", ["run_length_bp", "count"], "positive-runs"),
        ("freq.tsv", ["dataset", "score", "count"], "peak-score-frequency"),
        ("states.tsv", ["percentile_threshold", "state", "percentage_of_assigned_peaks"], "peak-states"),
        ("state_overlay.tsv", ["state", "rgb", "distance_bp", "raw_percent", "smoothed_percent"], "distance-state-overlay"),
        ("percentile_curves.tsv", ["percentile_threshold", "distance_bp", "raw_count", "raw_percent", "order", "scope"], "distance-percentile-curves"),
        ("percentile_counts.tsv", ["retained_peak_count", "percentile_threshold", "order", "scope"], "distance-percentile-peak-counts"),
        ("pairs.tsv", ["a_score", "b_score", "absolute_distance"], "compare-positions-score"),
        ("distance_bins.tsv", ["bin_start_inclusive", "bin_end_exclusive", "pair_count"], "compare-positions-histogram"),
        ("correlations.tsv", ["distance_bin", "pair_count", "spearman_score_correlation"], "compare-positions-correlation"),
        ("percentiles.tsv", ["percentile_group", "absolute_distance"], "compare-positions-percentile-boxplot"),
        ("spacing.tsv", ["sample", "region", "subset", "profile", "correlation"], "gene-expression-spacing"),
        ("spacing_scatter.tsv", ["sample", "profile", "body_median_spacing_bp", "transformed_expression"], "gene-expression-spacing-scatter"),
        ("fft.tsv", ["sample", "profile", "period_bp", "correlation"], "gene-expression-fft-trajectory"),
        ("ranking.tsv", ["sample", "rank", "profile", "correlation", "ranking_periods"], "gene-expression-ranking"),
        ("tss.tsv", ["quintile", "relative_position", "mean_signal"], "tss-expression"),
        ("dinuc.tsv", ["position", "AA_pct", "WW_pct", "SS_pct"], "dinucleotide-profile"),
        ("ww_lengths.tsv", ["fragment_length", "type1_percent_of_classified", "type2_percent_of_classified"], "ww-type-by-length"),
        ("genes.tsv", ["gene_id", "candidate_sets"], "gene-sets-venn"),
        ("relocation.tsv", ["relocation_bp", "count"], "count-profile"),
        ("overlay.tsv", ["relative_position", "sample_a", "sample_b"], "profile-overlay"),
    ]
    for filename, headers, expected in cases:
        assert detect_plot_type(tmp_path / filename, headers) == expected


def test_plot_parser_exposes_independent_major_minor_ticks_and_grids() -> None:
    parser = build_parser()
    options = {flag for action in parser._actions for flag in action.option_strings}
    for option in (
        "--x-major-tick", "--x-minor-tick", "--y-major-tick", "--y-minor-tick",
        "--x-major-grid", "--x-minor-grid", "--y-major-grid", "--y-minor-grid",
        "--major-grid-style", "--minor-grid-style", "--vmin", "--vmax",
        "--mpl-rc", "--mpl-kw", "--nrl-inset",
        "--normalization",
    ):
        assert option in options


def test_dac_replot_writes_figure_with_auto_detection(tmp_path: Path) -> None:
    path = tmp_path / "sample_dac.tsv"
    rows = []
    for distance in range(1, 801):
        value = 2.0 + 0.3 * np.cos(2 * np.pi * distance / 190.0)
        rows.append([distance, value, 0, value, 1, 0])
    _write(
        path,
        ["Distance", "DAC Value", "DAC Value Percent", "Raw DAC Value", "Opportunities", "DAC per million signal-pairs"],
        rows,
    )
    output = tmp_path / "dac.png"
    assert main([str(path), "--output", str(output), "--nrl-inset", "on"]) == 0
    assert output.is_file() and output.stat().st_size > 1000

def test_dac_replot_default_is_raw_only_and_detection_uses_nrl_layers(tmp_path: Path) -> None:
    from nucleosuite.replot import _plot_dac, _read_table

    path = tmp_path / "sample_dac.tsv"
    rows = []
    for distance in range(1, 1301):
        value = 2.0 + 0.25 * np.cos(2 * np.pi * distance / 190.0)
        rows.append([distance, value, 0, value, 1, 0])
    _write(
        path,
        ["Distance", "DAC Value", "DAC Value Percent", "Raw DAC Value", "Opportunities", "DAC per million signal-pairs"],
        rows,
    )
    headers, table_rows = _read_table(path)

    raw_args = build_parser().parse_args([str(path), "--output", str(tmp_path / "raw.png")])
    (_, raw_fig) = _plot_dac(path, headers, table_rows, raw_args, tmp_path / "raw.png", {})
    raw_ax = raw_fig.axes[0]
    assert len(raw_ax.lines) == 1
    assert not raw_ax.collections
    assert not raw_ax.texts

    detect_args = build_parser().parse_args([
        str(path), "--output", str(tmp_path / "detected.png"),
        "--detect-peaks", "--peak-resolution", "160", "--nrl-inset", "off",
    ])
    (_, detected_fig) = _plot_dac(path, headers, table_rows, detect_args, tmp_path / "detected.png", {})
    detected_ax = detected_fig.axes[0]
    labels = [line.get_label() for line in detected_ax.lines]
    assert "Unsmoothed" in labels
    assert "Local maxima (21 bp)" in labels
    assert "Peak detection (51 bp)" in labels
    assert detected_ax.collections
    assert detected_ax.texts
    for annotation in detected_ax.texts:
        if hasattr(annotation, "xy"):
            assert annotation.get_ha() == "center"
            assert annotation.get_va() == "bottom"


def test_dcc_replot_accepts_explicit_tick_grid_controls(tmp_path: Path) -> None:
    path = tmp_path / "sample_dcc.tsv"
    rows = [[lag, 100 + lag * lag, 0, 100 + lag * lag, 1, 0] for lag in range(-20, 21)]
    _write(
        path,
        ["Lag", "DCC Value", "DCC Value Percent", "Raw DCC Value", "Opportunities", "DCC per million signal-pairs"],
        rows,
    )
    output = tmp_path / "dcc.svg"
    assert main([
        str(path), "--output", str(output), "--format", "svg",
        "--x-major-tick", "10", "--x-minor-tick", "5",
        "--x-major-grid", "--x-minor-grid",
    ]) == 0
    assert output.is_file() and "<svg" in output.read_text(encoding="utf-8")[:500]


def test_heatmap_replot_respects_numeric_matrix_and_saturation_options(tmp_path: Path) -> None:
    path = tmp_path / "sample_heatmap_matrix.tsv"
    _write(path, ["row_index", "-10", "0", "10"], [[1, -2, 0, 2], [2, 2, 0, -2]])
    output = tmp_path / "heatmap.png"
    assert main([
        str(path), "--output", str(output), "--vmin", "-1", "--vmax", "1",
        "--x-major-tick", "10", "--x-minor-tick", "5",
    ]) == 0
    assert output.is_file() and output.stat().st_size > 1000


def test_nrl_peak_table_can_recreate_regression(tmp_path: Path) -> None:
    path = tmp_path / "sample_nrl_peaks.tsv"
    _write(
        path,
        ["peak_number", "distance_or_lag", "unsmoothed_value", "local_max_smoothed_value", "detection_peak_distance", "detection_smoothed_value"],
        [[1, 185, 1, 1, 185, 1], [2, 370, 1, 1, 370, 1], [3, 555, 1, 1, 555, 1]],
    )
    output = tmp_path / "nrl.png"
    assert main([str(path), "--output", str(output)]) == 0
    assert output.is_file() and output.stat().st_size > 1000


def test_distance_replot_uses_solid_coloured_order_lines_and_legend(tmp_path: Path) -> None:
    from nucleosuite.replot import _plot_distances, _read_table

    path = tmp_path / "sample_distances.tsv"
    _write(
        path,
        ["order", "scope", "state", "distance_bp", "count", "smoothed_count"],
        [
            [1, "combined_chromosomes", "All", 180, 10, 9.0],
            [1, "combined_chromosomes", "All", 185, 15, 14.0],
            [2, "combined_chromosomes", "All", 180, 8, 7.0],
            [2, "combined_chromosomes", "All", 185, 12, 11.0],
        ],
    )
    headers, rows = _read_table(path)
    args = build_parser().parse_args([str(path), "--output", str(tmp_path / "distances.png")])
    (_, figure) = _plot_distances(path, headers, rows, args, tmp_path / "distances.png", {})
    lines = figure.axes[0].lines
    coloured = [line for line in lines if line.get_label() in {"+1", "+2"}]
    assert len(coloured) == 2
    assert all(line.get_linestyle() == "-" for line in coloured)
    assert all(line.get_marker() == "None" for line in coloured)
    assert coloured[0].get_color() != coloured[1].get_color()
    legend = figure.axes[0].get_legend()
    assert legend is not None
    assert legend.get_title().get_text() == "Neighbour order"
    assert tuple(figure.get_size_inches()) == (10.0, 5.5)
    assert not figure.axes[0].texts


def test_distance_replot_artist_options_override_default_style(tmp_path: Path) -> None:
    from nucleosuite.replot import _plot_distances, _read_table

    path = tmp_path / "sample_distances.tsv"
    _write(path, ["distance_bp", "count"], [[180, 10], [185, 15], [190, 8]])
    headers, rows = _read_table(path)
    args = build_parser().parse_args([str(path), "--output", str(tmp_path / "distances.png")])
    (_, figure) = _plot_distances(
        path,
        headers,
        rows,
        args,
        tmp_path / "distances.png",
        {"line": {"linestyle": "--", "marker": "s"}},
    )
    line = figure.axes[0].lines[0]
    assert line.get_linestyle() == "--"
    assert line.get_marker() == "s"


def test_regression_replot_is_square_by_default_and_size_is_overridable(tmp_path: Path) -> None:
    from nucleosuite.replot import _plot_nrl_regression, _read_table

    path = tmp_path / "sample_nrl_regression.tsv"
    _write(
        path,
        ["order", "peak_distance_bp", "fitted_distance_bp"],
        [[1, 185, 185], [2, 370, 370], [3, 555, 555]],
    )
    headers, rows = _read_table(path)

    args = build_parser().parse_args([str(path), "--output", str(tmp_path / "regression.png")])
    (_, figure) = _plot_nrl_regression(
        path, headers, rows, args, tmp_path / "regression.png", {}
    )
    assert tuple(figure.get_size_inches()) == (6.5, 6.5)
    assert figure.axes[0].lines[0].get_linestyle() == ":"
    assert figure.axes[0].collections[0].get_facecolors().size == 0

    custom_args = build_parser().parse_args([
        str(path), "--output", str(tmp_path / "custom.png"),
        "--width", "8", "--height", "4",
    ])
    (_, custom_figure) = _plot_nrl_regression(
        path, headers, rows, custom_args, tmp_path / "custom.png", {}
    )
    assert tuple(custom_figure.get_size_inches()) == (8.0, 4.0)


def test_aggregate_nrl_replots_central_order_and_exclusion_zone(tmp_path: Path) -> None:
    from nucleosuite.replot import (
        _plot_aggregate_nrl_profile,
        _plot_aggregate_nrl_regression,
        _read_table,
    )

    profile_path = tmp_path / "sample_aggregate_nrl_profile.tsv"
    _write(
        profile_path,
        [
            "relative_position",
            "unsmoothed_value",
            "local_max_smoothed_value",
            "detection_smoothed_value",
            "regression_exclusion_start_bp",
            "regression_exclusion_end_bp",
            "is_peak",
        ],
        [
            [-100, 0, 0, 0, -50, 50, 0],
            [0, 1, 1, 1, -50, 50, 1],
            [100, 0, 0, 0, -50, 50, 0],
        ],
    )
    headers, rows = _read_table(profile_path)
    args = build_parser().parse_args(
        [str(profile_path), "--output", str(tmp_path / "profile.png")]
    )
    (_, profile_figure) = _plot_aggregate_nrl_profile(
        profile_path, headers, rows, args, tmp_path / "profile.png", {}
    )
    assert len(profile_figure.axes[0].patches) == 1
    span = profile_figure.axes[0].patches[0]
    assert span.get_x() == -50
    assert span.get_x() + span.get_width() == 50

    regression_path = tmp_path / "sample_aggregate_nrl_positive_regression.tsv"
    _write(
        regression_path,
        [
            "direction",
            "peak_number",
            "signed_position_bp",
            "distance_from_zero_bp",
            "fitted_distance_from_zero_bp",
        ],
        [
            ["positive", 0, 0, 0, 0],
            ["positive", 1, 185, 185, 185],
            ["positive", 4, 740, 740, 740],
        ],
    )
    headers, rows = _read_table(regression_path)
    args = build_parser().parse_args(
        [str(regression_path), "--output", str(tmp_path / "regression.png")]
    )
    (_, regression_figure) = _plot_aggregate_nrl_regression(
        regression_path, headers, rows, args, tmp_path / "regression.png", {}
    )
    plotted_orders = regression_figure.axes[0].collections[0].get_offsets()[:, 0]
    assert list(plotted_orders) == [0, 1, 4]
    assert tuple(regression_figure.get_size_inches()) == (6.5, 6.5)


def test_fragment_length_replot_recovers_source_command_normalisation(tmp_path: Path) -> None:
    from nucleosuite.replot import _plot_fragment_lengths, _read_table

    density_path = tmp_path / "sample_fragment_lengths.tsv"
    _write(density_path, ["fragment_length", "count"], [[100, 2], [101, 6]])
    headers, rows = _read_table(density_path)
    density_args = build_parser().parse_args([str(density_path), "--output", str(tmp_path / "density.png")])
    (_, density_figure) = _plot_fragment_lengths(
        density_path, headers, rows, density_args, tmp_path / "density.png", {}
    )
    density_axis = density_figure.axes[0]
    assert density_axis.get_ylabel() == "Density"
    assert np.sum(density_axis.lines[0].get_ydata()) == 1.0
    assert density_axis.lines[0].get_xdata()[0] == 0
    assert density_axis.lines[0].get_xdata()[-1] == 101

    count_path = tmp_path / "sample.fragment_length_counts.tsv"
    _write(count_path, ["fragment_length", "count"], [[100, 2], [101, 6]])
    headers, rows = _read_table(count_path)
    count_args = build_parser().parse_args([str(count_path), "--output", str(tmp_path / "count.png")])
    (_, count_figure) = _plot_fragment_lengths(
        count_path, headers, rows, count_args, tmp_path / "count.png", {}
    )
    count_axis = count_figure.axes[0]
    assert count_axis.get_ylabel() == "Fragment count"
    assert count_axis.get_title() == "Fragment-length distribution"
    assert tuple(count_figure.get_size_inches()) == (10.0, 5.0)


def test_fragment_heatmap_replot_uses_layout_sidecars(tmp_path: Path) -> None:
    from nucleosuite.replot import _plot_heatmap, _read_table

    matrix_path = tmp_path / "comparison_normalised_matrix.tsv"
    _write(matrix_path, ["profile", "100", "150", "200"], [["B", -1, 0, 1], ["A", 1, 0, -1]])
    _write(
        tmp_path / "comparison_clustered_profiles.tsv",
        ["cluster_order", "original_index", "profile", "sample", "label", "category", "source_file"],
        [[1, 1, "B", "", "", "Cancer", "b.tsv"], [2, 0, "A", "", "", "Healthy", "a.tsv"]],
    )
    _write(
        tmp_path / "comparison_heatmap_plot_metadata.tsv",
        ["setting", "value"],
        [
            ["colourbar_label", "Fragment-length z-score across profiles"],
            ["heatmap_centre", "0"], ["colour_percentile", "99"],
            ["low_colour", "#2166AC"], ["mid_colour", "#FFFFFF"], ["high_colour", "#F58518"],
            ["label_gutter", "1.2"], ["max_yticks", "80"], ["title", "Comparison"],
            ["dendrogram_colour", "black"], ["category_colour:Cancer", "#F58518"],
            ["category_colour:Healthy", "#54A24B"],
        ],
    )
    _write(
        tmp_path / "comparison_heatmap_linkage.tsv",
        ["left_child", "right_child", "distance", "member_count"],
        [[0, 1, 2.0, 2]],
    )
    headers, rows = _read_table(matrix_path)
    args = build_parser().parse_args([str(matrix_path), "--output", str(tmp_path / "heatmap.png")])
    (_, figure) = _plot_heatmap(matrix_path, headers, rows, args, tmp_path / "heatmap.png", {})
    assert len(figure.axes) == 5
    heatmap_axis = next(axis for axis in figure.axes if axis.get_xlabel() == "Fragment length (bp)")
    assert heatmap_axis.get_title() == "Comparison"
    assert heatmap_axis.get_legend() is not None


def test_compare_score_replot_retains_original_scatter_selection_and_method(tmp_path: Path) -> None:
    from nucleosuite.replot import _plot_compare_positions_score, _read_table

    path = tmp_path / "comparison_pairs.tsv"
    _write(
        path,
        [
            "a_score", "b_score", "a_score_z", "b_score_z", "absolute_distance",
            "plot_selected", "plot_score_normalization", "plot_correlation_method",
            "plot_label_a", "plot_label_b", "plot_score_z_limit",
        ],
        [
            [1, 2, -1, -1, 5, 1, "zscore", "spearman", "Method A", "Method B", 8],
            [2, 4, 0, 0, 10, 0, "zscore", "spearman", "Method A", "Method B", 8],
            [3, 6, 1, 1, 15, 1, "zscore", "spearman", "Method A", "Method B", 8],
        ],
    )
    headers, rows = _read_table(path)
    args = build_parser().parse_args([str(path), "--output", str(tmp_path / "score.png")])
    (_, figure) = _plot_compare_positions_score(path, headers, rows, args, tmp_path / "score.png", {})
    axis = figure.axes[0]
    assert axis.collections[0].get_offsets().shape[0] == 2
    assert axis.get_xlim() == (-8.0, 8.0)
    assert axis.get_xlabel() == "Method A score z-score"
    assert "Spearman" in axis.texts[0].get_text()
    assert "Pearson" not in axis.texts[0].get_text()


def test_dinucleotide_replot_writes_figure(tmp_path: Path) -> None:
    path = tmp_path / "sample_dinuc_profile.tsv"
    _write(path, ["position", "AA_pct", "AT_pct", "WW_pct", "SS_pct"], [[-1, 5, 6, 11, 8], [0, 7, 8, 15, 9], [1, 5, 6, 11, 8]])
    output = tmp_path / "dinuc.png"
    assert main([str(path), "--output", str(output)]) == 0
    assert output.is_file() and output.stat().st_size > 1000


def test_profile_overlay_replot_writes_all_numeric_series(tmp_path: Path) -> None:
    path = tmp_path / "overlay.tsv"
    _write(path, ["relative_position", "sample_a", "sample_b"], [[-10, 1, 2], [0, 3, 4], [10, 1, 2]])
    output = tmp_path / "overlay.png"
    assert main([str(path), "--output", str(output)]) == 0
    assert output.is_file() and output.stat().st_size > 1000


def test_plot_parser_makes_dac_peak_detection_opt_in_and_nrl_labels_default():
    args = build_parser().parse_args(["sample.tsv"])
    assert args.detect_peaks is False
    assert args.peak_resolution == 160.0
    assert args.label_peaks == "auto"
