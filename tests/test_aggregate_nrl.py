"""Unified aggregate peak calling and directional repeat-length tests."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pytest

from nucleosuite.align import (
    AlignmentConfig,
    analyse_aggregate_nrl,
    resolve_output_paths,
    write_aggregate_nrl_outputs,
)
from nucleosuite.cli.aggregate import add_aggregate_parser
from nucleosuite.nrl import moving_average_by_distance
from nucleosuite.replot import _read_table, detect_plot_type


def _aggregate_profile(central_position: int = 0) -> tuple[np.ndarray, np.ndarray]:
    positions = np.arange(-1200, 1201, dtype=float)
    values = np.zeros_like(positions)
    peaks = (
        (-950, 3.0),
        (-760, 4.0),
        (-570, 5.0),
        (-380, 6.0),
        (-190, 7.0),
        (central_position, 10.0),
        (185, 8.0),
        (370, 7.0),
        (555, 6.0),
        (740, 5.0),
        (925, 4.0),
    )
    for centre, amplitude in peaks:
        values += amplitude * np.exp(-0.5 * ((positions - centre) / 22.0) ** 2)
    return positions, values


def test_aggregate_nrl_calls_once_across_zero_then_regresses_directions() -> None:
    positions, values = _aggregate_profile()
    result = analyse_aggregate_nrl(
        values,
        positions=positions,
        peak_resolution=160,
        regression_min=100,
        regression_max=1000,
    )
    assert result.detection_window == 51
    assert result.local_max_window == 21
    assert np.allclose(
        result.local_values,
        moving_average_by_distance(positions, values, 21),
    )
    assert [peak.distance for peak in result.peaks] == [
        -950,
        -760,
        -570,
        -380,
        -190,
        0,
        185,
        370,
        555,
        740,
        925,
    ]
    assert [peak.distance for peak in result.positive_peaks] == [185, 370, 555, 740, 925]
    assert [peak.distance for peak in result.negative_peaks] == [190, 380, 570, 760, 950]
    assert result.positive_regression.slope == pytest.approx(185.0)
    assert result.negative_regression.slope == pytest.approx(190.0)


def test_central_peak_is_shared_as_order_zero_by_both_regressions() -> None:
    positions, values = _aggregate_profile()
    result = analyse_aggregate_nrl(
        values,
        positions=positions,
        regression_min=0,
        regression_max=1000,
    )
    assert result.central_peak is not None
    assert result.central_peak.distance == 0
    assert [peak.distance for peak in result.positive_peaks] == [
        0,
        185,
        370,
        555,
        740,
        925,
    ]
    assert [peak.distance for peak in result.negative_peaks] == [
        0,
        190,
        380,
        570,
        760,
        950,
    ]
    assert result.positive_peak_numbers == (0, 1, 2, 3, 4, 5)
    assert result.negative_peak_numbers == (0, 1, 2, 3, 4, 5)
    assert result.positive_regression.slope == pytest.approx(185.0)
    assert result.negative_regression.slope == pytest.approx(190.0)


def test_near_zero_peak_is_shared_within_half_resolution() -> None:
    positions, values = _aggregate_profile(central_position=7)
    result = analyse_aggregate_nrl(
        values,
        positions=positions,
        regression_max=1000,
    )
    assert result.central_peak is not None
    assert result.central_peak.distance == 7
    assert result.positive_peaks[0].distance == 7
    assert result.negative_peaks[0].distance == 7
    assert result.positive_peak_numbers[0] == 0
    assert result.negative_peak_numbers[0] == 0


def test_signed_exclusion_zone_preserves_peak_orders_and_unified_calls() -> None:
    positions, values = _aggregate_profile()
    result = analyse_aggregate_nrl(
        values,
        positions=positions,
        regression_max=1000,
        exclusion_start=-600,
        exclusion_end=-300,
    )
    assert [peak.distance for peak in result.peaks] == [
        -950,
        -760,
        -570,
        -380,
        -190,
        0,
        185,
        370,
        555,
        740,
        925,
    ]
    assert result.positive_peak_numbers == (0, 1, 2, 3, 4, 5)
    assert result.negative_peak_numbers == (0, 1, 4, 5)
    assert [peak.distance for peak in result.negative_peaks] == [0, 190, 760, 950]
    assert result.negative_regression.slope == pytest.approx(190.0)


def test_regression_range_does_not_limit_unified_peak_calling() -> None:
    positions, values = _aggregate_profile()
    result = analyse_aggregate_nrl(
        values,
        positions=positions,
        regression_min=300,
        regression_max=800,
    )
    assert any(peak.distance == 0 for peak in result.peaks)
    assert any(abs(peak.distance) > 800 for peak in result.peaks)
    assert [peak.distance for peak in result.positive_peaks] == [370, 555, 740]
    assert [peak.distance for peak in result.negative_peaks] == [380, 570, 760]


def test_aggregate_cli_nrl_defaults_and_range_options() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_aggregate_parser(subparsers)
    defaults = parser.parse_args(
        ["aggregate", "--bigwig", "signal.bw", "--region-bed", "regions.bed"]
    )
    assert defaults.nrl is True
    assert defaults.nrl_peak_resolution == 160
    assert defaults.nrl_regression_min == 0
    assert defaults.nrl_regression_max is None
    assert defaults.nrl_regression_exclusion_start is None
    assert defaults.nrl_regression_exclusion_end is None
    assert defaults.nrl_exclusion is True

    custom = parser.parse_args(
        [
            "aggregate",
            "--bigwig",
            "signal.bw",
            "--region-bed",
            "regions.bed",
            "--nrl-regression-min",
            "200",
            "--nrl-regression-max",
            "1200",
            "--nrl-exclusion-start",
            "-100",
            "--nrl-exclusion-end",
            "100",
            "--no-nrl",
        ]
    )
    assert custom.nrl is False
    assert custom.nrl_regression_min == 200
    assert custom.nrl_regression_max == 1200
    assert custom.nrl_regression_exclusion_start == -100
    assert custom.nrl_regression_exclusion_end == 100


def test_resolution_derived_exclusion_can_be_disabled_or_overridden() -> None:
    from nucleosuite.align import resolve_nrl_exclusion

    base = dict(bigwig=Path("signal.bw"), region_bed=Path("regions.bed"))
    assert resolve_nrl_exclusion(
        AlignmentConfig(**base, nrl_peak_resolution=200)
    ) == (-100.0, 100.0)
    assert resolve_nrl_exclusion(
        AlignmentConfig(**base, nrl_exclusion=False)
    ) == (None, None)
    assert resolve_nrl_exclusion(
        AlignmentConfig(
            **base,
            nrl_regression_exclusion_start=-25,
            nrl_regression_exclusion_end=40,
        )
    ) == (-25.0, 40.0)


def test_aggregate_nrl_exclusion_bounds_must_be_paired_and_ordered() -> None:
    positions, values = _aggregate_profile()
    with pytest.raises(ValueError, match="must be supplied together"):
        analyse_aggregate_nrl(
            values,
            positions=positions,
            exclusion_start=-100,
        )
    with pytest.raises(ValueError, match="greater than or equal"):
        analyse_aggregate_nrl(
            values,
            positions=positions,
            exclusion_start=100,
            exclusion_end=-100,
        )


def test_aggregate_nrl_outputs_are_unified_and_replottable(tmp_path: Path) -> None:
    positions, values = _aggregate_profile()
    result = analyse_aggregate_nrl(
        values,
        positions=positions,
        regression_min=100,
        regression_max=1000,
        exclusion_start=-600,
        exclusion_end=-300,
    )
    config = AlignmentConfig(
        bigwig=Path("signal.bw"),
        region_bed=Path("regions.bed"),
        output_dir=tmp_path,
        output_prefix="sample",
        window_half=1200,
        nrl_regression_min=100,
        nrl_regression_max=1000,
        nrl_regression_exclusion_start=-600,
        nrl_regression_exclusion_end=-300,
    )
    outputs = resolve_output_paths(config)
    write_aggregate_nrl_outputs(result, config, outputs)
    assert all(path.is_file() for key, path in outputs.items() if key.startswith("nrl_"))

    profile_headers, _ = _read_table(outputs["nrl_profile"])
    positive_headers, _ = _read_table(outputs["nrl_positive_regression"])
    negative_headers, _ = _read_table(outputs["nrl_negative_regression"])
    assert detect_plot_type(outputs["nrl_profile"], profile_headers) == "aggregate-nrl-profile"
    assert detect_plot_type(outputs["nrl_positive_regression"], positive_headers) == "aggregate-nrl-regression"
    assert detect_plot_type(outputs["nrl_negative_regression"], negative_headers) == "aggregate-nrl-regression"

    with outputs["nrl_summary"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["direction"] for row in rows] == ["positive", "negative"]
    assert float(rows[0]["repeat_length_bp"]) == pytest.approx(185.0)
    assert float(rows[1]["repeat_length_bp"]) == pytest.approx(190.0)
    assert all(row["peak_calling_scope"] == "complete_aggregate_alignment" for row in rows)
    assert all(float(row["regression_exclusion_start_bp"]) == -600 for row in rows)
    assert all(float(row["regression_exclusion_end_bp"]) == -300 for row in rows)

    with outputs["nrl_negative_regression"].open() as handle:
        negative_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["peak_number"]) for row in negative_rows] == [1, 4, 5]
