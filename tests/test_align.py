from pathlib import Path

import numpy as np

from nucleosuite.align import (
    AlignmentConfig,
    ChromosomeIntervals,
    find_relative_nucleosome_center,
    has_consecutive_zeros,
    make_output_prefix,
    overlaps_any,
    resolve_output_paths,
    sort_matrix,
)


def test_relative_nucleosome_positive_strand():
    centers = [100, 200, 300, 400]
    assert find_relative_nucleosome_center(centers, 250, "+", 1) == 300
    assert find_relative_nucleosome_center(centers, 250, "+", -1) == 200


def test_relative_nucleosome_negative_strand():
    centers = [100, 200, 300, 400]
    assert find_relative_nucleosome_center(centers, 250, "-", 1) == 200
    assert find_relative_nucleosome_center(centers, 250, "-", -1) == 300


def test_relative_nucleosome_is_strict():
    centers = [100, 200, 300]
    assert find_relative_nucleosome_center(centers, 200, "+", 1) == 300
    assert find_relative_nucleosome_center(centers, 200, "+", -1) == 100


def test_no_clamping():
    assert find_relative_nucleosome_center([100], 100, "+", 1) is None


def test_half_open_state_overlap():
    data = ChromosomeIntervals(starts=[100], intervals=[(100, 200)])
    assert overlaps_any(data, 150, 250)
    assert not overlaps_any(data, 200, 300)
    assert not overlaps_any(data, 0, 100)


def test_consecutive_zeros():
    assert has_consecutive_zeros(np.array([1, 0, 0, 0, 2.0]), 3)
    assert not has_consecutive_zeros(np.array([0, 0, 1, 0, 0.0]), 3)
    assert not has_consecutive_zeros(np.array([0, 0, 0.0]), 0)


def test_default_prefix_uses_all_inputs():
    config = AlignmentConfig(
        bigwig=Path("signal.bw"),
        region_bed=Path("regions.bed"),
        nucleosome_bed=Path("nucleosomes.bed.gz"),
        nucleosome_offset=-2,
        state_bed=Path("active.bed"),
    )
    assert make_output_prefix(config) == "regions_signal_nucleosomes_minus2_active"


def test_alignment_outputs_include_exact_heatmap_matrix(tmp_path: Path):
    config = AlignmentConfig(
        bigwig=Path("signal.bw"),
        region_bed=Path("regions.bed"),
        output_dir=tmp_path,
        output_prefix="aligned",
    )
    outputs = resolve_output_paths(config)
    assert outputs["heatmap"] == tmp_path / "aligned_heatmap.png"
    assert outputs["heatmap_matrix"] == tmp_path / "aligned_heatmap_matrix.tsv.gz"
    assert outputs["plotted_mean"] == tmp_path / "aligned_heatmap_mean.tsv"


def test_mean_absolute_sort_places_highest_signal_first():
    matrix = np.array([[1.0, -1.0, 1.0], [4.0, -4.0, 4.0], [2.0, 2.0, 2.0]])
    sorted_matrix, order, scores = sort_matrix(matrix, "mean_absolute")
    assert order.tolist() == [1, 2, 0]
    assert scores.tolist() == [4.0, 2.0, 1.0]
    assert np.array_equal(sorted_matrix[0], matrix[1])


def test_alignment_defaults_use_2500_mean_absolute_and_nan_to_zero():
    config = AlignmentConfig(bigwig=Path("signal.bw"), region_bed=Path("regions.bed"))
    assert config.window_half == 2500
    assert config.sort_mode == "mean_absolute"
    assert config.nan_to_zero is True


def test_aggregate_cli_can_disable_default_nan_to_zero():
    import argparse
    from nucleosuite.cli.aggregate import add_aggregate_parser

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_aggregate_parser(subparsers)
    default_args = parser.parse_args([
        "aggregate", "--bigwig", "signal.bw", "--region-bed", "regions.bed"
    ])
    disabled_args = parser.parse_args([
        "aggregate", "--bigwig", "signal.bw", "--region-bed", "regions.bed",
        "--no-nan-to-zero",
    ])
    assert default_args.nan_to_zero is True
    assert disabled_args.nan_to_zero is False


def test_sparse_bigwig_nan_values_are_zero_filled_by_default(tmp_path, monkeypatch):
    import sys
    from types import SimpleNamespace
    import pytest
    from nucleosuite.align import run_alignment

    class FakeBigWig:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return False

        def chroms(self):
            return {"chr1": 100}

        def values(self, chromosome, start, end, numpy=False):
            assert chromosome == "chr1"
            result = np.full(end - start, np.nan, dtype=float)
            signal_index = 49 - start
            if 0 <= signal_index < result.size:
                result[signal_index] = 2.0
            return result

    monkeypatch.setitem(sys.modules, "pyBigWig", SimpleNamespace(open=lambda path: FakeBigWig()))
    bigwig_path = tmp_path / "signal.bw"
    bigwig_path.touch()
    region_bed = tmp_path / "regions.bed"
    region_bed.write_text("chr1\t49\t50\tregion\t0\t+\n")

    outputs = run_alignment(AlignmentConfig(
        bigwig=bigwig_path,
        region_bed=region_bed,
        output_dir=tmp_path / "default",
        output_prefix="default",
        window_half=2,
        zero_thresh=0,
        max_score=None,
    ))
    assert outputs["summary"].exists()

    with pytest.raises(RuntimeError, match="No valid records remained"):
        run_alignment(AlignmentConfig(
            bigwig=bigwig_path,
            region_bed=region_bed,
            output_dir=tmp_path / "strict",
            output_prefix="strict",
            window_half=2,
            zero_thresh=0,
            max_score=None,
            nan_to_zero=False,
        ))


def test_minus_strand_regions_are_reversed_into_feature_orientation(tmp_path, monkeypatch):
    import gzip
    import sys
    from types import SimpleNamespace
    from nucleosuite.align import run_alignment

    signal = np.zeros(100, dtype=float)
    signal[18:23] = np.array([1, 2, 3, 4, 5], dtype=float)
    signal[38:43] = np.array([5, 4, 3, 2, 1], dtype=float)

    class FakeBigWig:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return False

        def chroms(self):
            return {"chr1": signal.size}

        def values(self, chromosome, start, end, numpy=False):
            assert chromosome == "chr1"
            return signal[start:end].copy()

    monkeypatch.setitem(sys.modules, "pyBigWig", SimpleNamespace(open=lambda path: FakeBigWig()))
    bigwig_path = tmp_path / "signal.bw"
    bigwig_path.touch()
    regions = tmp_path / "ctcf.bed"
    regions.write_text(
        "chr1\t20\t21\tplus\t0\t+\n"
        "chr1\t40\t41\tminus\t0\t-\n"
    )

    outputs = run_alignment(AlignmentConfig(
        bigwig=bigwig_path,
        region_bed=regions,
        output_dir=tmp_path / "out",
        output_prefix="ctcf",
        window_half=2,
        zero_thresh=0,
        max_score=None,
        sort_mode="unsorted",
        missing_strand="error",
    ))

    with gzip.open(outputs["heatmap_matrix"], "rt") as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle]
    rows = np.asarray([[float(value) for value in line[1:]] for line in lines[1:]])
    assert rows.shape == (2, 5)
    assert np.array_equal(rows[0], np.array([1, 2, 3, 4, 5], dtype=float))
    assert np.array_equal(rows[1], np.array([1, 2, 3, 4, 5], dtype=float))
