from pathlib import Path

import numpy as np
import pytest

from nucleosuite.mean_scale import (
    _bigwig_nonzero_mean,
    _default_output,
    _iter_scaled_intervals,
    _region_score_mean,
    _write_scaled_intervals,
    build_parser,
    main,
)


def test_region_score_mean_uses_column_five_by_default_semantics(tmp_path: Path):
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr1\t0\t10\tp1\t10\t.\n"
        "chr1\t20\t30\tp2\t20\t.\n"
        "chr1\t40\t50\tp3\tnan\t.\n",
        encoding="utf-8",
    )
    mean, count, nonfinite = _region_score_mean(bed, 5)
    assert mean == 15.0
    assert count == 2
    assert nonfinite == 1


class _ValuesBigWig:
    def chroms(self):
        return {"chr1": 5, "chr2": 4}

    def values(self, chrom, start, end, numpy=False):
        values = {
            "chr1": np.array([0.0, 2.0, np.nan, 4.0, 0.0]),
            "chr2": np.array([-2.0, 0.0, 6.0, np.inf]),
        }[chrom][start:end]
        return values if numpy else values.tolist()


def test_bigwig_mean_uses_only_finite_nonzero_bases():
    mean, count = _bigwig_nonzero_mean(_ValuesBigWig())
    assert count == 4
    assert mean == pytest.approx((2 + 4 - 2 + 6) / 4)


class _IntervalsBigWig:
    def intervals(self, chrom, start, end):
        assert chrom == "chr1"
        # Simulates one interval crossing the artificial query chunk boundary.
        rows = [(0, 7, 2.0), (7, 10, -1.0)]
        return [row for row in rows if row[0] < end and row[1] > start]


def test_scaled_interval_iterator_clips_chunks_without_duplication(monkeypatch):
    import nucleosuite.mean_scale as module

    monkeypatch.setattr(module, "_CHUNK_BP", 5)
    rows = list(_iter_scaled_intervals(_IntervalsBigWig(), "chr1", 10, 3.0))
    assert rows == [
        (0, 5, 6.0),
        (5, 7, 6.0),
        (7, 10, -3.0),
    ]


def test_default_output_encodes_reference_mode_and_scale(tmp_path: Path):
    source = tmp_path / "PNS.bw"
    region = tmp_path / "nucleosome_peaks.bed"
    region_path = _default_output(
        source,
        mode="region-score-mean",
        reference_mean=16.7644,
        scale=100,
        regions=region,
        score_column=5,
    )
    assert region_path.name == "PNS_meanscale_regions-nucleosome-peaks-col5_x100.bw"
    explicit_path = _default_output(
        source,
        mode="supplied-reference-mean",
        reference_mean=16.7644,
        scale=100,
        regions=None,
        score_column=5,
    )
    assert explicit_path.name == "PNS_meanscale_mean-16p7644_x100.bw"


def test_mean_scale_parser_supports_reference_mean_alias_and_region_mode():
    parser = build_parser()
    args = parser.parse_args(["signal.bw", "--normalization-mean", "16.7644"])
    assert args.reference_mean == pytest.approx(16.7644)
    with pytest.raises(SystemExit):
        parser.parse_args([
            "signal.bw", "--reference-mean", "16.7644", "--regions", "peaks.bed"
        ])


def test_mean_scale_parser_accepts_a_reference_bigwig_exclusively():
    parser = build_parser()
    args = parser.parse_args(["tns.bw", "--reference-bigwig", "posTNS.bw", "--scale", "1"])
    assert args.reference_bigwig == "posTNS.bw"
    assert args.scale == pytest.approx(1.0)
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["tns.bw", "--reference-bigwig", "posTNS.bw", "--regions", "peaks.bed"]
        )


def test_mean_scale_bed_defaults_to_own_score_mean_and_same_format(tmp_path: Path):
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr1\t0\t10\tp1\t1\t+\n"
        "chr1\t20\t30\tp2\t3\t-\n",
        encoding="utf-8",
    )
    assert main([str(bed), "--scale", "100"]) == 0
    output = tmp_path / "peaks_meanscale_scores-col5_x100.bed"
    assert output.is_file()
    rows = [line.split("\t") for line in output.read_text().splitlines()]
    assert [float(row[4]) for row in rows] == pytest.approx([50.0, 150.0])
    assert rows[0][:4] == ["chr1", "0", "10", "p1"]
    assert rows[1][5] == "-"


def test_mean_scale_bed_can_integer_round_and_clamp(tmp_path: Path):
    bed = tmp_path / "scores.bed"
    bed.write_text(
        "chr1\t0\t10\ta\t-1\t.\n"
        "chr1\t10\t20\tb\t1\t.\n"
        "chr1\t20\t30\tc\t4\t.\n",
        encoding="utf-8",
    )
    output = tmp_path / "scaled.bed"
    assert main([
        str(bed), "--reference-mean", "1", "--scale", "1000",
        "--integer-scores", "--clamp-min", "0", "--clamp-max", "1000",
        "--output", str(output),
    ]) == 0
    values = [int(line.split("\t")[4]) for line in output.read_text().splitlines()]
    assert values == [0, 1000, 1000]


def test_mean_scale_bigbed_controls_are_automatic():
    from nucleosuite.mean_scale import _effective_interval_controls
    integer, lower, upper = _effective_interval_controls(
        output_format="bigbed", integer_scores=False, clamp_min=None, clamp_max=None
    )
    assert integer is True
    assert lower == 0
    assert upper == 1000
    integer, lower, upper = _effective_interval_controls(
        output_format="bigbed", integer_scores=False, clamp_min=100, clamp_max=900
    )
    assert integer is True
    assert lower == 100
    assert upper == 900


def test_mean_scale_parser_exposes_interval_output_controls():
    args = build_parser().parse_args([
        "peaks.bed", "--integer-scores", "--clamp-min", "0", "--clamp-max", "1000",
        "--output-format", "bed.gz",
    ])
    assert args.input == "peaks.bed"
    assert args.integer_scores is True
    assert args.clamp_min == 0
    assert args.clamp_max == 1000
    assert args.output_format == "bed.gz"
