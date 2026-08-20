from pathlib import Path

import numpy as np
import pytest

from nucleosuite.mean_scale import (
    _bigwig_nonzero_mean,
    _default_output,
    _iter_scaled_intervals,
    _region_score_mean,
    build_parser,
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
