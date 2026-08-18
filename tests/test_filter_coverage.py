from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np

from nucleosuite.filter_coverage import (
    BigWigCoverageReader,
    bed_position,
    default_output_path,
    filter_bed_by_coverage,
)


def test_default_output_name_contains_filter_threshold(tmp_path):
    path = tmp_path / "sample_nucleosome_regions.bed"
    assert default_output_path(path, 2) == tmp_path / "sample_nucleosome_regions_coverage_ge2.bed"
    assert default_output_path(path, 2.5) == tmp_path / "sample_nucleosome_regions_coverage_ge2.5.bed"


def test_bed_position_defaults_to_midpoint_or_uses_one_based_column():
    fields = ["chr1", "10", "21", "peak", "1", ".", "13", "14"]
    assert bed_position(fields, None) == 15
    assert bed_position(fields, 7) == 13


class FakeBigWig:
    def __init__(self, values):
        self.array = np.asarray(values, dtype=float)
        self.closed = False
        self.calls = 0

    def chroms(self):
        return {"chr1": len(self.array)}

    def values(self, chrom, start, end, numpy=False):
        assert chrom == "chr1"
        self.calls += 1
        result = self.array[start:end]
        return result.copy() if numpy else result.tolist()

    def close(self):
        self.closed = True


def test_bigwig_reader_resolves_chr_alias_and_caches_chunks():
    values = np.zeros(100, dtype=float)
    values[12] = 3.0
    values[13] = 4.0
    handle = FakeBigWig(values)
    reader = BigWigCoverageReader(handle, chunk_size=20)
    assert reader.value("1", 12) == (3.0, False)
    assert reader.value("1", 13) == (4.0, False)
    assert handle.calls == 1


def test_filter_coverage_retains_threshold_or_greater_and_preserves_rows(
    tmp_path, monkeypatch
):
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "# example\n"
        "chr1\t10\t20\tp1\t5\t.\t12\t13\n"
        "chr1\t30\t40\tp2\t6\t.\t35\t36\n"
        "chr1\t50\t60\tp3\t7\t.\t55\t56\n",
        encoding="utf-8",
    )
    bigwig = tmp_path / "coverage.bw"
    bigwig.touch()
    values = np.zeros(100, dtype=float)
    values[12] = 1.99
    values[35] = 2.0
    values[55] = 3.0
    fake = FakeBigWig(values)
    monkeypatch.setitem(sys.modules, "pyBigWig", types.SimpleNamespace(open=lambda _: fake))

    output, summary_path, summary = filter_bed_by_coverage(
        bed,
        bigwig,
        coverage_threshold=2,
        position_column=7,
        chunk_size=100,
    )

    assert output.name == "peaks_coverage_ge2.bed"
    assert output.read_text(encoding="utf-8") == (
        "# example\n"
        "chr1\t30\t40\tp2\t6\t.\t35\t36\n"
        "chr1\t50\t60\tp3\t7\t.\t55\t56\n"
    )
    assert summary.total_peaks == 3
    assert summary.retained_peaks == 2
    assert summary.filtered_peaks == 1
    assert "coverage_threshold\t2" in summary_path.read_text(encoding="utf-8")
    assert fake.closed


def test_filter_coverage_treats_missing_bigwig_value_as_zero(tmp_path, monkeypatch):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t10\t20\n", encoding="utf-8")
    bigwig = tmp_path / "coverage.bw"
    bigwig.touch()
    values = np.zeros(100, dtype=float)
    values[15] = np.nan
    fake = FakeBigWig(values)
    monkeypatch.setitem(sys.modules, "pyBigWig", types.SimpleNamespace(open=lambda _: fake))

    output, _, summary = filter_bed_by_coverage(
        bed, bigwig, coverage_threshold=1, chunk_size=100
    )
    assert output.read_text(encoding="utf-8") == ""
    assert summary.missing_values == 1
