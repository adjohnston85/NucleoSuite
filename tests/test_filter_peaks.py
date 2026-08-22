from __future__ import annotations

import gzip
import os
import sys
import types
from pathlib import Path

import numpy as np
import pytest

from nucleosuite import filter_peaks


def _write_fake_tool(path: Path, body: str) -> None:
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)


def test_filter_peaks_absolute_score_and_length_preserves_bed_records(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr1\t0\t100\ta\t-5\t+\textraA\n"
        "chr1\t100\t250\tb\t10\t-\textraB\n"
        "chr1\t250\t450\tc\t20\t.\textraC\n"
        "chr1\t450\t700\td\t30\t.\textraD\n"
    )
    output, summary, result = filter_peaks.filter_peaks(
        bed, min_score=10, max_score=20, min_length=120, max_length=220
    )
    assert output.suffix == ".bed"
    assert output.name == "peaks_filtered_scoremin10_scoremax20_lenmin120_lenmax220.bed"
    assert output.read_text().splitlines() == [
        "chr1\t100\t250\tb\t10\t-\textraB",
        "chr1\t250\t450\tc\t20\t.\textraC",
    ]
    assert summary.is_file()
    assert result.valid_records == 4
    assert result.length_filtered_records == 2
    assert result.retained_records == 2


def test_filter_peaks_percentile_is_calculated_after_length_filter(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr1\t0\t50\tshort_high\t1000\n"
        "chr1\t100\t200\ta\t1\n"
        "chr1\t300\t400\tb\t2\n"
        "chr1\t500\t600\tc\t3\n"
    )
    output, _summary, result = filter_peaks.filter_peaks(
        bed, score_percentile=50, min_length=100, max_length=100
    )
    assert result.percentile_threshold == pytest.approx(2.0)
    assert [line.split("\t")[3] for line in output.read_text().splitlines()] == ["b", "c"]


def test_filter_peaks_abs_score_is_opt_in_and_scaling_changes_bed_score(tmp_path):
    bed = tmp_path / "negative.bed"
    bed.write_text(
        "chr1\t0\t100\tneg\t-4\n"
        "chr1\t100\t200\tpos\t3\n"
    )
    sparse, _summary, _result = filter_peaks.filter_peaks(bed, min_score=3.5)
    assert [line.split("\t")[3] for line in sparse.read_text().splitlines()] == []

    absolute, _summary, _result = filter_peaks.filter_peaks(
        bed, min_score=3.5, abs_score=True, score_scale=10
    )
    fields = absolute.read_text().splitlines()[0].split("\t")
    assert fields[3] == "neg"
    assert fields[4] == "40"


def test_filter_peaks_bed_gz_defaults_to_bed_gz(tmp_path):
    source = tmp_path / "peaks.bed.gz"
    with gzip.open(source, "wt") as handle:
        handle.write("chr1\t0\t100\ta\t5\nchr1\t100\t200\tb\t10\n")
    output, _summary, _result = filter_peaks.filter_peaks(source, min_score=10)
    assert output.name.endswith(".bed.gz")
    with gzip.open(output, "rt") as handle:
        assert handle.read().splitlines() == ["chr1\t100\t200\tb\t10"]


def test_filter_peaks_bigbed_output_scales_rounds_and_clamps_scores(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bedToBigBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[-3], sys.argv[-1])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr1\t0\t100\tneg\t-2\n"
        "chr1\t100\t200\tmid\t3\n"
        "chr1\t200\t300\thigh\t20\n"
    )
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t1000\n")
    output, _summary, result = filter_peaks.filter_peaks(
        bed,
        output_format="bigbed",
        chrom_sizes=sizes,
        abs_score=True,
        score_scale=100,
    )
    assert output.suffix == ".bb"
    assert result.retained_records == 3
    scores = [line.split("\t")[4] for line in output.read_text().splitlines()]
    assert scores == ["200", "300", "1000"]


def test_filter_peaks_bigbed_input_defaults_back_to_bigbed(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bigBedToBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[1], sys.argv[2])\n",
    )
    _write_fake_tool(
        tools / "bedToBigBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[-3], sys.argv[-1])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")
    monkeypatch.setattr(filter_peaks, "infer_bigbed_chrom_sizes", lambda _path: {"chr1": 1000})
    source = tmp_path / "peaks.bb"
    source.write_text("chr1\t0\t100\ta\t5\nchr1\t100\t200\tb\t10\n")
    output, _summary, _result = filter_peaks.filter_peaks(source, min_score=10)
    assert output.suffix == ".bb"
    assert output.read_text().splitlines() == ["chr1\t100\t200\tb\t10"]


def test_filter_peaks_rejects_absolute_and_percentile_score_filters_together(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t0\t100\ta\t5\n")
    with pytest.raises(ValueError, match="cannot be combined"):
        filter_peaks.filter_peaks(bed, min_score=2, score_percentile=50)


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


def _install_fake_bigwig(monkeypatch, *handles):
    queue = list(handles)
    monkeypatch.setitem(
        sys.modules,
        "pyBigWig",
        types.SimpleNamespace(open=lambda _path: queue.pop(0)),
    )


def test_filter_peaks_bed_position_defaults_to_midpoint_or_selected_column():
    fields = ["chr1", "10", "21", "peak", "1", ".", "13", "14"]
    assert filter_peaks.bed_position(fields, None) == 15
    assert filter_peaks.bed_position(fields, 7) == 13


def test_filter_peaks_bigwig_reader_resolves_chr_alias_and_caches_chunks():
    values = np.zeros(100, dtype=float)
    values[12] = 3.0
    values[13] = 4.0
    handle = FakeBigWig(values)
    reader = filter_peaks.BigWigCoverageReader(handle, chunk_size=20)
    assert reader.value("1", 12) == (3.0, False)
    assert reader.value("1", 13) == (4.0, False)
    assert handle.calls == 1


def test_filter_peaks_coverage_filter_supports_bed3_and_preserves_rows(tmp_path, monkeypatch):
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
    _install_fake_bigwig(monkeypatch, FakeBigWig(values), FakeBigWig(values))

    output, summary_path, summary = filter_peaks.filter_peaks(
        bed,
        coverage_bigwig=bigwig,
        min_coverage=2,
        coverage_position_column=7,
        coverage_chunk_size=100,
    )

    assert output.name == "peaks_filtered_covcoverage_covmin2_covposcol7.bed"
    assert output.read_text(encoding="utf-8") == (
        "# example\n"
        "chr1\t30\t40\tp2\t6\t.\t35\t36\n"
        "chr1\t50\t60\tp3\t7\t.\t55\t56\n"
    )
    assert summary.coverage_filtered_records == 1
    assert summary.retained_records == 2
    summary_text = summary_path.read_text(encoding="utf-8")
    assert "min_coverage\t2" in summary_text
    assert "coverage_position_source\tbed_column_7" in summary_text


def test_filter_peaks_treats_missing_coverage_as_zero(tmp_path, monkeypatch):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t10\t20\n", encoding="utf-8")
    bigwig = tmp_path / "coverage.bw"
    bigwig.touch()
    values = np.zeros(100, dtype=float)
    values[15] = np.nan
    _install_fake_bigwig(monkeypatch, FakeBigWig(values), FakeBigWig(values))

    output, _summary_path, summary = filter_peaks.filter_peaks(
        bed,
        coverage_bigwig=bigwig,
        min_coverage=1,
        coverage_chunk_size=100,
    )
    assert output.read_text(encoding="utf-8") == ""
    assert summary.missing_coverage_values == 1
    assert summary.coverage_filtered_records == 1


def test_filter_peaks_percentile_is_after_length_and_coverage_filters(tmp_path, monkeypatch):
    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr1\t0\t100\ta\t1\n"
        "chr1\t100\t200\tb\t2\n"
        "chr1\t200\t300\tcovered_outlier\t1000\n"
        "chr1\t300\t400\tc\t3\n"
        "chr1\t400\t450\tshort\t5000\n",
        encoding="utf-8",
    )
    bigwig = tmp_path / "coverage.bw"
    bigwig.touch()
    values = np.ones(500, dtype=float) * 10
    values[250] = 0
    _install_fake_bigwig(monkeypatch, FakeBigWig(values), FakeBigWig(values))

    output, _summary_path, summary = filter_peaks.filter_peaks(
        bed,
        min_length=100,
        max_length=100,
        coverage_bigwig=bigwig,
        min_coverage=5,
        score_percentile=50,
        coverage_chunk_size=500,
    )
    assert summary.percentile_threshold == pytest.approx(2.0)
    assert summary.length_filtered_records == 1
    assert summary.coverage_filtered_records == 1
    assert [line.split("\t")[3] for line in output.read_text().splitlines()] == ["b", "c"]


def test_filter_peaks_requires_coverage_track_and_threshold_together(tmp_path):
    bed = tmp_path / "peaks.bed"
    bed.write_text("chr1\t0\t100\ta\t5\n")
    with pytest.raises(ValueError, match="supplied together"):
        filter_peaks.filter_peaks(bed, min_coverage=2)
