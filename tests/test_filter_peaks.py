from __future__ import annotations

import gzip
import os
from pathlib import Path

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
