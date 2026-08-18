from __future__ import annotations

import os
from pathlib import Path

from nucleosuite.io import open_text
from nucleosuite.io.intervals import convert_bed_to_bigbed, finalise_interval_files


def _write_fake_tool(path: Path, body: str) -> None:
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)


def test_bed_to_bigbed_sorts_and_clamps_scores(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bedToBigBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[-3], sys.argv[-1])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")

    bed = tmp_path / "peaks.bed"
    bed.write_text(
        "chr2\t5\t10\tb\t2000\t.\t7\t8\n"
        "chr1\t20\t30\ta\t-4\t.\t25\t26\n"
        "chr1\t2\t8\tc\t50\t.\t4\t5\n"
        "chr1\t40\t50\td\t49.6\t.\t45\t46\n"
    )
    output = convert_bed_to_bigbed(
        bed,
        [("chr1", 100), ("chr2", 100)],
    )
    assert output == tmp_path / "peaks.bb"
    lines = output.read_text().splitlines()
    assert [line.split("\t")[:3] for line in lines] == [
        ["chr1", "2", "8"],
        ["chr1", "20", "30"],
        ["chr1", "40", "50"],
        ["chr2", "5", "10"],
    ]
    assert [line.split("\t")[4] for line in lines] == ["50", "0", "50", "1000"]


def test_bigbed_input_materialises_to_text(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bigBedToBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[1], sys.argv[2])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")
    bigbed = tmp_path / "regions.bb"
    bigbed.write_text("chr1\t1\t5\tregion\n")
    with open_text(bigbed) as handle:
        assert handle.read() == "chr1\t1\t5\tregion\n"


def test_bigbed_only_removes_bed_after_conversion(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bedToBigBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[-3], sys.argv[-1])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t1\t5\tGENE1\n")
    outputs = finalise_interval_files([bed], "bigbed", [("chr1", 100)])
    assert outputs == [tmp_path / "genes.bb"]
    assert not bed.exists()
    assert (tmp_path / "genes.bb").exists()


def test_bed_to_bigbed_rewrites_chr_alias_to_chrom_sizes_namespace(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bedToBigBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[-3], sys.argv[-1])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")
    bed = tmp_path / "alias.bed"
    bed.write_text("20\t2\t8\tpeak\t10\n")
    output = convert_bed_to_bigbed(bed, [("chr20", 100)])
    assert output.read_text().splitlines()[0].startswith("chr20\t2\t8\t")


def test_bigbed_score_multiplier_applies_before_rounding(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bedToBigBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[-3], sys.argv[-1])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")
    bed = tmp_path / "pns_peaks.bed"
    bed.write_text(
        "chr1\t1\t5\tpeak1\t0.123456\t.\t3\t4\n"
        "chr1\t10\t20\tpeak2\t12.000000\t.\t15\t16\n"
    )
    output = convert_bed_to_bigbed(
        bed,
        [("chr1", 100)],
        bigbed_score_multiplier=100.0,
    )
    scores = [line.split("\t")[4] for line in output.read_text().splitlines()]
    assert scores == ["12", "1000"]
    # Text BED retains its six-decimal floating scores.
    assert "0.123456" in bed.read_text()
    assert "12.000000" in bed.read_text()


def test_pns_default_1000_scale_maps_fractional_scores_to_bed_range(tmp_path, monkeypatch):
    tools = tmp_path / "tools_default_scale"
    tools.mkdir()
    _write_fake_tool(
        tools / "bedToBigBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[-3], sys.argv[-1])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH','')}")
    bed = tmp_path / "pns_default_scale.bed"
    bed.write_text(
        "chr1\t1\t5\tpeak1\t0.123456\t.\t3\t4\n"
        "chr1\t10\t20\tpeak2\t1.000000\t.\t15\t16\n"
    )
    output = convert_bed_to_bigbed(
        bed,
        [("chr1", 100)],
        bigbed_score_multiplier=1000.0,
    )
    scores = [line.split("\t")[4] for line in output.read_text().splitlines()]
    assert scores == ["123", "1000"]

