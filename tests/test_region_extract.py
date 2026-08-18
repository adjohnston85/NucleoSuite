from __future__ import annotations

import gzip
import os
from pathlib import Path

from nucleosuite.pns_region_extractor import parse_bed


def _write_fake_tool(path: Path, body: str) -> None:
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)


def test_region_extract_parses_plain_bed(tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t10\t20\tregion_a\nchr1\t30\t40\tregion_b\n")

    regions, skipped, max_columns = parse_bed(str(bed))

    assert [region.region_id for region in regions] == ["region_a", "region_b"]
    assert skipped == []
    assert max_columns == 4


def test_region_extract_parses_gzipped_bed(tmp_path):
    bed = tmp_path / "regions.bed.gz"
    with gzip.open(bed, "wt") as handle:
        handle.write("chr1\t10\t20\tregion_a\n")

    regions, skipped, max_columns = parse_bed(str(bed))

    assert len(regions) == 1
    assert regions[0].fields == ("chr1", "10", "20", "region_a")
    assert skipped == []
    assert max_columns == 4


def test_region_extract_parses_bigbed_via_shared_interval_reader(tmp_path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_fake_tool(
        tools / "bigBedToBed",
        "import shutil, sys\nshutil.copyfile(sys.argv[1], sys.argv[2])\n",
    )
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")

    bigbed = tmp_path / "regions.bb"
    bigbed.write_text("chr1\t10\t20\tregion_a\n")

    regions, skipped, max_columns = parse_bed(str(bigbed))

    assert len(regions) == 1
    assert regions[0].region_id == "region_a"
    assert skipped == []
    assert max_columns == 4
