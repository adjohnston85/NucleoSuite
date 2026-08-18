from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys

from nucleosuite.chrom_sizes_command import main as chrom_sizes_main
from nucleosuite.core.chrom_sizes import read_chrom_sizes_source


class _FakeAlignmentFile:
    def __init__(self, path: str, mode: str, **kwargs: str) -> None:
        self.path = path
        self.mode = mode
        self.kwargs = kwargs
        self.references = ("chr2", "chr1", "chrX")
        self.lengths = (2000, 1000, 500)

    def close(self) -> None:
        return None

    def __enter__(self) -> "_FakeAlignmentFile":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


def _install_fake_pysam(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "pysam", SimpleNamespace(AlignmentFile=_FakeAlignmentFile))


def test_read_chrom_sizes_from_bam_preserves_header_order(tmp_path: Path, monkeypatch) -> None:
    bam = tmp_path / "empty.bam"
    bam.touch()
    _install_fake_pysam(monkeypatch)
    assert read_chrom_sizes_source(bam) == [
        ("chr2", 2000),
        ("chr1", 1000),
        ("chrX", 500),
    ]


def test_chrom_sizes_command_filters_contigs(tmp_path: Path, monkeypatch) -> None:
    bam = tmp_path / "empty.bam"
    output = tmp_path / "selected.chrom.sizes"
    bam.touch()
    _install_fake_pysam(monkeypatch)
    assert chrom_sizes_main([
        "--bam", str(bam),
        "--output", str(output),
        "--contigs", "chr1,chrX",
    ]) == 0
    assert output.read_text() == "chr1\t1000\nchrX\t500\n"


def test_text_chrom_sizes_remain_supported(tmp_path: Path) -> None:
    source = tmp_path / "input.chrom.sizes"
    source.write_text("chr1\t1000\nchr2\t2000\n")
    assert read_chrom_sizes_source(source) == [("chr1", 1000), ("chr2", 2000)]
