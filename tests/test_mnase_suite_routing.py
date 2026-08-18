"""Tests for contig-aware BAM routing in multicontig mnase-suite runs."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys

from nucleosuite.cli.mnase_suite import _derive_sample_name, _route_bams_by_contig


class _FakeAlignmentFile:
    mapped_by_path: dict[str, dict[str, int]] = {}

    def __init__(self, path: str, mode: str):
        self.path = path

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def get_index_statistics(self):
        return [
            SimpleNamespace(contig=contig, mapped=mapped)
            for contig, mapped in self.mapped_by_path[self.path].items()
        ]


def test_routes_chromosome_specific_and_whole_genome_bams_by_index_statistics(tmp_path, monkeypatch):
    chr1 = tmp_path / "merged_chr1.bam"
    chr2 = tmp_path / "merged_chr2.bam"
    whole = tmp_path / "replicate.bam"
    for path in (chr1, chr2, whole):
        path.touch()
        Path(str(path) + ".bai").touch()

    _FakeAlignmentFile.mapped_by_path = {
        str(chr1): {"chr1": 10, "chr2": 0},
        str(chr2): {"chr1": 0, "chr2": 12},
        str(whole): {"chr1": 20, "chr2": 18},
    }
    fake_pysam = SimpleNamespace(AlignmentFile=_FakeAlignmentFile, index=lambda path: None)
    monkeypatch.setitem(sys.modules, "pysam", fake_pysam)

    routed = _route_bams_by_contig([str(chr1), str(chr2), str(whole)], ["chr1", "chr2"])

    assert routed["chr1"] == [str(chr1), str(whole)]
    assert routed["chr2"] == [str(chr2), str(whole)]


def test_sample_name_is_derived_before_contig_suffix_is_added():
    paths = ["/data/merged_chr1.bam", "/data/merged_chr2.bam", "/data/merged_chr10.bam"]
    assert _derive_sample_name(paths) == "merged"
    assert _derive_sample_name(paths, "Gaffney 32") == "Gaffney_32"


def test_sample_name_strips_fragment_extensions_and_contigs():
    paths = ["/data/merged_chr1.bed.gz", "/data/merged_chr2.bed.gz"]
    assert _derive_sample_name(paths, fallback="multi_fragments") == "merged"


def test_routes_chr_prefixed_fasta_contigs_to_unprefixed_bam_contigs(tmp_path, monkeypatch):
    bam_path = tmp_path / "BH01.bam"
    bam_path.touch()
    Path(str(bam_path) + ".bai").touch()

    _FakeAlignmentFile.mapped_by_path = {
        str(bam_path): {"1": 100, "2": 80, "X": 25, "MT": 10},
    }
    fake_pysam = SimpleNamespace(AlignmentFile=_FakeAlignmentFile, index=lambda path: None)
    monkeypatch.setitem(sys.modules, "pysam", fake_pysam)

    routed = _route_bams_by_contig(
        [str(bam_path)],
        ["chr1", "chr2", "chrX", "chrM"],
    )

    assert routed == {
        "chr1": [str(bam_path)],
        "chr2": [str(bam_path)],
        "chrX": [str(bam_path)],
        "chrM": [str(bam_path)],
    }


def test_exact_bam_contig_name_is_preferred_over_alias(tmp_path, monkeypatch):
    bam_path = tmp_path / "mixed_names.bam"
    bam_path.touch()
    Path(str(bam_path) + ".bai").touch()

    _FakeAlignmentFile.mapped_by_path = {
        str(bam_path): {"1": 100, "chr1": 0},
    }
    fake_pysam = SimpleNamespace(AlignmentFile=_FakeAlignmentFile, index=lambda path: None)
    monkeypatch.setitem(sys.modules, "pysam", fake_pysam)

    routed = _route_bams_by_contig([str(bam_path)], ["chr1"])

    assert routed["chr1"] == []


def test_validation_consumes_cores_before_bash_preflight(monkeypatch):
    import nucleosuite.cli.mnase_suite as suite

    observed: dict[str, list[str]] = {}

    def fake_run(command, **kwargs):
        observed["command"] = list(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(suite.subprocess, "run", fake_run)
    suite.validate_argv(["--bam", "sample.bam", "--cores", "4", "--fasta", "genome.fa"])

    assert "--cores" not in observed["command"]
    assert "4" not in observed["command"]
    assert "--validate-only" in observed["command"]


def test_validation_consumes_equals_cores_before_bash_preflight(monkeypatch):
    import nucleosuite.cli.mnase_suite as suite

    observed: dict[str, list[str]] = {}

    def fake_run(command, **kwargs):
        observed["command"] = list(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(suite.subprocess, "run", fake_run)
    suite.validate_argv(["--bam", "sample.bam", "--cores=4", "--fasta", "genome.fa"])

    assert not any(token.startswith("--cores") for token in observed["command"])
    assert "--validate-only" in observed["command"]


def test_validation_consumes_combine_bigwig_method_before_bash_preflight(monkeypatch):
    import nucleosuite.cli.mnase_suite as suite

    observed: dict[str, list[str]] = {}

    def fake_run(command, **kwargs):
        observed["command"] = list(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(suite.subprocess, "run", fake_run)
    suite.validate_argv(
        [
            "--bam", "sample.bam",
            "--fasta", "genome.fa",
            "--combine-bigwig-method", "bedgraph",
        ]
    )

    assert "--combine-bigwig-method" not in observed["command"]
    assert "bedgraph" not in observed["command"]
    assert "--validate-only" in observed["command"]


def test_analysis_scope_defaults_to_combined_only_and_is_consumed():
    import nucleosuite.cli.mnase_suite as suite

    scope, remaining = suite._extract_analysis_scope(
        ["--bam", "sample.bam", "--analysis-scope", "per-contig-and-combined"]
    )
    assert scope == "per-contig-and-combined"
    assert "--analysis-scope" not in remaining
    assert suite._extract_analysis_scope(["--bam", "sample.bam"])[0] == "combined-only"


def test_cfdna_analysis_scope_defaults_to_combined_only():
    import nucleosuite.cli.cfdna_suite as suite

    scope, remaining = suite._extract_analysis_scope(
        ["--analysis-scope=combined-only", "--bam", "sample.bam"]
    )
    assert scope == "combined-only"
    assert not any(token.startswith("--analysis-scope") for token in remaining)
