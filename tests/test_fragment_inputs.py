from __future__ import annotations

import gzip
from pathlib import Path

from nucleosuite.core.fragment_inputs import IntervalFragmentSource
from nucleosuite.fragments_command import main as fragments_main
from nucleosuite.cli.main import build_parser
from nucleosuite.dcc import build_parser as build_dcc_parser, validate_args


def test_interval_fragment_source_accepts_extra_columns_and_gzip(tmp_path: Path):
    first = tmp_path / "a.bed"
    first.write_text(
        "chr1\t10\t20\tname\t999\t+\n"
        "chr1\t10\t20\tduplicate\n"
        "chr1\t30\t45\n"
    )
    second = tmp_path / "b.bed.gz"
    with gzip.open(second, "wt") as handle:
        handle.write("chr1\t10\t20\tmore\tcolumns\tare\tignored\n")
        handle.write("chr2\t5\t15\n")

    source = IntervalFragmentSource([str(first), str(second)])
    try:
        assert source.references == ["chr1", "chr2"]
        assert source.fetch(
            "chr1", 0, 100,
            max_per_coordinate=1,
            subsample=None,
            dedup_scope="all_bams",
        ) == [(10, 20), (30, 45)]
        assert source.fetch(
            "chr1", 0, 100,
            max_per_coordinate=1,
            subsample=None,
            dedup_scope="per_bam",
        ) == [(10, 20), (10, 20), (30, 45)]
    finally:
        source.close()


def test_fragments_command_combines_interval_inputs(tmp_path: Path):
    a = tmp_path / "a.bed"
    b = tmp_path / "b.bed"
    a.write_text("chr1\t10\t20\tA\nchr1\t10\t20\tA2\n")
    b.write_text("chr1\t10\t20\tB\nchr1\t30\t45\tB2\n")
    prefix = tmp_path / "combined"
    assert fragments_main([
        "--fragments", str(a), str(b),
        "--output-prefix", str(prefix),
        "--max-duplicates", "1",
        "--output-format", "bed",
    ]) == 0
    prefix = tmp_path / "combined_fragmin1_fragmax1000_mapq0_maxdup1_dedupall-bams_subsamplenone_seednone"
    assert Path(f"{prefix}.fragments.bed").read_text().splitlines() == [
        "chr1\t10\t20",
        "chr1\t30\t45",
    ]
    assert Path(f"{prefix}.fragment_length_distribution.png").exists()


def test_pns_parser_accepts_fragment_bed_input():
    parser = build_parser()
    args = parser.parse_args(["pns", "--fragments", "sample.bed.gz"])
    assert args.fragment_files == ["sample.bed.gz"]
    assert args.bamfiles is None


def test_dcc_bam_mode_accepts_fragment_inputs():
    parser = build_dcc_parser()
    args = parser.parse_args([
        "bam",
        "--fragments-a", "a.bed",
        "--fragments-b", "b.bed.gz",
        "--length-a", "147",
        "--length-b", "167",
        "--chrom-sizes", "genome.sizes",
    ])
    assert args.fragments_a == ["a.bed"]
    assert args.fragments_b == ["b.bed.gz"]


def test_dcc_fragment_mode_can_mix_bam_and_interval_inputs():
    parser = build_dcc_parser()
    args = parser.parse_args([
        "bam",
        "--bam-a", "a.bam",
        "--fragments-b", "b.bed.gz",
        "--length-a", "147",
        "--length-b", "167",
        "--chrom-sizes", "genome.sizes",
    ])
    validate_args(args, parser)
    assert args.bam_a == ["a.bam"]
    assert args.fragments_b == ["b.bed.gz"]


def test_fragments_command_streams_multiple_contigs_to_combined_bed(tmp_path: Path):
    source = tmp_path / "fragments.bed"
    source.write_text("chr1\t10\t20\textra\nchr2\t5\t18\tignored\n")
    prefix = tmp_path / "out"
    assert fragments_main([
        "--fragments", str(source),
        "--output-prefix", str(prefix),
        "--output-format", "bed",
    ]) == 0
    prefix = tmp_path / "out_fragmin1_fragmax1000_mapq0_maxdup0_dedupall-bams_subsamplenone_seednone"
    assert Path(f"{prefix}.fragments.bed").read_text().splitlines() == [
        "chr1\t10\t20",
        "chr2\t5\t18",
    ]


def test_interval_fragment_source_normalizes_sqlite_contigs_to_supplied_namespace(tmp_path: Path):
    sizes = tmp_path / "analysis.chrom.sizes"
    sizes.write_text("chr20\t1000\nchrX\t500\n")
    a = tmp_path / "a.bed"
    b = tmp_path / "b.bed"
    a.write_text("20\t10\t155\n")
    b.write_text("chr20\t10\t155\n")

    source = IntervalFragmentSource([str(a), str(b)], chrom_sizes=str(sizes))
    try:
        assert source.references == ["chr20", "chrX"]
        assert source.fetch(
            "chr20", 0, 1000,
            max_per_coordinate=1,
            subsample=None,
            dedup_scope="all_bams",
        ) == [(10, 155)]
        assert source.fetch(
            "20", 0, 1000,
            max_per_coordinate=1,
            subsample=None,
            dedup_scope="per_bam",
        ) == [(10, 155), (10, 155)]
        aliases = source._connection.execute(
            "SELECT source_chrom, canonical_chrom FROM contig_aliases ORDER BY source_id"
        ).fetchall()
        assert aliases == [("20", "chr20"), ("chr20", "chr20")]
    finally:
        source.close()


def test_bam_fragment_source_honours_supplied_global_chr_namespace(tmp_path: Path, monkeypatch):
    from types import SimpleNamespace
    import nucleosuite.core.fragment_inputs as fragment_inputs

    sizes = tmp_path / "selected.chrom.sizes"
    sizes.write_text("chr20\t1000\n")

    class FakeBam:
        references = ("20",)
        lengths = (1000,)

        def close(self):
            return None

    fake_pysam = SimpleNamespace(AlignmentFile=lambda _path, _mode: FakeBam())
    monkeypatch.setattr(fragment_inputs, "pysam", fake_pysam)
    monkeypatch.setattr(fragment_inputs, "require_bam_indexes", lambda _paths: None)

    source = fragment_inputs.BamFragmentSource(
        ["worker_without_chr.bam"], chrom_sizes=str(sizes)
    )
    try:
        assert source.references == ["chr20"]
        assert source.lengths == [1000]
        assert source.source_contigs == [{"chr20": "20"}]
    finally:
        source.close()


def test_open_fragment_source_uses_tabix_random_access(tmp_path: Path):
    import pytest
    pysam = pytest.importorskip("pysam")
    from nucleosuite.core.fragment_inputs import (
        IndexedIntervalFragmentSource,
        open_fragment_source,
    )

    plain = tmp_path / "fragments.bed"
    plain.write_text(
        "20\t10\t155\n"
        "20\t200\t365\n"
        "21\t5\t150\n",
        encoding="utf-8",
    )
    compressed = Path(
        pysam.tabix_index(
            str(plain),
            preset="bed",
            force=True,
            keep_original=True,
        )
    )
    sizes = tmp_path / "analysis.chrom.sizes"
    sizes.write_text("chr20\t1000\nchr21\t800\n", encoding="utf-8")

    source = open_fragment_source(
        fragment_paths=[str(compressed)],
        chrom_sizes=str(sizes),
    )
    try:
        assert isinstance(source, IndexedIntervalFragmentSource)
        assert source.references == ["chr20", "chr21"]
        assert source.fetch(
            "chr20",
            0,
            180,
            max_per_coordinate=0,
            subsample=None,
            dedup_scope="all_bams",
        ) == [(10, 155)]
        assert source.fetch(
            "21",
            0,
            800,
            max_per_coordinate=0,
            subsample=None,
            dedup_scope="all_bams",
        ) == [(5, 150)]
    finally:
        source.close()


def test_tabix_fragment_source_requires_reference_lengths(tmp_path: Path):
    import pytest
    pysam = pytest.importorskip("pysam")
    from nucleosuite.core.fragment_inputs import open_fragment_source

    plain = tmp_path / "fragments.bed"
    plain.write_text("chr1\t1\t10\n", encoding="utf-8")
    compressed = Path(
        pysam.tabix_index(
            str(plain), preset="bed", force=True, keep_original=True
        )
    )
    with pytest.raises(ValueError, match="requires --chrom-sizes or --fasta"):
        open_fragment_source(fragment_paths=[str(compressed)])
