"""Tests for compatible multi-BAM collections."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

import nucleosuite.fragment_lengths as fragment_lengths
from nucleosuite.fragment_lengths import build_parser, default_output_path
from nucleosuite.core.bam_headers import merge_bam_reference_headers


def test_merge_bam_reference_headers_accepts_disjoint_chromosome_bams():
    bam1 = SimpleNamespace(references=("chr1",), lengths=(249250621,))
    bam2 = SimpleNamespace(references=("chr2",), lengths=(243199373,))
    references, lengths = merge_bam_reference_headers([bam1, bam2])
    assert references == ["chr1", "chr2"]
    assert lengths == [249250621, 243199373]


def test_merge_bam_reference_headers_accepts_shared_compatible_headers():
    bam1 = SimpleNamespace(references=("chr1", "chr2"), lengths=(100, 200))
    bam2 = SimpleNamespace(references=("chr1", "chr2"), lengths=(100, 200))
    references, lengths = merge_bam_reference_headers([bam1, bam2])
    assert references == ["chr1", "chr2"]
    assert lengths == [100, 200]


def test_merge_bam_reference_headers_rejects_conflicting_lengths():
    bam1 = SimpleNamespace(references=("chr1",), lengths=(100,))
    bam2 = SimpleNamespace(references=("chr1",), lengths=(101,))
    with pytest.raises(ValueError, match="conflicting lengths"):
        merge_bam_reference_headers([bam1, bam2])


def test_fragment_lengths_parser_accepts_multiple_bams():
    args = build_parser().parse_args(["--bam", "chr1.bam", "chr2.bam"])
    assert args.bamfiles == ["chr1.bam", "chr2.bam"]


def test_fragment_lengths_default_name_for_multiple_bams():
    assert str(default_output_path(["chr1.bam", "chr2.bam"], None)) == "combined_bams_fragment_lengths.tsv"


def test_fragment_length_counts_are_aggregated_across_disjoint_bams(monkeypatch):
    class Read:
        is_paired = True
        is_read1 = True
        is_unmapped = False
        mate_is_unmapped = False
        is_secondary = False
        is_supplementary = False
        is_qcfail = False
        mapping_quality = 60
        is_duplicate = False
        is_proper_pair = True
        template_length = 145
        next_reference_id = 0
        reference_id = 0
        reference_start = 100
        next_reference_start = 200

        def __init__(self, reference_name):
            self.reference_name = reference_name

    class Bam:
        def __init__(self, path):
            self.path = path
            self.references = ("chr1",) if path == "chr1.bam" else ("chr2",)

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def fetch(self, until_eof=False):
            return iter([Read(self.references[0])])

    fake_pysam = SimpleNamespace(AlignmentFile=lambda path, mode: Bam(path))
    monkeypatch.setattr(fragment_lengths, "pysam", fake_pysam)

    counts, summary = fragment_lengths.count_fragment_lengths(
        ["chr1.bam", "chr2.bam"],
        contigs={"chr1", "chr2"},
        min_length=100,
        max_length=200,
    )
    assert counts["all"][145] == 2
    assert summary.fragments_counted == 2


def test_merge_bam_headers_prefers_chr_namespace_and_records_source_aliases():
    from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases

    bam1 = SimpleNamespace(references=("1", "20", "X"), lengths=(100, 200, 300))
    bam2 = SimpleNamespace(references=("chr1", "chr20", "chrY"), lengths=(100, 200, 400))
    merged = merge_bam_reference_headers_with_aliases([bam1, bam2])
    assert merged.references == ["chr1", "chr20", "X", "chrY"]
    assert merged.lengths == [100, 200, 300, 400]
    assert merged.source_contigs[0] == {"chr1": "1", "chr20": "20", "X": "X"}
    assert merged.source_contigs[1] == {"chr1": "chr1", "chr20": "chr20", "chrY": "chrY"}


def test_merge_bam_headers_rejects_aliases_within_one_header():
    from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases

    bam = SimpleNamespace(references=("20", "chr20"), lengths=(200, 200))
    with pytest.raises(ValueError, match="ambiguous equivalent contigs"):
        merge_bam_reference_headers_with_aliases([bam])
