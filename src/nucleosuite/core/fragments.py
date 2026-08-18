"""Shared paired-end fragment extraction and BAM validation."""

from __future__ import annotations

import os
import random
from collections import defaultdict
from typing import DefaultDict, Dict, Iterable, Iterator, Optional, Tuple

try:
    import pysam
except ImportError:  # allow command help/tests without runtime dependency
    pysam = None

Fragment = Tuple[int, int]
FragmentKey = Tuple[str, int, int]


def is_softclipped_or_padded(cigartuples) -> bool:
    """Return ``True`` when a CIGAR contains S, H or P operations."""
    if not cigartuples:
        return False
    return any(op in (4, 5, 6) for op, _length in cigartuples)


def generate_paired_reads(
    bamfile: pysam.AlignmentFile,
    contig: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    max_per_coordinate: int = 1,
    subsample: Optional[float] = None,
    fragment_counts: Optional[DefaultDict[FragmentKey, int]] = None,
    min_mapq: int = 0,
    require_proper_pair: bool = False,
    include_marked_duplicates: bool = False,
    include_qcfail: bool = False,
    allow_softclipped: bool = False,
    include_secondary: bool = False,
    include_supplementary: bool = False,
    dedup_contig: Optional[str] = None,
):
    """Yield accepted read pairs as ``(forward_read, reverse_read)``.

    Coordinate-based deduplication is applied to complete fragment coordinates.
    ``max_per_coordinate=1`` keeps one fragment per coordinate and ``0`` disables
    coordinate-based deduplication. A shared ``fragment_counts`` mapping applies
    the limit across multiple BAM inputs; omitting it applies the limit only to
    this generator invocation.
    """
    if max_per_coordinate < 0:
        raise ValueError("max_per_coordinate must be 0 or greater")
    if subsample is not None and not 0.0 <= subsample <= 1.0:
        raise ValueError("subsample must be between 0 and 1")
    if min_mapq < 0:
        raise ValueError("min_mapq must be 0 or greater")

    unpaired: Dict[str, pysam.AlignedSegment] = {}
    counts = fragment_counts if fragment_counts is not None else defaultdict(int)

    try:
        iterator = bamfile.fetch(contig, start, end, multiple_iterators=True)
    except (ValueError, OSError):
        return

    for read in iterator:
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_secondary and not include_secondary:
            continue
        if read.is_supplementary and not include_supplementary:
            continue
        if read.is_duplicate and not include_marked_duplicates:
            continue
        if read.is_qcfail and not include_qcfail:
            continue
        if require_proper_pair and not read.is_proper_pair:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if not allow_softclipped and is_softclipped_or_padded(read.cigartuples):
            continue
        if read.reference_end is None or read.next_reference_start is None:
            continue

        query_name = read.query_name
        mate = unpaired.pop(query_name, None)
        if mate is None:
            unpaired[query_name] = read
            continue

        if mate.is_unmapped or mate.mate_is_unmapped:
            continue
        if mate.is_secondary and not include_secondary:
            continue
        if mate.is_supplementary and not include_supplementary:
            continue
        if mate.is_duplicate and not include_marked_duplicates:
            continue
        if mate.is_qcfail and not include_qcfail:
            continue
        if require_proper_pair and not mate.is_proper_pair:
            continue
        if mate.mapping_quality < min_mapq:
            continue
        if not allow_softclipped and is_softclipped_or_padded(mate.cigartuples):
            continue
        if mate.reference_end is None or mate.next_reference_start is None:
            continue
        if read.reference_name != mate.reference_name:
            continue
        if read.is_reverse == mate.is_reverse:
            continue

        fragment_start = min(read.reference_start, mate.reference_start)
        fragment_end = max(read.reference_end, mate.reference_end)
        if fragment_end <= fragment_start:
            continue

        key = (dedup_contig or read.reference_name, fragment_start, fragment_end)
        if max_per_coordinate > 0 and counts[key] >= max_per_coordinate:
            continue
        counts[key] += 1

        if subsample is not None and random.random() > subsample:
            continue

        if read.is_reverse:
            yield mate, read
        else:
            yield read, mate


def generate_fragment_ranges(
    bamfile: pysam.AlignmentFile,
    contig: str,
    fetch_start: int,
    fetch_end: int,
    max_per_coordinate: int,
    subsample: Optional[float],
    fragment_counts: Optional[DefaultDict[FragmentKey, int]] = None,
    min_mapq: int = 0,
    require_proper_pair: bool = False,
    include_marked_duplicates: bool = False,
    include_qcfail: bool = False,
    allow_softclipped: bool = False,
    include_secondary: bool = False,
    include_supplementary: bool = False,
    dedup_contig: Optional[str] = None,
) -> Iterator[Fragment]:
    """Convert accepted read pairs into half-open genomic fragments."""
    for forward, reverse in generate_paired_reads(
        bamfile=bamfile,
        contig=contig,
        start=fetch_start,
        end=fetch_end,
        max_per_coordinate=max_per_coordinate,
        subsample=subsample,
        fragment_counts=fragment_counts,
        min_mapq=min_mapq,
        require_proper_pair=require_proper_pair,
        include_marked_duplicates=include_marked_duplicates,
        include_qcfail=include_qcfail,
        allow_softclipped=allow_softclipped,
        include_secondary=include_secondary,
        include_supplementary=include_supplementary,
        dedup_contig=dedup_contig,
    ):
        if forward.is_reverse:
            forward, reverse = reverse, forward
        if forward.reference_start > reverse.reference_start:
            continue
        if reverse.reference_end < forward.reference_end:
            continue
        yield forward.reference_start, reverse.reference_end


def collect_fragments(
    bamfiles: Iterable[pysam.AlignmentFile],
    contig: str,
    start: int,
    end: int,
    max_per_coordinate: int,
    subsample: Optional[float],
    dedup_scope: str = "all_bams",
    source_contigs: Optional[Iterable[Optional[str]]] = None,
) -> list[Fragment]:
    """Collect fragments from BAMs using a common deduplication policy."""
    if dedup_scope not in {"all_bams", "per_bam"}:
        raise ValueError("dedup_scope must be 'all_bams' or 'per_bam'")

    bam_list = list(bamfiles)
    if source_contigs is None:
        source_names: list[Optional[str]] = [contig] * len(bam_list)
    else:
        source_names = list(source_contigs)
        if len(source_names) != len(bam_list):
            raise ValueError("source_contigs must contain one entry per BAM")

    shared_counts = defaultdict(int) if dedup_scope == "all_bams" else None
    fragments: list[Fragment] = []
    for bamfile, source_contig in zip(bam_list, source_names):
        if source_contig is None:
            continue
        counts = shared_counts if dedup_scope == "all_bams" else None
        fragments.extend(
            generate_fragment_ranges(
                bamfile=bamfile,
                contig=source_contig,
                fetch_start=start,
                fetch_end=end,
                max_per_coordinate=max_per_coordinate,
                subsample=subsample,
                fragment_counts=counts,
                dedup_contig=contig,
            )
        )
    return fragments


def bam_index_candidates(bam_path: str) -> tuple[str, str, str]:
    absolute = os.path.abspath(bam_path)
    directory = os.path.dirname(absolute)
    basename = os.path.splitext(os.path.basename(absolute))[0]
    return (
        absolute + ".bai",
        os.path.join(directory, basename + ".bai"),
        absolute + ".csi",
    )


def require_bam_indexes(bam_paths: Iterable[str]) -> None:
    """Raise a detailed error if any BAM index is absent."""
    missing = []
    for path in bam_paths:
        candidate1, candidate2, candidate3 = bam_index_candidates(path)
        if not any(os.path.exists(candidate) for candidate in (candidate1, candidate2, candidate3)):
            missing.append((os.path.abspath(path), candidate1, candidate2, candidate3))

    if not missing:
        return

    lines = [
        "Missing BAM index (.bai) for the following BAM file(s):",
        "",
    ]
    for bam, candidate1, candidate2, candidate3 in missing:
        lines.extend(
            [
                f"  BAM: {bam}",
                f"    expected: {candidate1}",
                f"         or : {candidate2}",
                f"         or : {candidate3}",
                "",
            ]
        )
    lines.extend(["Create each index with:", "  samtools index <file.bam>"])
    raise FileNotFoundError("\n".join(lines))
