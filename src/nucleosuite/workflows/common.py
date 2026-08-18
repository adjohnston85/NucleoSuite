"""Shared workflow setup and fragment randomisation."""

from __future__ import annotations

import os
import random
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
try:
    import pysam
except ImportError:  # interval-only commands can run without BAM/FASTA support
    pysam = None

from nucleosuite.core.fragment_inputs import FragmentSource, open_fragment_source
from nucleosuite.core.randomization import (
    build_dinuc_index,
    dinuc_anchor_randomize_fragments,
    uniform_randomize_fragments,
)
from nucleosuite.core.reference import ReferenceContext, prepare_reference_context
from nucleosuite.core.regions import ProcessingRegion, build_processing_regions, expand_contig_tokens
from nucleosuite.io.tracks import build_bigwig_header


@dataclass
class FragmentRunContext:
    source: FragmentSource
    fasta: Optional[pysam.FastaFile]
    regions: list[ProcessingRegion]
    selected_contigs: list[str]
    bigwig_header: list[tuple[str, int]]

    def collect(
        self,
        *,
        contig: str,
        start: int,
        end: int,
        max_duplicates: int,
        subsample: float | None,
        dedup_scope: str,
    ) -> list[tuple[int, int]]:
        return self.source.fetch(
            contig,
            start,
            end,
            max_per_coordinate=max_duplicates,
            subsample=subsample,
            dedup_scope=dedup_scope,
        )

    def close(self) -> None:
        try:
            self.source.close()
        finally:
            if self.fasta is not None:
                try:
                    self.fasta.close()
                except Exception:
                    pass


def ensure_coordinate_order(records):
    """Sort a fragment chunk only when its coordinate order contains an inversion.

    Fragment pairs recovered from a coordinate-sorted BAM can be emitted slightly
    out of start-coordinate order because a pair is not complete until both mates
    have been observed.  Checking first avoids an unnecessary sort for chunks that
    are already ordered.
    """

    if any(
        (records[index][0], records[index][1])
        > (records[index + 1][0], records[index + 1][1])
        for index in range(len(records) - 1)
    ):
        records.sort(key=lambda record: (record[0], record[1]))
    return records


def set_random_seed(seed: Optional[int]) -> None:
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)


def default_output_prefix(input_paths: Sequence[str], contigs=None) -> str:
    basenames = []
    for path in input_paths:
        name = os.path.basename(path)
        for suffix in (".bed.gz", ".bigbed", ".bigBed", ".bam", ".bed", ".bb", ".gz"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        basenames.append(name)
    prefix = "_".join(basenames)
    if contigs and len(contigs) == 1:
        prefix += "_" + contigs[0].replace(":", "_").replace(",", "-")
    return prefix


def input_paths_from_args(args) -> list[str]:
    return list(getattr(args, "bamfiles", None) or getattr(args, "fragment_files", None) or [])


def ensure_output_parent(output_prefix: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)


def prepare_fragment_run(
    *,
    bam_paths: Sequence[str] | None,
    fragment_paths: Sequence[str] | None,
    fasta_path: Optional[str],
    chrom_sizes_path: Optional[str],
    contig_tokens,
    chunk_bp: int,
    overlap_bp: int,
    blacklist_path: str | None = None,
) -> FragmentRunContext:
    if fasta_path and pysam is None:
        raise RuntimeError("pysam is required when --fasta is supplied")
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None
    source: FragmentSource | None = None
    try:
        source = open_fragment_source(
            bam_paths=bam_paths,
            fragment_paths=fragment_paths,
            chrom_sizes=chrom_sizes_path,
            fasta=fasta,
            blacklist_path=blacklist_path,
        )
        references = source.references
        lengths = source.lengths
        selected_specs = expand_contig_tokens(contig_tokens, references)
        regions, selected_names = build_processing_regions(
            selected_specs=selected_specs,
            references=references,
            lengths=lengths,
            chunk_bp=chunk_bp,
            overlap_bp=overlap_bp,
        )
        header = build_bigwig_header(references, lengths, selected_names)
        return FragmentRunContext(
            source=source,
            fasta=fasta,
            regions=regions,
            selected_contigs=selected_names,
            bigwig_header=header,
        )
    except Exception:
        if source is not None:
            source.close()
        if fasta is not None:
            try:
                fasta.close()
            except Exception:
                pass
        raise


def prepare_bam_run(
    bam_paths: Sequence[str],
    fasta_path: Optional[str],
    contig_tokens,
    chunk_bp: int,
    overlap_bp: int,
    blacklist_path: str | None = None,
) -> FragmentRunContext:
    """Compatibility wrapper retained for code outside the package."""
    return prepare_fragment_run(
        bam_paths=bam_paths,
        fragment_paths=None,
        fasta_path=fasta_path,
        chrom_sizes_path=None,
        contig_tokens=contig_tokens,
        chunk_bp=chunk_bp,
        overlap_bp=overlap_bp,
        blacklist_path=blacklist_path,
    )


def prepare_reference_if_needed(
    fasta: Optional[pysam.FastaFile],
    contig: str,
    start: int,
    end: int,
    required: bool,
) -> Optional[ReferenceContext]:
    if not required:
        return None
    if fasta is None:
        raise ValueError("This operation requires --fasta")
    return prepare_reference_context(fasta, contig, start, end)


def randomize_fragments(
    fragments,
    mode: str,
    start: int,
    end: int,
    reference_context: Optional[ReferenceContext],
    anchor_prob_start: float,
    max_anchor_tries: int,
    fallback: str,
):
    if mode == "none":
        return fragments
    if mode == "uniform":
        return uniform_randomize_fragments(fragments, start, end)
    if mode == "dinuc_anchor":
        if reference_context is None:
            raise ValueError("dinuc_anchor randomisation requires --fasta")
        positions = build_dinuc_index(reference_context["seq"])
        return dinuc_anchor_randomize_fragments(
            fragments=fragments,
            start=start,
            end=end,
            window_sequence=reference_context["seq"],
            dinuc_positions=positions,
            anchor_prob_start=anchor_prob_start,
            max_anchor_tries=max_anchor_tries,
            fallback=fallback,
        )
    raise ValueError(f"Unknown randomisation mode: {mode}")


def print_progress(
    index: int,
    total: int,
    region: ProcessingRegion,
    *,
    previous_contig: str | None = None,
) -> None:
    """Report contig transitions plus occasional chunk milestones."""
    contig_changed = region.contig != previous_contig
    if contig_changed or index == total or index % 1_000 == 0:
        print(
            f"[INFO] Processing {region.contig}: chunk {index:,}/{total:,}; "
            f"region {region.original_start:,}-{region.original_end:,}.",
            flush=True,
        )
