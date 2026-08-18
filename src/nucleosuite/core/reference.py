"""Reference FASTA helpers shared by sequence-aware commands."""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, TypedDict

from nucleosuite.core.regions import resolve_contig_name

if TYPE_CHECKING:
    import pysam

VALID_BASES = frozenset("ACGT")


class ReferenceContext(TypedDict):
    contig: str
    length: int
    start: int
    end: int
    seq: str


def resolve_fasta_contig(fasta: pysam.FastaFile, contig: str) -> str:
    return resolve_contig_name(contig, list(fasta.references), source_label="FASTA index")



def prepare_reference_context(
    fasta: pysam.FastaFile,
    contig: str,
    start: int,
    end: int,
) -> ReferenceContext:
    fasta_contig = resolve_fasta_contig(fasta, contig)
    fasta_length = fasta.get_reference_length(fasta_contig)
    sequence = fasta.fetch(fasta_contig, start, end).upper()
    expected = end - start
    if len(sequence) != expected:
        sequence = sequence[:expected].ljust(expected, "N")
    return {
        "contig": fasta_contig,
        "length": fasta_length,
        "start": start,
        "end": end,
        "seq": sequence,
    }


def extract_reference_sequence(
    fasta: pysam.FastaFile,
    context: ReferenceContext,
    seq_start: int,
    seq_end: int,
) -> Optional[str]:
    if seq_start < 0 or seq_end <= seq_start or seq_end > context["length"]:
        return None

    if seq_start >= context["start"] and seq_end <= context["end"]:
        relative_start = seq_start - context["start"]
        relative_end = seq_end - context["start"]
        sequence = context["seq"][relative_start:relative_end]
    else:
        sequence = fasta.fetch(context["contig"], seq_start, seq_end).upper()

    if len(sequence) != seq_end - seq_start:
        return None
    return sequence


def sequence_is_acgt(sequence: str) -> bool:
    return bool(sequence) and all(base in VALID_BASES for base in sequence)
