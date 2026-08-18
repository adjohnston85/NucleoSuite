"""Utilities for combining reference headers across BAM collections."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from nucleosuite.core.regions import canonical_contig_key


@dataclass(frozen=True)
class MergedBamHeader:
    """Canonical union of one or more BAM reference headers.

    ``references`` and ``lengths`` define the output namespace.  Each entry in
    ``source_contigs`` maps a canonical output name to the corresponding name
    in one BAM header.  Missing keys mean that BAM does not contain the contig.
    """

    references: list[str]
    lengths: list[int]
    source_contigs: list[dict[str, str]]


def _choose_output_name(names: Sequence[str]) -> str:
    """Prefer a ``chr`` spelling when equivalent BAM headers are mixed."""
    for name in names:
        if name.lower().startswith("chr"):
            return name
    return names[0]


def merge_bam_reference_headers_with_aliases(
    bamfiles: Sequence[object],
) -> MergedBamHeader:
    """Return a BAM-derived canonical contig union and per-BAM aliases.

    Equivalent names that differ only by an optional ``chr`` prefix (plus the
    package's conservative mitochondrial aliases) are merged.  If any BAM uses
    a ``chr`` spelling, that spelling becomes the output name.  Equivalent
    contigs must have identical lengths.  A single BAM containing two aliases
    for the same canonical contig is rejected as ambiguous.
    """
    if not bamfiles:
        raise ValueError("At least one BAM header is required")

    key_order: list[str] = []
    names_by_key: dict[str, list[str]] = {}
    lengths_by_key: dict[str, int] = {}
    source_by_key: list[dict[str, str]] = []

    for bam_index, bamfile in enumerate(bamfiles):
        current: dict[str, str] = {}
        references = tuple(bamfile.references)
        raw_lengths = getattr(bamfile, "lengths", (0,) * len(references))
        for reference, raw_length in zip(references, raw_lengths):
            reference = str(reference)
            length = int(raw_length)
            key = canonical_contig_key(reference)
            previous_name = current.get(key)
            if previous_name is not None and previous_name != reference:
                raise ValueError(
                    f"BAM input {bam_index + 1} contains ambiguous equivalent "
                    f"contigs {previous_name!r} and {reference!r}."
                )
            current[key] = reference
            if key not in lengths_by_key:
                key_order.append(key)
                names_by_key[key] = [reference]
                lengths_by_key[key] = length
            else:
                if lengths_by_key[key] != length:
                    names = names_by_key[key] + [reference]
                    raise ValueError(
                        "Equivalent BAM contigs "
                        + ", ".join(repr(name) for name in dict.fromkeys(names))
                        + f" report conflicting lengths: {lengths_by_key[key]} and {length}."
                    )
                if reference not in names_by_key[key]:
                    names_by_key[key].append(reference)
        source_by_key.append(current)

    if not key_order:
        raise ValueError("BAM inputs do not contain any reference sequences")

    canonical_by_key = {key: _choose_output_name(names_by_key[key]) for key in key_order}
    references = [canonical_by_key[key] for key in key_order]
    lengths = [lengths_by_key[key] for key in key_order]
    source_contigs = [
        {canonical_by_key[key]: source_name for key, source_name in current.items()}
        for current in source_by_key
    ]
    return MergedBamHeader(references, lengths, source_contigs)


def merge_bam_reference_headers(bamfiles: Sequence[object]) -> tuple[list[str], list[int]]:
    """Compatibility wrapper returning only the canonical union and lengths."""
    merged = merge_bam_reference_headers_with_aliases(bamfiles)
    return merged.references, merged.lengths
