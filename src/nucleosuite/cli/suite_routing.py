"""Shared BAM-index routing for multicontig suite workflows."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from nucleosuite.core.regions import candidate_contig_aliases


def ensure_bam_index(path: str) -> None:
    """Create a BAM index when neither BAI nor CSI is present."""
    import pysam

    bam_path = Path(path)
    candidates = (
        Path(str(bam_path) + ".bai"),
        bam_path.with_suffix(".bai"),
        Path(str(bam_path) + ".csi"),
    )
    if not any(candidate.exists() for candidate in candidates):
        pysam.index(str(bam_path))


def _resolve_index_contig(requested: str, available: Sequence[str], *, bam_path: str) -> str | None:
    """Resolve a requested FASTA contig against BAM-index names.

    Exact matches are preferred. Otherwise, only an optional ``chr`` prefix and
    common mitochondrial aliases are considered. ``None`` means that the BAM
    does not contain the requested contig.
    """
    available_set = set(available)
    if requested in available_set:
        return requested

    matches = [alias for alias in candidate_contig_aliases(requested) if alias in available_set]
    matches = list(dict.fromkeys(matches))
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(
            f"Contig '{requested}' is ambiguous in BAM index for {bam_path}; "
            f"matching names: {', '.join(matches)}"
        )
    return matches[0]


def route_bams_by_contig(bam_paths: Sequence[str], contigs: Sequence[str]) -> dict[str, list[str]]:
    """Select BAMs with mapped records for each requested contig.

    FASTA/BAM naming differences such as ``chr1`` versus ``1`` and ``chrM``
    versus ``MT`` are resolved conservatively without changing either file.
    """
    import pysam

    routed = {contig: [] for contig in contigs}
    for path in bam_paths:
        ensure_bam_index(path)
        with pysam.AlignmentFile(path, "rb") as bam:
            statistics = {row.contig: int(row.mapped) for row in bam.get_index_statistics()}

        for contig in contigs:
            bam_contig = _resolve_index_contig(contig, tuple(statistics), bam_path=path)
            if bam_contig is not None and statistics[bam_contig] > 0:
                routed[contig].append(path)
    return routed
