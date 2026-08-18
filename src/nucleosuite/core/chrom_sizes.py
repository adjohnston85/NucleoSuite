"""Chromosome-size input handling for text tables, BAM files and CRAM files."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence
import gzip

_ALIGNMENT_SUFFIXES = {".bam", ".cram"}


def is_alignment_path(path: str | Path) -> bool:
    """Return ``True`` when *path* is named as a BAM or CRAM alignment file."""

    return Path(path).suffix.lower() in _ALIGNMENT_SUFFIXES


def _read_alignment_header(
    path: str | Path,
    *,
    reference_fasta: str | Path | None = None,
) -> list[tuple[str, int]]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - runtime dependency
        raise RuntimeError("pysam is required to read chromosome sizes from BAM/CRAM") from exc

    source = Path(path)
    mode = "rc" if source.suffix.lower() == ".cram" else "rb"
    kwargs: dict[str, str] = {}
    if reference_fasta is not None:
        kwargs["reference_filename"] = str(reference_fasta)
    try:
        handle = pysam.AlignmentFile(str(source), mode, **kwargs)
    except (OSError, ValueError) as exc:
        if source.suffix.lower() == ".cram" and reference_fasta is None:
            raise OSError(
                f"Unable to read CRAM header from {source}. Provide the matching reference FASTA."
            ) from exc
        raise
    try:
        rows = [(str(name), int(length)) for name, length in zip(handle.references, handle.lengths)]
    finally:
        handle.close()
    if not rows:
        raise ValueError(f"No reference sequences were found in the alignment header: {source}")
    return rows


def _read_text_table(path: str | Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    seen: dict[str, int] = {}
    source = Path(path)
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    with opener(source, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected chromosome and length")
            name = fields[0]
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: chromosome length must be an integer"
                ) from exc
            if length < 1:
                raise ValueError(f"{path}:{line_number}: chromosome length must be positive")
            if name in seen:
                if seen[name] != length:
                    raise ValueError(f"{path}:{line_number}: conflicting length for {name}")
                continue
            seen[name] = length
            rows.append((name, length))
    if not rows:
        raise ValueError(f"No chromosome sizes were found in {path}")
    return rows


def read_chrom_sizes_source(
    source: str | Path,
    *,
    reference_fasta: str | Path | None = None,
) -> list[tuple[str, int]]:
    """Read chromosome sizes from a two-column table or alignment header.

    BAM and CRAM inputs are detected by their file suffix. Reference order is
    preserved for alignment files and for text tables.
    """

    path = Path(source)
    if not path.exists():
        raise FileNotFoundError(f"Chromosome-size source was not found: {path}")
    if is_alignment_path(path):
        return _read_alignment_header(path, reference_fasta=reference_fasta)
    return _read_text_table(path)


def filter_chrom_sizes(
    rows: Sequence[tuple[str, int]],
    contig_tokens: Sequence[str] | None,
) -> list[tuple[str, int]]:
    """Filter chromosome sizes using NucleoSuite contig selectors."""

    if not contig_tokens:
        return list(rows)
    from nucleosuite.core.regions import expand_contig_tokens

    references = [name for name, _length in rows]
    selected = expand_contig_tokens(list(contig_tokens), references)
    if any(":" in item for item in selected):
        raise ValueError("chrom-sizes accepts whole contigs, not coordinate intervals")
    selected_set = set(selected)
    return [(name, length) for name, length in rows if name in selected_set]


def write_chrom_sizes_table(
    rows: Iterable[tuple[str, int]],
    output: str | Path,
) -> Path:
    """Write a two-column chromosome-size table."""

    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    materialised = [(str(name), int(length)) for name, length in rows]
    if not materialised:
        raise ValueError("No chromosome sizes remained after filtering")
    with destination.open("wt", encoding="utf-8") as handle:
        for name, length in materialised:
            handle.write(f"{name}\t{length}\n")
    return destination


def chromosome_size_dict(
    source: str | Path | None,
    *,
    reference_fasta: str | Path | None = None,
) -> dict[str, int]:
    """Read *source* and return an insertion-ordered chromosome-size mapping."""

    if source is None:
        return {}
    return dict(read_chrom_sizes_source(source, reference_fasta=reference_fasta))
