"""Contig token parsing and genomic chunk construction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Iterator, Sequence


@dataclass(frozen=True)
class ProcessingRegion:
    contig: str
    adjusted_start: int
    adjusted_end: int
    original_start: int
    original_end: int


def strip_chr_prefix(name: str) -> str:
    return name[3:] if name.lower().startswith("chr") else name


def candidate_contig_aliases(name: str) -> list[str]:
    """Return conservative aliases for a genomic contig name.

    Exact names are always preferred.  The only automatic substitutions are an
    optional ``chr`` prefix and common mitochondrial aliases.
    """
    aliases = [name]
    stripped = strip_chr_prefix(name)
    aliases.extend([stripped, f"chr{stripped}"])
    mito = {
        "M": ("MT", "chrM", "chrMT"),
        "MT": ("M", "chrM", "chrMT"),
        "chrM": ("M", "MT", "chrMT"),
        "chrMT": ("M", "MT", "chrM"),
    }
    aliases.extend(mito.get(name, ()))
    aliases.extend(mito.get(stripped, ()))
    output: list[str] = []
    seen: set[str] = set()
    for alias in aliases:
        if alias and alias not in seen:
            seen.add(alias)
            output.append(alias)
    return output


def resolve_contig_name(
    requested: str,
    references: Sequence[str],
    *,
    source_label: str = "input",
) -> str:
    """Resolve a contig against an input while rejecting ambiguous aliases."""
    reference_set = set(references)
    if requested in reference_set:
        return requested
    matches = [alias for alias in candidate_contig_aliases(requested) if alias in reference_set]
    matches = list(dict.fromkeys(matches))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise KeyError(
            f"Contig '{requested}' is ambiguous in {source_label}; matching names: "
            + ", ".join(matches)
        )
    raise KeyError(f"Contig '{requested}' was not found in {source_label}.")


def resolve_mapping_contig(requested: str, mapping, *, source_label: str = "input") -> str:
    """Resolve a contig against mapping keys."""
    return resolve_contig_name(requested, list(mapping.keys()), source_label=source_label)

def canonical_contig_key(name: str) -> str:
    """Return a comparison key for optional chr prefixes and mitochondrial aliases."""
    stripped = strip_chr_prefix(name)
    if stripped.upper() in {"M", "MT"}:
        return "MT"
    return stripped


def _numeric_contigs(references: Sequence[str]) -> list[str]:
    """Return numeric contigs in header order; used for generic ``autosomes``."""
    return [name for name in references if strip_chr_prefix(name).isdigit()]


def _expand_numeric_range(token: str) -> list[str] | None:
    if ":" in token or "-" not in token:
        return None
    left, right = token.split("-", 1)
    left_value = strip_chr_prefix(left)
    right_value = strip_chr_prefix(right)
    if not (left_value.isdigit() and right_value.isdigit()):
        return None
    start, end = int(left_value), int(right_value)
    step = 1 if start <= end else -1
    return [str(value) for value in range(start, end + step, step)]


def expand_contig_tokens(
    raw_tokens: Sequence[str] | None,
    references: Sequence[str],
) -> list[str]:
    """Expand comma lists, numeric ranges, ``autosomes`` and ``all``.

    Unlike a human-specific implementation, ``autosomes`` means every BAM
    contig whose name becomes entirely numeric after removing an optional
    ``chr`` prefix. No fixed chromosome count is assumed.
    """
    if not raw_tokens:
        return list(references)

    tokens: list[str] = []
    for raw in raw_tokens:
        tokens.extend(part.strip() for part in raw.split(",") if part.strip())

    requested: list[str] = []
    for token in tokens:
        lowered = token.lower()
        if lowered == "all":
            requested.extend(references)
            continue
        if lowered == "autosomes":
            requested.extend(_numeric_contigs(references))
            continue
        expanded = _expand_numeric_range(token)
        if expanded is not None:
            requested.extend(expanded)
            continue
        requested.append(token)

    resolved: list[str] = []
    seen: set[str] = set()
    for item in requested:
        if ":" in item:
            contig_part, coordinates = item.split(":", 1)
            contig = resolve_contig_name(contig_part, references)
            value = f"{contig}:{coordinates}"
        else:
            value = resolve_contig_name(item, references)
        if value not in seen:
            seen.add(value)
            resolved.append(value)
    return resolved


def parse_contig_spec(spec: str, reference_length: int) -> tuple[int, int]:
    if ":" not in spec:
        return 0, reference_length
    _contig, coordinates = spec.split(":", 1)
    start_string, end_string = coordinates.replace(",", "").split("-", 1)
    start, end = int(start_string), int(end_string)
    if not 0 <= start < end <= reference_length:
        raise ValueError(
            f"Invalid interval {spec}; contig length is {reference_length:,}."
        )
    return start, end


def split_into_regions(
    contig: str,
    start: int,
    end: int,
    contig_length: int,
    max_length: int = 100_000,
    overlap: int = 1_000,
) -> list[ProcessingRegion]:
    if max_length < 1:
        raise ValueError("max_length must be positive")
    if overlap < 0:
        raise ValueError("overlap must be non-negative")

    regions: list[ProcessingRegion] = []
    current = start
    while current < end:
        core_start = current
        core_end = min(current + max_length, end)
        regions.append(
            ProcessingRegion(
                contig=contig,
                adjusted_start=max(0, core_start - overlap),
                adjusted_end=min(contig_length, core_end + overlap),
                original_start=core_start,
                original_end=core_end,
            )
        )
        current = core_end
    return regions


def build_processing_regions(
    selected_specs: Sequence[str],
    references: Sequence[str],
    lengths: Sequence[int],
    chunk_bp: int,
    overlap_bp: int,
) -> tuple[list[ProcessingRegion], list[str]]:
    """Build regions in reference-header and coordinate order."""
    length_by_name = dict(zip(references, lengths))
    reference_order = {name: index for index, name in enumerate(references)}

    def sort_key(spec: str):
        contig = spec.split(":", 1)[0]
        start = 0
        if ":" in spec:
            start = int(spec.split(":", 1)[1].replace(",", "").split("-", 1)[0])
        return reference_order[contig], start

    ordered_specs = sorted(selected_specs, key=sort_key)
    regions: list[ProcessingRegion] = []
    selected_names: list[str] = []
    for spec in ordered_specs:
        contig = spec.split(":", 1)[0]
        contig_length = int(length_by_name[contig])
        start, end = parse_contig_spec(spec, contig_length)
        selected_names.append(contig)
        regions.extend(
            split_into_regions(
                contig=contig,
                start=start,
                end=end,
                contig_length=contig_length,
                max_length=chunk_bp,
                overlap=overlap_bp,
            )
        )
    return regions, selected_names
