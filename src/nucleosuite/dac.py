#!/usr/bin/env python3
"""
Calculate multiplicative distance autocorrelation (DAC) from one or more
bigWig signal tracks.

The DAC at distance d is calculated as:

    DAC[d] = sum(signal[i] * signal[i + d])

for all valid position pairs within each selected region. Pair products are
never calculated across region boundaries.

Region input modes
------------------
1. BED intervals

   --regions-bed regions.bed

   The first three columns must be chromosome, start and end. Column 4 is used
   as the state/group name by default when present. Regions can be used exactly
   as supplied or converted to strand-aware upstream/downstream windows.

2. Whole-genome or chromosome windows

   --chrom-sizes genome.chrom.sizes --scope combined_chromosomes
   --chrom-sizes genome.chrom.sizes --scope chromosome --chromosome chr1

   Windows are generated from the chromosome sizes file and DAC is calculated
   separately for the assigned state name.

State categorisation
--------------------
BED state/group names are preserved by default. Repeatable ``--category``
rules can combine labels using exact, prefix, or regular-expression matching.

Multiple bigWigs
----------------
By default, matching bigWigs are combined into one DAC result. This is useful
for chromosome-specific bigWigs such as sample_chr*.bw. If multiple bigWigs
contain the same chromosome, each track contributes independently to the sum.
Use --separate-bigwigs to write one set of outputs per bigWig instead.

Chromosome names
----------------
Chromosome names are matched exactly first. If no exact match is available,
common ``chr``-prefix aliases such as ``1`` and ``chr1`` are checked.

Dependencies
------------
    numpy
    pyBigWig

Example
-------
    nucleosuite dac \
        --bigwig sample_chr*.bw \
        --chrom-sizes genome.chrom.sizes \
        --scope combined_chromosomes \
        --dmax 2000 \
        --out-prefix sample
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

from nucleosuite.core.chrom_sizes import read_chrom_sizes_source
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.io import open_text
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.progress import ProgressReporter

import numpy as np


@dataclass(frozen=True)
class Region:
    """A genomic interval assigned to a DAC group/state."""

    chrom: str
    start: int
    end: int
    state: str
    strand: str = "+"
    anchor_start: int | None = None
    anchor_end: int | None = None

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Track:
    """An open bigWig track and its chromosome metadata."""

    path: str
    handle: object
    chrom_sizes: Mapping[str, int]
    chrom_cache: MutableMapping[str, Optional[str]]


@dataclass
class DacStats:
    """Processing statistics for one DAC output."""

    regions_seen: int = 0
    region_track_pairs_used: int = 0
    skipped_missing_chromosome: int = 0
    skipped_short_or_empty: int = 0
    clipped_regions: int = 0
    signal_positions: int = 0
    nonzero_signal_positions: int = 0
    total_signal: float = 0.0
    blacklisted_signal_positions: int = 0
    missing_signal_positions: int = 0


def sanitize_filename(value: str) -> str:
    """Return a filesystem-friendly name while preserving useful detail."""

    value = str(value).strip()
    value = re.sub(r"[\\/\s]+", "_", value)
    value = re.sub(r"[^A-Za-z0-9._+-]+", "_", value)
    value = value.strip("._")
    return value or "unnamed"


def split_comma_values(values: Optional[Sequence[str]]) -> Optional[set[str]]:
    """Expand repeated and comma-separated command-line values."""

    if not values:
        return None

    expanded: set[str] = set()
    for value in values:
        expanded.update(item.strip() for item in value.split(",") if item.strip())

    return expanded or None


def expand_bigwig_inputs(inputs: Sequence[str]) -> List[str]:
    """Expand shell-style glob patterns and remove duplicate paths."""

    files: List[str] = []
    seen: set[str] = set()

    for item in inputs:
        matches = sorted(glob.glob(item))
        candidates = matches if matches else [item]

        for candidate in candidates:
            candidate = os.path.abspath(os.path.expanduser(candidate))
            if not os.path.isfile(candidate):
                continue
            if candidate not in seen:
                seen.add(candidate)
                files.append(candidate)

    if not files:
        raise FileNotFoundError(
            "No bigWig files were found for: " + " ".join(map(str, inputs))
        )

    return files


def read_chrom_sizes(path: str) -> List[Tuple[str, int]]:
    """Read chromosome sizes from a table, BAM header or CRAM header."""

    return list(read_chrom_sizes_source(path))


def make_windows_from_chrom_sizes(
    chrom_sizes_path: str,
    scope: str,
    selected_chromosomes: Optional[set[str]],
    window_size: int,
    state_name: str,
    min_region_length: int,
) -> List[Region]:
    """Generate non-overlapping windows from chromosome sizes."""

    if window_size <= 0:
        raise ValueError("--window-size must be greater than zero.")

    chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    if selected_chromosomes is not None:
        available_names = [chrom for chrom, _size in chrom_sizes]
        resolved_selected: set[str] = set()
        for requested in selected_chromosomes:
            try:
                resolved_selected.add(
                    resolve_contig_name(
                        requested, available_names, source_label="chromosome sizes"
                    )
                )
            except KeyError as exc:
                raise ValueError(str(exc)) from exc
        selected_chromosomes = resolved_selected

    if scope == "chromosome" and not selected_chromosomes:
        raise ValueError(
            "At least one --chromosome is required when --scope chromosome."
        )

    regions: List[Region] = []
    found: set[str] = set()

    for chrom, size in chrom_sizes:
        if selected_chromosomes is not None and chrom not in selected_chromosomes:
            continue

        found.add(chrom)
        for start in range(0, size, window_size):
            end = min(start + window_size, size)
            if end - start < min_region_length:
                continue
            regions.append(Region(chrom, start, end, state_name, "+"))

    if selected_chromosomes:
        missing = selected_chromosomes - found
        if missing:
            raise ValueError(
                "Chromosome(s) not found in chromosome sizes file: "
                + ", ".join(sorted(missing))
            )

    if not regions:
        raise ValueError("No windows remained after chromosome and length filtering.")

    return regions


StateMatcher = Tuple[str, object]
StateCategoryRule = Tuple[str, Tuple[StateMatcher, ...]]


def parse_state_category_rules(specifications: Optional[Sequence[str]]) -> List[StateCategoryRule]:
    """Parse repeatable CATEGORY=MATCHER[,MATCHER...] categorisation rules.

    Matcher forms:
      exact:STATE     exact state-name match
      prefix:PREFIX   state starts with PREFIX
      regex:PATTERN   Python regular-expression search
      STATE           shorthand for exact:STATE

    Rules are applied in command-line order and the first matching category wins.
    States that match no rule retain their original names.
    """

    rules: List[StateCategoryRule] = []

    for specification in specifications or []:
        if "=" not in specification:
            raise ValueError(
                "Invalid --category rule {!r}. Expected "
                "CATEGORY=MATCHER[,MATCHER...].".format(specification)
            )

        category, matcher_text = specification.split("=", 1)
        category = category.strip()
        if not category:
            raise ValueError(
                f"Invalid --category rule {specification!r}: category is empty."
            )

        raw_matchers = [item.strip() for item in matcher_text.split(",") if item.strip()]
        if not raw_matchers:
            raise ValueError(
                f"Invalid --category rule {specification!r}: no matchers supplied."
            )

        parsed_matchers: List[StateMatcher] = []
        for matcher in raw_matchers:
            if matcher.startswith("exact:"):
                value = matcher[len("exact:") :]
                if not value:
                    raise ValueError(
                        f"Invalid empty exact matcher in --category {specification!r}."
                    )
                parsed_matchers.append(("exact", value))

            elif matcher.startswith("prefix:"):
                value = matcher[len("prefix:") :]
                if not value:
                    raise ValueError(
                        f"Invalid empty prefix matcher in --category {specification!r}."
                    )
                parsed_matchers.append(("prefix", value))

            elif matcher.startswith("regex:"):
                pattern = matcher[len("regex:") :]
                if not pattern:
                    raise ValueError(
                        f"Invalid empty regex matcher in --category {specification!r}."
                    )
                try:
                    compiled = re.compile(pattern)
                except re.error as exc:
                    raise ValueError(
                        f"Invalid regex in --category {specification!r}: {exc}"
                    ) from exc
                parsed_matchers.append(("regex", compiled))

            else:
                parsed_matchers.append(("exact", matcher))

        rules.append((category, tuple(parsed_matchers)))

    return rules


def categorize_state_name(state: str, rules: Sequence[StateCategoryRule]) -> str:
    """Return the first user-defined category matching a state name."""

    for category, matchers in rules:
        for match_type, matcher in matchers:
            if match_type == "exact" and state == matcher:
                return category
            if match_type == "prefix" and state.startswith(str(matcher)):
                return category
            if match_type == "regex" and matcher.search(state):
                return category

    return state


def read_regions_bed(
    path: str,
    state_column: Optional[int],
    strand_column: Optional[int],
    state_name: str,
    region_mode: str,
    extend: int,
    strands: str,
    selected_chromosomes: Optional[set[str]],
    min_region_length: int,
    category_rules: Sequence[StateCategoryRule],
) -> List[Region]:
    """Read BED regions and optionally convert them to strand-aware windows."""

    if state_column is not None and state_column < 1:
        raise ValueError("--state-column must be a 1-based column number.")
    if strand_column is not None and strand_column < 1:
        raise ValueError("--strand-column must be a 1-based column number.")
    if region_mode != "interval" and extend <= 0:
        raise ValueError("--extend must be greater than zero for anchored regions.")

    state_index = state_column - 1 if state_column is not None else None
    strand_index = strand_column - 1 if strand_column is not None else None
    regions: List[Region] = []

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}: line {line_number} has fewer than three BED columns."
                )

            chrom = fields[0]
            if selected_chromosomes is not None:
                try:
                    resolve_contig_name(
                        chrom,
                        list(selected_chromosomes),
                        source_label="selected chromosomes",
                    )
                except KeyError:
                    continue

            try:
                feature_start = int(fields[1])
                feature_end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}: line {line_number} has non-integer coordinates."
                ) from exc

            if feature_start < 0 or feature_end <= feature_start:
                continue

            state = state_name
            if state_index is not None and state_index < len(fields):
                state = fields[state_index]
            state = categorize_state_name(str(state), category_rules)

            strand = "+"
            if strand_index is not None and strand_index < len(fields):
                candidate = fields[strand_index]
                if candidate in {"+", "-"}:
                    strand = candidate
                elif region_mode != "interval":
                    continue
            elif region_mode != "interval":
                raise ValueError(
                    f"{path}: line {line_number} does not contain the requested "
                    f"strand column {strand_column}."
                )

            if strands == "plus" and strand != "+":
                continue
            if strands == "minus" and strand != "-":
                continue

            if region_mode == "interval":
                start, end = feature_start, feature_end
            elif region_mode == "downstream":
                if strand == "+":
                    start, end = feature_start, feature_start + extend
                else:
                    start, end = feature_end - extend, feature_end
            elif region_mode == "upstream":
                if strand == "+":
                    start, end = feature_start - extend, feature_start
                else:
                    start, end = feature_end, feature_end + extend
            else:
                raise ValueError(f"Unsupported region mode: {region_mode}")

            start = max(0, start)
            if end - start < min_region_length:
                continue

            regions.append(
                Region(
                    chrom, start, end, state, strand,
                    anchor_start=feature_start,
                    anchor_end=feature_end,
                )
            )

    if not regions:
        raise ValueError("No BED regions remained after filtering.")

    return regions


def _strip_ensembl_version(value: str) -> str:
    """Remove a numeric version suffix from an Ensembl-style identifier."""

    value = value.strip()
    if value.upper().startswith("ENS") and "." in value:
        head, suffix = value.rsplit(".", 1)
        if suffix.isdigit():
            return head
    return value


def expand_gene_queries(
    values: Optional[Sequence[str]],
    list_paths: Optional[Sequence[str]],
) -> List[str]:
    """Read repeated, comma-separated and file-based gene queries."""

    queries: List[str] = []
    seen: set[str] = set()

    def add(value: str) -> None:
        value = value.strip()
        if value and value not in seen:
            seen.add(value)
            queries.append(value)

    for raw in values or []:
        for value in raw.split(","):
            add(value)

    for path in list_paths or []:
        with open_text(path) as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                add(text.split()[0])

    return queries


def _individual_gene_label(
    gene_id: str,
    gene_name: str,
    label_mode: str,
) -> str:
    if label_mode == "id":
        return gene_id
    if label_mode == "name":
        return gene_name or gene_id
    if gene_name and gene_name != gene_id:
        return f"{gene_name}__{gene_id}"
    return gene_id


def read_selected_gene_regions(
    path: str,
    queries: Sequence[str],
    gene_id_column: int,
    gene_name_column: int,
    mode: str,
    pool_name: str,
    label_mode: str,
    selected_chromosomes: Optional[set[str]],
    min_region_length: int,
) -> Tuple[List[Region], List[Mapping[str, object]]]:
    """Select genes by Ensembl identifier or gene name and create DAC groups."""

    if gene_id_column < 1 or gene_name_column < 1:
        raise ValueError("Gene ID and name columns must be 1-based positive integers.")
    if not queries:
        raise ValueError("--genes-bed requires at least one --gene or --gene-list value.")

    id_index = gene_id_column - 1
    name_index = gene_name_column - 1
    query_by_id: Dict[str, List[str]] = defaultdict(list)
    query_by_name: Dict[str, List[str]] = defaultdict(list)
    for query in queries:
        query_by_id[_strip_ensembl_version(query)].append(query)
        query_by_name[query.casefold()].append(query)

    selected: List[Tuple[str, int, int, str, str, str, List[str]]] = []
    matched_queries: set[str] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t") if "\t" in text else text.split()
            if len(fields) < 3:
                raise ValueError(f"{path}: line {line_number} has fewer than three BED columns.")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number} has non-integer coordinates.") from exc
            chrom = fields[0]
            if start < 0 or end <= start or end - start < min_region_length:
                continue
            if selected_chromosomes is not None and chrom not in selected_chromosomes:
                continue

            gene_id = fields[id_index] if id_index < len(fields) else ""
            gene_name = fields[name_index] if name_index < len(fields) else ""
            matched = []
            matched.extend(query_by_id.get(_strip_ensembl_version(gene_id), []))
            matched.extend(query_by_name.get(gene_name.casefold(), []))
            matched = list(dict.fromkeys(matched))
            if not matched:
                continue
            matched_queries.update(matched)
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else "+"
            selected.append((chrom, start, end, gene_id, gene_name, strand, matched))

    unmatched = [query for query in queries if query not in matched_queries]
    if unmatched:
        raise ValueError(
            "Gene query or queries did not match the selected gene BED after filtering: "
            + ", ".join(unmatched)
        )
    if not selected:
        raise ValueError("No selected gene intervals remained after filtering.")

    regions: List[Region] = []
    selection_rows: List[Mapping[str, object]] = []
    individual_labels: set[str] = set()
    for chrom, start, end, gene_id, gene_name, strand, matched in selected:
        label = _individual_gene_label(gene_id, gene_name, label_mode)
        if mode in {"individual", "both"}:
            if label in individual_labels:
                raise ValueError(
                    f"Individual DAC label {label!r} is not unique. Use "
                    "--gene-output-label name-id or id."
                )
            individual_labels.add(label)
        if mode in {"pooled", "both"}:
            regions.append(
                Region(chrom, start, end, pool_name, strand, start, end)
            )
        if mode in {"individual", "both"}:
            regions.append(
                Region(chrom, start, end, label, strand, start, end)
            )
        selection_rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "strand": strand,
                "matched_queries": ",".join(matched),
                "individual_label": label,
                "included_in_pool": "yes" if mode in {"pooled", "both"} else "no",
            }
        )

    return regions, selection_rows


def write_gene_selection_tsv(
    output_path: str,
    rows: Sequence[Mapping[str, object]],
) -> None:
    if not rows:
        return
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def group_regions_by_state(regions: Iterable[Region]) -> Dict[str, List[Region]]:
    """Group regions while retaining deterministic state ordering."""

    grouped: Dict[str, List[Region]] = defaultdict(list)
    for region in regions:
        grouped[region.state].append(region)

    return dict(sorted(grouped.items(), key=lambda item: item[0]))


def candidate_chromosome_aliases(chrom: str) -> List[str]:
    """Return conservative chromosome aliases without assuming a genome."""

    aliases = [chrom]

    if chrom.startswith("chr") and len(chrom) > 3:
        aliases.append(chrom[3:])
    else:
        aliases.append(f"chr{chrom}")

    mitochondrial_aliases = {
        "M": ["MT", "chrM", "chrMT"],
        "MT": ["M", "chrM", "chrMT"],
        "chrM": ["M", "MT", "chrMT"],
        "chrMT": ["M", "MT", "chrM"],
    }
    aliases.extend(mitochondrial_aliases.get(chrom, []))

    unique: List[str] = []
    seen: set[str] = set()
    for alias in aliases:
        if alias not in seen:
            seen.add(alias)
            unique.append(alias)

    return unique


def resolve_bigwig_chromosome(track: Track, chrom: str) -> Optional[str]:
    """Resolve a region chromosome against one bigWig chromosome dictionary."""

    if chrom in track.chrom_cache:
        return track.chrom_cache[chrom]

    resolved: Optional[str] = None
    for candidate in candidate_chromosome_aliases(chrom):
        if candidate in track.chrom_sizes:
            resolved = candidate
            break

    track.chrom_cache[chrom] = resolved
    return resolved


def open_tracks(paths: Sequence[str]) -> List[Track]:
    """Open bigWigs with pyBigWig and return track metadata."""

    try:
        import pyBigWig  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "pyBigWig is required. Install it with 'conda install -c bioconda "
            "pybigwig' or 'pip install pyBigWig'."
        ) from exc

    tracks: List[Track] = []
    try:
        for path in paths:
            handle = pyBigWig.open(path)
            if handle is None or not handle.isBigWig():
                raise ValueError(f"Not a readable bigWig file: {path}")
            tracks.append(
                Track(
                    path=path,
                    handle=handle,
                    chrom_sizes=handle.chroms(),
                    chrom_cache={},
                )
            )
    except Exception:
        close_tracks(tracks)
        raise

    return tracks


def close_tracks(tracks: Iterable[Track]) -> None:
    """Close all open bigWig handles."""

    for track in tracks:
        try:
            track.handle.close()
        except Exception:
            pass


def read_bigwig_region(
    track: Track,
    region: Region,
    value_limit: Optional[float],
) -> Tuple[Optional[np.ndarray], bool]:
    """Read one region, clipping it to the bigWig chromosome boundary."""

    bigwig_chrom = resolve_bigwig_chromosome(track, region.chrom)
    if bigwig_chrom is None:
        return None, False

    chrom_size = int(track.chrom_sizes[bigwig_chrom])
    start = max(0, region.start)
    end = min(region.end, chrom_size)
    clipped = start != region.start or end != region.end

    if end <= start:
        return np.empty(0, dtype=np.float64), clipped

    values = track.handle.values(bigwig_chrom, start, end, numpy=True)
    if values is None:
        values = np.zeros(end - start, dtype=np.float64)
    else:
        values = np.asarray(values, dtype=np.float64)
        np.nan_to_num(values, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

    if value_limit is not None:
        np.clip(values, -value_limit, value_limit, out=values)

    return values, clipped


def opportunity_vector(region_length: int, dmax: int) -> np.ndarray:
    """Return the number of possible pairs at each lag for one region."""

    opportunities = np.zeros(dmax + 1, dtype=np.float64)
    max_lag = min(dmax, region_length - 1)
    if max_lag > 0:
        lags = np.arange(1, max_lag + 1, dtype=np.int64)
        opportunities[1 : max_lag + 1] = region_length - lags
    return opportunities


def opportunity_vector_from_mask(valid: np.ndarray, dmax: int) -> np.ndarray:
    """Count valid base pairs per lag without treating masked bases as zero."""
    mask = np.asarray(valid, dtype=np.float64)
    result = np.zeros(dmax + 1, dtype=np.float64)
    if mask.size < 2 or np.count_nonzero(mask) < 2:
        return result
    max_lag = min(dmax, mask.size - 1)
    fft_length = 1 << (2 * mask.size - 1).bit_length()
    spectrum = np.fft.rfft(mask, n=fft_length)
    correlation = np.fft.irfft(
        spectrum * np.conjugate(spectrum), n=fft_length
    )
    result[1 : max_lag + 1] = np.rint(
        np.maximum(0.0, correlation[1 : max_lag + 1])
    )
    return result


def update_dac_sparse(dac: np.ndarray, values: np.ndarray, dmax: int) -> None:
    """Update DAC using only non-zero signal positions."""

    positions = np.flatnonzero(values)
    if positions.size < 2:
        return

    nonzero_values = values[positions]

    for i in range(positions.size - 1):
        stop = int(
            np.searchsorted(
                positions,
                positions[i] + dmax,
                side="right",
            )
        )
        if stop <= i + 1:
            continue

        distances = positions[i + 1 : stop] - positions[i]
        products = nonzero_values[i] * nonzero_values[i + 1 : stop]
        np.add.at(dac, distances, products)


def update_dac_fft(dac: np.ndarray, values: np.ndarray, dmax: int) -> None:
    """Update DAC using a zero-padded FFT autocorrelation."""

    n = values.size
    max_lag = min(dmax, n - 1)
    if max_lag < 1:
        return

    fft_length = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(values, n=fft_length)
    correlation = np.fft.irfft(
        spectrum * np.conjugate(spectrum),
        n=fft_length,
    )
    dac[1 : max_lag + 1] += correlation[1 : max_lag + 1]


def update_dac(
    dac: np.ndarray,
    values: np.ndarray,
    dmax: int,
    algorithm: str,
    sparse_threshold: float,
) -> str:
    """Select and run the sparse or FFT DAC implementation."""

    if values.size < 2:
        return "none"

    selected = algorithm
    if algorithm == "auto":
        density = float(np.count_nonzero(values)) / float(values.size)
        selected = "sparse" if density <= sparse_threshold else "fft"

    if selected == "sparse":
        update_dac_sparse(dac, values, dmax)
    elif selected == "fft":
        update_dac_fft(dac, values, dmax)
    else:
        raise ValueError(f"Unsupported DAC algorithm: {selected}")

    return selected


def calculate_state_dac(
    tracks: Sequence[Track],
    regions: Sequence[Region],
    dmax: int,
    value_limit: Optional[float],
    min_region_length: int,
    algorithm: str,
    sparse_threshold: float,
    progress_every: int,
    quiet: bool,
    blacklist: BlacklistIndex | None = None,
) -> Tuple[np.ndarray, np.ndarray, DacStats, Dict[str, int]]:
    """Calculate raw DAC and opportunities for one state across tracks."""

    raw_dac = np.zeros(dmax + 1, dtype=np.float64)
    opportunities = np.zeros(dmax + 1, dtype=np.float64)
    stats = DacStats()
    algorithm_counts = {"sparse": 0, "fft": 0, "none": 0}
    opportunity_cache: Dict[int, np.ndarray] = {}
    matching_track_cache: Dict[str, List[Track]] = {}

    for region_index, region in enumerate(regions, start=1):
        stats.regions_seen += 1

        if region.chrom not in matching_track_cache:
            matching_track_cache[region.chrom] = [
                track
                for track in tracks
                if resolve_bigwig_chromosome(track, region.chrom) is not None
            ]

        matching_tracks = matching_track_cache[region.chrom]
        if not matching_tracks:
            stats.skipped_missing_chromosome += 1
            continue

        for track in matching_tracks:
            values, clipped = read_bigwig_region(track, region, value_limit)

            # The chromosome was already resolved when matching_tracks was built.
            if values is None:
                continue

            if clipped:
                stats.clipped_regions += 1

            region_start = max(0, region.start)
            blacklisted = 0
            ordinary_missing = 0
            if blacklist is not None:
                blacklisted = blacklist.mask_values(
                    region.chrom, region_start, values
                )
                stats.blacklisted_signal_positions += blacklisted

            region_length = int(values.size)
            if region_length < min_region_length:
                stats.skipped_short_or_empty += 1
                continue

            valid = np.isfinite(values)
            valid_count = int(np.count_nonzero(valid))
            if valid_count < min_region_length:
                stats.skipped_short_or_empty += 1
                continue
            stats.region_track_pairs_used += 1
            stats.signal_positions += valid_count
            stats.missing_signal_positions += ordinary_missing
            working = np.where(valid, values, 0.0)
            stats.nonzero_signal_positions += int(np.count_nonzero(working))
            stats.total_signal += float(np.sum(working))

            if valid_count == region_length:
                if region_length not in opportunity_cache:
                    opportunity_cache[region_length] = opportunity_vector(
                        region_length, dmax
                    )
                opportunities += opportunity_cache[region_length]
            else:
                opportunities += opportunity_vector_from_mask(valid, dmax)

            selected_algorithm = update_dac(
                raw_dac,
                working,
                dmax,
                algorithm,
                sparse_threshold,
            )
            algorithm_counts[selected_algorithm] += 1

        if (
            not quiet
            and progress_every > 0
            and region_index % progress_every == 0
        ):
            print(
                f"  Processed {region_index:,}/{len(regions):,} regions; "
                f"used {stats.region_track_pairs_used:,} region-track pairs",
                flush=True,
            )

    return raw_dac, opportunities, stats, algorithm_counts


def build_reported_dac(
    raw_dac: np.ndarray,
    opportunities: np.ndarray,
    normalize_dac: bool,
) -> np.ndarray:
    """Return either raw or opportunity-normalized DAC values."""

    if not normalize_dac:
        return raw_dac.copy()

    reported = np.zeros_like(raw_dac, dtype=np.float64)
    valid = opportunities > 0
    reported[valid] = raw_dac[valid] / opportunities[valid]
    return reported


def write_dac_tsv(
    output_path: str,
    raw_dac: np.ndarray,
    opportunities: np.ndarray,
    normalize_dac: bool,
    total_signal: float,
    cpm_scale: float,
) -> None:
    """Write raw, opportunity-normalized, percentage and depth-scaled DAC."""

    reported = build_reported_dac(raw_dac, opportunities, normalize_dac)
    reported_total = float(np.sum(reported[1:]))

    if not math.isclose(reported_total, 0.0, abs_tol=0.0):
        percent = (reported / reported_total) * 100.0
    else:
        percent = np.zeros_like(reported)

    signal_pair_denominator = total_signal * total_signal
    if not math.isclose(signal_pair_denominator, 0.0, abs_tol=0.0):
        per_million_pairs = (raw_dac / signal_pair_denominator) * cpm_scale
    else:
        per_million_pairs = np.zeros_like(raw_dac)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Distance",
                "DAC Value",
                "DAC Value Percent",
                "Raw DAC Value",
                "Opportunities",
                "DAC per million signal-pairs",
            ]
        )

        for distance in range(1, len(raw_dac)):
            writer.writerow(
                [
                    distance,
                    f"{reported[distance]:.12g}",
                    f"{percent[distance]:.12g}",
                    f"{raw_dac[distance]:.12g}",
                    f"{opportunities[distance]:.12g}",
                    f"{per_million_pairs[distance]:.12g}",
                ]
            )


def plot_dac_tsv(
    tsv_path: str,
    output_path: str,
    title: str | None = None,
    *,
    peak_resolution: float = 160.0,
) -> Path | None:
    """Plot the primary DAC values written to a TSV.

    DAC plots show the raw profile only by default. Peak detection is performed
    only when peak labels are explicitly requested with ``--plot-label-points
    peaks``. In that case the plot uses exactly the same resolution-derived
    smoothing and peak caller as :mod:`nucleosuite.nrl`: the raw profile is
    retained in grey behind the local-maximum and detection-smoothed profiles.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    data = np.genfromtxt(tsv_path, names=True, delimiter="\t", dtype=None, encoding="utf-8")
    if getattr(data, "size", 0) == 0:
        return
    names = list(data.dtype.names or ())
    x = np.atleast_1d(data[names[0]]).astype(float)
    y_name = next((name for name in names if name.replace("_", " ") == "DAC Value"), names[1])
    y = np.atleast_1d(data[y_name]).astype(float)
    plot_mask = x != 1
    x_plot = x[plot_mask]
    y_plot = y[plot_mask]
    if x_plot.size == 0:
        return

    from nucleosuite.plotting import annotate_points, get_plot_options, save_figure
    options = get_plot_options()
    detect_and_label = options.label_points == "peaks" and x_plot.size >= 3

    fig, ax = plt.subplots(figsize=(10, 5))
    if detect_and_label:
        from nucleosuite.nrl import (
            call_resolution_peaks,
            moving_average_by_distance,
            resolution_smoothing_windows,
        )
        detection_window, local_window = resolution_smoothing_windows(float(peak_resolution))
        local_values = moving_average_by_distance(x_plot, y_plot, local_window)
        detection_values = moving_average_by_distance(x_plot, y_plot, detection_window)
        called = call_resolution_peaks(
            x_plot, y_plot, local_values, detection_values, float(peak_resolution)
        )

        ax.plot(x_plot, y_plot, color="0.72", linewidth=0.9, label="Unsmoothed")
        ax.plot(
            x_plot, local_values, color="black", linewidth=1.5,
            label=(f"Local maxima ({local_window} bp)" if local_window > 1 else "Local maxima signal"),
        )
        if detection_window != local_window:
            ax.plot(
                x_plot, detection_values, color="0.4", linewidth=1.2, linestyle="--",
                label=(f"Peak detection ({detection_window} bp)" if detection_window > 1 else "Peak detection signal"),
            )
        if called:
            ax.scatter(
                [peak.distance for peak in called],
                [peak.smoothed_value for peak in called],
                s=28, facecolors="white", edgecolors="black", linewidths=1.0,
                zorder=4, label="Called peaks",
            )
            annotate_points(
                ax,
                [peak.distance for peak in called],
                [peak.smoothed_value for peak in called],
                points_are_peaks=True,
                options=options,
            )
        ax.legend(frameon=False)
    else:
        ax.plot(
            x_plot, y_plot, linewidth=1.2, marker="o", markersize=2.0,
            markeredgewidth=0,
        )

    from nucleosuite.plotting import apply_base_pair_x_axis
    apply_base_pair_x_axis(ax, x_plot)
    ax.set_xlabel("Distance (bp)")
    ax.set_ylabel("DAC value")
    if title:
        ax.set_title(title)
    ax.grid(axis="x", alpha=0.5)
    fig.tight_layout()
    saved = save_figure(fig, output_path, default_dpi=220, bbox_inches="tight")
    plt.close(fig)
    return saved


def write_summary_tsv(
    output_path: str,
    rows: Sequence[Mapping[str, object]],
) -> None:
    """Write a compact run summary for all completed state outputs."""

    if not rows:
        return

    fieldnames = list(rows[0].keys())
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def default_output_prefix(
    bigwig_files: Sequence[str],
    region_source_path: Optional[str],
    scope: str,
    selected_chromosomes: Optional[set[str]],
) -> str:
    """Construct an informative default output prefix."""

    if len(bigwig_files) == 1:
        bigwig_part = Path(bigwig_files[0]).stem
    else:
        stems = [Path(path).stem for path in bigwig_files]
        common = os.path.commonprefix(stems).rstrip("._-")
        bigwig_part = common if common else f"combined_{len(bigwig_files)}_bigwigs"

    if region_source_path:
        region_part = Path(region_source_path).stem
    elif scope == "chromosome" and selected_chromosomes:
        region_part = "_".join(sorted(selected_chromosomes))
    else:
        region_part = "combined_chromosomes"

    return sanitize_filename(f"{bigwig_part}_{region_part}")


def process_track_group(
    bigwig_files: Sequence[str],
    regions_by_state: Mapping[str, Sequence[Region]],
    output_prefix: str,
    output_dir: str,
    args: argparse.Namespace,
) -> List[Mapping[str, object]]:
    """Calculate all state outputs for one combined or single-track group."""

    tracks = open_tracks(bigwig_files)
    blacklist = load_blacklist_unbounded(getattr(args, "blacklist_bed", None))
    summary_rows: List[Mapping[str, object]] = []

    try:
        for state, state_regions in regions_by_state.items():
            if not args.quiet:
                print(
                    f"State/group: {state} ({len(state_regions):,} regions; "
                    f"{len(tracks):,} bigWig file(s))",
                    flush=True,
                )

            raw_dac, opportunities, stats, algorithm_counts = calculate_state_dac(
                tracks=tracks,
                regions=state_regions,
                dmax=args.dmax,
                value_limit=args.value_limit,
                min_region_length=args.min_region_length,
                algorithm=args.algorithm,
                sparse_threshold=args.sparse_threshold,
                progress_every=args.progress_every,
                quiet=args.quiet,
                blacklist=blacklist,
            )

            state_part = sanitize_filename(state)
            normalization_suffix = (
                "opportunity_normalized"
                if not args.no_normalize_dac
                else "raw"
            )
            output_name = (
                f"{sanitize_filename(output_prefix)}_{state_part}_"
                f"DAC_{normalization_suffix}.tsv"
            )
            output_path = os.path.join(output_dir, output_name)

            write_dac_tsv(
                output_path=output_path,
                raw_dac=raw_dac,
                opportunities=opportunities,
                normalize_dac=not args.no_normalize_dac,
                total_signal=stats.total_signal,
                cpm_scale=args.cpm_scale,
            )
            plot_dac_tsv(
                output_path,
                os.path.splitext(output_path)[0] + ".png",
                title=f"{output_prefix}: {state}",
                peak_resolution=float(getattr(args, "peak_resolution", 160.0)),
            )

            summary_rows.append(
                {
                    "State": state,
                    "Output": output_path,
                    "BigWig files": len(bigwig_files),
                    "Regions": stats.regions_seen,
                    "Region-track pairs used": stats.region_track_pairs_used,
                    "Missing-chromosome skips": stats.skipped_missing_chromosome,
                    "Short-or-empty skips": stats.skipped_short_or_empty,
                    "Clipped regions": stats.clipped_regions,
                    "Signal positions": stats.signal_positions,
                    "Non-zero signal positions": stats.nonzero_signal_positions,
                    "Total signal": f"{stats.total_signal:.12g}",
                    "Blacklisted signal positions": stats.blacklisted_signal_positions,
                    "Missing signal positions": stats.missing_signal_positions,
                    "Blacklisted anchors excluded": getattr(args, "_blacklisted_anchor_exclusions", 0),
                    "Sparse calculations": algorithm_counts["sparse"],
                    "FFT calculations": algorithm_counts["fft"],
                }
            )

            if not args.quiet:
                print(f"  Wrote: {output_path}")
                print(
                    f"  Region-track pairs used: "
                    f"{stats.region_track_pairs_used:,}"
                )
                print(
                    f"  Non-zero signal positions (overlaps counted): "
                    f"{stats.nonzero_signal_positions:,}"
                )
                print(
                    f"  Algorithms: sparse={algorithm_counts['sparse']:,}, "
                    f"FFT={algorithm_counts['fft']:,}",
                    flush=True,
                )

    finally:
        close_tracks(tracks)

    return summary_rows


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Validate combinations not expressible through argparse alone."""

    if args.dmax < 1:
        parser.error("--dmax must be at least 1.")
    if args.min_region_length < 2:
        parser.error("--min-region-length must be at least 2.")
    if args.value_limit is not None and args.value_limit <= 0:
        parser.error("--value-limit must be greater than zero.")
    if args.cpm_scale <= 0:
        parser.error("--cpm-scale must be greater than zero.")
    if not 0.0 <= args.sparse_threshold <= 1.0:
        parser.error("--sparse-threshold must be between 0 and 1.")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative.")
    if float(getattr(args, "peak_resolution", 160.0)) < 0:
        parser.error("--peak-resolution must be 0 or greater.")

    gene_options_used = bool(args.gene or args.gene_list)
    if args.regions_bed:
        if args.scope not in {"genome", "combined_chromosomes"}:
            parser.error("--scope is only used with --chrom-sizes.")
        if gene_options_used:
            parser.error("--gene and --gene-list require --genes-bed.")
    elif args.genes_bed:
        if args.scope not in {"genome", "combined_chromosomes"}:
            parser.error("--scope is only used with --chrom-sizes.")
        if args.region_mode != "interval":
            parser.error("--region-mode is only used with --regions-bed.")
        if args.category:
            parser.error("--category is only used with --regions-bed.")
        if not gene_options_used:
            parser.error("--genes-bed requires at least one --gene or --gene-list value.")
    else:
        if not args.chrom_sizes:
            parser.error("Provide --regions-bed, --genes-bed or --chrom-sizes.")
        if args.region_mode != "interval":
            parser.error("--region-mode is only used with --regions-bed.")
        if args.category:
            parser.error("--category is only used with --regions-bed.")
        if gene_options_used:
            parser.error("--gene and --gene-list require --genes-bed.")


def build_parser() -> argparse.ArgumentParser:
    """Create the DAC command-line parser."""

    parser = argparse.ArgumentParser(
        prog="nucleosuite dac",
        description=(
            "Calculate multiplicative distance autocorrelation from bigWig "
            "signal within BED regions or genome/chromosome windows."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Whole genome from chromosome-specific bigWigs
  nucleosuite dac --bigwig 'sample_chr*.bw' \\
      --chrom-sizes genome.chrom.sizes --scope combined_chromosomes --dmax 2000

  # One chromosome
  nucleosuite dac --bigwig sample.bw \\
      --chrom-sizes genome.chrom.sizes --scope chromosome \\
      --chromosome chr1

  # DAC separately for states in BED column 4
  nucleosuite dac --bigwig sample.bw --regions-bed states.bed \\
      --state-column 4

  # Optionally combine arbitrary state labels into user-defined categories
  nucleosuite dac --bigwig sample.bw --regions-bed states.bed \\
      --category 'Open=exact:Promoter,exact:Enhancer' \\
      --category 'Repressed=prefix:13_,prefix:14_'

  # Pool selected genes into one DAC
  nucleosuite dac --bigwig sample.bw --genes-bed genes.bed \
      --gene TP53 --gene ENSG00000146648 --gene-dac-mode pooled

  # Produce pooled and individual DACs from a gene list
  nucleosuite dac --bigwig sample.bw --genes-bed genes.bed \
      --gene-list genes.txt --gene-dac-mode both

  # Strand-aware downstream windows
  nucleosuite dac --bigwig sample.bw --regions-bed features.bed \\
      --region-mode downstream --extend 2000 --strand-column 6
""",
    )

    parser.add_argument(
        "--bigwig",
        nargs="+",
        required=True,
        help=(
            "One or more bigWig paths or glob patterns. Quote globs when you "
            "want the script to expand them."
        ),
    )
    parser.add_argument(
        "--blacklist-bed",
        help=(
            "BED blacklist. Overlapping BED/gene anchors are excluded; "
            "blacklisted bases in generated windows are missing opportunities."
        ),
    )

    region_source = parser.add_mutually_exclusive_group(required=True)
    region_source.add_argument(
        "--regions-bed",
        help=(
            "BED file containing analysis regions. Columns 1-3 are required; "
            "column 4 is used as the state by default when present."
        ),
    )
    region_source.add_argument(
        "--genes-bed",
        help=(
            "BED3+ gene intervals selected by Ensembl gene ID or gene name. "
            "Use with --gene and/or --gene-list."
        ),
    )
    region_source.add_argument(
        "--chrom-sizes",
        help=(
            "Two-column chromosome-size table, BAM, or CRAM used to generate "
            "genome or chromosome windows."
        ),
    )

    parser.add_argument(
        "--scope",
        choices=["combined_chromosomes", "genome", "chromosome"],
        default="combined_chromosomes",
        help=(
            "Window scope when --chrom-sizes is used. Default: combined_chromosomes."
        ),
    )
    parser.add_argument(
        "--chromosome",
        action="append",
        default=None,
        help=(
            "Chromosome to include. May be repeated or comma-separated. "
            "Required with --scope chromosome; also filters BED input."
        ),
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=100_000,
        help="Window size for --chrom-sizes mode. Default: 100000.",
    )
    parser.add_argument(
        "--state-name",
        default="Combined chromosomes",
        help=(
            "State assigned to generated windows or BED rows lacking the "
            "requested state column. Default: Combined chromosomes."
        ),
    )
    parser.add_argument(
        "--state-column",
        type=int,
        default=4,
        help=(
            "1-based BED column containing the state/group name. If the "
            "column is absent, --state-name is used. Default: 4."
        ),
    )
    parser.add_argument(
        "--strand-column",
        type=int,
        default=6,
        help=(
            "1-based BED strand column. Required to exist for upstream or "
            "downstream region modes. Default: 6."
        ),
    )
    parser.add_argument(
        "--region-mode",
        choices=["interval", "downstream", "upstream"],
        default="interval",
        help=(
            "Use BED intervals directly or create strand-aware windows from "
            "their boundaries. Default: interval."
        ),
    )
    parser.add_argument(
        "--extend",
        type=int,
        default=2000,
        help=(
            "Window length for upstream/downstream region modes. Default: 2000."
        ),
    )
    parser.add_argument(
        "--strands",
        choices=["plus", "minus", "both"],
        default="both",
        help="BED strands to retain. Default: both.",
    )
    parser.add_argument(
        "--category",
        action="append",
        default=None,
        metavar="CATEGORY=MATCHER[,MATCHER...]",
        help=(
            "Optionally combine BED state names into user-defined categories. "
            "May be repeated. Matchers can be exact:STATE, prefix:PREFIX, "
            "regex:PATTERN, or an unqualified exact state name. Rules are "
            "applied in order; unmatched states retain their original names. "
            "Example: --category "
            "'Open=exact:Promoter,exact:Enhancer' --category "
            "'Repressed=prefix:13_,prefix:14_'."
        ),
    )

    parser.add_argument(
        "--gene",
        action="append",
        default=None,
        help=(
            "Gene to select from --genes-bed by Ensembl ID or gene name. May be "
            "repeated or comma-separated."
        ),
    )
    parser.add_argument(
        "--gene-list",
        action="append",
        default=None,
        help=(
            "Text file containing one Ensembl ID or gene name per line. May be repeated."
        ),
    )
    parser.add_argument(
        "--gene-id-column",
        type=int,
        default=4,
        help="1-based Ensembl gene ID column in --genes-bed. Default: 4.",
    )
    parser.add_argument(
        "--gene-name-column",
        type=int,
        default=5,
        help="1-based gene-name column in --genes-bed. Default: 5.",
    )
    parser.add_argument(
        "--gene-dac-mode",
        choices=["pooled", "individual", "both"],
        default="pooled",
        help=(
            "Pool all selected genes, calculate one DAC per selected gene, or do both. "
            "Default: pooled."
        ),
    )
    parser.add_argument(
        "--gene-pool-name",
        default="selected_genes",
        help="State/group label for pooled selected genes. Default: selected_genes.",
    )
    parser.add_argument(
        "--gene-output-label",
        choices=["name", "id", "name-id"],
        default="name-id",
        help="Label used for individual gene DAC outputs. Default: name-id.",
    )

    parser.add_argument(
        "--dmax",
        type=int,
        default=2000,
        help=(
            "Maximum DAC distance, inclusive. Output distances are 1 through "
            "dmax. Default: 2000."
        ),
    )
    parser.add_argument(
        "--min-region-length",
        type=int,
        default=2,
        help="Minimum clipped region length to analyse. Default: 2.",
    )
    parser.add_argument(
        "--value-limit",
        type=float,
        default=None,
        help="Optional absolute cap applied to bigWig values before DAC.",
    )
    parser.add_argument(
        "--no-normalize-dac",
        action="store_true",
        help=(
            "Report raw multiplicative DAC as DAC Value instead of dividing "
            "by the number of possible base-pair opportunities."
        ),
    )
    parser.add_argument(
        "--cpm-scale",
        type=float,
        default=1_000_000.0,
        help=(
            "Scale used for the DAC per million signal-pairs column. "
            "Default: 1000000."
        ),
    )
    parser.add_argument(
        "--algorithm",
        choices=["auto", "sparse", "fft"],
        default="auto",
        help=(
            "DAC calculation method. Auto uses sparse pair enumeration for "
            "sparse tracks and FFT autocorrelation for denser tracks. "
            "Default: auto."
        ),
    )
    parser.add_argument(
        "--sparse-threshold",
        type=float,
        default=0.10,
        help=(
            "Maximum non-zero fraction at which --algorithm auto uses the "
            "sparse method. Default: 0.10."
        ),
    )

    parser.add_argument(
        "--separate-bigwigs",
        action="store_true",
        help="Write independent outputs for each bigWig instead of combining them.",
    )
    parser.add_argument(
        "--out-prefix",
        default=None,
        help="Output prefix. Default is derived from the bigWig and region input.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for DAC outputs. Default: current directory.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help=(
            "Print progress after this many regions; use 0 to disable. "
            "Default: 1000."
        ),
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress and per-output summaries.",
    )
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive signal analysis. "
            "This budget defaults independently to 1."
        ),
    )

    parser.add_argument(
        "--peak-resolution",
        type=float,
        default=160.0,
        help=(
            "Resolution in bp used only when DAC peak labels are explicitly enabled with "
            "--plot-label-points peaks. Peak detection then uses the same resolution-derived "
            "smoothing and refinement method as `nucleosuite nrl` (default: 160 bp)."
        ),
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser, label_points_default="none")
    return parser


def run(args: argparse.Namespace) -> int:
    parser = build_parser()
    validate_args(args, parser)
    from nucleosuite.parallel import run_region_per_contig
    if not getattr(args, "_per_contig_worker", False) and int(getattr(args, "cores", 1) or 1) > 1:
        return run_region_per_contig("dac", args, run)

    reporter = ProgressReporter("dac")
    selected_chromosomes = split_comma_values(args.chromosome)
    try:
        category_rules = parse_state_category_rules(args.category)
    except ValueError as exc:
        parser.error(str(exc))

    reporter.stage("Loading analysis regions")
    bigwig_files = expand_bigwig_inputs(args.bigwig)
    os.makedirs(args.output_dir, exist_ok=True)

    gene_selection_rows: List[Mapping[str, object]] = []
    try:
        if args.regions_bed:
            regions = read_regions_bed(
                path=args.regions_bed,
                state_column=args.state_column,
                strand_column=args.strand_column,
                state_name=args.state_name,
                region_mode=args.region_mode,
                extend=args.extend,
                strands=args.strands,
                selected_chromosomes=selected_chromosomes,
                min_region_length=args.min_region_length,
                category_rules=category_rules,
            )
        elif args.genes_bed:
            gene_queries = expand_gene_queries(args.gene, args.gene_list)
            regions, gene_selection_rows = read_selected_gene_regions(
                path=args.genes_bed,
                queries=gene_queries,
                gene_id_column=args.gene_id_column,
                gene_name_column=args.gene_name_column,
                mode=args.gene_dac_mode,
                pool_name=args.gene_pool_name,
                label_mode=args.gene_output_label,
                selected_chromosomes=selected_chromosomes,
                min_region_length=args.min_region_length,
            )
        else:
            regions = make_windows_from_chrom_sizes(
                chrom_sizes_path=args.chrom_sizes,
                scope=args.scope,
                selected_chromosomes=selected_chromosomes,
                window_size=args.window_size,
                state_name=args.state_name,
                min_region_length=args.min_region_length,
            )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    for chromosome in dict.fromkeys(region.chrom for region in regions):
        reporter.reading_contig("regions", chromosome)

    blacklist = load_blacklist_unbounded(args.blacklist_bed)
    args._blacklisted_anchor_exclusions = 0
    if blacklist is not None and (args.regions_bed or args.genes_bed):
        retained_regions = []
        for region in regions:
            anchor_start = region.anchor_start if region.anchor_start is not None else region.start
            anchor_end = region.anchor_end if region.anchor_end is not None else region.end
            if blacklist.overlaps(region.chrom, anchor_start, anchor_end):
                args._blacklisted_anchor_exclusions += 1
            else:
                retained_regions.append(region)
        regions = retained_regions
        if not regions:
            parser.error("No regions remained after blacklist filtering")

    regions_by_state = group_regions_by_state(regions)
    reporter.stage(
        f"Prepared {len(regions):,} regions in {len(regions_by_state):,} groups; "
        f"processing {len(bigwig_files):,} signal track(s)"
    )

    base_prefix = args.out_prefix or default_output_prefix(
        bigwig_files=bigwig_files,
        region_source_path=args.regions_bed or args.genes_bed,
        scope=args.scope,
        selected_chromosomes=selected_chromosomes,
    )

    if gene_selection_rows:
        gene_selection_path = os.path.join(
            args.output_dir,
            f"{sanitize_filename(base_prefix)}_selected_genes.tsv",
        )
        write_gene_selection_tsv(gene_selection_path, gene_selection_rows)
    else:
        gene_selection_path = None

    if not args.quiet:
        print(f"bigWig files: {len(bigwig_files):,}")
        print(f"Regions: {len(regions):,}")
        print(f"States/groups: {len(regions_by_state):,}")
        print(f"Maximum DAC distance: {args.dmax:,} bp")
        print(
            "Opportunity normalization: "
            + ("off" if args.no_normalize_dac else "on")
        )

    all_summary_rows: List[Mapping[str, object]] = []

    if args.separate_bigwigs:
        for bigwig_file in bigwig_files:
            track_prefix = sanitize_filename(
                f"{base_prefix}_{Path(bigwig_file).stem}"
            )
            all_summary_rows.extend(
                process_track_group(
                    bigwig_files=[bigwig_file],
                    regions_by_state=regions_by_state,
                    output_prefix=track_prefix,
                    output_dir=args.output_dir,
                    args=args,
                )
            )
    else:
        all_summary_rows.extend(
            process_track_group(
                bigwig_files=bigwig_files,
                regions_by_state=regions_by_state,
                output_prefix=base_prefix,
                output_dir=args.output_dir,
                args=args,
            )
        )

    reporter.stage("Writing DAC summary")
    summary_path = os.path.join(
        args.output_dir,
        f"{sanitize_filename(base_prefix)}_DAC_summary.tsv",
    )
    write_summary_tsv(summary_path, all_summary_rows)

    if not args.quiet:
        if gene_selection_path:
            print(f"Selected genes: {gene_selection_path}")
        print(f"Summary: {summary_path}")

    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
