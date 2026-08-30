"""Bootstrap-stabilized fragment-mode estimation for fragment collections."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from nucleosuite.core.fragments import is_softclipped_or_padded
from nucleosuite.core.fragment_inputs import BamFragmentSource


@dataclass(frozen=True)
class ModeEstimate:
    mode: int
    ci_low: float
    ci_high: float
    sampled_fragments: int
    mode_search_fragments: int
    converged: bool
    checkpoints: int
    histogram: tuple[int, ...]
    search_lower: int
    search_upper: int


def mode_estimate_message(label: str, estimate: ModeEstimate) -> str:
    """Return a concise user-facing summary of a resolved estimate."""

    return (
        f"{label}: {estimate.mode} bp "
        f"(bootstrap 95% CI {estimate.ci_low:.6g}-{estimate.ci_high:.6g} bp; "
        f"{estimate.sampled_fragments:,} accepted fragments sampled; "
        f"{estimate.mode_search_fragments:,} within "
        f"{estimate.search_lower}-{estimate.search_upper} bp; "
        f"converged={str(estimate.converged).lower()})"
    )


def write_single_mode_report(
    path: str | Path,
    *,
    estimate: ModeEstimate | None,
    resolved_mode: int,
    mode_source: str,
    seed: int,
) -> Path:
    """Write the resolved mode and the calculation settings used to obtain it."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"mode_source\t{mode_source}\n")
        handle.write(f"resolved_mode_bp\t{resolved_mode}\n")
        handle.write(f"seed\t{seed}\n")
        if estimate is None:
            return output
        handle.write(f"bootstrap_ci_low_bp\t{estimate.ci_low:.12g}\n")
        handle.write(f"bootstrap_ci_high_bp\t{estimate.ci_high:.12g}\n")
        handle.write(f"sampled_fragments\t{estimate.sampled_fragments}\n")
        handle.write(
            f"mode_search_fragments\t{estimate.mode_search_fragments}\n"
        )
        handle.write(f"mode_search_lower_bp\t{estimate.search_lower}\n")
        handle.write(f"mode_search_upper_bp\t{estimate.search_upper}\n")
        handle.write(f"converged\t{str(estimate.converged).lower()}\n")
        handle.write(f"checkpoints\t{estimate.checkpoints}\n")
        histogram = ";".join(
            f"{estimate.search_lower + index}:{count}"
            for index, count in enumerate(estimate.histogram)
        )
        handle.write(f"histogram_length_bp_count\t{histogram}\n")
    return output


def _point_mode(counts: np.ndarray, lower: int) -> int:
    """Return the lowest integer length with the largest observed count."""
    if counts.size == 0 or int(np.sum(counts)) <= 0:
        raise ValueError("No fragments occurred in the requested mode-search range")
    # argmax selects the first maximum, giving stable lower-length tie breaking.
    return int(lower + int(np.argmax(counts)))


def bootstrap_histogram_mode(
    counts: Sequence[int] | np.ndarray,
    *,
    lower: int,
    replicates: int = 200,
    seed: int = 12345,
) -> tuple[int, float, float, np.ndarray]:
    """Estimate an integer mode and percentile bootstrap interval."""

    values = np.asarray(counts, dtype=np.int64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("counts must be a non-empty one-dimensional histogram")
    if np.any(values < 0):
        raise ValueError("Histogram counts cannot be negative")
    total = int(np.sum(values))
    if total < 1:
        raise ValueError("No fragments occurred in the requested mode-search range")
    if replicates < 1:
        raise ValueError("Bootstrap replicates must be at least 1")

    point = _point_mode(values, lower)
    probabilities = values.astype(np.float64) / float(total)
    rng = np.random.default_rng(seed)
    modes = np.empty(replicates, dtype=np.int64)
    for index in range(replicates):
        draw = rng.multinomial(total, probabilities)
        modes[index] = _point_mode(draw, lower)
    ci_low, ci_high = np.percentile(modes, [2.5, 97.5])
    return point, float(ci_low), float(ci_high), modes


def estimate_mode_from_lengths(
    lengths: Sequence[int] | np.ndarray,
    *,
    search_lower: int = 120,
    search_upper: int = 250,
    bootstrap_replicates: int = 200,
    seed: int = 12345,
) -> ModeEstimate:
    """Estimate one mode from an already sampled fragment-length collection."""

    if search_lower < 1 or search_upper < search_lower:
        raise ValueError("Require 1 <= mode-search-lower <= mode-search-upper")
    array = np.asarray(lengths, dtype=np.int64)
    selected = array[(array >= search_lower) & (array <= search_upper)]
    counts = np.bincount(
        selected - search_lower,
        minlength=search_upper - search_lower + 1,
    )
    mode, low, high, _modes = bootstrap_histogram_mode(
        counts,
        lower=search_lower,
        replicates=bootstrap_replicates,
        seed=seed,
    )
    return ModeEstimate(
        mode=mode,
        ci_low=low,
        ci_high=high,
        sampled_fragments=int(array.size),
        mode_search_fragments=int(selected.size),
        converged=False,
        checkpoints=1,
        histogram=tuple(int(value) for value in counts),
        search_lower=search_lower,
        search_upper=search_upper,
    )


def _bam_blocks(references: Sequence[str], lengths: Sequence[int], block_bp: int):
    blocks: list[tuple[str, int, int]] = []
    for chrom, length in zip(references, lengths):
        for start in range(0, int(length), block_bp):
            blocks.append((str(chrom), start, min(int(length), start + block_bp)))
    return blocks


def estimate_bam_fragment_mode(
    bam_paths: Sequence[str | Path],
    *,
    frag_lower: int = 120,
    frag_upper: int = 500,
    search_lower: int = 120,
    search_upper: int = 250,
    minimum_fragments: int = 100_000,
    batch_fragments: int = 25_000,
    maximum_fragments: int = 1_000_000,
    stable_checkpoints: int = 3,
    maximum_mode_change: int = 1,
    maximum_ci_width: float = 4.0,
    bootstrap_replicates: int = 200,
    block_bp: int = 1_000_000,
    seed: int = 12345,
    blacklist_bed: str | Path | None = None,
    max_duplicates: int = 1,
    dedup_scope: str = "all_bams",
    contig_tokens: Sequence[str] | None = None,
) -> ModeEstimate:
    """Randomly visit indexed genomic blocks until the bootstrap mode stabilizes."""

    if not bam_paths:
        raise ValueError("At least one BAM is required for automatic mode estimation")
    if frag_lower < 1 or frag_upper < frag_lower:
        raise ValueError("Invalid fragment-length range")
    if search_lower < frag_lower or search_upper > frag_upper or search_upper < search_lower:
        raise ValueError("Mode-search range must lie within the accepted fragment range")
    if min(minimum_fragments, batch_fragments, maximum_fragments, stable_checkpoints, block_bp) < 1:
        raise ValueError("Mode-sampling counts and block size must be positive")
    if maximum_fragments < minimum_fragments:
        raise ValueError("maximum_fragments must be at least minimum_fragments")
    if max_duplicates < 0:
        raise ValueError("max_duplicates must be non-negative")
    if dedup_scope not in {"all_bams", "per_bam"}:
        raise ValueError("dedup_scope must be 'all_bams' or 'per_bam'")

    paths = [Path(path) for path in bam_paths]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    source = BamFragmentSource([str(path) for path in paths])
    try:
        references = tuple(source.references)
        lengths = tuple(int(value) for value in source.lengths)
        handles = source.handles

        blacklist = None
        if blacklist_bed is not None:
            from nucleosuite.core.blacklist import load_blacklist
            blacklist = load_blacklist(blacklist_bed, references, lengths)

        if contig_tokens:
            from nucleosuite.core.regions import (
                build_processing_regions,
                expand_contig_tokens,
            )

            selected_specs = expand_contig_tokens(contig_tokens, references)
            regions, _selected_names = build_processing_regions(
                selected_specs=selected_specs,
                references=references,
                lengths=lengths,
                chunk_bp=block_bp,
                overlap_bp=0,
            )
            blocks = [
                (region.contig, region.original_start, region.original_end)
                for region in regions
            ]
        else:
            blocks = _bam_blocks(references, lengths, block_bp)
        random.Random(seed).shuffle(blocks)
        sampled: list[int] = []
        search_counts = np.zeros(search_upper - search_lower + 1, dtype=np.int64)
        coordinate_counts: dict[tuple[object, ...], int] = {}
        checkpoint_modes: list[int] = []
        last_checkpoint = 0
        latest_ci = (math.nan, math.nan)
        converged = False

        for chrom, block_start, block_end in blocks:
            for bam_index, handle in enumerate(handles):
                source_chrom = source.source_contigs[bam_index].get(chrom)
                if source_chrom is None:
                    continue
                for read in handle.fetch(source_chrom, block_start, block_end):
                    if (
                        read.is_unmapped
                        or read.mate_is_unmapped
                        or read.is_secondary
                        or read.is_supplementary
                        or read.is_qcfail
                        or read.is_duplicate
                        or not read.is_paired
                        or read.reference_id != read.next_reference_id
                        or read.is_reverse == read.mate_is_reverse
                        or int(read.template_length) <= 0
                    ):
                        continue
                    if is_softclipped_or_padded(read.cigartuples):
                        continue
                    start = int(read.reference_start)
                    if not block_start <= start < block_end:
                        continue
                    end = start + int(read.template_length)
                    length = end - start
                    if length < frag_lower or length > frag_upper:
                        continue
                    if blacklist is not None and blacklist.overlaps(chrom, start, end):
                        continue
                    key = (
                        (chrom, start, end)
                        if dedup_scope == "all_bams"
                        else (bam_index, chrom, start, end)
                    )
                    if max_duplicates > 0:
                        count = coordinate_counts.get(key, 0)
                        if count >= max_duplicates:
                            continue
                        coordinate_counts[key] = count + 1
                    sampled.append(length)
                    if search_lower <= length <= search_upper:
                        search_counts[length - search_lower] += 1

                    if len(sampled) >= maximum_fragments:
                        break
                    if (
                        len(sampled) >= minimum_fragments
                        and len(sampled) - last_checkpoint >= batch_fragments
                        and int(np.sum(search_counts)) > 0
                    ):
                        point, low, high, _ = bootstrap_histogram_mode(
                            search_counts,
                            lower=search_lower,
                            replicates=bootstrap_replicates,
                            seed=seed + len(checkpoint_modes),
                        )
                        checkpoint_modes.append(point)
                        latest_ci = (low, high)
                        last_checkpoint = len(sampled)
                        if len(checkpoint_modes) >= stable_checkpoints:
                            recent = checkpoint_modes[-stable_checkpoints:]
                            converged = (
                                max(recent) - min(recent) <= maximum_mode_change
                                and high - low <= maximum_ci_width
                            )
                            if converged:
                                break
                if converged or len(sampled) >= maximum_fragments:
                    break
            if converged or len(sampled) >= maximum_fragments:
                break

        if int(np.sum(search_counts)) == 0:
            raise ValueError("No accepted fragments occurred in the mode-search range")
        final_mode, low, high, _ = bootstrap_histogram_mode(
            search_counts,
            lower=search_lower,
            replicates=bootstrap_replicates,
            seed=seed + 100_000,
        )
        if not checkpoint_modes:
            checkpoint_modes.append(final_mode)
        latest_ci = (low, high)
        return ModeEstimate(
            mode=final_mode,
            ci_low=latest_ci[0],
            ci_high=latest_ci[1],
            sampled_fragments=len(sampled),
            mode_search_fragments=int(np.sum(search_counts)),
            converged=converged,
            checkpoints=len(checkpoint_modes),
            histogram=tuple(int(value) for value in search_counts),
            search_lower=search_lower,
            search_upper=search_upper,
        )
    finally:
        source.close()


def estimate_fragment_file_mode(
    fragment_paths: Sequence[str | Path],
    *,
    frag_lower: int,
    frag_upper: int,
    search_lower: int,
    search_upper: int,
    minimum_fragments: int = 100_000,
    batch_fragments: int = 25_000,
    maximum_fragments: int = 1_000_000,
    stable_checkpoints: int = 3,
    maximum_mode_change: int = 1,
    maximum_ci_width: float = 4.0,
    bootstrap_replicates: int = 200,
    block_bp: int = 1_000_000,
    seed: int = 12345,
    blacklist_bed: str | Path | None = None,
    max_duplicates: int = 1,
    dedup_scope: str = "all_bams",
    contig_tokens: Sequence[str] | None = None,
    chrom_sizes: str | Path | None = None,
    fasta_path: str | Path | None = None,
) -> ModeEstimate:
    """Estimate a mode from BED, BED.gz, or bigBed fragment coordinates.

    Indexed genomic blocks are shuffled just as they are for BAM estimation.
    A fragment is owned by the block containing its start, which prevents a
    boundary-spanning interval from entering the histogram more than once.
    """

    if not fragment_paths:
        raise ValueError("At least one fragment interval file is required")
    if frag_lower < 1 or frag_upper < frag_lower:
        raise ValueError("Invalid fragment-length range")
    if search_lower < frag_lower or search_upper > frag_upper or search_upper < search_lower:
        raise ValueError("Mode-search range must lie within the accepted fragment range")
    if min(minimum_fragments, batch_fragments, maximum_fragments, stable_checkpoints, block_bp) < 1:
        raise ValueError("Mode-sampling counts and block size must be positive")
    if maximum_fragments < minimum_fragments:
        raise ValueError("maximum_fragments must be at least minimum_fragments")

    from nucleosuite.core.fragment_inputs import open_fragment_source
    from nucleosuite.core.regions import build_processing_regions, expand_contig_tokens

    fasta = None
    source = None
    try:
        if fasta_path is not None:
            import pysam

            fasta = pysam.FastaFile(str(fasta_path))
        source = open_fragment_source(
            fragment_paths=[str(Path(path)) for path in fragment_paths],
            chrom_sizes=str(chrom_sizes) if chrom_sizes is not None else None,
            fasta=fasta,
            blacklist_path=str(blacklist_bed) if blacklist_bed is not None else None,
        )
        selected_specs = expand_contig_tokens(contig_tokens, source.references)
        regions, _selected_names = build_processing_regions(
            selected_specs=selected_specs,
            references=source.references,
            lengths=source.lengths,
            chunk_bp=block_bp,
            overlap_bp=0,
        )
        blocks = [
            (region.contig, region.original_start, region.original_end)
            for region in regions
        ]
        random.Random(seed).shuffle(blocks)

        sampled_fragments = 0
        search_counts = np.zeros(search_upper - search_lower + 1, dtype=np.int64)
        checkpoint_modes: list[int] = []
        last_checkpoint = 0
        converged = False

        for chrom, block_start, block_end in blocks:
            fragments = source.fetch(
                chrom,
                block_start,
                block_end,
                max_per_coordinate=max_duplicates,
                subsample=None,
                dedup_scope=dedup_scope,
            )
            for start, end in fragments:
                if not block_start <= int(start) < block_end:
                    continue
                length = int(end) - int(start)
                if length < frag_lower or length > frag_upper:
                    continue
                sampled_fragments += 1
                if search_lower <= length <= search_upper:
                    search_counts[length - search_lower] += 1
                if sampled_fragments >= maximum_fragments:
                    break
                if (
                    sampled_fragments >= minimum_fragments
                    and sampled_fragments - last_checkpoint >= batch_fragments
                    and int(np.sum(search_counts)) > 0
                ):
                    point, low, high, _ = bootstrap_histogram_mode(
                        search_counts,
                        lower=search_lower,
                        replicates=bootstrap_replicates,
                        seed=seed + len(checkpoint_modes),
                    )
                    checkpoint_modes.append(point)
                    last_checkpoint = sampled_fragments
                    if len(checkpoint_modes) >= stable_checkpoints:
                        recent = checkpoint_modes[-stable_checkpoints:]
                        converged = (
                            max(recent) - min(recent) <= maximum_mode_change
                            and high - low <= maximum_ci_width
                        )
                        if converged:
                            break
            if converged or sampled_fragments >= maximum_fragments:
                break

        if int(np.sum(search_counts)) == 0:
            raise ValueError("No accepted fragments occurred in the mode-search range")
        mode, low, high, _ = bootstrap_histogram_mode(
            search_counts,
            lower=search_lower,
            replicates=bootstrap_replicates,
            seed=seed + 100_000,
        )
        if not checkpoint_modes:
            checkpoint_modes.append(mode)
        return ModeEstimate(
            mode=mode,
            ci_low=low,
            ci_high=high,
            sampled_fragments=sampled_fragments,
            mode_search_fragments=int(np.sum(search_counts)),
            converged=converged,
            checkpoints=len(checkpoint_modes),
            histogram=tuple(int(value) for value in search_counts),
            search_lower=search_lower,
            search_upper=search_upper,
        )
    finally:
        if source is not None:
            source.close()
        if fasta is not None:
            fasta.close()


def pooled_mode_estimate(
    target: ModeEstimate,
    control: ModeEstimate,
    *,
    bootstrap_replicates: int = 200,
    seed: int = 12345,
) -> ModeEstimate:
    """Return an equal-sample-weighted pooled mode from two mode histograms."""

    if (target.search_lower, target.search_upper) != (
        control.search_lower,
        control.search_upper,
    ):
        raise ValueError("Target and control mode-search ranges differ")
    left = np.asarray(target.histogram, dtype=np.float64)
    right = np.asarray(control.histogram, dtype=np.float64)
    if np.sum(left) <= 0 or np.sum(right) <= 0:
        raise ValueError("Target and control histograms must both be non-empty")
    probabilities = 0.5 * (left / np.sum(left) + right / np.sum(right))
    effective_n = max(1, min(target.mode_search_fragments, control.mode_search_fragments))
    pooled_counts = np.rint(probabilities * effective_n).astype(np.int64)
    if int(np.sum(pooled_counts)) == 0:
        pooled_counts[int(np.argmax(probabilities))] = 1
    mode, low, high, _ = bootstrap_histogram_mode(
        pooled_counts,
        lower=target.search_lower,
        replicates=bootstrap_replicates,
        seed=seed,
    )
    return ModeEstimate(
        mode=mode,
        ci_low=low,
        ci_high=high,
        sampled_fragments=target.sampled_fragments + control.sampled_fragments,
        mode_search_fragments=int(np.sum(pooled_counts)),
        converged=target.converged and control.converged,
        checkpoints=min(target.checkpoints, control.checkpoints),
        histogram=tuple(int(value) for value in pooled_counts),
        search_lower=target.search_lower,
        search_upper=target.search_upper,
    )
