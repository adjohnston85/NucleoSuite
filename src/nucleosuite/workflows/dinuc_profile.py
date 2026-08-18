"""Standalone observed dinucleotide-profile workflow."""

from __future__ import annotations

import os
from collections import Counter

import numpy as np

from nucleosuite.core.reference import extract_reference_sequence
from nucleosuite.io.summaries import write_fragment_outputs
from nucleosuite.sequence import dinucleotide
from nucleosuite.profile_plots import plot_dinucleotide_profile, plot_ww_ss_profile
from nucleosuite.workflows.common import (
    default_output_prefix,
    ensure_output_parent,
    prepare_fragment_run,
    input_paths_from_args,
    prepare_reference_if_needed,
    print_progress,
    randomize_fragments,
    set_random_seed,
)


def run(args) -> int:
    set_random_seed(args.seed)
    if args.out_prefix is None:
        args.out_prefix = default_output_prefix(input_paths_from_args(args), args.contigs)
    args.out_prefix = (
        f"{args.out_prefix}_dinuc_lower{args.frag_lower}_upper{args.frag_upper}"
    )
    ensure_output_parent(args.out_prefix)

    context = prepare_fragment_run(
        bam_paths=args.bamfiles,
        fragment_paths=args.fragment_files,
        fasta_path=args.fasta,
        chrom_sizes_path=args.chrom_sizes,
        contig_tokens=args.contigs,
        chunk_bp=args.chunk_bp,
        overlap_bp=args.overlap_bp,
        blacklist_path=getattr(args, "blacklist_bed", None),
    )
    fragment_length_set = set(range(args.frag_lower, args.frag_upper + 1))
    profile_positions = dinucleotide.expected_profile_positions(
        args.frag_lower, args.frag_upper
    )
    accumulator = dinucleotide.new_accumulator()
    length_counts = Counter()
    total_filtered = 0
    total_used = 0
    unique_bases = 0

    for suffix in (
        "_dinuc_profile.tsv",
        "_fragment_summary.tsv",
        "_fragment_length_counts.tsv",
    ):
        path = args.out_prefix + suffix
        if os.path.exists(path):
            os.remove(path)

    try:
        total_chunks = len(context.regions)
        previous_contig = None
        for chunk_index, region in enumerate(context.regions, start=1):
            print_progress(
                chunk_index, total_chunks, region, previous_contig=previous_contig
            )
            previous_contig = region.contig
            fragments = context.collect(
                contig=region.contig,
                start=region.adjusted_start,
                end=region.adjusted_end,
                max_duplicates=args.max_duplicates,
                subsample=args.subsample,
                dedup_scope=args.dedup_scope,
            )
            reference_context = prepare_reference_if_needed(
                fasta=context.fasta,
                contig=region.contig,
                start=region.adjusted_start,
                end=region.adjusted_end,
                required=True,
            )
            fragments = randomize_fragments(
                fragments=fragments,
                mode=args.randomize_mode,
                start=region.adjusted_start,
                end=region.adjusted_end,
                reference_context=reference_context,
                anchor_prob_start=args.anchor_prob_start,
                max_anchor_tries=args.max_anchor_tries,
                fallback=args.randomize_fallback,
            )

            owned = [
                fragment
                for fragment in fragments
                if region.original_start <= fragment[0] < region.original_end
            ]
            total_filtered += len(owned)
            covered = np.zeros(region.original_end - region.original_start, dtype=bool)
            for fragment_start, fragment_end in owned:
                fragment_length = fragment_end - fragment_start
                if fragment_length not in fragment_length_set:
                    continue
                total_used += 1
                length_counts[fragment_length] += 1
                sequence = extract_reference_sequence(
                    fasta=context.fasta,
                    context=reference_context,
                    seq_start=fragment_start,
                    seq_end=fragment_end,
                )
                dinucleotide.add_fragment(
                    accumulator=accumulator,
                    fragment_sequence=sequence,
                    fragment_start=fragment_start,
                    fragment_end=fragment_end,
                )
                overlap_start = max(fragment_start, region.original_start)
                overlap_end = min(fragment_end, region.original_end)
                if overlap_end > overlap_start:
                    covered[
                        overlap_start - region.original_start :
                        overlap_end - region.original_start
                    ] = True
            unique_bases += int(covered.sum())

        dinucleotide.write_profile(
            output_path=f"{args.out_prefix}_dinuc_profile.tsv",
            accumulator=accumulator,
            positions=profile_positions,
            fraction=args.dinuc_fraction,
        )
        plot_dinucleotide_profile(
            f"{args.out_prefix}_dinuc_profile.tsv",
            f"{args.out_prefix}_dinuc_profile.png",
            title="Dinucleotide profile",
        )
        plot_ww_ss_profile(
            f"{args.out_prefix}_dinuc_profile.tsv",
            f"{args.out_prefix}_ww_ss_profile.png",
            title="WW/SS dinucleotide profile",
        )
        write_fragment_outputs(
            output_prefix=args.out_prefix,
            total_fragments_filtered=total_filtered,
            total_fragments_used=total_used,
            unique_bases_covered=unique_bases,
            length_counts=length_counts,
            dedup_scope=args.dedup_scope,
            max_duplicates=args.max_duplicates,
            max_per_coordinate=getattr(args, "max_per_coordinate", 0),
        )
        print(
            f"[INFO] Dinucleotide fragments used: {accumulator['fragments_used']:,}; "
            f"skipped: {accumulator['fragments_skipped']:,}."
        )
    finally:
        context.close()
    return 0
