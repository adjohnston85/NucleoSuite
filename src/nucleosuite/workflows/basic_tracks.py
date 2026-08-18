"""Standalone fragment-coordinate-to-track workflows for coverage, dyads and fragment ends."""

from __future__ import annotations

import os
from collections import Counter

import numpy as np

from nucleosuite.io.summaries import write_fragment_outputs
from nucleosuite.io.tracks import (
    close_track_handles,
    open_track_handles,
    remove_stale_track_outputs,
    write_tracks,
)
from nucleosuite.scoring import basic_tracks
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


def run(args, command_label: str, tracks: list[str]) -> int:
    set_random_seed(args.seed)
    if args.out_prefix is None:
        args.out_prefix = default_output_prefix(input_paths_from_args(args), args.contigs)
    args.out_prefix = (
        f"{args.out_prefix}_{command_label}_lower{args.frag_lower}"
        f"_upper{args.frag_upper}"
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
    remove_stale_track_outputs(args.out_prefix, tracks)
    for suffix in ("_fragment_summary.tsv", "_fragment_length_counts.tsv"):
        path = args.out_prefix + suffix
        if os.path.exists(path):
            os.remove(path)

    bigwig_handles, wig_handles = open_track_handles(
        output_prefix=args.out_prefix,
        tracks=tracks,
        output_format=args.output_format,
        bigwig_header=context.bigwig_header,
    )

    total_filtered = 0
    total_used = 0
    unique_bases = 0
    length_counts = Counter()
    reference_required = args.randomize_mode == "dinuc_anchor"

    completed = False
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
                required=reference_required,
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

            arrays = basic_tracks.new_arrays(
                region.adjusted_end - region.adjusted_start
            )
            for fragment_start, fragment_end in fragments:
                if fragment_end - fragment_start not in fragment_length_set:
                    continue
                basic_tracks.add_fragment(
                    arrays=arrays,
                    fragment_start=fragment_start,
                    fragment_end=fragment_end,
                    window_start=region.adjusted_start,
                    window_end=region.adjusted_end,
                    even_dyad=getattr(args, "even_dyad", "split"),
                )

            basic_tracks.cap_sparse_arrays(arrays, args.max_per_coordinate)
            scores = basic_tracks.to_scores(
                arrays, region.contig, region.adjusted_start
            )
            write_tracks(
                scores=scores,
                contig=region.contig,
                adjusted_start=region.adjusted_start,
                original_start=region.original_start,
                original_end=region.original_end,
                tracks=tracks,
                bigwig_handles=bigwig_handles,
                wig_handles=wig_handles,
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
                overlap_start = max(fragment_start, region.original_start)
                overlap_end = min(fragment_end, region.original_end)
                if overlap_end > overlap_start:
                    covered[
                        overlap_start - region.original_start :
                        overlap_end - region.original_start
                    ] = True
            unique_bases += int(covered.sum())

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
        completed = True
    finally:
        close_track_handles(bigwig_handles, wig_handles, commit=completed)
        context.close()

    return 0
