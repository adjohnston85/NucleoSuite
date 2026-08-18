"""Full WPS workflow used by ``nucleosuite wps``."""

from __future__ import annotations

import os
from collections import Counter

import numpy as np

from nucleosuite.io.summaries import write_fragment_outputs
from nucleosuite.io.intervals import finalise_interval_files
from nucleosuite.io.tracks import (
    close_track_handles,
    open_track_handles,
    remove_stale_track_outputs,
    write_tracks,
)
from nucleosuite.peaks.common import write_bed8_records
from nucleosuite.peaks.pns import call_records as call_pns_records
from nucleosuite.peaks.pns import write_records as write_pns_records
from nucleosuite.peaks.wps import call_records as call_wps_records
from nucleosuite.scoring import basic_tracks
from nucleosuite.scoring import wps as wps_scoring
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
        f"{args.out_prefix}_prot{args.protection}_lower{args.frag_lower}"
        f"_upper{args.frag_upper}"
    )
    if args.randomize_mode != "none":
        args.out_prefix += f"_rand{args.randomize_mode}"
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
    fragment_lengths = range(args.frag_lower, args.frag_upper + 1)
    fragment_length_set = set(fragment_lengths)
    distributions = wps_scoring.precompute_distributions(
        fragment_lengths, args.protection
    )

    all_tracks = set(wps_scoring.WPS_TRACKS) | set(basic_tracks.BASIC_TRACKS)
    remove_stale_track_outputs(args.out_prefix, all_tracks)
    for suffix in (
        "_nucleosome_regions.bed",
        "_breakpoint_peaks.bed",
        "_nucleosome_regions.bb",
        "_breakpoint_peaks.bb",
        "_fragment_summary.tsv",
        "_fragment_length_counts.tsv",
    ):
        path = args.out_prefix + suffix
        if os.path.exists(path):
            os.remove(path)

    bigwig_handles, wig_handles = open_track_handles(
        output_prefix=args.out_prefix,
        tracks=args.score_tracks,
        output_format=args.score_format,
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

            reference_length = region.adjusted_end - region.adjusted_start
            basic_arrays = basic_tracks.new_arrays(reference_length)
            wps_array = wps_scoring.new_array(reference_length)
            for fragment_start, fragment_end in fragments:
                if fragment_end - fragment_start in fragment_length_set:
                    basic_tracks.add_fragment(
                        arrays=basic_arrays,
                        fragment_start=fragment_start,
                        fragment_end=fragment_end,
                        window_start=region.adjusted_start,
                        window_end=region.adjusted_end,
                        even_dyad=args.even_dyad,
                    )
                    wps_scoring.add_fragment(
                        wps_array=wps_array,
                        fragment_start=fragment_start,
                        fragment_end=fragment_end,
                        window_start=region.adjusted_start,
                        protection=args.protection,
                        distributions=distributions,
                    )

            basic_tracks.cap_sparse_arrays(basic_arrays, args.max_per_coordinate)
            scores = basic_tracks.to_scores(
                basic_arrays, region.contig, region.adjusted_start
            )
            scores.update(
                wps_scoring.to_scores(
                    wps_array=wps_array,
                    contig=region.contig,
                    start=region.adjusted_start,
                    baseline_window=args.baseline_window,
                    smooth_window=args.sg_window,
                    smooth_order=args.sg_order,
                )
            )

            owned_fragments = [
                fragment
                for fragment in fragments
                if region.original_start <= fragment[0] < region.original_end
            ]
            total_filtered += len(owned_fragments)
            covered = np.zeros(region.original_end - region.original_start, dtype=bool)
            for fragment_start, fragment_end in owned_fragments:
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

            selected_track = scores[args.peak_track][0][2]
            mode = "w" if chunk_index == 1 else "a"
            if args.peak_caller == "wps":
                nucleosome_records = call_wps_records(
                    scores=selected_track,
                    chrom=region.contig,
                    adjusted_start=region.adjusted_start,
                    core_start=region.original_start,
                    core_end=region.original_end,
                    merge_gap_bp=args.peak_merge_gap,
                    min_length=args.peak_minlen,
                    max_length=args.peak_maxlen,
                    max_region_length=args.peak_maxregion,
                    score_cutoff=args.peak_varicutoff,
                    flip_scores=False,
                )
                breakpoint_records = call_wps_records(
                    scores=-1.0 * selected_track,
                    chrom=region.contig,
                    adjusted_start=region.adjusted_start,
                    core_start=region.original_start,
                    core_end=region.original_end,
                    merge_gap_bp=args.peak_merge_gap,
                    min_length=args.peak_minlen,
                    max_length=args.peak_maxlen,
                    max_region_length=args.peak_maxregion,
                    score_cutoff=args.peak_varicutoff,
                    flip_scores=True,
                )
                write_bed8_records(
                    f"{args.out_prefix}_nucleosome_regions.bed",
                    nucleosome_records,
                    "nuc",
                    args.peak_score_scale,
                    mode,
                )
                write_bed8_records(
                    f"{args.out_prefix}_breakpoint_peaks.bed",
                    breakpoint_records,
                    "brk",
                    args.peak_score_scale,
                    mode,
                )
            elif args.peak_caller == "pns":
                coverage = scores["coverage"][0][2]
                nucleosome_records = call_pns_records(
                    scores=selected_track,
                    chrom=region.contig,
                    adjusted_start=region.adjusted_start,
                    core_start=region.original_start,
                    core_end=region.original_end,
                    min_length=args.peak_minlen,
                    max_nonpositive_run=args.peak_merge_gap,
                    coverage_scores=coverage,
                )
                breakpoint_records = call_pns_records(
                    scores=-1.0 * selected_track,
                    chrom=region.contig,
                    adjusted_start=region.adjusted_start,
                    core_start=region.original_start,
                    core_end=region.original_end,
                    min_length=args.peak_minlen,
                    max_nonpositive_run=args.peak_merge_gap,
                    flip_scores=True,
                    coverage_scores=coverage,
                )
                write_pns_records(
                    f"{args.out_prefix}_nucleosome_regions.bed",
                    nucleosome_records,
                    "nuc",
                    args.peak_score_scale,
                    mode,
                )
                write_pns_records(
                    f"{args.out_prefix}_breakpoint_peaks.bed",
                    breakpoint_records,
                    "brk",
                    args.peak_score_scale,
                    mode,
                )

            write_tracks(
                scores=scores,
                contig=region.contig,
                adjusted_start=region.adjusted_start,
                original_start=region.original_start,
                original_end=region.original_end,
                tracks=args.score_tracks,
                bigwig_handles=bigwig_handles,
                wig_handles=wig_handles,
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
        if args.peak_caller != "none":
            finalise_interval_files(
                [
                    f"{args.out_prefix}_nucleosome_regions.bed",
                    f"{args.out_prefix}_breakpoint_peaks.bed",
                ],
                args.interval_format,
                context.bigwig_header,
            )
        completed = True
    finally:
        close_track_handles(bigwig_handles, wig_handles, commit=completed)
        context.close()

    return 0
