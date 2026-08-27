"""Nucleosome-scoring workflow used by ``nucleosuite nuc-score``."""

from __future__ import annotations

import os
from collections import Counter

import numpy as np

from nucleosuite.core.reference import extract_reference_sequence
from nucleosuite.io.summaries import write_fragment_outputs
from nucleosuite.io.intervals import finalise_interval_files
from nucleosuite.io.tracks import (
    TrackHandleMap,
    close_track_handles,
    open_track_handles,
    remove_stale_track_outputs,
    write_tracks,
)
from nucleosuite.peaks.pns import (
    call_records,
    filter_records_by_coverage,
    write_records,
)
from nucleosuite.scoring import basic_tracks
from nucleosuite.scoring import pns as pns_scoring
from nucleosuite.sequence import dinucleotide
from nucleosuite.sequence.ww_types import (
    ALL_OUTPUT_GROUPS,
    WW_TYPE_GROUPS,
    classify_fragment,
    group_output_prefix,
    write_summary as write_ww_summary,
)
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


def _remove_stale_outputs(output_prefix: str) -> None:
    all_tracks = set(pns_scoring.PNS_TRACKS) | set(pns_scoring.SNS_TRACKS) | set(pns_scoring.BNS_TRACKS) | set(pns_scoring.TNS_TRACKS) | set(basic_tracks.BASIC_TRACKS)
    for group in ALL_OUTPUT_GROUPS:
        prefix = group_output_prefix(output_prefix, group)
        remove_stale_track_outputs(prefix, all_tracks)
        for suffix in (
            "_fragment_summary.tsv",
            "_fragment_length_counts.tsv",
            "_nucleosome_regions.bed",
            "_breakpoint_peaks.bed",
            "_nucleosome_regions.bb",
            "_breakpoint_peaks.bb",
            "_dinuc_profile.tsv",
        ):
            path = prefix + suffix
            if os.path.exists(path):
                os.remove(path)
    summary = f"{output_prefix}_ww_type_summary.tsv"
    if os.path.exists(summary):
        os.remove(summary)


def run(args) -> int:
    set_random_seed(args.seed)
    if args.out_prefix is None:
        args.out_prefix = default_output_prefix(input_paths_from_args(args), args.contigs)
    args.out_prefix = (
        f"{args.out_prefix}_method{args.scoring_method}_mode{args.mode_length}_lower{args.frag_lower}"
        f"_upper{args.frag_upper}_smooth{args.smooth_window}x{args.smooth_order}"
    )
    ensure_output_parent(args.out_prefix)

    reference_required = (
        args.randomize_mode == "dinuc_anchor"
        or args.dinuc_profile
        or args.split_ww_types
    )
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
    if args.score_mode == "on":
        centred_distributions, positive_distributions = (
            pns_scoring.precompute_distributions(
                fragment_lengths, args.mode_length, args.scoring_method
            )
        )
    else:
        centred_distributions, positive_distributions = {}, {}

    output_groups = ALL_OUTPUT_GROUPS if args.split_ww_types else ("all",)
    _remove_stale_outputs(args.out_prefix)

    stats = {
        group: {
            "total_filtered": 0,
            "total_used": 0,
            "unique_bases": 0,
            "length_counts": Counter(),
        }
        for group in output_groups
    }
    type_counts = Counter()
    peak_filter_stats = {
        group: {"called": 0, "retained": 0}
        for group in output_groups
    }

    profile_positions = dinucleotide.expected_profile_positions(
        args.frag_lower, args.frag_upper
    )
    dinuc_accumulators = (
        {group: dinucleotide.new_accumulator() for group in output_groups}
        if args.dinuc_profile
        else {}
    )

    bigwig_handles = {group: TrackHandleMap() for group in output_groups}
    wig_handles = {group: TrackHandleMap() for group in output_groups}
    completed = False
    try:
        for group in output_groups:
            prefix = group_output_prefix(args.out_prefix, group)
            pns_bigwig, pns_wig = open_track_handles(
                output_prefix=prefix,
                tracks=args.score_tracks,
                output_format=args.score_format,
                bigwig_header=context.bigwig_header,
            )
            other_bigwig, other_wig = open_track_handles(
                output_prefix=prefix,
                tracks=args.other_tracks,
                output_format=args.other_format,
                bigwig_header=context.bigwig_header,
            )
            bigwig_handles[group].update(pns_bigwig)
            bigwig_handles[group].update(other_bigwig)
            bigwig_handles[group].paths.update(pns_bigwig.paths)
            bigwig_handles[group].paths.update(other_bigwig.paths)
            wig_handles[group].update(pns_wig)
            wig_handles[group].update(other_wig)
            wig_handles[group].paths.update(pns_wig.paths)
            wig_handles[group].paths.update(other_wig.paths)

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
            basic_by_group = {
                group: basic_tracks.new_arrays(reference_length)
                for group in output_groups
            }
            pns_by_group = (
                {
                    group: pns_scoring.new_arrays(reference_length, args.scoring_method)
                    for group in output_groups
                }
                if args.score_mode == "on"
                else {}
            )

            fragment_records = []
            for fragment_start, fragment_end in fragments:
                fragment_length = fragment_end - fragment_start
                ww_type = None
                if args.split_ww_types and fragment_length in fragment_length_set:
                    ww_type = classify_fragment(
                        fasta=context.fasta,
                        reference_context=reference_context,
                        fragment_start=fragment_start,
                        fragment_end=fragment_end,
                    )
                fragment_records.append((fragment_start, fragment_end, ww_type))

                if fragment_length not in fragment_length_set:
                    continue
                target_groups = ["all"]
                if args.split_ww_types and ww_type in WW_TYPE_GROUPS:
                    target_groups.append(ww_type)

                for group in target_groups:
                    basic_tracks.add_fragment(
                        arrays=basic_by_group[group],
                        fragment_start=fragment_start,
                        fragment_end=fragment_end,
                        window_start=region.adjusted_start,
                        window_end=region.adjusted_end,
                        even_dyad=args.even_dyad,
                    )
                    if args.score_mode == "on":
                        pns_scoring.add_fragment(
                            arrays=pns_by_group[group],
                            fragment_start=fragment_start,
                            fragment_end=fragment_end,
                            window_start=region.adjusted_start,
                            window_end=region.adjusted_end,
                            mode_dna_length=args.mode_length,
                            centred_distributions=centred_distributions,
                            positive_distributions=positive_distributions,
                            scoring_method=args.scoring_method,
                        )

            for arrays in basic_by_group.values():
                basic_tracks.cap_sparse_arrays(arrays, args.max_per_coordinate)

            scores_by_group = {}
            for group in output_groups:
                scores = basic_tracks.to_scores(
                    basic_by_group[group], region.contig, region.adjusted_start
                )
                if args.score_mode == "on" and args.peak_calling:
                    scores.update(
                        pns_scoring.to_scores(
                            pns_by_group[group],
                            region.contig,
                            region.adjusted_start,
                            smooth_window=args.smooth_window,
                            smooth_order=args.smooth_order,
                            scoring_method=args.scoring_method,
                        )
                    )
                scores_by_group[group] = scores

            owned_records = [
                record
                for record in fragment_records
                if region.original_start <= record[0] < region.original_end
            ]
            stats["all"]["total_filtered"] += len(owned_records)
            covered_by_group = {
                group: np.zeros(
                    region.original_end - region.original_start, dtype=bool
                )
                for group in output_groups
            }

            for fragment_start, fragment_end, ww_type in owned_records:
                fragment_length = fragment_end - fragment_start
                if fragment_length not in fragment_length_set:
                    continue

                stats["all"]["total_used"] += 1
                stats["all"]["length_counts"][fragment_length] += 1
                targets = ["all"]
                if args.split_ww_types:
                    if ww_type in WW_TYPE_GROUPS:
                        type_counts[ww_type] += 1
                        targets.append(ww_type)
                        stats[ww_type]["total_filtered"] += 1
                        stats[ww_type]["total_used"] += 1
                        stats[ww_type]["length_counts"][fragment_length] += 1
                    else:
                        type_counts["unclassified"] += 1

                overlap_start = max(fragment_start, region.original_start)
                overlap_end = min(fragment_end, region.original_end)
                if overlap_end > overlap_start:
                    for group in targets:
                        covered_by_group[group][
                            overlap_start - region.original_start :
                            overlap_end - region.original_start
                        ] = True

                if args.dinuc_profile:
                    fragment_sequence = extract_reference_sequence(
                        fasta=context.fasta,
                        context=reference_context,
                        seq_start=fragment_start,
                        seq_end=fragment_end,
                    )
                    for group in targets:
                        dinucleotide.add_fragment(
                            accumulator=dinuc_accumulators[group],
                            fragment_sequence=fragment_sequence,
                            fragment_start=fragment_start,
                            fragment_end=fragment_end,
                        )

            for group in output_groups:
                stats[group]["unique_bases"] += int(covered_by_group[group].sum())
                prefix = group_output_prefix(args.out_prefix, group)
                scores = scores_by_group[group]

                if args.score_mode == "on":
                    smoothed_track, score_track, _ = pns_scoring.scoring_track_names(
                        args.scoring_method
                    )
                    peak_track = smoothed_track if args.smooth_window > 0 else score_track
                    score_values = scores[peak_track][0][2]
                    coverage_values = scores["coverage"][0][2]
                    nucleosome_records = call_records(
                        scores=score_values,
                        chrom=region.contig,
                        adjusted_start=region.adjusted_start,
                        core_start=region.original_start,
                        core_end=region.original_end,
                        min_length=args.min_region_length,
                        max_nonpositive_run=args.max_neg_run,
                        flip_scores=False,
                        coverage_scores=coverage_values,
                    )
                    if args.peak_coverage_threshold is not None:
                        called_count = len(nucleosome_records)
                        nucleosome_records, _filtered_count = filter_records_by_coverage(
                            nucleosome_records,
                            coverage_values,
                            region.adjusted_start,
                            args.peak_coverage_threshold,
                        )
                        peak_filter_stats[group]["called"] += called_count
                        peak_filter_stats[group]["retained"] += len(nucleosome_records)
                    breakpoint_records = call_records(
                        scores=-1.0 * score_values,
                        chrom=region.contig,
                        adjusted_start=region.adjusted_start,
                        core_start=region.original_start,
                        core_end=region.original_end,
                        min_length=args.min_region_length,
                        max_nonpositive_run=args.max_neg_run,
                        flip_scores=True,
                        coverage_scores=coverage_values,
                    )
                    mode = "w" if chunk_index == 1 else "a"
                    write_records(
                        f"{prefix}_nucleosome_regions.bed",
                        nucleosome_records,
                        "nuc",
                        args.peak_score_scale,
                        mode,
                    )
                    write_records(
                        f"{prefix}_breakpoint_peaks.bed",
                        breakpoint_records,
                        "brk",
                        args.peak_score_scale,
                        mode,
                    )

                tracks = list(dict.fromkeys(args.score_tracks + args.other_tracks))
                write_tracks(
                    scores=scores,
                    contig=region.contig,
                    adjusted_start=region.adjusted_start,
                    original_start=region.original_start,
                    original_end=region.original_end,
                    tracks=tracks,
                    bigwig_handles=bigwig_handles[group],
                    wig_handles=wig_handles[group],
                )

        for group in output_groups:
            prefix = group_output_prefix(args.out_prefix, group)
            write_fragment_outputs(
                output_prefix=prefix,
                total_fragments_filtered=stats[group]["total_filtered"],
                total_fragments_used=stats[group]["total_used"],
                unique_bases_covered=stats[group]["unique_bases"],
                length_counts=stats[group]["length_counts"],
                dedup_scope=args.dedup_scope,
                max_duplicates=args.max_duplicates,
                max_per_coordinate=getattr(args, "max_per_coordinate", 0),
            )
            if args.dinuc_profile:
                dinucleotide.write_profile(
                    output_path=f"{prefix}_dinuc_profile.tsv",
                    accumulator=dinuc_accumulators[group],
                    positions=profile_positions,
                    fraction=args.dinuc_fraction,
                )

        if args.peak_coverage_threshold is not None:
            method_name = args.scoring_method.upper()
            threshold = format(args.peak_coverage_threshold, ".12g")
            for group in output_groups:
                called = peak_filter_stats[group]["called"]
                retained = peak_filter_stats[group]["retained"]
                filtered = called - retained
                label = "all" if group == "all" else group
                print(
                    f"[INFO] {method_name} nucleosome peak coverage filter ({label}) "
                    f">= {threshold}: retained {retained:,}/{called:,}; "
                    f"filtered {filtered:,}."
                )

        if args.split_ww_types:
            write_ww_summary(
                output_prefix=args.out_prefix,
                type_counts=type_counts,
                total_in_range=stats["all"]["total_used"],
            )

        if args.dinuc_profile:
            for group, accumulator in dinuc_accumulators.items():
                print(
                    f"[INFO] {group} dinucleotide fragments used: "
                    f"{accumulator['fragments_used']:,}; "
                    f"skipped: {accumulator['fragments_skipped']:,}."
                )

        if args.score_mode == "on" and args.peak_calling:
            peak_beds = []
            for group in output_groups:
                prefix = group_output_prefix(args.out_prefix, group)
                peak_beds.extend([
                    f"{prefix}_nucleosome_regions.bed",
                    f"{prefix}_breakpoint_peaks.bed",
                ])
            finalise_interval_files(
                peak_beds,
                args.interval_format,
                context.bigwig_header,
                bigbed_score_multiplier=args.bigbed_score_scale,
            )
        completed = True
    finally:
        for group in output_groups:
            close_track_handles(
                bigwig_handles[group], wig_handles[group], commit=completed
            )
        context.close()

    return 0
