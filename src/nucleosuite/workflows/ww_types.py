"""Standalone WW/SS type-classification workflow."""

from __future__ import annotations

import os
from collections import Counter, defaultdict

import numpy as np

from nucleosuite.core.reference import extract_reference_sequence
from nucleosuite.io.summaries import write_fragment_bed_rows, write_fragment_outputs
from nucleosuite.io.intervals import finalise_interval_files
from nucleosuite.io.tracks import (
    close_track_handles,
    open_track_handles,
    remove_stale_track_outputs,
    write_tracks,
)
from nucleosuite.scoring import basic_tracks
from nucleosuite.sequence import dinucleotide
from nucleosuite.plotting import remove_plot_variants
from nucleosuite.profile_plots import (
    plot_category_counts,
    plot_dinucleotide_profile,
    plot_ww_ss_profile,
    plot_ww_type_length_stacked,
)
from nucleosuite.sequence.ww_types import (
    ALL_OUTPUT_GROUPS,
    WW_TYPE_GROUPS,
    classify_fragment,
    group_output_prefix,
    write_summary,
    write_length_summary,
)
from nucleosuite.workflows.common import (
    default_output_prefix,
    ensure_coordinate_order,
    ensure_output_parent,
    prepare_fragment_run,
    input_paths_from_args,
    prepare_reference_if_needed,
    print_progress,
    randomize_fragments,
    set_random_seed,
)


# Private alias used by internal callers.
_ensure_coordinate_order = ensure_coordinate_order


def run(args) -> int:
    set_random_seed(args.seed)
    if args.out_prefix is None:
        args.out_prefix = default_output_prefix(input_paths_from_args(args), args.contigs)
    args.out_prefix = (
        f"{args.out_prefix}_wwtypes_lower{args.frag_lower}"
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
    output_groups = ALL_OUTPUT_GROUPS
    type_counts = Counter()
    total_filtered = 0
    total_used = 0
    unique_bases = 0
    length_counts = Counter()
    type_counts_by_length = defaultdict(Counter)

    combined_bed = f"{args.out_prefix}_ww_types.bed"
    stale = [
        combined_bed,
        f"{args.out_prefix}_ww_type_summary.tsv",
        f"{args.out_prefix}_fragment_summary.tsv",
        f"{args.out_prefix}_fragment_length_counts.tsv",
        f"{args.out_prefix}_ww_type_by_length.tsv",
    ]
    if args.split_beds:
        stale.extend(
            f"{args.out_prefix}_{group}.bed"
            for group in (*WW_TYPE_GROUPS, "unclassified")
        )
    stale.extend([combined_bed[:-4] + ".bb"])
    if args.split_beds:
        stale.extend(
            f"{args.out_prefix}_{group}.bb"
            for group in (*WW_TYPE_GROUPS, "unclassified")
        )
    for group in output_groups:
        prefix = group_output_prefix(args.out_prefix, group)
        stale.append(f"{prefix}_dinuc_profile.tsv")
        remove_stale_track_outputs(prefix, ["dyad"])
    for path in stale:
        if os.path.exists(path):
            os.remove(path)
    for plot_name in (
        f"{args.out_prefix}_ww_type_summary.png",
        f"{args.out_prefix}_ww_type_by_length_stacked.png",
    ):
        remove_plot_variants(plot_name)
    for group in output_groups:
        group_prefix = group_output_prefix(args.out_prefix, group)
        remove_plot_variants(f"{group_prefix}_dinuc_profile.png")
        remove_plot_variants(f"{group_prefix}_ww_ss_profile.png")

    profile_positions = dinucleotide.expected_profile_positions(
        args.frag_lower, args.frag_upper
    )
    accumulators = (
        {group: dinucleotide.new_accumulator() for group in output_groups}
        if args.dinuc_profile
        else {}
    )

    bigwig_handles = {group: {} for group in output_groups}
    wig_handles = {group: {} for group in output_groups}
    if args.dyad_tracks:
        for group in output_groups:
            prefix = group_output_prefix(args.out_prefix, group)
            bigwig, wig = open_track_handles(
                output_prefix=prefix,
                tracks=["dyad"],
                output_format=args.output_format,
                bigwig_header=context.bigwig_header,
            )
            bigwig_handles[group] = bigwig
            wig_handles[group] = wig

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

            arrays_by_group = (
                {
                    group: basic_tracks.new_arrays(
                        region.adjusted_end - region.adjusted_start
                    )
                    for group in output_groups
                }
                if args.dyad_tracks
                else {}
            )
            records = []
            for fragment_start, fragment_end in fragments:
                fragment_length = fragment_end - fragment_start
                if fragment_length not in fragment_length_set:
                    continue
                ww_type = classify_fragment(
                    fasta=context.fasta,
                    reference_context=reference_context,
                    fragment_start=fragment_start,
                    fragment_end=fragment_end,
                )
                records.append((fragment_start, fragment_end, ww_type))
                if args.dyad_tracks:
                    basic_tracks.add_fragment(
                        arrays=arrays_by_group["all"],
                        fragment_start=fragment_start,
                        fragment_end=fragment_end,
                        window_start=region.adjusted_start,
                        window_end=region.adjusted_end,
                        even_dyad=args.even_dyad,
                    )
                    if ww_type in WW_TYPE_GROUPS:
                        basic_tracks.add_fragment(
                            arrays=arrays_by_group[ww_type],
                            fragment_start=fragment_start,
                            fragment_end=fragment_end,
                            window_start=region.adjusted_start,
                            window_end=region.adjusted_end,
                            even_dyad=args.even_dyad,
                        )

            owned = _ensure_coordinate_order([
                record
                for record in records
                if region.original_start <= record[0] < region.original_end
            ])
            total_filtered += len(owned)
            total_used += len(owned)
            covered = np.zeros(region.original_end - region.original_start, dtype=bool)
            mode = "w" if chunk_index == 1 else "a"
            write_fragment_bed_rows(
                output_path=combined_bed,
                contig=region.contig,
                records=owned,
                mode=mode,
                include_type=True,
            )

            if args.split_beds:
                for group in (*WW_TYPE_GROUPS, "unclassified"):
                    selected = [
                        record
                        for record in owned
                        if (record[2] or "unclassified") == group
                    ]
                    write_fragment_bed_rows(
                        output_path=f"{args.out_prefix}_{group}.bed",
                        contig=region.contig,
                        records=selected,
                        mode=mode,
                        include_type=False,
                    )

            for fragment_start, fragment_end, ww_type in owned:
                fragment_length = fragment_end - fragment_start
                length_counts[fragment_length] += 1
                type_name = ww_type or "unclassified"
                type_counts[type_name] += 1
                type_counts_by_length[fragment_length][type_name] += 1
                overlap_start = max(fragment_start, region.original_start)
                overlap_end = min(fragment_end, region.original_end)
                if overlap_end > overlap_start:
                    covered[
                        overlap_start - region.original_start :
                        overlap_end - region.original_start
                    ] = True

                if args.dinuc_profile:
                    sequence = extract_reference_sequence(
                        fasta=context.fasta,
                        context=reference_context,
                        seq_start=fragment_start,
                        seq_end=fragment_end,
                    )
                    dinucleotide.add_fragment(
                        accumulators["all"], sequence, fragment_start, fragment_end
                    )
                    if ww_type in WW_TYPE_GROUPS:
                        dinucleotide.add_fragment(
                            accumulators[ww_type],
                            sequence,
                            fragment_start,
                            fragment_end,
                        )
            unique_bases += int(covered.sum())

            if args.dyad_tracks:
                for arrays in arrays_by_group.values():
                    basic_tracks.cap_sparse_arrays(arrays, args.max_per_coordinate)
                for group in output_groups:
                    scores = basic_tracks.to_scores(
                        arrays_by_group[group], region.contig, region.adjusted_start
                    )
                    write_tracks(
                        scores=scores,
                        contig=region.contig,
                        adjusted_start=region.adjusted_start,
                        original_start=region.original_start,
                        original_end=region.original_end,
                        tracks=["dyad"],
                        bigwig_handles=bigwig_handles[group],
                        wig_handles=wig_handles[group],
                    )

        write_summary(args.out_prefix, type_counts, total_used)
        plot_category_counts(
            f"{args.out_prefix}_ww_type_summary.tsv",
            f"{args.out_prefix}_ww_type_summary.png",
            title="WW/SS fragment types",
        )
        length_summary = write_length_summary(args.out_prefix, type_counts_by_length)
        plot_ww_type_length_stacked(
            length_summary,
            f"{args.out_prefix}_ww_type_by_length_stacked.png",
            title="WW/SS type frequencies by fragment length",
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
        if args.dinuc_profile:
            for group in output_groups:
                dinucleotide.write_profile(
                    output_path=(
                        f"{group_output_prefix(args.out_prefix, group)}"
                        "_dinuc_profile.tsv"
                    ),
                    accumulator=accumulators[group],
                    positions=profile_positions,
                    fraction=args.dinuc_fraction,
                )
                group_prefix = group_output_prefix(args.out_prefix, group)
                plot_dinucleotide_profile(
                    f"{group_prefix}_dinuc_profile.tsv",
                    f"{group_prefix}_dinuc_profile.png",
                    title="Dinucleotide profile",
                )
                plot_ww_ss_profile(
                    f"{group_prefix}_dinuc_profile.tsv",
                    f"{group_prefix}_ww_ss_profile.png",
                    title="WW/SS dinucleotide profile",
                )

        interval_beds = [combined_bed]
        if args.split_beds:
            interval_beds.extend(
                f"{args.out_prefix}_{group}.bed"
                for group in (*WW_TYPE_GROUPS, "unclassified")
            )
        finalise_interval_files(
            interval_beds, args.interval_format, context.bigwig_header
        )
        completed = True
    finally:
        if args.dyad_tracks:
            for group in output_groups:
                close_track_handles(
                    bigwig_handles[group], wig_handles[group], commit=completed
                )
        context.close()

    return 0
