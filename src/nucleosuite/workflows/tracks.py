"""Combined multi-range track generation.

The workflow reads each accepted fragment once per genomic chunk, assigns it to
all matching fragment-length ranges, and writes any requested PNS, WPS, basic
coordinate tracks and peak calls from the shared in-memory range accumulators.
"""

from __future__ import annotations

import csv
import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

from nucleosuite.io.bedgraph import (
    close_staged_bedgraph_handles,
    open_staged_bedgraph_handles,
    write_staged_bedgraph_tracks,
)
from nucleosuite.io.intervals import finalise_interval_files
from nucleosuite.io.summaries import write_fragment_bed_rows, write_fragment_outputs
from nucleosuite.io.tracks import (
    close_track_handles,
    open_track_handles,
    remove_stale_track_outputs,
    write_tracks,
)
from nucleosuite.plotting import remove_plot_variants
from nucleosuite.peaks.common import write_bed8_records
from nucleosuite.peaks.pns import call_records as call_pns_records
from nucleosuite.peaks.pns import write_records as write_pns_records
from nucleosuite.peaks.wps import call_records as call_wps_records
from nucleosuite.scoring import basic_tracks
from nucleosuite.scoring import pns as pns_scoring
from nucleosuite.scoring import wps as wps_scoring
from nucleosuite.core.reference import extract_reference_sequence
from nucleosuite.sequence import dinucleotide
from nucleosuite.sequence.ww_types import (
    ALL_OUTPUT_GROUPS,
    WW_TYPE_GROUPS,
    classify_core_sequence,
    group_output_prefix,
    write_length_summary,
    write_summary,
)
from nucleosuite.profile_plots import (
    plot_category_counts,
    plot_dinucleotide_profile,
    plot_ww_ss_profile,
    plot_ww_type_length_stacked,
)
from nucleosuite.workflows.common import prepare_reference_if_needed
from nucleosuite.workflows.common import (
    default_output_prefix,
    ensure_coordinate_order,
    input_paths_from_args,
    prepare_fragment_run,
    print_progress,
    set_random_seed,
)


PEAK_TOKENS = frozenset({"pns_peaks", "wps_peaks"})
SEQUENCE_TOKENS = frozenset({"dinuc_profile", "ww_types", "type_dyads"})
PNS_TOKENS = frozenset({"pns", "posPNS", "pns_smoothed", "pns_peaks"})
WPS_TOKENS = frozenset({"wps", "wps_smoothed", "mWPS", "sm_mWPS", "wps_peaks"})
BASIC_TOKENS = frozenset(basic_tracks.BASIC_TRACKS)
OUTPUT_TOKENS = BASIC_TOKENS | (PNS_TOKENS - PEAK_TOKENS) | (WPS_TOKENS - PEAK_TOKENS)
VALID_TOKENS = OUTPUT_TOKENS | PEAK_TOKENS | SEQUENCE_TOKENS
ALIASES = {
    "pospns": "posPNS",
    "positive-pns": "posPNS",
    "pns-smoothed": "pns_smoothed",
    "pns-peaks": "pns_peaks",
    "wps-smoothed": "wps_smoothed",
    "mwps": "mWPS",
    "sm-mwps": "sm_mWPS",
    "sm_mwps": "sm_mWPS",
    "wps-peaks": "wps_peaks",
    "fragment-ends": "fragment_ends",
    "left-end": "fragment_left_ends",
    "left-ends": "fragment_left_ends",
    "right-end": "fragment_right_ends",
    "right-ends": "fragment_right_ends",
    "dinuc-profile": "dinuc_profile",
    "ww-types": "ww_types",
    "type-dyads": "type_dyads",
}


@dataclass(frozen=True, order=True)
class FragmentRange:
    lower: int
    upper: int

    @property
    def label(self) -> str:
        return str(self.lower) if self.lower == self.upper else f"{self.lower}-{self.upper}"

    @property
    def directory_label(self) -> str:
        return f"exact_{self.lower}" if self.lower == self.upper else f"range_{self.lower}_{self.upper}"

    def contains(self, length: int) -> bool:
        return self.lower <= length <= self.upper


@dataclass
class OutputSpec:
    fragment_range: FragmentRange
    output_prefix: str
    tracks: tuple[str, ...]
    basic_scope: str = "range"
    bigwig_handles: dict = field(default_factory=dict)
    wig_handles: dict = field(default_factory=dict)
    type_bigwig_handles: dict = field(default_factory=dict)
    type_wig_handles: dict = field(default_factory=dict)
    bedgraph_handles: dict = field(default_factory=dict)
    type_bedgraph_handles: dict = field(default_factory=dict)

    @property
    def output_tracks(self) -> list[str]:
        return [track for track in self.tracks if track in OUTPUT_TOKENS]


@dataclass
class RangeState:
    fragment_range: FragmentRange
    need_basic_range: bool = False
    need_pns: bool = False
    need_wps: bool = False
    length_counts: Counter = field(default_factory=Counter)
    total_used: int = 0
    unique_bases: int = 0
    pns_distributions: tuple[dict, dict] | None = None
    wps_distributions: dict | None = None
    need_dinuc: bool = False
    need_ww_types: bool = False
    need_type_dyads: bool = False
    dinuc_accumulators: dict = field(default_factory=dict)
    type_counts: Counter = field(default_factory=Counter)
    type_counts_by_length: dict = field(default_factory=lambda: defaultdict(Counter))


def parse_fragment_range(text: str) -> FragmentRange:
    value = str(text).strip()
    if not value:
        raise ValueError("Fragment range cannot be empty")
    if "-" in value:
        left, right = value.split("-", 1)
    else:
        left = right = value
    try:
        lower, upper = int(left), int(right)
    except ValueError as error:
        raise ValueError(f"Invalid fragment range: {text!r}") from error
    if lower < 1 or upper < lower:
        raise ValueError(f"Require 1 <= lower <= upper in fragment range {text!r}")
    return FragmentRange(lower, upper)


def normalise_tracks(values: Iterable[str]) -> tuple[str, ...]:
    output: list[str] = []
    for raw in values:
        key = str(raw).strip()
        if not key:
            continue
        canonical = ALIASES.get(key.lower(), key)
        if canonical not in VALID_TOKENS:
            raise ValueError(
                f"Unknown track {raw!r}. Valid tracks: {', '.join(sorted(VALID_TOKENS))}"
            )
        if canonical not in output:
            output.append(canonical)
    if not output:
        raise ValueError("Each fragment range must request at least one track")
    return tuple(output)


def _parse_inline_spec(text: str, args) -> OutputSpec:
    if "=" not in text:
        raise ValueError(
            "--fragment-range requires RANGE=TRACKS, for example "
            "145=dyad,fragment_left_ends,fragment_right_ends"
        )
    range_text, tracks_text = text.split("=", 1)
    fragment_range = parse_fragment_range(range_text)
    tracks = normalise_tracks(tracks_text.split(","))
    directory = Path(args.output_dir) / fragment_range.directory_label
    prefix = directory / f"{args.output_prefix}_{fragment_range.directory_label}"
    return OutputSpec(fragment_range, str(prefix), tracks, "range")


def _read_spec_file(path: str) -> list[OutputSpec]:
    specs: list[OutputSpec] = []
    with open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        if not first:
            return specs
        handle.seek(0)
        first_fields = first.rstrip("\n").split("\t")
        has_header = "fragment_range" in first_fields and "output_prefix" in first_fields
        if has_header:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = reader
        else:
            rows = (
                {
                    "fragment_range": fields[0],
                    "output_prefix": fields[1],
                    "tracks": fields[2],
                    "basic_scope": fields[3] if len(fields) > 3 else "range",
                }
                for raw in handle
                if (fields := raw.rstrip("\n").split("\t")) and len(fields) >= 3
            )
        for row in rows:
            fragment_range = parse_fragment_range(row.get("fragment_range", ""))
            output_prefix = str(row.get("output_prefix", "")).strip()
            if not output_prefix:
                raise ValueError(f"{path}: output_prefix cannot be empty")
            tracks = normalise_tracks(str(row.get("tracks", "")).split(","))
            basic_scope = str(row.get("basic_scope", "range") or "range").strip().lower()
            if basic_scope not in {"range", "all"}:
                raise ValueError(f"{path}: basic_scope must be range or all")
            specs.append(OutputSpec(fragment_range, output_prefix, tracks, basic_scope))
    return specs


def load_specs(args) -> list[OutputSpec]:
    if args.output_prefix is None:
        args.output_prefix = default_output_prefix(input_paths_from_args(args), args.contigs)
    specs = [_parse_inline_spec(value, args) for value in args.fragment_range]
    for path in args.spec_file:
        specs.extend(_read_spec_file(path))
    seen_prefixes: set[str] = set()
    for spec in specs:
        if spec.output_prefix in seen_prefixes:
            raise ValueError(f"Duplicate output prefix: {spec.output_prefix}")
        seen_prefixes.add(spec.output_prefix)
        if "pns_peaks" in spec.tracks and "wps_peaks" in spec.tracks:
            raise ValueError(
                f"{spec.output_prefix}: PNS and WPS peaks cannot share one output prefix"
            )
        if "pns_smoothed" in spec.tracks and args.pns_smooth_window == 0:
            raise ValueError(
                f"{spec.output_prefix}: pns_smoothed requires --pns-smooth-window"
            )
    return specs


def _range_states(specs: list[OutputSpec], args) -> dict[FragmentRange, RangeState]:
    states: dict[FragmentRange, RangeState] = {}
    for spec in specs:
        state = states.setdefault(spec.fragment_range, RangeState(spec.fragment_range))
        state.need_basic_range |= bool(set(spec.tracks) & BASIC_TOKENS) and spec.basic_scope == "range"
        state.need_pns |= bool(set(spec.tracks) & PNS_TOKENS)
        state.need_wps |= bool(set(spec.tracks) & WPS_TOKENS)
        state.need_dinuc |= "dinuc_profile" in spec.tracks or "ww_types" in spec.tracks
        state.need_ww_types |= "ww_types" in spec.tracks or "type_dyads" in spec.tracks
        state.need_type_dyads |= "type_dyads" in spec.tracks
    for state in states.values():
        lengths = range(state.fragment_range.lower, state.fragment_range.upper + 1)
        if state.need_pns:
            state.pns_distributions = pns_scoring.precompute_distributions(
                lengths, args.pns_mode_length
            )
        if state.need_dinuc:
            groups = ALL_OUTPUT_GROUPS if state.need_ww_types else ("all",)
            state.dinuc_accumulators = {group: dinucleotide.new_accumulator() for group in groups}
        if state.need_wps:
            state.wps_distributions = wps_scoring.precompute_distributions(
                lengths, args.wps_protection
            )
    return states


def _remove_outputs(spec: OutputSpec) -> None:
    Path(spec.output_prefix).parent.mkdir(parents=True, exist_ok=True)
    remove_stale_track_outputs(spec.output_prefix, spec.output_tracks)
    for suffix in (
        "_fragment_summary.tsv",
        "_fragment_length_counts.tsv",
        "_nucleosome_regions.bed",
        "_breakpoint_peaks.bed",
        "_nucleosome_regions.bb",
        "_breakpoint_peaks.bb",
        "_dinuc_profile.tsv", "_dinuc_profile_counts.tsv",
        "_ww_type_summary.tsv",
        "_ww_type_by_length.tsv",
        "_ww_types.bed", "_ww_types.bb",
    ):
        try:
            os.remove(spec.output_prefix + suffix)
        except FileNotFoundError:
            pass
    for suffix in (
        "_fragment_length_distribution.png",
        "_dinuc_profile.png", "_ww_ss_profile.png",
        "_ww_type_summary.png", "_ww_type_by_length_stacked.png",
    ):
        remove_plot_variants(spec.output_prefix + suffix)


def _write_pns_peaks(spec, scores, region, args, mode):
    peak_track = "pns_smoothed" if args.pns_smooth_window > 0 else "pns"
    values = scores[peak_track][0][2]
    coverage = scores["coverage"][0][2]
    nuc = call_pns_records(
        scores=values,
        chrom=region.contig,
        adjusted_start=region.adjusted_start,
        core_start=region.original_start,
        core_end=region.original_end,
        min_length=args.pns_min_region_length,
        max_nonpositive_run=args.pns_max_neg_run,
        coverage_scores=coverage,
    )
    brk = call_pns_records(
        scores=-1.0 * values,
        chrom=region.contig,
        adjusted_start=region.adjusted_start,
        core_start=region.original_start,
        core_end=region.original_end,
        min_length=args.pns_min_region_length,
        max_nonpositive_run=args.pns_max_neg_run,
        flip_scores=True,
        coverage_scores=coverage,
    )
    blacklist = getattr(args, "_blacklist_index", None)
    if blacklist is not None:
        nuc = [r for r in nuc if not blacklist.overlaps(r["chrom"], r["region_start"], r["region_end"])]
        brk = [r for r in brk if not blacklist.overlaps(r["chrom"], r["region_start"], r["region_end"])]
    write_pns_records(
        spec.output_prefix + "_nucleosome_regions.bed",
        nuc,
        "nuc",
        args.pns_peak_score_scale,
        mode,
    )
    write_pns_records(
        spec.output_prefix + "_breakpoint_peaks.bed",
        brk,
        "brk",
        args.pns_peak_score_scale,
        mode,
    )


def _write_wps_peaks(spec, scores, region, args, mode):
    values = scores[args.wps_peak_track][0][2]
    nuc = call_wps_records(
        scores=values,
        chrom=region.contig,
        adjusted_start=region.adjusted_start,
        core_start=region.original_start,
        core_end=region.original_end,
        merge_gap_bp=args.wps_peak_merge_gap,
        min_length=args.wps_peak_minlen,
        max_length=args.wps_peak_maxlen,
        max_region_length=args.wps_peak_maxregion,
        score_cutoff=args.wps_peak_varicutoff,
    )
    brk = call_wps_records(
        scores=-1.0 * values,
        chrom=region.contig,
        adjusted_start=region.adjusted_start,
        core_start=region.original_start,
        core_end=region.original_end,
        merge_gap_bp=args.wps_peak_merge_gap,
        min_length=args.wps_peak_minlen,
        max_length=args.wps_peak_maxlen,
        max_region_length=args.wps_peak_maxregion,
        score_cutoff=args.wps_peak_varicutoff,
        flip_scores=True,
    )
    blacklist = getattr(args, "_blacklist_index", None)
    if blacklist is not None:
        nuc = [r for r in nuc if not blacklist.overlaps(r["chrom"], r["region_start"], r["region_end"])]
        brk = [r for r in brk if not blacklist.overlaps(r["chrom"], r["region_start"], r["region_end"])]
    write_bed8_records(
        spec.output_prefix + "_nucleosome_regions.bed",
        nuc,
        "nuc",
        args.wps_peak_score_scale,
        mode,
    )
    write_bed8_records(
        spec.output_prefix + "_breakpoint_peaks.bed",
        brk,
        "brk",
        args.wps_peak_score_scale,
        mode,
    )



def _write_sequence_outputs(spec, state, args, context):
    positions = dinucleotide.expected_profile_positions(
        state.fragment_range.lower, state.fragment_range.upper
    )
    groups = ALL_OUTPUT_GROUPS if state.need_ww_types else ("all",)
    if state.need_dinuc:
        for group in groups:
            prefix = group_output_prefix(spec.output_prefix, group)
            path = f"{prefix}_dinuc_profile.tsv"
            dinucleotide.write_profile(path, state.dinuc_accumulators[group], positions)
            plot_dinucleotide_profile(path, f"{prefix}_dinuc_profile.png", title=f"{group} dinucleotide profile")
            plot_ww_ss_profile(path, f"{prefix}_ww_ss_profile.png", title=f"{group} WW/SS profile")
    if state.need_ww_types:
        write_summary(spec.output_prefix, state.type_counts, state.total_used)
        plot_category_counts(
            f"{spec.output_prefix}_ww_type_summary.tsv",
            f"{spec.output_prefix}_ww_type_summary.png",
            title="WW/SS fragment types",
        )
        table = write_length_summary(spec.output_prefix, state.type_counts_by_length)
        plot_ww_type_length_stacked(
            table,
            f"{spec.output_prefix}_ww_type_by_length_stacked.png",
            title="WW/SS type frequencies by fragment length",
        )
        beds = [f"{spec.output_prefix}_ww_types.bed"]
        beds.extend(f"{spec.output_prefix}_{group}.bed" for group in (*WW_TYPE_GROUPS, "unclassified"))
        finalise_interval_files(beds, args.interval_format, context.bigwig_header)



def _sequence_features_for_fragment(
    *,
    fasta,
    reference_context,
    fragment_start: int,
    fragment_end: int,
    need_ww_type: bool,
) -> tuple[str | None, str | None]:
    """Fetch one sequence window and derive all sequence-dependent features.

    The fetched window covers both the complete fragment and, when needed, the
    dyad-centred 147-bp WW/SS classification core.  The caller caches this tuple
    for the lifetime of one genomic chunk.
    """
    if not need_ww_type:
        sequence = extract_reference_sequence(
            fasta=fasta,
            context=reference_context,
            seq_start=fragment_start,
            seq_end=fragment_end,
        )
        return sequence, None

    dyad = dinucleotide.fragment_dyad(fragment_start, fragment_end)
    core_start = dyad - 73
    core_end = dyad + 74
    fetch_start = min(fragment_start, core_start)
    fetch_end = max(fragment_end, core_end)
    window = extract_reference_sequence(
        fasta=fasta,
        context=reference_context,
        seq_start=fetch_start,
        seq_end=fetch_end,
    )
    if window is None:
        return None, None
    fragment_sequence = window[fragment_start - fetch_start : fragment_end - fetch_start]
    core = window[core_start - fetch_start : core_end - fetch_start]
    if len(fragment_sequence) != fragment_end - fragment_start:
        fragment_sequence = None
    return fragment_sequence, classify_core_sequence(core)

def run(args) -> int:
    set_random_seed(args.seed)
    specs = load_specs(args)
    if any(set(spec.tracks) & SEQUENCE_TOKENS for spec in specs) and not args.fasta:
        raise ValueError("dinuc_profile, ww_types and type_dyads require --fasta")
    if any("wps_peaks" in spec.tracks for spec in specs) and args.wps_peak_track == "sm_mWPS":
        required_overlap = (
            args.wps_peak_maxregion
            + args.wps_baseline_window // 2
            + (args.wps_sg_window // 2 if args.wps_sg_window else 0)
        )
        if args.overlap_bp < required_overlap:
            raise ValueError(
                "WPS peak calling requires --overlap-bp >= "
                f"{required_overlap} for the selected peak and preprocessing windows"
            )
    states = _range_states(specs, args)
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
    source = getattr(context, "source", None)
    args._blacklist_index = getattr(source, "blacklist", None)
    total_filtered = 0
    any_basic_all = any(
        bool(set(spec.tracks) & BASIC_TOKENS) and spec.basic_scope == "all"
        for spec in specs
    )
    if getattr(args, "staged_bedgraph_root", None):
        print(
            "Staging validated bedGraphs during track generation under "
            f"{Path(args.staged_bedgraph_root).resolve()}"
        )

    completed = False
    try:
        for spec in specs:
            _remove_outputs(spec)
            spec.bigwig_handles, spec.wig_handles = open_track_handles(
                output_prefix=spec.output_prefix,
                tracks=spec.output_tracks,
                output_format=args.output_format,
                bigwig_header=context.bigwig_header,
            )
            spec.bedgraph_handles = open_staged_bedgraph_handles(
                output_prefix=spec.output_prefix,
                tracks=spec.output_tracks,
                staging_root=getattr(args, "staged_bedgraph_root", None),
                source_root=getattr(args, "staged_bedgraph_source_root", None),
                source_id=getattr(args, "staged_bedgraph_source_id", None),
                chrom_order=[chrom for chrom, _length in context.bigwig_header],
            )
            if "type_dyads" in spec.tracks:
                for group in ALL_OUTPUT_GROUPS:
                    bw, wig = open_track_handles(
                        output_prefix=group_output_prefix(spec.output_prefix, group),
                        tracks=["dyad"], output_format=args.output_format,
                        bigwig_header=context.bigwig_header,
                    )
                    spec.type_bigwig_handles[group] = bw
                    spec.type_wig_handles[group] = wig
                    spec.type_bedgraph_handles[group] = open_staged_bedgraph_handles(
                        output_prefix=group_output_prefix(spec.output_prefix, group),
                        tracks=["dyad"],
                        staging_root=getattr(args, "staged_bedgraph_root", None),
                        source_root=getattr(args, "staged_bedgraph_source_root", None),
                        source_id=getattr(args, "staged_bedgraph_source_id", None),
                        chrom_order=[chrom for chrom, _length in context.bigwig_header],
                    )

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
            need_sequence = any(state.need_dinuc or state.need_ww_types for state in states.values())
            reference_context = prepare_reference_if_needed(
                fasta=getattr(context, "fasta", None), contig=region.contig,
                start=region.adjusted_start, end=region.adjusted_end,
                required=need_sequence,
            )
            reference_length = region.adjusted_end - region.adjusted_start
            basic_all = basic_tracks.new_arrays(reference_length) if any_basic_all else None
            basic_by_range: dict[FragmentRange, dict] = {}
            pns_by_range: dict[FragmentRange, dict] = {}
            wps_by_range: dict[FragmentRange, np.ndarray] = {}
            type_arrays_by_range: dict[FragmentRange, dict] = {}
            covered_by_range = {
                key: np.zeros(region.original_end - region.original_start, dtype=bool)
                for key in states
            }
            for key, state in states.items():
                if state.need_basic_range or state.need_pns:
                    # PNS peak calling requires range-specific coverage even when it
                    # is not selected as an output track.
                    basic_by_range[key] = basic_tracks.new_arrays(reference_length)
                if state.need_pns:
                    pns_by_range[key] = pns_scoring.new_arrays(reference_length)
                if state.need_wps:
                    wps_by_range[key] = wps_scoring.new_array(reference_length)
                if state.need_type_dyads:
                    type_arrays_by_range[key] = {group: basic_tracks.new_arrays(reference_length) for group in ALL_OUTPUT_GROUPS}

            sequence_features: dict[tuple[int, int], tuple[str | None, str | None]] = {}
            if need_sequence:
                for fragment_start, fragment_end in fragments:
                    length = fragment_end - fragment_start
                    matching_states = [
                        state for key, state in states.items() if key.contains(length)
                    ]
                    if not any(state.need_dinuc or state.need_ww_types for state in matching_states):
                        continue
                    need_ww_type = any(state.need_ww_types for state in matching_states)
                    fragment_key = (fragment_start, fragment_end)
                    if fragment_key not in sequence_features:
                        sequence_features[fragment_key] = _sequence_features_for_fragment(
                            fasta=context.fasta,
                            reference_context=reference_context,
                            fragment_start=fragment_start,
                            fragment_end=fragment_end,
                            need_ww_type=need_ww_type,
                        )

            for fragment_start, fragment_end in fragments:
                length = fragment_end - fragment_start
                if basic_all is not None:
                    basic_tracks.add_fragment(
                        basic_all,
                        fragment_start,
                        fragment_end,
                        region.adjusted_start,
                        region.adjusted_end,
                        args.even_dyad,
                    )
                for key, state in states.items():
                    if not key.contains(length):
                        continue
                    if key in basic_by_range:
                        basic_tracks.add_fragment(
                            basic_by_range[key],
                            fragment_start,
                            fragment_end,
                            region.adjusted_start,
                            region.adjusted_end,
                            args.even_dyad,
                        )
                    if state.need_pns:
                        centred, positive = state.pns_distributions or ({}, {})
                        pns_scoring.add_fragment(
                            pns_by_range[key],
                            fragment_start,
                            fragment_end,
                            region.adjusted_start,
                            region.adjusted_end,
                            args.pns_mode_length,
                            centred,
                            positive,
                        )
                    if state.need_wps:
                        wps_scoring.add_fragment(
                            wps_by_range[key],
                            fragment_start,
                            fragment_end,
                            region.adjusted_start,
                            args.wps_protection,
                            state.wps_distributions or {},
                        )
                    if state.need_type_dyads:
                        _sequence, ww_type = sequence_features[(fragment_start, fragment_end)]
                        basic_tracks.add_fragment(
                            type_arrays_by_range[key]["all"], fragment_start, fragment_end,
                            region.adjusted_start, region.adjusted_end, args.even_dyad,
                        )
                        if ww_type in WW_TYPE_GROUPS:
                            basic_tracks.add_fragment(
                                type_arrays_by_range[key][ww_type], fragment_start, fragment_end,
                                region.adjusted_start, region.adjusted_end, args.even_dyad,
                            )

            if basic_all is not None:
                basic_tracks.cap_sparse_arrays(basic_all, args.max_per_coordinate)
            for arrays in basic_by_range.values():
                basic_tracks.cap_sparse_arrays(arrays, args.max_per_coordinate)
            for grouped in type_arrays_by_range.values():
                for arrays in grouped.values():
                    basic_tracks.cap_sparse_arrays(arrays, args.max_per_coordinate)

            scores_by_range: dict[FragmentRange, dict] = {}
            all_basic_scores = (
                basic_tracks.to_scores(basic_all, region.contig, region.adjusted_start)
                if basic_all is not None
                else {}
            )
            for key, state in states.items():
                scores = (
                    basic_tracks.to_scores(
                        basic_by_range[key], region.contig, region.adjusted_start
                    )
                    if key in basic_by_range
                    else {}
                )
                if state.need_pns:
                    scores.update(
                        pns_scoring.to_scores(
                            pns_by_range[key],
                            region.contig,
                            region.adjusted_start,
                            args.pns_smooth_window,
                            args.pns_smooth_order,
                        )
                    )
                if state.need_wps:
                    scores.update(
                        wps_scoring.to_scores(
                            wps_by_range[key],
                            region.contig,
                            region.adjusted_start,
                            args.wps_baseline_window,
                            args.wps_sg_window,
                            args.wps_sg_order,
                        )
                    )
                scores_by_range[key] = scores

            owned = ensure_coordinate_order([
                fragment
                for fragment in fragments
                if region.original_start <= fragment[0] < region.original_end
            ])
            total_filtered += len(owned)
            for fragment_start, fragment_end in owned:
                length = fragment_end - fragment_start
                for key, state in states.items():
                    if not key.contains(length):
                        continue
                    state.total_used += 1
                    state.length_counts[length] += 1
                    if state.need_dinuc or state.need_ww_types:
                        sequence, cached_ww_type = sequence_features[(fragment_start, fragment_end)]
                        ww_type = cached_ww_type if state.need_ww_types else None
                        if state.need_ww_types:
                            type_name = ww_type or "unclassified"
                            state.type_counts[type_name] += 1
                            state.type_counts_by_length[length][type_name] += 1
                        if state.need_dinuc:
                            dinucleotide.add_fragment(
                                state.dinuc_accumulators["all"], sequence,
                                fragment_start, fragment_end,
                            )
                            if state.need_ww_types and ww_type in WW_TYPE_GROUPS:
                                dinucleotide.add_fragment(
                                    state.dinuc_accumulators[ww_type], sequence,
                                    fragment_start, fragment_end,
                                )
                    overlap_start = max(fragment_start, region.original_start)
                    overlap_end = min(fragment_end, region.original_end)
                    if overlap_end > overlap_start:
                        covered_by_range[key][
                            overlap_start - region.original_start :
                            overlap_end - region.original_start
                        ] = True
            for key, state in states.items():
                state.unique_bases += int(covered_by_range[key].sum())

            mode = "w" if chunk_index == 1 else "a"
            for spec in specs:
                state = states[spec.fragment_range]
                if "ww_types" in spec.tracks:
                    typed_records = []
                    for fragment_start, fragment_end in owned:
                        if not spec.fragment_range.contains(fragment_end-fragment_start):
                            continue
                        _sequence, ww_type = sequence_features[(fragment_start, fragment_end)]
                        typed_records.append((fragment_start, fragment_end, ww_type))
                    write_fragment_bed_rows(
                        f"{spec.output_prefix}_ww_types.bed", region.contig,
                        typed_records, mode=mode, include_type=True,
                    )
                    for group in (*WW_TYPE_GROUPS, "unclassified"):
                        selected = [r for r in typed_records if (r[2] or "unclassified") == group]
                        write_fragment_bed_rows(
                            f"{spec.output_prefix}_{group}.bed", region.contig,
                            selected, mode=mode, include_type=False,
                        )
                scores = dict(scores_by_range[spec.fragment_range])
                if spec.basic_scope == "all":
                    scores.update(all_basic_scores)
                if "pns_peaks" in spec.tracks:
                    _write_pns_peaks(spec, scores, region, args, mode)
                if "wps_peaks" in spec.tracks:
                    _write_wps_peaks(spec, scores, region, args, mode)
                if "type_dyads" in spec.tracks:
                    for group in ALL_OUTPUT_GROUPS:
                        group_scores = basic_tracks.to_scores(
                            type_arrays_by_range[spec.fragment_range][group],
                            region.contig, region.adjusted_start,
                        )
                        write_tracks(
                            scores=group_scores, contig=region.contig,
                            adjusted_start=region.adjusted_start,
                            original_start=region.original_start,
                            original_end=region.original_end, tracks=["dyad"],
                            bigwig_handles=spec.type_bigwig_handles[group],
                            wig_handles=spec.type_wig_handles[group],
                        )
                        write_staged_bedgraph_tracks(
                            scores=group_scores,
                            contig=region.contig,
                            adjusted_start=region.adjusted_start,
                            original_start=region.original_start,
                            original_end=region.original_end,
                            tracks=["dyad"],
                            handles=spec.type_bedgraph_handles.get(group, {}),
                        )
                write_tracks(
                    scores=scores,
                    contig=region.contig,
                    adjusted_start=region.adjusted_start,
                    original_start=region.original_start,
                    original_end=region.original_end,
                    tracks=spec.output_tracks,
                    bigwig_handles=spec.bigwig_handles,
                    wig_handles=spec.wig_handles,
                )
                write_staged_bedgraph_tracks(
                    scores=scores,
                    contig=region.contig,
                    adjusted_start=region.adjusted_start,
                    original_start=region.original_start,
                    original_end=region.original_end,
                    tracks=spec.output_tracks,
                    handles=spec.bedgraph_handles,
                )

        for spec in specs:
            state = states[spec.fragment_range]
            write_fragment_outputs(
                output_prefix=spec.output_prefix,
                total_fragments_filtered=total_filtered,
                total_fragments_used=state.total_used,
                unique_bases_covered=state.unique_bases,
                length_counts=state.length_counts,
                dedup_scope=args.dedup_scope,
                max_duplicates=args.max_duplicates,
                max_per_coordinate=args.max_per_coordinate,
            )
            if set(spec.tracks) & SEQUENCE_TOKENS:
                _write_sequence_outputs(spec, state, args, context)
            if set(spec.tracks) & PEAK_TOKENS:
                finalise_interval_files(
                    [
                        spec.output_prefix + "_nucleosome_regions.bed",
                        spec.output_prefix + "_breakpoint_peaks.bed",
                    ],
                    args.interval_format,
                    context.bigwig_header,
                    bigbed_score_multiplier=(
                        args.bigbed_score_scale if "pns_peaks" in spec.tracks else 1.0
                    ),
                )

        if args.report:
            report = Path(args.report)
            report.parent.mkdir(parents=True, exist_ok=True)
            with report.open("wt", encoding="utf-8") as handle:
                source = getattr(context, "source", None)
                blacklist = getattr(source, "blacklist", None)
                blacklist_path = getattr(blacklist, "path", "") if blacklist else ""
                blacklist_intervals = (
                    blacklist.summary.interval_count if blacklist is not None else 0
                )
                blacklisted_bases = (
                    blacklist.summary.blacklisted_bases if blacklist is not None else 0
                )
                fragments_excluded = getattr(
                    source, "fragments_excluded", 0
                )
                handle.write(
                    "fragment_range\toutput_prefix\ttracks\tbasic_scope\t"
                    "max_duplicates\tmax_per_coordinate\tblacklist_file\t"
                    "blacklist_intervals\tblacklisted_bases\t"
                    "fragments_excluded_by_blacklist\n"
                )
                for spec in specs:
                    handle.write(
                        f"{spec.fragment_range.label}\t{spec.output_prefix}\t"
                        f"{','.join(spec.tracks)}\t{spec.basic_scope}\t"
                        f"{args.max_duplicates}\t{args.max_per_coordinate}\t"
                        f"{blacklist_path}\t{blacklist_intervals}\t"
                        f"{blacklisted_bases}\t{fragments_excluded}\n"
                    )
        completed = True
    finally:
        close_error = None
        for spec in specs:
            try:
                close_track_handles(
                    spec.bigwig_handles, spec.wig_handles, commit=completed
                )
            except Exception as error:
                close_error = close_error or error
            try:
                close_staged_bedgraph_handles(
                    spec.bedgraph_handles, commit=completed
                )
            except Exception as error:
                close_error = close_error or error
            for group in spec.type_bigwig_handles:
                try:
                    close_track_handles(
                        spec.type_bigwig_handles[group],
                        spec.type_wig_handles[group],
                        commit=completed,
                    )
                except Exception as error:
                    close_error = close_error or error
                try:
                    close_staged_bedgraph_handles(
                        spec.type_bedgraph_handles.get(group, {}), commit=completed
                    )
                except Exception as error:
                    close_error = close_error or error
        context.close()
        if close_error is not None:
            raise close_error
    return 0
