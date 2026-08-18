"""Call PNS-style or WPS-style peaks directly from a BigWig signal."""

from __future__ import annotations

import os

import numpy as np
from scipy.signal import savgol_filter

try:
    import pyBigWig
except ImportError:  # pragma: no cover
    pyBigWig = None

from nucleosuite.core.regions import build_processing_regions, expand_contig_tokens
from nucleosuite.io.intervals import finalise_interval_files
from nucleosuite.peaks.common import prepare_output, write_bed8_records
from nucleosuite.peaks.pns import call_records as call_pns_records
from nucleosuite.peaks.pns import write_records as write_pns_records
from nucleosuite.peaks.wps import call_records as call_wps_records
from nucleosuite.scoring import wps as wps_scoring
from nucleosuite.workflows.common import ensure_output_parent, print_progress


def smooth_signal(values: np.ndarray, window: int, order: int) -> np.ndarray:
    if window <= 0 or values.size == 0:
        return values.copy()
    if window % 2 == 0:
        window += 1
    if window < 3 or values.size < window:
        return values.copy()
    if order >= window:
        raise ValueError("--smooth-order must be smaller than --smooth-window")
    return savgol_filter(values, window_length=window, polyorder=order)


def get_values(bigwig, contig: str, start: int, end: int) -> np.ndarray:
    values = bigwig.values(contig, start, end, numpy=True)
    if values is None:
        return np.zeros(end - start, dtype=float)
    values = np.asarray(values, dtype=float)
    values[~np.isfinite(values)] = 0.0
    return values


def run(args) -> int:
    if pyBigWig is None:
        raise RuntimeError(
            "The call-peaks command requires pyBigWig. Install it with "
            "'conda install -c bioconda pybigwig'."
        )
    if args.out_prefix is None:
        args.out_prefix = os.path.splitext(os.path.basename(args.input_bigwig))[0]
    ensure_output_parent(args.out_prefix)

    bigwig = pyBigWig.open(args.input_bigwig)
    if bigwig is None:
        raise OSError(f"Could not open BigWig: {args.input_bigwig}")

    nucleosome_path = f"{args.out_prefix}_nucleosome_regions.bed"
    breakpoint_path = f"{args.out_prefix}_breakpoint_peaks.bed"
    if args.signal in {"both", "nucleosome"}:
        prepare_output(nucleosome_path)
    if args.signal in {"both", "breakpoint"}:
        prepare_output(breakpoint_path)

    excluded_by_blacklist = 0
    try:
        chrom_sizes = bigwig.chroms()
        references = list(chrom_sizes)
        lengths = [int(chrom_sizes[name]) for name in references]
        from nucleosuite.core.blacklist import load_blacklist
        blacklist = load_blacklist(
            getattr(args, "blacklist_bed", None), references, lengths
        )
        selected_specs = expand_contig_tokens(args.regions, references)
        regions, _selected_names = build_processing_regions(
            selected_specs=selected_specs,
            references=references,
            lengths=lengths,
            chunk_bp=args.chunk_bp,
            overlap_bp=args.overlap_bp,
        )

        total_chunks = len(regions)
        previous_contig = None
        for chunk_index, region in enumerate(regions, start=1):
            print_progress(
                chunk_index, total_chunks, region, previous_contig=previous_contig
            )
            previous_contig = region.contig
            raw_signal = get_values(
                bigwig,
                region.contig,
                region.adjusted_start,
                region.adjusted_end,
            )
            if args.method == "pns":
                signal = smooth_signal(
                    raw_signal,
                    args.smooth_window,
                    args.smooth_order,
                )
            elif args.wps_input_mode == "raw":
                signal = wps_scoring.kircher_adjusted_signal(
                    raw_signal,
                    baseline_window=args.wps_baseline_window,
                    smooth_window=args.smooth_window,
                    smooth_order=args.smooth_order,
                )
            else:
                # sm_mWPS or another already adjusted track is evaluated without
                # a second smoothing or median-subtraction pass.
                signal = raw_signal
            if blacklist is not None:
                blacklist.mask_values(
                    region.contig,
                    region.adjusted_start,
                    signal,
                )
            mode = "w" if chunk_index == 1 else "a"

            def keep_allowed(records):
                nonlocal excluded_by_blacklist
                if blacklist is None:
                    return records
                kept = []
                for record in records:
                    if blacklist.overlaps(
                        record["chrom"],
                        int(record["region_start"]),
                        int(record["region_end"]),
                    ):
                        excluded_by_blacklist += 1
                    else:
                        kept.append(record)
                return kept

            if args.method == "pns":
                if args.signal in {"both", "nucleosome"}:
                    records = call_pns_records(
                        scores=signal,
                        chrom=region.contig,
                        adjusted_start=region.adjusted_start,
                        core_start=region.original_start,
                        core_end=region.original_end,
                        min_length=args.min_region_length,
                        max_nonpositive_run=args.max_neg_run,
                    )
                    records = keep_allowed(records)
                    write_pns_records(
                        nucleosome_path,
                        records,
                        "nuc",
                        args.score_scale,
                        mode,
                    )
                if args.signal in {"both", "breakpoint"}:
                    records = call_pns_records(
                        scores=-1.0 * signal,
                        chrom=region.contig,
                        adjusted_start=region.adjusted_start,
                        core_start=region.original_start,
                        core_end=region.original_end,
                        min_length=args.min_region_length,
                        max_nonpositive_run=args.max_neg_run,
                        flip_scores=True,
                    )
                    records = keep_allowed(records)
                    write_pns_records(
                        breakpoint_path,
                        records,
                        "brk",
                        args.score_scale,
                        mode,
                    )
            else:
                if args.signal in {"both", "nucleosome"}:
                    records = call_wps_records(
                        scores=signal,
                        chrom=region.contig,
                        adjusted_start=region.adjusted_start,
                        core_start=region.original_start,
                        core_end=region.original_end,
                        merge_gap_bp=args.wps_merge_gap,
                        min_length=args.wps_min_length,
                        max_length=args.wps_max_length,
                        max_region_length=args.wps_max_region,
                        score_cutoff=args.wps_score_cutoff,
                    )
                    records = keep_allowed(records)
                    write_bed8_records(
                        nucleosome_path,
                        records,
                        "nuc",
                        args.score_scale,
                        mode,
                    )
                if args.signal in {"both", "breakpoint"}:
                    records = call_wps_records(
                        scores=-1.0 * signal,
                        chrom=region.contig,
                        adjusted_start=region.adjusted_start,
                        core_start=region.original_start,
                        core_end=region.original_end,
                        merge_gap_bp=args.wps_merge_gap,
                        min_length=args.wps_min_length,
                        max_length=args.wps_max_length,
                        max_region_length=args.wps_max_region,
                        score_cutoff=args.wps_score_cutoff,
                        flip_scores=True,
                    )
                    records = keep_allowed(records)
                    write_bed8_records(
                        breakpoint_path,
                        records,
                        "brk",
                        args.score_scale,
                        mode,
                    )
        interval_paths = []
        if args.signal in {"both", "nucleosome"}:
            interval_paths.append(nucleosome_path)
        if args.signal in {"both", "breakpoint"}:
            interval_paths.append(breakpoint_path)
        finalise_interval_files(
            interval_paths,
            args.interval_format,
            chrom_sizes,
            bigbed_score_multiplier=args.bigbed_score_scale if args.method == "pns" else 1.0,
        )
        if blacklist is not None:
            print(f"Blacklist-overlapping calls discarded: {excluded_by_blacklist:,}")
    finally:
        bigwig.close()
    return 0
