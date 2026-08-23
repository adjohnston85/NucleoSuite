#!/usr/bin/env python3
"""Extract, combine and filter genomic fragment intervals."""

from __future__ import annotations

import argparse
import glob
import gzip
import os
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

try:
    import pysam
except ImportError:  # allow --help without runtime dependency
    pysam = None

from nucleosuite.core.bam_headers import merge_bam_reference_headers_with_aliases
from nucleosuite.core.fragment_inputs import IntervalFragmentSource
from nucleosuite.core.fragments import generate_fragment_ranges, require_bam_indexes
from nucleosuite.core.regions import expand_contig_tokens
from nucleosuite.io.intervals import convert_bed_to_bigbed
from nucleosuite.profile_plots import plot_count_profile
from nucleosuite.progress import ProgressReporter


def expand_inputs(values: Sequence[str]) -> list[str]:
    output: list[str] = []
    for value in values:
        matches = sorted(glob.glob(value))
        output.extend(matches or [value])
    output = list(dict.fromkeys(output))
    missing = [path for path in output if not Path(path).exists()]
    if missing:
        raise FileNotFoundError("Input file(s) not found: " + ", ".join(missing))
    return output


def _open_output(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return path.open("wt", encoding="utf-8")


def _output_suffix(output_format: str) -> str:
    return ".fragments.bed.gz" if output_format == "bed.gz" else ".fragments.bed"


def _write_records(path: Path, records: Iterable[tuple[str, int, int]]) -> int:
    count = 0
    with _open_output(path) as handle:
        for chrom, start, end in records:
            handle.write(f"{chrom}\t{start}\t{end}\n")
            count += 1
    return count


def _bam_fragments_by_contig(args, paths: list[str]):
    reporter = ProgressReporter("fragments")
    require_bam_indexes(paths)
    handles = [pysam.AlignmentFile(path, "rb") for path in paths]
    try:
        merged = merge_bam_reference_headers_with_aliases(handles)
        references, lengths = merged.references, merged.lengths
        selected_specs = expand_contig_tokens(args.contigs, references)
        if any(":" in spec for spec in selected_specs):
            raise ValueError("The fragments command accepts whole-contig selectors, not coordinate intervals")
        selected = selected_specs
        length_by_contig = dict(zip(references, lengths))
        for index, contig in enumerate(selected, start=1):
            reporter.contig("Loading fragments", contig, index, len(selected))
            shared_counts = defaultdict(int) if args.dedup_scope == "all_bams" else None
            rows: list[tuple[int, int]] = []
            for handle, source_mapping in zip(handles, merged.source_contigs):
                source_contig = source_mapping.get(contig)
                if source_contig is None:
                    continue
                counts = shared_counts if args.dedup_scope == "all_bams" else None
                rows.extend(
                    generate_fragment_ranges(
                        bamfile=handle,
                        contig=source_contig,
                        fetch_start=0,
                        fetch_end=int(length_by_contig[contig]),
                        max_per_coordinate=args.max_duplicates,
                        subsample=args.subsample,
                        fragment_counts=counts,
                        min_mapq=args.min_mapq,
                        require_proper_pair=not args.allow_improper_pairs,
                        include_marked_duplicates=args.include_marked_duplicates,
                        include_qcfail=args.include_qcfail,
                        allow_softclipped=args.allow_softclipped,
                        include_secondary=args.include_secondary,
                        include_supplementary=args.include_supplementary,
                        dedup_contig=contig,
                    )
                )
            rows = [
                (start, end)
                for start, end in rows
                if args.frag_lower <= end - start <= args.frag_upper
            ]
            rows.sort()
            yield contig, int(length_by_contig[contig]), rows
    finally:
        for handle in handles:
            handle.close()


def _interval_fragments_by_contig(args, paths: list[str]):
    reporter = ProgressReporter("fragments")
    source = IntervalFragmentSource(paths, chrom_sizes=args.chrom_sizes)
    try:
        selected_specs = expand_contig_tokens(args.contigs, source.references)
        if any(":" in spec for spec in selected_specs):
            raise ValueError("The fragments command accepts whole-contig selectors, not coordinate intervals")
        lengths = dict(zip(source.references, source.lengths))
        for index, contig in enumerate(selected_specs, start=1):
            reporter.contig(
                "Loading fragments", contig, index, len(selected_specs)
            )
            rows = source.fetch(
                contig,
                0,
                lengths[contig],
                max_per_coordinate=args.max_duplicates,
                subsample=args.subsample,
                dedup_scope=args.dedup_scope,
            )
            rows = [
                (start, end)
                for start, end in rows
                if args.frag_lower <= end - start <= args.frag_upper
            ]
            rows.sort()
            yield contig, lengths[contig], rows
    finally:
        source.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite fragments",
        description=(
            "Extract fragment BED3 coordinates from paired-end BAM files, or combine, "
            "filter and deduplicate existing BED/BED.gz/bigBed fragment files."
        ),
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--bam", "--bamfiles", dest="bamfiles", nargs="+", help="Coordinate-sorted paired-end BAM input(s).")
    group.add_argument("--fragments", "--fragment-bed", dest="fragment_files", nargs="+", help="BED, BED.gz, bigBed, or .bb fragment interval input(s).")
    parser.add_argument("-o", "--output-prefix", help="Path prefix for fragment intervals and summary outputs. Default: primary input basename.")
    parser.add_argument("-c", "--contigs", nargs="+", default=["all"], help="Contigs to retain; supports names, comma lists, numeric ranges, autosomes, and all (default: all).")
    parser.add_argument("--frag-lower", type=int, default=1, help="Inclusive minimum fragment length in bp (default: 1).")
    parser.add_argument("--frag-upper", type=int, default=1000, help="Inclusive maximum fragment length in bp (default: 1000).")
    parser.add_argument("--max-duplicates", dest="max_duplicates", type=int, default=0, help="Maximum copies of an identical complete fragment retained; 0 permits all copies (default: 0).")
    parser.add_argument("--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams", help="Apply coordinate limits across all inputs or within each input (default: all_bams).")
    parser.add_argument("--subsample", type=float, help="Optional independent fragment-retention probability from 0 to 1.")
    parser.add_argument("--seed", type=int, help="Random seed used for reproducible subsampling.")
    parser.add_argument("--split-contigs", action="store_true", help="Write one fragment interval file per selected contig.")
    parser.add_argument("--output-format", choices=("bed", "bed.gz", "bigbed", "both"), default="bed.gz", help="Fragment interval format (default: bed.gz).")
    parser.add_argument("--chrom-sizes", help="Chromosome-size table, BAM or CRAM for BED input and bigBed output.")
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality for BAM reads (default: 0).")
    parser.add_argument("--allow-improper-pairs", action="store_true", help="Accept paired reads without the proper-pair SAM flag.")
    parser.add_argument("--include-marked-duplicates", action="store_true", help="Accept reads carrying the duplicate SAM flag.")
    parser.add_argument("--include-qcfail", action="store_true", help="Accept reads carrying the QC-fail SAM flag.")
    parser.add_argument("--allow-softclipped", action="store_true", help="Accept read pairs containing soft-clipped alignments.")
    parser.add_argument("--include-secondary", action="store_true", help="Accept secondary alignments.")
    parser.add_argument("--include-supplementary", action="store_true", help="Accept supplementary alignments.")
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser)
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def run(args: argparse.Namespace) -> int:
    parser = build_parser()
    reporter = ProgressReporter("fragments")
    if pysam is None and (args.bamfiles):
        parser.error("pysam is required for this operation")
    if args.frag_lower < 1 or args.frag_upper < args.frag_lower:
        parser.error("Require 1 <= --frag-lower <= --frag-upper")
    if args.max_duplicates < 0:
        parser.error("--max-duplicates must be 0 or greater")
    if args.min_mapq < 0:
        parser.error("--min-mapq must be 0 or greater")
    if args.subsample is not None and not 0 <= args.subsample <= 1:
        parser.error("--subsample must be between 0 and 1")
    if args.seed is not None:
        random.seed(args.seed)

    from nucleosuite.output_naming import automatic_prefix, parameterized_prefix

    if not args.output_prefix:
        primary = (args.bamfiles or args.fragment_files)[0]
        args.output_prefix = str(automatic_prefix(primary))
    args.output_prefix = str(
        parameterized_prefix(
            args.output_prefix,
            (
                ("fragmin", args.frag_lower),
                ("fragmax", args.frag_upper),
                ("mapq", args.min_mapq),
                ("maxdup", args.max_duplicates),
                ("dedup", args.dedup_scope),
                ("subsample", args.subsample),
                ("seed", args.seed),
            ),
        )
    )

    from nucleosuite.parallel import run_native_per_contig
    if not getattr(args, "_per_contig_worker", False) and int(getattr(args, "cores", 1) or 1) > 1:
        return run_native_per_contig("fragments", args, run)

    reporter.stage("Resolving fragment inputs")
    paths = expand_inputs(args.bamfiles or args.fragment_files)
    iterator = (
        _bam_fragments_by_contig(args, paths)
        if args.bamfiles
        else _interval_fragments_by_contig(args, paths)
    )
    prefix = Path(args.output_prefix)
    chrom_sizes: list[tuple[str, int]] = []
    length_counts: Counter[int] = Counter()
    output_paths: list[Path] = []

    combined_bed_path: Path | None = None
    combined_handle = None
    if not args.split_contigs:
        if args.output_format == "bed.gz":
            combined_bed_path = Path(f"{prefix}.fragments.bed.gz")
        else:
            # BigBed conversion requires plain BED. ``both`` intentionally
            # retains this file alongside the resulting .bb file.
            combined_bed_path = Path(f"{prefix}.fragments.bed")
        combined_handle = _open_output(combined_bed_path)

    try:
        for contig, contig_length, rows in iterator:
            chrom_sizes.append((contig, contig_length))
            for start, end in rows:
                length_counts[end - start] += 1
            if args.split_contigs:
                safe = contig.replace("/", "_")
                if args.output_format == "bed.gz":
                    bed_path = Path(f"{prefix}.{safe}.fragments.bed.gz")
                else:
                    bed_path = Path(f"{prefix}.{safe}.fragments.bed")
                _write_records(bed_path, ((contig, start, end) for start, end in rows))

                if args.output_format in {"bed", "bed.gz", "both"}:
                    output_paths.append(bed_path)
                if args.output_format in {"bigbed", "both"}:
                    bb = convert_bed_to_bigbed(bed_path, [(contig, contig_length)])
                    if bb is not None:
                        output_paths.append(bb)
                    elif args.output_format == "bigbed":
                        # Empty interval collections cannot be represented as
                        # bigBed; retain the empty BED as an explicit output.
                        output_paths.append(bed_path)
                    if args.output_format == "bigbed" and bb is not None:
                        bed_path.unlink(missing_ok=True)
            else:
                assert combined_handle is not None
                for start, end in rows:
                    combined_handle.write(f"{contig}\t{start}\t{end}\n")
    finally:
        if combined_handle is not None:
            combined_handle.close()

    if not args.split_contigs:
        assert combined_bed_path is not None
        if args.output_format in {"bed", "bed.gz", "both"}:
            output_paths.append(combined_bed_path)
        if args.output_format in {"bigbed", "both"}:
            bb = convert_bed_to_bigbed(combined_bed_path, chrom_sizes)
            if bb is not None:
                output_paths.append(bb)
            elif args.output_format == "bigbed":
                output_paths.append(combined_bed_path)
            if args.output_format == "bigbed" and bb is not None:
                combined_bed_path.unlink(missing_ok=True)

    summary_path = Path(f"{prefix}.fragments.summary.tsv")
    reporter.stage("Writing fragment intervals, summary, and length distribution")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("wt", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"input_type\t{'bam' if args.bamfiles else 'fragments'}\n")
        handle.write(f"input_files\t{len(paths)}\n")
        handle.write(f"fragments_written\t{sum(length_counts.values())}\n")
        handle.write(f"contigs_written\t{len(chrom_sizes)}\n")
        handle.write(f"fragment_length_min\t{args.frag_lower}\n")
        handle.write(f"fragment_length_max\t{args.frag_upper}\n")
        handle.write(f"max_per_coordinate\t{args.max_duplicates}\n")
        handle.write(f"dedup_scope\t{args.dedup_scope}\n")
    counts_path = Path(f"{prefix}.fragment_length_counts.tsv")
    with counts_path.open("wt", encoding="utf-8") as handle:
        handle.write("fragment_length\tcount\n")
        for length in sorted(length_counts):
            handle.write(f"{length}\t{length_counts[length]}\n")
    from nucleosuite.plotting import plot_path
    counts_plot = plot_path(Path(f"{prefix}.fragment_length_distribution.png"))
    saved_counts_plot = plot_count_profile(
        str(counts_path),
        str(counts_plot),
        x_column="fragment_length",
        y_column="count",
        xlabel="Fragment length (bp)",
        ylabel="Fragment count",
        title="Fragment-length distribution",
    )
    if saved_counts_plot is not None:
        counts_plot = Path(saved_counts_plot)

    for path in output_paths:
        print(f"Wrote: {path}")
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {counts_path}")
    if counts_plot.exists():
        print(f"Wrote: {counts_plot}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
