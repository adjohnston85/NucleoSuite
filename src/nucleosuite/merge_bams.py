#!/usr/bin/env python3
"""Merge, filter, coordinate-deduplicate and optionally split paired-end BAMs."""

from __future__ import annotations

import argparse
import json
import glob
import os
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence

try:
    import pysam
except ImportError:  # allow --help without runtime dependency
    pysam = None

from nucleosuite.core.fragments import generate_paired_reads, require_bam_indexes
from nucleosuite.core.regions import expand_contig_tokens
from nucleosuite.progress import ProgressReporter


def expand_inputs(values: Sequence[str]) -> list[str]:
    paths: list[str] = []
    for value in values:
        matches = sorted(glob.glob(value))
        paths.extend(matches or [value])
    paths = list(dict.fromkeys(paths))
    missing = [path for path in paths if not Path(path).exists()]
    if missing:
        raise FileNotFoundError("BAM file(s) not found: " + ", ".join(missing))
    return paths


def union_header(handles: list[pysam.AlignmentFile], selected_contigs: list[str] | None = None) -> dict:
    first = handles[0].header.to_dict()
    output: dict = {"HD": dict(first.get("HD", {}))}
    output["HD"]["SO"] = "coordinate"
    lengths: dict[str, int] = {}
    order: list[str] = []
    for handle in handles:
        for name, length in zip(handle.references, handle.lengths):
            length = int(length)
            if name not in lengths:
                lengths[name] = length
                order.append(name)
            elif lengths[name] != length:
                raise ValueError(f"Conflicting BAM reference length for {name}: {lengths[name]} vs {length}")
    if selected_contigs is not None:
        order = [name for name in order if name in selected_contigs]
    output["SQ"] = [{"SN": name, "LN": lengths[name]} for name in order]

    for key, identity in (("RG", "ID"), ("PG", "ID")):
        records = []
        seen = set()
        for handle in handles:
            for record in handle.header.to_dict().get(key, []):
                ident = record.get(identity)
                if ident in seen:
                    continue
                seen.add(ident)
                records.append(dict(record))
        if records:
            output[key] = records
    comments = []
    for handle in handles:
        comments.extend(handle.header.to_dict().get("CO", []))
    if comments:
        output["CO"] = list(dict.fromkeys(comments))
    return output


def _write_output(
    *,
    handles: list[pysam.AlignmentFile],
    contigs: list[str],
    output_path: Path,
    args,
    per_contig: bool,
    progress: ProgressReporter | None = None,
) -> tuple[int, int]:
    header_dict = union_header(handles, contigs if per_contig else None)
    output_header = pysam.AlignmentHeader.from_dict(header_dict)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_fd, temp_name = tempfile.mkstemp(prefix="nucleosuite_merge_", suffix=".bam", dir=str(output_path.parent))
    os.close(temp_fd)
    pairs_written = 0
    reads_written = 0
    try:
        with pysam.AlignmentFile(temp_name, "wb", header=output_header) as writer:
            for index, contig in enumerate(contigs, start=1):
                if progress is not None:
                    progress.contig("Merging BAM records", contig, index, len(contigs))
                shared_counts = defaultdict(int) if args.dedup_scope == "all_bams" else None
                for handle in handles:
                    if contig not in handle.references:
                        continue
                    counts = shared_counts if args.dedup_scope == "all_bams" else None
                    for forward, reverse in generate_paired_reads(
                        bamfile=handle,
                        contig=contig,
                        start=0,
                        end=handle.get_reference_length(contig),
                        max_per_coordinate=args.max_duplicates,
                        subsample=args.subsample,
                        fragment_counts=counts,
                        min_mapq=args.min_mapq,
                        require_proper_pair=not args.allow_improper_pairs,
                        include_marked_duplicates=not args.exclude_marked_duplicates,
                        include_qcfail=args.include_qcfail,
                        allow_softclipped=args.allow_softclipped,
                        include_secondary=args.include_secondary,
                        include_supplementary=args.include_supplementary,
                    ):
                        fragment_start = min(forward.reference_start, reverse.reference_start)
                        fragment_end = max(forward.reference_end or 0, reverse.reference_end or 0)
                        length = fragment_end - fragment_start
                        if length < args.frag_lower or length > args.frag_upper:
                            continue
                        for read in (forward, reverse):
                            converted = pysam.AlignedSegment.fromstring(read.to_string(), output_header)
                            writer.write(converted)
                            reads_written += 1
                        pairs_written += 1
        pysam.sort("-o", str(output_path), temp_name)
        pysam.index(str(output_path))
    finally:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
    return pairs_written, reads_written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite merge-bams",
        description=(
            "Combine paired-end BAM files, with optional fragment-coordinate "
            "deduplication, filtering and one-output-per-contig splitting."
        ),
    )
    parser.add_argument("-b", "--bam", "--bamfiles", dest="bamfiles", nargs="+", required=True, help="Coordinate-sorted, indexed paired-end BAM input(s).")
    parser.add_argument("-o", "--output", help="Merged BAM path when not using --split-contigs.")
    parser.add_argument("--output-prefix", help="Output prefix for --split-contigs; default is derived from the first BAM basename.")
    parser.add_argument("-c", "--contigs", nargs="+", default=["all"], help="Whole contigs to retain; supports names, comma lists, numeric ranges, autosomes, and all (default: all).")
    parser.add_argument("--split-contigs", action="store_true", help="Write one sorted and indexed BAM per selected contig.")
    parser.add_argument("--frag-lower", type=int, default=1, help="Inclusive minimum template length in bp (default: 1).")
    parser.add_argument("--frag-upper", type=int, default=1_000_000, help="Inclusive maximum template length in bp (default: 1000000).")
    parser.add_argument("--max-duplicates", dest="max_duplicates", type=int, default=0, help="Maximum read pairs retained per identical fragment coordinate; 0 permits all (default: 0).")
    parser.add_argument("--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams", help="Apply the coordinate limit across all BAMs or within each BAM (default: all_bams).")
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality (default: 0).")
    parser.add_argument("--allow-improper-pairs", action="store_true", help="Accept paired reads without the proper-pair SAM flag.")
    parser.add_argument("--exclude-marked-duplicates", action="store_true", help="Exclude reads carrying the duplicate SAM flag.")
    parser.add_argument("--include-qcfail", action="store_true", help="Accept reads carrying the QC-fail SAM flag.")
    parser.add_argument("--allow-softclipped", action="store_true", help="Accept read pairs containing soft-clipped alignments.")
    parser.add_argument("--include-secondary", action="store_true", help="Accept secondary alignments.")
    parser.add_argument("--include-supplementary", action="store_true", help="Accept supplementary alignments.")
    parser.add_argument("--subsample", type=float, help="Optional independent read-pair retention probability from 0 to 1.")
    parser.add_argument("--seed", type=int, help="Random seed used for reproducible subsampling.")
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser)
    return parser


def _run_serial(args: argparse.Namespace, parser: argparse.ArgumentParser | None = None) -> int:
    parser = parser or build_parser()
    reporter = ProgressReporter("merge-bams")
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
    if args.split_contigs and not args.output_prefix:
        from nucleosuite.output_naming import automatic_prefix
        args.output_prefix = str(automatic_prefix(args.bamfiles[0], "merged"))
    if not args.split_contigs and not args.output:
        from nucleosuite.output_naming import automatic_prefix
        args.output = str(automatic_prefix(args.bamfiles[0], "merged")) + ".bam"
    if args.seed is not None:
        import random
        random.seed(args.seed)

    reporter.stage("Resolving and validating BAM inputs")
    paths = expand_inputs(args.bamfiles)
    require_bam_indexes(paths)
    handles = [pysam.AlignmentFile(path, "rb") for path in paths]
    summary = Counter()
    outputs: list[Path] = []
    try:
        references = []
        for handle in handles:
            for reference in handle.references:
                if reference not in references:
                    references.append(reference)
        selected = expand_contig_tokens(args.contigs, references)
        if any(":" in spec for spec in selected):
            parser.error("merge-bams accepts whole-contig selectors, not coordinate intervals")
        if args.split_contigs:
            for contig in selected:
                output = Path(f"{args.output_prefix}.{contig.replace('/', '_')}.bam")
                pairs, reads = _write_output(
                    handles=handles,
                    contigs=[contig],
                    output_path=output,
                    args=args,
                    per_contig=True,
                    progress=reporter,
                )
                summary["pairs"] += pairs
                summary["reads"] += reads
                outputs.append(output)
        else:
            output = Path(args.output)
            pairs, reads = _write_output(
                handles=handles,
                contigs=selected,
                output_path=output,
                args=args,
                per_contig=False,
                progress=reporter,
            )
            summary["pairs"] = pairs
            summary["reads"] = reads
            outputs.append(output)
    finally:
        for handle in handles:
            handle.close()

    prefix = Path(args.output_prefix) if args.split_contigs else Path(args.output).with_suffix("")
    qc = Path(f"{prefix}.merge_summary.tsv")
    reporter.stage("Indexing outputs and writing merge summary")
    qc.parent.mkdir(parents=True, exist_ok=True)
    with qc.open("wt", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"input_bams\t{len(paths)}\n")
        handle.write(f"output_bams\t{len(outputs)}\n")
        handle.write(f"read_pairs_written\t{summary['pairs']}\n")
        handle.write(f"reads_written\t{summary['reads']}\n")
        handle.write(f"max_per_coordinate\t{args.max_duplicates}\n")
        handle.write(f"dedup_scope\t{args.dedup_scope}\n")
    for output in outputs:
        print(f"Wrote: {output}")
        print(f"Wrote: {output}.bai")
    print(f"Wrote: {qc}")
    return 0


def _merge_contig_worker(paths: list[str], contig: str, output_path: str, args_data: dict) -> tuple[int, int, str]:
    local_args = argparse.Namespace(**args_data)
    handles = [pysam.AlignmentFile(path, "rb") for path in paths]
    try:
        pairs, reads = _write_output(
            handles=handles,
            contigs=[contig],
            output_path=Path(output_path),
            args=local_args,
            per_contig=True,
        )
    finally:
        for handle in handles:
            handle.close()
    prefix = str(Path(output_path).with_suffix(""))
    qc = Path(prefix + ".merge_summary.tsv")
    with qc.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"input_bams\t{len(paths)}\n")
        handle.write("output_bams\t1\n")
        handle.write(f"read_pairs_written\t{pairs}\n")
        handle.write(f"reads_written\t{reads}\n")
        handle.write(f"max_per_coordinate\t{local_args.max_duplicates}\n")
        handle.write(f"dedup_scope\t{local_args.dedup_scope}\n")
    return pairs, reads, prefix


def _run_parallel(args: argparse.Namespace) -> int:
    if pysam is None:
        raise RuntimeError("pysam is required for this operation")
    paths = expand_inputs(args.bamfiles)
    require_bam_indexes(paths)
    handles = [pysam.AlignmentFile(path, "rb") for path in paths]
    try:
        references: list[str] = []
        lengths: dict[str, int] = {}
        for handle in handles:
            for reference, length in zip(handle.references, handle.lengths):
                if reference not in references:
                    references.append(reference)
                    lengths[reference] = int(length)
        selected = expand_contig_tokens(args.contigs, references)
    finally:
        for handle in handles:
            handle.close()
    if any(":" in spec for spec in selected):
        raise ValueError("merge-bams accepts whole-contig selectors, not coordinate intervals")
    if len(selected) <= 1:
        return _run_serial(args)

    from nucleosuite.parallel import MANIFEST_NAME, _safe_contig
    requested = Path(args.output or (str(args.output_prefix) + ".bam"))
    base_name = requested.stem
    root = (
        Path(args.parallel_dir)
        if args.parallel_dir
        else requested.parent / f"{base_name}_multicontig"
    ).resolve()
    per_root = root / "per_contig"
    combined_dir = root / "combined"
    per_root.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    args_data = {
        key: value for key, value in vars(args).items()
        if key not in {"command", "command_function", "command_runner"}
    }
    entries: list[dict[str, object]] = []
    failures: list[str] = []
    from nucleosuite.parallel import _worker_initializer
    with ProcessPoolExecutor(
        max_workers=min(int(args.cores), len(selected)),
        initializer=_worker_initializer,
    ) as executor:
        futures = {}
        for contig in selected:
            safe = _safe_contig(contig)
            directory = per_root / safe
            directory.mkdir(parents=True, exist_ok=True)
            output = directory / f"{base_name}_{safe}.bam"
            future = executor.submit(_merge_contig_worker, paths, contig, str(output), args_data)
            futures[future] = (contig, output.with_suffix(""))
        for future in as_completed(futures):
            contig, expected = futures[future]
            try:
                _pairs, _reads, prefix = future.result()
                entries.append({"contig": contig, "prefix": prefix, "exit_code": 0})
                print(f"Completed {contig}")
            except Exception as error:
                failures.append(f"{contig}: {error}")
                entries.append({"contig": contig, "prefix": str(expected), "exit_code": 2, "error": str(error)})
    order = {contig: index for index, contig in enumerate(selected)}
    entries.sort(key=lambda item: order[str(item["contig"])])
    manifest = {
        "schema_version": 1,
        "command": "merge-bams",
        "base_name": base_name,
        "combined_name": base_name,
        "root_dir": str(root),
        "per_contig_dir": str(per_root),
        "combined_dir": str(combined_dir),
        "chrom_sizes": [{"chrom": chrom, "size": lengths[chrom]} for chrom in selected],
        "per_contig": entries,
        "options": {"cores": int(args.cores)},
    }
    with (root / MANIFEST_NAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, default=str)
        handle.write("\n")
    if failures:
        raise RuntimeError("Per-contig BAM jobs failed: " + "; ".join(failures))
    if not args.skip_combine:
        from nucleosuite.combine import combine_run
        result = combine_run(root, combine_tracks=False, cores=int(args.cores))
        print(f"Combined outputs: {result['output_dir']}")
        for warning in result["warnings"]:
            print(f"WARNING: {warning}")
    else:
        print(f"Per-contig outputs: {per_root}")
    return 0


def run(args: argparse.Namespace) -> int:
    cores = int(getattr(args, "cores", 1) or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    if cores == 1 or getattr(args, "_per_contig_worker", False):
        return _run_serial(args)
    return _run_parallel(args)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
