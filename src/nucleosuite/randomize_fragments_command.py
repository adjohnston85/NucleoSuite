#!/usr/bin/env python3
"""Materialise reproducible randomised fragment BED files."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import random
import tempfile
from collections import Counter
from pathlib import Path
from typing import Sequence

try:
    import pysam
except ImportError:  # allow --help without runtime dependency
    pysam = None

from nucleosuite.core.fragment_inputs import open_fragment_source
from nucleosuite.core.blacklist import load_blacklist
from nucleosuite.core.randomization import (
    RandomizationBlock, place_dinucleotide_matched, uniform_randomize_fragment,
)
from nucleosuite.core.regions import expand_contig_tokens
from nucleosuite.core.reference import resolve_fasta_contig
from nucleosuite.profile_plots import plot_count_profile
from nucleosuite.progress import ProgressReporter


def _fallback(fragment, region_start, region_end, rng, fallback):
    """Uniform fallback helper for callers that supply direct region bounds.

    The command workflow uses :class:`RandomizationBlock` so blacklist and
    canonical-sequence constraints are enforced before sampling.
    """
    if fallback == "skip":
        return None
    if fallback != "uniform":
        raise ValueError(f"Unknown fallback mode: {fallback}")
    length = fragment[1] - fragment[0]
    maximum = region_end - length
    if maximum < region_start:
        return None
    number_of_starts = maximum - region_start + 1
    original_start = fragment[0]
    original_is_valid = region_start <= original_start <= maximum
    allowed = number_of_starts - (1 if original_is_valid else 0)
    if allowed <= 0:
        return None
    draw = rng.randrange(allowed)
    candidate = region_start + draw
    if original_is_valid and candidate >= original_start:
        candidate += 1
    return candidate, candidate + length


def _output_path(output_prefix: Path, *, chrom: str | None, compress: bool) -> Path:
    suffix = ".randomized.fragments.bed.gz" if compress else ".randomized.fragments.bed"
    if chrom is None:
        return Path(f"{output_prefix}{suffix}")
    return Path(f"{output_prefix}.{chrom.replace('/', '_')}{suffix}")


def _open_text_output(path: Path, compress: bool):
    path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(path, "wt", encoding="utf-8") if compress else path.open("wt", encoding="utf-8")


def _temporary_output_path(path: Path) -> Path:
    """Create a same-directory temporary path suitable for atomic replacement."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent
    )
    os.close(descriptor)
    return Path(name)


def _validate_randomized_bed(path: Path, *, compress: bool, expected: int) -> int:
    """Validate a completed temporary BED before it becomes user-visible."""
    opener = gzip.open if compress else open
    observed = 0
    previous_by_chrom: dict[str, tuple[int, int]] = {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise RuntimeError(
                    f"Randomized BED validation failed at line {line_number}: expected BED3"
                )
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise RuntimeError(
                    f"Randomized BED validation failed at line {line_number}: invalid coordinates"
                ) from exc
            if start < 0 or end <= start:
                raise RuntimeError(
                    f"Randomized BED validation failed at line {line_number}: require 0 <= start < end"
                )
            previous = previous_by_chrom.get(chrom)
            if previous is not None and (start, end) < previous:
                raise RuntimeError(
                    f"Randomized BED validation failed at line {line_number}: coordinates are unsorted"
                )
            previous_by_chrom[chrom] = (start, end)
            observed += 1
    if observed != expected:
        raise RuntimeError(
            f"Randomized BED validation counted {observed:,} records; expected {expected:,}"
        )
    return observed


def _atomic_text_path(path: Path, writer) -> None:
    temporary = _temporary_output_path(path)
    try:
        writer(temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite randomize-fragments",
        description=(
            "Create one materialised random fragment set for reuse by PNS, WPS and "
            "other analyses. Dinucleotide mode preserves a randomly chosen start- "
            "or end-boundary dinucleotide for each fragment."
        ),
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--bam", "--bamfiles", dest="bamfiles", nargs="+", help="Coordinate-sorted paired-end BAM input(s).")
    group.add_argument("--fragments", "--fragment-bed", dest="fragment_files", nargs="+", help="BED, BED.gz, bigBed, or .bb fragment interval input(s).")
    parser.add_argument("--fasta", required=True, help="Indexed reference FASTA used to validate placement and match terminal dinucleotides.")
    parser.add_argument("-o", "--output-prefix", required=True, help="Path prefix for the randomized fragment set and QC outputs.")
    parser.add_argument("-c", "--contigs", nargs="+", default=["all"], help="Contigs to randomize; supports names, comma lists, ranges, autosomes, and all (default: all).")
    parser.add_argument("--frag-lower", type=int, default=120, help="Inclusive minimum source-fragment length in bp (default: 120).")
    parser.add_argument("--frag-upper", type=int, default=180, help="Inclusive maximum source-fragment length in bp (default: 180).")
    parser.add_argument("--max-duplicates", dest="max_duplicates", type=int, default=1, help="Maximum source fragments retained with identical complete coordinates; 0 permits all (default: 1).")
    parser.add_argument("--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams", help="Apply source-coordinate limits across all inputs or within each input (default: all_bams).")
    parser.add_argument("--subsample", type=float, help="Optional independent source-fragment retention probability from 0 to 1.")
    parser.add_argument("--method", choices=("dinucleotide", "uniform"), default="dinucleotide", help="Relocation method: preserve one terminal dinucleotide or sample uniformly (default: dinucleotide).")
    parser.add_argument("--seed", type=int, default=12345, help="Randomization seed (default: 12345).")
    parser.add_argument("--search-window", type=int, default=100_000, help="Maximum local reference block used for dinucleotide candidate indexing (default: 100000 bp).")
    parser.add_argument("--anchor-prob-start", type=float, default=0.5, help="Probability of trying the source start dinucleotide before the end dinucleotide (default: 0.5).")
    parser.add_argument(
        "--max-anchor-tries", type=int, default=30,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--max-randomized-per-coordinate", type=int, default=0,
        help=(
            "Optional cap on randomized fragments with identical complete coordinates; "
            "0 permits independent sampling with replacement (default)."
        ),
    )
    parser.add_argument(
        "--blacklist-bed",
        help=(
            "Optional BED. Source fragments overlapping it are excluded and randomized "
            "fragments cannot be placed across it."
        ),
    )
    parser.add_argument(
        "--fallback",
        choices=("uniform", "skip"),
        default="uniform",
        help=(
            "Action when dinucleotide matching fails. 'uniform' chooses a different "
            "valid position; 'skip' omits the fragment. Default: uniform."
        ),
    )
    parser.add_argument("--chrom-sizes", help="Chromosome-size table, BAM, or CRAM used to select and validate fragment-input contigs.")
    parser.add_argument("--split-contigs", action="store_true", help="Write one randomized fragment file per selected contig.")
    parser.add_argument("--uncompressed", action="store_true", help="Write BED instead of BED.gz.")
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser)
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def _stable_seed(seed: int, contig: str) -> int:
    digest = hashlib.blake2b(contig.encode("utf-8"), digest_size=4).digest()
    return (int(seed) + int.from_bytes(digest, "big")) % (2**32)


def _run_serial(args: argparse.Namespace, parser: argparse.ArgumentParser | None = None) -> int:
    parser = parser or build_parser()
    if pysam is None:
        parser.error("pysam is required for randomize-fragments")
    if args.frag_lower < 1 or args.frag_upper < args.frag_lower:
        parser.error("Require 1 <= --frag-lower <= --frag-upper")
    if args.max_duplicates < 0:
        parser.error("--max-duplicates must be 0 or greater")
    if args.search_window < args.frag_upper:
        parser.error("--search-window must be at least --frag-upper")
    if not 0 <= args.anchor_prob_start <= 1:
        parser.error("--anchor-prob-start must be between 0 and 1")
    if args.max_anchor_tries < 1:
        parser.error("--max-anchor-tries must be positive")
    if args.max_randomized_per_coordinate < 0:
        parser.error("--max-randomized-per-coordinate must be 0 or greater")

    reporter = ProgressReporter("randomize-fragments")
    reporter.stage("Opening reference FASTA and fragment source")
    fasta = pysam.FastaFile(args.fasta)
    source = open_fragment_source(
        bam_paths=args.bamfiles,
        fragment_paths=args.fragment_files,
        chrom_sizes=args.chrom_sizes,
        fasta=fasta,
        blacklist_path=None,
    )
    blacklist = load_blacklist(args.blacklist_bed, source.references, source.lengths)
    stats = Counter()
    reasons = Counter()
    relocation = Counter()
    unique_randomized_coordinates = 0
    duplicate_randomized_fragments = 0
    maximum_randomized_multiplicity = 0
    total_randomized = 0
    output_prefix = Path(args.output_prefix)
    compress = not args.uncompressed
    paths: list[Path] = []
    temporary_paths: dict[Path, Path] = {}
    expected_by_path: Counter[Path] = Counter()
    combined_handle = None
    try:
        selected = expand_contig_tokens(args.contigs, source.references)
        if any(":" in spec for spec in selected):
            parser.error("randomize-fragments accepts whole-contig selectors, not coordinate intervals")
        lengths = dict(zip(source.references, source.lengths))
        if not args.split_contigs:
            path = _output_path(output_prefix, chrom=None, compress=compress)
            temporary_paths[path] = _temporary_output_path(path)
            combined_handle = _open_text_output(temporary_paths[path], compress)
            paths.append(path)
        for contig in selected:
            reporter.reading_contig("fragments", contig)
            contig_length = lengths[contig]
            try:
                fasta_contig = resolve_fasta_contig(fasta, contig)
            except KeyError as exc:
                parser.error(str(exc))
            fasta_contig_length = fasta.get_reference_length(fasta_contig)
            if fasta_contig_length != contig_length:
                parser.error(
                    f"Contig length mismatch for {contig!r} (fragment source: "
                    f"{contig_length:,} bp) and FASTA contig {fasta_contig!r} "
                    f"({fasta_contig_length:,} bp)."
                )
            contig_seed = (
                int(args.seed)
                if getattr(args, "_per_contig_worker", False)
                else _stable_seed(args.seed, contig)
            )
            rng = random.Random(contig_seed)
            # Fragment-source subsampling uses Python's module RNG. Seed it per
            # contig so results are invariant to serial/parallel routing.
            random.seed(contig_seed)
            if args.split_contigs:
                path = _output_path(output_prefix, chrom=contig, compress=compress)
                temporary_paths[path] = _temporary_output_path(path)
                handle = _open_text_output(temporary_paths[path], compress)
                paths.append(path)
            else:
                handle = combined_handle
            try:
                core_start = 0
                while core_start < contig_length:
                    core_end = min(core_start + args.search_window, contig_length)
                    # The right fetch extension recovers complete source fragments whose
                    # left boundary belongs to the core. It is never indexed as a candidate.
                    source_fetch_end = min(contig_length, core_end + args.frag_upper)
                    fragments = source.fetch(
                        contig, core_start, source_fetch_end,
                        max_per_coordinate=args.max_duplicates,
                        subsample=args.subsample, dedup_scope=args.dedup_scope,
                    )
                    fragments = [
                        fragment for fragment in fragments
                        if core_start <= fragment[0] < core_end
                        and fragment[1] <= source_fetch_end
                        and args.frag_lower <= fragment[1] - fragment[0] <= args.frag_upper
                    ]
                    if blacklist is not None:
                        retained = []
                        for fragment in fragments:
                            if blacklist.overlaps(contig, fragment[0], fragment[1]):
                                stats["fragments_excluded_by_blacklist"] += 1
                            else:
                                retained.append(fragment)
                        fragments = retained
                    randomized_chunk: list[tuple[int, int]] = []
                    block_coordinate_counts: Counter[tuple[int, int]] = Counter()
                    if fragments:
                        original_sequence = fasta.fetch(
                            fasta_contig, core_start, source_fetch_end
                        ).upper()
                        candidate_sequence = original_sequence[: core_end - core_start]
                        block = RandomizationBlock(
                            contig,
                            core_start,
                            core_end,
                            candidate_sequence,
                            blacklist=blacklist,
                            coordinate_counts=block_coordinate_counts,
                            max_per_coordinate=args.max_randomized_per_coordinate,
                        )
                        for fragment in fragments:
                            stats["input"] += 1
                            length = fragment[1] - fragment[0]
                            relative_start = fragment[0] - core_start
                            relative_end = fragment[1] - core_start
                            if not (
                                0 <= relative_start <= len(original_sequence) - 2
                                and 2 <= relative_end <= len(original_sequence)
                            ):
                                randomized = None
                                status = "skipped"
                                selected_anchor = matched_anchor = None
                                reason = "source_boundary_unavailable"
                            elif args.method == "uniform":
                                randomized = uniform_randomize_fragment(fragment, block, rng)
                                status = "uniform" if randomized is not None else "skipped"
                                selected_anchor = matched_anchor = None
                                reason = None if randomized is not None else "no_valid_uniform_candidate"
                            else:
                                randomized, status, selected_anchor, matched_anchor, reason = (
                                    place_dinucleotide_matched(
                                        fragment,
                                        start_dinuc=original_sequence[
                                            relative_start : relative_start + 2
                                        ],
                                        end_dinuc=original_sequence[
                                            relative_end - 2 : relative_end
                                        ],
                                        block=block,
                                        rng=rng,
                                        anchor_prob_start=args.anchor_prob_start,
                                        fallback=args.fallback,
                                    )
                                )
                            if selected_anchor:
                                stats[f"anchor_{selected_anchor}_selected"] += 1
                            if matched_anchor:
                                stats[f"anchor_{matched_anchor}_matched"] += 1
                            if reason:
                                reasons[reason] += 1
                            if randomized is None:
                                stats["skipped"] += 1
                                continue
                            stats[status] += 1
                            relocation[randomized[0] - fragment[0]] += 1
                            randomized_chunk.append(randomized)
                        stats["random_candidates_rejected_by_blacklist"] += (
                            block.blacklist_candidate_rejections
                        )
                        stats["random_candidates_rejected_non_acgt"] += (
                            block.non_acgt_candidate_rejections
                        )
                    randomized_chunk.sort()
                    for start, end in randomized_chunk:
                        handle.write(f"{contig}\t{start}\t{end}\n")
                    handle.flush()
                    expected_by_path[path] += len(randomized_chunk)
                    total_randomized += len(randomized_chunk)
                    multiplicities = tuple(block_coordinate_counts.values())
                    unique_randomized_coordinates += len(multiplicities)
                    duplicate_randomized_fragments += sum(
                        max(0, count - 1) for count in multiplicities
                    )
                    maximum_randomized_multiplicity = max(
                        maximum_randomized_multiplicity,
                        max(multiplicities, default=0),
                    )
                    core_start = core_end
            finally:
                if args.split_contigs:
                    handle.close()
        if combined_handle is not None:
            combined_handle.close()
            combined_handle = None
        if stats["input"] == 0:
            raise RuntimeError("No eligible source fragments remained for randomization")
        if total_randomized == 0:
            raise RuntimeError("Randomization produced no usable fragments")
        reporter.stage("Validating randomized fragments and writing QC outputs")
        for path in paths:
            _validate_randomized_bed(
                temporary_paths[path],
                compress=compress,
                expected=int(expected_by_path[path]),
            )
        qc_path = Path(f"{args.output_prefix}.randomization_qc.tsv")
        qc_path.parent.mkdir(parents=True, exist_ok=True)
        with qc_path.open("wt", encoding="utf-8") as handle:
            handle.write("metric\tvalue\n")
            for key in (
                "input", "matched", "uniform", "fallback", "skipped",
                "fragments_excluded_by_blacklist",
                "random_candidates_rejected_by_blacklist",
                "random_candidates_rejected_non_acgt",
                "anchor_start_selected", "anchor_end_selected",
                "anchor_start_matched", "anchor_end_matched",
            ):
                handle.write(f"{key}\t{stats[key]}\n")
            for reason in sorted(reasons):
                handle.write(f"reason_{reason}\t{reasons[reason]}\n")
            handle.write(f"unique_randomized_coordinates\t{unique_randomized_coordinates}\n")
            handle.write(f"duplicate_randomized_fragments\t{duplicate_randomized_fragments}\n")
            handle.write(
                f"maximum_randomized_multiplicity\t{maximum_randomized_multiplicity}\n"
            )
            collision_fraction = duplicate_randomized_fragments / total_randomized if total_randomized else 0.0
            handle.write(f"collision_fraction\t{collision_fraction:.12g}\n")
            handle.write(f"seed\t{args.seed}\n")
            handle.write(f"seed_derivation\tstable_seed_plus_contig\n")
            handle.write(f"method\t{args.method}\n")
            handle.write(f"search_window\t{args.search_window}\n")
            handle.write(f"candidate_flank_bp\t0\n")
            handle.write(f"fallback_mode\t{args.fallback}\n")
            handle.write(
                f"max_randomized_per_coordinate\t{args.max_randomized_per_coordinate}\n"
            )
            handle.write(f"blacklist_file\t{args.blacklist_bed or ''}\n")
            if blacklist is not None:
                selected_blacklist_intervals = sum(
                    len(blacklist.intervals.get(contig, ())) for contig in selected
                )
                selected_blacklisted_bases = sum(
                    end - start
                    for contig in selected
                    for start, end in blacklist.intervals.get(contig, ())
                )
                handle.write(
                    f"blacklist_intervals\t{selected_blacklist_intervals}\n"
                )
                handle.write(
                    f"blacklisted_bases\t{selected_blacklisted_bases}\n"
                )
            handle.write("original_position_allowed\tfalse\n")
            handle.write("candidate_sequence_requires_acgt\ttrue\n")
        relocation_path = Path(f"{args.output_prefix}.relocation_distances.tsv")
        with relocation_path.open("wt", encoding="utf-8") as handle:
            handle.write("relocation_bp\tcount\n")
            for distance in sorted(relocation):
                handle.write(f"{distance}\t{relocation[distance]}\n")
        from nucleosuite.plotting import plot_path
        relocation_plot = plot_path(Path(f"{args.output_prefix}.relocation_distances.png"))
        saved_relocation_plot = plot_count_profile(
            str(relocation_path), str(relocation_plot),
            x_column="relocation_bp", y_column="count",
            xlabel="Fragment relocation (bp)", ylabel="Fragment count",
            title="Randomized fragment relocation distances", vertical_zero=True,
        )
        if saved_relocation_plot is not None:
            relocation_plot = Path(saved_relocation_plot)
        # The validated randomized BED is the transactional suite input. Make
        # it visible only after every companion QC output has completed.
        for path in paths:
            os.replace(temporary_paths[path], path)
        for path in paths:
            print(f"Wrote: {path}")
        print(f"Wrote: {qc_path}")
        print(f"Wrote: {relocation_path}")
        if relocation_plot.exists():
            print(f"Wrote: {relocation_plot}")
    finally:
        if combined_handle is not None:
            combined_handle.close()
        for temporary in temporary_paths.values():
            temporary.unlink(missing_ok=True)
        source.close()
        fasta.close()
    return 0

def run(args: argparse.Namespace) -> int:
    from nucleosuite.parallel import run_native_per_contig
    return run_native_per_contig(
        "randomize-fragments", args, lambda namespace: _run_serial(namespace)
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
