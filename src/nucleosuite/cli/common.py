"""Argument helpers shared by fragment-coordinate commands."""

from __future__ import annotations

import argparse


def auto_or_integer_mode(text: str) -> str | int:
    """Parse a protected-DNA mode supplied as ``auto`` or an integer."""

    value = str(text).strip().lower()
    if value == "auto":
        return "auto"
    try:
        integer = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("mode must be auto or an integer") from exc
    if integer < 3:
        raise argparse.ArgumentTypeError("mode must be at least 3")
    return integer


def add_mode_estimation_arguments(
    parser: argparse.ArgumentParser,
    *,
    default_search_lower: int | None = None,
    default_search_upper: int | None = None,
) -> None:
    """Add options shared by commands that resolve an automatic fragment mode."""

    parser.add_argument(
        "--mode-search-lower",
        type=int,
        default=default_search_lower,
        help=(
            "Lower fragment length considered during automatic mode estimation. "
            + (
                f"Default: {default_search_lower}."
                if default_search_lower is not None
                else "Default: --frag-lower."
            )
        ),
    )
    parser.add_argument(
        "--mode-search-upper",
        type=int,
        default=default_search_upper,
        help=(
            "Upper fragment length considered during automatic mode estimation. "
            + (
                f"Default: {default_search_upper}."
                if default_search_upper is not None
                else "Default: --frag-upper."
            )
        ),
    )
    parser.add_argument("--mode-min-fragments", type=int, default=100_000, help="Minimum accepted fragments before convergence checks (default: 100000).")
    parser.add_argument("--mode-batch-fragments", type=int, default=25_000, help="Additional accepted fragments between convergence checks (default: 25000).")
    parser.add_argument("--mode-max-fragments", type=int, default=1_000_000, help="Maximum accepted fragments sampled for mode estimation (default: 1000000).")
    parser.add_argument("--mode-bootstrap", type=int, default=200, help="Bootstrap mode replicates per checkpoint (default: 200).")
    parser.add_argument("--mode-stable-checkpoints", type=int, default=3, help="Consecutive stable checkpoints required for early stopping (default: 3).")
    parser.add_argument("--mode-max-change", type=int, default=1, help="Largest mode range allowed across stable checkpoints (default: 1 bp).")
    parser.add_argument("--mode-max-ci-width", type=float, default=4.0, help="Largest bootstrap 95%% interval width allowed for convergence (default: 4 bp).")
    parser.add_argument("--mode-block-bp", type=int, default=1_000_000, help="Genomic block size used for seeded random sampling (default: 1000000 bp).")


def resolve_fragment_mode(args, value: str | int, *, command: str):
    """Resolve an explicit or automatically estimated mode for a fragment command."""

    from nucleosuite.mode_estimation import (
        estimate_bam_fragment_mode,
        estimate_fragment_file_mode,
        mode_estimate_message,
    )

    seed = 12345 if getattr(args, "seed", None) is None else int(args.seed)
    if value != "auto":
        mode = int(value)
        print(f"[{command}] Fragment mode: {mode} bp (explicit)", flush=True)
        return mode, None, "explicit", seed

    search_lower = (
        int(args.frag_lower)
        if args.mode_search_lower is None
        else int(args.mode_search_lower)
    )
    search_upper = (
        int(args.frag_upper)
        if args.mode_search_upper is None
        else int(args.mode_search_upper)
    )
    if search_lower < args.frag_lower or search_upper > args.frag_upper:
        raise ValueError("Automatic mode-search bounds must lie within the fragment range")
    if search_upper < search_lower:
        raise ValueError("--mode-search-upper must be at least --mode-search-lower")
    if min(
        args.mode_min_fragments,
        args.mode_batch_fragments,
        args.mode_max_fragments,
        args.mode_bootstrap,
        args.mode_stable_checkpoints,
        args.mode_block_bp,
    ) < 1:
        raise ValueError("Automatic mode-sampling counts must be positive")
    if args.mode_max_fragments < args.mode_min_fragments:
        raise ValueError("--mode-max-fragments must be at least --mode-min-fragments")
    if args.mode_max_change < 0 or args.mode_max_ci_width < 0:
        raise ValueError("Automatic mode-stability limits must be non-negative")

    print(
        f"[{command}] Estimating fragment mode from raw integer counts",
        flush=True,
    )
    common = dict(
        frag_lower=args.frag_lower,
        frag_upper=args.frag_upper,
        search_lower=search_lower,
        search_upper=search_upper,
        minimum_fragments=args.mode_min_fragments,
        batch_fragments=args.mode_batch_fragments,
        maximum_fragments=args.mode_max_fragments,
        stable_checkpoints=args.mode_stable_checkpoints,
        maximum_mode_change=args.mode_max_change,
        maximum_ci_width=args.mode_max_ci_width,
        bootstrap_replicates=args.mode_bootstrap,
        block_bp=args.mode_block_bp,
        seed=seed,
        blacklist_bed=getattr(args, "blacklist_bed", None),
        max_duplicates=args.max_duplicates,
        dedup_scope=args.dedup_scope,
        contig_tokens=args.contigs,
    )
    if getattr(args, "bamfiles", None):
        estimate = estimate_bam_fragment_mode(args.bamfiles, **common)
    else:
        estimate = estimate_fragment_file_mode(
            args.fragment_files,
            chrom_sizes=getattr(args, "chrom_sizes", None),
            fasta_path=getattr(args, "fasta", None),
            **common,
        )
    print(f"[{command}] {mode_estimate_message('Fragment mode estimate', estimate)}", flush=True)
    return estimate.mode, estimate, "automatic", seed


def add_fragment_input_arguments(
    parser: argparse.ArgumentParser,
    *,
    default_lower: int | None,
    default_upper: int | None,
    default_max_duplicates: int = 1,
    default_even_dyad: str = "split",
    include_even_dyad: bool = True,
) -> None:
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-b", "--bamfiles", "--bam", dest="bamfiles", nargs="+",
        help="Input coordinate-sorted paired-end BAM file(s).",
    )
    inputs.add_argument(
        "--fragments", "--fragment-bed", dest="fragment_files", nargs="+",
        help=(
            "Input fragment interval file(s): BED, BED.gz, .bb or .bigBed. "
            "Only chromosome, start and end (the first three columns) are required."
        ),
    )
    parser.add_argument(
        "--chrom-sizes",
        help=(
            "Optional chromosome-size table, BAM or CRAM for fragment BED input. For plain BED "
            "without this option or --fasta, each contig length is inferred from "
            "the largest fragment end."
        ),
    )
    parser.add_argument(
        "--blacklist-bed",
        help=(
            "Optional BED of excluded regions. Complete fragments overlapping any "
            "blacklisted base are removed before PNS, WPS, coverage, dyad, "
            "fragment-end and sequence-derived tracks are calculated."
        ),
    )
    parser.add_argument(
        "-o", "--out-prefix", dest="out_prefix", default=None,
        help="Output prefix. Default: derived from input names and selected contigs.",
    )
    parser.add_argument(
        "-c", "--contigs", nargs="+", default=None,
        help=(
            "Contigs or intervals. Supports comma lists, numeric ranges, "
            "autosomes, all, and forms such as chr2:100000-200000."
        ),
    )
    lower_help = (
        f"Inclusive minimum accepted fragment length in bp (default: {default_lower})."
        if default_lower is not None
        else "Inclusive minimum accepted fragment length in bp; when omitted the command resolves its automatic lower bound."
    )
    upper_help = (
        f"Inclusive maximum accepted fragment length in bp (default: {default_upper})."
        if default_upper is not None
        else "Inclusive maximum accepted fragment length in bp; when omitted the command resolves its automatic upper bound."
    )
    parser.add_argument(
        "--frag-lower",
        type=int,
        default=default_lower,
        help=lower_help,
    )
    parser.add_argument(
        "--frag-upper",
        type=int,
        default=default_upper,
        help=upper_help,
    )
    parser.add_argument(
        "--max-duplicates",
        dest="max_duplicates", type=int, default=default_max_duplicates,
        help=(
            "Maximum fragments retained with identical complete coordinates. "
            f"Use 1 for coordinate deduplication and 0 to disable it (default: {default_max_duplicates})."
        ),
    )
    parser.add_argument(
        "--max-per-coordinate",
        dest="max_per_coordinate", type=int, default=0,
        help=(
            "Optional cap on dyad and fragment-end signal at one output "
            "coordinate. Use 0 for no cap (default: 0)."
        ),
    )
    parser.add_argument(
        "--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams",
        help=(
            "Apply the identical-fragment coordinate limit across the combined input "
            "collection (all_bams) or independently within each input (per_bam). "
            "Default: all_bams."
        ),
    )
    parser.add_argument(
        "--chunk-bp", type=int, default=100_000,
        help="Core genomic chunk length in bp (default: 100000).",
    )
    parser.add_argument(
        "--overlap-bp", type=int, default=1_000,
        help="Context added to each side of a chunk before core-owned output is retained (default: 1000).",
    )
    parser.add_argument(
        "--subsample", type=float, default=None,
        help="Optional independent fragment-retention probability from 0 to 1 (default: retain all).",
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed used for reproducible fragment subsampling and randomization.",
    )
    if include_even_dyad:
        parser.add_argument(
            "--even-dyad", choices=("split", "left", "right"),
            default=default_even_dyad,
            help=(
                "Representation of even-length fragment centres in dyad outputs: "
                "split 0.5/0.5 between central bases, or assign one count to the "
                f"left or right central base (default: {default_even_dyad})."
            ),
        )
    from nucleosuite.parallel import add_parallel_arguments
    add_parallel_arguments(parser, combine_resources=True, resumable=True)


# Internal alias used by command modules.
add_bam_fragment_arguments = add_fragment_input_arguments


def add_randomization_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--randomize-mode", choices=("none", "uniform", "dinuc_anchor"),
        default="none",
        help=(
            "Fragment-coordinate mode: observed coordinates (none), uniform relocation "
            "within the selected contig (uniform), or terminal-dinucleotide-matched "
            "relocation (dinuc_anchor). Default: none."
        ),
    )
    parser.add_argument(
        "--fasta", default=None,
        help="Indexed reference FASTA; required by sequence-aware operations.",
    )
    parser.add_argument(
        "--anchor-prob-start", type=float, default=0.5,
        help="Probability of trying the fragment-start dinucleotide before the end dinucleotide (default: 0.5).",
    )
    parser.add_argument(
        "--max-anchor-tries", type=int, default=30, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--randomize-fallback", choices=("uniform", "skip"),
        default="uniform",
        help="Action when dinucleotide-matched placement is unavailable: uniform relocation or skip (default: uniform).",
    )


def validate_bam_arguments(args) -> None:
    if bool(getattr(args, "bamfiles", None)) == bool(getattr(args, "fragment_files", None)):
        raise ValueError("Provide exactly one of --bamfiles/--bam or --fragments")
    if args.max_duplicates < 0:
        raise ValueError("--max-duplicates must be 0 or greater")
    if args.max_per_coordinate < 0:
        raise ValueError("--max-per-coordinate must be 0 or greater")
    if args.frag_lower < 1 or args.frag_upper < args.frag_lower:
        raise ValueError("Require 1 <= --frag-lower <= --frag-upper")
    if args.chunk_bp < 1:
        raise ValueError("--chunk-bp must be positive")
    if args.overlap_bp < 0:
        raise ValueError("--overlap-bp must be non-negative")
    if args.subsample is not None and not 0.0 <= args.subsample <= 1.0:
        raise ValueError("--subsample must be between 0 and 1")
    if not 0.0 <= args.anchor_prob_start <= 1.0:
        raise ValueError("--anchor-prob-start must be between 0 and 1")
    if args.max_anchor_tries < 1:
        raise ValueError("--max-anchor-tries must be positive")
    if args.randomize_mode == "dinuc_anchor" and not args.fasta:
        raise ValueError("--randomize-mode dinuc_anchor requires --fasta")


def normalise_track_list(values, valid_tracks, option_name: str):
    if len(values) == 1 and values[0].lower() == "none":
        return []
    unknown = [value for value in values if value not in valid_tracks]
    if unknown:
        raise ValueError(
            f"Unknown {option_name} values: {unknown}. Valid: {sorted(valid_tracks)}"
        )
    return list(dict.fromkeys(values))


def add_interval_output_arguments(
    parser: argparse.ArgumentParser,
    *,
    default: str = "bed",
    include_chrom_sizes: bool = False,
) -> None:
    parser.add_argument(
        "--interval-format",
        choices=("bed", "bigbed", "both"),
        default=default,
        help=(
            "Interval-file output: BED text, bigBed only, or both. bigBed "
            f"requires the UCSC bedToBigBed executable (default: {default})."
        ),
    )
    if include_chrom_sizes:
        parser.add_argument(
            "--interval-chrom-sizes",
            help="Chromosome sizes required for bigBed output when they cannot be inferred.",
        )
