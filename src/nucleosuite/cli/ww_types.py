"""CLI registration for WW/SS classification."""

from __future__ import annotations

from nucleosuite.cli.common import (
    add_bam_fragment_arguments,
    add_randomization_arguments,
    add_interval_output_arguments,
    validate_bam_arguments,
)


def register(subparsers):
    parser = subparsers.add_parser(
        "ww-types",
        help="Group fragments by nucleosome-associated WW/SS sequence patterns.",
    )
    add_bam_fragment_arguments(
        parser, default_lower=137, default_upper=197,
        default_max_duplicates=1, default_even_dyad="split",
    )
    add_randomization_arguments(parser)
    add_interval_output_arguments(parser)
    parser.add_argument(
        "--split-beds", action="store_true",
        help="Also write one three-column BED file per type and unclassified group.",
    )
    parser.add_argument(
        "--dyad-tracks", action="store_true",
        help="Write all, type1, type2, type3 and type4 dyad tracks.",
    )
    parser.add_argument(
        "--output-format", choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="Format for optional --dyad-tracks (default: bigwig).",
    )
    parser.add_argument(
        "--no-dinuc-profile", dest="dinuc_profile", action="store_false",
        help="Do not write all-fragment and type-specific dinucleotide profiles.",
    )
    parser.set_defaults(dinuc_profile=True)
    parser.add_argument(
        "--dinuc-fraction", action="store_true",
        help="Report dinucleotide frequencies as fractions from 0 to 1 instead of percentages.",
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    parser.set_defaults(command_runner=run)


def run(args):
    validate_bam_arguments(args)
    if not args.fasta:
        raise ValueError("ww-types requires --fasta")
    if args.dinuc_fraction and not args.dinuc_profile:
        raise ValueError("--dinuc-fraction cannot be used with --no-dinuc-profile")
    if args.dyad_tracks and args.output_format == "none":
        raise ValueError("--dyad-tracks requires a non-'none' --output-format")
    from nucleosuite.workflows import ww_types as workflow
    from nucleosuite.parallel import run_native_per_contig
    return run_native_per_contig("ww-types", args, workflow.run)
