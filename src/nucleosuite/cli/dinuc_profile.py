"""CLI registration for dinucleotide profiles."""

from __future__ import annotations

from nucleosuite.cli.common import (
    add_bam_fragment_arguments,
    add_randomization_arguments,
    validate_bam_arguments,
)


def register(subparsers):
    parser = subparsers.add_parser(
        "dinuc-profile",
        help="Measure two-base sequence patterns around fragment centres.",
    )
    add_bam_fragment_arguments(
        parser, default_lower=137, default_upper=197,
        default_max_duplicates=1, include_even_dyad=False,
    )
    add_randomization_arguments(parser)
    parser.add_argument(
        "--dinuc-fraction", action="store_true",
        help="Write fractions from 0 to 1 instead of percentages.",
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    parser.set_defaults(command_runner=run)


def run(args):
    validate_bam_arguments(args)
    if args.frag_lower < 2:
        raise ValueError("dinuc-profile requires --frag-lower >= 2")
    if not args.fasta:
        raise ValueError("dinuc-profile requires --fasta")
    from nucleosuite.workflows import dinuc_profile as workflow
    from nucleosuite.parallel import run_native_per_contig
    return run_native_per_contig("dinuc-profile", args, workflow.run)
