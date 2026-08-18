"""CLI for ``nucleosuite dyads`` (alias: ``dyad``)."""

from nucleosuite.cli.basic_common import validate
from nucleosuite.cli.common import add_bam_fragment_arguments, add_randomization_arguments


def register(subparsers):
    parser = subparsers.add_parser(
        "dyads",
        aliases=["dyad"],
        help="Create a track of fragment centres, often used as nucleosome centres.",
    )
    add_bam_fragment_arguments(
        parser,
        default_lower=1,
        default_upper=1000,
        default_max_duplicates=1,
        default_even_dyad="split",
    )
    add_randomization_arguments(parser)
    parser.add_argument(
        "--output-format",
        choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="Dyad-track format: BigWig, compressed WIG, both, or no signal file (default: bigwig).",
    )
    parser.set_defaults(command_runner=run)


def run(args):
    validate(args)
    from nucleosuite.workflows import basic_tracks as workflow
    from nucleosuite.parallel import run_native_per_contig

    return run_native_per_contig(
        "dyads", args, lambda namespace: workflow.run(namespace, "dyads", ["dyad"])
    )
