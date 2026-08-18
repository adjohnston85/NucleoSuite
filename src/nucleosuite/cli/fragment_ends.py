"""CLI for ``nucleosuite fragment-ends``."""

from nucleosuite.cli.basic_common import register_basic_parser, validate


def register(subparsers):
    parser = register_basic_parser(
        subparsers,
        "fragment-ends",
        "Create tracks showing where fragments begin and end.",
    )
    parser.add_argument(
        "--tracks",
        nargs="+",
        choices=("combined", "left", "right"),
        default=("combined", "left", "right"),
        help="Boundary tracks to write: both ends combined, genomic-left ends, and/or genomic-right ends (default: all three).",
    )
    parser.set_defaults(command_runner=run)


def run(args):
    validate(args)
    mapping = {
        "combined": "fragment_ends",
        "left": "fragment_left_ends",
        "right": "fragment_right_ends",
    }
    tracks = [mapping[value] for value in dict.fromkeys(args.tracks)]
    from nucleosuite.workflows import basic_tracks as workflow
    from nucleosuite.parallel import run_native_per_contig

    return run_native_per_contig(
        "fragment-ends", args, lambda namespace: workflow.run(namespace, "fragment_ends", tracks)
    )
