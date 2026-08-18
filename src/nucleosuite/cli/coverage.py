"""CLI for ``nucleosuite coverage``."""

from nucleosuite.cli.basic_common import register_basic_parser, validate


def register(subparsers):
    parser = register_basic_parser(
        subparsers,
        "coverage",
        "Create a per-base track showing how many fragments cover each position.",
    )
    parser.set_defaults(command_runner=run)


def run(args):
    validate(args)
    from nucleosuite.workflows import basic_tracks as workflow
    from nucleosuite.parallel import run_native_per_contig

    return run_native_per_contig(
        "coverage", args, lambda namespace: workflow.run(namespace, "coverage", ["coverage"])
    )
