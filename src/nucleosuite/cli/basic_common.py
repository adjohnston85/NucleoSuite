"""Parser setup for coverage, dyad and fragment-end commands."""

from __future__ import annotations

from nucleosuite.cli.common import (
    add_bam_fragment_arguments,
    add_randomization_arguments,
    validate_bam_arguments,
)


def register_basic_parser(subparsers, name: str, help_text: str, *, include_even_dyad: bool = False):
    parser = subparsers.add_parser(name, help=help_text, description=help_text)
    add_bam_fragment_arguments(
        parser,
        default_lower=1,
        default_upper=1000,
        default_max_duplicates=1,
        default_even_dyad="split",
        include_even_dyad=include_even_dyad,
    )
    add_randomization_arguments(parser)
    parser.add_argument(
        "--output-format",
        choices=("bigwig", "wiggz", "both", "none"),
        default="bigwig",
        help="Signal-track format: BigWig, compressed WIG, both, or no signal file (default: bigwig).",
    )
    return parser


def validate(args) -> None:
    validate_bam_arguments(args)
