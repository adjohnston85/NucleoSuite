"""Generate chromosome-size tables from BAM or CRAM headers."""

from __future__ import annotations

import argparse
from pathlib import Path
from collections.abc import Sequence

from nucleosuite.core.chrom_sizes import (
    filter_chrom_sizes,
    read_chrom_sizes_source,
    write_chrom_sizes_table,
)
from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.progress import ProgressReporter


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite chrom-sizes",
        description="Write chromosome names and lengths from a BAM or CRAM header.",
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--bam", "--alignment", dest="alignment", required=True, help="BAM or CRAM whose reference header supplies names, lengths, and order.")
    parser.add_argument("--output", "-o", help="Destination two-column chromosome-size table. Default: <alignment-basename>.chrom.sizes.")
    parser.add_argument(
        "--contigs",
        nargs="+",
        default=None,
        help="Optional contig selectors such as chr1-22,chrX, autosomes, or all.",
    )
    parser.add_argument(
        "--fasta",
        default=None,
        help="Matching reference FASTA when required to open a CRAM file.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.output:
        from nucleosuite.output_naming import automatic_prefix
        args.output = str(automatic_prefix(args.alignment)) + ".chrom.sizes"
    reporter = ProgressReporter("chrom-sizes")
    reporter.stage(f"Reading reference header: {args.alignment}")
    rows = read_chrom_sizes_source(args.alignment, reference_fasta=args.fasta)
    rows = filter_chrom_sizes(rows, args.contigs)
    reporter.stage(f"Writing {len(rows):,} selected contig sizes")
    output = write_chrom_sizes_table(rows, args.output)
    print(f"Wrote: {output}")
    return 0
