"""Top-level ``nucleosuite`` command."""

from __future__ import annotations

import argparse
import importlib
import sys
from collections.abc import Callable, Sequence

from nucleosuite import __version__
from nucleosuite.console import message_pair, print_completion, print_startup
from nucleosuite.logo import render_logo
from nucleosuite.progress import ProgressReporter
from nucleosuite.cli import (
    call_peaks,
    coverage,
    dinuc_profile,
    dyads,
    fragment_ends,
    tracks,
    pns,
    wps,
    ww_types,
)

CommandMain = Callable[[Sequence[str] | None], int]


class TopLevelHelpParser(argparse.ArgumentParser):
    """Argument parser that prefixes only top-level help with the logo."""

    def format_help(self) -> str:
        return f"{render_logo()}\n\n{super().format_help()}"


def _aggregate_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.cli.aggregate import main
    return main(argv)





def _chrom_sizes_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.chrom_sizes_command import main
    return main(argv)

def _combine_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.combine import main
    return main(argv)

def _compare_positions_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.compare_positions import main
    return main(argv)


def _dac_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.dac import main
    return main(argv)


def _dcc_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.dcc import main
    return main(argv)


def _distances_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.distances import main
    return main(argv)


def _filter_coverage_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.filter_coverage import main
    return main(argv)



def _flank_spacing_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.flank_spacing import main
    return main(argv)

def _fragment_lengths_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.fragment_lengths import main
    return main(argv)


def _fragment_heatmap_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.fragment_heatmap import main
    return main(argv)





def _fragments_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.fragments_command import main
    return main(argv)




def _mean_scale_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.mean_scale import main
    return main(argv)

def _merge_bams_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.merge_bams import main
    return main(argv)


def _randomize_fragments_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.randomize_fragments_command import main
    return main(argv)

def _gene_sets_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.gene_sets import main
    return main(argv)


def _gene_expression_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.gene_expression import main
    return main(argv)


def _tss_expression_quintiles_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.tss_expression_quintiles import main
    return main(argv)


def _resources_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.resource_files import main
    return main(argv)


def _validate_inputs_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.validate_inputs import main
    return main(argv)



def _cfdna_suite_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.cli.cfdna_suite import main
    return main(argv)

def _mnase_suite_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.cli.mnase_suite import main
    return main(argv)


def _nrl_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.nrl import main
    return main(argv)




def _plot_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.replot import main
    return main(argv)

def _positive_runs_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.positive_runs import main
    return main(argv)


def _peak_score_frequency_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.peak_score_frequency import main
    return main(argv)


def _peak_states_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.peak_states import main
    return main(argv)


def _region_extract_main(argv: Sequence[str] | None = None) -> int:
    from nucleosuite.pns_region_extractor import main
    return main(argv)


DELEGATED_COMMANDS: dict[str, tuple[CommandMain, str]] = {
    "fragments": (_fragments_main, "Write paired-end fragments as BED intervals."),
    "merge-bams": (_merge_bams_main, "Combine BAM files while keeping full alignment records."),
    "mean-scale": (_mean_scale_main, "Scale a BigWig relative to a supplied or calculated reference mean."),
    "randomize-fragments": (_randomize_fragments_main, "Create a reproducible control fragment set for comparison."),
    "fragment-lengths": (_fragment_lengths_main, "Count how many fragments occur at each length."),
    "flank-spacing": (_flank_spacing_main, "Compare nucleosome spacing around categorized reference sites."),
    "filter-coverage": (_filter_coverage_main, "Filter BED peaks by coverage at their summit or interval midpoint."),
    "fragment-heatmap": (_fragment_heatmap_main, "Compare fragment-length patterns across samples or region groups."),
    "aggregate": (_aggregate_main, "Aggregate BigWig signal around genomic features."),
    "region-extract": (_region_extract_main, "Export signal values and nearby peaks for each region in a BED file."),
    "compare-positions": (_compare_positions_main, "Compare one main callset with one or more called-position sets."),
    "combine": (_combine_main, "Combine per-contig outputs using raw counts and denominators."),
    "chrom-sizes": (_chrom_sizes_main, "Write chromosome names and lengths from a BAM or CRAM header."),
    "distances": (_distances_main, "Measure spacing between neighbouring or higher-order peaks."),
    "dac": (_dac_main, "Find repeating spacing patterns within one signal track."),
    "dcc": (_dcc_main, "Measure the positional offset and similarity between two signals."),
    "nrl": (_nrl_main, "Estimate nucleosome repeat length from DAC or DCC peaks."),
    "plot": (_plot_main, "Recreate and customise figures from NucleoSuite output tables."),
    "positive-runs": (_positive_runs_main, "Measure the lengths of continuous positive-signal regions."),
    "peak-score-frequency": (_peak_score_frequency_main, "Compare the score distributions of peak callsets."),
    "peak-states": (_peak_states_main, "Measure peak abundance and enrichment by chromatin state."),
    "gene-sets": (_gene_sets_main, "Group genes by the chromatin states that overlap them."),
    "gene-expression": (_gene_expression_main, "Compare gene expression with nucleosome spacing or periodicity."),
    "tss-expression-quintiles": (_tss_expression_quintiles_main, "Aggregate TSS signal after splitting genes into expression quintiles."),
    "mnase-suite": (_mnase_suite_main, "Run the configurable MNase analysis workflow."),
    "cfdna-suite": (_cfdna_suite_main, "Run the cfDNA fragmentomics and nucleosome workflow."),
    "resources": (_resources_main, "List, check or copy reference files bundled with NucleoSuite."),
    "validate-inputs": (_validate_inputs_main, "Validate suite inputs and reference compatibility."),
    "region-peak-extractor": (_region_extract_main, "Alias of region-extract."),
    "pns-region-extractor": (_region_extract_main, "Alias of region-extract."),
}


DELEGATED_MODULES: dict[str, str] = {
    "fragments": "nucleosuite.fragments_command",
    "merge-bams": "nucleosuite.merge_bams",
    "mean-scale": "nucleosuite.mean_scale",
    "randomize-fragments": "nucleosuite.randomize_fragments_command",
    "fragment-lengths": "nucleosuite.fragment_lengths",
    "flank-spacing": "nucleosuite.flank_spacing",
    "filter-coverage": "nucleosuite.filter_coverage",
    "fragment-heatmap": "nucleosuite.fragment_heatmap",
    "aggregate": "nucleosuite.cli.aggregate",
    "region-extract": "nucleosuite.pns_region_extractor",
    "compare-positions": "nucleosuite.compare_positions",
    "combine": "nucleosuite.combine",
    "chrom-sizes": "nucleosuite.chrom_sizes_command",
    "distances": "nucleosuite.distances",
    "dac": "nucleosuite.dac",
    "dcc": "nucleosuite.dcc",
    "nrl": "nucleosuite.nrl",
    "plot": "nucleosuite.replot",
    "positive-runs": "nucleosuite.positive_runs",
    "peak-score-frequency": "nucleosuite.peak_score_frequency",
    "peak-states": "nucleosuite.peak_states",
    "gene-sets": "nucleosuite.gene_sets",
    "gene-expression": "nucleosuite.gene_expression",
    "tss-expression-quintiles": "nucleosuite.tss_expression_quintiles",
    "resources": "nucleosuite.resource_files",
    "validate-inputs": "nucleosuite.validate_inputs",
    "mnase-suite": "nucleosuite.cli.mnase_suite",
    "cfdna-suite": "nucleosuite.cli.cfdna_suite",
}


def _is_informational_invocation(args: Sequence[str]) -> bool:
    if not args or args[0].startswith("-"):
        return True
    if args[0] == "resources":
        return True
    return any(token in {"-h", "--help", "--help-plotting", "--version"} for token in args[1:])


def _preparse_delegated(
    command: str, argv: Sequence[str]
) -> argparse.Namespace | None:
    """Run argparse validation before printing a job-start message."""
    module_name = DELEGATED_MODULES.get(command)
    if module_name is None:
        return None
    module = importlib.import_module(module_name)
    validator = getattr(module, "validate_argv", None)
    if validator is not None:
        validator(list(argv))
    custom_parser = getattr(module, "parse_cli_args", None)
    if custom_parser is not None:
        return custom_parser(list(argv))
    builder = getattr(module, "build_parser", None)
    if builder is not None:
        return builder().parse_args(list(argv))
    return None


def build_parser() -> argparse.ArgumentParser:
    """Build the central parser for native commands and top-level help."""
    parser = TopLevelHelpParser(
        prog="nucleosuite",
        description=(
            "Turn paired-end sequencing fragments into nucleosome signal tracks, "
            "peak calls, spacing measurements, sequence profiles and comparison plots."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Documentation: README.md and docs/QUICKSTART.md.\n"
            "Run 'nucleosuite COMMAND --help' for command-specific options and defaults."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"{render_logo()}\n\n%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        metavar="COMMAND",
        parser_class=argparse.ArgumentParser,
    )

    # BAM and track-generation commands.
    pns.register(subparsers)
    wps.register(subparsers)
    coverage.register(subparsers)
    dyads.register(subparsers)
    fragment_ends.register(subparsers)
    tracks.register(subparsers)
    dinuc_profile.register(subparsers)
    ww_types.register(subparsers)
    call_peaks.register(subparsers)

    # Analysis commands with dedicated parsers.
    # Placeholder parsers make them visible in top-level help; dispatch occurs
    # before this parser handles command-specific arguments.
    visible = (
        "fragments",
        "merge-bams",
        "mean-scale",
        "randomize-fragments",
        "fragment-lengths",
        "flank-spacing",
        "filter-coverage",
        "fragment-heatmap",
        "aggregate",
        "region-extract",
        "compare-positions",
        "combine",
        "chrom-sizes",
        "distances",
        "dac",
        "dcc",
        "nrl",
        "plot",
        "positive-runs",
        "peak-score-frequency",
        "peak-states",
        "gene-sets",
        "gene-expression",
        "tss-expression-quintiles",
        "mnase-suite",
        "cfdna-suite",
        "resources",
        "validate-inputs",
    )
    for name in visible:
        _, description = DELEGATED_COMMANDS[name]
        subparsers.add_parser(name, help=description, add_help=False)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args_list = list(sys.argv[1:] if argv is None else argv)
    informational = _is_informational_invocation(args_list)
    try:
        if args_list and args_list[0] in DELEGATED_COMMANDS:
            command = args_list[0]
            parsed = (
                _preparse_delegated(command, args_list[1:])
                if not informational
                else None
            )
            runner, _ = DELEGATED_COMMANDS[command]
            execute = lambda: int(runner(args_list[1:]) or 0)
        else:
            # Parse first so top-level help, version output, unknown commands and
            # argparse errors never receive a Warcraft startup message.
            parser = build_parser()
            args = parser.parse_args(args_list)
            command = str(args.command)
            parsed = args
            execute = lambda: int(args.command_runner(args) or 0)
    except (FileNotFoundError, OSError, ValueError, KeyError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    from nucleosuite.plotting import configure_plot_metadata
    configure_plot_metadata(command, args_list, vars(parsed) if parsed is not None else None)

    if parsed is not None and hasattr(parsed, "plot_format"):
        from nucleosuite.plotting import configure_plot_options
        configure_plot_options(parsed)

    if informational:
        try:
            return execute()
        except (FileNotFoundError, OSError, ValueError, KeyError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2

    from nucleosuite.command_logging import CommandLog

    parameters = vars(parsed) if parsed is not None else None
    with CommandLog(
        command,
        args_list,
        version=__version__,
        parameters=parameters,
    ) as command_log:
        selected_messages = message_pair()
        print_startup(selected_messages.startup)
        progress = ProgressReporter(command)
        progress.stage("Starting analysis")
        try:
            exit_code = execute()
        except (FileNotFoundError, OSError, ValueError, KeyError, RuntimeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            exit_code = 2
        if exit_code == 0:
            progress.complete("Analysis completed")
            print_completion(selected_messages.completion)
        command_log.finish(exit_code)
        return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
