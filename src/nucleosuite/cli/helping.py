"""Two-level command help for NucleoSuite.

Ordinary ``--help`` output is intentionally compact and shows only the inputs
and analysis controls most users need routinely. ``--help-all`` expands the
remaining command-specific controls. Shared plotting controls remain isolated
under ``--help-plotting`` where the command supports them.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping, Sequence


# Core options are selected deliberately per command. Positionals and required
# actions are always shown even when they are not listed here. Aliases are
# canonicalised below so one profile covers the synonymous command names.
CORE_DESTINATIONS: dict[str, set[str]] = {
    "pns": {
        "bamfiles", "fragment_files", "chrom_sizes", "blacklist_bed", "out_prefix",
        "contigs", "frag_lower", "frag_upper", "mode_length",
        "score_tracks", "interval_format", "cores",
    },
    "wps": {
        "bamfiles", "fragment_files", "chrom_sizes", "blacklist_bed", "out_prefix",
        "contigs", "frag_lower", "frag_upper", "protection", "baseline_window",
        "sg_window", "score_tracks", "peak_caller", "peak_track", "cores",
    },
    "coverage": {
        "bamfiles", "fragment_files", "chrom_sizes", "blacklist_bed", "out_prefix",
        "contigs", "frag_lower", "frag_upper", "output_format", "cores",
    },
    "dyads": {
        "bamfiles", "fragment_files", "chrom_sizes", "blacklist_bed", "out_prefix",
        "contigs", "frag_lower", "frag_upper", "even_dyad", "output_format", "cores",
    },
    "fragment-ends": {
        "bamfiles", "fragment_files", "chrom_sizes", "blacklist_bed", "out_prefix",
        "contigs", "frag_lower", "frag_upper", "tracks", "output_format", "cores",
    },
    "tracks": {
        "bamfiles", "fragment_files", "fasta", "chrom_sizes", "blacklist_bed",
        "contigs", "fragment_range", "spec_file", "output_dir", "output_prefix",
        "output_format", "cores",
    },
    "dinuc-profile": {
        "bamfiles", "fragment_files", "fasta", "chrom_sizes", "blacklist_bed",
        "out_prefix", "contigs", "frag_lower", "frag_upper", "dinuc_fraction", "cores",
    },
    "ww-types": {
        "bamfiles", "fragment_files", "fasta", "chrom_sizes", "blacklist_bed",
        "out_prefix", "contigs", "frag_lower", "frag_upper", "split_beds",
        "dyad_tracks", "output_format", "cores",
    },
    "call-peaks": {
        "input_bigwig", "blacklist_bed", "out_prefix", "method", "signal", "regions",
        "interval_format", "smooth_window", "min_region_length", "score_scale", "cores",
    },
    "fragments": {
        "bamfiles", "fragment_files", "output_prefix", "contigs", "frag_lower", "frag_upper",
        "max_duplicates", "dedup_scope", "output_format", "chrom_sizes", "cores",
    },
    "merge-bams": {
        "bamfiles", "output", "output_prefix", "contigs", "split_contigs", "frag_lower",
        "frag_upper", "max_duplicates", "dedup_scope", "cores",
    },
    "mean-scale": {
        "input", "reference_mean", "regions", "reference_bigwig", "score_column", "scale",
        "integer_scores", "clamp_min", "clamp_max", "output_format", "chrom_sizes", "output",
    },
    "randomize-fragments": {
        "bamfiles", "fragment_files", "fasta", "output_prefix", "contigs", "frag_lower",
        "frag_upper", "method", "seed", "search_window", "blacklist_bed", "fallback", "cores",
    },
    "fragment-lengths": {
        "bamfiles", "fragment_files", "blacklist_bed", "bed", "contigs", "min_length",
        "max_length", "output", "separate_files", "plot", "plot_min", "plot_max",
        "nrl_min_length", "nrl_max_length", "no_fragment_size_nrl", "cores",
    },
    "flank-spacing": {
        "nucleosome_bed", "region_bed", "category_col", "point_col", "nucleosome_center_col",
        "distribution", "ratio_x1", "ratio_x2", "top_categories", "x_min", "x_max",
        "output_dir", "output_prefix", "write_detail_tables",
    },
    "filter-peaks": {
        "input", "score_column", "min_score", "max_score", "score_percentile",
        "min_length", "max_length", "abs_score", "score_scale", "coverage_bigwig",
        "min_coverage", "coverage_position_column", "coverage_chunk_size", "output",
        "output_format", "chrom_sizes", "summary_output", "strict",
    },
    "fragment-heatmap": {
        "input", "out_prefix", "min_frag", "max_frag", "normalisation", "downsample_to",
        "cluster_method", "cluster_metric", "no_cluster", "metadata", "require_metadata",
        "write_detail_tables",
    },
    "aggregate": {
        "bigwig", "region_bed", "blacklist_bed", "nucleosome_bed", "state_bed", "output_dir",
        "output_prefix", "write_detail_tables", "window_half", "point_col", "strand_col",
        "category_col", "nrl", "nrl_peak_resolution", "contigs", "cores",
    },
    "region-extract": {
        "bed", "blacklist_bed", "coverage_bw", "pns_bw", "signal_track", "nucleosome_peaks",
        "breakpoint_peaks", "peak_track", "peak_flank_bp", "out_prefix", "cores",
    },
    "compare-positions": {
        "main_bed", "compare_beds", "score_bigwigs", "blacklist_bed", "max_distance",
        "percentile_interval", "distance_bins", "score_normalization", "score_correlation",
        "write_detail_tables", "output_prefix", "cores",
    },
    "combine": {
        "input_dir", "output_dir", "chrom_sizes", "skip_tracks", "force", "cores",
        "bigwig_method",
    },
    "chrom-sizes": {"alignment", "output", "contigs", "fasta"},
    "distances": {
        "input", "blacklist_bed", "state_bed", "position_column", "score_column",
        "score_percentile", "target_peaks", "min_distance", "max_distance", "max_order",
        "nrl_mode", "scope", "regression_scope", "label_peaks", "output_prefix", "cores",
    },
    "dac": {
        "bigwig", "blacklist_bed", "regions_bed", "genes_bed", "chrom_sizes", "scope",
        "chromosome", "window_size", "category", "max_distance", "output_prefix", "cores",
    },
    "dcc": {"mode", "cores"},
    "dcc:bigwig": {
        "bigwig_a", "bigwig_b", "regions_bed", "chrom_sizes", "scope", "chromosome",
        "dmax", "category", "output_prefix", "algorithm",
    },
    "dcc:bam": {
        "bam_a", "fragments_a", "bam_b", "fragments_b", "length_a", "min_length_a",
        "max_length_a", "length_b", "min_length_b", "max_length_b", "position_a", "position_b",
        "regions_bed", "chrom_sizes", "scope", "chromosome", "dmax", "output_prefix",
    },
    "nrl": {
        "input", "distance_column", "value_column", "min_distance", "max_distance",
        "peak_resolution", "output_prefix", "title",
    },
    "plot": {
        "input", "plot_type", "from_command", "output", "format", "width", "height", "dpi",
        "title", "no_title", "x_label", "y_label", "x_min", "x_max", "y_min", "y_max",
        # Common plot-family controls are shown only when that family exposes them.
        "label_peaks", "peak_label_value", "show_boxplot_outliers", "stats",
        "percentile_boxplot_y_max", "score_normalization", "score_correlation", "vmin", "vmax",
        "normalization", "bar_gap", "detect_peaks",
    },
    "positive-runs": {
        "bigwig", "blacklist_bed", "output_prefix", "contigs", "threshold", "min_run_length",
        "max_run_length", "plot_x_max", "normalization", "cores",
    },
    "peak-score-frequency": {
        "peaks", "output_prefix", "blacklist_bed", "score_column", "score_scale", "integer_bins",
        "bins", "bin_width", "score_min", "score_max", "normalization", "plot_x_min", "plot_x_max",
        "log_y",
    },
    "empirical-peak-fdr": {
        "sample_peaks", "randomized_peaks", "score_column", "fdr",
        "output_prefix", "output",
    },
    "cutn-suite": {
        "treatment1_bam", "control1_bam", "treatment2_bam", "control2_bam",
        "outdir", "sample_name", "condition1_name", "condition2_name",
        "inspect_run", "rerun_from", "exclude_sample", "bam_mode",
        "mode", "mode_strategy", "frag_mode_padding",
        "score_frag_lower", "score_frag_upper", "coverage_frag_lower",
        "coverage_frag_upper", "blacklist_bed", "contigs", "cores",
        "peak_min_region_length", "peak_max_neg_run", "stage1_coverage_statistic",
        "cluster_seed_mode", "cluster_seed_gate_mode", "stage1_gate_mode",
        "cluster_seed_p_value", "cluster_member_mode", "cluster_max_non_member_gap",
        "max_cluster_gap", "min_cluster_members", "differential_fdr",
    },
    "cutn-compare": {
        "condition1_results", "condition2_results", "outdir", "feature_level",
        "peak_match_distance", "fdr",
    },
    "peak-states": {
        "peaks", "state_bed", "blacklist_bed", "position_column", "score_column",
        "state_label_column", "state_color_column", "contigs", "overlap_policy", "score_percentile",
        "pct_values", "pct_bin_size", "output_prefix", "cores",
    },
    "gene-sets": {
        "genes_bed", "states_bed", "config", "gene_set", "chrom_sizes", "blacklist_bed",
        "output_dir", "output_prefix", "venn_sets", "cores",
    },
    "gene-expression": {
        "expression", "genes_bed", "resource_set", "peaks", "signal", "signal_type", "analysis",
        "output_prefix", "write_detail_tables", "focus_profile", "correlation", "cores",
    },
    "tss-expression-quintiles": {
        "signal", "sample", "signal_label", "expression", "tissue", "genes_bed", "resource_set",
        "window", "output_prefix",
    },
    "resources": {"action"},
    "validate-inputs": {
        "bam", "fragments", "bed", "blacklist_bed", "fasta", "chrom_sizes",
        "require_bam_index", "require_sorted_fragments", "max_records", "report",
    },
}


ALIASES: dict[str, str] = {
    "dyad": "dyads",
    "peak-call": "call-peaks",
    "peakcall": "call-peaks",
    "region-peak-extractor": "region-extract",
    "pns-region-extractor": "region-extract",
}


def canonical_help_key(command: str, nested: str | None = None) -> str:
    command = ALIASES.get(command, command)
    if nested:
        nested_key = f"{command}:{nested}"
        if nested_key in CORE_DESTINATIONS:
            return nested_key
    return command


def parser_has_option(parser: argparse.ArgumentParser, option: str) -> bool:
    return any(option in action.option_strings for action in parser._actions)


def visible_actions(parser: argparse.ArgumentParser) -> list[argparse.Action]:
    return [
        action
        for action in parser._actions
        if action.help is not argparse.SUPPRESS
        and action.help != argparse.SUPPRESS
        and not str(action.dest).startswith("_")
    ]


def _action_required_by_group(parser: argparse.ArgumentParser, action: argparse.Action) -> bool:
    for group in parser._mutually_exclusive_groups:
        if group.required and action in group._group_actions:
            return True
    return False


def core_actions(parser: argparse.ArgumentParser, key: str) -> list[argparse.Action]:
    selected = CORE_DESTINATIONS.get(key, set())
    actions: list[argparse.Action] = []
    fallback_budget = 12
    for action in visible_actions(parser):
        dest = str(action.dest)
        positional = not action.option_strings
        required = bool(getattr(action, "required", False)) or _action_required_by_group(parser, action)
        help_switch = dest in {"help", "help_plotting"}
        if positional or required or help_switch or dest in selected:
            actions.append(action)
    # Commands added in the future should still get a useful compact help page
    # before a profile is curated. Keep the first visible options rather than
    # falling back to the original enormous help output.
    if not selected:
        existing = {id(action) for action in actions}
        for action in visible_actions(parser):
            if id(action) in existing:
                continue
            actions.append(action)
            if len(actions) >= fallback_budget:
                break
    return actions


def _format_with_actions(
    parser: argparse.ArgumentParser,
    actions: Sequence[argparse.Action],
    *,
    command_display: str,
    extended: bool,
) -> str:
    formatter = parser._get_formatter()
    # Deliberately omit the full mutually-exclusive group structure from the
    # compact usage line; the option descriptions still explain valid choices.
    formatter.add_usage(parser.usage, list(actions), [])
    formatter.add_text(parser.description)

    action_ids = {id(action) for action in actions}
    for group in parser._action_groups:
        group_actions = [action for action in group._group_actions if id(action) in action_ids]
        if not group_actions:
            continue
        formatter.start_section(group.title)
        formatter.add_text(group.description)
        formatter.add_arguments(group_actions)
        formatter.end_section()

    if extended and parser.epilog:
        formatter.add_text(parser.epilog)

    plotting_note = ""
    if parser_has_option(parser, "--help-plotting"):
        plotting_note = (
            f"\nRun 'nucleosuite {command_display} --help-plotting' for shared plot customization options."
        )
    if extended:
        formatter.add_text(
            f"Full command-specific help shown.{plotting_note}"
        )
    else:
        formatter.add_text(
            "Core options shown. "
            f"Run 'nucleosuite {command_display} --help-all' for all command-specific options."
            f"{plotting_note}"
        )
    return formatter.format_help()


def format_layered_help(
    parser: argparse.ArgumentParser,
    command: str,
    *,
    nested: str | None = None,
    extended: bool = False,
) -> str:
    """Format core or extended help without mutating the parser."""

    key = canonical_help_key(command, nested)
    actions = visible_actions(parser) if extended else core_actions(parser, key)
    display = command if nested is None else f"{command} {nested}"
    return _format_with_actions(
        parser,
        actions,
        command_display=display,
        extended=extended,
    )


def resolve_nested_parser(
    parser: argparse.ArgumentParser,
    tokens: Sequence[str],
) -> tuple[argparse.ArgumentParser, str | None]:
    """Resolve a directly named argparse subcommand, such as ``dcc bigwig``."""

    current = parser
    nested: str | None = None
    for token in tokens:
        if token.startswith("-"):
            continue
        subparsers = [
            action for action in current._actions
            if isinstance(action, argparse._SubParsersAction)
        ]
        if not subparsers:
            break
        chosen = None
        for subparser_action in subparsers:
            if token in subparser_action.choices:
                chosen = subparser_action.choices[token]
                nested = token
                break
        if chosen is None:
            break
        current = chosen
    return current, nested


SUITE_CORE_HELP: Mapping[str, str] = {
    "cfdna-suite": """\
usage: nucleosuite cfdna-suite --bam BAM [BAM ...] --fasta FASTA --outdir DIR [core options]

Run the cfDNA fragmentomics and nucleosome workflow.

core options:
  --bam FILE [MORE ...]         Coordinate-sorted paired-end BAM input; mutually exclusive with --fragments.
  --fragments FILE [MORE ...]   Fragment BED/BED.gz/bigBed input; mutually exclusive with --bam.
  --fasta FILE                  Matching reference FASTA.
  --resource-set NAME           Bundled resource collection, for example hg19-gm12878.
  --sample-name NAME            Output sample name; otherwise derived from inputs.
  --outdir DIR                  Output directory.
  --cores N                     Maximum concurrent contig workers (default: 1).
  --analysis-scope VALUE        combined-only (default) or per-contig-and-combined.
  --score-frag-lower N            PNS scoring lower fragment length (default: 137).
  --score-frag-upper N            PNS scoring upper fragment length (default: 197).
  --score-mode-length N           PNS scoring modal fragment length (default: 167).
  --with-randomized-control     Run complete observed and randomized workflows, then annotate combined peak BEDs with empirical FDR.
  --fdr N                       In paired mode, also write combined peak BEDs filtered at FDR N.
  --max-duplicates N            Identical-fragment copy limit (default: 1).
  --blacklist-bed FILE          Override the assembly-specific blacklist.
  --no-blacklist                Disable blacklist filtering.
  --help-plotting               Show shared plot customization options and exit.

Core options shown. Run 'nucleosuite cfdna-suite --help-all' for all command-specific options.
""",
    "mnase-suite": """\
usage: nucleosuite mnase-suite --bam BAM [BAM ...] --fasta FASTA --outdir DIR [core options]

Run the configurable MNase analysis workflow.

core options:
  --bam FILE [MORE ...]         Coordinate-sorted paired-end BAM input; mutually exclusive with --fragments.
  --fragments FILE [MORE ...]   Fragment BED/BED.gz/bigBed input; mutually exclusive with --bam.
  --fasta FILE                  Matching reference FASTA.
  --resource-set NAME           Bundled resource collection, for example hg19-gm12878.
  --ctcf-bed FILE               CTCF coordinates for aggregate analyses.
  --sample-name NAME            Output sample name; otherwise derived from inputs.
  --outdir DIR                  Output directory.
  --cores N                     Maximum concurrent contig workers (default: 1).
  --analysis-scope VALUE        combined-only (default) or per-contig-and-combined.
  --score-frag-lower N            PNS scoring lower fragment length (default: 120).
  --score-frag-upper N            PNS scoring upper fragment length (default: 180).
  --score-mode-length N           PNS scoring modal fragment length (default: 147).
  --fine-frag-lower/upper N     Ranged dyad/WW class (default: 146-148).
  --exact-size N                Exact dyad/fragment-end length (default: 147).
  --with-randomized-control     Run complete observed and randomized workflows, then annotate combined peak BEDs with empirical FDR.
  --fdr N                       In paired mode, also write combined peak BEDs filtered at FDR N.
  --max-duplicates N            Identical-fragment copy limit.
  --blacklist-bed FILE          Override the assembly-specific blacklist.
  --no-blacklist                Disable blacklist filtering.
  --help-plotting               Show shared plot customization options and exit.

Core options shown. Run 'nucleosuite mnase-suite --help-all' for all command-specific options.
""",
}
