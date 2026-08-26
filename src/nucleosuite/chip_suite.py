#!/usr/bin/env python3
"""Stage 1 and optional Stage 2 ChIP/CUT&RUN/CUT&Tag workflow."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from nucleosuite import __version__
from nucleosuite.bigwig_ops import average_bigwigs
from nucleosuite.chip_compare import (
    MANIFEST_NAME,
    MANIFEST_SCHEMA,
    MANIFEST_SCHEMA_VERSION,
    compare_stage1,
)
from nucleosuite.chip_aggregate import run_cluster_aggregate, write_cluster_anchor_bed
from nucleosuite.chip_peaks import analyze_chip_replicate_peaks
from nucleosuite.mean_scale import scale_bigwig_by_reference
from nucleosuite.mode_estimation import (
    ModeEstimate,
    estimate_bam_fragment_mode,
    mode_estimate_message,
    pooled_mode_estimate,
)
from nucleosuite.progress import ProgressReporter


_TRACK_NAMES = {
    "tns": ("tns", "posTNS"),
    "bns": ("bns", "posBNS"),
    "pns": ("pns", "posPNS"),
}
_MULTICONTIG_MANIFEST = "nucleosuite_multicontig_manifest.json"


def _mode_value(text: str) -> str | int:
    value = str(text).strip().lower()
    if value == "auto":
        return "auto"
    try:
        integer = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--mode must be auto or an integer") from exc
    if integer < 3:
        raise argparse.ArgumentTypeError("--mode must be at least 3")
    return integer


def _run_nucleosuite(arguments: Sequence[str]) -> None:
    command = [sys.executable, "-m", "nucleosuite.cli.main", *map(str, arguments)]
    completed = subprocess.run(command, check=False)
    if completed.returncode:
        raise RuntimeError(
            f"NucleoSuite subcommand failed with exit code {completed.returncode}: "
            + " ".join(map(str, arguments))
        )


def _resolved_prefix(
    base: Path,
    method: str,
    mode: int,
    frag_lower: int,
    frag_upper: int,
) -> Path:
    return Path(
        f"{base}_method{method}_mode{mode}_lower{frag_lower}_upper{frag_upper}_smooth0x2"
    )


def _scoring_fragment_range(args: argparse.Namespace, mode: int) -> tuple[int, int]:
    if args.score_frag_lower is not None and args.score_frag_upper is not None:
        return int(args.score_frag_lower), int(args.score_frag_upper)
    return max(1, int(mode) - int(args.score_fragment_flank)), int(mode) + int(
        args.score_fragment_flank
    )


def _resolved_coverage_prefix(base: Path, lower: int, upper: int) -> Path:
    return Path(f"{base}_coverage_lower{lower}_upper{upper}")


def _locate_completed_prefix(
    requested_base: Path,
    direct_prefix: Path,
    required_suffixes: Sequence[str],
) -> Path:
    """Locate serial or manifest-declared combined multicontig outputs."""

    candidates = [direct_prefix]
    multicontig_root = requested_base.parent / f"{requested_base.name}_multicontig"
    manifest_path = multicontig_root / _MULTICONTIG_MANIFEST
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(
                f"Could not read multicontig output manifest: {manifest_path}"
            ) from exc
        combined_dir = Path(
            str(manifest.get("combined_dir") or (multicontig_root / "combined"))
        )
        if not combined_dir.is_absolute():
            combined_dir = multicontig_root / combined_dir
        combined_name = str(manifest.get("combined_name") or direct_prefix.name)
        candidates.append(combined_dir / combined_name)
    else:
        candidates.append(multicontig_root / "combined" / direct_prefix.name)

    for prefix in candidates:
        if all(Path(f"{prefix}{suffix}").is_file() for suffix in required_suffixes):
            return prefix

    expected = ", ".join(
        str(Path(f"{prefix}{suffix}"))
        for prefix in candidates
        for suffix in required_suffixes
    )
    raise RuntimeError(
        f"Expected completed outputs were not created; searched: {expected}"
    )


def _write_mode_report(
    path: Path,
    *,
    mode_source: str,
    target: ModeEstimate | None,
    control: ModeEstimate | None,
    pooled: ModeEstimate | None,
    target_mode: int,
    control_mode: int,
    seeds: tuple[int, int, int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "sample\tmode\tbootstrap_ci_low\tbootstrap_ci_high\t"
            "sampled_fragments\tmode_search_fragments\tmode_search_lower\t"
            "mode_search_upper\thistogram_smoothing\tconverged\tcheckpoints\t"
            "seed\tmode_source\thistogram_length_bp_count\n"
        )
        for label, estimate, analysis_mode, seed in (
            ("treatment", target, target_mode, seeds[0]),
            ("control", control, control_mode, seeds[1]),
            ("pooled", pooled, target_mode if target_mode == control_mode else -1, seeds[2]),
        ):
            if estimate is None:
                if label != "pooled":
                    handle.write(
                        f"{label}\t{analysis_mode}\t\t\t\t\t\t\t"
                        f"not_applicable\t\t\t{seed}\t{mode_source}\t\n"
                    )
                continue
            histogram = ";".join(
                f"{estimate.search_lower + index}:{count}"
                for index, count in enumerate(estimate.histogram)
            )
            handle.write(
                f"{label}\t{estimate.mode}\t{estimate.ci_low:.6g}\t"
                f"{estimate.ci_high:.6g}\t{estimate.sampled_fragments}\t"
                f"{estimate.mode_search_fragments}\t{estimate.search_lower}\t"
                f"{estimate.search_upper}\t{estimate.histogram_smoothing}\t"
                f"{str(estimate.converged).lower()}\t"
                f"{estimate.checkpoints}\t{seed}\t{mode_source}\t{histogram}\n"
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite chip-suite",
        description=(
            "Run control-aware Stage 1 ChIP/CUT&RUN/CUT&Tag scoring for one "
            "condition and optionally compare it with a second condition."
        ),
    )
    parser.add_argument(
        "--treatment1-bam", "--target-bam",
        dest="treatment1_bam", nargs="+", required=True,
        help="Condition 1 treatment BAM file(s); --target-bam is a compatibility alias.",
    )
    parser.add_argument(
        "--control1-bam", "--control-bam",
        dest="control1_bam", nargs="+", required=True,
        help="Condition 1 control BAM file(s); --control-bam is a compatibility alias.",
    )
    parser.add_argument(
        "--treatment2-bam", nargs="+",
        help="Optional condition 2 treatment BAM file(s).",
    )
    parser.add_argument(
        "--control2-bam", nargs="+",
        help="Optional condition 2 control BAM file(s).",
    )
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--sample-name", default="chip", help="Output prefix (default: chip).")
    parser.add_argument("--condition1-name", help="Condition 1 label; defaults to --sample-name.")
    parser.add_argument("--condition2-name", default="condition2", help="Condition 2 label (default: condition2).")
    parser.add_argument(
        "--bam-mode", choices=("replicates", "merged"), default="replicates",
        help=(
            "Treat BAMs as independent biological replicates, or pool each BAM "
            "group as one logical merged sample (default: replicates)."
        ),
    )
    parser.add_argument(
        "--scoring-method", choices=("pns", "bns", "tns"), default="pns",
        help="Nucleosome scoring kernel (default: pns).",
    )
    parser.add_argument(
        "--mode", type=_mode_value, default="auto",
        help=(
            "Scoring mode: auto uses bootstrap-stabilized random-block sampling; "
            "an integer bypasses estimation (default: auto)."
        ),
    )
    parser.add_argument(
        "--mode-strategy", choices=("pooled", "separate", "target", "control"),
        default="pooled",
        help="How group estimates define treatment and control analysis modes (default: pooled).",
    )
    parser.add_argument(
        "--score-fragment-flank", type=int, default=30,
        help="Default scoring range extends this many bp below and above the resolved mode (default: 30).",
    )
    parser.add_argument(
        "--score-frag-lower", "--frag-lower", dest="score_frag_lower", type=int,
        help="Explicit scoring-fragment lower bound; supply with --score-frag-upper to override mode +/- flank.",
    )
    parser.add_argument(
        "--score-frag-upper", "--frag-upper", dest="score_frag_upper", type=int,
        help="Explicit scoring-fragment upper bound; supply with --score-frag-lower to override mode +/- flank.",
    )
    parser.add_argument(
        "--coverage-frag-lower", type=int, default=1,
        help="Minimum fragment length contributing to coverage measurements (default: 1).",
    )
    parser.add_argument(
        "--coverage-frag-upper", type=int, default=1000,
        help="Maximum fragment length contributing to coverage measurements (default: 1000).",
    )
    parser.add_argument("--mode-search-lower", type=int, default=120, help="Automatic mode-search lower bound (default: 120).")
    parser.add_argument("--mode-search-upper", type=int, default=250, help="Automatic mode-search upper bound (default: 250).")
    parser.add_argument("--mode-min-fragments", type=int, default=100_000, help="Minimum sampled fragments before convergence checks (default: 100000).")
    parser.add_argument("--mode-batch-fragments", type=int, default=25_000, help="Fragments between convergence checks (default: 25000).")
    parser.add_argument("--mode-max-fragments", type=int, default=1_000_000, help="Maximum sampled fragments per BAM group (default: 1000000).")
    parser.add_argument("--mode-bootstrap", type=int, default=200, help="Bootstrap mode replicates per checkpoint (default: 200).")
    parser.add_argument("--mode-stable-checkpoints", type=int, default=3, help="Required stable checkpoints (default: 3).")
    parser.add_argument("--mode-max-change", type=int, default=1, help="Maximum mode range across stable checkpoints (default: 1 bp).")
    parser.add_argument("--mode-max-ci-width", type=float, default=4.0, help="Maximum bootstrap 95%% interval width (default: 4 bp).")
    parser.add_argument("--mode-block-bp", type=int, default=1_000_000, help="Random genomic-block size (default: 1000000 bp).")
    parser.add_argument(
        "--mode-histogram-smoothing",
        choices=("none", "binomial"),
        default="none",
        help=(
            "Histogram processing for automatic mode estimation: none uses raw "
            "integer counts; binomial explicitly applies the optional 1,4,6,4,1 "
            "kernel (default: none)."
        ),
    )
    parser.add_argument("--seed", type=int, default=12345, help="Random seed (default: 12345).")
    parser.add_argument("--blacklist-bed", help="BED blacklist applied to every BAM group.")
    parser.add_argument("-c", "--contigs", nargs="+", default=["autosomes"], help="Contigs analysed in every sample (default: autosomes).")
    parser.add_argument("--max-duplicates", type=int, default=1, help="Maximum identical fragments retained (default: 1).")
    parser.add_argument(
        "--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams",
        help="Coordinate duplicate scope within a logical sample (default: all_bams).",
    )
    parser.add_argument("--cores", type=int, default=1, help="Parallel contig workers (default: 1).")
    parser.add_argument(
        "--stage1-p-value", type=float,
        help=(
            "Optional exploratory one-sided Stage 1 p-value cutoff. By default, "
            "the all-controls gate alone selects peaks."
        ),
    )
    parser.add_argument(
        "--peak-fdr", type=float,
        help=(
            "Optional Stage 1 peak FDR cutoff. No Stage 1 FDR cutoff is applied "
            "by default."
        ),
    )
    parser.add_argument(
        "--cluster-fdr", type=float,
        help=(
            "Optional Stage 1 maximum-seed FDR cutoff. No cluster FDR "
            "cutoff is applied by default."
        ),
    )
    parser.add_argument(
        "--cluster-seed-p-value", type=float, default=0.05,
        help=(
            "One-sided p-value required, together with the all-controls gate, "
            "for a treatment peak to seed a cluster (default: 0.05)."
        ),
    )
    parser.add_argument(
        "--cluster-member-p-value", dest="cluster_seed_p_value", type=float,
        default=argparse.SUPPRESS, help=argparse.SUPPRESS,
    )
    parser.add_argument("--differential-fdr", type=float, default=0.05, help="Stage 2 differential FDR cutoff (default: 0.05).")
    parser.add_argument(
        "--cluster-max-non-gated-gap", type=int, default=1,
        help="Consecutive non-gated peaks that may bridge gated cluster members (default: 1).",
    )
    parser.add_argument(
        "--cluster-break", dest="cluster_max_non_gated_gap", type=int,
        default=argparse.SUPPRESS, help=argparse.SUPPRESS,
    )
    parser.add_argument("--max-cluster-gap", type=int, default=1000, help="Maximum distance between adjacent gated-member summits (default: 1000 bp).")
    parser.add_argument(
        "--min-cluster-gated-peaks", type=int, default=2,
        help="Minimum total gate-passing members per seeded cluster (default: 2).",
    )
    parser.add_argument(
        "--min-significant-peaks", dest="min_cluster_gated_peaks", type=int,
        default=argparse.SUPPRESS, help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--cluster-aggregate-window-half", type=int, default=1000,
        help="Bases plotted on either side of the strongest peak in each cluster (default: 1000).",
    )
    parser.add_argument(
        "--cluster-aggregate-max-heatmap-rows", type=int, default=5000,
        help="Maximum randomly sampled cluster rows in each heatmap (default: 5000).",
    )
    parser.add_argument(
        "--cluster-aggregate-bootstrap", type=int, default=200,
        help="Cluster bootstrap replicates for aggregate 95%% confidence bands (default: 200; 0 disables).",
    )
    parser.add_argument(
        "--cluster-aggregate-nrl-resolution", type=float, default=140.0,
        help="Peak resolution for cluster-aligned directional NRLs (default: 140 bp).",
    )
    parser.add_argument(
        "--cluster-aggregate-nrl-min-order", type=int, default=0,
        help="First peak order used for cluster-aligned directional NRLs (default: 0).",
    )
    parser.add_argument(
        "--cluster-aggregate-nrl-max-order", type=int, default=3,
        help="Last peak order used for cluster-aligned directional NRLs (default: 3).",
    )
    parser.add_argument(
        "--skip-cluster-aggregate", action="store_true",
        help="Skip cluster-aligned PNS profiles, heatmaps and directional NRLs.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Validate and report planned stages without processing BAMs.")
    return parser


def _validate(args: argparse.Namespace) -> None:
    if (
        args.coverage_frag_lower < 1
        or args.coverage_frag_upper < args.coverage_frag_lower
    ):
        raise ValueError(
            "Require 1 <= --coverage-frag-lower <= --coverage-frag-upper"
        )
    if args.score_fragment_flank < 0:
        raise ValueError("--score-fragment-flank must be non-negative")
    if (args.score_frag_lower is None) != (args.score_frag_upper is None):
        raise ValueError(
            "--score-frag-lower and --score-frag-upper must be supplied together"
        )
    if args.score_frag_lower is not None and (
        args.score_frag_lower < 1 or args.score_frag_upper < args.score_frag_lower
    ):
        raise ValueError("Require 1 <= --score-frag-lower <= --score-frag-upper")
    if args.mode == "auto":
        if (
            args.mode_search_lower < args.coverage_frag_lower
            or args.mode_search_upper > args.coverage_frag_upper
        ):
            raise ValueError(
                "Automatic mode-search bounds must lie within the coverage fragment range"
            )
        if args.mode_search_upper < args.mode_search_lower:
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
    if args.cores < 1 or args.max_duplicates < 0:
        raise ValueError("--cores must be positive and --max-duplicates non-negative")
    if args.cluster_max_non_gated_gap < 0:
        raise ValueError("--cluster-max-non-gated-gap must be non-negative")
    if min(args.max_cluster_gap, args.min_cluster_gated_peaks) < 1:
        raise ValueError("Cluster gap and minimum gated-peak count must be positive")
    if args.cluster_aggregate_window_half < 1:
        raise ValueError("--cluster-aggregate-window-half must be positive")
    if args.cluster_aggregate_max_heatmap_rows < 1:
        raise ValueError("--cluster-aggregate-max-heatmap-rows must be positive")
    if args.cluster_aggregate_bootstrap < 0:
        raise ValueError("--cluster-aggregate-bootstrap must be non-negative")
    if args.cluster_aggregate_nrl_resolution < 0:
        raise ValueError("--cluster-aggregate-nrl-resolution must be non-negative")
    if (
        args.cluster_aggregate_nrl_min_order < 0
        or args.cluster_aggregate_nrl_max_order < args.cluster_aggregate_nrl_min_order
    ):
        raise ValueError("Cluster aggregate NRL orders must satisfy 0 <= min <= max")
    for name in (
        "stage1_p_value",
        "peak_fdr",
        "cluster_seed_p_value",
        "cluster_fdr",
        "differential_fdr",
    ):
        raw_value = getattr(args, name)
        if raw_value is None:
            continue
        value = float(raw_value)
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be between 0 and 1")
    if (args.treatment2_bam is None) != (args.control2_bam is None):
        raise ValueError("--treatment2-bam and --control2-bam must be supplied together")
    conditions = [(args.treatment1_bam, args.control1_bam)]
    if args.treatment2_bam is not None:
        conditions.append((args.treatment2_bam, args.control2_bam))
    for name in (args.sample_name, args.condition1_name, args.condition2_name):
        if name is not None and (not str(name).strip() or Path(str(name)).name != str(name)):
            raise ValueError("Sample and condition names must be non-empty path-safe names")
    paths = [*args.treatment1_bam, *args.control1_bam]
    if args.treatment2_bam is not None:
        paths.extend([*args.treatment2_bam, *args.control2_bam])
    for path in paths:
        if not Path(path).is_file():
            raise FileNotFoundError(path)


def _estimate_modes(
    args: argparse.Namespace,
    conditions: Sequence[tuple[Sequence[str], Sequence[str]]],
    reporter: ProgressReporter,
) -> tuple[int, int, str, list[tuple[ModeEstimate | None, ModeEstimate | None, ModeEstimate | None]]]:
    if args.mode != "auto":
        mode = int(args.mode)
        return mode, mode, "explicit", [(None, None, None) for _ in conditions]

    common = dict(
        frag_lower=args.coverage_frag_lower,
        frag_upper=args.coverage_frag_upper,
        search_lower=args.mode_search_lower,
        search_upper=args.mode_search_upper,
        minimum_fragments=args.mode_min_fragments,
        batch_fragments=args.mode_batch_fragments,
        maximum_fragments=args.mode_max_fragments,
        stable_checkpoints=args.mode_stable_checkpoints,
        maximum_mode_change=args.mode_max_change,
        maximum_ci_width=args.mode_max_ci_width,
        bootstrap_replicates=args.mode_bootstrap,
        block_bp=args.mode_block_bp,
        blacklist_bed=args.blacklist_bed,
        max_duplicates=args.max_duplicates,
        dedup_scope=args.dedup_scope,
        contig_tokens=args.contigs,
        histogram_smoothing=args.mode_histogram_smoothing,
    )
    estimates: list[tuple[ModeEstimate, ModeEstimate, ModeEstimate]] = []
    for index, (treatment_bams, control_bams) in enumerate(conditions):
        reporter.stage(f"Estimating condition {index + 1} treatment fragment mode")
        target = estimate_bam_fragment_mode(
            treatment_bams, **{**common, "seed": args.seed + index * 10}
        )
        reporter.stage(
            mode_estimate_message(
                f"Condition {index + 1} treatment fragment mode", target
            )
        )
        reporter.stage(f"Estimating condition {index + 1} control fragment mode")
        control = estimate_bam_fragment_mode(
            control_bams, **{**common, "seed": args.seed + index * 10 + 1}
        )
        reporter.stage(
            mode_estimate_message(
                f"Condition {index + 1} control fragment mode", control
            )
        )
        pooled = pooled_mode_estimate(
            target,
            control,
            bootstrap_replicates=args.mode_bootstrap,
            seed=args.seed + index * 10 + 2,
            histogram_smoothing=args.mode_histogram_smoothing,
        )
        reporter.stage(
            mode_estimate_message(
                f"Condition {index + 1} pooled fragment mode", pooled
            )
        )
        estimates.append((target, control, pooled))
        if abs(target.mode - control.mode) > 5:
            print(
                f"WARNING: condition {index + 1} treatment and control fragment modes "
                "differ by more than 5 bp; inspect the mode report",
                file=sys.stderr,
            )

    if len(estimates) == 1:
        target_pooled, control_pooled, all_pooled = estimates[0]
    else:
        target_pooled = pooled_mode_estimate(
            estimates[0][0], estimates[1][0],
            bootstrap_replicates=args.mode_bootstrap,
            seed=args.seed + 100,
            histogram_smoothing=args.mode_histogram_smoothing,
        )
        control_pooled = pooled_mode_estimate(
            estimates[0][1], estimates[1][1],
            bootstrap_replicates=args.mode_bootstrap,
            seed=args.seed + 101,
            histogram_smoothing=args.mode_histogram_smoothing,
        )
        all_pooled = pooled_mode_estimate(
            estimates[0][2], estimates[1][2],
            bootstrap_replicates=args.mode_bootstrap,
            seed=args.seed + 102,
            histogram_smoothing=args.mode_histogram_smoothing,
        )
    if args.mode_strategy == "pooled":
        target_mode = control_mode = all_pooled.mode
    elif args.mode_strategy == "separate":
        target_mode, control_mode = target_pooled.mode, control_pooled.mode
    elif args.mode_strategy == "target":
        target_mode = control_mode = target_pooled.mode
    else:
        target_mode = control_mode = control_pooled.mode
    reporter.stage(
        f"Resolved analysis modes: treatment={target_mode} bp; "
        f"control={control_mode} bp; strategy={args.mode_strategy}"
    )
    return target_mode, control_mode, f"automatic_{args.mode_strategy}", list(estimates)


def _generate_and_scale(
    *,
    args: argparse.Namespace,
    reporter: ProgressReporter,
    bam_paths: Sequence[str],
    label: str,
    mode: int,
    tracks_dir: Path,
    scaled_dir: Path,
    setup_dir: Path,
    score_track: str,
    positive_track: str,
    seed: int,
) -> dict[str, object]:
    score_lower, score_upper = _scoring_fragment_range(args, mode)
    reporter.stage(
        f"Generating {label} score tracks from {score_lower}-{score_upper} bp fragments"
    )
    base = tracks_dir / f"{args.sample_name}_{label}_{args.scoring_method}_discovery"
    command = [
        "pns", "--bam", *bam_paths,
        "--scoring-method", args.scoring_method,
        "--mode-length", str(mode),
        "--frag-lower", str(score_lower),
        "--frag-upper", str(score_upper),
        "--score-tracks", score_track, positive_track,
        "--other-tracks", "none",
        "--pns-format", "bigwig",
        "--other-format", "none",
        "--interval-format", "bed",
        "--no-peak-calling",
        "--out-prefix", str(base),
        "--contigs", *args.contigs,
        "--max-duplicates", str(args.max_duplicates),
        "--dedup-scope", args.dedup_scope,
        "--cores", str(args.cores),
        "--seed", str(seed),
    ]
    if args.blacklist_bed:
        command.extend(["--blacklist-bed", args.blacklist_bed])
    _run_nucleosuite(command)

    direct_prefix = _resolved_prefix(
        base, args.scoring_method, mode, score_lower, score_upper
    )
    prefix = _locate_completed_prefix(
        base,
        direct_prefix,
        (f"_{score_track}.bw", f"_{positive_track}.bw"),
    )
    score_path = Path(f"{prefix}_{score_track}.bw").resolve()
    positive_path = Path(f"{prefix}_{positive_track}.bw").resolve()
    coverage_base = tracks_dir / f"{args.sample_name}_{label}_coverage_all_fragments"
    reporter.stage(
        f"Generating {label} raw coverage from {args.coverage_frag_lower}-"
        f"{args.coverage_frag_upper} bp fragments"
    )
    coverage_command = [
        "coverage", "--bam", *bam_paths,
        "--frag-lower", str(args.coverage_frag_lower),
        "--frag-upper", str(args.coverage_frag_upper),
        "--output-format", "bigwig",
        "--out-prefix", str(coverage_base),
        "--contigs", *args.contigs,
        "--max-duplicates", str(args.max_duplicates),
        "--dedup-scope", args.dedup_scope,
        "--cores", str(args.cores),
        "--seed", str(seed + 1),
    ]
    if args.blacklist_bed:
        coverage_command.extend(["--blacklist-bed", args.blacklist_bed])
    _run_nucleosuite(coverage_command)
    direct_coverage_prefix = _resolved_coverage_prefix(
        coverage_base, args.coverage_frag_lower, args.coverage_frag_upper
    )
    coverage_prefix = _locate_completed_prefix(
        coverage_base, direct_coverage_prefix, ("_coverage.bw",)
    )
    coverage_path = Path(f"{coverage_prefix}_coverage.bw").resolve()
    scaled_path = (
        scaled_dir / f"{args.sample_name}_{label}_{score_track}_divided_by_mean_{positive_track}.bw"
    ).resolve()
    reporter.stage(f"Scaling {label} {score_track} by mean {positive_track}")
    _scaled, reference_mean, reference_count = scale_bigwig_by_reference(
        score_path, positive_path, scaled_path, scale=1.0, reporter=reporter
    )
    scaled_coverage_path = (
        scaled_dir / f"{args.sample_name}_{label}_coverage_divided_by_mean_x100.bw"
    ).resolve()
    reporter.stage(f"Scaling {label} coverage to a non-zero mean of 100")
    _coverage_scaled, coverage_mean, coverage_count = scale_bigwig_by_reference(
        coverage_path,
        coverage_path,
        scaled_coverage_path,
        scale=100.0,
        reporter=reporter,
    )
    with (setup_dir / f"{args.sample_name}_{label}_score_scaling.tsv").open(
        "wt", encoding="utf-8"
    ) as handle:
        handle.write("field\tvalue\n")
        handle.write(f"score_track\t{score_path}\n")
        handle.write(f"score_fragment_range\t{score_lower}-{score_upper}\n")
        handle.write(f"positive_reference_track\t{positive_path}\n")
        handle.write(f"raw_coverage_track\t{coverage_path}\n")
        handle.write(
            f"coverage_fragment_range\t{args.coverage_frag_lower}-"
            f"{args.coverage_frag_upper}\n"
        )
        handle.write(f"positive_reference_mean\t{reference_mean:.12g}\n")
        handle.write(f"finite_nonzero_reference_bases\t{reference_count}\n")
        handle.write(f"scaled_track\t{scaled_path}\n")
        handle.write(f"coverage_nonzero_mean\t{coverage_mean:.12g}\n")
        handle.write(f"coverage_finite_nonzero_bases\t{coverage_count}\n")
        handle.write(f"coverage_scaled_to_mean_100\t{scaled_coverage_path}\n")
    return {
        "bams": [str(Path(path).resolve()) for path in bam_paths],
        "score": str(score_path),
        "positive_score": str(positive_path),
        "coverage": str(coverage_path),
        "scaled_score": str(scaled_path),
        "scaled_coverage": str(scaled_coverage_path),
        "positive_score_mean": float(reference_mean),
        "positive_score_finite_nonzero_bases": int(reference_count),
        "coverage_nonzero_mean": float(coverage_mean),
        "coverage_finite_nonzero_bases": int(coverage_count),
        "score_frag_lower": score_lower,
        "score_frag_upper": score_upper,
        "coverage_frag_lower": int(args.coverage_frag_lower),
        "coverage_frag_upper": int(args.coverage_frag_upper),
    }


def _generate_pns_for_aggregate(
    *,
    args: argparse.Namespace,
    reporter: ProgressReporter,
    bam_paths: Sequence[str],
    label: str,
    mode: int,
    tracks_dir: Path,
    scaled_dir: Path,
    setup_dir: Path,
    seed: int,
) -> dict[str, object]:
    """Generate the PNS/posPNS pair used only for cluster-centred aggregates."""

    reporter.stage(f"Generating {label} PNS and posPNS aggregate tracks")
    score_lower, score_upper = _scoring_fragment_range(args, mode)
    base = tracks_dir / f"{args.sample_name}_{label}_pns_aggregate"
    command = [
        "pns", "--bam", *bam_paths,
        "--scoring-method", "pns",
        "--mode-length", str(mode),
        "--frag-lower", str(score_lower),
        "--frag-upper", str(score_upper),
        "--score-tracks", "pns", "posPNS",
        "--other-tracks", "none",
        "--pns-format", "bigwig",
        "--other-format", "none",
        "--interval-format", "bed",
        "--no-peak-calling",
        "--out-prefix", str(base),
        "--contigs", *args.contigs,
        "--max-duplicates", str(args.max_duplicates),
        "--dedup-scope", args.dedup_scope,
        "--cores", str(args.cores),
        "--seed", str(seed),
    ]
    if args.blacklist_bed:
        command.extend(["--blacklist-bed", args.blacklist_bed])
    _run_nucleosuite(command)
    direct_prefix = _resolved_prefix(
        base, "pns", mode, score_lower, score_upper
    )
    prefix = _locate_completed_prefix(
        base, direct_prefix, ("_pns.bw", "_posPNS.bw")
    )
    pns_path = Path(f"{prefix}_pns.bw").resolve()
    pos_pns_path = Path(f"{prefix}_posPNS.bw").resolve()
    scaled_path = (
        scaled_dir / f"{args.sample_name}_{label}_pns_divided_by_mean_posPNS.bw"
    ).resolve()
    reporter.stage(f"Scaling {label} PNS by mean posPNS for cluster aggregates")
    _scaled, reference_mean, reference_count = scale_bigwig_by_reference(
        pns_path, pos_pns_path, scaled_path, scale=1.0, reporter=reporter
    )
    report = setup_dir / f"{args.sample_name}_{label}_aggregate_pns_scaling.tsv"
    with report.open("wt", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"pns_track\t{pns_path}\n")
        handle.write(f"score_fragment_range\t{score_lower}-{score_upper}\n")
        handle.write(f"posPNS_reference_track\t{pos_pns_path}\n")
        handle.write(f"posPNS_reference_mean\t{reference_mean:.12g}\n")
        handle.write(f"finite_nonzero_posPNS_bases\t{reference_count}\n")
        handle.write(f"scaled_pns_track\t{scaled_path}\n")
    return {
        "pns": str(pns_path),
        "pos_pns": str(pos_pns_path),
        "scaled_pns": str(scaled_path),
        "pos_pns_mean": float(reference_mean),
        "pos_pns_finite_nonzero_bases": int(reference_count),
        "pns_scaling_report": str(report.resolve()),
        "score_frag_lower": score_lower,
        "score_frag_upper": score_upper,
    }


def _call_peaks(
    args: argparse.Namespace,
    reporter: ProgressReporter,
    *,
    scaled_path: Path,
    label: str,
    peaks_dir: Path,
    score_track: str,
) -> Path:
    peak_prefix = peaks_dir / f"{args.sample_name}_{label}_{score_track}_mean_scaled"
    command = [
        "call-peaks",
        "--input-bigwig", str(scaled_path),
        "--peak-caller", "pns",
        "--call-type", "nucleosome",
        "--out-prefix", str(peak_prefix),
        "--regions", *args.contigs,
        "--cores", str(args.cores),
        "--interval-format", "bed",
    ]
    if args.blacklist_bed:
        command.extend(["--blacklist-bed", args.blacklist_bed])
    reporter.stage(f"Calling {label} peaks on the mean scaled {score_track} track")
    _run_nucleosuite(command)
    completed = _locate_completed_prefix(
        peak_prefix, peak_prefix, ("_nucleosome_regions.bed",)
    )
    return Path(f"{completed}_nucleosome_regions.bed").resolve()


def _run_stage1(
    args: argparse.Namespace,
    *,
    outdir: Path,
    condition_name: str,
    treatment_bams: Sequence[str],
    control_bams: Sequence[str],
    target_mode: int,
    control_mode: int,
    mode_source: str,
    estimates: tuple[ModeEstimate | None, ModeEstimate | None, ModeEstimate | None],
    reporter: ProgressReporter,
    seed_offset: int,
    mode_condition_index: int,
) -> Path:
    tracks_dir = outdir / "01_score_tracks"
    scaled_dir = outdir / "02_mean_scaled_tracks"
    peaks_dir = outdir / "03_peak_calls"
    fdr_dir = outdir / "04_peak_fdr"
    aggregate_dir = outdir / "05_cluster_aggregate"
    setup_dir = outdir / "00_setup"
    for directory in (
        tracks_dir, scaled_dir, peaks_dir, fdr_dir, aggregate_dir, setup_dir
    ):
        directory.mkdir(parents=True, exist_ok=True)

    mode_report = setup_dir / f"{args.sample_name}_fragment_mode_estimation.tsv"
    _write_mode_report(
        mode_report,
        mode_source=mode_source,
        target=estimates[0],
        control=estimates[1],
        pooled=estimates[2],
        target_mode=target_mode,
        control_mode=control_mode,
        seeds=(
            args.seed + mode_condition_index * 10,
            args.seed + mode_condition_index * 10 + 1,
            args.seed + mode_condition_index * 10 + 2,
        ),
    )
    score_track, positive_track = _TRACK_NAMES[args.scoring_method]
    target_score_lower, target_score_upper = _scoring_fragment_range(
        args, target_mode
    )
    control_score_lower, control_score_upper = _scoring_fragment_range(
        args, control_mode
    )
    if args.bam_mode == "replicates":
        target_groups = [[path] for path in treatment_bams]
        control_groups = [[path] for path in control_bams]
    else:
        target_groups = [list(treatment_bams)]
        control_groups = [list(control_bams)]

    treatment_records: list[dict[str, object]] = []
    control_records: list[dict[str, object]] = []
    target_scaled: list[Path] = []
    control_scaled: list[Path] = []
    target_scaled_coverage: list[Path] = []
    control_scaled_coverage: list[Path] = []
    target_scaled_pns: list[Path] = []

    def generate_group(
        groups: Sequence[Sequence[str]],
        *,
        role: str,
        mode: int,
        seed_start: int,
    ) -> tuple[list[dict[str, object]], list[Path], list[Path], list[Path]]:
        records: list[dict[str, object]] = []
        scaled_scores: list[Path] = []
        scaled_coverages: list[Path] = []
        scaled_pns_tracks: list[Path] = []
        for index, bam_group in enumerate(groups, 1):
            suffix = "" if len(groups) == 1 else f"_rep{index}"
            generated = _generate_and_scale(
                args=args,
                reporter=reporter,
                bam_paths=bam_group,
                label=f"{role}{suffix}",
                mode=mode,
                tracks_dir=tracks_dir,
                scaled_dir=scaled_dir,
                setup_dir=setup_dir,
                score_track=score_track,
                positive_track=positive_track,
                seed=seed_start + index,
            )
            if role == "target":
                if args.scoring_method == "pns":
                    generated.update(
                        {
                            "pns": generated["score"],
                            "pos_pns": generated["positive_score"],
                            "scaled_pns": generated["scaled_score"],
                            "pos_pns_mean": generated["positive_score_mean"],
                            "pos_pns_finite_nonzero_bases": generated[
                                "positive_score_finite_nonzero_bases"
                            ],
                        }
                    )
                else:
                    generated.update(
                        _generate_pns_for_aggregate(
                            args=args,
                            reporter=reporter,
                            bam_paths=bam_group,
                            label=f"{role}{suffix}",
                            mode=mode,
                            tracks_dir=tracks_dir,
                            scaled_dir=scaled_dir,
                            setup_dir=setup_dir,
                            seed=seed_start + index + 10_000,
                        )
                    )
                scaled_pns_tracks.append(Path(str(generated["scaled_pns"])))
            scaled_scores.append(Path(str(generated["scaled_score"])))
            scaled_coverages.append(Path(str(generated["scaled_coverage"])))
            records.append(
                {
                    "replicate": index,
                    "bams": generated["bams"],
                    "score": generated["score"],
                    "positive_score": generated["positive_score"],
                    "coverage": generated["coverage"],
                    "scaled_score": generated["scaled_score"],
                    "scaled_coverage": generated["scaled_coverage"],
                    "positive_score_mean": generated["positive_score_mean"],
                    "positive_score_finite_nonzero_bases": generated[
                        "positive_score_finite_nonzero_bases"
                    ],
                    "coverage_nonzero_mean": generated["coverage_nonzero_mean"],
                    "coverage_finite_nonzero_bases": generated[
                        "coverage_finite_nonzero_bases"
                    ],
                    "score_frag_lower": generated["score_frag_lower"],
                    "score_frag_upper": generated["score_frag_upper"],
                    "coverage_frag_lower": generated["coverage_frag_lower"],
                    "coverage_frag_upper": generated["coverage_frag_upper"],
                    **(
                        {
                            "pns": generated["pns"],
                            "pos_pns": generated["pos_pns"],
                            "scaled_pns": generated["scaled_pns"],
                            "pos_pns_mean": generated["pos_pns_mean"],
                            "pos_pns_finite_nonzero_bases": generated[
                                "pos_pns_finite_nonzero_bases"
                            ],
                        }
                        if role == "target"
                        else {}
                    ),
                }
            )
        return records, scaled_scores, scaled_coverages, scaled_pns_tracks

    (
        treatment_records,
        target_scaled,
        target_scaled_coverage,
        target_scaled_pns,
    ) = generate_group(
        target_groups,
        role="target",
        mode=target_mode,
        seed_start=args.seed + seed_offset + 100,
    )
    control_records, control_scaled, control_scaled_coverage, _ = generate_group(
        control_groups,
        role="control",
        mode=control_mode,
        seed_start=args.seed + seed_offset + 200,
    )

    if len(target_scaled) == 1:
        mean_target = target_scaled[0]
        mean_target_coverage = target_scaled_coverage[0]
        mean_target_pns = target_scaled_pns[0]
    else:
        reporter.stage(f"Averaging {condition_name} replicate-scaled treatment tracks")
        mean_target = average_bigwigs(
            target_scaled,
            scaled_dir / f"{args.sample_name}_condition_mean_target_{score_track}.bw",
        )
        reporter.stage(f"Averaging {condition_name} replicate-scaled treatment coverage")
        mean_target_coverage = average_bigwigs(
            target_scaled_coverage,
            scaled_dir / f"{args.sample_name}_condition_mean_target_coverage_x100.bw",
        )
        if args.scoring_method == "pns":
            mean_target_pns = mean_target
        else:
            reporter.stage(
                f"Averaging {condition_name} replicate PNS tracks after independent posPNS scaling"
            )
            mean_target_pns = average_bigwigs(
                target_scaled_pns,
                scaled_dir
                / f"{args.sample_name}_condition_mean_target_pns_divided_by_mean_posPNS.bw",
            )

    if len(control_scaled) == 1:
        mean_control = control_scaled[0]
        mean_control_coverage = control_scaled_coverage[0]
    else:
        reporter.stage(f"Averaging {condition_name} replicate-scaled control tracks")
        mean_control = average_bigwigs(
            control_scaled,
            scaled_dir / f"{args.sample_name}_condition_mean_control_{score_track}.bw",
        )
        reporter.stage(f"Averaging {condition_name} replicate-scaled control coverage")
        mean_control_coverage = average_bigwigs(
            control_scaled_coverage,
            scaled_dir / f"{args.sample_name}_condition_mean_control_coverage_x100.bw",
        )

    target_peaks = _call_peaks(
        args, reporter, scaled_path=mean_target, label="target", peaks_dir=peaks_dir,
        score_track=score_track,
    )
    reporter.stage(
        f"Testing {condition_name} treatment candidates across all scaled-coverage replicates"
    )
    if args.bam_mode == "replicates" and (
        len(target_scaled_coverage) < 3 or len(control_scaled_coverage) < 3
    ):
        print(
            "WARNING: Stage 1 has fewer than three biological replicates in at "
            "least one group; Welch p-values and FDR are exploratory. Peak "
            "selection defaults to the all-controls gate.",
            file=sys.stderr,
        )
    outputs = analyze_chip_replicate_peaks(
        target_peaks,
        output_dir=fdr_dir,
        target_replicate_bigwigs=target_scaled_coverage,
        control_replicate_bigwigs=control_scaled_coverage,
        target_mean_bigwig=mean_target_coverage,
        peak_pvalue=args.stage1_p_value,
        peak_fdr=args.peak_fdr,
        cluster_seed_pvalue=args.cluster_seed_p_value,
        cluster_fdr=args.cluster_fdr,
        cluster_max_non_gated_gap=args.cluster_max_non_gated_gap,
        max_cluster_gap=args.max_cluster_gap,
        minimum_gate_peaks=args.min_cluster_gated_peaks,
    )

    anchor_bed = write_cluster_anchor_bed(
        outputs["selected_clusters"],
        aggregate_dir / f"{args.sample_name}_{condition_name}_cluster_anchors.bed",
    )
    if args.skip_cluster_aggregate:
        aggregate_outputs: dict[str, object] = {
            "status": "skipped_by_user",
            "anchor_bed": str(anchor_bed),
        }
    else:
        aggregate_outputs = run_cluster_aggregate(
            mean_scaled_pns=mean_target_pns,
            replicate_scaled_pns=target_scaled_pns,
            anchor_bed=anchor_bed,
            output_dir=aggregate_dir,
            label=condition_name,
            window_half=args.cluster_aggregate_window_half,
            maximum_heatmap_rows=args.cluster_aggregate_max_heatmap_rows,
            bootstrap_replicates=args.cluster_aggregate_bootstrap,
            nrl_peak_resolution=args.cluster_aggregate_nrl_resolution,
            nrl_min_order=args.cluster_aggregate_nrl_min_order,
            nrl_max_order=args.cluster_aggregate_nrl_max_order,
            seed=args.seed + seed_offset + 300,
            reporter=reporter,
        )

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "nucleosuite_version": __version__,
        "condition_name": condition_name,
        "bam_mode": args.bam_mode,
        "scoring_method": args.scoring_method,
        "score_track": score_track,
        "positive_track": positive_track,
        "mode_source": mode_source,
        "mode_histogram_smoothing": args.mode_histogram_smoothing,
        "mode_estimation_report": str(mode_report.resolve()),
        "target_mode": target_mode,
        "control_mode": control_mode,
        "frag_lower": target_score_lower,
        "frag_upper": target_score_upper,
        "score_fragment_flank": args.score_fragment_flank,
        "target_score_frag_lower": target_score_lower,
        "target_score_frag_upper": target_score_upper,
        "control_score_frag_lower": control_score_lower,
        "control_score_frag_upper": control_score_upper,
        "coverage_frag_lower": args.coverage_frag_lower,
        "coverage_frag_upper": args.coverage_frag_upper,
        "contigs": list(args.contigs),
        "blacklist_bed": str(Path(args.blacklist_bed).resolve()) if args.blacklist_bed else None,
        "max_duplicates": args.max_duplicates,
        "dedup_scope": args.dedup_scope,
        "stage1_p_value": args.stage1_p_value,
        "peak_fdr": args.peak_fdr,
        "cluster_seed_p_value": args.cluster_seed_p_value,
        "cluster_max_non_gated_gap": args.cluster_max_non_gated_gap,
        "minimum_cluster_gate_peaks": args.min_cluster_gated_peaks,
        "cluster_fdr": args.cluster_fdr,
        "stage1_selection": "all_treatments_exceed_all_controls",
        "stage1_statistics": "exploratory_one_sided_welch_bh_all_candidates",
        "treatment_replicates": treatment_records,
        "control_replicates": control_records,
        "condition_mean_treatment_score": str(mean_target.resolve()),
        "condition_mean_control_score": str(mean_control.resolve()),
        "condition_mean_treatment_coverage": str(mean_target_coverage.resolve()),
        "condition_mean_control_coverage": str(mean_control_coverage.resolve()),
        "condition_mean_treatment_pns_divided_by_mean_posPNS": str(
            mean_target_pns.resolve()
        ),
        "cluster_anchor_bed": str(anchor_bed),
        "cluster_aggregate": aggregate_outputs,
        "cluster_aggregate_parameters": {
            "window_half": args.cluster_aggregate_window_half,
            "maximum_heatmap_rows": args.cluster_aggregate_max_heatmap_rows,
            "bootstrap_replicates": args.cluster_aggregate_bootstrap,
            "nrl_peak_resolution": args.cluster_aggregate_nrl_resolution,
            "nrl_min_order": args.cluster_aggregate_nrl_min_order,
            "nrl_max_order": args.cluster_aggregate_nrl_max_order,
            "nrl_regression_exclusion": False,
        },
        "peak_discovery_track": score_track,
        "peak_measurement_track": "coverage_divided_by_nonzero_mean_x100",
        "target_candidate_peaks": str(target_peaks),
        "control_candidate_peaks": None,
        **{name: str(path.resolve()) for name, path in outputs.items()},
    }
    manifest_path = outdir / MANIFEST_NAME
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    summary = outdir / f"{args.sample_name}_chip_suite_summary.tsv"
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"condition_name\t{condition_name}\n")
        handle.write(f"bam_mode\t{args.bam_mode}\n")
        handle.write(f"scoring_method\t{args.scoring_method}\n")
        handle.write(f"mode_source\t{mode_source}\n")
        handle.write(f"mode_histogram_smoothing\t{args.mode_histogram_smoothing}\n")
        handle.write(f"mode_estimation_report\t{mode_report}\n")
        handle.write(f"target_analysis_mode\t{target_mode}\n")
        handle.write(f"control_analysis_mode\t{control_mode}\n")
        handle.write(
            f"treatment_score_fragment_range\t{target_score_lower}-"
            f"{target_score_upper}\n"
        )
        handle.write(
            f"control_score_fragment_range\t{control_score_lower}-"
            f"{control_score_upper}\n"
        )
        handle.write(
            f"coverage_fragment_range\t{args.coverage_frag_lower}-"
            f"{args.coverage_frag_upper}\n"
        )
        handle.write(
            "target_raw_coverage_track\t"
            + ";".join(str(record["coverage"]) for record in treatment_records)
            + "\n"
        )
        handle.write(
            "control_raw_coverage_track\t"
            + ";".join(str(record["coverage"]) for record in control_records)
            + "\n"
        )
        handle.write(f"condition_mean_treatment_score\t{mean_target}\n")
        handle.write(f"condition_mean_control_score\t{mean_control}\n")
        handle.write(f"condition_mean_treatment_coverage\t{mean_target_coverage}\n")
        handle.write(f"condition_mean_control_coverage\t{mean_control_coverage}\n")
        handle.write(f"condition_mean_treatment_scaled_pns\t{mean_target_pns}\n")
        handle.write(f"cluster_anchor_bed\t{anchor_bed}\n")
        handle.write(
            f"cluster_aggregate_status\t{aggregate_outputs.get('status')}\n"
        )
        handle.write(f"stage1_p_value\t{args.stage1_p_value}\n")
        handle.write(f"peak_fdr\t{args.peak_fdr}\n")
        handle.write(f"cluster_seed_p_value\t{args.cluster_seed_p_value}\n")
        handle.write(
            f"cluster_max_non_gated_gap\t{args.cluster_max_non_gated_gap}\n"
        )
        handle.write(
            f"minimum_cluster_gate_peaks\t{args.min_cluster_gated_peaks}\n"
        )
        handle.write(f"cluster_fdr\t{args.cluster_fdr}\n")
        for name, path in outputs.items():
            handle.write(f"{name}\t{path}\n")
        handle.write(f"stage1_manifest\t{manifest_path}\n")
    return manifest_path


def run(args: argparse.Namespace) -> int:
    _validate(args)
    outdir = Path(args.outdir).resolve()
    condition1_name = args.condition1_name or args.sample_name
    has_second = args.treatment2_bam is not None
    if args.dry_run:
        print(f"score_type\t{args.scoring_method}")
        print(f"mode\t{args.mode}")
        print(f"mode_histogram_smoothing\t{args.mode_histogram_smoothing}")
        print(f"bam_mode\t{args.bam_mode}")
        print(f"conditions\t{2 if has_second else 1}")
        print(
            "score_fragment_range\t"
            + (
                f"{args.score_frag_lower}-{args.score_frag_upper}"
                if args.score_frag_lower is not None
                else f"mode_plus_minus_{args.score_fragment_flank}"
            )
        )
        print(
            f"coverage_fragment_range\t{args.coverage_frag_lower}-"
            f"{args.coverage_frag_upper}"
        )
        print(
            "stages\tstage1-treatment-candidates-all-controls-gate"
            + (",stage2-log2-moderated-four-group-interaction-bh" if has_second else "")
        )
        return 0

    reporter = ProgressReporter("chip-suite")
    conditions: list[tuple[Sequence[str], Sequence[str]]] = [
        (args.treatment1_bam, args.control1_bam)
    ]
    if has_second:
        conditions.append((args.treatment2_bam, args.control2_bam))
    target_mode, control_mode, mode_source, estimates = _estimate_modes(
        args, conditions, reporter
    )

    first_outdir = outdir if not has_second else outdir / "01_condition1_stage1"
    first_manifest = _run_stage1(
        args,
        outdir=first_outdir,
        condition_name=condition1_name,
        treatment_bams=args.treatment1_bam,
        control_bams=args.control1_bam,
        target_mode=target_mode,
        control_mode=control_mode,
        mode_source=mode_source,
        estimates=estimates[0],
        reporter=reporter,
        seed_offset=1000,
        mode_condition_index=0,
    )
    if not has_second:
        print(f"chip_stage1_manifest\t{first_manifest}")
        return 0

    second_outdir = outdir / "02_condition2_stage1"
    second_manifest = _run_stage1(
        args,
        outdir=second_outdir,
        condition_name=args.condition2_name,
        treatment_bams=args.treatment2_bam,
        control_bams=args.control2_bam,
        target_mode=target_mode,
        control_mode=control_mode,
        mode_source=mode_source,
        estimates=estimates[1],
        reporter=reporter,
        seed_offset=2000,
        mode_condition_index=1,
    )
    reporter.stage("Comparing Stage 1 enrichments between conditions")
    comparison_manifest = compare_stage1(
        first_manifest,
        second_manifest,
        outdir=outdir / "03_condition_comparison",
        fdr=args.differential_fdr,
    )
    summary = outdir / f"{args.sample_name}_chip_suite_summary.tsv"
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"condition1_manifest\t{first_manifest}\n")
        handle.write(f"condition2_manifest\t{second_manifest}\n")
        handle.write(f"comparison_manifest\t{comparison_manifest}\n")
        handle.write(f"bam_mode\t{args.bam_mode}\n")
        handle.write(f"scoring_method\t{args.scoring_method}\n")
        handle.write(f"mode_histogram_smoothing\t{args.mode_histogram_smoothing}\n")
        handle.write(f"target_analysis_mode\t{target_mode}\n")
        handle.write(f"control_analysis_mode\t{control_mode}\n")
    print(f"chip_suite_summary\t{summary}")
    print(f"chip_comparison_manifest\t{comparison_manifest}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
