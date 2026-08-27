#!/usr/bin/env python3
"""Stage 1 and optional Stage 2 CUT&RUN/CUT&Tag workflow."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from nucleosuite import __version__
from nucleosuite.bigwig_ops import average_bigwigs
from nucleosuite.cutn_compare import (
    MANIFEST_NAME,
    MANIFEST_SCHEMA,
    MANIFEST_SCHEMA_VERSION,
    compare_stage1,
)
from nucleosuite.cutn_aggregate import run_cluster_aggregate, write_cluster_anchor_bed
from nucleosuite.cutn_peaks import analyze_cutn_replicate_peaks
from nucleosuite.mean_scale import scale_bigwig_by_reference
from nucleosuite.mode_estimation import (
    ModeEstimate,
    estimate_bam_fragment_mode,
    mode_estimate_message,
    pooled_mode_estimate,
)
from nucleosuite.progress import ProgressReporter


_TRACK_NAMES = {
    "sns": ("sns", "posSNS"),
    "tns": ("tns", "posTNS"),
    "bns": ("bns", "posBNS"),
    "pns": ("pns", "posPNS"),
}
_MULTICONTIG_MANIFEST = "nucleosuite_multicontig_manifest.json"
_RUN_MANIFEST = "cutn_suite_run_manifest.json"
_RUN_MANIFEST_SCHEMA = "nucleosuite-cutn-suite-run"
_RUN_MANIFEST_VERSION = 1


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
    lower = (
        int(args.score_frag_lower)
        if args.score_frag_lower is not None
        else max(1, int(mode) - int(args.frag_mode_padding))
    )
    upper = (
        int(args.score_frag_upper)
        if args.score_frag_upper is not None
        else int(mode) + int(args.frag_mode_padding)
    )
    return lower, upper


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


def _fragment_statistics(
    length_counts_path: Path,
    summary_path: Path | None = None,
    *,
    search_lower: int = 120,
    search_upper: int = 250,
) -> dict[str, int | None]:
    """Read lightweight per-sample fragment statistics from retained track outputs."""

    counts: dict[int, int] = {}
    if length_counts_path.is_file():
        try:
            with length_counts_path.open("rt", encoding="utf-8") as handle:
                next(handle, None)
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2:
                        continue
                    counts[int(fields[0])] = int(fields[1])
        except (OSError, ValueError):
            counts = {}

    mode = None
    mode_count = None
    candidates = {
        length: count
        for length, count in counts.items()
        if search_lower <= length <= search_upper
    }
    if candidates:
        mode = min(
            (length for length, count in candidates.items() if count == max(candidates.values())),
            default=None,
        )
        mode_count = candidates.get(mode) if mode is not None else None

    total_used = sum(counts.values()) if counts else None
    if summary_path is not None and summary_path.is_file():
        try:
            with summary_path.open("rt", encoding="utf-8") as handle:
                next(handle, None)
                for line in handle:
                    key, _, value = line.rstrip("\n").partition("\t")
                    if key == "total_fragments_used_in_range":
                        total_used = int(value)
                        break
        except (OSError, ValueError):
            pass

    return {
        "fragment_mode": mode,
        "fragment_mode_count": mode_count,
        "total_fragments_used": total_used,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite cutn-suite",
        description=(
            "Run control-aware Stage 1 CUT&RUN/CUT&Tag scoring for one "
            "condition and optionally compare it with a second condition."
        ),
    )
    parser.add_argument(
        "--treatment1-bam", dest="treatment1_bam", nargs="+",
        help="Condition 1 treatment BAM file(s). Required for a fresh run.",
    )
    parser.add_argument(
        "--control1-bam", dest="control1_bam", nargs="+",
        help="Condition 1 control BAM file(s). Required for a fresh run.",
    )
    parser.add_argument(
        "--treatment2-bam", nargs="+",
        help="Optional condition 2 treatment BAM file(s).",
    )
    parser.add_argument(
        "--control2-bam", nargs="+",
        help="Optional condition 2 control BAM file(s).",
    )
    parser.add_argument("--outdir", help="Output directory for a fresh run.")
    parser.add_argument("--sample-name", default="cutn", help="Output prefix (default: cutn).")
    parser.add_argument("--condition1-name", default="condition1", help="Condition 1 label (default: condition1).")
    parser.add_argument("--condition2-name", default="condition2", help="Condition 2 label (default: condition2).")
    management = parser.add_mutually_exclusive_group()
    management.add_argument(
        "--inspect-run",
        metavar="RUN_DIR",
        help="Inspect a completed cutn-suite run without reprocessing data.",
    )
    management.add_argument(
        "--rerun-from",
        metavar="RUN_DIR",
        help=(
            "Rerun downstream cutn-suite analysis from retained per-sample BigWigs. "
            "BAM-derived track settings are inherited and cannot be changed."
        ),
    )
    parser.add_argument(
        "--exclude-sample",
        action="append",
        default=[],
        metavar="BAM",
        help=(
            "With --rerun-from, exclude one biological replicate by BAM path, file "
            "name, or unambiguous stem. Repeat to exclude multiple samples."
        ),
    )
    parser.add_argument(
        "--bam-mode", choices=("replicates", "merged"), default="replicates",
        help=(
            "Treat BAMs as independent biological replicates, or pool each BAM "
            "group as one logical merged sample (default: replicates)."
        ),
    )
    parser.add_argument(
        "--scoring-method", choices=("sns", "pns", "bns", "tns"), default="sns",
        help="Nucleosome scoring kernel: sns, pns, bns, or tns (default: sns).",
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
        "--frag-mode-padding", dest="frag_mode_padding", type=int, default=30,
        help=(
            "When a scoring-fragment bound is omitted, derive it from the resolved "
            "mode plus or minus this many bp (default: 30). Explicit --score-frag-lower "
            "and --score-frag-upper values override their corresponding automatic "
            "bounds independently."
        ),
    )
    parser.add_argument(
        "--score-frag-lower", dest="score_frag_lower", type=int,
        help="Explicit scoring-fragment lower bound; overrides the automatic mode-minus-padding bound.",
    )
    parser.add_argument(
        "--score-frag-upper", dest="score_frag_upper", type=int,
        help="Explicit scoring-fragment upper bound; overrides the automatic mode-plus-padding bound.",
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
        "--peak-min-region-length", type=int, default=50,
        help="Minimum positive score-region length for treatment peak calls (default: 50 bp).",
    )
    parser.add_argument(
        "--peak-max-neg-run", type=int, default=0,
        help="Maximum zero-or-negative run bridged inside a treatment peak (default: 0 bp).",
    )
    parser.add_argument(
        "--peak-smooth-window", type=int, default=0,
        help="Savitzky-Golay window used only for treatment peak calling; 0 disables (default: 0).",
    )
    parser.add_argument(
        "--peak-smooth-order", type=int, default=2,
        help="Savitzky-Golay polynomial order for treatment peak calling (default: 2).",
    )
    parser.add_argument(
        "--stage1-p-value", type=float,
        help=(
            "Optional exploratory one-sided Stage 1 p-value cutoff. By default, "
            "the selected treatment-control gate alone selects peaks."
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
            "One-sided p-value required, together with the selected treatment-control gate, "
            "for a treatment peak to seed a cluster (default: 0.05)."
        ),
    )
    parser.add_argument("--differential-fdr", type=float, default=0.05, help="Stage 2 differential FDR cutoff (default: 0.05).")
    parser.add_argument(
        "--stage1-gate-mode", choices=("mean", "all-controls"), default="all-controls",
        help=(
            "Treatment-control gate: all-controls requires every treatment replicate to "
            "exceed every control replicate; mean compares group means (default: all-controls)."
        ),
    )
    parser.add_argument(
        "--cluster-member-mode", choices=("seed-and-gated", "significant-only"),
        default="seed-and-gated",
        help=(
            "Cluster S+G gate-passing members or only statistically significant S peaks "
            "(default: seed-and-gated)."
        ),
    )
    parser.add_argument(
        "--cluster-max-non-member-gap", type=int, default=1,
        help="Consecutive non-member candidates that may bridge included members (default: 1).",
    )
    parser.add_argument("--max-cluster-gap", type=int, default=1000, help="Maximum distance between adjacent gated-member summits (default: 1000 bp).")
    parser.add_argument(
        "--min-cluster-members", type=int, default=2,
        help="Minimum included members per seeded cluster (default: 2).",
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
        "--cluster-aggregate-nrl-resolution", type=float, default=130.0,
        help="Peak resolution for cluster-aligned directional NRLs (default: 130 bp).",
    )
    parser.add_argument(
        "--cluster-aggregate-nrl-min-order", type=int, default=0,
        help="First peak order used for cluster-aligned directional NRLs (default: 0).",
    )
    parser.add_argument(
        "--cluster-aggregate-nrl-max-order", type=int, default=3,
        help="Last peak order used for cluster-aligned directional NRLs (default: 3).",
    )
    aggregate_switch = parser.add_mutually_exclusive_group()
    aggregate_switch.add_argument(
        "--skip-cluster-aggregate", dest="skip_cluster_aggregate", action="store_true",
        help="Skip cluster-aligned score profiles, heatmaps and directional NRLs.",
    )
    aggregate_switch.add_argument(
        "--run-cluster-aggregate", dest="skip_cluster_aggregate", action="store_false",
        help=(
            "Run cluster aggregates explicitly. This is useful in a rerun when the "
            "source run skipped aggregate analysis."
        ),
    )
    parser.set_defaults(skip_cluster_aggregate=False)
    parser.add_argument("--dry-run", action="store_true", help="Validate and report planned stages without processing BAMs.")
    return parser


def _validate(args: argparse.Namespace) -> None:
    if args.inspect_run:
        if args.exclude_sample:
            raise ValueError("--exclude-sample is only valid with --rerun-from")
        if not Path(args.inspect_run).exists():
            raise FileNotFoundError(args.inspect_run)
        return

    if args.rerun_from and not Path(args.rerun_from).exists():
        raise FileNotFoundError(args.rerun_from)

    if (
        args.coverage_frag_lower < 1
        or args.coverage_frag_upper < args.coverage_frag_lower
    ):
        raise ValueError(
            "Require 1 <= --coverage-frag-lower <= --coverage-frag-upper"
        )
    if args.frag_mode_padding < 0:
        raise ValueError("--frag-mode-padding must be non-negative")
    if args.score_frag_lower is not None and args.score_frag_lower < 1:
        raise ValueError("--score-frag-lower must be at least 1")
    if args.score_frag_upper is not None and args.score_frag_upper < 1:
        raise ValueError("--score-frag-upper must be at least 1")
    if (
        args.score_frag_lower is not None
        and args.score_frag_upper is not None
        and args.score_frag_upper < args.score_frag_lower
    ):
        raise ValueError("Require --score-frag-upper >= --score-frag-lower")
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
    if args.peak_min_region_length < 1 or args.peak_max_neg_run < 0:
        raise ValueError("Peak-region length must be positive and negative-run allowance non-negative")
    if args.peak_smooth_window < 0 or (
        args.peak_smooth_window
        and (args.peak_smooth_window < 3 or args.peak_smooth_window % 2 == 0)
    ):
        raise ValueError("--peak-smooth-window must be 0 or an odd integer of at least 3")
    if args.peak_smooth_order < 0 or (
        args.peak_smooth_window and args.peak_smooth_order >= args.peak_smooth_window
    ):
        raise ValueError("--peak-smooth-order must be non-negative and smaller than the active window")
    if args.cluster_max_non_member_gap < 0:
        raise ValueError("--cluster-max-non-member-gap must be non-negative")
    if min(args.max_cluster_gap, args.min_cluster_members) < 1:
        raise ValueError("Cluster gap and minimum member count must be positive")
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

    if args.rerun_from:
        if args.outdir is not None:
            raise ValueError("--outdir is not used with --rerun-from; reruns are written inside the source run")
        return

    if args.exclude_sample:
        raise ValueError("--exclude-sample requires --rerun-from")
    if args.treatment1_bam is None or args.control1_bam is None or args.outdir is None:
        raise ValueError(
            "A fresh cutn-suite run requires --treatment1-bam, --control1-bam, and --outdir"
        )
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
        f"Generating {label} {score_track}/{positive_track} from {score_lower}-{score_upper} bp "
        f"and coverage from {args.coverage_frag_lower}-{args.coverage_frag_upper} bp in one tracks pass"
    )
    score_base = (tracks_dir / f"{args.sample_name}_{label}_{args.scoring_method}_discovery").resolve()
    coverage_base = (tracks_dir / f"{args.sample_name}_{label}_coverage_all_fragments").resolve()
    spec_path = setup_dir / f"{args.sample_name}_{label}_tracks_spec.tsv"
    with spec_path.open("wt", encoding="utf-8") as handle:
        handle.write("fragment_range\toutput_prefix\ttracks\tbasic_scope\n")
        handle.write(
            f"{score_lower}-{score_upper}\t{score_base}\t"
            f"{score_track},{positive_track}\trange\n"
        )
        handle.write(
            f"{args.coverage_frag_lower}-{args.coverage_frag_upper}\t"
            f"{coverage_base}\tcoverage\trange\n"
        )
    command = [
        "tracks", "--bam", *bam_paths,
        "--spec-file", str(spec_path),
        "--output-dir", str(tracks_dir.resolve()),
        "--scoring-method", args.scoring_method,
        "--score-mode-length", str(mode),
        "--output-format", "bigwig",
        "--interval-format", "bed",
        "--contigs", *args.contigs,
        "--max-duplicates", str(args.max_duplicates),
        "--dedup-scope", args.dedup_scope,
        "--cores", str(args.cores),
        "--seed", str(seed),
    ]
    if args.blacklist_bed:
        command.extend(["--blacklist-bed", args.blacklist_bed])
    _run_nucleosuite(command)
    score_prefix = _locate_completed_prefix(
        score_base, score_base, (f"_{score_track}.bw", f"_{positive_track}.bw")
    )
    coverage_prefix = _locate_completed_prefix(
        coverage_base, coverage_base, ("_coverage.bw",)
    )
    score_path = Path(f"{score_prefix}_{score_track}.bw").resolve()
    positive_path = Path(f"{score_prefix}_{positive_track}.bw").resolve()
    coverage_path = Path(f"{coverage_prefix}_coverage.bw").resolve()
    scaled_path = (scaled_dir / f"{args.sample_name}_{label}_{score_track}_divided_by_mean_{positive_track}.bw").resolve()
    reporter.stage(f"Scaling {label} {score_track} by mean {positive_track}")
    _scaled, reference_mean, reference_count = scale_bigwig_by_reference(
        score_path, positive_path, scaled_path, scale=1.0, reporter=reporter
    )
    scaled_coverage_path = (scaled_dir / f"{args.sample_name}_{label}_coverage_divided_by_mean_x100.bw").resolve()
    reporter.stage(f"Scaling {label} coverage to a non-zero mean of 100")
    _coverage_scaled, coverage_mean, coverage_count = scale_bigwig_by_reference(
        coverage_path, coverage_path, scaled_coverage_path, scale=100.0, reporter=reporter
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
    fragment_summary = Path(f"{coverage_prefix}_fragment_summary.tsv").resolve()
    fragment_length_counts = Path(f"{coverage_prefix}_fragment_length_counts.tsv").resolve()
    fragment_stats = _fragment_statistics(
        fragment_length_counts, fragment_summary,
        search_lower=args.mode_search_lower, search_upper=args.mode_search_upper,
    )
    return {
        "bams": [str(Path(path).resolve()) for path in bam_paths],
        "score": str(score_path), "positive_score": str(positive_path), "coverage": str(coverage_path),
        "scaled_score": str(scaled_path), "scaled_coverage": str(scaled_coverage_path),
        "positive_score_mean": float(reference_mean),
        "positive_score_finite_nonzero_bases": int(reference_count),
        "coverage_nonzero_mean": float(coverage_mean), "coverage_finite_nonzero_bases": int(coverage_count),
        "score_frag_lower": score_lower, "score_frag_upper": score_upper,
        "coverage_frag_lower": int(args.coverage_frag_lower), "coverage_frag_upper": int(args.coverage_frag_upper),
        "fragment_summary": str(fragment_summary),
        "fragment_length_counts": str(fragment_length_counts),
        "sample_fragment_mode": fragment_stats.get("fragment_mode"),
        "sample_fragment_mode_count": fragment_stats.get("fragment_mode_count"),
        "total_fragments_used": fragment_stats.get("total_fragments_used"),
        "tracks_spec": str(spec_path.resolve()),
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
        "--min-region-length", str(args.peak_min_region_length),
        "--max-neg-run", str(args.peak_max_neg_run),
        "--smooth-window", str(args.peak_smooth_window),
        "--smooth-order", str(args.peak_smooth_order),
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

    def generate_group(
        groups: Sequence[Sequence[str]],
        *,
        role: str,
        mode: int,
        seed_start: int,
    ) -> tuple[list[dict[str, object]], list[Path], list[Path]]:
        records: list[dict[str, object]] = []
        scaled_scores: list[Path] = []
        scaled_coverages: list[Path] = []
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
                    "fragment_summary": generated.get("fragment_summary"),
                    "fragment_length_counts": generated.get("fragment_length_counts"),
                    "sample_fragment_mode": generated.get("sample_fragment_mode"),
                    "sample_fragment_mode_count": generated.get("sample_fragment_mode_count"),
                    "total_fragments_used": generated.get("total_fragments_used"),
                }
            )
        return records, scaled_scores, scaled_coverages

    (
        treatment_records,
        target_scaled,
        target_scaled_coverage,
    ) = generate_group(
        target_groups,
        role="target",
        mode=target_mode,
        seed_start=args.seed + seed_offset + 100,
    )
    control_records, control_scaled, control_scaled_coverage = generate_group(
        control_groups,
        role="control",
        mode=control_mode,
        seed_start=args.seed + seed_offset + 200,
    )

    if len(target_scaled) == 1:
        mean_target = target_scaled[0]
        mean_target_coverage = target_scaled_coverage[0]
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
        gate_description = (
            "every treatment replicate > every control replicate"
            if args.stage1_gate_mode == "all-controls"
            else "mean treatment > mean control"
        )
        print(
            "WARNING: Stage 1 has fewer than three biological replicates in at "
            "least one group; Welch p-values and FDR are exploratory. Peak "
            f"selection uses the {args.stage1_gate_mode} gate ({gate_description}).",
            file=sys.stderr,
        )
    outputs = analyze_cutn_replicate_peaks(
        target_peaks,
        output_dir=fdr_dir,
        target_replicate_bigwigs=target_scaled_coverage,
        control_replicate_bigwigs=control_scaled_coverage,
        target_mean_bigwig=mean_target_coverage,
        peak_pvalue=args.stage1_p_value,
        peak_fdr=args.peak_fdr,
        cluster_seed_pvalue=args.cluster_seed_p_value,
        cluster_fdr=args.cluster_fdr,
        gate_mode=args.stage1_gate_mode,
        cluster_member_mode=args.cluster_member_mode,
        cluster_max_non_member_gap=args.cluster_max_non_member_gap,
        max_cluster_gap=args.max_cluster_gap,
        minimum_cluster_members=args.min_cluster_members,
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
            mean_scaled_score=mean_target,
            replicate_scaled_scores=target_scaled,
            anchor_bed=anchor_bed,
            output_dir=aggregate_dir,
            label=condition_name,
            scoring_method=args.scoring_method,
            positive_track=positive_track,
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
        "sample_name": args.sample_name,
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
        "frag_mode_padding": args.frag_mode_padding,
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
        "peak_min_region_length": args.peak_min_region_length,
        "peak_max_neg_run": args.peak_max_neg_run,
        "peak_smooth_window": args.peak_smooth_window,
        "peak_smooth_order": args.peak_smooth_order,
        "stage1_p_value": args.stage1_p_value,
        "peak_fdr": args.peak_fdr,
        "cluster_seed_p_value": args.cluster_seed_p_value,
        "cluster_max_non_member_gap": args.cluster_max_non_member_gap,
        "max_cluster_gap": args.max_cluster_gap,
        "minimum_cluster_members": args.min_cluster_members,
        "cluster_fdr": args.cluster_fdr,
        "stage1_gate_mode": args.stage1_gate_mode,
        "cluster_member_mode": args.cluster_member_mode,
        "stage1_selection": ("mean_treatment_exceeds_mean_control" if args.stage1_gate_mode == "mean" else "all_treatments_exceed_all_controls"),
        "stage1_statistics": "exploratory_one_sided_welch_bh_all_candidates",
        "treatment_replicates": treatment_records,
        "control_replicates": control_records,
        "condition_mean_treatment_score": str(mean_target.resolve()),
        "condition_mean_control_score": str(mean_control.resolve()),
        "condition_mean_treatment_coverage": str(mean_target_coverage.resolve()),
        "condition_mean_control_coverage": str(mean_control_coverage.resolve()),
        "condition_mean_treatment_cluster_aggregate_score": str(mean_target.resolve()),
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

    summary = outdir / f"{args.sample_name}_cutn_suite_summary.tsv"
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
        handle.write(f"condition_mean_treatment_cluster_aggregate_score\t{mean_target}\n")
        handle.write(f"cluster_anchor_bed\t{anchor_bed}\n")
        handle.write(
            f"cluster_aggregate_status\t{aggregate_outputs.get('status')}\n"
        )
        handle.write(f"peak_min_region_length\t{args.peak_min_region_length}\n")
        handle.write(f"peak_max_neg_run\t{args.peak_max_neg_run}\n")
        handle.write(f"peak_smooth_window\t{args.peak_smooth_window}\n")
        handle.write(f"peak_smooth_order\t{args.peak_smooth_order}\n")
        handle.write(f"stage1_p_value\t{args.stage1_p_value}\n")
        handle.write(f"peak_fdr\t{args.peak_fdr}\n")
        handle.write(f"cluster_seed_p_value\t{args.cluster_seed_p_value}\n")
        handle.write(f"stage1_gate_mode\t{args.stage1_gate_mode}\n")
        handle.write(f"cluster_member_mode\t{args.cluster_member_mode}\n")
        handle.write(
            f"cluster_max_non_member_gap\t{args.cluster_max_non_member_gap}\n"
        )
        handle.write(f"max_cluster_gap\t{args.max_cluster_gap}\n")
        handle.write(
            f"minimum_cluster_members\t{args.min_cluster_members}\n"
        )
        handle.write(f"cluster_fdr\t{args.cluster_fdr}\n")
        for name, path in outputs.items():
            handle.write(f"{name}\t{path}\n")
        handle.write(f"stage1_manifest\t{manifest_path}\n")
    return manifest_path


def _read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read cutn-suite manifest: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Invalid cutn-suite manifest: {path}")
    return value


def _resolve_run_reference(root: Path, value: object) -> Path:
    path = Path(str(value))
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def _discover_run(value: str | Path) -> tuple[Path, list[tuple[Path, dict[str, object]]], Path | None]:
    """Locate one or two Stage 1 manifests in a cutn-suite run directory."""

    supplied = Path(value).resolve()
    if supplied.is_file():
        if supplied.name == _RUN_MANIFEST:
            root = supplied.parent
            run_manifest_path = supplied
        elif supplied.name == MANIFEST_NAME:
            return supplied.parent, [(supplied, _read_json(supplied))], None
        else:
            raise ValueError(
                f"Expected {_RUN_MANIFEST}, {MANIFEST_NAME}, or a cutn-suite run directory: {supplied}"
            )
    else:
        root = supplied
        run_manifest_path = root / _RUN_MANIFEST

    if run_manifest_path.is_file():
        run_manifest = _read_json(run_manifest_path)
        refs = run_manifest.get("condition_manifests")
        if not isinstance(refs, list) or not refs:
            conditions = run_manifest.get("conditions")
            if isinstance(conditions, list):
                refs = [
                    item.get("manifest")
                    for item in conditions
                    if isinstance(item, dict) and item.get("manifest")
                ]
        if not isinstance(refs, list) or not refs:
            raise ValueError(f"Run manifest has no condition manifests: {run_manifest_path}")
        found: list[tuple[Path, dict[str, object]]] = []
        for ref in refs:
            manifest_path = _resolve_run_reference(root, ref)
            if not manifest_path.is_file():
                raise FileNotFoundError(manifest_path)
            found.append((manifest_path, _read_json(manifest_path)))
        return root, found, run_manifest_path

    direct = root / MANIFEST_NAME
    if direct.is_file():
        return root, [(direct, _read_json(direct))], None

    known = [
        root / "01_condition1_stage1" / MANIFEST_NAME,
        root / "02_condition2_stage1" / MANIFEST_NAME,
    ]
    found_paths = [path for path in known if path.is_file()]
    if not found_paths:
        found_paths = sorted(
            path
            for path in root.glob(f"*_condition*_stage1/{MANIFEST_NAME}")
            if path.is_file()
        )
    if not found_paths:
        raise FileNotFoundError(
            f"No cutn-suite Stage 1 manifest was found under {root}"
        )
    return root, [(path, _read_json(path)) for path in found_paths], None


def _infer_sample_name(manifest: dict[str, object]) -> str:
    value = manifest.get("sample_name")
    if value:
        return str(value)
    for key, token in (
        ("condition_mean_treatment_score", "_condition_mean_target_"),
        ("target_candidate_peaks", "_target_"),
    ):
        raw = manifest.get(key)
        if raw:
            name = Path(str(raw)).name
            if token in name:
                prefix = name.split(token, 1)[0]
                if prefix:
                    return prefix
    return "cutn"


def _mode_search_bounds_from_manifest(manifest: dict[str, object]) -> tuple[int, int]:
    report = manifest.get("mode_estimation_report")
    if report:
        path = Path(str(report))
        if path.is_file():
            try:
                with path.open("rt", encoding="utf-8") as handle:
                    header = next(handle).rstrip("\n").split("\t")
                    lower_i = header.index("mode_search_lower")
                    upper_i = header.index("mode_search_upper")
                    for line in handle:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) > max(lower_i, upper_i) and fields[lower_i] and fields[upper_i]:
                            return int(fields[lower_i]), int(fields[upper_i])
            except (OSError, ValueError, StopIteration):
                pass
    return 120, 250


def _record_fragment_statistics(
    record: dict[str, object], manifest: dict[str, object]
) -> dict[str, int | None]:
    stored = {
        "fragment_mode": record.get("sample_fragment_mode"),
        "fragment_mode_count": record.get("sample_fragment_mode_count"),
        "total_fragments_used": record.get("total_fragments_used"),
    }
    if all(value is not None for value in stored.values()):
        return {key: int(value) if value is not None else None for key, value in stored.items()}

    counts_raw = record.get("fragment_length_counts")
    summary_raw = record.get("fragment_summary")
    if not counts_raw:
        coverage_raw = record.get("coverage")
        if coverage_raw and str(coverage_raw).endswith("_coverage.bw"):
            prefix = str(coverage_raw)[: -len("_coverage.bw")]
            counts_raw = prefix + "_fragment_length_counts.tsv"
            summary_raw = prefix + "_fragment_summary.tsv"
    lower, upper = _mode_search_bounds_from_manifest(manifest)
    inferred = _fragment_statistics(
        Path(str(counts_raw)) if counts_raw else Path("/__missing_fragment_counts__"),
        Path(str(summary_raw)) if summary_raw else None,
        search_lower=lower,
        search_upper=upper,
    )
    for key, value in stored.items():
        if value is not None:
            inferred[key] = int(value)
    return inferred


def _mode_report_rows(manifest: dict[str, object]) -> dict[str, dict[str, str]]:
    path_raw = manifest.get("mode_estimation_report")
    if not path_raw:
        return {}
    path = Path(str(path_raw))
    if not path.is_file():
        return {}
    try:
        with path.open("rt", encoding="utf-8") as handle:
            header = next(handle).rstrip("\n").split("\t")
            rows = {}
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                fields += [""] * max(0, len(header) - len(fields))
                row = dict(zip(header, fields))
                if row.get("sample"):
                    rows[row["sample"]] = row
            return rows
    except (OSError, StopIteration):
        return {}


def inspect_run(value: str | Path) -> int:
    root, conditions, run_manifest_path = _discover_run(value)
    run_meta = _read_json(run_manifest_path) if run_manifest_path else {}
    version = run_meta.get("nucleosuite_version")
    if not version and conditions:
        version = conditions[0][1].get("nucleosuite_version")
    print("NucleoSuite cutn-suite run")
    print(f"run_directory\t{root}")
    print(f"nucleosuite_version\t{version or 'unknown'}")
    print(f"conditions\t{len(conditions)}")
    if run_manifest_path:
        print(f"run_manifest\t{run_manifest_path}")
    if run_meta.get("source_run"):
        print(f"source_run\t{run_meta['source_run']}")
    excluded = run_meta.get("excluded_bams")
    if isinstance(excluded, list) and excluded:
        print("excluded_bams\t" + ";".join(map(str, excluded)))

    for condition_index, (manifest_path, manifest) in enumerate(conditions, 1):
        name = str(manifest.get("condition_name") or f"condition{condition_index}")
        print(f"\nCondition {condition_index}: {name}")
        print(f"  manifest: {manifest_path}")
        print(f"  scoring method: {manifest.get('scoring_method', 'sns')}")
        print(
            "  analysis modes: "
            f"treatment={manifest.get('target_mode', 'unknown')} bp; "
            f"control={manifest.get('control_mode', 'unknown')} bp"
        )
        print(
            "  scoring ranges: "
            f"treatment={manifest.get('target_score_frag_lower', manifest.get('frag_lower', '?'))}-"
            f"{manifest.get('target_score_frag_upper', manifest.get('frag_upper', '?'))} bp; "
            f"control={manifest.get('control_score_frag_lower', '?')}-"
            f"{manifest.get('control_score_frag_upper', '?')} bp"
        )
        print(
            "  coverage range: "
            f"{manifest.get('coverage_frag_lower', '?')}-{manifest.get('coverage_frag_upper', '?')} bp"
        )
        print(
            "  Stage 1: "
            f"gate={manifest.get('stage1_gate_mode', 'all-controls')}; "
            f"cluster members={manifest.get('cluster_member_mode', 'seed-and-gated')}; "
            f"seed p<{manifest.get('cluster_seed_p_value', 0.05)}"
        )
        report_rows = _mode_report_rows(manifest)
        for role_key, role_label, report_key in (
            ("treatment_replicates", "Treatment", "treatment"),
            ("control_replicates", "Control", "control"),
        ):
            records = manifest.get(role_key)
            records = records if isinstance(records, list) else []
            group_report = report_rows.get(report_key, {})
            group_text = ""
            if group_report.get("mode"):
                group_text = f"; group mode={group_report['mode']} bp"
                if group_report.get("bootstrap_ci_low") and group_report.get("bootstrap_ci_high"):
                    group_text += (
                        f" (95% bootstrap CI {group_report['bootstrap_ci_low']}-"
                        f"{group_report['bootstrap_ci_high']})"
                    )
            print(f"  {role_label} ({len(records)} retained replicate track(s)){group_text}:")
            for record_index, raw_record in enumerate(records, 1):
                if not isinstance(raw_record, dict):
                    continue
                record = raw_record
                bams = record.get("bams")
                bams = bams if isinstance(bams, list) else []
                stats = _record_fragment_statistics(record, manifest)
                mode = stats.get("fragment_mode")
                used = stats.get("total_fragments_used")
                print(f"    replicate {record_index}:")
                for bam in bams:
                    print(f"      BAM: {bam}")
                print(f"      fragment mode: {mode if mode is not None else 'unavailable'}" + (" bp" if mode is not None else ""))
                print(f"      fragments used (coverage range): {used if used is not None else 'unavailable'}")
                if record.get("positive_score_mean") is not None:
                    print(f"      positive-score mean: {record['positive_score_mean']}")
                if record.get("coverage_nonzero_mean") is not None:
                    print(f"      coverage non-zero mean before x100 scaling: {record['coverage_nonzero_mean']}")
    return 0


def _explicit(args: argparse.Namespace, *flags: str) -> bool:
    supplied = getattr(args, "_explicit_options", set())
    return any(flag in supplied for flag in flags)


def _inherit_rerun_parameters(args: argparse.Namespace, manifest: dict[str, object]) -> dict[str, object]:
    """Populate downstream parameters from the source run unless explicitly overridden."""

    changed: dict[str, object] = {}
    scalar_map = {
        "peak_min_region_length": ("--peak-min-region-length", "peak_min_region_length", 50),
        "peak_max_neg_run": ("--peak-max-neg-run", "peak_max_neg_run", 0),
        "peak_smooth_window": ("--peak-smooth-window", "peak_smooth_window", 0),
        "peak_smooth_order": ("--peak-smooth-order", "peak_smooth_order", 2),
        "stage1_p_value": ("--stage1-p-value", "stage1_p_value", None),
        "peak_fdr": ("--peak-fdr", "peak_fdr", None),
        "cluster_fdr": ("--cluster-fdr", "cluster_fdr", None),
        "cluster_seed_p_value": ("--cluster-seed-p-value", "cluster_seed_p_value", 0.05),
        "stage1_gate_mode": ("--stage1-gate-mode", "stage1_gate_mode", "all-controls"),
        "cluster_member_mode": ("--cluster-member-mode", "cluster_member_mode", "seed-and-gated"),
        "cluster_max_non_member_gap": ("--cluster-max-non-member-gap", "cluster_max_non_member_gap", 1),
        "max_cluster_gap": ("--max-cluster-gap", "max_cluster_gap", 1000),
        "min_cluster_members": ("--min-cluster-members", "minimum_cluster_members", 2),
    }
    for attr, (flag, field, fallback) in scalar_map.items():
        source_value = manifest.get(field, fallback)
        if _explicit(args, flag):
            changed[attr] = getattr(args, attr)
        else:
            setattr(args, attr, source_value)

    aggregate = manifest.get("cluster_aggregate_parameters")
    aggregate = aggregate if isinstance(aggregate, dict) else {}
    aggregate_map = {
        "cluster_aggregate_window_half": ("--cluster-aggregate-window-half", "window_half", 1000),
        "cluster_aggregate_max_heatmap_rows": ("--cluster-aggregate-max-heatmap-rows", "maximum_heatmap_rows", 5000),
        "cluster_aggregate_bootstrap": ("--cluster-aggregate-bootstrap", "bootstrap_replicates", 200),
        "cluster_aggregate_nrl_resolution": ("--cluster-aggregate-nrl-resolution", "nrl_peak_resolution", 130.0),
        "cluster_aggregate_nrl_min_order": ("--cluster-aggregate-nrl-min-order", "nrl_min_order", 0),
        "cluster_aggregate_nrl_max_order": ("--cluster-aggregate-nrl-max-order", "nrl_max_order", 3),
    }
    for attr, (flag, field, fallback) in aggregate_map.items():
        source_value = aggregate.get(field, fallback)
        if _explicit(args, flag):
            changed[attr] = getattr(args, attr)
        else:
            setattr(args, attr, source_value)

    if not _explicit(args, "--skip-cluster-aggregate", "--run-cluster-aggregate"):
        aggregate_result = manifest.get("cluster_aggregate")
        args.skip_cluster_aggregate = bool(
            isinstance(aggregate_result, dict)
            and aggregate_result.get("status") == "skipped_by_user"
        )
    else:
        changed["skip_cluster_aggregate"] = args.skip_cluster_aggregate

    if not _explicit(args, "--sample-name"):
        args.sample_name = _infer_sample_name(manifest)
    args.bam_mode = str(manifest.get("bam_mode") or "replicates")
    args.scoring_method = str(manifest.get("scoring_method") or "sns")
    args.frag_mode_padding = int(manifest.get("frag_mode_padding", manifest.get("score_fragment_flank", 30)))
    args.score_frag_lower = manifest.get("target_score_frag_lower", manifest.get("frag_lower"))
    args.score_frag_upper = manifest.get("target_score_frag_upper", manifest.get("frag_upper"))
    args.coverage_frag_lower = int(manifest.get("coverage_frag_lower", 1))
    args.coverage_frag_upper = int(manifest.get("coverage_frag_upper", 1000))
    args.contigs = list(manifest.get("contigs") or ["autosomes"])
    args.blacklist_bed = manifest.get("blacklist_bed")
    args.max_duplicates = int(manifest.get("max_duplicates", 1))
    args.dedup_scope = str(manifest.get("dedup_scope") or "all_bams")
    args.mode_histogram_smoothing = str(
        manifest.get("mode_histogram_smoothing") or "none"
    )
    return changed


def _validate_rerun_immutable_options(args: argparse.Namespace) -> None:
    immutable = {
        "--treatment1-bam", "--control1-bam", "--treatment2-bam", "--control2-bam",
        "--bam-mode", "--scoring-method", "--mode", "--mode-strategy",
        "--frag-mode-padding", "--score-frag-lower", "--score-frag-upper",
        "--coverage-frag-lower", "--coverage-frag-upper", "--mode-search-lower",
        "--mode-search-upper", "--mode-min-fragments", "--mode-batch-fragments",
        "--mode-max-fragments", "--mode-bootstrap", "--mode-stable-checkpoints",
        "--mode-max-change", "--mode-max-ci-width", "--mode-block-bp",
        "--mode-histogram-smoothing", "--seed", "--blacklist-bed", "--contigs", "-c",
        "--max-duplicates", "--dedup-scope",
    }
    supplied = getattr(args, "_explicit_options", set())
    blocked = sorted(immutable.intersection(supplied))
    if blocked:
        raise ValueError(
            "These options change BAM-derived per-sample BigWigs and cannot be changed "
            "by --rerun-from: " + ", ".join(blocked)
        )


def _all_replicate_entries(
    conditions: list[tuple[Path, dict[str, object]]]
) -> list[tuple[int, str, int, dict[str, object]]]:
    entries = []
    for condition_index, (_path, manifest) in enumerate(conditions, 1):
        for role_key, role in (("treatment_replicates", "treatment"), ("control_replicates", "control")):
            records = manifest.get(role_key)
            records = records if isinstance(records, list) else []
            for record_index, record in enumerate(records):
                if isinstance(record, dict):
                    entries.append((condition_index, role, record_index, record))
    return entries


def _match_excluded_records(
    conditions: list[tuple[Path, dict[str, object]]], tokens: Sequence[str]
) -> tuple[set[tuple[int, str, int]], list[str]]:
    excluded: set[tuple[int, str, int]] = set()
    matched_bams: list[str] = []
    entries = _all_replicate_entries(conditions)
    for token in tokens:
        token_path = Path(token)
        token_name = token_path.name
        token_stem = token_path.stem
        exact_path = str(token_path.resolve()) if token_path.is_absolute() else None
        matches: list[tuple[int, str, int, dict[str, object], str]] = []
        for condition_index, role, record_index, record in entries:
            bams = record.get("bams")
            bams = bams if isinstance(bams, list) else []
            for bam in bams:
                bam_str = str(bam)
                bam_path = Path(bam_str)
                if (
                    (exact_path is not None and str(bam_path.resolve()) == exact_path)
                    or bam_str == token
                    or bam_path.name == token_name
                    or bam_path.stem == token_stem
                ):
                    matches.append((condition_index, role, record_index, record, bam_str))
        unique = {(c, r, i): (rec, bam) for c, r, i, rec, bam in matches}
        if not unique:
            raise ValueError(f"Excluded sample did not match any retained BAM: {token}")
        if len(unique) > 1:
            raise ValueError(
                f"Excluded sample name is ambiguous across retained replicates: {token}. "
                "Use the full BAM path."
            )
        (condition_index, role, record_index), (record, bam_str) = next(iter(unique.items()))
        bams = record.get("bams")
        bams = bams if isinstance(bams, list) else []
        if len(bams) != 1:
            raise ValueError(
                "Cannot exclude one BAM from a retained merged replicate track. "
                "Rerun the original BAM processing or use --bam-mode replicates in the source run."
            )
        excluded.add((condition_index, role, record_index))
        matched_bams.append(bam_str)
    return excluded, matched_bams


def _filtered_records(
    manifest: dict[str, object], condition_index: int, excluded: set[tuple[int, str, int]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    result = []
    for role_key, role in (("treatment_replicates", "treatment"), ("control_replicates", "control")):
        records = manifest.get(role_key)
        records = records if isinstance(records, list) else []
        kept = [
            dict(record)
            for index, record in enumerate(records)
            if isinstance(record, dict) and (condition_index, role, index) not in excluded
        ]
        if not kept:
            raise ValueError(
                f"Exclusions leave condition {condition_index} with no {role} replicates"
            )
        result.append(kept)
    return result[0], result[1]


def _safe_rerun_token(value: str) -> str:
    stem = Path(value).stem or Path(value).name or "sample"
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._-") or "sample"
    return cleaned[:60]


def _next_rerun_directory(root: Path, exclusions: Sequence[str]) -> Path:
    if len(exclusions) == 1:
        base = f"rerun_excluding_{_safe_rerun_token(exclusions[0])}"
    elif len(exclusions) > 1:
        base = f"rerun_excluding_{len(exclusions)}_samples"
    else:
        base = "rerun"
    pattern = re.compile(rf"^{re.escape(base)}_(\d+)$")
    existing = []
    for child in root.iterdir():
        if not child.is_dir():
            continue
        match = pattern.match(child.name)
        if match:
            existing.append(int(match.group(1)))
    number = max(existing, default=0) + 1
    return root / f"{base}_{number:02d}"


def _warn_small_replicate_groups(
    args: argparse.Namespace, target_scaled_coverage: Sequence[Path], control_scaled_coverage: Sequence[Path]
) -> None:
    if args.bam_mode == "replicates" and (
        len(target_scaled_coverage) < 3 or len(control_scaled_coverage) < 3
    ):
        gate_description = (
            "every treatment replicate > every control replicate"
            if args.stage1_gate_mode == "all-controls"
            else "mean treatment > mean control"
        )
        print(
            "WARNING: Stage 1 has fewer than three biological replicates in at "
            "least one group; Welch p-values and FDR are exploratory. Peak "
            f"selection uses the {args.stage1_gate_mode} gate ({gate_description}).",
            file=sys.stderr,
        )


def _run_stage1_reuse(
    args: argparse.Namespace,
    *,
    source_manifest_path: Path,
    source_manifest: dict[str, object],
    outdir: Path,
    condition_name: str,
    treatment_records: list[dict[str, object]],
    control_records: list[dict[str, object]],
    reporter: ProgressReporter,
    seed_offset: int,
    excluded_bams: Sequence[str],
) -> Path:
    """Recompute Stage 1 downstream of retained per-replicate BigWigs."""

    scaled_dir = outdir / "02_mean_scaled_tracks"
    peaks_dir = outdir / "03_peak_calls"
    fdr_dir = outdir / "04_peak_fdr"
    aggregate_dir = outdir / "05_cluster_aggregate"
    setup_dir = outdir / "00_setup"
    for directory in (scaled_dir, peaks_dir, fdr_dir, aggregate_dir, setup_dir):
        directory.mkdir(parents=True, exist_ok=True)

    score_track, positive_track = _TRACK_NAMES[args.scoring_method]
    target_scaled = [Path(str(record["scaled_score"])) for record in treatment_records]
    control_scaled = [Path(str(record["scaled_score"])) for record in control_records]
    target_scaled_coverage = [Path(str(record["scaled_coverage"])) for record in treatment_records]
    control_scaled_coverage = [Path(str(record["scaled_coverage"])) for record in control_records]
    for path in [*target_scaled, *control_scaled, *target_scaled_coverage, *control_scaled_coverage]:
        if not path.is_file():
            raise FileNotFoundError(
                f"A retained per-sample BigWig required for rerun is missing: {path}"
            )

    if len(target_scaled) == 1:
        mean_target = target_scaled[0]
        mean_target_coverage = target_scaled_coverage[0]
    else:
        reporter.stage(f"Re-averaging {condition_name} retained treatment score tracks")
        mean_target = average_bigwigs(
            target_scaled,
            scaled_dir / f"{args.sample_name}_condition_mean_target_{score_track}.bw",
        )
        reporter.stage(f"Re-averaging {condition_name} retained treatment coverage tracks")
        mean_target_coverage = average_bigwigs(
            target_scaled_coverage,
            scaled_dir / f"{args.sample_name}_condition_mean_target_coverage_x100.bw",
        )

    if len(control_scaled) == 1:
        mean_control = control_scaled[0]
        mean_control_coverage = control_scaled_coverage[0]
    else:
        reporter.stage(f"Re-averaging {condition_name} retained control score tracks")
        mean_control = average_bigwigs(
            control_scaled,
            scaled_dir / f"{args.sample_name}_condition_mean_control_{score_track}.bw",
        )
        reporter.stage(f"Re-averaging {condition_name} retained control coverage tracks")
        mean_control_coverage = average_bigwigs(
            control_scaled_coverage,
            scaled_dir / f"{args.sample_name}_condition_mean_control_coverage_x100.bw",
        )

    target_peaks = _call_peaks(
        args,
        reporter,
        scaled_path=mean_target,
        label="target",
        peaks_dir=peaks_dir,
        score_track=score_track,
    )
    reporter.stage(
        f"Retesting {condition_name} treatment candidates across retained coverage replicates"
    )
    _warn_small_replicate_groups(args, target_scaled_coverage, control_scaled_coverage)
    outputs = analyze_cutn_replicate_peaks(
        target_peaks,
        output_dir=fdr_dir,
        target_replicate_bigwigs=target_scaled_coverage,
        control_replicate_bigwigs=control_scaled_coverage,
        target_mean_bigwig=mean_target_coverage,
        peak_pvalue=args.stage1_p_value,
        peak_fdr=args.peak_fdr,
        cluster_seed_pvalue=args.cluster_seed_p_value,
        cluster_fdr=args.cluster_fdr,
        gate_mode=args.stage1_gate_mode,
        cluster_member_mode=args.cluster_member_mode,
        cluster_max_non_member_gap=args.cluster_max_non_member_gap,
        max_cluster_gap=args.max_cluster_gap,
        minimum_cluster_members=args.min_cluster_members,
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
            mean_scaled_score=mean_target,
            replicate_scaled_scores=target_scaled,
            anchor_bed=anchor_bed,
            output_dir=aggregate_dir,
            label=condition_name,
            scoring_method=args.scoring_method,
            positive_track=positive_track,
            window_half=args.cluster_aggregate_window_half,
            maximum_heatmap_rows=args.cluster_aggregate_max_heatmap_rows,
            bootstrap_replicates=args.cluster_aggregate_bootstrap,
            nrl_peak_resolution=args.cluster_aggregate_nrl_resolution,
            nrl_min_order=args.cluster_aggregate_nrl_min_order,
            nrl_max_order=args.cluster_aggregate_nrl_max_order,
            seed=args.seed + seed_offset + 300,
            reporter=reporter,
        )

    manifest = dict(source_manifest)
    manifest.update(
        {
            "nucleosuite_version": __version__,
            "condition_name": condition_name,
            "sample_name": args.sample_name,
            "peak_min_region_length": args.peak_min_region_length,
            "peak_max_neg_run": args.peak_max_neg_run,
            "peak_smooth_window": args.peak_smooth_window,
            "peak_smooth_order": args.peak_smooth_order,
            "stage1_p_value": args.stage1_p_value,
            "peak_fdr": args.peak_fdr,
            "cluster_seed_p_value": args.cluster_seed_p_value,
            "cluster_max_non_member_gap": args.cluster_max_non_member_gap,
            "max_cluster_gap": args.max_cluster_gap,
            "minimum_cluster_members": args.min_cluster_members,
            "cluster_fdr": args.cluster_fdr,
            "stage1_gate_mode": args.stage1_gate_mode,
            "cluster_member_mode": args.cluster_member_mode,
            "stage1_selection": (
                "mean_treatment_exceeds_mean_control"
                if args.stage1_gate_mode == "mean"
                else "all_treatments_exceed_all_controls"
            ),
            "treatment_replicates": treatment_records,
            "control_replicates": control_records,
            "condition_mean_treatment_score": str(Path(mean_target).resolve()),
            "condition_mean_control_score": str(Path(mean_control).resolve()),
            "condition_mean_treatment_coverage": str(Path(mean_target_coverage).resolve()),
            "condition_mean_control_coverage": str(Path(mean_control_coverage).resolve()),
            "condition_mean_treatment_cluster_aggregate_score": str(Path(mean_target).resolve()),
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
            "target_candidate_peaks": str(target_peaks),
            "rerun_source_manifest": str(source_manifest_path.resolve()),
            "rerun_excluded_bams": list(excluded_bams),
            **{name: str(path.resolve()) for name, path in outputs.items()},
        }
    )
    manifest_path = outdir / MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    summary = outdir / f"{args.sample_name}_cutn_suite_summary.tsv"
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"condition_name\t{condition_name}\n")
        handle.write(f"rerun_source_manifest\t{source_manifest_path}\n")
        handle.write(f"excluded_bams\t{';'.join(excluded_bams)}\n")
        handle.write(f"retained_treatment_replicates\t{len(treatment_records)}\n")
        handle.write(f"retained_control_replicates\t{len(control_records)}\n")
        handle.write(f"peak_min_region_length\t{args.peak_min_region_length}\n")
        handle.write(f"peak_max_neg_run\t{args.peak_max_neg_run}\n")
        handle.write(f"peak_smooth_window\t{args.peak_smooth_window}\n")
        handle.write(f"peak_smooth_order\t{args.peak_smooth_order}\n")
        handle.write(f"stage1_gate_mode\t{args.stage1_gate_mode}\n")
        handle.write(f"cluster_member_mode\t{args.cluster_member_mode}\n")
        handle.write(f"cluster_seed_p_value\t{args.cluster_seed_p_value}\n")
        handle.write(f"cluster_max_non_member_gap\t{args.cluster_max_non_member_gap}\n")
        handle.write(f"max_cluster_gap\t{args.max_cluster_gap}\n")
        handle.write(f"minimum_cluster_members\t{args.min_cluster_members}\n")
        handle.write(f"stage1_manifest\t{manifest_path}\n")
    return manifest_path


def _relative_reference(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _write_run_manifest(
    root: Path,
    *,
    condition_manifests: Sequence[Path],
    comparison_manifest: Path | None,
    args: argparse.Namespace,
    source_run: Path | None = None,
    excluded_bams: Sequence[str] = (),
    changed_parameters: dict[str, object] | None = None,
) -> Path:
    conditions = []
    for path in condition_manifests:
        manifest = _read_json(path)
        conditions.append(
            {
                "condition_name": manifest.get("condition_name"),
                "manifest": _relative_reference(root, path),
                "treatment_bams": [
                    bam
                    for record in manifest.get("treatment_replicates", [])
                    if isinstance(record, dict)
                    for bam in record.get("bams", [])
                ],
                "control_bams": [
                    bam
                    for record in manifest.get("control_replicates", [])
                    if isinstance(record, dict)
                    for bam in record.get("bams", [])
                ],
            }
        )
    payload = {
        "schema": _RUN_MANIFEST_SCHEMA,
        "schema_version": _RUN_MANIFEST_VERSION,
        "nucleosuite_version": __version__,
        "sample_name": args.sample_name,
        "run_directory": str(root.resolve()),
        "condition_manifests": [_relative_reference(root, path) for path in condition_manifests],
        "conditions": conditions,
        "comparison_manifest": (
            _relative_reference(root, comparison_manifest) if comparison_manifest else None
        ),
        "source_run": str(source_run.resolve()) if source_run else None,
        "excluded_bams": list(excluded_bams),
        "changed_parameters": changed_parameters or {},
        "parameters": {
            "bam_mode": args.bam_mode,
            "scoring_method": args.scoring_method,
            "coverage_frag_lower": args.coverage_frag_lower,
            "coverage_frag_upper": args.coverage_frag_upper,
            "frag_mode_padding": args.frag_mode_padding,
            "score_frag_lower": args.score_frag_lower,
            "score_frag_upper": args.score_frag_upper,
            "contigs": list(args.contigs),
            "blacklist_bed": args.blacklist_bed,
            "max_duplicates": args.max_duplicates,
            "dedup_scope": args.dedup_scope,
            "peak_min_region_length": args.peak_min_region_length,
            "peak_max_neg_run": args.peak_max_neg_run,
            "peak_smooth_window": args.peak_smooth_window,
            "peak_smooth_order": args.peak_smooth_order,
            "stage1_p_value": args.stage1_p_value,
            "peak_fdr": args.peak_fdr,
            "cluster_fdr": args.cluster_fdr,
            "cluster_seed_p_value": args.cluster_seed_p_value,
            "stage1_gate_mode": args.stage1_gate_mode,
            "cluster_member_mode": args.cluster_member_mode,
            "cluster_max_non_member_gap": args.cluster_max_non_member_gap,
            "max_cluster_gap": args.max_cluster_gap,
            "min_cluster_members": args.min_cluster_members,
            "differential_fdr": args.differential_fdr,
            "cluster_aggregate_window_half": args.cluster_aggregate_window_half,
            "cluster_aggregate_max_heatmap_rows": args.cluster_aggregate_max_heatmap_rows,
            "cluster_aggregate_bootstrap": args.cluster_aggregate_bootstrap,
            "cluster_aggregate_nrl_resolution": args.cluster_aggregate_nrl_resolution,
            "cluster_aggregate_nrl_min_order": args.cluster_aggregate_nrl_min_order,
            "cluster_aggregate_nrl_max_order": args.cluster_aggregate_nrl_max_order,
            "skip_cluster_aggregate": args.skip_cluster_aggregate,
            "seed": args.seed,
        },
    }
    path = root / _RUN_MANIFEST
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _run_rerun(args: argparse.Namespace) -> int:
    _validate_rerun_immutable_options(args)
    source_root, source_conditions, source_run_manifest = _discover_run(args.rerun_from)
    if len(source_conditions) > 2:
        raise ValueError("cutn-suite reruns support at most two conditions")
    changed = _inherit_rerun_parameters(args, source_conditions[0][1])
    if source_run_manifest:
        source_run_meta = _read_json(source_run_manifest)
        source_parameters = source_run_meta.get("parameters")
        source_parameters = source_parameters if isinstance(source_parameters, dict) else {}
        if not _explicit(args, "--differential-fdr"):
            args.differential_fdr = float(source_parameters.get("differential_fdr", 0.05))
        else:
            changed["differential_fdr"] = args.differential_fdr
        if source_parameters.get("seed") is not None:
            args.seed = int(source_parameters["seed"])
    else:
        if _explicit(args, "--differential-fdr"):
            changed["differential_fdr"] = args.differential_fdr
        else:
            old_comparison = source_root / "03_condition_comparison" / "cutn_comparison_manifest.json"
            if old_comparison.is_file():
                old_comparison_meta = _read_json(old_comparison)
                if old_comparison_meta.get("fdr") is not None:
                    args.differential_fdr = float(old_comparison_meta["fdr"])
    excluded_records, matched_bams = _match_excluded_records(source_conditions, args.exclude_sample)
    rerun_dir = _next_rerun_directory(source_root, matched_bams)
    rerun_dir.mkdir(parents=True, exist_ok=False)

    reporter = ProgressReporter("cutn-suite")
    reporter.stage(f"Reusing per-sample BigWigs from {source_root}")
    if matched_bams:
        reporter.stage("Excluding: " + ", ".join(matched_bams))

    rerun_manifests: list[Path] = []
    for condition_index, (source_manifest_path, source_manifest) in enumerate(source_conditions, 1):
        treatment_records, control_records = _filtered_records(
            source_manifest, condition_index, excluded_records
        )
        if condition_index == 1:
            condition_name = (
                args.condition1_name
                if _explicit(args, "--condition1-name")
                else str(source_manifest.get("condition_name") or "condition1")
            )
        else:
            condition_name = (
                args.condition2_name
                if _explicit(args, "--condition2-name")
                else str(source_manifest.get("condition_name") or "condition2")
            )
        condition_outdir = (
            rerun_dir
            if len(source_conditions) == 1
            else rerun_dir / f"0{condition_index}_condition{condition_index}_stage1"
        )
        rerun_manifests.append(
            _run_stage1_reuse(
                args,
                source_manifest_path=source_manifest_path,
                source_manifest=source_manifest,
                outdir=condition_outdir,
                condition_name=condition_name,
                treatment_records=treatment_records,
                control_records=control_records,
                reporter=reporter,
                seed_offset=condition_index * 1000,
                excluded_bams=matched_bams,
            )
        )

    comparison_manifest = None
    if len(rerun_manifests) == 2:
        reporter.stage("Recomputing Stage 2 comparison from rerun Stage 1 outputs")
        comparison_manifest = compare_stage1(
            rerun_manifests[0],
            rerun_manifests[1],
            outdir=rerun_dir / "03_condition_comparison",
            fdr=args.differential_fdr,
        )
    run_manifest = _write_run_manifest(
        rerun_dir,
        condition_manifests=rerun_manifests,
        comparison_manifest=comparison_manifest,
        args=args,
        source_run=source_root,
        excluded_bams=matched_bams,
        changed_parameters=changed,
    )
    print(f"cutn_rerun_directory\t{rerun_dir}")
    print(f"cutn_run_manifest\t{run_manifest}")
    if comparison_manifest:
        print(f"cutn_comparison_manifest\t{comparison_manifest}")
    return 0

def run(args: argparse.Namespace) -> int:
    _validate(args)
    if args.inspect_run:
        return inspect_run(args.inspect_run)
    if args.rerun_from:
        return _run_rerun(args)

    outdir = Path(args.outdir).resolve()
    condition1_name = args.condition1_name
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
                else f"mode_plus_minus_{args.frag_mode_padding}"
            )
        )
        print(
            f"coverage_fragment_range\t{args.coverage_frag_lower}-"
            f"{args.coverage_frag_upper}"
        )
        print(
            f"stages\tstage1-treatment-candidates-{args.stage1_gate_mode}-gate"
            + (",stage2-log2-moderated-four-group-interaction-bh" if has_second else "")
        )
        return 0

    reporter = ProgressReporter("cutn-suite")
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
        run_manifest = _write_run_manifest(
            outdir,
            condition_manifests=[first_manifest],
            comparison_manifest=None,
            args=args,
        )
        print(f"cutn_stage1_manifest\t{first_manifest}")
        print(f"cutn_run_manifest\t{run_manifest}")
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
    summary = outdir / f"{args.sample_name}_cutn_suite_summary.tsv"
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
    run_manifest = _write_run_manifest(
        outdir,
        condition_manifests=[first_manifest, second_manifest],
        comparison_manifest=comparison_manifest,
        args=args,
    )
    print(f"cutn_suite_summary\t{summary}")
    print(f"cutn_comparison_manifest\t{comparison_manifest}")
    print(f"cutn_run_manifest\t{run_manifest}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    raw = list(argv) if argv is not None else sys.argv[1:]
    args = build_parser().parse_args(raw)
    args._explicit_options = {
        token.split("=", 1)[0]
        for token in raw
        if token.startswith("-")
    }
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
