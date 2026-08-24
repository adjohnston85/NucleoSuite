#!/usr/bin/env python3
"""Matched ChIP-seq/CUT&RUN/CUT&Tag nucleosome-scoring workflow."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from nucleosuite.chip_peaks import analyze_chip_peaks
from nucleosuite.mean_scale import scale_bigwig_by_reference
from nucleosuite.mode_estimation import (
    ModeEstimate,
    estimate_bam_fragment_mode,
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
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "sample\tmode\tbootstrap_ci_low\tbootstrap_ci_high\t"
            "sampled_fragments\tmode_search_fragments\tconverged\tcheckpoints\tmode_source\n"
        )
        for label, estimate, analysis_mode in (
            ("target", target, target_mode),
            ("control", control, control_mode),
            ("pooled", pooled, target_mode if target_mode == control_mode else -1),
        ):
            if estimate is None:
                if label == "target":
                    handle.write(
                        f"target\t{target_mode}\t\t\t\t\t\t\t{mode_source}\n"
                    )
                elif label == "control":
                    handle.write(
                        f"control\t{control_mode}\t\t\t\t\t\t\t{mode_source}\n"
                    )
                continue
            handle.write(
                f"{label}\t{estimate.mode}\t{estimate.ci_low:.6g}\t"
                f"{estimate.ci_high:.6g}\t{estimate.sampled_fragments}\t"
                f"{estimate.mode_search_fragments}\t{str(estimate.converged).lower()}\t"
                f"{estimate.checkpoints}\t{mode_source}\n"
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite chip-suite",
        description=(
            "Generate matched, positive-score-mean-scaled TNS/BNS/PNS tracks from "
            "target and control ChIP-seq, CUT&RUN or CUT&Tag BAMs; call peaks "
            "independently; and estimate target peak and cluster FDR from the control."
        ),
    )
    parser.add_argument(
        "--target-bam", nargs="+", required=True,
        help="Coordinate-sorted, indexed target BAM file(s).",
    )
    parser.add_argument(
        "--control-bam", nargs="+", required=True,
        help="Coordinate-sorted, indexed matched-control BAM file(s).",
    )
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--sample-name", default="chip", help="Output sample prefix (default: chip).")
    parser.add_argument(
        "--scoring-method", choices=("tns", "bns", "pns"), default="tns",
        help="Nucleosome scoring kernel (default: tns).",
    )
    parser.add_argument(
        "--mode", type=_mode_value, default="auto",
        help=(
            "Scoring mode: auto estimates target and control modes by random-block "
            "bootstrap sampling; an integer bypasses estimation and is used exactly "
            "for both tracks (default: auto)."
        ),
    )
    parser.add_argument(
        "--mode-strategy", choices=("pooled", "separate", "target", "control"),
        default="pooled",
        help="How independently estimated modes define analysis modes (default: pooled).",
    )
    parser.add_argument("--frag-lower", type=int, default=120, help="Minimum fragment length (default: 120).")
    parser.add_argument("--frag-upper", type=int, default=500, help="Maximum fragment length (default: 500).")
    parser.add_argument("--mode-search-lower", type=int, default=120, help="Automatic mode-search lower bound (default: 120).")
    parser.add_argument("--mode-search-upper", type=int, default=250, help="Automatic mode-search upper bound (default: 250).")
    parser.add_argument("--mode-min-fragments", type=int, default=100_000, help="Minimum sampled fragments before convergence checks (default: 100000).")
    parser.add_argument("--mode-batch-fragments", type=int, default=25_000, help="Fragments between convergence checks (default: 25000).")
    parser.add_argument("--mode-max-fragments", type=int, default=1_000_000, help="Maximum sampled fragments per BAM collection (default: 1000000).")
    parser.add_argument("--mode-bootstrap", type=int, default=200, help="Bootstrap mode replicates per checkpoint (default: 200).")
    parser.add_argument("--mode-stable-checkpoints", type=int, default=3, help="Required stable checkpoints (default: 3).")
    parser.add_argument("--mode-max-change", type=int, default=1, help="Maximum mode range across stable checkpoints (default: 1 bp).")
    parser.add_argument("--mode-max-ci-width", type=float, default=4.0, help="Maximum bootstrap 95%% interval width (default: 4 bp).")
    parser.add_argument("--mode-block-bp", type=int, default=1_000_000, help="Random genomic-block size (default: 1000000 bp).")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed (default: 12345).")
    parser.add_argument("--blacklist-bed", help="BED blacklist applied to both samples.")
    parser.add_argument("-c", "--contigs", nargs="+", default=["autosomes"], help="Contigs analysed in both samples (default: autosomes).")
    parser.add_argument("--max-duplicates", type=int, default=1, help="Maximum identical fragments retained (default: 1).")
    parser.add_argument(
        "--dedup-scope", choices=("all_bams", "per_bam"), default="all_bams",
        help="Apply coordinate duplicate limits across all BAMs or within each BAM (default: all_bams).",
    )
    parser.add_argument("--cores", type=int, default=1, help="Parallel contig workers (default: 1).")
    parser.add_argument("--peak-fdr", type=float, default=0.05, help="Peak empirical-FDR cutoff (default: 0.05).")
    parser.add_argument("--cluster-fdr", type=float, default=0.05, help="Cluster empirical-FDR cutoff (default: 0.05).")
    parser.add_argument("--peak-match-distance", type=int, help="Target-control summit matching distance; default: half the analysis mode.")
    parser.add_argument("--cluster-break", type=int, default=5, help="Consecutive nonsignificant peaks ending a cluster (default: 5).")
    parser.add_argument("--max-cluster-gap", type=int, default=1000, help="Maximum gap between significant peak summits (default: 1000 bp).")
    parser.add_argument("--min-significant-peaks", type=int, default=2, help="Minimum significant peaks per cluster (default: 2).")
    parser.add_argument("--dry-run", action="store_true", help="Validate and report planned stages without processing BAMs.")
    return parser


def _validate(args: argparse.Namespace) -> None:
    if args.frag_lower < 1 or args.frag_upper < args.frag_lower:
        raise ValueError("Require 1 <= --frag-lower <= --frag-upper")
    if args.mode == "auto":
        if (
            args.mode_search_lower < args.frag_lower
            or args.mode_search_upper > args.frag_upper
        ):
            raise ValueError("Automatic mode-search bounds must lie within the fragment range")
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
    if args.peak_match_distance is not None and args.peak_match_distance < 0:
        raise ValueError("--peak-match-distance must be non-negative")
    if min(args.cluster_break, args.max_cluster_gap, args.min_significant_peaks) < 1:
        raise ValueError("Cluster break, gap and significant-peak count must be positive")
    for name in ("peak_fdr", "cluster_fdr"):
        value = float(getattr(args, name))
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be between 0 and 1")
    for path in [*args.target_bam, *args.control_bam]:
        if not Path(path).is_file():
            raise FileNotFoundError(path)


def run(args: argparse.Namespace) -> int:
    _validate(args)
    outdir = Path(args.outdir).resolve()
    tracks_dir = outdir / "01_score_tracks"
    scaled_dir = outdir / "02_mean_scaled_tracks"
    peaks_dir = outdir / "03_peak_calls"
    fdr_dir = outdir / "04_peak_fdr"
    setup_dir = outdir / "00_setup"
    if args.dry_run:
        print(f"score_type\t{args.scoring_method}")
        print(f"mode\t{args.mode}")
        print(f"fragment_range\t{args.frag_lower}-{args.frag_upper}")
        print("stages\tmode-estimation,score-tracks,positive-track-mean-scaling,peak-calling,peak-fdr,cluster-fdr")
        return 0

    for directory in (tracks_dir, scaled_dir, peaks_dir, fdr_dir, setup_dir):
        directory.mkdir(parents=True, exist_ok=True)

    reporter = ProgressReporter("chip-suite")
    target_estimate = control_estimate = pooled_estimate = None
    if args.mode == "auto":
        reporter.stage("Estimating target fragment mode")
        common = dict(
            frag_lower=args.frag_lower,
            frag_upper=args.frag_upper,
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
            seed=args.seed,
            blacklist_bed=args.blacklist_bed,
            max_duplicates=args.max_duplicates,
            dedup_scope=args.dedup_scope,
            contig_tokens=args.contigs,
        )
        target_estimate = estimate_bam_fragment_mode(args.target_bam, **common)
        reporter.stage("Estimating control fragment mode")
        control_estimate = estimate_bam_fragment_mode(
            args.control_bam, **{**common, "seed": args.seed + 1}
        )
        pooled_estimate = pooled_mode_estimate(
            target_estimate,
            control_estimate,
            bootstrap_replicates=args.mode_bootstrap,
            seed=args.seed + 2,
        )
        if args.mode_strategy == "pooled":
            target_mode = control_mode = pooled_estimate.mode
        elif args.mode_strategy == "separate":
            target_mode, control_mode = target_estimate.mode, control_estimate.mode
        elif args.mode_strategy == "target":
            target_mode = control_mode = target_estimate.mode
        else:
            target_mode = control_mode = control_estimate.mode
        if abs(target_estimate.mode - control_estimate.mode) > 5:
            print(
                "WARNING: target and control fragment modes differ by more than 5 bp; "
                "inspect the mode-estimation report before interpretation",
                file=sys.stderr,
            )
        mode_source = f"automatic_{args.mode_strategy}"
    else:
        target_mode = control_mode = int(args.mode)
        mode_source = "explicit"

    _write_mode_report(
        setup_dir / f"{args.sample_name}_fragment_mode_estimation.tsv",
        mode_source=mode_source,
        target=target_estimate,
        control=control_estimate,
        pooled=pooled_estimate,
        target_mode=target_mode,
        control_mode=control_mode,
    )

    score_track, positive_track = _TRACK_NAMES[args.scoring_method]
    peak_paths: dict[str, Path] = {}
    for label, bam_paths, mode in (
        ("target", args.target_bam, target_mode),
        ("control", args.control_bam, control_mode),
    ):
        reporter.stage(f"Generating {label} {args.scoring_method.upper()} and {positive_track}")
        base = tracks_dir / f"{args.sample_name}_{label}"
        command = [
            "pns", "--bam", *bam_paths,
            "--scoring-method", args.scoring_method,
            "--mode-length", str(mode),
            "--frag-lower", str(args.frag_lower),
            "--frag-upper", str(args.frag_upper),
            "--score-tracks", score_track, positive_track,
            "--other-tracks", "none",
            "--pns-format", "bigwig",
            "--other-format", "none",
            "--interval-format", "bed",
            "--out-prefix", str(base),
            "--contigs", *args.contigs,
            "--max-duplicates", str(args.max_duplicates),
            "--dedup-scope", args.dedup_scope,
            "--cores", str(args.cores),
            "--seed", str(args.seed),
        ]
        if args.blacklist_bed:
            command.extend(["--blacklist-bed", args.blacklist_bed])
        _run_nucleosuite(command)

        direct_prefix = _resolved_prefix(
            base, args.scoring_method, mode, args.frag_lower, args.frag_upper
        )
        prefix = _locate_completed_prefix(
            base,
            direct_prefix,
            (f"_{score_track}.bw", f"_{positive_track}.bw"),
        )
        score_path = Path(f"{prefix}_{score_track}.bw")
        positive_path = Path(f"{prefix}_{positive_track}.bw")
        scaled_path = scaled_dir / (
            f"{args.sample_name}_{label}_{score_track}_divided_by_mean_{positive_track}.bw"
        )
        reporter.stage(f"Scaling {label} {score_track} by mean {positive_track}")
        _scaled, reference_mean, reference_count = scale_bigwig_by_reference(
            score_path, positive_path, scaled_path, scale=1.0, reporter=reporter
        )
        with (setup_dir / f"{args.sample_name}_{label}_score_scaling.tsv").open(
            "wt", encoding="utf-8"
        ) as handle:
            handle.write("field\tvalue\n")
            handle.write(f"score_track\t{score_path}\n")
            handle.write(f"positive_reference_track\t{positive_path}\n")
            handle.write(f"positive_reference_mean\t{reference_mean:.12g}\n")
            handle.write(f"finite_nonzero_reference_bases\t{reference_count}\n")
            handle.write(f"scaled_track\t{scaled_path}\n")

        peak_prefix = peaks_dir / f"{args.sample_name}_{label}_{score_track}_mean_scaled"
        peak_command = [
            "call-peaks",
            "--input-bigwig", str(scaled_path),
            "--peak-caller", "pns",
            "--call-type", "both",
            "--out-prefix", str(peak_prefix),
            "--regions", *args.contigs,
            "--cores", str(args.cores),
            "--interval-format", "bed",
        ]
        if args.blacklist_bed:
            peak_command.extend(["--blacklist-bed", args.blacklist_bed])
        reporter.stage(f"Calling {label} peaks on mean-scaled {score_track}")
        _run_nucleosuite(peak_command)
        completed_peak_prefix = _locate_completed_prefix(
            peak_prefix,
            peak_prefix,
            ("_nucleosome_regions.bed",),
        )
        peak_paths[label] = Path(
            f"{completed_peak_prefix}_nucleosome_regions.bed"
        )

    match_distance = (
        args.peak_match_distance
        if args.peak_match_distance is not None
        else max(0, int(round((target_mode + control_mode) / 4.0)))
    )
    reporter.stage("Estimating target peak and cluster FDR from the matched control")
    outputs = analyze_chip_peaks(
        peak_paths["target"],
        peak_paths["control"],
        output_dir=fdr_dir,
        match_distance=match_distance,
        peak_fdr=args.peak_fdr,
        cluster_fdr=args.cluster_fdr,
        cluster_break=args.cluster_break,
        max_cluster_gap=args.max_cluster_gap,
        minimum_significant_peaks=args.min_significant_peaks,
    )

    summary = outdir / f"{args.sample_name}_chip_suite_summary.tsv"
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"scoring_method\t{args.scoring_method}\n")
        handle.write(f"mode_source\t{mode_source}\n")
        handle.write(f"target_analysis_mode\t{target_mode}\n")
        handle.write(f"control_analysis_mode\t{control_mode}\n")
        handle.write(f"fragment_range\t{args.frag_lower}-{args.frag_upper}\n")
        handle.write(f"peak_match_distance\t{match_distance}\n")
        handle.write(f"peak_fdr\t{args.peak_fdr}\n")
        handle.write(f"cluster_fdr\t{args.cluster_fdr}\n")
        for name, path in outputs.items():
            handle.write(f"{name}\t{path}\n")
    print(f"chip_suite_summary\t{summary}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
