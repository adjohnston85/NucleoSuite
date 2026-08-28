"""Shared paired observed/randomized execution for cfDNA and MNase suites."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence

from nucleosuite.peak_fdr import PeakFdrResult, annotate_peak_fdr


def extract_paired_options(
    argv: Sequence[str],
) -> tuple[bool, float | None, list[str]]:
    """Consume wrapper-only paired-run and optional FDR-filter options."""
    paired = False
    fdr: float | None = None
    output: list[str] = []
    index = 0
    while index < len(argv):
        token = argv[index]
        if token == "--with-randomized-control":
            paired = True
            index += 1
            continue
        if token == "--fdr":
            if index + 1 >= len(argv):
                raise ValueError("--fdr requires a value")
            fdr = float(argv[index + 1])
            index += 2
            continue
        if token.startswith("--fdr="):
            fdr = float(token.split("=", 1)[1])
            index += 1
            continue
        output.append(token)
        index += 1
    if fdr is not None and (not math.isfinite(fdr) or not 0 <= fdr <= 1):
        raise ValueError("--fdr must be between 0 and 1")
    if fdr is not None and not paired:
        raise ValueError("Suite-level --fdr requires --with-randomized-control")
    return paired, fdr, output


def _single(paths: list[Path], label: str) -> Path:
    if len(paths) != 1:
        raise RuntimeError(
            f"Expected one {label} combined peak BED, found {len(paths)}: "
            + ", ".join(map(str, paths))
        )
    return paths[0]


def _paired_peak_paths(scaled_dir: Path, call_type: str) -> tuple[Path, Path]:
    matches = sorted(scaled_dir.glob(f"*_{call_type}_mean_scaled.bed"))
    observed = [path for path in matches if "_randomized_control_" not in path.name]
    randomized = [path for path in matches if "_randomized_control_" in path.name]
    return (
        _single(observed, f"observed {call_type}"),
        _single(randomized, f"randomized {call_type}"),
    )


def annotate_suite_combined_peaks(
    outdir: str | Path,
    *,
    suite_name: str,
    fdr_threshold: float | None,
) -> dict[str, PeakFdrResult]:
    """Annotate observed combined suite peaks using matched randomized outputs."""
    root = Path(outdir).resolve()
    combined = root / "combined" if (root / "combined" / "01_combined_tracks").is_dir() else root
    scaled = combined / "01_combined_tracks" / "scaled"
    if not scaled.is_dir():
        raise RuntimeError(f"Combined scaled peak directory was not created: {scaled}")
    output_dir = combined / "13_peak_analysis" / "pns" / "empirical_fdr"
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs: dict[str, PeakFdrResult] = {}
    for label, suffix in (
        ("nucleosome", "nucleosome_regions"),
        ("breakpoint", "breakpoint_peaks"),
    ):
        observed, randomized = _paired_peak_paths(scaled, suffix)
        prefix = output_dir / f"{observed.stem}_{suite_name}"
        outputs[label] = annotate_peak_fdr(
            observed,
            [randomized],
            score_column=5,
            fdr_threshold=fdr_threshold,
            output_prefix=prefix,
        )

    manifest = output_dir / f"{suite_name}_combined_peak_fdr_outputs.tsv"
    with manifest.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("peak_type\tannotated_bed\tsignificant_bed\tsummary\n")
        for label, result in outputs.items():
            handle.write(
                f"{label}\t{result.annotated_path}\t"
                f"{result.significant_path or ''}\t{result.summary_path}\n"
            )
    return outputs
