"""Cluster-centred aggregate outputs for the ChIP workflow."""

from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from nucleosuite.align import AlignmentConfig, run_alignment
from nucleosuite.progress import ProgressReporter


def common_symmetric_bigwig_limit(paths: Sequence[str | Path]) -> float | None:
    """Estimate one symmetric colour limit shared by compatible BigWigs."""

    try:
        import pyBigWig
    except ImportError:  # pragma: no cover
        return None
    absolute = 0.0
    for path in paths:
        handle = pyBigWig.open(str(path))
        if handle is None:
            continue
        try:
            for chrom, length in handle.chroms().items():
                for statistic in ("min", "max"):
                    value = handle.stats(chrom, 0, int(length), type=statistic)[0]
                    if value is not None and math.isfinite(float(value)):
                        absolute = max(absolute, abs(float(value)))
        finally:
            handle.close()
    return absolute if absolute > 0 else None


def write_cluster_anchor_bed(cluster_bed: str | Path, output_path: str | Path) -> Path:
    """Write one-base anchors at the strongest member summit of each cluster."""

    source = Path(cluster_bed)
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rt", encoding="utf-8") as handle, output.open(
        "wt", encoding="utf-8"
    ) as destination:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                raise ValueError(
                    f"Cluster BED line {line_number} in {source} lacks summit column 7"
                )
            try:
                summit = int(fields[6])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid cluster summit at line {line_number} in {source}"
                ) from exc
            cluster_id = fields[3] if len(fields) >= 4 and fields[3] else f"cluster_{line_number}"
            score = fields[4] if len(fields) >= 5 else "0"
            destination.write(
                f"{fields[0]}\t{summit}\t{summit + 1}\t{cluster_id}\t{score}\t.\t"
                f"{summit}\t{summit + 1}\n"
            )
    return output


def _read_profile(path: Path) -> tuple[np.ndarray, np.ndarray]:
    positions: list[float] = []
    values: list[float] = []
    with path.open("rt", encoding="utf-8") as handle:
        next(handle, None)
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            positions.append(float(fields[0]))
            values.append(float(fields[1]))
    return np.asarray(positions), np.asarray(values)


def _read_matrix(path: Path) -> tuple[np.ndarray, np.ndarray]:
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    rows: list[list[float]] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n").split("\t")
        positions = np.asarray([float(value) for value in header[1:]])
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) == len(header):
                rows.append([float(value) for value in fields[1:]])
    return positions, np.asarray(rows, dtype=float)


def _write_bootstrap_confidence(
    matrix_path: Path,
    output_dir: Path,
    *,
    prefix: str,
    replicates: int,
    seed: int,
    maximum_rows: int = 1000,
) -> dict[str, str]:
    if replicates < 1:
        return {}
    positions, matrix = _read_matrix(matrix_path)
    if matrix.size == 0:
        return {}
    rng = np.random.default_rng(seed)
    if matrix.shape[0] > maximum_rows:
        matrix = matrix[
            np.sort(rng.choice(matrix.shape[0], size=maximum_rows, replace=False))
        ]
    estimates = np.empty((replicates, matrix.shape[1]), dtype=np.float32)
    for index in range(replicates):
        selected = rng.integers(0, matrix.shape[0], size=matrix.shape[0])
        estimates[index] = np.nanmean(matrix[selected], axis=0)
    lower, upper = np.nanquantile(estimates, (0.025, 0.975), axis=0)
    mean = np.nanmean(matrix, axis=0)
    table_path = output_dir / f"{prefix}_bootstrap95_profile.tsv"
    with table_path.open("wt", encoding="utf-8") as handle:
        handle.write(
            "relative_position\tmean_scaled_pns\tbootstrap_95_ci_lower\t"
            "bootstrap_95_ci_upper\theatmap_rows_used\tbootstrap_replicates\n"
        )
        for position, observed, low, high in zip(positions, mean, lower, upper):
            handle.write(
                f"{position:g}\t{observed:.12g}\t{low:.12g}\t{high:.12g}\t"
                f"{matrix.shape[0]}\t{replicates}\n"
            )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_path = output_dir / f"{prefix}_bootstrap95_profile.png"
    figure, axis = plt.subplots(figsize=(11, 5))
    axis.fill_between(positions, lower, upper, color="0.75", linewidth=0)
    axis.plot(positions, mean, color="black", linewidth=1.4)
    axis.axvline(0, color="0.5", linewidth=0.8)
    axis.set_xlabel("Distance from strongest cluster peak (bp)")
    axis.set_ylabel("Mean PNS divided by mean posPNS")
    figure.tight_layout()
    figure.savefig(plot_path, dpi=300)
    plt.close(figure)
    return {
        "bootstrap_profile": str(table_path.resolve()),
        "bootstrap_profile_plot": str(plot_path.resolve()),
    }


def _write_replicate_overlay(
    profiles: Sequence[tuple[str, Path]], output_path: Path
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(11, 5))
    for label, path in profiles:
        positions, values = _read_profile(path)
        combined = label == "replicate_combined"
        axis.plot(
            positions,
            values,
            color="black" if combined else None,
            linewidth=2.0 if combined else 0.9,
            alpha=1.0 if combined else 0.75,
            label=label,
        )
    axis.axvline(0, color="0.5", linewidth=0.8)
    axis.set_xlabel("Distance from strongest cluster peak (bp)")
    axis.set_ylabel("Mean PNS divided by mean posPNS")
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    plt.close(figure)
    return output_path.resolve()


def run_cluster_aggregate(
    *,
    mean_scaled_pns: str | Path,
    replicate_scaled_pns: Sequence[str | Path],
    anchor_bed: str | Path,
    output_dir: str | Path,
    label: str,
    window_half: int = 1000,
    maximum_heatmap_rows: int = 5000,
    bootstrap_replicates: int = 200,
    nrl_peak_resolution: float = 140.0,
    nrl_min_order: int = 0,
    nrl_max_order: int = 3,
    seed: int = 12345,
    reporter: ProgressReporter | None = None,
    vlim: float | None = None,
) -> dict[str, object]:
    """Align scaled PNS to cluster anchors and write profiles, heatmaps and NRLs."""

    anchors = Path(anchor_bed).resolve()
    directory = Path(output_dir).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if not anchors.is_file():
        raise FileNotFoundError(anchors)
    if not anchors.read_text(encoding="utf-8").strip():
        status = directory / f"{label}_cluster_aggregate_status.tsv"
        status.write_text("status\tno_selected_clusters\n", encoding="utf-8")
        return {"status": "no_selected_clusters", "status_table": str(status)}

    if reporter is not None:
        reporter.stage(f"Aggregating {label} mean-scaled PNS at strongest cluster peaks")
    common = dict(
        region_bed=anchors,
        output_dir=directory,
        window_half=window_half,
        point_col=7,
        strand_col=6,
        missing_strand="forward",
        zero_thresh=0,
        max_score=float("inf"),
        nan_to_zero=True,
        max_heatmap_rows=maximum_heatmap_rows,
        subsample_mode="random",
        seed=seed,
        sort_mode="unsorted",
        axis_label="Distance from strongest cluster peak",
        colorbar_label="PNS divided by mean posPNS",
        mean_ylabel="Mean PNS divided by mean posPNS",
        vmin=None if vlim is None else -abs(vlim),
        vmax=None if vlim is None else abs(vlim),
    )
    combined_config = AlignmentConfig(
        bigwig=Path(mean_scaled_pns).resolve(),
        output_prefix=f"{label}_replicate_combined_scaled_pns_cluster_aligned",
        write_detail_tables=True,
        nrl=True,
        nrl_peak_resolution=nrl_peak_resolution,
        nrl_regression_min=0.0,
        nrl_regression_max=None,
        nrl_regression_min_order=nrl_min_order,
        nrl_regression_max_order=nrl_max_order,
        nrl_exclusion=False,
        **common,
    )
    combined_outputs = run_alignment(combined_config, progress=reporter)

    profile_paths: list[tuple[str, Path]] = [
        ("replicate_combined", combined_outputs["aggregate"])
    ]
    replicate_outputs: list[dict[str, str]] = []
    for index, path in enumerate(replicate_scaled_pns, 1):
        if reporter is not None:
            reporter.stage(f"Aggregating {label} PNS replicate {index}")
        config = AlignmentConfig(
            bigwig=Path(path).resolve(),
            output_prefix=f"{label}_replicate_{index}_scaled_pns_cluster_aligned",
            write_detail_tables=False,
            nrl=False,
            **common,
        )
        outputs = run_alignment(config, progress=None)
        profile_paths.append((f"replicate_{index}", outputs["aggregate"]))
        replicate_outputs.append({name: str(value) for name, value in outputs.items()})

    overlay = _write_replicate_overlay(
        profile_paths, directory / f"{label}_replicate_and_combined_profiles.png"
    )
    bootstrap = _write_bootstrap_confidence(
        combined_outputs["heatmap_matrix"],
        directory,
        prefix=f"{label}_replicate_combined_scaled_pns_cluster_aligned",
        replicates=bootstrap_replicates,
        seed=seed + 1,
    )
    return {
        "status": "complete",
        "anchor_bed": str(anchors),
        "combined": {name: str(value) for name, value in combined_outputs.items()},
        "replicates": replicate_outputs,
        "replicate_overlay_plot": str(overlay),
        **bootstrap,
    }
