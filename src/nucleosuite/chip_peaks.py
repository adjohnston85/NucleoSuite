"""Target-control peak competition and empirical cluster FDR for chip-suite."""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

import numpy as np

from nucleosuite.peak_fdr import PeakRow, _read_peak_rows, empirical_peak_qvalues


@dataclass(frozen=True)
class CompetitivePeak:
    row: PeakRow
    chrom: str
    start: int
    end: int
    summit: int
    source: str
    winner: bool = False
    matched_score: float = 0.0
    qvalue: float = 1.0


@dataclass(frozen=True)
class PeakCluster:
    chrom: str
    start: int
    end: int
    significant_peaks: tuple[CompetitivePeak, ...]
    score: float
    summit: int
    max_peak_score: float
    qvalue: float = 1.0


def _to_competitive(
    rows: Sequence[PeakRow], source: str, summit_column: int
) -> list[CompetitivePeak]:
    if summit_column < 1:
        raise ValueError("summit column must be at least 1")
    index = summit_column - 1
    output: list[CompetitivePeak] = []
    for row in rows:
        if index >= len(row.fields):
            raise ValueError(
                f"Peak record at input line {row.source_line} lacks summit column {summit_column}"
            )
        try:
            start, end, summit = int(row.fields[1]), int(row.fields[2]), int(row.fields[index])
        except ValueError as exc:
            raise ValueError(
                f"Invalid peak coordinates at input line {row.source_line}"
            ) from exc
        if start < 0 or end <= start or summit < start or summit >= end:
            raise ValueError(
                f"Invalid BED interval or summit at input line {row.source_line}"
            )
        output.append(
            CompetitivePeak(row, row.fields[0], start, end, summit, source)
        )
    output.sort(key=lambda peak: (peak.chrom, peak.summit, peak.start, peak.end))
    return output


def compete_peaks(
    target_rows: Sequence[PeakRow],
    control_rows: Sequence[PeakRow],
    *,
    match_distance: int,
    summit_column: int = 7,
) -> tuple[list[CompetitivePeak], list[CompetitivePeak]]:
    """Match nearby peaks one-to-one and label the higher-scoring member."""

    if match_distance < 0:
        raise ValueError("match_distance must be non-negative")
    target = _to_competitive(target_rows, "target", summit_column)
    control = _to_competitive(control_rows, "control", summit_column)
    candidates: list[tuple[int, float, int, int]] = []
    by_chrom_control: dict[str, list[tuple[int, CompetitivePeak]]] = {}
    for index, peak in enumerate(control):
        by_chrom_control.setdefault(peak.chrom, []).append((index, peak))
    control_summits = {
        chrom: [peak.summit for _index, peak in entries]
        for chrom, entries in by_chrom_control.items()
    }
    for target_index, target_peak in enumerate(target):
        chrom_controls = by_chrom_control.get(target_peak.chrom, [])
        summits = control_summits.get(target_peak.chrom, [])
        left = bisect_left(summits, target_peak.summit - match_distance)
        right = bisect_right(summits, target_peak.summit + match_distance)
        for control_index, control_peak in chrom_controls[left:right]:
            distance = abs(target_peak.summit - control_peak.summit)
            if distance <= match_distance:
                candidates.append(
                    (distance, -max(target_peak.row.score, control_peak.row.score), target_index, control_index)
                )
    candidates.sort()
    target_used: set[int] = set()
    control_used: set[int] = set()
    target_matches: dict[int, int] = {}
    control_matches: dict[int, int] = {}
    for _distance, _score, target_index, control_index in candidates:
        if target_index in target_used or control_index in control_used:
            continue
        target_used.add(target_index)
        control_used.add(control_index)
        target_matches[target_index] = control_index
        control_matches[control_index] = target_index

    updated_target: list[CompetitivePeak] = []
    for index, peak in enumerate(target):
        match = target_matches.get(index)
        if match is None:
            updated_target.append(replace(peak, winner=True, matched_score=0.0))
            continue
        other = control[match]
        updated_target.append(
            replace(
                peak,
                winner=peak.row.score > other.row.score,
                matched_score=other.row.score,
            )
        )
    updated_control: list[CompetitivePeak] = []
    for index, peak in enumerate(control):
        match = control_matches.get(index)
        if match is None:
            updated_control.append(replace(peak, winner=True, matched_score=0.0))
            continue
        other = target[match]
        # Ties conservatively belong to the control.
        updated_control.append(
            replace(
                peak,
                winner=peak.row.score >= other.row.score,
                matched_score=other.row.score,
            )
        )
    return updated_target, updated_control


def assign_competition_qvalues(
    target: Sequence[CompetitivePeak], control: Sequence[CompetitivePeak]
) -> tuple[list[CompetitivePeak], float | None]:
    """Assign target-decoy q-values to target-winning peaks."""

    target_indices = [index for index, peak in enumerate(target) if peak.winner]
    control_scores = [peak.row.score for peak in control if peak.winner]
    if not target_indices:
        return list(target), None
    qvalues, _sample_counts, _control_counts = empirical_peak_qvalues(
        [target[index].row.score for index in target_indices],
        [control_scores],
    )
    output = list(target)
    for index, qvalue in zip(target_indices, qvalues):
        output[index] = replace(output[index], qvalue=float(qvalue))
    return output, min((output[index].row.score for index in target_indices), default=None)


def cluster_peaks(
    peaks: Sequence[CompetitivePeak],
    *,
    significance_score: float,
    require_qvalue: float | None,
    cluster_break: int = 5,
    max_cluster_gap: int = 1000,
    minimum_significant_peaks: int = 2,
) -> list[PeakCluster]:
    """Cluster significant peaks while allowing a bounded number of bridges."""

    if cluster_break < 1 or max_cluster_gap < 1 or minimum_significant_peaks < 1:
        raise ValueError("Cluster break, maximum gap and minimum peak count must be positive")
    ordered = sorted(peaks, key=lambda peak: (peak.chrom, peak.summit, peak.start, peak.end))
    clusters: list[PeakCluster] = []
    current_significant: list[CompetitivePeak] = []
    nonsignificant_run = 0
    current_chrom: str | None = None

    def is_significant(peak: CompetitivePeak) -> bool:
        if not peak.winner or peak.row.score < significance_score:
            return False
        return require_qvalue is None or peak.qvalue <= require_qvalue

    def finish() -> None:
        nonlocal current_significant, nonsignificant_run
        if len(current_significant) >= minimum_significant_peaks:
            best = max(current_significant, key=lambda peak: (peak.row.score, -peak.summit))
            score = float(
                sum(max(peak.row.score - significance_score, 0.0) for peak in current_significant)
            )
            clusters.append(
                PeakCluster(
                    chrom=current_significant[0].chrom,
                    start=min(peak.start for peak in current_significant),
                    end=max(peak.end for peak in current_significant),
                    significant_peaks=tuple(current_significant),
                    score=score,
                    summit=best.summit,
                    max_peak_score=best.row.score,
                )
            )
        current_significant = []
        nonsignificant_run = 0

    for peak in ordered:
        if current_chrom is not None and peak.chrom != current_chrom:
            finish()
        current_chrom = peak.chrom
        if is_significant(peak):
            if (
                current_significant
                and peak.summit - current_significant[-1].summit > max_cluster_gap
            ):
                finish()
            current_significant.append(peak)
            nonsignificant_run = 0
        elif current_significant:
            nonsignificant_run += 1
            if nonsignificant_run >= cluster_break:
                finish()
    finish()
    return clusters


def assign_cluster_qvalues(
    target_clusters: Sequence[PeakCluster],
    control_clusters: Sequence[PeakCluster],
) -> list[PeakCluster]:
    if not target_clusters:
        return []
    qvalues, _sample_counts, _control_counts = empirical_peak_qvalues(
        [cluster.score for cluster in target_clusters],
        [[cluster.score for cluster in control_clusters]],
    )
    return [
        replace(cluster, qvalue=float(qvalue))
        for cluster, qvalue in zip(target_clusters, qvalues)
    ]


def analyze_chip_peaks(
    target_path: str | Path,
    control_path: str | Path,
    *,
    output_dir: str | Path,
    score_column: int = 5,
    summit_column: int = 7,
    match_distance: int = 84,
    peak_fdr: float = 0.05,
    cluster_fdr: float = 0.05,
    cluster_break: int = 5,
    max_cluster_gap: int = 1000,
    minimum_significant_peaks: int = 2,
) -> dict[str, Path]:
    """Run competition, peak FDR and cluster FDR and write BED/TSV outputs."""

    for name, value in (("peak_fdr", peak_fdr), ("cluster_fdr", cluster_fdr)):
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{name} must be between 0 and 1")
    target_rows = _read_peak_rows(target_path, score_column, allow_empty=True)
    control_rows = _read_peak_rows(control_path, score_column, allow_empty=True)
    target, control = compete_peaks(
        target_rows, control_rows,
        match_distance=match_distance,
        summit_column=summit_column,
    )
    target, _ = assign_competition_qvalues(target, control)
    passing = [peak for peak in target if peak.winner and peak.qvalue <= peak_fdr]
    threshold = min((peak.row.score for peak in passing), default=math.inf)

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    annotated = directory / "target_peaks_empirical_fdr.bed"
    significant = directory / f"target_peaks_fdr{peak_fdr:g}_significant.bed"
    with annotated.open("wt", encoding="utf-8") as all_handle, significant.open(
        "wt", encoding="utf-8"
    ) as sig_handle:
        for peak in sorted(target, key=lambda value: (value.chrom, value.start, value.end)):
            text = "\t".join((*peak.row.fields, f"{peak.qvalue:.12g}")) + "\n"
            all_handle.write(text)
            if peak.winner and peak.qvalue <= peak_fdr:
                sig_handle.write(text)

    target_clusters: list[PeakCluster] = []
    control_clusters: list[PeakCluster] = []
    if math.isfinite(threshold):
        target_clusters = cluster_peaks(
            target,
            significance_score=threshold,
            require_qvalue=peak_fdr,
            cluster_break=cluster_break,
            max_cluster_gap=max_cluster_gap,
            minimum_significant_peaks=minimum_significant_peaks,
        )
        control_clusters = cluster_peaks(
            control,
            significance_score=threshold,
            require_qvalue=None,
            cluster_break=cluster_break,
            max_cluster_gap=max_cluster_gap,
            minimum_significant_peaks=minimum_significant_peaks,
        )
        target_clusters = assign_cluster_qvalues(target_clusters, control_clusters)

    cluster_table = directory / "target_clusters_empirical_fdr.tsv"
    significant_clusters = directory / f"target_clusters_fdr{cluster_fdr:g}_significant.bed"
    with cluster_table.open("wt", encoding="utf-8") as table, significant_clusters.open(
        "wt", encoding="utf-8"
    ) as bed:
        table.write(
            "chromosome\tstart\tend\tsignificant_peak_count\tcluster_score\t"
            "max_peak_score\tsummit\tempirical_fdr\n"
        )
        for index, cluster in enumerate(target_clusters, 1):
            table.write(
                f"{cluster.chrom}\t{cluster.start}\t{cluster.end}\t"
                f"{len(cluster.significant_peaks)}\t{cluster.score:.12g}\t"
                f"{cluster.max_peak_score:.12g}\t{cluster.summit}\t{cluster.qvalue:.12g}\n"
            )
            if cluster.qvalue <= cluster_fdr:
                bed.write(
                    f"{cluster.chrom}\t{cluster.start}\t{cluster.end}\t"
                    f"chip_cluster_{index}\t{cluster.score:.6f}\t.\t"
                    f"{cluster.summit}\t{cluster.summit + 1}\t{cluster.qvalue:.12g}\n"
                )

    return {
        "annotated_peaks": annotated,
        "significant_peaks": significant,
        "cluster_table": cluster_table,
        "significant_clusters": significant_clusters,
    }
