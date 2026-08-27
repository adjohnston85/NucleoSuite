"""Target-control peak competition and empirical cluster FDR for cutn-suite."""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import ttest_ind

from nucleosuite.bigwig_ops import interval_max, open_bigwigs
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
    signal_score: float | None = None
    competition_score: float | None = None
    treatment_replicate_scores: tuple[float, ...] = ()
    control_replicate_scores: tuple[float, ...] = ()
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
    seed_peak_count: int = 0
    bridged_non_member_peak_count: int = 0
    minimum_seed_pvalue: float = math.nan


@dataclass(frozen=True)
class ReplicatePeakStatistics:
    peak: CompetitivePeak
    treatment_scores: tuple[float, ...]
    control_scores: tuple[float, ...]
    treatment_mean: float
    control_mean: float
    mean_difference: float
    minimum_treatment: float
    maximum_control: float
    conservative_excess: float
    conservative_fold_enrichment: float
    conservative_log2_enrichment: float
    all_controls_gate: bool
    pvalue: float
    qvalue: float = math.nan


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
            updated_target.append(
                replace(
                    peak,
                    winner=True,
                    matched_score=0.0,
                    signal_score=peak.row.score,
                    competition_score=peak.row.score,
                )
            )
            continue
        other = control[match]
        updated_target.append(
            replace(
                peak,
                winner=peak.row.score > other.row.score,
                matched_score=other.row.score,
                signal_score=peak.row.score,
                competition_score=max(peak.row.score - other.row.score, 0.0),
            )
        )
    updated_control: list[CompetitivePeak] = []
    for index, peak in enumerate(control):
        match = control_matches.get(index)
        if match is None:
            updated_control.append(
                replace(
                    peak,
                    winner=True,
                    matched_score=0.0,
                    signal_score=peak.row.score,
                    competition_score=peak.row.score,
                )
            )
            continue
        other = target[match]
        # Ties conservatively belong to the control.
        updated_control.append(
            replace(
                peak,
                winner=peak.row.score >= other.row.score,
                matched_score=other.row.score,
                signal_score=peak.row.score,
                competition_score=max(peak.row.score - other.row.score, 0.0),
            )
        )
    return updated_target, updated_control


def compete_peaks_with_bigwigs(
    target_rows: Sequence[PeakRow],
    control_rows: Sequence[PeakRow],
    *,
    target_bigwig: str | Path,
    control_bigwig: str | Path,
    summit_column: int = 7,
) -> tuple[list[CompetitivePeak], list[CompetitivePeak]]:
    """Compare both mean-scaled coverage tracks in each candidate interval.

    Target candidates use the maximum target and control values in the target
    peak interval. Control candidates use the reciprocal comparison in the
    control peak interval. This avoids treating the absence of a separately
    called overlapping control peak as zero control signal.
    """

    target = _to_competitive(target_rows, "target", summit_column)
    control = _to_competitive(control_rows, "control", summit_column)
    target_handle, control_handle = open_bigwigs([target_bigwig, control_bigwig])
    try:
        updated_target: list[CompetitivePeak] = []
        for peak in target:
            target_score = max(
                interval_max(target_handle, peak.chrom, peak.start, peak.end), 0.0
            )
            control_score = max(
                interval_max(control_handle, peak.chrom, peak.start, peak.end), 0.0
            )
            updated_target.append(
                replace(
                    peak,
                    winner=target_score > control_score,
                    matched_score=control_score,
                    signal_score=target_score,
                    competition_score=max(target_score - control_score, 0.0),
                )
            )

        updated_control: list[CompetitivePeak] = []
        for peak in control:
            target_score = max(
                interval_max(target_handle, peak.chrom, peak.start, peak.end), 0.0
            )
            control_score = max(
                interval_max(control_handle, peak.chrom, peak.start, peak.end), 0.0
            )
            updated_control.append(
                replace(
                    peak,
                    winner=control_score >= target_score,
                    matched_score=target_score,
                    signal_score=control_score,
                    competition_score=max(control_score - target_score, 0.0),
                )
            )
        return updated_target, updated_control
    finally:
        target_handle.close()
        control_handle.close()


def compete_peaks_all_controls(
    target_rows: Sequence[PeakRow],
    control_rows: Sequence[PeakRow],
    *,
    target_bigwigs: Sequence[str | Path],
    control_bigwigs: Sequence[str | Path],
    target_mean_bigwig: str | Path,
    control_mean_bigwig: str | Path,
    summit_column: int = 7,
) -> tuple[list[CompetitivePeak], list[CompetitivePeak]]:
    """Require every treatment replicate to exceed every control replicate.

    Replicate maxima are measured over the candidate coordinates.  Therefore a
    target candidate wins exactly when the minimum treatment-replicate maximum
    is greater than the maximum control-replicate maximum.  Control candidates
    are evaluated reciprocally to provide empirical decoys.  The condition-mean
    coverage maximum remains the reported BED score.
    """

    if not target_bigwigs or not control_bigwigs:
        raise ValueError("All-controls competition requires treatment and control tracks")
    target = _to_competitive(target_rows, "target", summit_column)
    control = _to_competitive(control_rows, "control", summit_column)
    paths = [
        *target_bigwigs,
        *control_bigwigs,
        target_mean_bigwig,
        control_mean_bigwig,
    ]
    handles = open_bigwigs(paths)
    target_handles = handles[: len(target_bigwigs)]
    control_handles = handles[
        len(target_bigwigs) : len(target_bigwigs) + len(control_bigwigs)
    ]
    target_mean_handle = handles[-2]
    control_mean_handle = handles[-1]

    def maxima(peak: CompetitivePeak) -> tuple[tuple[float, ...], tuple[float, ...]]:
        treatment_scores = tuple(
            max(interval_max(handle, peak.chrom, peak.start, peak.end), 0.0)
            for handle in target_handles
        )
        control_scores = tuple(
            max(interval_max(handle, peak.chrom, peak.start, peak.end), 0.0)
            for handle in control_handles
        )
        return treatment_scores, control_scores

    try:
        updated_target: list[CompetitivePeak] = []
        for peak in target:
            treatment_scores, control_scores = maxima(peak)
            minimum_treatment = min(treatment_scores)
            maximum_control = max(control_scores)
            updated_target.append(
                replace(
                    peak,
                    winner=minimum_treatment > maximum_control,
                    matched_score=maximum_control,
                    signal_score=max(
                        interval_max(
                            target_mean_handle, peak.chrom, peak.start, peak.end
                        ),
                        0.0,
                    ),
                    competition_score=max(minimum_treatment - maximum_control, 0.0),
                    treatment_replicate_scores=treatment_scores,
                    control_replicate_scores=control_scores,
                )
            )

        updated_control: list[CompetitivePeak] = []
        for peak in control:
            treatment_scores, control_scores = maxima(peak)
            maximum_treatment = max(treatment_scores)
            minimum_control = min(control_scores)
            updated_control.append(
                replace(
                    peak,
                    winner=minimum_control >= maximum_treatment,
                    matched_score=maximum_treatment,
                    signal_score=max(
                        interval_max(
                            control_mean_handle, peak.chrom, peak.start, peak.end
                        ),
                        0.0,
                    ),
                    competition_score=max(minimum_control - maximum_treatment, 0.0),
                    treatment_replicate_scores=treatment_scores,
                    control_replicate_scores=control_scores,
                )
            )
        return updated_target, updated_control
    finally:
        for handle in handles:
            handle.close()


def _competition_score(peak: CompetitivePeak) -> float:
    return (
        peak.row.score
        if peak.competition_score is None
        else float(peak.competition_score)
    )


def _signal_score(peak: CompetitivePeak) -> float:
    return peak.row.score if peak.signal_score is None else float(peak.signal_score)


def assign_competition_qvalues(
    target: Sequence[CompetitivePeak], control: Sequence[CompetitivePeak]
) -> tuple[list[CompetitivePeak], float | None]:
    """Assign target-decoy q-values to target-winning peaks."""

    target_indices = [index for index, peak in enumerate(target) if peak.winner]
    control_scores = [_competition_score(peak) for peak in control if peak.winner]
    if not target_indices:
        return list(target), None
    qvalues, _sample_counts, _control_counts = empirical_peak_qvalues(
        [_competition_score(target[index]) for index in target_indices],
        [control_scores],
    )
    output = list(target)
    for index, qvalue in zip(target_indices, qvalues):
        output[index] = replace(output[index], qvalue=float(qvalue))
    return output, min(
        (_competition_score(output[index]) for index in target_indices), default=None
    )


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
        if not peak.winner or _competition_score(peak) < significance_score:
            return False
        return require_qvalue is None or peak.qvalue <= require_qvalue

    def finish() -> None:
        nonlocal current_significant, nonsignificant_run
        if len(current_significant) >= minimum_significant_peaks:
            best = max(
                current_significant,
                key=lambda peak: (_competition_score(peak), -peak.summit),
            )
            score = float(
                sum(
                    max(_competition_score(peak) - significance_score, 0.0)
                    for peak in current_significant
                )
            )
            clusters.append(
                PeakCluster(
                    chrom=current_significant[0].chrom,
                    start=min(peak.start for peak in current_significant),
                    end=max(peak.end for peak in current_significant),
                    significant_peaks=tuple(current_significant),
                    score=score,
                    summit=best.summit,
                    max_peak_score=max(_signal_score(peak) for peak in current_significant),
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


def cluster_seeded_gate_peaks(
    statistics: Sequence[ReplicatePeakStatistics],
    *,
    seed_pvalue: float = 0.05,
    gate_mode: str = "all-controls",
    member_mode: str = "seed-and-gated",
    maximum_non_member_gap: int = 1,
    max_cluster_gap: int = 1000,
    minimum_member_peaks: int = 2,
) -> list[PeakCluster]:
    """Build seeded clusters using the selected treatment-control gate.

    ``gate_mode='all-controls'`` uses the conservative minimum-treatment >
    maximum-control rule and is the default. ``mean`` uses mean treatment >
    mean control and scores each member by mean treatment minus mean control. In
    ``seed-and-gated`` mode, significant seeds and other gate-passing peaks are
    members. In ``significant-only`` mode, only significant seeds are members.
    Up to ``maximum_non_member_gap`` consecutive non-members may bridge two
    included members; a longer run ends the current cluster.
    """
    if not math.isfinite(seed_pvalue) or not 0 <= seed_pvalue <= 1:
        raise ValueError("Seed p-value must be between 0 and 1")
    if gate_mode not in {"mean", "all-controls"}:
        raise ValueError("gate_mode must be mean or all-controls")
    if member_mode not in {"seed-and-gated", "significant-only"}:
        raise ValueError("member_mode must be seed-and-gated or significant-only")
    if maximum_non_member_gap < 0:
        raise ValueError("Maximum non-member gap must be non-negative")
    if max_cluster_gap < 1 or minimum_member_peaks < 1:
        raise ValueError("Maximum cluster gap and minimum member count must be positive")

    def passes_gate(record: ReplicatePeakStatistics) -> bool:
        return record.mean_difference > 0 if gate_mode == "mean" else record.all_controls_gate

    def excess(record: ReplicatePeakStatistics) -> float:
        return record.mean_difference if gate_mode == "mean" else record.conservative_excess

    def is_seed(record: ReplicatePeakStatistics) -> bool:
        return passes_gate(record) and math.isfinite(record.pvalue) and record.pvalue < seed_pvalue

    def is_member(record: ReplicatePeakStatistics) -> bool:
        return is_seed(record) if member_mode == "significant-only" else passes_gate(record)

    ordered = sorted(
        statistics,
        key=lambda record: (record.peak.chrom, record.peak.summit, record.peak.start, record.peak.end),
    )
    clusters: list[PeakCluster] = []
    current: list[ReplicatePeakStatistics] = []
    pending_non_members = 0
    bridged_non_members = 0
    current_chrom: str | None = None

    def finish() -> None:
        nonlocal current, pending_non_members, bridged_non_members
        seeds = [record for record in current if is_seed(record)]
        if len(current) >= minimum_member_peaks and seeds:
            best = max(
                current,
                key=lambda record: (_signal_score(record.peak), excess(record), -record.peak.summit),
            )
            finite_seed_qvalues = [record.qvalue for record in seeds if math.isfinite(record.qvalue)]
            clusters.append(
                PeakCluster(
                    chrom=current[0].peak.chrom,
                    start=min(record.peak.start for record in current),
                    end=max(record.peak.end for record in current),
                    significant_peaks=tuple(record.peak for record in current),
                    score=float(sum(excess(record) for record in current)),
                    summit=best.peak.summit,
                    max_peak_score=max(_signal_score(record.peak) for record in current),
                    qvalue=max(finite_seed_qvalues) if finite_seed_qvalues else math.nan,
                    seed_peak_count=len(seeds),
                    bridged_non_member_peak_count=bridged_non_members,
                    minimum_seed_pvalue=min(record.pvalue for record in seeds),
                )
            )
        current = []
        pending_non_members = 0
        bridged_non_members = 0

    for record in ordered:
        peak = record.peak
        if current_chrom is not None and peak.chrom != current_chrom:
            finish()
        current_chrom = peak.chrom
        if is_member(record):
            if current and peak.summit - current[-1].peak.summit > max_cluster_gap:
                finish()
            if current:
                bridged_non_members += pending_non_members
                current.append(record)
            else:
                current = [record]
            pending_non_members = 0
        elif current:
            pending_non_members += 1
            if pending_non_members > maximum_non_member_gap:
                finish()
    finish()
    return clusters


def analyze_cutn_peaks(
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
    target_bigwig: str | Path | None = None,
    control_bigwig: str | Path | None = None,
    target_replicate_bigwigs: Sequence[str | Path] | None = None,
    control_replicate_bigwigs: Sequence[str | Path] | None = None,
    control_mode: str = "condition-mean",
) -> dict[str, Path]:
    """Run competition, peak FDR and cluster FDR and write BED/TSV outputs."""

    for name, value in (("peak_fdr", peak_fdr), ("cluster_fdr", cluster_fdr)):
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{name} must be between 0 and 1")
    target_rows = _read_peak_rows(target_path, score_column, allow_empty=True)
    control_rows = _read_peak_rows(control_path, score_column, allow_empty=True)
    if control_mode not in {"all-controls", "condition-mean"}:
        raise ValueError("control_mode must be all-controls or condition-mean")
    if (target_bigwig is None) != (control_bigwig is None):
        raise ValueError("target_bigwig and control_bigwig must be supplied together")
    if control_mode == "all-controls":
        if target_bigwig is None or control_bigwig is None:
            raise ValueError("all-controls mode requires condition-mean BigWigs")
        if not target_replicate_bigwigs or not control_replicate_bigwigs:
            raise ValueError("all-controls mode requires replicate BigWigs")
        target, control = compete_peaks_all_controls(
            target_rows,
            control_rows,
            target_bigwigs=target_replicate_bigwigs,
            control_bigwigs=control_replicate_bigwigs,
            target_mean_bigwig=target_bigwig,
            control_mean_bigwig=control_bigwig,
            summit_column=summit_column,
        )
    elif target_bigwig is not None and control_bigwig is not None:
        target, control = compete_peaks_with_bigwigs(
            target_rows,
            control_rows,
            target_bigwig=target_bigwig,
            control_bigwig=control_bigwig,
            summit_column=summit_column,
        )
    else:
        target, control = compete_peaks(
            target_rows, control_rows,
            match_distance=match_distance,
            summit_column=summit_column,
        )
    target, _ = assign_competition_qvalues(target, control)
    passing = [peak for peak in target if peak.winner and peak.qvalue <= peak_fdr]
    threshold = min((_competition_score(peak) for peak in passing), default=math.inf)

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    annotated = directory / "target_peaks_empirical_fdr.bed"
    significant = directory / f"target_peaks_fdr{peak_fdr:g}_significant.bed"
    competition_table = directory / "target_control_interval_competition.tsv"
    with annotated.open("wt", encoding="utf-8") as all_handle, significant.open(
        "wt", encoding="utf-8"
    ) as sig_handle:
        for peak in sorted(target, key=lambda value: (value.chrom, value.start, value.end)):
            fields = list(peak.row.fields)
            fields[score_column - 1] = f"{_signal_score(peak):.12g}"
            text = "\t".join((*fields, f"{peak.qvalue:.12g}")) + "\n"
            all_handle.write(text)
            if peak.winner and peak.qvalue <= peak_fdr:
                sig_handle.write(text)
    with competition_table.open("wt", encoding="utf-8") as handle:
        handle.write(
            "source\tchromosome\tstart\tend\tsummit\tdiscovery_score\t"
            "scaled_coverage_max\t"
            "other_track_interval_max\tlocal_excess\twinner\tempirical_fdr\t"
            "treatment_replicate_maxima\tcontrol_replicate_maxima\n"
        )
        for peak in [*target, *control]:
            treatment_scores = ";".join(
                f"{value:.12g}" for value in peak.treatment_replicate_scores
            )
            control_scores = ";".join(
                f"{value:.12g}" for value in peak.control_replicate_scores
            )
            handle.write(
                f"{peak.source}\t{peak.chrom}\t{peak.start}\t{peak.end}\t"
                f"{peak.summit}\t{peak.row.score:.12g}\t{_signal_score(peak):.12g}\t"
                f"{peak.matched_score:.12g}\t"
                f"{_competition_score(peak):.12g}\t{str(peak.winner).lower()}\t"
                f"{peak.qvalue:.12g}\t{treatment_scores}\t{control_scores}\n"
            )

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
            "chromosome\tstart\tend\tqualifying_peak_count\tcluster_score\t"
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
                    f"cutn_cluster_{index}\t{cluster.score:.6f}\t.\t"
                    f"{cluster.summit}\t{cluster.summit + 1}\t{cluster.qvalue:.12g}\n"
                )

    return {
        "annotated_peaks": annotated,
        "significant_peaks": significant,
        "competition_table": competition_table,
        "cluster_table": cluster_table,
        "significant_clusters": significant_clusters,
    }


def _one_sided_welch_greater(
    treatment: Sequence[float], control: Sequence[float]
) -> float:
    """Test whether the treatment mean exceeds the control mean."""

    if len(treatment) < 2 or len(control) < 2:
        return math.nan
    treatment_variance = float(np.var(treatment, ddof=1))
    control_variance = float(np.var(control, ddof=1))
    treatment_mean = float(np.mean(treatment))
    control_mean = float(np.mean(control))
    if treatment_variance == 0 and control_variance == 0:
        return 0.0 if treatment_mean > control_mean else 1.0
    result = ttest_ind(
        treatment,
        control,
        equal_var=False,
        nan_policy="omit",
        alternative="greater",
    )
    value = float(result.pvalue)
    return value if math.isfinite(value) else 1.0


def _bh_qvalues(pvalues: Sequence[float]) -> list[float]:
    qvalues = [math.nan] * len(pvalues)
    finite = [(index, value) for index, value in enumerate(pvalues) if math.isfinite(value)]
    if not finite:
        return qvalues
    ordered = sorted(finite, key=lambda item: item[1])
    running = 1.0
    total = len(ordered)
    for rank in range(total, 0, -1):
        index, pvalue = ordered[rank - 1]
        running = min(running, max(0.0, min(1.0, pvalue * total / rank)))
        qvalues[index] = running
    return qvalues


def _format_optional(value: float) -> str:
    return "." if not math.isfinite(value) else f"{value:.12g}"


def analyze_cutn_replicate_peaks(
    target_path: str | Path,
    *,
    output_dir: str | Path,
    target_replicate_bigwigs: Sequence[str | Path],
    control_replicate_bigwigs: Sequence[str | Path],
    target_mean_bigwig: str | Path,
    score_column: int = 5,
    summit_column: int = 7,
    peak_pvalue: float | None = None,
    peak_fdr: float | None = None,
    cluster_seed_pvalue: float = 0.05,
    cluster_fdr: float | None = None,
    gate_mode: str = "all-controls",
    cluster_member_mode: str = "seed-and-gated",
    cluster_max_non_member_gap: int = 1,
    max_cluster_gap: int = 1000,
    minimum_cluster_members: int = 2,
) -> dict[str, Path]:
    """Test treatment-defined candidates from replicate normalized coverage.

    Candidate maxima are measured in every treatment and control replicate.
    The default all-controls gate requires every treatment replicate to exceed
    every control replicate. The optional mean gate requires mean treatment >
    mean control. Welch p-values and BH FDR are reported independently of the gate.
    """
    for name, value in (
        ("peak_pvalue", peak_pvalue), ("peak_fdr", peak_fdr),
        ("cluster_seed_pvalue", cluster_seed_pvalue), ("cluster_fdr", cluster_fdr),
    ):
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f"{name} must be between 0 and 1")
    if gate_mode not in {"mean", "all-controls"}:
        raise ValueError("gate_mode must be mean or all-controls")
    if cluster_member_mode not in {"seed-and-gated", "significant-only"}:
        raise ValueError("cluster_member_mode must be seed-and-gated or significant-only")
    if cluster_max_non_member_gap < 0:
        raise ValueError("cluster_max_non_member_gap must be non-negative")
    if minimum_cluster_members < 1:
        raise ValueError("minimum_cluster_members must be positive")
    if not target_replicate_bigwigs or not control_replicate_bigwigs:
        raise ValueError("Treatment and control replicate BigWigs are required")

    target_rows = _read_peak_rows(target_path, score_column, allow_empty=True)
    peaks = _to_competitive(target_rows, "target", summit_column)
    paths = [*target_replicate_bigwigs, *control_replicate_bigwigs, target_mean_bigwig]
    handles = open_bigwigs(paths)
    treatment_handles = handles[:len(target_replicate_bigwigs)]
    control_handles = handles[len(target_replicate_bigwigs):len(target_replicate_bigwigs)+len(control_replicate_bigwigs)]
    mean_handle = handles[-1]
    statistics: list[ReplicatePeakStatistics] = []
    try:
        for peak in peaks:
            treatment_scores = tuple(max(interval_max(h, peak.chrom, peak.start, peak.end), 0.0) for h in treatment_handles)
            control_scores = tuple(max(interval_max(h, peak.chrom, peak.start, peak.end), 0.0) for h in control_handles)
            treatment_mean = float(np.mean(treatment_scores))
            control_mean = float(np.mean(control_scores))
            mean_difference = treatment_mean - control_mean
            minimum_treatment = min(treatment_scores)
            maximum_control = max(control_scores)
            conservative_difference = minimum_treatment - maximum_control
            conservative_excess = max(conservative_difference, 0.0)
            fold_enrichment = (minimum_treatment + 1.0) / (maximum_control + 1.0)
            all_controls_gate = minimum_treatment > maximum_control
            selected_gate = mean_difference > 0 if gate_mode == "mean" else all_controls_gate
            selected_excess = mean_difference if gate_mode == "mean" else conservative_difference
            mean_score = max(interval_max(mean_handle, peak.chrom, peak.start, peak.end), 0.0)
            measured_peak = replace(
                peak,
                winner=selected_gate,
                matched_score=control_mean if gate_mode == "mean" else maximum_control,
                signal_score=mean_score,
                competition_score=max(selected_excess, 0.0),
                treatment_replicate_scores=treatment_scores,
                control_replicate_scores=control_scores,
            )
            statistics.append(ReplicatePeakStatistics(
                peak=measured_peak,
                treatment_scores=treatment_scores,
                control_scores=control_scores,
                treatment_mean=treatment_mean,
                control_mean=control_mean,
                mean_difference=mean_difference,
                minimum_treatment=minimum_treatment,
                maximum_control=maximum_control,
                conservative_excess=conservative_excess,
                conservative_fold_enrichment=fold_enrichment,
                conservative_log2_enrichment=math.log2(fold_enrichment),
                all_controls_gate=all_controls_gate,
                pvalue=_one_sided_welch_greater(treatment_scores, control_scores),
            ))
    finally:
        for handle in handles:
            handle.close()

    qvalues = _bh_qvalues([record.pvalue for record in statistics])
    statistics = [replace(record, qvalue=q, peak=replace(record.peak, qvalue=q)) for record, q in zip(statistics, qvalues)]

    def selected_gate(record: ReplicatePeakStatistics) -> bool:
        return record.mean_difference > 0 if gate_mode == "mean" else record.all_controls_gate

    def selected_excess(record: ReplicatePeakStatistics) -> float:
        return record.mean_difference if gate_mode == "mean" else record.conservative_excess

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    annotated = directory / f"target_peaks_replicate_statistics_gate_{gate_mode}.bed"
    selection_tokens = [f"gate_{gate_mode}"]
    if peak_pvalue is not None: selection_tokens.append(f"p{peak_pvalue:g}")
    if peak_fdr is not None: selection_tokens.append(f"fdr{peak_fdr:g}")
    significant = directory / f"target_peaks_{'_'.join(selection_tokens)}.bed"
    seed_tokens = [f"gate_{gate_mode}", f"seed_p{cluster_seed_pvalue:g}"]
    seed_peaks = directory / f"target_seed_peaks_{'_'.join(seed_tokens)}.bed"
    statistics_table = directory / f"target_peak_replicate_statistics_gate_{gate_mode}.tsv"
    with (
        annotated.open("wt", encoding="utf-8") as all_handle,
        significant.open("wt", encoding="utf-8") as sig_handle,
        seed_peaks.open("wt", encoding="utf-8") as seed_handle,
    ):
        for record in statistics:
            fields = list(record.peak.row.fields)
            fields[score_column - 1] = f"{_signal_score(record.peak):.12g}"
            text = "\t".join((
                *fields,
                _format_optional(record.pvalue),
                _format_optional(record.qvalue),
            )) + "\n"
            all_handle.write(text)
            passes_p = peak_pvalue is None or (math.isfinite(record.pvalue) and record.pvalue <= peak_pvalue)
            passes_q = peak_fdr is None or (math.isfinite(record.qvalue) and record.qvalue <= peak_fdr)
            if selected_gate(record) and passes_p and passes_q:
                sig_handle.write(text)
            if (
                selected_gate(record)
                and math.isfinite(record.pvalue)
                and record.pvalue < cluster_seed_pvalue
            ):
                seed_handle.write(text)

    with statistics_table.open("wt", encoding="utf-8") as handle:
        handle.write(
            "chromosome	start	end	summit	discovery_score	condition_mean_scaled_coverage_max	"
            "treatment_replicate_maxima	control_replicate_maxima	treatment_mean	control_mean	"
            "treatment_minus_control	minimum_treatment	maximum_control	conservative_excess	"
            "conservative_fold_enrichment_pseudocount1	conservative_log2_enrichment_pseudocount1	"
            "mean_treatment_exceeds_mean_control	all_treatments_exceed_all_controls	gate_mode	"
            "selected_gate\tselected_excess\tselected_for_stage2\tp_value\tfdr\n"
        )
        for record in statistics:
            passes_p = peak_pvalue is None or (math.isfinite(record.pvalue) and record.pvalue <= peak_pvalue)
            passes_q = peak_fdr is None or (math.isfinite(record.qvalue) and record.qvalue <= peak_fdr)
            selected = selected_gate(record) and passes_p and passes_q
            peak = record.peak
            handle.write(
                f"{peak.chrom}	{peak.start}	{peak.end}	{peak.summit}	{peak.row.score:.12g}	{_signal_score(peak):.12g}	"
                + ";".join(f"{v:.12g}" for v in record.treatment_scores) + "	"
                + ";".join(f"{v:.12g}" for v in record.control_scores)
                + f"	{record.treatment_mean:.12g}	{record.control_mean:.12g}	{record.mean_difference:.12g}	"
                f"{record.minimum_treatment:.12g}	{record.maximum_control:.12g}	{record.conservative_excess:.12g}	"
                f"{record.conservative_fold_enrichment:.12g}	{record.conservative_log2_enrichment:.12g}	"
                f"{str(record.mean_difference > 0).lower()}	{str(record.all_controls_gate).lower()}	{gate_mode}	"
                f"{str(selected_gate(record)).lower()}	{selected_excess(record):.12g}	{str(selected).lower()}	"
                f"{_format_optional(record.pvalue)}\t{_format_optional(record.qvalue)}\n"
            )

    clusters = cluster_seeded_gate_peaks(
        statistics,
        seed_pvalue=cluster_seed_pvalue,
        gate_mode=gate_mode,
        member_mode=cluster_member_mode,
        maximum_non_member_gap=cluster_max_non_member_gap,
        max_cluster_gap=max_cluster_gap,
        minimum_member_peaks=minimum_cluster_members,
    )
    cluster_table = directory / "target_clusters_seeded.tsv"
    cluster_tokens = [f"gate_{gate_mode}", cluster_member_mode, f"seed_p{cluster_seed_pvalue:g}", f"gap{cluster_max_non_member_gap}", f"min{minimum_cluster_members}"]
    if cluster_fdr is not None: cluster_tokens.append(f"maximum_seed_fdr{cluster_fdr:g}")
    significant_clusters = directory / f"target_clusters_{'_'.join(cluster_tokens)}.bed"
    with cluster_table.open("wt", encoding="utf-8") as table, significant_clusters.open("wt", encoding="utf-8") as bed:
        table.write(
            "cluster_id	chromosome	start	end	seed_peak_count	member_count	bridged_non_member_peak_count	"
            "cluster_score\tmax_peak_score\tstrongest_peak_summit\tminimum_seed_p_value\tmaximum_seed_fdr\n"
        )
        for index, cluster in enumerate(clusters, 1):
            cluster_id = f"cutn_cluster_{index}"
            table.write(
                f"{cluster_id}	{cluster.chrom}	{cluster.start}	{cluster.end}	{cluster.seed_peak_count}	"
                f"{len(cluster.significant_peaks)}	{cluster.bridged_non_member_peak_count}	{cluster.score:.12g}	"
                f"{cluster.max_peak_score:.12g}\t{cluster.summit}\t"
                f"{_format_optional(cluster.minimum_seed_pvalue)}\t{_format_optional(cluster.qvalue)}\n"
            )
            passes_cluster_fdr = cluster_fdr is None or (math.isfinite(cluster.qvalue) and cluster.qvalue <= cluster_fdr)
            if passes_cluster_fdr:
                bed.write(
                    f"{cluster.chrom}	{cluster.start}	{cluster.end}	{cluster_id}	{cluster.score:.6f}	.	"
                    f"{cluster.summit}\t{cluster.summit + 1}\t"
                    f"{_format_optional(cluster.qvalue)}\t{cluster.max_peak_score:.12g}\n"
                )
    return {
        "annotated_peaks": annotated, "selected_peaks": significant, "significant_peaks": significant,
        "seed_peaks": seed_peaks,
        "competition_table": statistics_table, "cluster_table": cluster_table,
        "selected_clusters": significant_clusters, "significant_clusters": significant_clusters,
    }

