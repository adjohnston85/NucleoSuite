"""Target-control peak competition and empirical cluster FDR for cutn-suite."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import ttest_ind

from nucleosuite.bigwig_ops import interval_max, interval_mean, open_bigwigs
from nucleosuite.peak_fdr import PeakRow, _read_peak_rows


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










def _signal_score(peak: CompetitivePeak) -> float:
    return peak.row.score if peak.signal_score is None else float(peak.signal_score)








def cluster_seeded_gate_peaks(
    statistics: Sequence[ReplicatePeakStatistics],
    *,
    seed_pvalue: float = 0.05,
    seed_mode: str = "pvalue",
    seed_gate_mode: str = "mean",
    member_gate_mode: str = "all-controls",
    member_mode: str = "seed-and-gated",
    maximum_non_member_gap: int = 1,
    max_cluster_gap: int = 1000,
    minimum_member_peaks: int = 2,
) -> list[PeakCluster]:
    """Build clusters from independently defined seed (S) and member (G) rules."""

    if not math.isfinite(seed_pvalue) or not 0 <= seed_pvalue <= 1:
        raise ValueError("Seed p-value must be between 0 and 1")
    if seed_mode not in {"pvalue", "gated"}:
        raise ValueError("seed_mode must be pvalue or gated")
    if seed_gate_mode not in {"mean", "all-controls"}:
        raise ValueError("seed_gate_mode must be mean or all-controls")
    if member_gate_mode not in {"mean", "all-controls"}:
        raise ValueError("member_gate_mode must be mean or all-controls")
    if member_mode not in {"seed-and-gated", "significant-only"}:
        raise ValueError("member_mode must be seed-and-gated or significant-only")
    if maximum_non_member_gap < 0:
        raise ValueError("Maximum non-member gap must be non-negative")
    if max_cluster_gap < 1 or minimum_member_peaks < 1:
        raise ValueError("Maximum cluster gap and minimum member count must be positive")

    def passes_gate(record: ReplicatePeakStatistics, mode: str) -> bool:
        return record.mean_difference > 0 if mode == "mean" else record.all_controls_gate

    def gate_excess(record: ReplicatePeakStatistics, mode: str) -> float:
        return record.mean_difference if mode == "mean" else record.conservative_excess

    def is_seed(record: ReplicatePeakStatistics) -> bool:
        if not passes_gate(record, seed_gate_mode):
            return False
        if seed_mode == "gated":
            return True
        return math.isfinite(record.pvalue) and record.pvalue < seed_pvalue

    def is_member(record: ReplicatePeakStatistics) -> bool:
        if member_mode == "significant-only":
            return is_seed(record)
        return is_seed(record) or passes_gate(record, member_gate_mode)

    def member_excess(record: ReplicatePeakStatistics) -> float:
        if is_seed(record) and not passes_gate(record, member_gate_mode):
            return max(gate_excess(record, seed_gate_mode), 0.0)
        return max(gate_excess(record, member_gate_mode), 0.0)

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
                key=lambda record: (_signal_score(record.peak), member_excess(record), -record.peak.summit),
            )
            finite_seed_p = [record.pvalue for record in seeds if math.isfinite(record.pvalue)]
            clusters.append(
                PeakCluster(
                    chrom=current[0].peak.chrom,
                    start=min(record.peak.start for record in current),
                    end=max(record.peak.end for record in current),
                    significant_peaks=tuple(record.peak for record in current),
                    score=float(sum(member_excess(record) for record in current)),
                    summit=best.peak.summit,
                    max_peak_score=max(_signal_score(record.peak) for record in current),
                    qvalue=math.nan,
                    seed_peak_count=len(seeds),
                    bridged_non_member_peak_count=bridged_non_members,
                    minimum_seed_pvalue=min(finite_seed_p) if finite_seed_p else math.nan,
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
    cluster_seed_pvalue: float = 0.05,
    seed_mode: str = "pvalue",
    seed_gate_mode: str = "mean",
    member_gate_mode: str = "all-controls",
    compute_pvalues: bool = True,
    coverage_statistic: str = "mean",
    cluster_member_mode: str = "seed-and-gated",
    cluster_max_non_member_gap: int = 1,
    max_cluster_gap: int = 1000,
    minimum_cluster_members: int = 2,
) -> dict[str, Path]:
    """Measure treatment-defined candidates across replicate coverage tracks.

    Mean interval coverage is the default measurement. Seed (S) and member (G)
    gates are independent. Raw one-sided Welch p-values are calculated when
    requested and are used only by p-value seed mode.
    """
    if not math.isfinite(cluster_seed_pvalue) or not 0 <= cluster_seed_pvalue <= 1:
        raise ValueError("cluster_seed_pvalue must be between 0 and 1")
    if seed_mode not in {"pvalue", "gated"}:
        raise ValueError("seed_mode must be pvalue or gated")
    if seed_gate_mode not in {"mean", "all-controls"}:
        raise ValueError("seed_gate_mode must be mean or all-controls")
    if member_gate_mode not in {"mean", "all-controls"}:
        raise ValueError("member_gate_mode must be mean or all-controls")
    if coverage_statistic not in {"mean", "max"}:
        raise ValueError("coverage_statistic must be mean or max")
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
    measure = interval_mean if coverage_statistic == "mean" else interval_max
    statistics: list[ReplicatePeakStatistics] = []
    try:
        for peak in peaks:
            treatment_scores = tuple(max(measure(h, peak.chrom, peak.start, peak.end), 0.0) for h in treatment_handles)
            control_scores = tuple(max(measure(h, peak.chrom, peak.start, peak.end), 0.0) for h in control_handles)
            treatment_mean = float(np.mean(treatment_scores))
            control_mean = float(np.mean(control_scores))
            mean_difference = treatment_mean - control_mean
            minimum_treatment = min(treatment_scores)
            maximum_control = max(control_scores)
            conservative_difference = minimum_treatment - maximum_control
            conservative_excess = max(conservative_difference, 0.0)
            fold_enrichment = (minimum_treatment + 1.0) / (maximum_control + 1.0)
            all_controls_gate = minimum_treatment > maximum_control
            mean_score = max(measure(mean_handle, peak.chrom, peak.start, peak.end), 0.0)
            pvalue = (
                _one_sided_welch_greater(treatment_scores, control_scores)
                if compute_pvalues else math.nan
            )
            measured_peak = replace(
                peak,
                signal_score=mean_score,
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
                pvalue=pvalue,
                qvalue=math.nan,
            ))
    finally:
        for handle in handles:
            handle.close()

    def passes_gate(record: ReplicatePeakStatistics, mode: str) -> bool:
        return record.mean_difference > 0 if mode == "mean" else record.all_controls_gate

    def is_seed(record: ReplicatePeakStatistics) -> bool:
        if not passes_gate(record, seed_gate_mode):
            return False
        if seed_mode == "gated":
            return True
        return math.isfinite(record.pvalue) and record.pvalue < cluster_seed_pvalue

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    annotated = directory / "target_peaks_replicate_statistics.bed"
    selected = directory / f"target_cluster_eligible_peaks_S-{seed_mode}-{seed_gate_mode}_G-{member_gate_mode}.bed"
    seed_peaks = directory / f"target_seed_peaks_S-{seed_mode}-{seed_gate_mode}.bed"
    statistics_table = directory / "target_peak_replicate_statistics.tsv"
    with annotated.open("wt", encoding="utf-8") as all_handle, selected.open("wt", encoding="utf-8") as selected_handle, seed_peaks.open("wt", encoding="utf-8") as seed_handle:
        for record in statistics:
            fields = list(record.peak.row.fields)
            fields[score_column - 1] = f"{_signal_score(record.peak):.12g}"
            text = "\t".join((*fields, _format_optional(record.pvalue))) + "\n"
            all_handle.write(text)
            seed = is_seed(record)
            member = passes_gate(record, member_gate_mode)
            if seed or member:
                selected_handle.write(text)
            if seed:
                seed_handle.write(text)

    with statistics_table.open("wt", encoding="utf-8") as handle:
        handle.write(
            "chromosome\tstart\tend\tsummit\tdiscovery_score\tcondition_mean_coverage_"
            f"{coverage_statistic}\ttreatment_replicate_{coverage_statistic}s\tcontrol_replicate_{coverage_statistic}s\t"
            "treatment_mean\tcontrol_mean\ttreatment_minus_control\tminimum_treatment\tmaximum_control\t"
            "conservative_excess\tconservative_fold_enrichment_pseudocount1\tconservative_log2_enrichment_pseudocount1\t"
            "mean_treatment_exceeds_mean_control\tall_treatments_exceed_all_controls\tseed_mode\tseed_gate_mode\tmember_gate_mode\t"
            "is_seed\tis_gated_member\tp_value\n"
        )
        for record in statistics:
            peak = record.peak
            handle.write(
                f"{peak.chrom}\t{peak.start}\t{peak.end}\t{peak.summit}\t{peak.row.score:.12g}\t{_signal_score(peak):.12g}\t"
                + ";".join(f"{v:.12g}" for v in record.treatment_scores) + "\t"
                + ";".join(f"{v:.12g}" for v in record.control_scores)
                + f"\t{record.treatment_mean:.12g}\t{record.control_mean:.12g}\t{record.mean_difference:.12g}\t"
                f"{record.minimum_treatment:.12g}\t{record.maximum_control:.12g}\t{record.conservative_excess:.12g}\t"
                f"{record.conservative_fold_enrichment:.12g}\t{record.conservative_log2_enrichment:.12g}\t"
                f"{str(record.mean_difference > 0).lower()}\t{str(record.all_controls_gate).lower()}\t"
                f"{seed_mode}\t{seed_gate_mode}\t{member_gate_mode}\t{str(is_seed(record)).lower()}\t"
                f"{str(passes_gate(record, member_gate_mode)).lower()}\t{_format_optional(record.pvalue)}\n"
            )

    clusters = cluster_seeded_gate_peaks(
        statistics,
        seed_pvalue=cluster_seed_pvalue,
        seed_mode=seed_mode,
        seed_gate_mode=seed_gate_mode,
        member_gate_mode=member_gate_mode,
        member_mode=cluster_member_mode,
        maximum_non_member_gap=cluster_max_non_member_gap,
        max_cluster_gap=max_cluster_gap,
        minimum_member_peaks=minimum_cluster_members,
    )
    cluster_table = directory / "target_clusters_seeded.tsv"
    cluster_tokens = [f"S-{seed_mode}-{seed_gate_mode}", f"G-{member_gate_mode}", cluster_member_mode, f"gap{cluster_max_non_member_gap}", f"min{minimum_cluster_members}"]
    if seed_mode == "pvalue":
        cluster_tokens.insert(1, f"p{cluster_seed_pvalue:g}")
    selected_clusters = directory / f"target_clusters_{'_'.join(cluster_tokens)}.bed"
    with cluster_table.open("wt", encoding="utf-8") as table, selected_clusters.open("wt", encoding="utf-8") as bed:
        table.write(
            "cluster_id\tchromosome\tstart\tend\tseed_peak_count\tmember_count\tbridged_non_member_peak_count\t"
            "cluster_score\tmax_peak_score\tstrongest_peak_summit\tminimum_seed_p_value\n"
        )
        for index, cluster in enumerate(clusters, 1):
            cluster_id = f"cutn_cluster_{index}"
            table.write(
                f"{cluster_id}\t{cluster.chrom}\t{cluster.start}\t{cluster.end}\t{cluster.seed_peak_count}\t"
                f"{len(cluster.significant_peaks)}\t{cluster.bridged_non_member_peak_count}\t{cluster.score:.12g}\t"
                f"{cluster.max_peak_score:.12g}\t{cluster.summit}\t{_format_optional(cluster.minimum_seed_pvalue)}\n"
            )
            bed.write(
                f"{cluster.chrom}\t{cluster.start}\t{cluster.end}\t{cluster_id}\t{cluster.score:.6f}\t.\t"
                f"{cluster.summit}\t{cluster.summit + 1}\t{cluster.max_peak_score:.12g}\n"
            )
    return {
        "annotated_peaks": annotated,
        "selected_peaks": selected,
        "significant_peaks": selected,
        "seed_peaks": seed_peaks,
        "competition_table": statistics_table,
        "cluster_table": cluster_table,
        "selected_clusters": selected_clusters,
        "significant_clusters": selected_clusters,
    }
