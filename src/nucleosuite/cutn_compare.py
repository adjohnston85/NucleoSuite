#!/usr/bin/env python3
"""Differential comparison of two completed cutn-suite Stage 1 analyses."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import brentq
from scipy.special import digamma, polygamma
from scipy.stats import t as student_t

from nucleosuite.bigwig_ops import (
    bigwig_chroms,
    interval_mean,
    open_bigwigs,
)
from nucleosuite.cutn_aggregate import (
    common_symmetric_bigwig_limit,
    run_cluster_aggregate,
)


MANIFEST_NAME = "cutn_stage1_manifest.json"
MANIFEST_SCHEMA = "nucleosuite_cutn_stage1"
MANIFEST_SCHEMA_VERSION = 6


@dataclass(frozen=True)
class ClusterRecord:
    cluster_id: str
    chrom: str
    start: int
    end: int
    summit: int
    score: float
    condition: int


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    summit: int
    condition1_support: bool
    condition2_support: bool
    origin: str
    condition1_cluster_ids: tuple[str, ...] = ()
    condition2_cluster_ids: tuple[str, ...] = ()
    measurement_intervals: tuple[tuple[int, int], ...] = ()


def _manifest_path(value: str | Path) -> Path:
    path = Path(value).resolve()
    if path.is_dir():
        path = path / MANIFEST_NAME
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def load_stage1_manifest(value: str | Path) -> tuple[Path, dict[str, object]]:
    path = _manifest_path(value)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read Stage 1 manifest: {path}") from exc
    if payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"Not a cutn-suite Stage 1 manifest: {path}")
    if int(payload.get("schema_version", -1)) not in {
        1, 2, 3, 4, 5, MANIFEST_SCHEMA_VERSION
    }:
        raise ValueError(f"Unsupported Stage 1 manifest schema version: {path}")
    return path, payload


def _require_manifest_path(manifest: dict[str, object], field: str) -> Path:
    value = manifest.get(field)
    if not isinstance(value, str) or not value:
        raise ValueError(f"Stage 1 manifest is missing {field}")
    path = Path(value).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def _require_feature_path(
    manifest: dict[str, object], preferred: str, legacy: str
) -> Path:
    """Read the gate-first schema field with a legacy manifest fallback."""

    field = preferred if manifest.get(preferred) else legacy
    return _require_manifest_path(manifest, field)


def _validate_compatibility(
    first: dict[str, object], second: dict[str, object]
) -> dict[str, int]:
    fields = [
        "scoring_method",
        "score_track",
        "positive_track",
        "peak_discovery_track",
        "target_mode",
        "control_mode",
        "contigs",
    ]
    if first.get("target_score_frag_lower") is not None or second.get(
        "target_score_frag_lower"
    ) is not None:
        fields.extend(
            [
                "target_score_frag_lower",
                "target_score_frag_upper",
                "control_score_frag_lower",
                "control_score_frag_upper",
                "coverage_frag_lower",
                "coverage_frag_upper",
            ]
        )
    else:
        fields.extend(["frag_lower", "frag_upper"])
    differences = [field for field in fields if first.get(field) != second.get(field)]
    if differences:
        joined = ", ".join(differences)
        raise ValueError(
            "Stage 1 results are not comparable; incompatible fields: " + joined
        )
    first_track = _require_manifest_path(first, "condition_mean_treatment_coverage")
    second_track = _require_manifest_path(second, "condition_mean_treatment_coverage")
    first_chroms = bigwig_chroms(first_track)
    second_chroms = bigwig_chroms(second_track)
    if first_chroms != second_chroms:
        raise ValueError("Stage 1 coverage BigWigs have different chromosome definitions")
    return first_chroms


def _read_bed_regions(path: Path, *, summit_column: int = 7) -> list[tuple[str, int, int, int]]:
    rows: list[tuple[str, int, int, int]] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED line {line_number} in {path} has fewer than 3 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
                summit = (
                    int(fields[summit_column - 1])
                    if len(fields) >= summit_column
                    else start + (end - start) // 2
                )
            except ValueError as exc:
                raise ValueError(f"Invalid BED coordinates at line {line_number} in {path}") from exc
            if start < 0 or end <= start or not start <= summit < end:
                raise ValueError(f"Invalid BED interval at line {line_number} in {path}")
            rows.append((fields[0], start, end, summit))
    rows.sort(key=lambda value: (value[0], value[3], value[1], value[2]))
    return rows


def _consensus_peak_regions(
    first: Sequence[tuple[str, int, int, int]],
    second: Sequence[tuple[str, int, int, int]],
    *,
    match_distance: int,
    chroms: dict[str, int],
) -> list[Region]:
    candidates: list[tuple[int, int, int, int, str]] = []
    for first_index, left in enumerate(first):
        for second_index, right in enumerate(second):
            if left[0] != right[0]:
                continue
            distance = abs(left[3] - right[3])
            overlaps = max(left[1], right[1]) < min(left[2], right[2])
            if overlaps or distance <= match_distance:
                origin = "overlap_union" if overlaps else "proximity_union"
                candidates.append(
                    (0 if overlaps else 1, distance, first_index, second_index, origin)
                )
    candidates.sort()
    used_first: set[int] = set()
    used_second: set[int] = set()
    matches: list[tuple[int, int, str]] = []
    for _priority, _distance, first_index, second_index, origin in candidates:
        if first_index in used_first or second_index in used_second:
            continue
        used_first.add(first_index)
        used_second.add(second_index)
        matches.append((first_index, second_index, origin))

    seeds: list[tuple[str, int, int, int, bool, bool, str]] = []
    for first_index, second_index, origin in matches:
        left, right = first[first_index], second[second_index]
        seeds.append(
            (
                left[0],
                min(left[1], right[1]),
                max(left[2], right[2]),
                int(round((left[3] + right[3]) / 2.0)),
                True,
                True,
                origin,
            )
        )
    for index, row in enumerate(first):
        if index not in used_first:
            seeds.append(
                (row[0], row[1], row[2], row[3], True, False, "condition1_only")
            )
    for index, row in enumerate(second):
        if index not in used_second:
            seeds.append(
                (row[0], row[1], row[2], row[3], False, True, "condition2_only")
            )

    regions: list[Region] = []
    for chrom, start, end, summit, support1, support2, origin in seeds:
        if chrom not in chroms:
            raise ValueError(f"Peak contig {chrom!r} is absent from the Stage 1 BigWigs")
        start = max(0, start)
        end = min(chroms[chrom], end)
        if end > start:
            regions.append(
                Region(chrom, start, end, summit, support1, support2, origin)
            )
    regions.sort(key=lambda value: (value.chrom, value.start, value.end, value.summit))
    return regions


def _read_cluster_records(path: Path, *, condition: int) -> list[ClusterRecord]:
    records: list[ClusterRecord] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"BED line {line_number} in {path} has fewer than 3 columns"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
                summit = int(fields[6]) if len(fields) >= 7 else start + (end - start) // 2
                score = (
                    float(fields[9])
                    if len(fields) >= 10
                    else float(fields[4])
                    if len(fields) >= 5
                    else 0.0
                )
            except ValueError as exc:
                raise ValueError(
                    f"Invalid cluster BED values at line {line_number} in {path}"
                ) from exc
            if start < 0 or end <= start or not start <= summit < end:
                raise ValueError(f"Invalid BED interval at line {line_number} in {path}")
            records.append(
                ClusterRecord(
                    cluster_id=(
                        fields[3]
                        if len(fields) >= 4 and fields[3]
                        else f"condition{condition}_cluster_{line_number}"
                    ),
                    chrom=fields[0],
                    start=start,
                    end=end,
                    summit=summit,
                    score=score if math.isfinite(score) else 0.0,
                    condition=condition,
                )
            )
    records.sort(key=lambda value: (value.chrom, value.start, value.end, value.summit))
    return records


def _merged_overlap_intervals(
    first: Sequence[ClusterRecord], second: Sequence[ClusterRecord]
) -> tuple[tuple[int, int], ...]:
    overlaps: list[tuple[int, int]] = []
    for left in first:
        for right in second:
            start = max(left.start, right.start)
            end = min(left.end, right.end)
            if start < end:
                overlaps.append((start, end))
    if not overlaps:
        return ()
    overlaps.sort()
    merged: list[list[int]] = []
    for start, end in overlaps:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def _consensus_cluster_regions(
    first: Sequence[ClusterRecord],
    second: Sequence[ClusterRecord],
    *,
    chroms: dict[str, int],
) -> list[Region]:
    """Create connected overlap components while preserving every source cluster."""

    tagged = sorted(
        [*first, *second],
        key=lambda value: (value.chrom, value.start, value.end, value.condition),
    )
    components: list[list[ClusterRecord]] = []
    current: list[ClusterRecord] = []
    current_end = -1
    current_chrom: str | None = None
    for record in tagged:
        if current and record.chrom == current_chrom and record.start < current_end:
            current.append(record)
            current_end = max(current_end, record.end)
        else:
            if current:
                components.append(current)
            current = [record]
            current_chrom = record.chrom
            current_end = record.end
    if current:
        components.append(current)

    output: list[Region] = []
    for component in components:
        chrom = component[0].chrom
        if chrom not in chroms:
            raise ValueError(f"Cluster contig {chrom!r} is absent from Stage 1 BigWigs")
        start = max(0, min(record.start for record in component))
        end = min(chroms[chrom], max(record.end for record in component))
        if end <= start:
            continue
        first_records = [record for record in component if record.condition == 1]
        second_records = [record for record in component if record.condition == 2]
        first_ids = tuple(record.cluster_id for record in first_records)
        second_ids = tuple(record.cluster_id for record in second_records)
        support1, support2 = bool(first_ids), bool(second_ids)
        measurement_intervals = (
            _merged_overlap_intervals(first_records, second_records)
            if support1 and support2
            else ((start, end),)
        )
        origin = (
            "overlap_union"
            if support1 and support2
            else "condition1_only"
            if support1
            else "condition2_only"
        )
        strongest = max(
            component, key=lambda record: (record.score, -record.summit)
        )
        output.append(
            Region(
                chrom=chrom,
                start=start,
                end=end,
                summit=strongest.summit,
                condition1_support=support1,
                condition2_support=support2,
                origin=origin,
                condition1_cluster_ids=first_ids,
                condition2_cluster_ids=second_ids,
                measurement_intervals=measurement_intervals or ((start, end),),
            )
        )
    output.sort(key=lambda value: (value.chrom, value.start, value.end))
    return output


def _merge_intervals(records: Sequence[ClusterRecord]) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for record in sorted(records, key=lambda value: (value.chrom, value.start, value.end)):
        if merged and merged[-1][0] == record.chrom and record.start <= merged[-1][2]:
            chrom, start, end = merged[-1]
            merged[-1] = (chrom, start, max(end, record.end))
        else:
            merged.append((record.chrom, record.start, record.end))
    return merged


def _interval_bp(intervals: Sequence[tuple[str, int, int]]) -> int:
    return sum(end - start for _chrom, start, end in intervals)


def _overlap_bp(
    first: Sequence[tuple[str, int, int]],
    second: Sequence[tuple[str, int, int]],
) -> int:
    by_chrom_first: dict[str, list[tuple[int, int]]] = {}
    by_chrom_second: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in first:
        by_chrom_first.setdefault(chrom, []).append((start, end))
    for chrom, start, end in second:
        by_chrom_second.setdefault(chrom, []).append((start, end))
    total = 0
    for chrom in set(by_chrom_first) & set(by_chrom_second):
        left = by_chrom_first[chrom]
        right = by_chrom_second[chrom]
        i = j = 0
        while i < len(left) and j < len(right):
            total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
            if left[i][1] <= right[j][1]:
                i += 1
            else:
                j += 1
    return total


def _write_cluster_overlap_outputs(
    first: Sequence[ClusterRecord],
    second: Sequence[ClusterRecord],
    regions: Sequence[Region],
    output_dir: Path,
    *,
    condition1_name: str,
    condition2_name: str,
) -> dict[str, object]:
    mapping = output_dir / "cluster_overlap_components.tsv"
    topology = {"1_to_1": 0, "1_to_many": 0, "many_to_1": 0, "many_to_many": 0}
    with mapping.open("wt", encoding="utf-8") as handle:
        handle.write(
            "cluster_locus_id\tchromosome\tstart\tend\tregion_origin\t"
            "condition1_cluster_count\tcondition2_cluster_count\t"
            "condition1_cluster_ids\tcondition2_cluster_ids\trelationship\t"
            "aggregate_anchor_summit\n"
        )
        for index, region in enumerate(regions, 1):
            first_count = len(region.condition1_cluster_ids)
            second_count = len(region.condition2_cluster_ids)
            if first_count and second_count:
                relationship = (
                    "1_to_1" if first_count == second_count == 1
                    else "1_to_many" if first_count == 1
                    else "many_to_1" if second_count == 1
                    else "many_to_many"
                )
                topology[relationship] += 1
            else:
                relationship = region.origin
            handle.write(
                f"cluster_locus_{index}\t{region.chrom}\t{region.start}\t{region.end}\t"
                f"{region.origin}\t{first_count}\t{second_count}\t"
                f"{';'.join(region.condition1_cluster_ids)}\t"
                f"{';'.join(region.condition2_cluster_ids)}\t{relationship}\t"
                f"{region.summit}\n"
            )

    first_intervals = _merge_intervals(first)
    second_intervals = _merge_intervals(second)
    first_bp = _interval_bp(first_intervals)
    second_bp = _interval_bp(second_intervals)
    overlap_bp = _overlap_bp(first_intervals, second_intervals)
    union_bp = first_bp + second_bp - overlap_bp
    first_only_loci = sum(region.origin == "condition1_only" for region in regions)
    second_only_loci = sum(region.origin == "condition2_only" for region in regions)
    shared_loci = sum(region.origin == "overlap_union" for region in regions)
    first_overlapping_ids = {
        cluster_id
        for region in regions
        if region.origin == "overlap_union"
        for cluster_id in region.condition1_cluster_ids
    }
    second_overlapping_ids = {
        cluster_id
        for region in regions
        if region.origin == "overlap_union"
        for cluster_id in region.condition2_cluster_ids
    }
    summary_values: list[tuple[str, object]] = [
        ("condition1_name", condition1_name),
        ("condition2_name", condition2_name),
        ("condition1_raw_cluster_count", len(first)),
        ("condition2_raw_cluster_count", len(second)),
        ("condition1_clusters_with_any_overlap", len(first_overlapping_ids)),
        ("condition2_clusters_with_any_overlap", len(second_overlapping_ids)),
        ("condition1_only_cluster_loci", first_only_loci),
        ("condition2_only_cluster_loci", second_only_loci),
        ("shared_cluster_loci", shared_loci),
        ("union_cluster_loci", len(regions)),
        ("shared_1_to_1_loci", topology["1_to_1"]),
        ("shared_1_to_many_loci", topology["1_to_many"]),
        ("shared_many_to_1_loci", topology["many_to_1"]),
        ("shared_many_to_many_loci", topology["many_to_many"]),
        ("condition1_cluster_bp", first_bp),
        ("condition2_cluster_bp", second_bp),
        ("overlapping_cluster_bp", overlap_bp),
        ("condition1_unique_cluster_bp", first_bp - overlap_bp),
        ("condition2_unique_cluster_bp", second_bp - overlap_bp),
        ("union_cluster_bp", union_bp),
        (
            "condition1_cluster_bp_overlapping_percent",
            100.0 * overlap_bp / first_bp if first_bp else math.nan,
        ),
        (
            "condition2_cluster_bp_overlapping_percent",
            100.0 * overlap_bp / second_bp if second_bp else math.nan,
        ),
        ("cluster_bp_jaccard_percent", 100.0 * overlap_bp / union_bp if union_bp else math.nan),
    ]
    summary = output_dir / "cluster_overlap_summary.tsv"
    with summary.open("wt", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for name, value in summary_values:
            handle.write(f"{name}\t{_format_number(value) if isinstance(value, float) else value}\n")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    venn_path = output_dir / "cluster_locus_venn.png"
    figure, axis = plt.subplots(figsize=(7, 6))
    venn2(
        subsets=(first_only_loci, second_only_loci, shared_loci),
        set_labels=(condition1_name, condition2_name),
        ax=axis,
    )
    axis.set_title("Overlap-connected Stage 1 cluster loci")
    figure.tight_layout()
    figure.savefig(venn_path, dpi=300)
    plt.close(figure)
    return {
        "component_mapping": str(mapping.resolve()),
        "summary": str(summary.resolve()),
        "venn_plot": str(venn_path.resolve()),
        "raw_condition1_clusters": len(first),
        "raw_condition2_clusters": len(second),
        "shared_cluster_loci": shared_loci,
        "condition1_only_cluster_loci": first_only_loci,
        "condition2_only_cluster_loci": second_only_loci,
        "overlapping_cluster_bp": overlap_bp,
        "union_cluster_bp": union_bp,
    }


def _write_union_anchor_bed(regions: Sequence[Region], output_path: Path) -> Path:
    with output_path.open("wt", encoding="utf-8") as handle:
        for index, region in enumerate(regions, 1):
            handle.write(
                f"{region.chrom}\t{region.summit}\t{region.summit + 1}\t"
                f"cluster_locus_{index}\t0\t.\t{region.summit}\t{region.summit + 1}\n"
            )
    return output_path.resolve()


def _group_tracks(manifest: dict[str, object]) -> tuple[list[Path], list[Path]]:
    """Return independent treatment and control raw-coverage tracks."""

    treatment_records = manifest.get("treatment_replicates")
    control_records = manifest.get("control_replicates")
    if isinstance(treatment_records, list) and isinstance(control_records, list):
        if not treatment_records or not control_records:
            raise ValueError("Stage 1 manifest has an empty replicate group")

        def paths(records: list[object], label: str) -> list[Path]:
            output: list[Path] = []
            for record in records:
                if not isinstance(record, dict):
                    raise ValueError(f"Invalid {label} replicate record")
                path = Path(str(record.get("coverage") or record.get("scaled_coverage", ""))).resolve()
                if not path.is_file():
                    raise FileNotFoundError(path)
                output.append(path)
            return output

        return paths(treatment_records, "treatment"), paths(control_records, "control")

    # Schema 1 compatibility: old manifests stored order-paired records.
    records = manifest.get("replicates")
    if not isinstance(records, list) or not records:
        raise ValueError("Stage 1 manifest has no replicate score tracks")
    treatment: list[Path] = []
    control: list[Path] = []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("Invalid replicate record in Stage 1 manifest")
        target_path = Path(str(record.get("treatment_scaled_coverage", ""))).resolve()
        control_path = Path(str(record.get("control_scaled_coverage", ""))).resolve()
        if not target_path.is_file():
            raise FileNotFoundError(target_path)
        if not control_path.is_file():
            raise FileNotFoundError(control_path)
        treatment.append(target_path)
        control.append(control_path)
    return treatment, control


def _bh_qvalues(pvalues: Sequence[float]) -> list[float]:
    qvalues = [math.nan] * len(pvalues)
    finite = [(index, value) for index, value in enumerate(pvalues) if math.isfinite(value)]
    if not finite:
        return qvalues
    ordered = sorted(finite, key=lambda value: value[1])
    running = 1.0
    total = len(ordered)
    for rank in range(total, 0, -1):
        index, pvalue = ordered[rank - 1]
        running = min(running, max(0.0, min(1.0, pvalue * total / rank)))
        qvalues[index] = running
    return qvalues


def _interaction_design(group_sizes: Sequence[int]) -> np.ndarray:
    """Build condition + treatment + condition:treatment factorial design."""

    if len(group_sizes) != 4:
        raise ValueError("Four group sizes are required")
    rows: list[list[float]] = []
    # Input order is condition1 treatment, condition1 control,
    # condition2 treatment, condition2 control.
    for group_index, size in enumerate(group_sizes):
        condition2 = 1.0 if group_index >= 2 else 0.0
        treatment = 1.0 if group_index in {0, 2} else 0.0
        rows.extend(
            [[1.0, condition2, treatment, condition2 * treatment]] * size
        )
    return np.asarray(rows, dtype=float)


def _estimate_variance_prior(
    residual_variances: Sequence[float], residual_df: int
) -> tuple[float, float]:
    """Estimate a scaled-inverse-chi-square variance prior across regions.

    This compact empirical-Bayes estimator follows the moment logic used for
    moderated t statistics.  Log variances are winsorized before fitting so a
    small number of extreme regions cannot determine the shared prior.
    """

    values = np.asarray(residual_variances, dtype=float)
    values = values[np.isfinite(values) & (values >= 0)]
    if residual_df <= 0 or values.size < 2:
        return 0.0, math.nan
    positive = values[values > 0]
    reference = float(np.median(positive)) if positive.size else 1.0
    floor = max(reference * 1e-8, np.finfo(float).tiny)
    logs = np.log(np.maximum(values, floor))
    if logs.size >= 20:
        lower, upper = np.quantile(logs, (0.05, 0.95))
        logs = np.clip(logs, lower, upper)
    observed_variance = float(np.var(logs, ddof=1))
    sampling_variance = float(polygamma(1, residual_df / 2.0))
    target = observed_variance - sampling_variance
    if not math.isfinite(target) or target <= 1e-8:
        prior_df = 1_000_000.0
    else:
        try:
            prior_df = 2.0 * brentq(
                lambda value: float(polygamma(1, value)) - target,
                1e-6,
                1e8,
            )
        except ValueError:
            prior_df = 1_000_000.0 if target < float(polygamma(1, 1e8)) else 0.001
    finite_df = max(prior_df, 1e-6)
    sampling_mean = float(digamma(residual_df / 2.0) - math.log(residual_df / 2.0))
    prior_mean = float(digamma(finite_df / 2.0) - math.log(finite_df / 2.0))
    prior_variance = math.exp(float(np.mean(logs)) - sampling_mean + prior_mean)
    return prior_df, max(prior_variance, floor)


def _moderated_interaction_statistics(
    log_groups_by_region: Sequence[Sequence[Sequence[float]]],
) -> tuple[list[dict[str, float]], dict[str, float | int | bool]]:
    """Fit and moderate the four-group interaction for every region."""

    if not log_groups_by_region:
        return [], {
            "available": False,
            "moderated": False,
            "residual_degrees_freedom": 0,
            "prior_degrees_freedom": math.nan,
            "prior_variance": math.nan,
        }
    group_sizes = [len(group) for group in log_groups_by_region[0]]
    if any(size < 2 for size in group_sizes):
        return [
            {
                "effect": math.nan,
                "ordinary_pvalue": math.nan,
                "moderated_pvalue": math.nan,
                "ordinary_standard_error": math.nan,
                "moderated_standard_error": math.nan,
                "ordinary_ci_lower": math.nan,
                "ordinary_ci_upper": math.nan,
                "moderated_ci_lower": math.nan,
                "moderated_ci_upper": math.nan,
                "residual_variance": math.nan,
                "posterior_variance": math.nan,
            }
            for _ in log_groups_by_region
        ], {
            "available": False,
            "moderated": False,
            "residual_degrees_freedom": 0,
            "prior_degrees_freedom": math.nan,
            "prior_variance": math.nan,
        }
    if any([len(group) for group in groups] != group_sizes for groups in log_groups_by_region):
        raise ValueError("Replicate-group sizes changed between regions")
    design = _interaction_design(group_sizes)
    residual_df = int(design.shape[0] - np.linalg.matrix_rank(design))
    if residual_df <= 0:
        raise ValueError("The four-group interaction model has no residual degrees of freedom")
    inverse = np.linalg.inv(design.T @ design)
    contrast_variance = float(inverse[3, 3])
    fitted: list[dict[str, float]] = []
    residual_variances: list[float] = []
    for groups in log_groups_by_region:
        response = np.asarray([value for group in groups for value in group], dtype=float)
        beta = inverse @ design.T @ response
        residuals = response - design @ beta
        residual_variance = max(float(residuals @ residuals) / residual_df, 0.0)
        effect = float(beta[3])
        if residual_variance == 0:
            ordinary_pvalue = 1.0 if effect == 0 else 0.0
        else:
            statistic = effect / math.sqrt(residual_variance * contrast_variance)
            ordinary_pvalue = 2.0 * float(student_t.sf(abs(statistic), residual_df))
        ordinary_standard_error = math.sqrt(residual_variance * contrast_variance)
        ordinary_critical = float(student_t.ppf(0.975, residual_df))
        fitted.append(
            {
                "effect": effect,
                "ordinary_pvalue": ordinary_pvalue,
                "ordinary_standard_error": ordinary_standard_error,
                "ordinary_ci_lower": effect - ordinary_critical * ordinary_standard_error,
                "ordinary_ci_upper": effect + ordinary_critical * ordinary_standard_error,
                "residual_variance": residual_variance,
            }
        )
        residual_variances.append(residual_variance)
    prior_df, prior_variance = _estimate_variance_prior(residual_variances, residual_df)
    moderated = prior_df > 0 and math.isfinite(prior_variance)
    for record in fitted:
        if moderated:
            posterior_variance = (
                prior_df * prior_variance
                + residual_df * record["residual_variance"]
            ) / (prior_df + residual_df)
            total_df = residual_df + prior_df
            if posterior_variance == 0:
                pvalue = 1.0 if record["effect"] == 0 else 0.0
            else:
                statistic = record["effect"] / math.sqrt(
                    posterior_variance * contrast_variance
                )
                pvalue = 2.0 * float(student_t.sf(abs(statistic), total_df))
        else:
            posterior_variance = record["residual_variance"]
            pvalue = record["ordinary_pvalue"]
            total_df = residual_df
        moderated_standard_error = math.sqrt(
            posterior_variance * contrast_variance
        )
        moderated_critical = float(student_t.ppf(0.975, total_df))
        record["posterior_variance"] = posterior_variance
        record["moderated_standard_error"] = moderated_standard_error
        record["moderated_ci_lower"] = (
            record["effect"] - moderated_critical * moderated_standard_error
        )
        record["moderated_ci_upper"] = (
            record["effect"] + moderated_critical * moderated_standard_error
        )
        record["moderated_pvalue"] = min(max(pvalue, 0.0), 1.0)
    return fitted, {
        "available": True,
        "moderated": moderated,
        "residual_degrees_freedom": residual_df,
        "prior_degrees_freedom": prior_df,
        "prior_variance": prior_variance,
    }


def _format_number(value: float) -> str:
    return "" if not math.isfinite(value) else f"{value:.12g}"


def _compare_regions(
    regions: Sequence[Region],
    first_manifest: dict[str, object],
    second_manifest: dict[str, object],
    *,
    statistic_name: str,
    statistic: Callable[[object, str, int, int], float],
    output_path: Path,
    fdr: float,
) -> dict[str, object]:
    first_treatment, first_control = _group_tracks(first_manifest)
    second_treatment, second_control = _group_tracks(second_manifest)
    groups = [first_treatment, first_control, second_treatment, second_control]
    all_paths = [path for group in groups for path in group]
    handles = open_bigwigs(all_paths)
    grouped_handles: list[list[object]] = []
    offset = 0
    for group in groups:
        grouped_handles.append(handles[offset : offset + len(group)])
        offset += len(group)
    first_treatment_handles, first_control_handles, second_treatment_handles, second_control_handles = grouped_handles
    rows: list[dict[str, object]] = []
    try:
        for region in regions:
            intervals = region.measurement_intervals or ((region.start, region.end),)
            total_bp = sum(end - start for start, end in intervals)
            def measure_handle(handle: object) -> float:
                if total_bp <= 0:
                    return 0.0
                weighted = sum(
                    float(statistic(handle, region.chrom, start, end)) * (end - start)
                    for start, end in intervals
                )
                return max(weighted / total_bp, 0.0)
            measured = [[measure_handle(handle) for handle in group] for group in grouped_handles]
            log_measured = [
                [math.log2(value + 1.0) for value in group] for group in measured
            ]
            (
                first_treatment_values,
                first_control_values,
                second_treatment_values,
                second_control_values,
            ) = measured
            (
                first_treatment_log,
                first_control_log,
                second_treatment_log,
                second_control_log,
            ) = log_measured
            first_treatment_mean = float(np.mean(first_treatment_values))
            first_control_mean = float(np.mean(first_control_values))
            second_treatment_mean = float(np.mean(second_treatment_values))
            second_control_mean = float(np.mean(second_control_values))
            first_enrichment = first_treatment_mean - first_control_mean
            second_enrichment = second_treatment_mean - second_control_mean
            raw_delta = second_enrichment - first_enrichment
            log_means = tuple(float(np.mean(group)) for group in log_measured)
            first_log_enrichment = log_means[0] - log_means[1]
            second_log_enrichment = log_means[2] - log_means[3]
            log_delta = second_log_enrichment - first_log_enrichment
            first_lower = min(first_treatment_log) - max(first_control_log)
            first_upper = max(first_treatment_log) - min(first_control_log)
            second_lower = min(second_treatment_log) - max(second_control_log)
            second_upper = max(second_treatment_log) - min(second_control_log)
            if second_lower > first_upper:
                consistency = "robust_gain"
            elif second_upper < first_lower:
                consistency = "robust_loss"
            else:
                consistency = "not_replicate_separated"
            rows.append(
                {
                    "region": region,
                    "groups": measured,
                    "log_groups": log_measured,
                    "means": (
                        first_treatment_mean,
                        first_control_mean,
                        second_treatment_mean,
                        second_control_mean,
                    ),
                    "log_means": log_means,
                    "first_enrichment": first_enrichment,
                    "second_enrichment": second_enrichment,
                    "raw_delta": raw_delta,
                    "first_log_enrichment": first_log_enrichment,
                    "second_log_enrichment": second_log_enrichment,
                    "log_delta": log_delta,
                    "first_bounds": (first_lower, first_upper),
                    "second_bounds": (second_lower, second_upper),
                    "replicate_consistency": consistency,
                }
            )
    finally:
        for handle in handles:
            handle.close()

    inferential_requested = (
        first_manifest.get("bam_mode") == "replicates"
        and second_manifest.get("bam_mode") == "replicates"
        and all(len(group) >= 2 for group in groups)
    )
    if inferential_requested:
        model_rows, model_metadata = _moderated_interaction_statistics(
            [row["log_groups"] for row in rows]  # type: ignore[list-item]
        )
    else:
        model_rows = [
            {
                "effect": float(row["log_delta"]),
                "ordinary_pvalue": math.nan,
                "moderated_pvalue": math.nan,
                "ordinary_standard_error": math.nan,
                "moderated_standard_error": math.nan,
                "ordinary_ci_lower": math.nan,
                "ordinary_ci_upper": math.nan,
                "moderated_ci_lower": math.nan,
                "moderated_ci_upper": math.nan,
                "residual_variance": math.nan,
                "posterior_variance": math.nan,
            }
            for row in rows
        ]
        model_metadata = {
            "available": False,
            "moderated": False,
            "residual_degrees_freedom": 0,
            "prior_degrees_freedom": math.nan,
            "prior_variance": math.nan,
        }
    for row, model in zip(rows, model_rows):
        row.pop("log_groups", None)
        row["model"] = model
    qvalues = _bh_qvalues(
        [float(model["moderated_pvalue"]) for model in model_rows]
    )
    inferential = bool(model_metadata["available"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gain_path = output_path.with_name(output_path.stem + f"_fdr{fdr:g}_gains.bed")
    loss_path = output_path.with_name(output_path.stem + f"_fdr{fdr:g}_losses.bed")
    all_gain_path = output_path.with_name(output_path.stem + "_all_gains.bed")
    all_loss_path = output_path.with_name(output_path.stem + "_all_losses.bed")
    robust_gain_path = output_path.with_name(output_path.stem + "_robust_gains.bed")
    robust_loss_path = output_path.with_name(output_path.stem + "_robust_losses.bed")
    with (
        output_path.open("wt", encoding="utf-8") as table,
        gain_path.open("wt", encoding="utf-8") as gains,
        loss_path.open("wt", encoding="utf-8") as losses,
        all_gain_path.open("wt", encoding="utf-8") as all_gains,
        all_loss_path.open("wt", encoding="utf-8") as all_losses,
        robust_gain_path.open("wt", encoding="utf-8") as robust_gains,
        robust_loss_path.open("wt", encoding="utf-8") as robust_losses,
    ):
        replicate_columns = [
            *[f"condition1_treatment_replicate_{index}" for index in range(1, len(first_treatment) + 1)],
            *[f"condition1_control_replicate_{index}" for index in range(1, len(first_control) + 1)],
            *[f"condition2_treatment_replicate_{index}" for index in range(1, len(second_treatment) + 1)],
            *[f"condition2_control_replicate_{index}" for index in range(1, len(second_control) + 1)],
        ]
        table.write(
            "chromosome\tstart\tend\tsummit\tcondition1_stage1_support\t"
            "condition2_stage1_support\tregion_origin\tmeasurement_intervals\tmeasurement_bp\t"
            "condition1_cluster_count\tcondition2_cluster_count\t"
            "condition1_cluster_ids\tcondition2_cluster_ids\tstatistic\t"
            + "\t".join(replicate_columns)
            + "\tcondition1_treatment_mean\tcondition1_control_mean\t"
            "condition2_treatment_mean\tcondition2_control_mean\t"
            "condition1_mean_enrichment\tcondition2_mean_enrichment\t"
            "raw_interaction_difference\tcondition1_log2_treatment_mean\t"
            "condition1_log2_control_mean\tcondition2_log2_treatment_mean\t"
            "condition2_log2_control_mean\tcondition1_log2_enrichment\t"
            "condition2_log2_enrichment\tlog2_interaction_difference\t"
            "condition1_log2_enrichment_lower\tcondition1_log2_enrichment_upper\t"
            "condition2_log2_enrichment_lower\tcondition2_log2_enrichment_upper\t"
            "effect_direction\treplicate_consistency\tordinary_p_value\t"
            "moderated_p_value\tdifferential_fdr\tordinary_standard_error\t"
            "moderated_standard_error\tordinary_95_ci_lower\tordinary_95_ci_upper\t"
            "moderated_95_ci_lower\tmoderated_95_ci_upper\tresidual_variance\t"
            "posterior_variance\tstatus\n"
        )
        for index, (row, qvalue) in enumerate(zip(rows, qvalues), 1):
            region = row["region"]
            delta = float(row["log_delta"])
            direction = "gain" if delta > 0 else "loss" if delta < 0 else "stable"
            if inferential and math.isfinite(qvalue) and qvalue <= fdr and delta > 0:
                status = "significant_gain"
            elif inferential and math.isfinite(qvalue) and qvalue <= fdr and delta < 0:
                status = "significant_loss"
            elif inferential:
                status = f"not_significant_{direction}"
            elif delta > 0:
                status = "descriptive_gain"
            elif delta < 0:
                status = "descriptive_loss"
            else:
                status = "stable"
            values = [
                region.chrom,
                str(region.start),
                str(region.end),
                str(region.summit),
                str(region.condition1_support).lower(),
                str(region.condition2_support).lower(),
                region.origin,
                ";".join(f"{start}-{end}" for start, end in (region.measurement_intervals or ((region.start, region.end),))),
                str(sum(end - start for start, end in (region.measurement_intervals or ((region.start, region.end),)))),
                str(len(region.condition1_cluster_ids)),
                str(len(region.condition2_cluster_ids)),
                ";".join(region.condition1_cluster_ids),
                ";".join(region.condition2_cluster_ids),
                statistic_name,
                *(
                    _format_number(float(value))
                    for group in row["groups"]
                    for value in group
                ),
                *(_format_number(float(value)) for value in row["means"]),
                _format_number(float(row["first_enrichment"])),
                _format_number(float(row["second_enrichment"])),
                _format_number(float(row["raw_delta"])),
                *(_format_number(float(value)) for value in row["log_means"]),
                _format_number(float(row["first_log_enrichment"])),
                _format_number(float(row["second_log_enrichment"])),
                _format_number(delta),
                *(_format_number(float(value)) for value in row["first_bounds"]),
                *(_format_number(float(value)) for value in row["second_bounds"]),
                direction,
                str(row["replicate_consistency"]),
                _format_number(float(row["model"]["ordinary_pvalue"])),
                _format_number(float(row["model"]["moderated_pvalue"])),
                _format_number(qvalue),
                _format_number(float(row["model"]["ordinary_standard_error"])),
                _format_number(float(row["model"]["moderated_standard_error"])),
                _format_number(float(row["model"]["ordinary_ci_lower"])),
                _format_number(float(row["model"]["ordinary_ci_upper"])),
                _format_number(float(row["model"]["moderated_ci_lower"])),
                _format_number(float(row["model"]["moderated_ci_upper"])),
                _format_number(float(row["model"]["residual_variance"])),
                _format_number(float(row["model"]["posterior_variance"])),
                status,
            ]
            table.write("\t".join(values) + "\n")
            bed_line = (
                f"{region.chrom}\t{region.start}\t{region.end}\t"
                f"{statistic_name}_{index}\t{abs(delta):.6f}\t.\t"
                f"{region.summit}\t{region.summit + 1}\t{_format_number(qvalue)}\t"
                f"{region.origin}\n"
            )
            if direction == "gain":
                all_gains.write(bed_line)
            elif direction == "loss":
                all_losses.write(bed_line)
            if status == "significant_gain":
                gains.write(bed_line)
            elif status == "significant_loss":
                losses.write(bed_line)
            if row["replicate_consistency"] == "robust_gain":
                robust_gains.write(bed_line)
            elif row["replicate_consistency"] == "robust_loss":
                robust_losses.write(bed_line)
    origin_counts: dict[str, int] = {}
    for region in regions:
        origin_counts[region.origin] = origin_counts.get(region.origin, 0) + 1
    return {
        "table": str(output_path),
        "significant_gains": str(gain_path),
        "significant_losses": str(loss_path),
        "all_gains": str(all_gain_path),
        "all_losses": str(all_loss_path),
        "robust_gains": str(robust_gain_path),
        "robust_losses": str(robust_loss_path),
        "regions": len(regions),
        "region_origin_counts": origin_counts,
        "inferential_fdr_available": inferential,
        "statistical_model": "empirical_bayes_log2_factorial_interaction",
        "moderation": model_metadata,
    }


def _manifest_score_tracks(
    manifest: dict[str, object],
) -> tuple[Path, list[Path], str, str] | None:
    """Return native treatment PNS tracks used for cluster aggregates."""

    method = str(manifest.get("scoring_method") or "pns").lower()
    if method != "pns":
        return None
    positive_track = "posPNS"
    mean_value = manifest.get("condition_mean_treatment_cluster_aggregate_score")
    records = manifest.get("treatment_replicates")
    replicate_field = "analysis_score"

    if not isinstance(mean_value, str) or not mean_value or not isinstance(records, list):
        return None
    mean_path = Path(mean_value).resolve()
    replicate_paths: list[Path] = []
    for record in records:
        if not isinstance(record, dict) or not record.get(replicate_field):
            return None
        replicate_paths.append(Path(str(record[replicate_field])).resolve())
    if not mean_path.is_file() or any(not path.is_file() for path in replicate_paths):
        return None
    return mean_path, replicate_paths, method, positive_track


def _run_shared_cluster_aggregates(
    first: dict[str, object],
    second: dict[str, object],
    regions: Sequence[Region],
    output_dir: Path,
) -> dict[str, object]:
    first_tracks = _manifest_score_tracks(first)
    second_tracks = _manifest_score_tracks(second)
    if first_tracks is None or second_tracks is None:
        return {
            "status": "unavailable",
            "reason": (
                "Stage 1 manifest lacks native replicate PNS tracks; "
                "rerun Stage 1 with NucleoSuite 0.10.11 or later."
            ),
        }
    aggregate_dir = output_dir / "cluster_aligned_aggregates"
    aggregate_dir.mkdir(parents=True, exist_ok=True)
    anchors = _write_union_anchor_bed(
        regions, aggregate_dir / "union_cluster_locus_anchors.bed"
    )
    parameters = first.get("cluster_aggregate_parameters")
    parameters = parameters if isinstance(parameters, dict) else {}
    limit = common_symmetric_bigwig_limit([first_tracks[0], second_tracks[0]])
    common = dict(
        anchor_bed=anchors,
        window_half=int(parameters.get("window_half", 1000)),
        maximum_heatmap_rows=int(parameters.get("maximum_heatmap_rows", 5000)),
        bootstrap_replicates=int(parameters.get("bootstrap_replicates", 200)),
        nrl_peak_resolution=float(parameters.get("nrl_peak_resolution", 130.0)),
        nrl_min_order=int(parameters.get("nrl_min_order", 0)),
        nrl_max_order=int(parameters.get("nrl_max_order", 3)),
        vlim=limit,
    )
    first_name = str(first.get("condition_name") or "condition1")
    second_name = str(second.get("condition_name") or "condition2")
    first_outputs = run_cluster_aggregate(
        mean_score=first_tracks[0],
        replicate_scores=first_tracks[1],
        output_dir=aggregate_dir / "condition1",
        label=first_name,
        scoring_method=first_tracks[2],
        positive_track=first_tracks[3],
        seed=12345,
        **common,
    )
    second_outputs = run_cluster_aggregate(
        mean_score=second_tracks[0],
        replicate_scores=second_tracks[1],
        output_dir=aggregate_dir / "condition2",
        label=second_name,
        scoring_method=second_tracks[2],
        positive_track=second_tracks[3],
        seed=12345,
        **common,
    )
    return {
        "status": "complete",
        "common_anchor_bed": str(anchors),
        "common_symmetric_heatmap_limit": limit,
        "condition1": first_outputs,
        "condition2": second_outputs,
    }


def compare_stage1(
    condition1_results: str | Path,
    condition2_results: str | Path,
    *,
    outdir: str | Path,
    fdr: float = 0.05,
) -> Path:
    if not math.isfinite(fdr) or not 0 <= fdr <= 1:
        raise ValueError("fdr must be between 0 and 1")
    first_path, first = load_stage1_manifest(condition1_results)
    second_path, second = load_stage1_manifest(condition2_results)
    chroms = _validate_compatibility(first, second)
    output_dir = Path(outdir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    first_clusters = _read_cluster_records(
        _require_feature_path(first, "selected_clusters", "significant_clusters"),
        condition=1,
    )
    second_clusters = _read_cluster_records(
        _require_feature_path(second, "selected_clusters", "significant_clusters"),
        condition=2,
    )
    cluster_regions = _consensus_cluster_regions(
        first_clusters, second_clusters, chroms=chroms
    )
    overlap = _write_cluster_overlap_outputs(
        first_clusters,
        second_clusters,
        cluster_regions,
        output_dir,
        condition1_name=str(first.get("condition_name") or "condition1"),
        condition2_name=str(second.get("condition_name") or "condition2"),
    )
    results: dict[str, object] = {
        "clusters": _compare_regions(
            cluster_regions,
            first,
            second,
            statistic_name="mean_raw_coverage_over_overlap",
            statistic=interval_mean,
            output_path=output_dir / "differential_clusters.tsv",
            fdr=fdr,
        )
    }
    aggregate_outputs = _run_shared_cluster_aggregates(
        first, second, cluster_regions, output_dir
    )

    manifest_path = output_dir / "cutn_comparison_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "nucleosuite_cutn_comparison",
                "schema_version": 4,
                "condition1_manifest": str(first_path),
                "condition2_manifest": str(second_path),
                "condition1_name": first.get("condition_name"),
                "condition2_name": second.get("condition_name"),
                "feature_level": "clusters",
                "peak_discovery_track": first.get("peak_discovery_track"),
                "measurement_track": first.get("peak_measurement_track"),
                "statistical_model": "empirical_bayes_log2_factorial_interaction",
                "coverage_transform": "log2(mean_raw_coverage_plus_1)",
                "fdr": fdr,
                "cluster_overlap": overlap,
                "cluster_aligned_aggregates": aggregate_outputs,
                "genomic_randomization_overlap_test": "not_performed",
                "results": results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite cutn-compare",
        description=(
            "Compare two completed cutn-suite Stage 1 analyses using their "
            "coverage BigWigs scaled to a non-zero mean of 100 and a log2 "
            "empirical-Bayes interaction model at overlap-connected cluster "
            "loci; BAM files are not revisited."
        ),
    )
    parser.add_argument(
        "--condition1-results",
        required=True,
        help=f"Condition 1 Stage 1 directory or {MANIFEST_NAME}.",
    )
    parser.add_argument(
        "--condition2-results",
        required=True,
        help=f"Condition 2 Stage 1 directory or {MANIFEST_NAME}.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument(
        "--fdr", type=float, default=0.05,
        help=(
            "Moderated differential FDR cutoff for the separate significant "
            "gain/loss BEDs; all regions remain in the complete tables (default: 0.05)."
        ),
    )
    return parser


def run(args: argparse.Namespace) -> int:
    output = compare_stage1(
        args.condition1_results,
        args.condition2_results,
        outdir=args.outdir,
        fdr=args.fdr,
    )
    print(f"cutn_comparison_manifest\t{output}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
