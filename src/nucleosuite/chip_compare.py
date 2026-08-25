#!/usr/bin/env python3
"""Differential comparison of two completed chip-suite Stage 1 analyses."""

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
    interval_max,
    interval_positive_area,
    open_bigwigs,
)


MANIFEST_NAME = "chip_stage1_manifest.json"
MANIFEST_SCHEMA = "nucleosuite_chip_stage1"
MANIFEST_SCHEMA_VERSION = 3


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    summit: int
    condition1_support: bool
    condition2_support: bool
    origin: str


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
        raise ValueError(f"Not a chip-suite Stage 1 manifest: {path}")
    if int(payload.get("schema_version", -1)) not in {1, 2, MANIFEST_SCHEMA_VERSION}:
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
    fields = (
        "scoring_method",
        "score_track",
        "positive_track",
        "peak_discovery_track",
        "peak_measurement_track",
        "target_mode",
        "control_mode",
        "frag_lower",
        "frag_upper",
        "contigs",
        "stage1_statistics",
    )
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
        raise ValueError("Stage 1 scaled-coverage BigWigs have different chromosome definitions")
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


def _consensus_cluster_regions(
    first: Sequence[tuple[str, int, int, int]],
    second: Sequence[tuple[str, int, int, int]],
) -> list[Region]:
    def origin(support1: bool, support2: bool) -> str:
        if support1 and support2:
            return "overlap_union"
        return "condition1_only" if support1 else "condition2_only"

    tagged = [(*row[:3], True, False) for row in first]
    tagged.extend((*row[:3], False, True) for row in second)
    tagged.sort(key=lambda value: (value[0], value[1], value[2]))
    output: list[Region] = []
    for chrom, start, end, support1, support2 in tagged:
        if output and output[-1].chrom == chrom and start <= output[-1].end:
            previous = output.pop()
            new_start, new_end = previous.start, max(previous.end, end)
            combined_support1 = previous.condition1_support or support1
            combined_support2 = previous.condition2_support or support2
            output.append(
                Region(
                    chrom,
                    new_start,
                    new_end,
                    new_start + (new_end - new_start) // 2,
                    combined_support1,
                    combined_support2,
                    origin(combined_support1, combined_support2),
                )
            )
        else:
            output.append(
                Region(
                    chrom,
                    start,
                    end,
                    start + (end - start) // 2,
                    support1,
                    support2,
                    origin(support1, support2),
                )
            )
    return output


def _group_tracks(manifest: dict[str, object]) -> tuple[list[Path], list[Path]]:
    """Return independent treatment and control scaled-coverage tracks."""

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
                path = Path(str(record.get("scaled_coverage", ""))).resolve()
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
        fitted.append(
            {
                "effect": effect,
                "ordinary_pvalue": ordinary_pvalue,
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
        record["posterior_variance"] = posterior_variance
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
            measured = [
                [
                    max(
                        float(statistic(handle, region.chrom, region.start, region.end)),
                        0.0,
                    )
                    for handle in group
                ]
                for group in grouped_handles
            ]
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
            "condition2_stage1_support\tregion_origin\tstatistic\t"
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
            "moderated_p_value\tdifferential_fdr\tresidual_variance\t"
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


def compare_stage1(
    condition1_results: str | Path,
    condition2_results: str | Path,
    *,
    outdir: str | Path,
    fdr: float = 0.05,
    feature_level: str = "both",
    peak_match_distance: int | None = None,
) -> Path:
    if not math.isfinite(fdr) or not 0 <= fdr <= 1:
        raise ValueError("fdr must be between 0 and 1")
    if feature_level not in {"peaks", "clusters", "both"}:
        raise ValueError("feature_level must be peaks, clusters or both")
    first_path, first = load_stage1_manifest(condition1_results)
    second_path, second = load_stage1_manifest(condition2_results)
    chroms = _validate_compatibility(first, second)
    output_dir = Path(outdir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, object] = {}
    if feature_level in {"peaks", "both"}:
        first_peaks = _read_bed_regions(
            _require_feature_path(first, "selected_peaks", "significant_peaks")
        )
        second_peaks = _read_bed_regions(
            _require_feature_path(second, "selected_peaks", "significant_peaks")
        )
        distance = (
            int(peak_match_distance)
            if peak_match_distance is not None
            else 0
        )
        peak_regions = _consensus_peak_regions(
            first_peaks,
            second_peaks,
            match_distance=distance,
            chroms=chroms,
        )
        results["peaks"] = _compare_regions(
            peak_regions,
            first,
            second,
            statistic_name="peak_max",
            statistic=interval_max,
            output_path=output_dir / "differential_peaks.tsv",
            fdr=fdr,
        )

    if feature_level in {"clusters", "both"}:
        first_clusters = _read_bed_regions(
            _require_feature_path(first, "selected_clusters", "significant_clusters")
        )
        second_clusters = _read_bed_regions(
            _require_feature_path(second, "selected_clusters", "significant_clusters")
        )
        cluster_regions = _consensus_cluster_regions(first_clusters, second_clusters)
        results["clusters"] = _compare_regions(
            cluster_regions,
            first,
            second,
            statistic_name="cluster_positive_area",
            statistic=interval_positive_area,
            output_path=output_dir / "differential_clusters.tsv",
            fdr=fdr,
        )

    manifest_path = output_dir / "chip_comparison_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "nucleosuite_chip_comparison",
                "schema_version": 3,
                "condition1_manifest": str(first_path),
                "condition2_manifest": str(second_path),
                "condition1_name": first.get("condition_name"),
                "condition2_name": second.get("condition_name"),
                "feature_level": feature_level,
                "peak_discovery_track": first.get("peak_discovery_track"),
                "measurement_track": first.get("peak_measurement_track"),
                "statistical_model": "empirical_bayes_log2_factorial_interaction",
                "coverage_transform": "log2(mean_scaled_coverage_plus_1)",
                "fdr": fdr,
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
        prog="nucleosuite chip-compare",
        description=(
            "Compare two completed chip-suite Stage 1 analyses using their "
            "coverage BigWigs scaled to a non-zero mean of 100 and a log2 "
            "empirical-Bayes interaction model; BAM files are not revisited."
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
        "--feature-level",
        choices=("peaks", "clusters", "both"),
        default="both",
        help="Compare Stage 1 peaks, clusters, or both (default: both).",
    )
    parser.add_argument(
        "--peak-match-distance",
        type=int,
        help=(
            "Also merge non-overlapping peaks whose summits are within this distance; "
            "by default only overlapping peaks are merged."
        ),
    )
    parser.add_argument(
        "--fdr", type=float, default=0.05,
        help=(
            "Moderated differential FDR cutoff for the separate significant "
            "gain/loss BEDs; all regions remain in the complete tables (default: 0.05)."
        ),
    )
    return parser


def run(args: argparse.Namespace) -> int:
    if args.peak_match_distance is not None and args.peak_match_distance < 0:
        raise ValueError("--peak-match-distance must be non-negative")
    output = compare_stage1(
        args.condition1_results,
        args.condition2_results,
        outdir=args.outdir,
        fdr=args.fdr,
        feature_level=args.feature_level,
        peak_match_distance=args.peak_match_distance,
    )
    print(f"chip_comparison_manifest\t{output}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
