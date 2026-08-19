#!/usr/bin/env python3
"""Relate nucleosome peak spacing and periodic signal to gene expression.

The command implements three complementary analyses for any continuous signal
track. PNS is the default signal label in NucleoSuite.

* Peak-spacing analysis: median adjacent peak-centre distance in each gene body
  and 10 kb flanks, correlated with expression profiles.
* Per-gene FFT analysis: periodograms from a strand-aware fixed window beginning
  at each transcription start site, with expression-correlation trajectories.
* Profile ranking: correlations between expression and the mean FFT intensity at
  selected periods (default 193, 196, and 199 bp), ranked from most negative.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nucleosuite.plotting import configure_unique_category_cycle
configure_unique_category_cycle()
import numpy as np
from scipy import fft as scipy_fft
from scipy import signal as scipy_signal
from scipy import stats as scipy_stats

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.io import open_text
from nucleosuite.parallel import add_parallel_arguments
from nucleosuite.partitioned import run_partitioned_command
from nucleosuite.resource_files import (
    materialized_resource_path,
    resolve_set_resource_name,
)
from nucleosuite.core.blacklist import BlacklistIndex, load_blacklist_unbounded
from nucleosuite.progress import ProgressReporter


HEADER_PREFIXES = ("#", "track", "browser")
ENSEMBL_VERSION_RE = re.compile(r"^(ENS[A-Z]*G\d+)(?:\.\d+)?$", re.IGNORECASE)


@dataclass(frozen=True)
class GeneRecord:
    chrom: str
    start: int
    end: int
    gene_id: str
    gene_name: str
    strand: str

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end - 1


@dataclass
class ExpressionData:
    values: dict[str, dict[str, float]]
    gene_names: dict[str, str]
    profiles: list[str]
    profile_types: dict[str, str]
    genes: set[str]
    duplicate_records: int
    invalid_records: int
    filtered_gene_count: int


@dataclass(frozen=True)
class NamedPath:
    name: str
    path: Path


@dataclass
class CorrelationResult:
    r: float
    p: float
    n: int


@dataclass
class FFTRow:
    sample: str
    gene: GeneRecord
    valid_fraction: float
    dominant_period: int
    dominant_intensity: float
    ranking_intensity: float
    intensities: np.ndarray


def normalize_gene_id(value: str) -> str:
    text = value.strip()
    match = ENSEMBL_VERSION_RE.match(text)
    return match.group(1).upper() if match else text



def open_expression_table(path: str | Path):
    """Open plain or gzip-compressed expression tables with BOM-safe decoding."""
    input_path = Path(path)
    if input_path.suffix.casefold() == ".gz":
        return gzip.open(input_path, "rt", newline="", encoding="utf-8-sig")
    return input_path.open("rt", newline="", encoding="utf-8-sig")

def split_fields(line: str) -> list[str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) == 1:
        fields = line.split()
    return fields


def _column_index(value: int, option: str) -> int:
    if value < 1:
        raise ValueError(f"{option} must be a positive one-based column")
    return value - 1



def _normalise_profile_metadata_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def _normalise_profile_type(value: str) -> str:
    text = re.sub(r"[^a-z0-9]+", "_", value.strip().casefold()).strip("_")
    if text in {"cell_line", "cellline", "cell_lines", "celllines"}:
        return "cell_line"
    if text in {"tissue", "tissues"}:
        return "tissue"
    return text


def read_profile_metadata(path: str | Path) -> dict[str, dict[str, str]]:
    """Read cell-line metadata using exact and conservative normalized names."""
    exact: dict[str, dict[str, str]] = {}
    normalized_candidates: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Cell line", "Cancer cell line", "Primary disease", "Disease subtype",
            "Primary/Metastasis", "Sample collection site", "Cellosaurus ID",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            available = ", ".join(reader.fieldnames or [])
            raise ValueError(f"Cell-line metadata has incompatible columns: {available}")
        for row in reader:
            name = (row.get("Cell line") or "").strip()
            if not name:
                continue
            record = {
                "cancer_cell_line": (row.get("Cancer cell line") or "").strip(),
                "primary_disease": (row.get("Primary disease") or "").strip(),
                "disease_subtype": (row.get("Disease subtype") or "").strip(),
                "primary_or_metastasis": (row.get("Primary/Metastasis") or "").strip(),
                "sample_collection_site": (row.get("Sample collection site") or "").strip(),
                "cellosaurus_id": (row.get("Cellosaurus ID") or "").strip(),
            }
            exact[name.casefold()] = record
            normalized_candidates[_normalise_profile_metadata_name(name)].append(record)
    result = dict(exact)
    for key, records in normalized_candidates.items():
        unique = {tuple(sorted(record.items())) for record in records}
        if len(unique) == 1:
            result[f"normalized:{key}"] = records[0]
    return result


def annotate_profile(
    profile: str,
    metadata: Mapping[str, Mapping[str, str]] | None,
    *,
    profile_type: str = "",
) -> dict[str, str]:
    empty = {
        "profile_type": profile_type,
        "cancer_cell_line": "",
        "primary_disease": "",
        "disease_subtype": "",
        "primary_or_metastasis": "",
        "sample_collection_site": "",
        "cellosaurus_id": "",
        "metadata_matched": "false",
        "metadata_status": "not_found",
    }
    if profile_type == "tissue":
        return {**empty, "metadata_status": "not_applicable"}
    if not metadata:
        return {**empty, "metadata_status": "unavailable"}
    record = metadata.get(profile.casefold())
    if record is None:
        record = metadata.get(f"normalized:{_normalise_profile_metadata_name(profile)}")
    if record is None:
        return empty
    return {
        **empty,
        **dict(record),
        "metadata_matched": "true",
        "metadata_status": "matched",
    }

def read_genes(
    path: str | Path,
    *,
    gene_id_column: int = 4,
    gene_name_column: int = 5,
    strand_column: int = 6,
) -> list[GeneRecord]:
    id_idx = _column_index(gene_id_column, "--gene-id-column")
    name_idx = _column_index(gene_name_column, "--gene-name-column")
    strand_idx = _column_index(strand_column, "--gene-strand-column")
    required = max(2, id_idx, name_idx, strand_idx)
    genes: list[GeneRecord] = []
    seen: set[str] = set()

    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(HEADER_PREFIXES):
                continue
            fields = split_fields(raw)
            if len(fields) <= required:
                raise ValueError(
                    f"{path}:{line_no}: expected at least {required + 1} columns"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: gene coordinates must be integers") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: require 0 <= start < end")
            strand = fields[strand_idx]
            if strand not in {"+", "-"}:
                raise ValueError(f"{path}:{line_no}: strand must be '+' or '-'")
            gene_id = normalize_gene_id(fields[id_idx])
            if not gene_id:
                raise ValueError(f"{path}:{line_no}: empty gene identifier")
            if gene_id in seen:
                raise ValueError(f"{path}:{line_no}: duplicate gene identifier {gene_id!r}")
            seen.add(gene_id)
            genes.append(
                GeneRecord(
                    chrom=fields[0],
                    start=start,
                    end=end,
                    gene_id=gene_id,
                    gene_name=fields[name_idx] or gene_id,
                    strand=strand,
                )
            )
    if not genes:
        raise ValueError(f"No valid genes were found in {path}")
    genes.sort(key=lambda g: (g.chrom, g.start, g.end, g.gene_id))
    return genes


def filter_blacklisted_gene_anchors(
    genes: Sequence[GeneRecord], blacklist: BlacklistIndex | None
) -> tuple[list[GeneRecord], int]:
    """Exclude only genes whose one-base strand-aware TSS anchor is blacklisted."""
    if blacklist is None:
        return list(genes), 0
    retained = [
        gene
        for gene in genes
        if not blacklist.overlaps(gene.chrom, gene.tss, gene.tss + 1)
    ]
    return retained, len(genes) - len(retained)


def read_expression(
    path: str | Path,
    *,
    gene_column: str,
    name_column: str,
    profile_column: str,
    profile_type_column: str | None,
    value_column: str,
    min_nonzero_profiles: int,
) -> ExpressionData:
    sums: dict[tuple[str, str], float] = defaultdict(float)
    counts: dict[tuple[str, str], int] = defaultdict(int)
    gene_names: dict[str, str] = {}
    profiles: set[str] = set()
    profile_type_candidates: dict[str, set[str]] = defaultdict(set)
    invalid = 0
    duplicates = 0

    with open_expression_table(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Expression file has no header: {path}")
        missing = [
            name
            for name in (gene_column, name_column, profile_column, value_column)
            if name not in reader.fieldnames
        ]
        if missing:
            raise ValueError(
                f"Expression file is missing required column(s): {', '.join(missing)}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        available_profile_type_column = (
            profile_type_column
            if profile_type_column and profile_type_column in reader.fieldnames
            else None
        )
        inferred_profile_type = ""
        if available_profile_type_column is None:
            if profile_column.casefold() == "cell line":
                inferred_profile_type = "cell_line"
            elif profile_column.casefold() == "tissue":
                inferred_profile_type = "tissue"
        for row in reader:
            gene_id = normalize_gene_id(row.get(gene_column, ""))
            profile = row.get(profile_column, "").strip()
            if not gene_id or not profile:
                invalid += 1
                continue
            try:
                value = float(row.get(value_column, ""))
            except (TypeError, ValueError):
                invalid += 1
                continue
            if not math.isfinite(value) or value < 0:
                invalid += 1
                continue
            key = (gene_id, profile)
            if counts[key]:
                duplicates += 1
            sums[key] += value
            counts[key] += 1
            profiles.add(profile)
            raw_profile_type = (
                row.get(available_profile_type_column, "")
                if available_profile_type_column
                else inferred_profile_type
            )
            profile_type = _normalise_profile_type(raw_profile_type or "")
            if profile_type:
                profile_type_candidates[profile].add(profile_type)
            name = row.get(name_column, "").strip()
            if name and gene_id not in gene_names:
                gene_names[gene_id] = name

    if not sums:
        raise ValueError(f"No valid expression values were found in {path}")

    profile_list = sorted(profiles, key=str.casefold)
    conflicting_types = {
        profile: sorted(values)
        for profile, values in profile_type_candidates.items()
        if len(values) > 1
    }
    if conflicting_types:
        profile, values = next(iter(conflicting_types.items()))
        raise ValueError(
            f"Expression profile {profile!r} has conflicting profile types: {', '.join(values)}"
        )
    profile_types = {
        profile: next(iter(profile_type_candidates.get(profile, set())), "")
        for profile in profile_list
    }
    values: dict[str, dict[str, float]] = {profile: {} for profile in profile_list}
    nonzero_by_gene: dict[str, int] = defaultdict(int)
    all_genes: set[str] = set()
    for (gene_id, profile), total in sums.items():
        value = total / counts[(gene_id, profile)]
        values[profile][gene_id] = value
        all_genes.add(gene_id)
        if value > 0:
            nonzero_by_gene[gene_id] += 1

    effective_min = min(max(1, min_nonzero_profiles), len(profile_list))
    keep_genes = {gene for gene in all_genes if nonzero_by_gene.get(gene, 0) >= effective_min}
    filtered = len(all_genes) - len(keep_genes)
    for profile in profile_list:
        values[profile] = {
            gene: value for gene, value in values[profile].items() if gene in keep_genes
        }

    return ExpressionData(
        values=values,
        gene_names=gene_names,
        profiles=profile_list,
        profile_types=profile_types,
        genes=keep_genes,
        duplicate_records=duplicates,
        invalid_records=invalid,
        filtered_gene_count=filtered,
    )


def parse_named_paths(values: Sequence[str], *, option: str) -> list[NamedPath]:
    parsed: list[NamedPath] = []
    seen: set[str] = set()
    for value in values:
        if "=" in value:
            name, path_text = value.split("=", 1)
            name = name.strip()
            path_text = path_text.strip()
        else:
            path_text = value.strip()
            name = Path(path_text).stem
        if not name or not path_text:
            raise ValueError(f"{option} entries must be NAME=PATH or PATH")
        if name in seen:
            raise ValueError(f"Duplicate sample name for {option}: {name}")
        path = Path(path_text)
        if not path.is_file():
            raise FileNotFoundError(f"{option} file not found: {path}")
        parsed.append(NamedPath(name, path))
        seen.add(name)
    return parsed


def read_peak_positions(
    path: str | Path,
    *,
    position_column: int,
    blacklist: BlacklistIndex | None = None,
) -> dict[str, np.ndarray]:
    pos_idx = _column_index(position_column, "--peak-position-column")
    by_chrom: dict[str, list[int]] = defaultdict(list)
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            text = raw.strip()
            if not text or text.startswith(HEADER_PREFIXES):
                continue
            fields = split_fields(raw)
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: expected BED3 or greater")
            try:
                start = int(fields[1])
                end = int(fields[2])
                position = int(fields[pos_idx]) if pos_idx < len(fields) else (start + end) // 2
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: non-integer peak coordinate") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}: invalid peak interval")
            if blacklist is not None and blacklist.overlaps(fields[0], start, end):
                continue
            by_chrom[fields[0]].append(position)
    if not by_chrom:
        raise ValueError(f"No peaks were found in {path}")
    return {chrom: np.asarray(sorted(set(values)), dtype=np.int64) for chrom, values in by_chrom.items()}


def contig_candidates(chrom: str) -> tuple[str, ...]:
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append(f"chr{chrom}")
    if chrom in {"chrM", "M"}:
        candidates.extend(["chrMT", "MT"])
    elif chrom in {"chrMT", "MT"}:
        candidates.extend(["chrM", "M"])
    return tuple(dict.fromkeys(candidates))


def resolve_contig(chrom: str, available: Mapping[str, object]) -> str | None:
    return next((candidate for candidate in contig_candidates(chrom) if candidate in available), None)


def positions_in_interval(positions: np.ndarray, start: int, end: int) -> np.ndarray:
    left = int(np.searchsorted(positions, start, side="left"))
    right = int(np.searchsorted(positions, end, side="left"))
    return positions[left:right]


def median_adjacent_distance(positions: np.ndarray) -> tuple[int, float]:
    count = int(positions.size)
    if count < 2:
        return count, math.nan
    return count, float(np.median(np.diff(positions)))


def gene_regions(gene: GeneRecord, flank: int) -> dict[str, tuple[int, int]]:
    if gene.strand == "+":
        upstream = (max(0, gene.start - flank), gene.start)
        downstream = (gene.end, gene.end + flank)
    else:
        upstream = (gene.end, gene.end + flank)
        downstream = (max(0, gene.start - flank), gene.start)
    return {
        "body": (gene.start, gene.end),
        "upstream": upstream,
        "downstream": downstream,
    }


def transform_expression(value: float, mode: str, floor: float) -> float:
    if mode == "none":
        return value
    if mode == "log2p1":
        return math.log2(value + 1.0)
    if mode == "log2-floor":
        return math.log2(max(value, floor))
    raise ValueError(f"Unknown expression transform: {mode}")


def correlate(x: np.ndarray, y: np.ndarray, *, method: str, min_genes: int) -> CorrelationResult:
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < min_genes:
        return CorrelationResult(math.nan, math.nan, n)
    xv = x[mask]
    yv = y[mask]
    if np.allclose(xv, xv[0]) or np.allclose(yv, yv[0]):
        return CorrelationResult(math.nan, math.nan, n)
    if method == "pearson":
        result = scipy_stats.pearsonr(xv, yv)
    else:
        result = scipy_stats.spearmanr(xv, yv)
    return CorrelationResult(float(result.statistic), float(result.pvalue), n)


def run_spacing_analysis(
    genes: Sequence[GeneRecord],
    expression: ExpressionData,
    peak_specs: Sequence[NamedPath],
    *,
    output_prefix: Path,
    peak_position_column: int,
    flank: int,
    high_confidence_peaks: int,
    expression_transform: str,
    expression_floor: float,
    correlation_method: str,
    min_correlation_genes: int,
    focus_profiles: Sequence[str],
    blacklist: BlacklistIndex | None = None,
) -> dict[str, Path]:
    per_gene_path = Path(f"{output_prefix}_gene_peak_spacing.tsv")
    correlation_path = Path(f"{output_prefix}_spacing_expression_correlations.tsv")
    from nucleosuite.plotting import plot_path as resolve_plot_path
    plot_path = resolve_plot_path(Path(f"{output_prefix}_spacing_expression_correlations.png"))
    scatter_path = resolve_plot_path(Path(f"{output_prefix}_spacing_expression_scatter.png"))
    scatter_data_path = Path(f"{output_prefix}_spacing_expression_scatter.tsv")
    per_gene_path.parent.mkdir(parents=True, exist_ok=True)

    gene_rows: list[dict[str, object]] = []
    correlations: list[dict[str, object]] = []

    for peak_spec in peak_specs:
        peaks_by_chrom = read_peak_positions(
            peak_spec.path,
            position_column=peak_position_column,
            blacklist=blacklist,
        )
        available = peaks_by_chrom
        for gene in genes:
            peak_chrom = resolve_contig(gene.chrom, available)
            row: dict[str, object] = {
                "sample": peak_spec.name,
                "chrom": gene.chrom,
                "start": gene.start,
                "end": gene.end,
                "gene_id": gene.gene_id,
                "gene_name": gene.gene_name,
                "strand": gene.strand,
            }
            for region_name, (start, end) in gene_regions(gene, flank).items():
                if peak_chrom is None or end <= start:
                    count, median = 0, math.nan
                else:
                    selected = positions_in_interval(peaks_by_chrom[peak_chrom], start, end)
                    count, median = median_adjacent_distance(selected)
                row[f"{region_name}_peak_count"] = count
                row[f"{region_name}_median_spacing_bp"] = median
            gene_rows.append(row)

    with per_gene_path.open("wt", newline="") as handle:
        fieldnames = [
            "sample", "chrom", "start", "end", "gene_id", "gene_name", "strand",
            "body_peak_count", "body_median_spacing_bp",
            "upstream_peak_count", "upstream_median_spacing_bp",
            "downstream_peak_count", "downstream_median_spacing_bp",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in gene_rows:
            writer.writerow(row)

    rows_by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in gene_rows:
        rows_by_sample[str(row["sample"])].append(row)

    for sample, rows in rows_by_sample.items():
        gene_ids = [str(row["gene_id"]) for row in rows]
        for region in ("body", "upstream", "downstream"):
            spacing = np.asarray([float(row[f"{region}_median_spacing_bp"]) for row in rows])
            peak_counts = np.asarray([int(row[f"{region}_peak_count"]) for row in rows])
            for subset, subset_mask in (
                ("all", np.ones(len(rows), dtype=bool)),
                (f"at_least_{high_confidence_peaks}_peaks", peak_counts >= high_confidence_peaks),
            ):
                for profile in expression.profiles:
                    expr = np.asarray(
                        [
                            transform_expression(
                                expression.values[profile].get(gene_id, math.nan),
                                expression_transform,
                                expression_floor,
                            )
                            if gene_id in expression.values[profile]
                            else math.nan
                            for gene_id in gene_ids
                        ],
                        dtype=float,
                    )
                    result = correlate(
                        np.where(subset_mask, spacing, np.nan),
                        expr,
                        method=correlation_method,
                        min_genes=min_correlation_genes,
                    )
                    correlations.append(
                        {
                            "sample": sample,
                            "region": region,
                            "subset": subset,
                            "profile": profile,
                            "correlation_method": correlation_method,
                            "correlation": result.r,
                            "p_value": result.p,
                            "matched_genes": result.n,
                            "expression_transform": expression_transform,
                        }
                    )

    with correlation_path.open("wt", newline="") as handle:
        fieldnames = [
            "sample", "region", "subset", "profile", "correlation_method",
            "correlation", "p_value", "matched_genes", "expression_transform",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(correlations)

    # Summary plot: body correlations across all profiles for each sample.
    fig, ax = plt.subplots(figsize=(11, 6))
    plotted = False
    for sample in sorted(rows_by_sample):
        selected = [
            row for row in correlations
            if row["sample"] == sample and row["region"] == "body" and row["subset"] == "all"
            and math.isfinite(float(row["correlation"]))
        ]
        selected.sort(key=lambda row: float(row["correlation"]))
        if selected:
            ax.plot(
                np.arange(1, len(selected) + 1),
                [float(row["correlation"]) for row in selected],
                label=sample,
                marker="o",
                markersize=2.2,
                markeredgewidth=0,
            )
            plotted = True
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Expression profile rank by correlation")
    ax.set_ylabel(f"{correlation_method.title()} correlation: median spacing vs expression")
    ax.set_title("Gene-body nucleosome spacing and expression")
    if plotted:
        ax.legend(frameon=False)
    from nucleosuite.plotting import save_figure
    fig.tight_layout()
    plot_path = save_figure(fig, plot_path, default_dpi=200)
    plt.close(fig)

    # Scatter requested profiles, or the strongest negative body result.
    chosen = list(dict.fromkeys(focus_profiles))
    if not chosen:
        body_rows = [
            row for row in correlations
            if row["region"] == "body" and row["subset"] == "all"
            and math.isfinite(float(row["correlation"]))
        ]
        if body_rows:
            chosen = [str(min(body_rows, key=lambda row: float(row["correlation"]))["profile"])]
    n_panels = max(1, len(rows_by_sample) * max(1, len(chosen)))
    columns = min(3, n_panels)
    rows_n = int(math.ceil(n_panels / columns))
    fig, axes = plt.subplots(rows_n, columns, figsize=(5.2 * columns, 4.5 * rows_n), squeeze=False)
    panel = 0
    scatter_rows: list[dict[str, object]] = []
    for sample, rows in sorted(rows_by_sample.items()):
        for profile in chosen:
            ax = axes.flat[panel]
            panel += 1
            if profile not in expression.values:
                ax.text(0.5, 0.5, f"Profile not found: {profile}", ha="center", va="center")
                ax.set_axis_off()
                continue
            x: list[float] = []
            y: list[float] = []
            for row in rows:
                spacing = float(row["body_median_spacing_bp"])
                gene_id = str(row["gene_id"])
                value = expression.values[profile].get(gene_id)
                if math.isfinite(spacing) and value is not None:
                    x.append(spacing)
                    y.append(transform_expression(value, expression_transform, expression_floor))
            result = correlate(np.asarray(x), np.asarray(y), method=correlation_method, min_genes=min_correlation_genes)
            scatter_rows.extend(
                {
                    "sample": sample,
                    "profile": profile,
                    "body_median_spacing_bp": x_value,
                    "transformed_expression": y_value,
                    "expression_transform": expression_transform,
                    "correlation_method": correlation_method,
                    "correlation": result.r,
                    "matched_genes": result.n,
                }
                for x_value, y_value in zip(x, y)
            )
            ax.scatter(x, y, s=8, alpha=0.35, linewidths=0)
            ax.set_xlabel("Median adjacent peak spacing in gene body (bp)")
            ax.set_ylabel(f"{expression_transform} expression")
            ax.set_title(f"{sample} vs {profile}\nr={result.r:.3f}, n={result.n:,}")
    for ax in axes.flat[panel:]:
        ax.set_axis_off()
    fig.tight_layout()
    scatter_path = save_figure(fig, scatter_path, default_dpi=200)
    plt.close(fig)

    with scatter_data_path.open("wt", newline="") as handle:
        fieldnames = [
            "sample", "profile", "body_median_spacing_bp",
            "transformed_expression", "expression_transform",
            "correlation_method", "correlation", "matched_genes",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(scatter_rows)

    return {
        "spacing_per_gene": per_gene_path,
        "spacing_correlations": correlation_path,
        "spacing_plot": plot_path,
        "spacing_scatter_data": scatter_data_path,
        "spacing_scatter": scatter_path,
    }


def recursive_filter(values: np.ndarray) -> np.ndarray:
    """Apply the 24-coefficient recursive filter used in the Snyder et al. (2016) workflow.

    R ``stats::filter(..., method="recursive")`` follows
    ``y[n] = x[n] + sum(a[k] * y[n-k])``.  The equivalent SciPy IIR
    denominator is therefore ``[1, -a1, -a2, ...]``.
    """
    coefficients = 1.0 / np.arange(5.0, 101.0, 4.0)
    warm_count = min(300, values.size)
    extended = np.concatenate([values[:warm_count], values]).astype(float, copy=False)
    denominator = np.concatenate(([1.0], -coefficients))
    filtered = scipy_signal.lfilter([1.0], denominator, extended)
    return filtered[warm_count:]


def periodogram_intensities(
    values: np.ndarray,
    periods: np.ndarray,
    *,
    trim_fraction: float,
    pad_fraction: float,
    taper_fraction: float,
    recursive: bool,
    valid_mask: np.ndarray | None = None,
) -> np.ndarray:
    trajectory = np.asarray(values, dtype=float)
    mask = None if valid_mask is None else np.asarray(valid_mask, dtype=bool)
    if mask is not None and mask.shape != trajectory.shape:
        raise ValueError("FFT valid mask must match the signal trajectory")
    if mask is not None and np.all(mask):
        mask = None
    if recursive:
        if mask is None:
            trajectory = recursive_filter(trajectory)
        else:
            # Filter each retained run independently so a blacklisted value can
            # neither contribute directly nor propagate through the IIR state.
            filtered = np.zeros_like(trajectory)
            padded = np.concatenate(([False], mask, [False])).astype(np.int8)
            boundaries = np.flatnonzero(np.diff(padded))
            for start, end in boundaries.reshape(-1, 2):
                filtered[start:end] = recursive_filter(trajectory[start:end])
            trajectory = filtered
    if mask is None:
        # Preserve the pre-blacklist numerical path exactly.
        trajectory = trajectory - float(
            scipy_stats.trim_mean(trajectory, trim_fraction)
        )
        trajectory = scipy_signal.detrend(trajectory, type="linear")
        trajectory = trajectory - float(np.mean(trajectory))
        normalization_count = trajectory.size
    else:
        if np.count_nonzero(mask) < 3:
            raise ValueError("At least three unmasked signal positions are required")
        coordinates = np.arange(trajectory.size, dtype=float)
        trajectory = trajectory - float(
            scipy_stats.trim_mean(trajectory[mask], trim_fraction)
        )
        slope, intercept = np.polyfit(coordinates[mask], trajectory[mask], 1)
        trajectory = trajectory - (slope * coordinates + intercept)
        trajectory = trajectory - float(np.mean(trajectory[mask]))
        normalization_count = int(np.count_nonzero(mask))
    alpha = min(1.0, max(0.0, 2.0 * taper_fraction))
    trajectory = trajectory * scipy_signal.windows.tukey(trajectory.size, alpha=alpha)
    if mask is not None:
        # Zero here means no weighted contribution after centring, not a zero
        # genomic signal value. Blacklisted positions are excluded from the DFT.
        trajectory = np.where(mask, trajectory, 0.0)
    nfft = scipy_fft.next_fast_len(int(math.ceil(trajectory.size * (1.0 + pad_fraction))))
    frequencies = np.fft.rfftfreq(nfft, d=1.0)
    spectrum = np.abs(np.fft.rfft(trajectory, n=nfft)) ** 2 / max(
        1, normalization_count
    )
    # Modified Daniell kernel with m=2, equivalent to R spec.pgram(span=2).
    kernel = np.asarray([0.125, 0.25, 0.25, 0.25, 0.125], dtype=float)
    spectrum = np.convolve(spectrum, kernel, mode="same")
    keep = frequencies > 0
    source_periods = 1.0 / frequencies[keep]
    source_spectrum = spectrum[keep]
    order = np.argsort(source_periods)
    return np.interp(periods.astype(float), source_periods[order], source_spectrum[order])


def extract_gene_signal(
    bw,
    gene: GeneRecord,
    *,
    window_size: int,
    missing_value: float,
    blacklist: BlacklistIndex | None = None,
) -> tuple[np.ndarray | None, float, np.ndarray | None]:
    chroms = bw.chroms()
    chrom = resolve_contig(gene.chrom, chroms)
    if chrom is None:
        return None, 0.0, None
    chrom_length = int(chroms[chrom])
    if gene.strand == "+":
        start = gene.start
        end = start + window_size
    else:
        end = gene.end
        start = end - window_size
    if start < 0 or end > chrom_length or end <= start:
        return None, 0.0, None
    values = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
    if values.size != window_size:
        return None, 0.0, None
    blacklist_valid = (
        blacklist.valid_mask(chrom, start, end)
        if blacklist is not None
        else None
    )
    finite = np.isfinite(values)
    valid = finite if blacklist_valid is None else finite & blacklist_valid
    valid_fraction = float(valid.mean())
    values[~finite] = missing_value
    if blacklist_valid is not None and not np.any(blacklist_valid):
        return None, 0.0, blacklist_valid
    if gene.strand == "-":
        values = values[::-1]
        if blacklist_valid is not None:
            blacklist_valid = blacklist_valid[::-1]
    return values, valid_fraction, blacklist_valid


def _expression_vector(
    expression: ExpressionData,
    profile: str,
    gene_ids: Sequence[str],
    *,
    transform: str,
    floor: float,
) -> np.ndarray:
    profile_values = expression.values[profile]
    return np.asarray(
        [
            transform_expression(profile_values[gene], transform, floor)
            if gene in profile_values
            else math.nan
            for gene in gene_ids
        ],
        dtype=float,
    )


def run_fft_analysis(
    genes: Sequence[GeneRecord],
    expression: ExpressionData,
    signal_specs: Sequence[NamedPath],
    *,
    output_prefix: Path,
    window_size: int,
    period_min: int,
    period_max: int,
    ranking_periods: Sequence[int],
    min_valid_fraction: float,
    missing_value: float,
    trim_fraction: float,
    pad_fraction: float,
    taper_fraction: float,
    recursive_filter_enabled: bool,
    expression_transform: str,
    expression_floor: float,
    correlation_method: str,
    min_correlation_genes: int,
    focus_profiles: Sequence[str],
    control_samples: Sequence[str],
    top_profiles: int,
    profile_metadata: Mapping[str, Mapping[str, str]] | None = None,
    blacklist: BlacklistIndex | None = None,
) -> dict[str, Path]:
    periods = np.arange(period_min, period_max + 1, dtype=int)
    ranking_indices = [int(np.where(periods == period)[0][0]) for period in ranking_periods]
    fft_path = Path(f"{output_prefix}_per_gene_fft.tsv.gz")
    trajectory_path = Path(f"{output_prefix}_fft_expression_correlations.tsv")
    ranking_path = Path(f"{output_prefix}_expression_profile_rankings.tsv")
    rank_change_path = Path(f"{output_prefix}_expression_rank_changes.tsv")
    from nucleosuite.plotting import plot_path as resolve_plot_path
    trajectory_plot_path = resolve_plot_path(Path(f"{output_prefix}_fft_expression_correlations.png"))
    ranking_plot_path = resolve_plot_path(Path(f"{output_prefix}_expression_profile_rankings.png"))
    qc_path = Path(f"{output_prefix}_fft_qc.tsv")
    fft_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("pyBigWig is required for gene-expression FFT analysis") from exc

    all_rows: list[FFTRow] = []
    qc_rows: list[dict[str, object]] = []
    for spec in signal_specs:
        used = 0
        skipped_missing_contig = 0
        skipped_boundary = 0
        skipped_valid = 0
        with pyBigWig.open(str(spec.path)) as bw:
            chroms = bw.chroms()
            for index, gene in enumerate(genes, 1):
                chrom = resolve_contig(gene.chrom, chroms)
                if chrom is None:
                    skipped_missing_contig += 1
                    continue
                values, valid_fraction, blacklist_valid = extract_gene_signal(
                    bw,
                    gene,
                    window_size=window_size,
                    missing_value=missing_value,
                    blacklist=blacklist,
                )
                if values is None:
                    skipped_boundary += 1
                    continue
                if valid_fraction < min_valid_fraction:
                    skipped_valid += 1
                    continue
                intensities = periodogram_intensities(
                    values,
                    periods,
                    trim_fraction=trim_fraction,
                    pad_fraction=pad_fraction,
                    taper_fraction=taper_fraction,
                    recursive=recursive_filter_enabled,
                    valid_mask=blacklist_valid,
                )
                dominant_index = int(np.argmax(intensities))
                all_rows.append(
                    FFTRow(
                        sample=spec.name,
                        gene=gene,
                        valid_fraction=valid_fraction,
                        dominant_period=int(periods[dominant_index]),
                        dominant_intensity=float(intensities[dominant_index]),
                        ranking_intensity=float(np.mean(intensities[ranking_indices])),
                        intensities=intensities,
                    )
                )
                used += 1
                if index % 1000 == 0:
                    print(f"FFT {spec.name}: processed {index:,}/{len(genes):,} genes", file=sys.stderr)
        qc_rows.append(
            {
                "sample": spec.name,
                "genes_total": len(genes),
                "genes_used": used,
                "skipped_missing_contig": skipped_missing_contig,
                "skipped_boundary": skipped_boundary,
                "skipped_low_valid_fraction": skipped_valid,
            }
        )

    if not all_rows:
        raise ValueError("No genes passed the FFT signal extraction filters")

    with gzip.open(fft_path, "wt", newline="") as handle:
        fieldnames = [
            "sample", "chrom", "start", "end", "gene_id", "gene_name", "strand",
            "valid_fraction", "dominant_period_bp", "dominant_intensity",
            "ranking_periods", "ranking_mean_intensity",
        ] + [f"period_{period}_bp" for period in periods]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in all_rows:
            record: dict[str, object] = {
                "sample": row.sample,
                "chrom": row.gene.chrom,
                "start": row.gene.start,
                "end": row.gene.end,
                "gene_id": row.gene.gene_id,
                "gene_name": row.gene.gene_name,
                "strand": row.gene.strand,
                "valid_fraction": f"{row.valid_fraction:.6f}",
                "dominant_period_bp": row.dominant_period,
                "dominant_intensity": f"{row.dominant_intensity:.12g}",
                "ranking_periods": ",".join(map(str, ranking_periods)),
                "ranking_mean_intensity": f"{row.ranking_intensity:.12g}",
            }
            record.update(
                {f"period_{period}_bp": f"{value:.12g}" for period, value in zip(periods, row.intensities)}
            )
            writer.writerow(record)

    rows_by_sample: dict[str, list[FFTRow]] = defaultdict(list)
    for row in all_rows:
        rows_by_sample[row.sample].append(row)

    trajectory_rows: list[dict[str, object]] = []
    ranking_rows: list[dict[str, object]] = []
    ranks_by_sample: dict[str, dict[str, int]] = {}

    for sample, rows in rows_by_sample.items():
        gene_ids = [row.gene.gene_id for row in rows]
        matrix = np.vstack([row.intensities for row in rows])
        band = np.asarray([row.ranking_intensity for row in rows], dtype=float)
        sample_rank_rows: list[dict[str, object]] = []
        for profile in expression.profiles:
            expr = _expression_vector(
                expression,
                profile,
                gene_ids,
                transform=expression_transform,
                floor=expression_floor,
            )
            for period_index, period in enumerate(periods):
                result = correlate(
                    matrix[:, period_index],
                    expr,
                    method=correlation_method,
                    min_genes=min_correlation_genes,
                )
                trajectory_rows.append(
                    {
                        "sample": sample,
                        "profile": profile,
                        "period_bp": int(period),
                        "correlation_method": correlation_method,
                        "correlation": result.r,
                        "p_value": result.p,
                        "matched_genes": result.n,
                    }
                )
            result = correlate(
                band,
                expr,
                method=correlation_method,
                min_genes=min_correlation_genes,
            )
            sample_rank_rows.append(
                {
                    "sample": sample,
                    "profile": profile,
                    "correlation": result.r,
                    "p_value": result.p,
                    "matched_genes": result.n,
                    "ranking_periods": ",".join(map(str, ranking_periods)),
                    **annotate_profile(
                        profile,
                        profile_metadata,
                        profile_type=expression.profile_types.get(profile, ""),
                    ),
                }
            )
        finite_rows = [row for row in sample_rank_rows if math.isfinite(float(row["correlation"]))]
        finite_rows.sort(key=lambda row: (float(row["correlation"]), str(row["profile"]).casefold()))
        ranks = {str(row["profile"]): rank for rank, row in enumerate(finite_rows, 1)}
        ranks_by_sample[sample] = ranks
        for row in sample_rank_rows:
            row["rank"] = ranks.get(str(row["profile"]), "")
            row["plot_selected"] = bool(
                row["rank"] != "" and int(row["rank"]) <= top_profiles
            )
        ranking_rows.extend(sample_rank_rows)

    highlighted_by_sample: dict[str, set[str]] = {}
    for sample in rows_by_sample:
        available = {
            str(row["profile"])
            for row in trajectory_rows
            if row["sample"] == sample
        }
        selected_profiles = [profile for profile in focus_profiles if profile in available]
        if not selected_profiles:
            sample_ranking = [
                row for row in ranking_rows
                if row["sample"] == sample and row["rank"] != ""
            ]
            sample_ranking.sort(key=lambda row: int(row["rank"]))
            selected_profiles = [str(row["profile"]) for row in sample_ranking[:1]]
        highlighted_by_sample[sample] = set(selected_profiles)
    for row in trajectory_rows:
        row["plot_highlight"] = str(row["profile"]) in highlighted_by_sample[str(row["sample"])]
        row["ranking_periods"] = ",".join(map(str, ranking_periods))

    with trajectory_path.open("wt", newline="") as handle:
        fieldnames = [
            "sample", "profile", "period_bp", "correlation_method",
            "correlation", "p_value", "matched_genes", "plot_highlight",
            "ranking_periods",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(trajectory_rows)

    with ranking_path.open("wt", newline="") as handle:
        fieldnames = [
            "sample", "rank", "profile", "correlation", "p_value", "matched_genes", "ranking_periods", "plot_selected",
            "profile_type",
            "cancer_cell_line", "primary_disease", "disease_subtype",
            "primary_or_metastasis", "sample_collection_site", "cellosaurus_id",
            "metadata_matched", "metadata_status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        ranking_rows.sort(
            key=lambda row: (
                str(row["sample"]),
                int(row["rank"]) if row["rank"] != "" else 10**9,
                str(row["profile"]).casefold(),
            )
        )
        writer.writerows(ranking_rows)

    controls = [sample for sample in control_samples if sample in ranks_by_sample]
    rank_change_rows: list[dict[str, object]] = []
    if controls:
        all_profiles = sorted(set().union(*(set(ranks_by_sample[s]) for s in controls)), key=str.casefold)
        control_average = {
            profile: float(np.mean([ranks_by_sample[s][profile] for s in controls if profile in ranks_by_sample[s]]))
            for profile in all_profiles
        }
        for sample, ranks in ranks_by_sample.items():
            for profile, rank in ranks.items():
                if profile not in control_average:
                    continue
                rank_change_rows.append(
                    {
                        "sample": sample,
                        "profile": profile,
                        "sample_rank": rank,
                        "control_mean_rank": control_average[profile],
                        "rank_increase": control_average[profile] - rank,
                        "control_samples": ",".join(controls),
                    }
                )
    with rank_change_path.open("wt", newline="") as handle:
        fieldnames = [
            "sample", "profile", "sample_rank", "control_mean_rank", "rank_increase", "control_samples",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rank_change_rows, key=lambda r: (str(r["sample"]), -float(r["rank_increase"]))))

    with qc_path.open("wt", newline="") as handle:
        fieldnames = list(qc_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)

    # Plot all trajectories in grey and highlight selected/focus profiles.
    sample_names = sorted(rows_by_sample)
    fig, axes = plt.subplots(len(sample_names), 1, figsize=(10, 5 * len(sample_names)), squeeze=False)
    for ax, sample in zip(axes.flat, sample_names):
        rows = [row for row in trajectory_rows if row["sample"] == sample]
        by_profile: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in rows:
            by_profile[str(row["profile"])].append(row)
        selected_profiles = [profile for profile in focus_profiles if profile in by_profile]
        if not selected_profiles:
            sample_ranking = [row for row in ranking_rows if row["sample"] == sample and row["rank"] != ""]
            sample_ranking.sort(key=lambda row: int(row["rank"]))
            selected_profiles = [str(row["profile"]) for row in sample_ranking[:1]]
        for profile, profile_rows in by_profile.items():
            profile_rows.sort(key=lambda row: int(row["period_bp"]))
            ax.plot(
                [int(row["period_bp"]) for row in profile_rows],
                [float(row["correlation"]) for row in profile_rows],
                linewidth=0.6,
                alpha=0.25,
                color="0.5",
                marker="o",
                markersize=1.2,
                markeredgewidth=0,
            )
        for profile in selected_profiles:
            profile_rows = sorted(by_profile[profile], key=lambda row: int(row["period_bp"]))
            ax.plot(
                [int(row["period_bp"]) for row in profile_rows],
                [float(row["correlation"]) for row in profile_rows],
                linewidth=2.0,
                marker="o",
                markersize=2.5,
                label=profile,
            )
        ax.axhline(0, linewidth=0.8)
        ax.axvspan(min(ranking_periods), max(ranking_periods), alpha=0.08)
        from nucleosuite.plotting import apply_base_pair_x_axis
        apply_base_pair_x_axis(ax, range(period_min, period_max + 1))
        ax.set_xlabel("Nucleosome period (bp)")
        ax.set_ylabel(f"{correlation_method.title()} correlation")
        ax.set_title(f"{sample}: FFT intensity versus expression")
        if selected_profiles:
            ax.legend(frameon=False)
    from nucleosuite.plotting import save_figure
    fig.tight_layout()
    trajectory_plot_path = save_figure(fig, trajectory_plot_path, default_dpi=200)
    plt.close(fig)

    # Plot profile rankings for each signal sample.
    fig, axes = plt.subplots(len(sample_names), 1, figsize=(10, max(5, 0.28 * top_profiles) * len(sample_names)), squeeze=False)
    for ax, sample in zip(axes.flat, sample_names):
        selected = [row for row in ranking_rows if row["sample"] == sample and row["rank"] != ""]
        selected.sort(key=lambda row: int(row["rank"]))
        selected = selected[:top_profiles]
        labels = [str(row["profile"]) for row in selected][::-1]
        values = [float(row["correlation"]) for row in selected][::-1]
        ax.barh(np.arange(len(selected)), values)
        ax.set_yticks(np.arange(len(selected)), labels=labels)
        ax.axvline(0, linewidth=0.8)
        ax.set_xlabel(
            f"{correlation_method.title()} correlation with mean FFT intensity at "
            + ", ".join(map(str, ranking_periods))
            + " bp"
        )
        ax.set_title(f"{sample}: expression profiles ranked from most negative correlation")
    fig.tight_layout()
    ranking_plot_path = save_figure(fig, ranking_plot_path, default_dpi=200)
    plt.close(fig)

    return {
        "fft_per_gene": fft_path,
        "fft_correlations": trajectory_path,
        "fft_rankings": ranking_path,
        "fft_rank_changes": rank_change_path,
        "fft_correlation_plot": trajectory_plot_path,
        "fft_ranking_plot": ranking_plot_path,
        "fft_qc": qc_path,
    }


def write_metadata(
    output_prefix: Path,
    *,
    args: argparse.Namespace,
    gene_path: Path,
    genes: Sequence[GeneRecord],
    expression: ExpressionData,
    outputs: Mapping[str, Path],
) -> Path:
    path = Path(f"{output_prefix}_metadata.tsv")
    rows = [
        ("analysis", args.analysis),
        ("signal_type", args.signal_type),
        ("gene_bed", gene_path),
        ("gene_count", len(genes)),
        ("blacklist_bed", args.blacklist_bed or ""),
        (
            "blacklist_overlapping_gene_anchors_excluded",
            getattr(args, "_blacklisted_gene_anchors", 0),
        ),
        ("expression", args.expression),
        ("expression_value_column", args.expression_value_column),
        ("expression_profile_column", args.expression_profile_column),
        ("expression_profile_type_column", args.expression_profile_type_column),
        ("expression_profile_count", len(expression.profiles)),
        ("expression_cell_line_profile_count", sum(value == "cell_line" for value in expression.profile_types.values())),
        ("expression_tissue_profile_count", sum(value == "tissue" for value in expression.profile_types.values())),
        ("expression_gene_count", len(expression.genes)),
        ("expression_duplicate_records_averaged", expression.duplicate_records),
        ("expression_invalid_records_skipped", expression.invalid_records),
        ("expression_genes_filtered_for_nonzero_profiles", expression.filtered_gene_count),
        ("spacing_expression_transform", args.spacing_expression_transform),
        ("fft_expression_transform", args.fft_expression_transform),
        ("fft_window", args.fft_window),
        ("fft_period_min", args.fft_period_min),
        ("fft_period_max", args.fft_period_max),
        ("fft_ranking_periods", args.fft_ranking_periods),
    ]
    for name, output in outputs.items():
        rows.append((f"output_{name}", output))
    with path.open("wt") as handle:
        handle.write("parameter\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite gene-expression",
        description=(
            "Relate PNS or WPS peak spacing and per-gene FFT intensity to long-format "
            "gene-expression profiles. PNS is the default signal type."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--expression", required=True, help="Long-format tab-delimited expression table.")
    parser.add_argument("--genes-bed", help="Gene BED. Default: genes from --resource-set.")
    parser.add_argument(
        "--blacklist-bed",
        help=(
            "BED blacklist. Genes whose one-base TSS anchor overlaps it and "
            "overlapping peaks are excluded; blacklisted BigWig positions are "
            "masked from FFT power calculations."
        ),
    )
    parser.add_argument("--resource-set", default="hg19-gm12878", help="Bundled resource set used when --genes-bed is omitted.")
    parser.add_argument("--peaks", "--peak-bed", action="append", default=[], metavar="NAME=BED", help="Peak BED for spacing analysis; repeat for multiple samples.")
    parser.add_argument("--signal", "--bigwig", action="append", default=[], metavar="NAME=BIGWIG", help="PNS/WPS BigWig for FFT analysis; repeat for multiple samples.")
    parser.add_argument("--signal-type", choices=("pns", "wps", "other"), default="pns", help="Signal label used in metadata and plot text.")
    parser.add_argument("--analysis", choices=("all", "spacing", "fft"), default="all", help="Run peak-spacing correlations, FFT correlations, or both.")
    parser.add_argument("--output-prefix", "--out-prefix", required=True, help="Path prefix for analysis tables, metadata, and plots.")

    parser.add_argument("--gene-id-column", type=int, default=4, help="One-based gene BED column containing the gene identifier.")
    parser.add_argument("--gene-name-column", type=int, default=5, help="One-based gene BED column containing the gene name.")
    parser.add_argument("--gene-strand-column", type=int, default=6, help="One-based gene BED column containing + or - strand.")
    parser.add_argument("--peak-position-column", type=int, default=7, help="One-based peak-centre column; midpoint is used when unavailable.")
    parser.add_argument("--gene-flank", type=int, default=10000, help="Upstream/downstream flank used for spacing comparisons.")
    parser.add_argument("--high-confidence-peaks", type=int, default=60, help="Peak-count threshold for the high-confidence spacing subset.")

    parser.add_argument("--expression-gene-column", default="Gene", help="Expression-table gene-ID column.")
    parser.add_argument("--expression-name-column", default="Gene name", help="Expression-table gene-name column.")
    parser.add_argument("--expression-profile-column", default="Cell line", help="Expression-table profile, tissue, or cell-line column.")
    parser.add_argument(
        "--expression-profile-type-column",
        default="Profile type",
        help=(
            "Optional column identifying profile classes such as cell_line or tissue. "
            "The default column is used when present and ignored when absent."
        ),
    )
    parser.add_argument("--expression-value-column", default="nTPM", help="Numeric expression column used in correlations (default: nTPM).")
    parser.add_argument("--profile-metadata", help="Optional cell-line metadata TSV. Cell-line analyses default to the bundled resource.")
    parser.add_argument("--no-profile-metadata", action="store_true", help="Do not annotate profile-ranking rows with bundled cell-line metadata.")
    parser.add_argument("--min-nonzero-profiles-per-gene", type=int, default=3, help="Minimum profiles with expression above zero required to retain a gene (default: 3).")
    parser.add_argument("--expression-floor", type=float, default=0.04, help="Positive floor used by the log2-floor transform (default: 0.04).")
    parser.add_argument("--spacing-expression-transform", choices=("none", "log2p1", "log2-floor"), default="log2p1", help="Expression transform for peak-spacing correlations (default: log2p1).")
    parser.add_argument("--fft-expression-transform", choices=("none", "log2p1", "log2-floor"), default="log2-floor", help="Expression transform for FFT correlations (default: log2-floor).")
    parser.add_argument("--focus-profile", action="append", default=[], help="Expression profile to highlight in plots; may be repeated.")
    parser.add_argument("--correlation", choices=("pearson", "spearman"), default="pearson", help="Correlation statistic (default: pearson).")
    parser.add_argument("--min-correlation-genes", type=int, default=30, help="Minimum matched genes required to report a correlation (default: 30).")

    parser.add_argument("--fft-window", type=int, default=10000, help="Fixed strand-aware window beginning at the TSS.")
    parser.add_argument("--fft-period-min", type=int, default=120, help="Minimum integer period evaluated in bp (default: 120).")
    parser.add_argument("--fft-period-max", type=int, default=280, help="Maximum integer period evaluated in bp (default: 280).")
    parser.add_argument("--fft-ranking-periods", default="193,196,199", help="Comma-separated periods averaged for profile ranking.")
    parser.add_argument("--fft-min-valid-fraction", type=float, default=0.9, help="Minimum finite BigWig fraction required within a gene FFT window (default: 0.9).")
    parser.add_argument("--fft-missing-value", type=float, default=0.0, help="Value used for ordinary missing bases after the valid-fraction filter (default: 0).")
    parser.add_argument("--fft-trim-fraction", type=float, default=0.1, help="Fraction trimmed from each tail when centring the FFT signal (default: 0.1).")
    parser.add_argument("--fft-pad-fraction", type=float, default=0.3, help="Zero-padding fraction appended for spectral evaluation (default: 0.3).")
    parser.add_argument("--fft-taper-fraction", type=float, default=0.3, help="Fraction of the signal covered by the split-cosine taper (default: 0.3).")
    parser.add_argument("--no-fft-recursive-filter", action="store_true", help="Use the detrended and tapered signal without recursive background filtering.")
    parser.add_argument("--control-sample", action="append", default=[], help="Signal sample used as a rank-change control; may be repeated.")
    parser.add_argument("--top-profiles", type=int, default=30, help="Number of ranked expression profiles displayed in summary plots (default: 30).")
    add_parallel_arguments(
        parser,
        cores_option="--memory-intensive-analysis-cores",
        cores_help=(
            "Concurrent contig workers for this memory-intensive gene analysis. "
            "This budget defaults independently to 1."
        ),
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def _run_serial(args: argparse.Namespace) -> int:
    if args.gene_flank < 0:
        raise ValueError("--gene-flank must be non-negative")
    if args.high_confidence_peaks < 2:
        raise ValueError("--high-confidence-peaks must be at least 2")
    if args.min_nonzero_profiles_per_gene < 1:
        raise ValueError("--min-nonzero-profiles-per-gene must be at least 1")
    if args.min_correlation_genes < 3:
        raise ValueError("--min-correlation-genes must be at least 3")
    if args.expression_floor <= 0:
        raise ValueError("--expression-floor must be greater than 0")
    if args.fft_window < 100:
        raise ValueError("--fft-window must be at least 100 bp")
    if args.fft_period_min < 2 or args.fft_period_max <= args.fft_period_min:
        raise ValueError("Require 2 <= --fft-period-min < --fft-period-max")
    if not 0 <= args.fft_min_valid_fraction <= 1:
        raise ValueError("--fft-min-valid-fraction must be between 0 and 1")
    if not 0 <= args.fft_trim_fraction < 0.5:
        raise ValueError("--fft-trim-fraction must be in [0, 0.5)")
    if args.fft_pad_fraction < 0:
        raise ValueError("--fft-pad-fraction must be non-negative")
    if not 0 <= args.fft_taper_fraction <= 0.5:
        raise ValueError("--fft-taper-fraction must be between 0 and 0.5")
    if args.top_profiles < 1:
        raise ValueError("--top-profiles must be at least 1")

    reporter = ProgressReporter("gene-expression")
    peak_specs = parse_named_paths(args.peaks, option="--peaks")
    signal_specs = parse_named_paths(args.signal, option="--signal")
    if args.analysis in {"all", "spacing"} and not peak_specs:
        raise ValueError("--peaks is required for spacing analysis")
    if args.analysis in {"all", "fft"} and not signal_specs:
        raise ValueError("--signal is required for FFT analysis")

    ranking_periods = [int(token) for token in args.fft_ranking_periods.split(",") if token.strip()]
    if not ranking_periods:
        raise ValueError("--fft-ranking-periods must contain at least one integer")
    if any(period < args.fft_period_min or period > args.fft_period_max for period in ranking_periods):
        raise ValueError("Every --fft-ranking-periods value must be inside the FFT period range")

    from nucleosuite.output_naming import parameter_range, parameterized_prefix

    output_prefix = parameterized_prefix(
        args.output_prefix,
        (
            ("analysis", args.analysis),
            ("flank", args.gene_flank),
            ("fftwin", args.fft_window),
            ("period", parameter_range(args.fft_period_min, args.fft_period_max)),
            ("corr", args.correlation),
        ),
    )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    if args.genes_bed:
        gene_path = Path(args.genes_bed)
        if not gene_path.is_file():
            raise FileNotFoundError(f"Gene BED not found: {gene_path}")
        gene_context = None
    else:
        resource_name = resolve_set_resource_name(args.resource_set, "genes")
        gene_context = materialized_resource_path(resource_name)
        gene_path = gene_context.__enter__()

    try:
        reporter.stage("Loading genes and expression profiles")
        genes = read_genes(
            gene_path,
            gene_id_column=args.gene_id_column,
            gene_name_column=args.gene_name_column,
            strand_column=args.gene_strand_column,
        )
        blacklist = load_blacklist_unbounded(args.blacklist_bed)
        genes, args._blacklisted_gene_anchors = filter_blacklisted_gene_anchors(
            genes, blacklist
        )
        if blacklist is not None:
            if not genes:
                raise ValueError(
                    "No genes remained after blacklist TSS-anchor filtering"
                )
        expression = read_expression(
            args.expression,
            gene_column=args.expression_gene_column,
            name_column=args.expression_name_column,
            profile_column=args.expression_profile_column,
            profile_type_column=args.expression_profile_type_column or None,
            value_column=args.expression_value_column,
            min_nonzero_profiles=args.min_nonzero_profiles_per_gene,
        )
        matched = sum(gene.gene_id in expression.genes for gene in genes)
        if matched < args.min_correlation_genes:
            raise ValueError(
                f"Only {matched} genes match between the gene BED and expression table; "
                f"at least {args.min_correlation_genes} are required"
            )
        print(
            f"Genes: {len(genes):,}; expression genes: {len(expression.genes):,}; "
            f"matched: {matched:,}; profiles: {len(expression.profiles):,}",
            file=sys.stderr,
        )

        profile_metadata = None
        metadata_context = None
        should_load_profile_metadata = (
            not args.no_profile_metadata
            and (
                bool(args.profile_metadata)
                or args.expression_profile_column.casefold() == "cell line"
                or any(value == "cell_line" for value in expression.profile_types.values())
            )
        )
        if should_load_profile_metadata:
            if args.profile_metadata:
                metadata_path = Path(args.profile_metadata)
                if not metadata_path.is_file():
                    raise FileNotFoundError(f"Profile metadata not found: {metadata_path}")
            else:
                metadata_name = resolve_set_resource_name(args.resource_set, "cell_line_metadata")
                metadata_context = materialized_resource_path(metadata_name)
                metadata_path = metadata_context.__enter__()
            profile_metadata = read_profile_metadata(metadata_path)

        outputs: dict[str, Path] = {}
        if args.analysis in {"all", "spacing"}:
            reporter.stage(
                f"Calculating peak spacing for {len(peak_specs):,} peak callset(s)"
            )
            outputs.update(
                run_spacing_analysis(
                    genes,
                    expression,
                    peak_specs,
                    output_prefix=output_prefix,
                    peak_position_column=args.peak_position_column,
                    flank=args.gene_flank,
                    high_confidence_peaks=args.high_confidence_peaks,
                    expression_transform=args.spacing_expression_transform,
                    expression_floor=args.expression_floor,
                    correlation_method=args.correlation,
                    min_correlation_genes=args.min_correlation_genes,
                    focus_profiles=args.focus_profile,
                    blacklist=blacklist,
                )
            )
        if args.analysis in {"all", "fft"}:
            reporter.stage(
                f"Calculating per-gene spectra for {len(signal_specs):,} signal track(s)"
            )
            outputs.update(
                run_fft_analysis(
                    genes,
                    expression,
                    signal_specs,
                    output_prefix=output_prefix,
                    window_size=args.fft_window,
                    period_min=args.fft_period_min,
                    period_max=args.fft_period_max,
                    ranking_periods=ranking_periods,
                    min_valid_fraction=args.fft_min_valid_fraction,
                    missing_value=args.fft_missing_value,
                    trim_fraction=args.fft_trim_fraction,
                    pad_fraction=args.fft_pad_fraction,
                    taper_fraction=args.fft_taper_fraction,
                    recursive_filter_enabled=not args.no_fft_recursive_filter,
                    expression_transform=args.fft_expression_transform,
                    expression_floor=args.expression_floor,
                    correlation_method=args.correlation,
                    min_correlation_genes=args.min_correlation_genes,
                    focus_profiles=args.focus_profile,
                    control_samples=args.control_sample,
                    top_profiles=args.top_profiles,
                    profile_metadata=profile_metadata,
                    blacklist=blacklist,
                )
            )
        reporter.stage("Writing gene-expression results and metadata")
        metadata = write_metadata(
            output_prefix,
            args=args,
            gene_path=Path(gene_path),
            genes=genes,
            expression=expression,
            outputs=outputs,
        )
        outputs["metadata"] = metadata
        for name, path in outputs.items():
            print(f"Wrote {name}: {path}", file=sys.stderr)
    finally:
        if 'metadata_context' in locals() and metadata_context is not None:
            metadata_context.__exit__(None, None, None)
        if gene_context is not None:
            gene_context.__exit__(None, None, None)
    return 0


def _run_partitioned_with_genes(args: argparse.Namespace) -> int:
    return run_partitioned_command(
        "gene-expression",
        args,
        _run_serial,
        runner_module="nucleosuite.gene_expression",
        runner_function="_run_serial",
        primary_attr="genes_bed",
        output_prefix_attr="output_prefix",
        named_path_list_attrs=("peaks",),
        path_attrs=("blacklist_bed",),
    )


def run(args: argparse.Namespace) -> int:
    from nucleosuite.output_naming import parameter_range, parameterized_prefix

    args.output_prefix = str(
        parameterized_prefix(
            args.output_prefix,
            (
                ("analysis", args.analysis),
                ("flank", args.gene_flank),
                ("fftwin", args.fft_window),
                ("period", parameter_range(args.fft_period_min, args.fft_period_max)),
                ("corr", args.correlation),
            ),
        )
    )
    cores = int(getattr(args, "cores", 1) or 1)
    if cores <= 1 or getattr(args, "_per_contig_worker", False):
        return _run_serial(args)
    if args.genes_bed:
        return _run_partitioned_with_genes(args)
    resource_name = resolve_set_resource_name(args.resource_set, "genes")
    with materialized_resource_path(resource_name) as gene_path:
        copied = argparse.Namespace(**vars(args))
        copied.genes_bed = str(gene_path)
        return _run_partitioned_with_genes(copied)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
