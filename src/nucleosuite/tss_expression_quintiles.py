#!/usr/bin/env python3
"""Aggregate genomic signals around TSSs after expression-quintile stratification."""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nucleosuite.plotting import configure_unique_category_cycle
configure_unique_category_cycle()
import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.gene_expression import GeneRecord, normalize_gene_id, open_expression_table, read_genes
from nucleosuite.io import open_text
from nucleosuite.resource_files import materialized_resource_path, resolve_set_resource_name
from nucleosuite.progress import ProgressReporter


QUINTILES = (
    ("Q1_lowest", "Low"),
    ("Q2_20_40_percent", "20–40%"),
    ("Q3_middle", "Medium"),
    ("Q4_60_80_percent", "60–80%"),
    ("Q5_highest", "High"),
)


@dataclass(frozen=True)
class ExpressionRecord:
    gene_id: str
    gene_name: str
    value: float


def _normalise_profile_name(value: str) -> str:
    return re.sub(r"[\s_]+", " ", value.strip()).casefold()


def _display_profile_name(value: str) -> str:
    return re.sub(r"[\s_]+", " ", value.strip())


def read_profile_expression(
    path: str | Path,
    *,
    profile: str,
    gene_column: str,
    name_column: str,
    profile_column: str,
    value_column: str,
) -> tuple[str, dict[str, ExpressionRecord], int, int]:
    target = _normalise_profile_name(profile)
    sums: dict[str, float] = {}
    counts: dict[str, int] = {}
    names: dict[str, str] = {}
    matched_profile_names: set[str] = set()
    invalid = 0
    duplicates = 0

    with open_expression_table(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Expression file has no header: {path}")
        missing = [
            column for column in (gene_column, name_column, profile_column, value_column)
            if column not in reader.fieldnames
        ]
        if missing:
            raise ValueError(
                f"Expression file is missing required column(s): {', '.join(missing)}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )
        for row in reader:
            raw_profile = (row.get(profile_column) or "").strip()
            if _normalise_profile_name(raw_profile) != target:
                continue
            matched_profile_names.add(raw_profile)
            gene_id = normalize_gene_id(row.get(gene_column, ""))
            if not gene_id:
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
            if gene_id in counts:
                duplicates += 1
            sums[gene_id] = sums.get(gene_id, 0.0) + value
            counts[gene_id] = counts.get(gene_id, 0) + 1
            name = (row.get(name_column) or "").strip()
            if name and gene_id not in names:
                names[gene_id] = name

    if not matched_profile_names:
        raise ValueError(
            f"Expression profile {profile!r} was not found in column {profile_column!r}. "
            "Use underscores in the CLI value for spaces, for example bone_marrow."
        )
    if len(matched_profile_names) > 1:
        raise ValueError(
            f"Profile selector {profile!r} matched multiple source labels: "
            + ", ".join(sorted(matched_profile_names))
        )
    records = {
        gene_id: ExpressionRecord(
            gene_id=gene_id,
            gene_name=names.get(gene_id, gene_id),
            value=total / counts[gene_id],
        )
        for gene_id, total in sums.items()
    }
    if not records:
        raise ValueError(f"No valid expression values were found for profile {profile!r}")
    return next(iter(matched_profile_names)), records, duplicates, invalid


def assign_quintiles(
    genes: Sequence[GeneRecord],
    expression: dict[str, ExpressionRecord],
) -> list[tuple[str, str, list[tuple[GeneRecord, ExpressionRecord]]]]:
    matched = [(gene, expression[gene.gene_id]) for gene in genes if gene.gene_id in expression]
    matched.sort(key=lambda item: (item[1].value, item[0].gene_id))
    if len(matched) < 5:
        raise ValueError("At least five genes are required to construct expression quintiles")
    indices = np.array_split(np.arange(len(matched), dtype=int), 5)
    return [
        (key, label, [matched[int(index)] for index in group])
        for (key, label), group in zip(QUINTILES, indices)
    ]


def _resolve_contig(name: str, chroms: dict[str, int]) -> str | None:
    if name in chroms:
        return name
    if name.startswith("chr") and name[3:] in chroms:
        return name[3:]
    prefixed = f"chr{name}"
    if prefixed in chroms:
        return prefixed
    return None


def aggregate_quintiles(
    bigwig_path: str | Path,
    quintiles: Sequence[tuple[str, str, list[tuple[GeneRecord, ExpressionRecord]]]],
    *,
    window: int,
    missing_to_zero: bool,
    blacklist_bed: str | Path | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("pyBigWig is required for TSS expression-quintile analysis") from exc

    offsets = np.arange(-window, window + 1, dtype=int)
    profile_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    with pyBigWig.open(str(bigwig_path)) as bw:
        chroms = {str(name): int(length) for name, length in bw.chroms().items()}
        from nucleosuite.core.blacklist import load_blacklist
        blacklist = load_blacklist(
            blacklist_bed, list(chroms), list(chroms.values())
        )
        for key, label, members in quintiles:
            total = np.zeros(offsets.size, dtype=float)
            valid_counts = np.zeros(offsets.size, dtype=np.int64)
            genes_used = 0
            skipped_contig = 0
            skipped_boundary = 0
            skipped_blacklisted_anchor = 0
            expression_values: list[float] = []

            for gene, expr in members:
                chrom = _resolve_contig(gene.chrom, chroms)
                if chrom is None:
                    skipped_contig += 1
                    continue
                if blacklist is not None and blacklist.overlaps(
                    chrom, gene.tss, gene.tss + 1
                ):
                    skipped_blacklisted_anchor += 1
                    continue
                start = gene.tss - window
                end = gene.tss + window + 1
                if start < 0 or end > chroms[chrom]:
                    skipped_boundary += 1
                    continue
                values = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
                if values.size != offsets.size:
                    skipped_boundary += 1
                    continue
                blacklist_mask = np.zeros(values.size, dtype=bool)
                if blacklist is not None:
                    blacklist_mask = ~blacklist.valid_mask(chrom, start, end)
                    values[blacklist_mask] = np.nan
                if gene.strand == "-":
                    values = values[::-1]
                    blacklist_mask = blacklist_mask[::-1]
                finite = np.isfinite(values)
                if missing_to_zero:
                    ordinary_missing = ~finite & ~blacklist_mask
                    values[ordinary_missing] = 0.0
                    usable = ~blacklist_mask
                    total[usable] += values[usable]
                    valid_counts[usable] += 1
                else:
                    total[finite] += values[finite]
                    valid_counts[finite] += 1
                genes_used += 1
                expression_values.append(expr.value)

            mean = np.divide(
                total,
                valid_counts,
                out=np.full(total.shape, np.nan, dtype=float),
                where=valid_counts > 0,
            )
            for offset, value, count in zip(offsets, mean, valid_counts):
                profile_rows.append(
                    {
                        "quintile": key,
                        "quintile_label": label,
                        "relative_position": int(offset),
                        "mean_signal": float(value),
                        "contributing_genes": int(count),
                        "genes_used": genes_used,
                    }
                )
            values_array = np.asarray(expression_values, dtype=float)
            summary_rows.append(
                {
                    "quintile": key,
                    "quintile_label": label,
                    "assigned_gene_count": len(members),
                    "genes_used": genes_used,
                    "skipped_missing_contig": skipped_contig,
                    "skipped_boundary": skipped_boundary,
                    "skipped_blacklisted_anchor": skipped_blacklisted_anchor,
                    "minimum_nTPM": float(np.min(values_array)) if values_array.size else math.nan,
                    "maximum_nTPM": float(np.max(values_array)) if values_array.size else math.nan,
                    "median_nTPM": float(np.median(values_array)) if values_array.size else math.nan,
                    "mean_nTPM": float(np.mean(values_array)) if values_array.size else math.nan,
                }
            )
    return profile_rows, summary_rows


def write_outputs(
    *,
    output_prefix: Path,
    sample: str,
    tissue: str,
    signal_label: str,
    profile_rows: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
    expression_path: Path,
    genes_path: Path,
    duplicate_records: int,
    invalid_records: int,
    window: int,
    blacklist_bed: str | Path | None = None,
) -> dict[str, Path]:
    profile_path = Path(f"{output_prefix}_tss_expression_quintiles.tsv")
    summary_path = Path(f"{output_prefix}_tss_expression_quintile_summary.tsv")
    plot_path = Path(f"{output_prefix}_tss_expression_quintiles.png")
    metadata_path = Path(f"{output_prefix}_tss_expression_quintiles_metadata.tsv")
    profile_path.parent.mkdir(parents=True, exist_ok=True)

    profile_fields = [
        "sample", "tissue", "signal", "quintile", "quintile_label",
        "relative_position", "mean_signal", "contributing_genes", "genes_used",
    ]
    with profile_path.open("wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=profile_fields, delimiter="\t")
        writer.writeheader()
        for row in profile_rows:
            writer.writerow({"sample": sample, "tissue": tissue, "signal": signal_label, **row})

    summary_fields = [
        "sample", "tissue", "signal", "quintile", "quintile_label",
        "assigned_gene_count", "genes_used", "skipped_missing_contig", "skipped_boundary",
        "skipped_blacklisted_anchor",
        "minimum_nTPM", "maximum_nTPM", "median_nTPM", "mean_nTPM",
    ]
    with summary_path.open("wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({"sample": sample, "tissue": tissue, "signal": signal_label, **row})

    by_quintile: dict[str, list[dict[str, object]]] = {}
    for row in profile_rows:
        by_quintile.setdefault(str(row["quintile"]), []).append(row)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for key, label in QUINTILES:
        rows = sorted(by_quintile.get(key, []), key=lambda row: int(row["relative_position"]))
        if not rows:
            continue
        ax.plot(
            [int(row["relative_position"]) for row in rows],
            [float(row["mean_signal"]) for row in rows],
            linewidth=1.4,
            label=label,
        )
    ax.axvline(0, linewidth=0.8, linestyle="--")
    ax.set_xlim(-window, window)
    ax.set_xlabel("Position relative to TSS (bp)")
    ax.set_ylabel(f"Mean {signal_label}")
    ax.set_title(f"{sample}: {signal_label} at TSSs by {tissue} nTPM quintile")
    ax.legend(frameon=False, title="Expression quintile")
    from nucleosuite.plotting import save_figure
    fig.tight_layout()
    plot_path = save_figure(fig, plot_path, default_dpi=200)
    plt.close(fig)

    with metadata_path.open("wt", encoding="utf-8") as handle:
        handle.write("parameter\tvalue\n")
        rows = (
            ("sample", sample),
            ("tissue", tissue),
            ("signal", signal_label),
            ("window", window),
            ("expression", expression_path),
            ("genes_bed", genes_path),
            ("blacklist_bed", blacklist_bed or ""),
            ("expression_duplicate_records_averaged", duplicate_records),
            ("expression_invalid_records_skipped", invalid_records),
            ("profile_tsv", profile_path),
            ("summary_tsv", summary_path),
            ("plot", plot_path),
        )
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")
    return {
        "profile": profile_path,
        "summary": summary_path,
        "plot": plot_path,
        "metadata": metadata_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite tss-expression-quintiles",
        description=(
            "Split genes into five expression quintiles for one tissue/profile and plot a "
            "strand-aware aggregate signal around transcription start sites."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--signal", "--bigwig", required=True, help="Input PNS, WPS, or other BigWig.")
    parser.add_argument("--sample", default="sample", help="Sample label written to outputs.")
    parser.add_argument("--signal-label", default="PNS", help="Signal label used in outputs and plots.")
    parser.add_argument("--expression", help="Expression TSV/TSV.gz. Default: bundled HPA tissue consensus resource.")
    parser.add_argument("--tissue", default="bone_marrow", help="Profile name; underscores are interpreted as spaces.")
    parser.add_argument("--genes-bed", help="Gene BED. Default: genes from --resource-set.")
    parser.add_argument(
        "--blacklist-bed",
        help="BED blacklist; overlapping TSS anchors are skipped and window bases are masked.",
    )
    parser.add_argument("--resource-set", default="hg19-gm12878", help="Bundled resource set supplying genes and tissue expression when explicit files are absent.")
    parser.add_argument("--expression-gene-column", default="Gene", help="Expression-table gene-ID column.")
    parser.add_argument("--expression-name-column", default="Gene name", help="Expression-table gene-name column.")
    parser.add_argument("--expression-profile-column", default="Tissue", help="Expression-table tissue or profile column.")
    parser.add_argument("--expression-value-column", default="nTPM", help="Numeric expression column used for quintile ranking (default: nTPM).")
    parser.add_argument("--gene-id-column", type=int, default=4, help="One-based gene BED column containing the gene identifier.")
    parser.add_argument("--gene-name-column", type=int, default=5, help="One-based gene BED column containing the gene name.")
    parser.add_argument("--gene-strand-column", type=int, default=6, help="One-based gene BED column containing + or - strand.")
    parser.add_argument("--window", type=int, default=2000, help="Bases on each side of the TSS.")
    parser.add_argument("--preserve-missing", action="store_true", help="Exclude missing BigWig bases rather than treating them as zero.")
    parser.add_argument("--output-prefix", required=True, help="Path prefix for quintile profiles, summary, metadata, and plot.")
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    reporter = ProgressReporter("tss-expression-quintiles")
    if args.window < 1:
        raise ValueError("--window must be positive")
    signal_path = Path(args.signal)
    if not signal_path.is_file():
        raise FileNotFoundError(f"Signal BigWig not found: {signal_path}")

    expression_context = None
    gene_context = None
    if args.expression:
        expression_path = Path(args.expression)
        if not expression_path.is_file():
            raise FileNotFoundError(f"Expression table not found: {expression_path}")
    else:
        resource_name = resolve_set_resource_name(args.resource_set, "tissue_expression")
        expression_context = materialized_resource_path(resource_name)
        expression_path = expression_context.__enter__()

    if args.genes_bed:
        gene_path = Path(args.genes_bed)
        if not gene_path.is_file():
            raise FileNotFoundError(f"Gene BED not found: {gene_path}")
    else:
        resource_name = resolve_set_resource_name(args.resource_set, "genes")
        gene_context = materialized_resource_path(resource_name)
        gene_path = gene_context.__enter__()

    try:
        reporter.stage("Loading genes and expression profile")
        genes = read_genes(
            gene_path,
            gene_id_column=args.gene_id_column,
            gene_name_column=args.gene_name_column,
            strand_column=args.gene_strand_column,
        )
        tissue, expression, duplicates, invalid = read_profile_expression(
            expression_path,
            profile=args.tissue,
            gene_column=args.expression_gene_column,
            name_column=args.expression_name_column,
            profile_column=args.expression_profile_column,
            value_column=args.expression_value_column,
        )
        quintiles = assign_quintiles(genes, expression)
        reporter.stage(
            f"Aggregating signal for {sum(len(members) for _, _, members in quintiles):,} "
            "genes across five expression quintiles"
        )
        profile_rows, summary_rows = aggregate_quintiles(
            signal_path,
            quintiles,
            window=args.window,
            missing_to_zero=not args.preserve_missing,
            blacklist_bed=args.blacklist_bed,
        )
        reporter.stage("Writing quintile profiles, summary, metadata, and plot")
        outputs = write_outputs(
            output_prefix=Path(args.output_prefix),
            sample=args.sample,
            tissue=tissue,
            signal_label=args.signal_label,
            profile_rows=profile_rows,
            summary_rows=summary_rows,
            expression_path=Path(expression_path),
            genes_path=Path(gene_path),
            duplicate_records=duplicates,
            invalid_records=invalid,
            window=args.window,
            blacklist_bed=args.blacklist_bed,
        )
        for name, path in outputs.items():
            print(f"{name}\t{path}")
        return 0
    finally:
        if gene_context is not None:
            gene_context.__exit__(None, None, None)
        if expression_context is not None:
            expression_context.__exit__(None, None, None)


if __name__ == "__main__":
    raise SystemExit(main())
