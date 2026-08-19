"""Category-aware orchestration for :mod:`nucleosuite aggregate`.

The core aggregate implementation remains deliberately category-agnostic. This
module partitions a BED file by a selected column, runs the ordinary aggregate
pipeline independently for every category, then writes an overlay profile and a
combined NRL summary. Reusing the ordinary aggregate runner keeps all filtering,
strand orientation, blacklist handling, per-contig combination and NRL behavior
identical to a standalone aggregate run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import tempfile
from dataclasses import replace
from pathlib import Path
from typing import Callable

from nucleosuite.align import AlignmentConfig, make_output_prefix, resolve_output_paths
from nucleosuite.profile_plots import plot_profile_overlay


_EXPLICIT_SINGLE_OUTPUTS = (
    "heatmap_output",
    "heatmap_matrix_output",
    "aggregate_output",
    "plotted_mean_output",
    "mean_plot_output",
    "summary_output",
)


def _safe_category_name(value: str) -> str:
    """Return a filesystem-safe category token while preserving readable names."""
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    return token or "category"


def _unique_category_names(categories: list[str]) -> dict[str, str]:
    """Map category labels to deterministic collision-free filesystem tokens."""
    used: set[str] = set()
    mapping: dict[str, str] = {}
    for category in categories:
        base = _safe_category_name(category)
        candidate = base
        suffix = 2
        while candidate.casefold() in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used.add(candidate.casefold())
        mapping[category] = candidate
    return mapping


def _open_region_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).lower().endswith(".gz") else path.open("rt", encoding="utf-8")


def split_region_bed_by_category(
    region_bed: Path,
    *,
    category_col: int,
    skip_header: bool,
    destination: Path,
) -> tuple[list[str], dict[str, Path], dict[str, int]]:
    """Split a BED/TSV-like region file by a one-based category column."""
    if category_col < 1:
        raise ValueError("category_col must be >= 1")
    destination.mkdir(parents=True, exist_ok=True)

    records: dict[str, list[str]] = {}
    with _open_region_text(Path(region_bed)) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if skip_header and line_number == 1:
                continue
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < category_col:
                raise ValueError(
                    f"Region line {line_number} has {len(fields)} columns but "
                    f"--category-col {category_col} was requested"
                )
            category = fields[category_col - 1].strip()
            if not category:
                raise ValueError(
                    f"Region line {line_number} has an empty category in column {category_col}"
                )
            records.setdefault(category, []).append(raw if raw.endswith("\n") else raw + "\n")

    if not records:
        raise ValueError("No category-bearing region records were found")

    categories = sorted(records, key=lambda value: (value.casefold(), value))
    safe_names = _unique_category_names(categories)
    paths: dict[str, Path] = {}
    counts: dict[str, int] = {}
    for category in categories:
        path = destination / f"{safe_names[category]}.bed"
        with path.open("wt", encoding="utf-8") as output:
            output.writelines(records[category])
        paths[category] = path
        counts[category] = len(records[category])
    return categories, paths, counts


def _base_output_name(args: argparse.Namespace) -> str:
    """Return the unparameterized user-facing base name for category outputs."""
    if getattr(args, "output_prefix", None):
        return str(args.output_prefix)
    from nucleosuite.io import strip_known_suffix

    parts = [strip_known_suffix(Path(args.region_bed)), strip_known_suffix(Path(args.bigwig))]
    if getattr(args, "nucleosome_bed", None) is not None:
        offset = int(getattr(args, "nucleosome_offset", 1))
        label = f"plus{offset}" if offset > 0 else f"minus{abs(offset)}"
        parts.extend([strip_known_suffix(Path(args.nucleosome_bed)), label])
    if getattr(args, "state_bed", None) is not None:
        parts.append(strip_known_suffix(Path(args.state_bed)))
    return "_".join(str(part) for part in parts)


def _category_overall_prefix(args: argparse.Namespace) -> str:
    """Build a parameterized prefix for combined category outputs."""
    from nucleosuite.output_naming import parameter_range, parameterized_prefix
    from nucleosuite.align import resolve_nrl_exclusion

    parameters: list[tuple[str, object]] = [
        ("catcol", int(args.category_col)),
        ("win", int(args.window_half)),
        ("zero", int(args.zero_thresh)),
        ("maxscore", args.max_score),
        ("missing", "zero" if bool(args.nan_to_zero) else "reject"),
        ("sort", str(args.sort_mode)),
    ]
    if bool(args.nrl):
        config_values = {
            key: value
            for key, value in vars(args).items()
            if key in AlignmentConfig.__dataclass_fields__
        }
        config = AlignmentConfig(**config_values)
        exclusion_start, exclusion_end = resolve_nrl_exclusion(config)
        regression_max = config.window_half if config.nrl_regression_max is None else config.nrl_regression_max
        parameters.extend(
            [
                ("nrlres", config.nrl_peak_resolution),
                ("nrlmin", config.nrl_regression_min),
                ("nrlmax", regression_max),
                (
                    "excl",
                    "none" if exclusion_start is None else parameter_range(exclusion_start, exclusion_end),
                ),
            ]
        )
    else:
        parameters.append(("nrl", "off"))
    return str(parameterized_prefix(_base_output_name(args), parameters))


def _category_result_paths(args: argparse.Namespace) -> dict[str, Path]:
    """Resolve final aggregate paths for one category, including multicontig mode."""
    from nucleosuite.aggregate_parallel import _resolve_contigs

    config_values = {
        key: value
        for key, value in vars(args).items()
        if key in AlignmentConfig.__dataclass_fields__
    }
    config = AlignmentConfig(**config_values)
    cores = int(getattr(args, "cores", 1) or 1)
    final_dir = Path(config.output_dir)
    if cores > 1:
        contigs, _sizes = _resolve_contigs(args)
        if len(contigs) > 1:
            prefix = make_output_prefix(config)
            root = (
                Path(args.parallel_dir)
                if getattr(args, "parallel_dir", None)
                else Path(config.output_dir) / f"{prefix}_multicontig"
            ).resolve()
            final_dir = root / "combined"
    return resolve_output_paths(replace(config, output_dir=final_dir))


def _read_valid_count(summary_path: Path) -> int:
    if not summary_path.is_file():
        return 0
    with summary_path.open("rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("section") == "statistics" and row.get("key") == "valid_total":
                try:
                    return int(float(row.get("value", "0") or 0))
                except ValueError:
                    return 0
    return 0


def write_category_nrl_summary(
    category_outputs: list[tuple[str, dict[str, Path]]],
    output_path: Path,
) -> Path:
    """Combine the two-direction NRL summaries from every category."""
    rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None
    for category, outputs in category_outputs:
        nrl_path = outputs.get("nrl_summary")
        if nrl_path is None or not Path(nrl_path).is_file():
            continue
        valid_count = _read_valid_count(Path(outputs["summary"]))
        with Path(nrl_path).open("rt", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            if fieldnames is None:
                fieldnames = list(reader.fieldnames)
            for row in reader:
                rows.append(
                    {
                        "category": category,
                        "valid_region_count": str(valid_count),
                        **{name: row.get(name, "") for name in reader.fieldnames},
                    }
                )

    if fieldnames is None:
        raise RuntimeError("No per-category NRL summaries were available to combine")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["category", "valid_region_count", *fieldnames],
        )
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def run_category_aggregate(
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
) -> int:
    """Run ordinary aggregate independently for each category in ``args.region_bed``."""
    if int(getattr(args, "category_col", 0) or 0) < 1:
        raise ValueError("--category-col must be at least 1 in category mode")

    conflicting = [
        f"--{name.replace('_', '-')}"
        for name in _EXPLICIT_SINGLE_OUTPUTS
        if getattr(args, name, None) is not None
    ]
    if conflicting:
        raise ValueError(
            "Category aggregation produces multiple per-category outputs, so these "
            "single-output overrides cannot be used with --category-col: "
            + ", ".join(conflicting)
        )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    overall_prefix = _category_overall_prefix(args)
    combined_profile = Path(
        getattr(args, "category_profile_output", None)
        or output_dir / f"{overall_prefix}_category_profiles.tsv"
    )
    from nucleosuite.plotting import plot_path

    combined_plot = plot_path(
        getattr(args, "category_plot_output", None)
        or output_dir / f"{overall_prefix}_category_profiles.png"
    )
    combined_nrl = Path(
        getattr(args, "category_nrl_summary_output", None)
        or output_dir / f"{overall_prefix}_category_nrl_summary.tsv"
    )

    original_region = Path(args.region_bed)
    category_outputs: list[tuple[str, dict[str, Path]]] = []
    profile_specs: list[tuple[str, Path]] = []

    with tempfile.TemporaryDirectory(prefix=".aggregate_categories_", dir=output_dir) as temp_name:
        categories, category_beds, input_counts = split_region_bed_by_category(
            original_region,
            category_col=int(args.category_col),
            skip_header=bool(args.skip_header),
            destination=Path(temp_name),
        )
        safe_names = _unique_category_names(categories)
        base_name = _base_output_name(args)

        for category in categories:
            safe = safe_names[category]
            category_dir = output_dir / safe
            category_args = argparse.Namespace(**vars(args))
            category_args.category_col = 0
            category_args.region_bed = category_beds[category]
            category_args.skip_header = False
            category_args.output_dir = category_dir
            category_args.output_prefix = f"{base_name}_{safe}"
            category_args.category_profile_output = None
            category_args.category_plot_output = None
            category_args.category_nrl_summary_output = None
            if getattr(args, "parallel_dir", None) is not None:
                category_args.parallel_dir = Path(args.parallel_dir) / safe

            from nucleosuite.aggregate_parallel import run_aggregate_per_contig

            exit_code = int(run_aggregate_per_contig(category_args, serial_runner) or 0)
            if exit_code:
                raise RuntimeError(f"Aggregate failed for category {category!r} with exit code {exit_code}")
            outputs = _category_result_paths(category_args)
            aggregate_path = Path(outputs["aggregate"])
            if not aggregate_path.is_file():
                raise RuntimeError(
                    f"Aggregate output was not produced for category {category!r}: {aggregate_path}"
                )
            profile_specs.append((category, aggregate_path))
            category_outputs.append((category, outputs))
            valid_count = _read_valid_count(Path(outputs["summary"]))
            print(
                f"Category {category}: input regions={input_counts[category]:,}; "
                f"valid aggregate regions={valid_count:,}"
            )

    plot_profile_overlay(
        profile_specs,
        combined_profile,
        combined_plot,
        xlabel=f"{args.axis_label} (bp)",
        ylabel=str(args.mean_ylabel),
        title=f"{_base_output_name(args)}: aggregate profiles by category",
    )
    print(f"Combined category profile: {combined_profile}")
    print(f"Combined category plot: {combined_plot}")

    if bool(args.nrl):
        write_category_nrl_summary(category_outputs, combined_nrl)
        print(f"Combined category NRL summary: {combined_nrl}")
    return 0
