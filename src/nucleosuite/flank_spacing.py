"""Category-wise spacing between nucleosomes flanking reference sites."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.output_naming import parameter_range, parameterized_prefix
from nucleosuite.plotting import add_plotting_arguments, category_colors, save_figure


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def _read_bed_rows(path: Path, *, skip_header: bool = False) -> Iterable[list[str]]:
    with _open_text(path) as handle:
        skipped = False
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            if skip_header and not skipped:
                skipped = True
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            yield fields


def _coordinate(fields: Sequence[str], column: int | None) -> int:
    if column is None:
        return (int(fields[1]) + int(fields[2])) // 2
    index = column - 1
    if index < 0 or index >= len(fields):
        raise ValueError(f"Requested BED column {column} is not present in a row with {len(fields)} columns")
    return int(float(fields[index]))


def load_nucleosome_centres(path: str | Path, *, centre_col: int | None = None) -> dict[str, list[int]]:
    centres: dict[str, list[int]] = defaultdict(list)
    for fields in _read_bed_rows(Path(path)):
        centres[fields[0]].append(_coordinate(fields, centre_col))
    if not centres:
        raise ValueError("No nucleosome calls were read from --nucleosome-bed")
    for values in centres.values():
        values.sort()
    return dict(centres)


def compute_flanking_spacings(
    region_bed: str | Path,
    nucleosome_centres: dict[str, list[int]],
    *,
    category_col: int = 4,
    point_col: int | None = None,
    skip_header: bool = False,
) -> tuple[list[dict[str, object]], dict[str, list[int]], dict[str, int]]:
    rows: list[dict[str, object]] = []
    by_category: dict[str, list[int]] = defaultdict(list)
    totals: dict[str, int] = Counter()
    references = list(nucleosome_centres)
    for fields in _read_bed_rows(Path(region_bed), skip_header=skip_header):
        if category_col < 1 or category_col > len(fields):
            raise ValueError(f"--category-col {category_col} is not present in a row with {len(fields)} columns")
        category = fields[category_col - 1]
        totals[category] += 1
        centre = _coordinate(fields, point_col)
        try:
            nuc_chrom = resolve_contig_name(fields[0], references, source_label="nucleosome BED")
        except KeyError:
            rows.append({
                "chrom": fields[0], "region_start": int(fields[1]), "region_end": int(fields[2]),
                "category": category, "region_center": centre,
                "upstream_nucleosome_center": "", "downstream_nucleosome_center": "",
                "flanking_spacing_bp": "", "status": "no_matching_contig",
            })
            continue
        values = nucleosome_centres[nuc_chrom]
        left_i = bisect_left(values, centre) - 1
        right_i = bisect_right(values, centre)
        if left_i < 0 or right_i >= len(values):
            rows.append({
                "chrom": fields[0], "region_start": int(fields[1]), "region_end": int(fields[2]),
                "category": category, "region_center": centre,
                "upstream_nucleosome_center": values[left_i] if left_i >= 0 else "",
                "downstream_nucleosome_center": values[right_i] if right_i < len(values) else "",
                "flanking_spacing_bp": "", "status": "missing_flank",
            })
            continue
        upstream = values[left_i]
        downstream = values[right_i]
        spacing = downstream - upstream
        by_category[category].append(spacing)
        rows.append({
            "chrom": fields[0], "region_start": int(fields[1]), "region_end": int(fields[2]),
            "category": category, "region_center": centre,
            "upstream_nucleosome_center": upstream, "downstream_nucleosome_center": downstream,
            "flanking_spacing_bp": spacing, "status": "ok",
        })
    if not totals:
        raise ValueError("No reference regions were read from --region-bed")
    return rows, dict(by_category), dict(totals)


def distribution_curve(values: Sequence[int], x: np.ndarray, mode: str) -> np.ndarray:
    data = np.asarray(values, dtype=float)
    if mode == "count":
        counts = Counter(int(v) for v in values)
        return np.asarray([counts.get(int(round(value)), 0) for value in x], dtype=float)
    if data.size < 2:
        return np.zeros_like(x, dtype=float)
    if np.allclose(data, data[0]):
        y = np.zeros_like(x, dtype=float)
        idx = int(np.argmin(np.abs(x - data[0])))
        if y.size:
            y[idx] = 1.0
        return y
    from scipy.stats import gaussian_kde
    return np.asarray(gaussian_kde(data)(x), dtype=float)


def _ratio(y1: float, y2: float) -> float:
    if not math.isfinite(y1) or not math.isfinite(y2):
        return float("nan")
    if y2 == 0:
        return float("inf") if y1 > 0 else float("nan")
    return y1 / y2


def rank_categories(
    by_category: dict[str, list[int]],
    *,
    mode: str,
    ratio_x1: int,
    ratio_x2: int,
    x_grid: np.ndarray,
) -> tuple[list[dict[str, object]], dict[str, np.ndarray]]:
    curves: dict[str, np.ndarray] = {}
    rows: list[dict[str, object]] = []
    for category in sorted(by_category):
        curve = distribution_curve(by_category[category], x_grid, mode)
        curves[category] = curve
        i1 = int(np.argmin(np.abs(x_grid - ratio_x1)))
        i2 = int(np.argmin(np.abs(x_grid - ratio_x2)))
        y1, y2 = float(curve[i1]), float(curve[i2])
        ratio = _ratio(y1, y2)
        rows.append({
            "category": category,
            "valid_flanking_pair_count": len(by_category[category]),
            f"y_at_{ratio_x1}": y1,
            f"y_at_{ratio_x2}": y2,
            f"ratio_{ratio_x1}_to_{ratio_x2}": ratio,
        })
    def key(row):
        value = float(row[f"ratio_{ratio_x1}_to_{ratio_x2}"])
        if math.isnan(value):
            return (2, float("inf"), str(row["category"]))
        if math.isinf(value):
            return (1, value, str(row["category"]))
        return (0, value, str(row["category"]))
    rows.sort(key=key)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows, curves


def _write_tsv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _plot(
    curves: dict[str, np.ndarray],
    rankings: Sequence[dict[str, object]],
    x_grid: np.ndarray,
    output: Path,
    *,
    mode: str,
    x_min: float,
    x_max: float,
    top_categories: int,
    ratio_x1: int,
    ratio_x2: int,
):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    selected = [str(row["category"]) for row in rankings[:top_categories]]
    selected_set = set(selected)
    for category in sorted(curves):
        if category not in selected_set:
            ax.plot(x_grid, curves[category], color="0.72", linewidth=1.0, alpha=0.75, zorder=1)
    colours = category_colors(len(selected))
    handle_by_category = {}
    # Draw lower-priority selected categories first; rank 1 is drawn last/on top.
    for reverse_index, category in enumerate(reversed(selected)):
        rank = selected.index(category) + 1
        color = colours[rank - 1]
        line, = ax.plot(
            x_grid, curves[category], color=color, linewidth=1.6,
            label=category, zorder=10 + (len(selected) - rank),
        )
        handle_by_category[category] = line
    if selected:
        ax.legend([handle_by_category[c] for c in selected], selected, frameon=False, title="Category")
    ax.set_xlim(float(x_min), float(x_max))
    ax.set_xlabel("Distance (bp) between flanking nucleosome centres")
    ax.set_ylabel("Density" if mode == "density" else "Count")
    ax.set_title(f"Flanking nucleosome spacing ({ratio_x1}/{ratio_x2} ratio ranking)")
    ax.spines[["top", "right"]].set_visible(False)
    return save_figure(fig, output, default_dpi=220), fig


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite flank-spacing",
        description=(
            "Measure the distance between the nearest upstream and downstream nucleosome "
            "centres around BED reference sites, compare distributions by category, and rank "
            "categories by the ratio of distribution heights at two spacing values."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--nucleosome-bed", required=True, type=Path, help="BED of nucleosome calls.")
    parser.add_argument("--region-bed", required=True, type=Path, help="BED of reference sites/regions.")
    parser.add_argument("--category-col", type=int, default=4, help="One-based category column in the reference BED.")
    parser.add_argument("--point-col", type=int, help="One-based exact reference coordinate column; otherwise interval midpoint.")
    parser.add_argument("--nucleosome-center-col", type=int, help="One-based exact nucleosome-centre column; otherwise interval midpoint.")
    parser.add_argument("--skip-header", action="store_true", help="Skip the first non-comment reference-BED line.")
    parser.add_argument("--distribution", choices=("density", "count"), default="density", help="Plot density curves or raw counts.")
    parser.add_argument("--ratio-x1", type=int, default=190, help="Numerator spacing position for category ranking.")
    parser.add_argument("--ratio-x2", type=int, default=260, help="Denominator spacing position for category ranking.")
    parser.add_argument("--top-categories", type=int, default=7, help="Number of best-ranked categories to colour and label.")
    parser.add_argument("--x-min", type=int, default=0, help="Minimum spacing displayed/evaluated in the distribution table.")
    parser.add_argument("--x-max", type=int, default=500, help="Maximum spacing displayed/evaluated in the distribution table.")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Output directory.")
    parser.add_argument("--output-prefix", help="Base output prefix; automatic parameter tokens are appended.")
    parser.add_argument(
        "--write-detail-tables", action="store_true",
        help="Write the per-reference-site flanking-pair table; omitted by default.",
    )
    add_plotting_arguments(parser)
    return parser


def validate_argv(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    for name in ("category_col", "point_col", "nucleosome_center_col"):
        value = getattr(args, name)
        if value is not None and value < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be at least 1")
    if args.top_categories < 0:
        raise ValueError("--top-categories must be zero or greater")
    if args.x_min < 0 or args.x_max <= args.x_min:
        raise ValueError("--x-max must be greater than --x-min and --x-min must be non-negative")
    for value, option in ((args.ratio_x1, "--ratio-x1"), (args.ratio_x2, "--ratio-x2")):
        if value < args.x_min or value > args.x_max:
            raise ValueError(f"{option} must lie between --x-min and --x-max")


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    validate_argv(argv)
    centres = load_nucleosome_centres(args.nucleosome_bed, centre_col=args.nucleosome_center_col)
    detail, by_category, totals = compute_flanking_spacings(
        args.region_bed, centres, category_col=args.category_col,
        point_col=args.point_col, skip_header=args.skip_header,
    )
    if not by_category:
        raise RuntimeError("No reference regions had both an upstream and downstream nucleosome call")
    x_grid = np.arange(args.x_min, args.x_max + 1, dtype=float)
    rankings, curves = rank_categories(
        by_category, mode=args.distribution, ratio_x1=args.ratio_x1,
        ratio_x2=args.ratio_x2, x_grid=x_grid,
    )
    for row in rankings:
        row["reference_site_count"] = totals.get(str(row["category"]), 0)

    base_name = args.output_prefix or f"{args.region_bed.stem}_flank_spacing"
    base = parameterized_prefix(
        args.output_dir / base_name,
        (("dist", args.distribution), ("ratio", parameter_range(args.ratio_x1, args.ratio_x2)), ("xmax", args.x_max)),
    )
    detail_path = Path(f"{base}_sites.tsv")
    distribution_path = Path(f"{base}_distributions.tsv")
    ranking_path = Path(f"{base}_ranking.tsv")
    plot_path = Path(f"{base}.png")

    if args.write_detail_tables:
        _write_tsv(detail_path, detail, [
            "chrom", "region_start", "region_end", "category", "region_center",
            "upstream_nucleosome_center", "downstream_nucleosome_center", "flanking_spacing_bp", "status",
        ])
    dist_rows = []
    rank_by_category = {str(row["category"]): int(row["rank"]) for row in rankings}
    for category in sorted(curves):
        rank = rank_by_category[category]
        for x, y in zip(x_grid, curves[category]):
            dist_rows.append({
                "category": category,
                "spacing_bp": int(x),
                "value": float(y),
                "distribution": args.distribution,
                "rank": rank,
                "highlighted": int(rank <= args.top_categories),
                "ratio_x1": args.ratio_x1,
                "ratio_x2": args.ratio_x2,
                "top_categories": args.top_categories,
                "x_min": args.x_min,
                "x_max": args.x_max,
            })
    _write_tsv(
        distribution_path, dist_rows,
        ["category", "spacing_bp", "value", "distribution", "rank", "highlighted",
         "ratio_x1", "ratio_x2", "top_categories", "x_min", "x_max"],
    )
    ratio_col = f"ratio_{args.ratio_x1}_to_{args.ratio_x2}"
    _write_tsv(ranking_path, rankings, [
        "rank", "category", "reference_site_count", "valid_flanking_pair_count",
        f"y_at_{args.ratio_x1}", f"y_at_{args.ratio_x2}", ratio_col,
    ])
    saved, fig = _plot(
        curves, rankings, x_grid, plot_path, mode=args.distribution,
        x_min=args.x_min, x_max=args.x_max, top_categories=args.top_categories,
        ratio_x1=args.ratio_x1, ratio_x2=args.ratio_x2,
    )
    from nucleosuite.plotting import write_plot_metadata
    write_plot_metadata(saved, extra={"source_table": str(distribution_path), "detected_plot_type": "flank-spacing"})
    import matplotlib.pyplot as plt
    plt.close(fig)
    if args.write_detail_tables:
        print(f"Sites: {detail_path}")
    print(f"Distributions: {distribution_path}")
    print(f"Ranking: {ranking_path}")
    print(f"Plot: {saved}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
