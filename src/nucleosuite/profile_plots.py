"""Small plotting helpers for tabular NucleoSuite profiles."""

from __future__ import annotations

from pathlib import Path
import csv
import math
import re


def plot_dinucleotide_profile(tsv_path: str, png_path: str, title: str | None = None) -> Path | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    with open(tsv_path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return
    position = [int(row["position"]) for row in rows]
    columns = [
        name for name in rows[0]
        if re.fullmatch(r"[ACGT]{2}_(?:pct|frac)", name)
    ]
    if not columns:
        columns = [
            name for name in rows[0]
            if name.endswith(("_pct", "_frac"))
            and not name.startswith(("WW_", "SS_"))
        ]
    fig, ax = plt.subplots(figsize=(12, 5))
    for column in columns:
        values = []
        for row in rows:
            try:
                values.append(float(row[column]))
            except (ValueError, TypeError):
                values.append(float("nan"))
        ax.plot(position, values, linewidth=1.1, label=column.rsplit("_", 1)[0], marker="o", markersize=1.8, markeredgewidth=0)
    ax.axvline(0, color="black", linewidth=1.0, alpha=1.0, zorder=3)
    from nucleosuite.plotting import apply_dyad_profile_x_axis
    apply_dyad_profile_x_axis(ax, position)
    ax.set_xlabel("Position relative to dyad (bp)")
    ax.set_ylabel("Dinucleotide fraction" if columns and columns[0].endswith("_frac") else "Dinucleotide percentage")
    from nucleosuite.plotting import automatic_plot_title
    ax.set_title(automatic_plot_title(tsv_path, title or "Dinucleotide profile"))
    if columns:
        ax.legend(frameon=False, ncol=4, fontsize="small")
    fig.tight_layout()
    from nucleosuite.plotting import save_figure
    saved_path = save_figure(fig, png_path, default_dpi=220)
    from nucleosuite.plotting import write_plot_metadata
    write_plot_metadata(
        saved_path,
        extra={"source_table": str(tsv_path), "resolved_title": ax.get_title()},
    )
    plt.close(fig)
    return saved_path


def plot_ww_ss_profile(tsv_path: str, png_path: str, title: str | None = None) -> Path | None:
    """Plot aggregate WW and SS profiles from a dinucleotide-profile TSV."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    with open(tsv_path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return
    position = [int(row["position"]) for row in rows]
    columns = [name for name in ("WW_pct", "SS_pct", "WW_frac", "SS_frac") if name in rows[0]]
    if not columns:
        return
    fig, ax = plt.subplots(figsize=(12, 5))
    for column in columns:
        values = []
        for row in rows:
            try:
                values.append(float(row[column]))
            except (ValueError, TypeError):
                values.append(float("nan"))
        ax.plot(position, values, linewidth=1.4, label=column.rsplit("_", 1)[0], marker="o", markersize=2.0, markeredgewidth=0)
    ax.axvline(0, color="black", linewidth=1.0, alpha=1.0, zorder=3)
    from nucleosuite.plotting import apply_dyad_profile_x_axis
    apply_dyad_profile_x_axis(ax, position)
    ax.set_xlabel("Position relative to dyad (bp)")
    ax.set_ylabel("Dinucleotide fraction" if columns[0].endswith("_frac") else "Dinucleotide percentage")
    from nucleosuite.plotting import automatic_plot_title
    ax.set_title(automatic_plot_title(tsv_path, title or "WW/SS dinucleotide profile"))
    ax.legend(frameon=False)
    fig.tight_layout()
    from nucleosuite.plotting import save_figure
    saved_path = save_figure(fig, png_path, default_dpi=220)
    plt.close(fig)
    return saved_path


def plot_category_counts(tsv_path: str, png_path: str, title: str | None = None) -> Path | None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()
    with open(tsv_path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return
    category_key = next(iter(rows[0]))
    numeric_keys = [k for k in rows[0] if k != category_key]
    count_key = next((k for k in numeric_keys if "count" in k.lower()), numeric_keys[0] if numeric_keys else None)
    if count_key is None:
        return
    labels=[]; values=[]
    for row in rows:
        try:
            value=float(row[count_key])
        except (ValueError, TypeError):
            continue
        labels.append(row[category_key]); values.append(value)
    if not values:
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    from nucleosuite.plotting import category_colors
    ax.bar(labels, values, color=category_colors(len(labels)))
    if "count" in count_key.lower():
        from nucleosuite.plotting import apply_integer_y_axis
        apply_integer_y_axis(ax)
    ax.set_xlabel(category_key.replace("_", " ").title())
    ax.set_ylabel(count_key.replace("_", " ").title())
    if title:
        ax.set_title(title)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    from nucleosuite.plotting import save_figure
    saved_path = save_figure(fig, png_path, default_dpi=220)
    plt.close(fig)
    return saved_path


def plot_count_profile(
    tsv_path: str,
    png_path: str,
    *,
    x_column: str,
    y_column: str = "count",
    xlabel: str | None = None,
    ylabel: str = "Count",
    title: str | None = None,
    vertical_zero: bool = False,
) -> Path | None:
    """Plot a simple numeric count profile from a TSV.

    This helper is intentionally small and is used for tabular outputs whose
    natural visualization is a one-dimensional frequency curve, such as
    fragment lengths and randomization relocation distances.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    with open(tsv_path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    x_values: list[float] = []
    y_values: list[float] = []
    for row in rows:
        try:
            x_value = float(row[x_column])
            y_value = float(row[y_column])
        except (KeyError, TypeError, ValueError):
            continue
        if not (math.isfinite(x_value) and math.isfinite(y_value)):
            continue
        x_values.append(x_value)
        y_values.append(y_value)
    if not x_values:
        return
    order = sorted(range(len(x_values)), key=x_values.__getitem__)
    x_values = [x_values[index] for index in order]
    y_values = [y_values[index] for index in order]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x_values, y_values, linewidth=1.3, marker="o", markersize=2.2, markeredgewidth=0)
    from nucleosuite.plotting import apply_base_pair_x_axis, apply_integer_y_axis
    apply_base_pair_x_axis(ax, x_values)
    if "count" in y_column.lower() or ylabel.strip().lower() == "count":
        apply_integer_y_axis(ax)
    if vertical_zero:
        ax.axvline(0, linewidth=0.8, alpha=0.5)
    ax.set_xlabel(xlabel or x_column.replace("_", " ").title())
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    from nucleosuite.plotting import save_figure
    saved_path = save_figure(fig, png_path, default_dpi=220)
    plt.close(fig)
    return saved_path


def plot_profile_overlay(
    inputs: list[tuple[str, str | Path]],
    output_tsv: str | Path,
    output_png: str | Path,
    *,
    xlabel: str = "Position relative to feature (bp)",
    ylabel: str = "Mean signal",
    title: str | None = None,
) -> Path | None:
    """Combine NucleoSuite two-column aggregate profiles into one TSV and PNG.

    Each input must contain ``relative_position`` and ``score`` columns. All
    profiles must cover the same coordinates so the resulting overlay is a
    direct point-for-point comparison.
    """

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    if not inputs:
        raise ValueError("At least one profile is required")

    positions: list[int] | None = None
    profiles: list[tuple[str, list[float]]] = []
    seen_labels: set[str] = set()
    for label, path_value in inputs:
        clean_label = str(label).strip()
        if not clean_label:
            raise ValueError("Profile labels must not be empty")
        if clean_label in seen_labels:
            raise ValueError(f"Duplicate profile label: {clean_label}")
        seen_labels.add(clean_label)
        path = Path(path_value)
        with path.open("rt", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            raise ValueError(f"Profile contains no rows: {path}")
        try:
            current_positions = [int(row["relative_position"]) for row in rows]
            current_scores = [float(row["score"]) for row in rows]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Profile must contain numeric relative_position and score columns: {path}"
            ) from exc
        if positions is None:
            positions = current_positions
        elif current_positions != positions:
            raise ValueError(
                f"Profile coordinates differ from the first input: {path}"
            )
        profiles.append((clean_label, current_scores))

    assert positions is not None
    tsv_path = Path(output_tsv)
    png_path = Path(output_png)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)

    with tsv_path.open("wt", encoding="utf-8") as output:
        output.write("relative_position\t" + "\t".join(label for label, _ in profiles) + "\n")
        for row_index, position in enumerate(positions):
            values = [f"{scores[row_index]:.10g}" for _, scores in profiles]
            output.write(f"{position}\t" + "\t".join(values) + "\n")

    figure, axis = plt.subplots(figsize=(12, 5))
    for label, scores in profiles:
        axis.plot(
            positions,
            scores,
            linewidth=1.4,
            marker="o",
            markersize=1.8,
            markeredgewidth=0,
            label=label,
        )
    axis.axvline(0, linewidth=0.8, alpha=0.5)
    from nucleosuite.plotting import apply_base_pair_x_axis

    apply_base_pair_x_axis(axis, positions)
    axis.set_xlim(float(positions[0]), float(positions[-1]))
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    if title:
        axis.set_title(title)
    axis.legend(frameon=False)
    figure.tight_layout()
    from nucleosuite.plotting import save_figure
    saved_path = save_figure(figure, png_path, default_dpi=300)
    plt.close(figure)
    return saved_path


def plot_ww_type_length_stacked(
    tsv_path: str | Path,
    png_path: str | Path,
    *,
    title: str | None = None,
) -> Path | None:
    """Plot type1-type4 relative frequencies as stacked bars by fragment length.

    Bars are normalized across classified fragments so the four WW/SS classes
    sum to 100% at each fragment length. Unclassified counts remain available
    in the TSV but are not drawn as a dinucleotide type.
    """

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_unique_category_cycle
    configure_unique_category_cycle()

    from nucleosuite.sequence.ww_types import WW_TYPE_GROUPS

    with open(tsv_path, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return

    lengths: list[int] = []
    values_by_group: dict[str, list[float]] = {group: [] for group in WW_TYPE_GROUPS}
    for row in rows:
        try:
            fragment_length = int(float(row["fragment_length"]))
        except (KeyError, TypeError, ValueError):
            continue
        current: dict[str, float] = {}
        valid = True
        for group in WW_TYPE_GROUPS:
            try:
                value = float(row[f"{group}_percent_of_classified"])
            except (KeyError, TypeError, ValueError):
                valid = False
                break
            if not math.isfinite(value):
                value = 0.0
            current[group] = value
        if not valid:
            continue
        lengths.append(fragment_length)
        for group in WW_TYPE_GROUPS:
            values_by_group[group].append(current[group])
    if not lengths:
        return

    order = sorted(range(len(lengths)), key=lengths.__getitem__)
    lengths = [lengths[index] for index in order]
    values_by_group = {
        group: [values[index] for index in order]
        for group, values in values_by_group.items()
    }

    positions = list(range(len(lengths)))
    bottom = [0.0] * len(lengths)
    figure, axis = plt.subplots(figsize=(max(7.0, 1.1 * len(lengths) + 4.0), 5.5))
    for group in WW_TYPE_GROUPS:
        values = values_by_group[group]
        axis.bar(positions, values, bottom=bottom, label=group)
        bottom = [current + value for current, value in zip(bottom, values)]

    axis.set_xticks(positions, [str(length) for length in lengths])
    axis.set_xlabel("Fragment length (bp)")
    axis.set_ylabel("Relative frequency among classified fragments (%)")
    axis.set_ylim(0.0, 100.0)
    if title:
        axis.set_title(title)
    axis.legend(frameon=False, ncol=4)
    axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    from nucleosuite.plotting import save_figure
    saved_path = save_figure(figure, png_path, default_dpi=300)
    plt.close(figure)
    return saved_path
