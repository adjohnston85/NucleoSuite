"""Per-contig execution and exact recombination for ``nucleosuite aggregate``."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np

from nucleosuite.align import (
    AlignmentConfig,
    AlignmentStats,
    add_heatmap_row,
    central_crop,
    make_output_prefix,
    no_valid_regions_message,
    plot_mean_profile,
    plot_outputs,
    resolve_output_paths,
    resolve_nrl_exclusion,
    sort_matrix,
    write_heatmap_matrix,
    write_heatmap_row_metadata,
    write_profile,
    write_summary,
)
from nucleosuite.parallel import MANIFEST_NAME, _safe_contig, _worker_initializer
from nucleosuite.core.regions import resolve_contig_name


def _stable_seed(seed: int | None, contig: str) -> int | None:
    if seed is None:
        return None
    digest = hashlib.blake2b(contig.encode("utf-8"), digest_size=4).digest()
    return (int(seed) + int.from_bytes(digest, "big")) % (2**32)


def _filter_bed(
    source: Path,
    destination: Path,
    contig: str,
    *,
    chrom_col: int,
    skip_header: bool,
) -> None:
    opener = gzip.open if str(source).lower().endswith(".gz") else open
    destination.parent.mkdir(parents=True, exist_ok=True)
    with opener(source, "rt", encoding="utf-8") as input_handle, destination.open("w", encoding="utf-8") as output:
        for line_number, raw in enumerate(input_handle, start=1):
            if skip_header and line_number == 1:
                output.write(raw)
                continue
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < chrom_col:
                continue
            try:
                matched = resolve_contig_name(
                    fields[chrom_col - 1], [contig], source_label="aggregate contig"
                )
            except KeyError:
                continue
            if matched == contig:
                output.write(raw if raw.endswith("\n") else raw + "\n")


def _region_contigs(region_bed: Path, chrom_col: int, skip_header: bool) -> list[str]:
    opener = gzip.open if str(region_bed).lower().endswith(".gz") else open
    seen: set[str] = set()
    contigs: list[str] = []
    with opener(region_bed, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if skip_header and line_number == 1:
                continue
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < chrom_col:
                continue
            contig = fields[chrom_col - 1]
            if contig not in seen:
                seen.add(contig)
                contigs.append(contig)
    return contigs


def _resolve_contigs(args: argparse.Namespace) -> tuple[list[str], list[tuple[str, int]]]:
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("aggregate multicontig processing requires pyBigWig") from exc
    from nucleosuite.core.regions import expand_contig_tokens

    handle = pyBigWig.open(str(args.bigwig))
    if handle is None:
        raise OSError(f"Could not open BigWig: {args.bigwig}")
    try:
        sizes = {str(chrom): int(size) for chrom, size in handle.chroms().items()}
    finally:
        handle.close()
    present = _region_contigs(Path(args.region_bed), int(args.chrom_col), bool(args.skip_header))
    available: list[str] = []
    for chrom in sizes:
        for region_chrom in present:
            try:
                if resolve_contig_name(
                    region_chrom, [chrom], source_label="aggregate BigWig"
                ) == chrom:
                    available.append(chrom)
                    break
            except KeyError:
                continue
    tokens = getattr(args, "contigs", None) or ["all"]
    selected_specs = expand_contig_tokens(tokens, available)
    selected = []
    seen: set[str] = set()
    for spec in selected_specs:
        chrom = spec.split(":", 1)[0]
        if chrom not in seen:
            seen.add(chrom)
            selected.append(chrom)
    return selected, [(chrom, sizes[chrom]) for chrom in selected]


def _worker(
    namespace_data: dict,
    contig: str,
    contig_dir: str,
    prefix: str,
    region_bed: str,
    seed: int | None,
) -> tuple[int, str, str]:
    from nucleosuite.cli.aggregate import run_from_args

    data = dict(namespace_data)
    requested_details = bool(data.get("write_detail_tables", False))
    requested_stop = data.get("stop_after_valid")
    needs_matrix_for_combine = requested_details or requested_stop is not None
    if requested_stop is not None:
        # The combiner needs at most the first K rows from each ordered contig
        # to reproduce a global --stop-after-valid K without scanning or retaining
        # the remainder of any contig.
        worker_max_rows = None
        worker_subsample_mode = "first"
    else:
        worker_max_rows = data.get("max_heatmap_rows") if requested_details else None
        worker_subsample_mode = (
            data.get("subsample_mode", "first") if requested_details else "first"
        )
    accumulator = Path(contig_dir) / f"{prefix}_aggregate_accumulator.tsv.gz"
    data.update(
        {
            "region_bed": Path(region_bed),
            "output_dir": Path(contig_dir),
            "output_prefix": prefix,
            "heatmap_output": None,
            "heatmap_matrix_output": None,
            "aggregate_output": None,
            "plotted_mean_output": None,
            "mean_plot_output": None,
            "summary_output": None,
            "accumulator_output": accumulator,
            "write_detail_tables": needs_matrix_for_combine,
            "max_heatmap_rows": worker_max_rows,
            "stop_after_valid": requested_stop,
            "subsample_mode": worker_subsample_mode,
            "breadth": 1.0,
            "sort_mode": "unsorted",
            "seed": seed,
            "cores": 1,
            "parallel_dir": None,
            "skip_combine": True,
            "skip_combined_tracks": True,
            "_per_contig_worker": True,
        }
    )
    return (
        int(run_from_args(argparse.Namespace(**data)) or 0),
        str(Path(contig_dir) / prefix),
        str(accumulator),
    )


def _iter_matrix(path: Path) -> tuple[np.ndarray, Iterator[np.ndarray]]:
    handle = gzip.open(path, "rt", encoding="utf-8") if str(path).lower().endswith(".gz") else path.open("r", encoding="utf-8")
    header = handle.readline().rstrip("\n").split("\t")
    x_values = np.asarray([int(float(value)) for value in header[1:]], dtype=int)

    def rows() -> Iterator[np.ndarray]:
        try:
            for raw in handle:
                fields = raw.rstrip("\n").split("\t")
                if len(fields) != len(header):
                    continue
                yield np.asarray([float(value) for value in fields[1:]], dtype=float)
        finally:
            handle.close()

    return x_values, rows()


def _read_accumulator(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    positions: list[int] = []
    totals: list[float] = []
    counts: list[int] = []
    accepted_regions: int | None = None
    with opener(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "relative_position", "signal_sum", "valid_count", "accepted_regions"
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid aggregate accumulator: {path}")
        for row in reader:
            positions.append(int(float(row["relative_position"])))
            totals.append(float(row["signal_sum"]))
            counts.append(int(float(row["valid_count"])))
            row_total = int(float(row["accepted_regions"]))
            if accepted_regions is None:
                accepted_regions = row_total
            elif accepted_regions != row_total:
                raise ValueError(
                    f"Inconsistent accepted-region total in {path}"
                )
    return (
        np.asarray(positions, dtype=int),
        np.asarray(totals, dtype=float),
        np.asarray(counts, dtype=np.int64),
        int(accepted_regions or 0),
    )


def _read_stats(path: Path) -> AlignmentStats:
    values = AlignmentStats()
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("section") != "statistics":
                continue
            key = row.get("key", "")
            if hasattr(values, key):
                setattr(values, key, int(float(row.get("value", 0) or 0)))
    return values


def _combine(
    args: argparse.Namespace,
    root: Path,
    entries: Sequence[dict[str, object]],
    base_name: str,
    output_dir: Path | None = None,
) -> dict[str, str]:
    combined_dir = output_dir or (root / "combined")
    combined_dir.mkdir(parents=True, exist_ok=True)
    config_values = {
        key: value
        for key, value in vars(args).items()
        if key in AlignmentConfig.__dataclass_fields__
    }
    config_values.update(
        {
            "output_dir": combined_dir,
            "output_prefix": base_name,
            "heatmap_output": None,
            "heatmap_matrix_output": None,
            "aggregate_output": None,
            "plotted_mean_output": None,
            "mean_plot_output": None,
            "summary_output": None,
        }
    )
    config = AlignmentConfig(**config_values)
    outputs = resolve_output_paths(config)

    stats = AlignmentStats()
    selected_rows: list[np.ndarray] = []
    rng = np.random.default_rng(config.seed)
    running_sum: np.ndarray | None = None
    running_count: np.ndarray | None = None
    valid_total = 0
    reference_x: np.ndarray | None = None
    stop = False
    entry_valid_totals: list[int] = []

    for entry in entries:
        prefix = Path(str(entry["prefix"]))
        summary_path = prefix.parent / f"{prefix.name}_summary.tsv"
        if summary_path.exists():
            part_stats = _read_stats(summary_path)
            for field in AlignmentStats.__dataclass_fields__:
                if field in {"valid_total", "selected_for_plot", "stopped_after_valid_limit"}:
                    continue
                setattr(stats, field, getattr(stats, field) + getattr(part_stats, field))

    accumulator_paths = [
        Path(str(entry["accumulator"]))
        for entry in entries
        if entry.get("accumulator")
    ]
    use_accumulators = (
        config.stop_after_valid is None
        and len(accumulator_paths) == len(entries)
        and all(path.is_file() for path in accumulator_paths)
    )

    if use_accumulators:
        for entry, accumulator_path in zip(entries, accumulator_paths):
            x_values, part_sum, part_count, part_valid_total = _read_accumulator(
                accumulator_path
            )
            if reference_x is None:
                reference_x = x_values
                running_sum = np.zeros_like(part_sum)
                running_count = np.zeros_like(part_count)
            elif not np.array_equal(reference_x, x_values):
                raise ValueError(
                    f"Aggregate accumulator positions differ in {accumulator_path}"
                )
            assert running_sum is not None and running_count is not None
            running_sum += part_sum
            running_count += part_count
            valid_total += part_valid_total
            entry_valid_totals.append(part_valid_total)

    if not use_accumulators or config.write_detail_tables:
        random_allocations: list[int] | None = None
        if (
            use_accumulators
            and config.write_detail_tables
            and config.max_heatmap_rows is not None
            and config.subsample_mode == "random"
        ):
            remaining_draw = min(config.max_heatmap_rows, valid_total)
            remaining_population = valid_total
            random_allocations = []
            for population in entry_valid_totals:
                if remaining_draw == 0:
                    selected = 0
                elif population == remaining_population:
                    selected = remaining_draw
                else:
                    selected = int(
                        rng.hypergeometric(
                            population,
                            remaining_population - population,
                            remaining_draw,
                        )
                    )
                random_allocations.append(selected)
                remaining_draw -= selected
                remaining_population -= population

        sampled_index = 0
        for entry_index, entry in enumerate(entries):
            prefix = Path(str(entry["prefix"]))
            matrix_path = prefix.parent / f"{prefix.name}_heatmap_matrix.tsv.gz"
            x_values, rows = _iter_matrix(matrix_path)
            if reference_x is None:
                reference_x = x_values
            elif not np.array_equal(reference_x, x_values):
                raise ValueError(f"Aggregate matrix positions differ in {matrix_path}")
            if random_allocations is not None:
                part_rows = list(rows)
                take = random_allocations[entry_index]
                if take > len(part_rows):
                    raise ValueError(
                        f"Aggregate reservoir in {matrix_path} has {len(part_rows)} "
                        f"rows but {take} are required"
                    )
                if take:
                    chosen = rng.choice(len(part_rows), size=take, replace=False)
                    selected_rows.extend(part_rows[int(index)].copy() for index in chosen)
                continue

            for row in rows:
                if not use_accumulators:
                    if running_sum is None:
                        running_sum = np.zeros_like(row)
                        running_count = np.zeros_like(row, dtype=np.int64)
                    finite = np.isfinite(row)
                    running_sum[finite] += row[finite]
                    assert running_count is not None
                    running_count[finite] += 1
                    valid_total += 1
                if config.write_detail_tables:
                    sampled_index += 1
                    add_heatmap_row(
                        selected_rows,
                        row,
                        sampled_index,
                        config.max_heatmap_rows,
                        config.subsample_mode,
                        rng,
                    )
                if (
                    not use_accumulators
                    and config.stop_after_valid is not None
                    and valid_total >= config.stop_after_valid
                ):
                    stop = True
                    break
            if stop:
                break

    if running_sum is None or running_count is None or valid_total == 0:
        raise RuntimeError(no_valid_regions_message(stats, config))
    if config.write_detail_tables and not selected_rows:
        raise RuntimeError("No vectors were retained for heatmap output")
    stats.valid_total = valid_total
    stats.selected_for_plot = len(selected_rows)
    stats.stopped_after_valid_limit = int(stop)

    full_mean = np.divide(
        running_sum,
        running_count,
        out=np.full_like(running_sum, np.nan),
        where=running_count > 0,
    )
    write_profile(full_mean, outputs["aggregate"])
    if config.nrl:
        from nucleosuite.align import analyse_aggregate_nrl, write_aggregate_nrl_outputs

        exclusion_start, exclusion_end = resolve_nrl_exclusion(config)
        nrl_result = analyse_aggregate_nrl(
            full_mean,
            positions=reference_x,
            peak_resolution=config.nrl_peak_resolution,
            regression_min=config.nrl_regression_min,
            regression_max=config.nrl_regression_max,
            regression_min_order=config.nrl_regression_min_order,
            regression_max_order=config.nrl_regression_max_order,
            exclusion_start=exclusion_start,
            exclusion_end=exclusion_end,
        )
        write_aggregate_nrl_outputs(nrl_result, config, outputs)
    if config.write_detail_tables:
        matrix = np.vstack(selected_rows)
        matrix, x_values = central_crop(matrix, config.breadth)
        matrix, original_order, sort_scores = sort_matrix(matrix, config.sort_mode)
        write_heatmap_matrix(matrix, x_values, outputs["heatmap_matrix"])
        write_heatmap_row_metadata(outputs["row_metadata"], original_order, sort_scores, config.sort_mode)
        plotted_mean = plot_outputs(matrix, x_values, config, outputs)
    else:
        if config.breadth < 1.0:
            cropped, x_values = central_crop(
                full_mean[np.newaxis, :], config.breadth
            )
            plotted_mean = cropped[0]
        else:
            assert reference_x is not None
            x_values = reference_x
            plotted_mean = full_mean
        plot_mean_profile(plotted_mean, x_values, config, outputs)
    write_profile(plotted_mean, outputs["plotted_mean"], x_values)
    write_summary(outputs["summary"], config, stats, outputs)
    return {name: str(path) for name, path in outputs.items()}



def combine_aggregate_manifest(
    root: Path,
    manifest: dict,
    *,
    output_dir: Path | None = None,
) -> dict[str, object]:
    """Rebuild a combined aggregate profile from a stored multicontig manifest."""
    data = dict(manifest.get("aggregate_args", {}))
    for field in (
        "bigwig", "region_bed", "blacklist_bed", "nucleosome_bed", "state_bed", "output_dir",
        "heatmap_output", "heatmap_matrix_output", "aggregate_output",
        "plotted_mean_output", "mean_plot_output", "summary_output",
    ):
        if data.get(field) not in (None, ""):
            data[field] = Path(data[field])
    args = argparse.Namespace(**data)
    entries = list(manifest["per_contig"])
    base_name = str(manifest.get("combined_name") or manifest.get("base_name") or "aggregate")
    outputs = _combine(args, root, entries, base_name, output_dir=output_dir)
    return {
        "input_dir": str(root),
        "output_dir": str(output_dir or (root / "combined")),
        "combined_prefix": str((output_dir or (root / "combined")) / base_name),
        "written": list(outputs.values()),
        "warnings": [],
        "bigwig_stage_last": True,
    }


def run_aggregate_per_contig(args: argparse.Namespace, serial_runner) -> int:
    cores = int(getattr(args, "cores", 1) or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    if getattr(args, "_per_contig_worker", False) or cores == 1:
        return int(serial_runner(args) or 0)
    contigs, chrom_sizes = _resolve_contigs(args)
    if len(contigs) <= 1:
        return int(serial_runner(args) or 0)

    config_values = {key: value for key, value in vars(args).items() if key in AlignmentConfig.__dataclass_fields__}
    base_name = make_output_prefix(AlignmentConfig(**config_values))
    root = (
        Path(args.parallel_dir)
        if getattr(args, "parallel_dir", None)
        else Path(args.output_dir) / f"{base_name}_multicontig"
    ).resolve()
    per_contig = root / "per_contig"
    per_contig.mkdir(parents=True, exist_ok=True)
    namespace_data = {
        key: value
        for key, value in vars(args).items()
        if key not in {"command", "command_function", "command_runner"}
    }

    entries: list[dict[str, object]] = []
    failures: list[str] = []
    with ProcessPoolExecutor(
        max_workers=min(cores, len(contigs)), initializer=_worker_initializer
    ) as executor:
        futures = {}
        for contig in contigs:
            safe = _safe_contig(contig)
            directory = per_contig / safe
            filtered = directory / "regions.bed"
            _filter_bed(
                Path(args.region_bed),
                filtered,
                contig,
                chrom_col=int(args.chrom_col),
                skip_header=bool(args.skip_header),
            )
            prefix = f"{base_name}_{safe}"
            future = executor.submit(
                _worker,
                namespace_data,
                contig,
                str(directory),
                prefix,
                str(filtered),
                _stable_seed(getattr(args, "seed", None), contig),
            )
            futures[future] = (contig, directory / prefix)
        for future in as_completed(futures):
            contig, expected = futures[future]
            try:
                exit_code, prefix, accumulator = future.result()
                if exit_code:
                    failures.append(f"{contig}: exit code {exit_code}")
                entries.append(
                    {
                        "contig": contig,
                        "prefix": prefix,
                        "accumulator": accumulator,
                        "exit_code": exit_code,
                    }
                )
                if exit_code == 0:
                    print(f"Completed {contig}")
            except Exception as error:
                failures.append(f"{contig}: {error}")
                entries.append({"contig": contig, "prefix": str(expected), "exit_code": 2, "error": str(error)})

    order = {contig: index for index, contig in enumerate(contigs)}
    entries.sort(key=lambda item: order[str(item["contig"])])
    manifest = {
        "schema_version": 1,
        "command": "aggregate",
        "combine_strategy": (
            "aggregate_matrix"
            if bool(getattr(args, "write_detail_tables", False))
            or getattr(args, "stop_after_valid", None) is not None
            else "aggregate_accumulator"
        ),
        "base_name": base_name,
        "combined_name": base_name,
        "root_dir": str(root),
        "per_contig_dir": str(per_contig),
        "combined_dir": str(root / "combined"),
        "chrom_sizes": [{"chrom": chrom, "size": size} for chrom, size in chrom_sizes],
        "per_contig": entries,
        "options": {"cores": cores},
        "aggregate_args": {
            key: value
            for key, value in namespace_data.items()
            if key not in {"command_function", "command_runner"}
        },
    }
    with (root / MANIFEST_NAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, default=str)
        handle.write("\n")
    if failures:
        raise RuntimeError("Per-contig aggregate jobs failed: " + "; ".join(failures))
    if not getattr(args, "skip_combine", False):
        outputs = _combine(args, root, entries, base_name)
        print(f"Combined outputs: {root / 'combined'}")
        for path in outputs.values():
            print(path)
        if not bool(getattr(args, "write_detail_tables", False)):
            for entry in entries:
                part_prefix = Path(str(entry["prefix"]))
                for suffix in (
                    "_heatmap_matrix.tsv.gz", "_heatmap_rows.tsv",
                    "_heatmap.png", "_heatmap.svg", "_heatmap_metadata.tsv",
                ):
                    detail_path = part_prefix.parent / f"{part_prefix.name}{suffix}"
                    if detail_path.exists():
                        detail_path.unlink()
    else:
        print(f"Per-contig outputs: {per_contig}")
    return 0
