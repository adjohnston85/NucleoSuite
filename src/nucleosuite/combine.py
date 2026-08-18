"""Combine per-contig NucleoSuite outputs using sufficient statistics.

The combiner is used automatically by multicontig runs and is also exposed as
``nucleosuite combine``.  Tabular outputs are combined before track files.
BigWig files can be combined either by streaming per-contig BigWigs directly
into a new BigWig or by consuming bedGraphs staged during per-contig track
generation. Detailed progress is logged, and incomplete outputs are never
published under the final filename.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import json
import math
import os
import shutil
import subprocess
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter
from nucleosuite.core.chrom_sizes import read_chrom_sizes_source
from nucleosuite.core.regions import resolve_contig_name
from nucleosuite.dac import plot_dac_tsv, write_dac_tsv
from nucleosuite.dcc import plot_dcc_tsv, write_dcc_tsv, write_shift_summary
from nucleosuite.profile_plots import (
    plot_category_counts,
    plot_dinucleotide_profile,
    plot_ww_ss_profile,
    plot_ww_type_length_stacked,
)
from nucleosuite.sequence import dinucleotide
from nucleosuite.sequence.ww_types import (
    WW_TYPE_GROUPS,
    write_length_summary_table,
    write_summary as write_ww_summary,
)

try:
    import pyBigWig
except ImportError:  # pragma: no cover
    pyBigWig = None

MANIFEST_NAME = "nucleosuite_multicontig_manifest.json"
BIGWIG_COMBINE_CHUNK_SIZE = 100_000
# Matplotlib is not thread-safe. Table aggregation remains parallel, while
# figure creation is serialized through this lock inside threaded combine jobs.
_PLOT_LOCK = threading.Lock()


class _CombineLogger:
    """Write combine progress to both the console and a persistent log."""

    def __init__(self, path: Path | None):
        self.path = path
        self.handle = None
        self._lock = threading.Lock()
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            self.handle = path.open("a", encoding="utf-8")

    def log(self, message: str) -> None:
        stamp = datetime.now().astimezone().isoformat(timespec="seconds")
        rendered = f"[COMBINE {stamp}] {message}"
        with self._lock:
            print(rendered, flush=True)
            if self.handle is not None:
                self.handle.write(rendered + "\n")
                self.handle.flush()

    def close(self) -> None:
        with self._lock:
            if self.handle is not None:
                self.handle.close()
                self.handle = None


def _rss_mb() -> float | None:
    """Return current resident memory on Linux/WSL when available."""
    status = Path("/proc/self/status")
    if not status.exists():
        return None
    try:
        for line in status.read_text(encoding="utf-8").splitlines():
            if line.startswith("VmRSS:"):
                return float(line.split()[1]) / 1024.0
    except (OSError, ValueError, IndexError):
        return None
    return None


def _format_bytes(value: int) -> str:
    amount = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if amount < 1024.0 or unit == "TiB":
            return f"{amount:.1f} {unit}"
        amount /= 1024.0
    return f"{amount:.1f} TiB"


def _run_logged_subprocess(command: Sequence[str], logger: _CombineLogger | None) -> None:
    if logger is not None:
        logger.log("Command: " + " ".join(str(value) for value in command))
    process = subprocess.Popen(
        [str(value) for value in command],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    for line in process.stdout:
        text = line.rstrip("\n")
        if logger is not None:
            logger.log(f"tool: {text}")
        else:
            print(text, flush=True)
    return_code = process.wait()
    if return_code != 0:
        if return_code < 0:
            reason = f"terminated by signal {-return_code}"
        else:
            reason = f"exit status {return_code}"
        raise subprocess.CalledProcessError(return_code, command, output=reason)


def _safe_contig(value: str) -> str:
    return value.replace("/", "_").replace("\\", "_").replace(":", "_")


def _open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return [], []
        return list(reader.fieldnames), list(reader)


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_manifest(input_dir: Path) -> dict:
    manifest_path = input_dir / MANIFEST_NAME
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"No {MANIFEST_NAME} was found in {input_dir}. "
            "Use --manifest to identify a multicontig run directory."
        )
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not manifest.get("per_contig"):
        raise ValueError(f"{manifest_path} contains no per-contig outputs")
    return manifest


def _chrom_sizes_from_manifest(manifest: Mapping[str, object]) -> list[tuple[str, int]]:
    raw = manifest.get("chrom_sizes", [])
    output: list[tuple[str, int]] = []
    for item in raw:
        if isinstance(item, dict):
            output.append((str(item["chrom"]), int(item["size"])))
        else:
            chrom, size = item
            output.append((str(chrom), int(size)))
    return output


def _write_chrom_sizes(path: Path, chrom_sizes: Sequence[tuple[str, int]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for chrom, size in chrom_sizes:
            handle.write(f"{chrom}\t{int(size)}\n")


def _output_groups(manifest: Mapping[str, object]) -> dict[str, list[Path]]:
    """Group files by the suffix following each per-contig output prefix."""
    groups: dict[str, list[Path]] = defaultdict(list)
    for entry in manifest["per_contig"]:
        prefix = Path(str(entry["prefix"]))
        parent = prefix.parent
        if not parent.exists():
            continue
        prefix_name = prefix.name
        for path in sorted(parent.glob(prefix_name + "*")):
            if not path.is_file():
                continue
            suffix = path.name[len(prefix_name):]
            if not suffix:
                continue
            if suffix.endswith((".bb.empty", ".bigBed.empty")):
                groups.setdefault(suffix[:-6], [])
                continue
            groups[suffix].append(path)
    return dict(groups)


def _combined_path(combined_prefix: Path, suffix: str) -> Path:
    return combined_prefix.parent / f"{combined_prefix.name}{suffix}"


def _concatenate_text(inputs: Sequence[Path], output: Path, *, keep_one_header: bool = False) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with _open_text(output, "wt") as destination:
        first_line: str | None = None
        for input_index, path in enumerate(inputs):
            with _open_text(path) as source:
                for line_index, line in enumerate(source):
                    if keep_one_header and line_index == 0:
                        if first_line is None:
                            first_line = line
                            destination.write(line)
                        elif line == first_line:
                            continue
                        else:
                            raise ValueError(f"Header mismatch while combining {path}")
                    else:
                        destination.write(line)


def _interval_fields_and_key(
    text: str,
    *,
    path: Path,
    line_number: int,
    canonical: Sequence[str],
    rank: Mapping[str, int],
) -> tuple[list[str], tuple[int, int, int]]:
    fields = text.split("\t") if "\t" in text else text.split()
    if len(fields) < 3:
        raise ValueError(f"{path}:{line_number}: expected at least three BED columns")
    try:
        fields[0] = resolve_contig_name(
            fields[0], canonical, source_label="combined chromosome sizes"
        )
    except KeyError as exc:
        raise ValueError(str(exc)) from exc
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_number}: BED start and end must be integers"
        ) from exc
    return fields, (rank[fields[0]], start, end)


def _write_interval_sort_chunk(
    records: list[tuple[tuple[int, int, int], str]],
    path: Path,
) -> None:
    records.sort(key=lambda item: item[0])
    with path.open("wt", encoding="utf-8") as handle:
        for (chrom_rank, start, end), line in records:
            handle.write(f"{chrom_rank}\t{start}\t{end}\t{line}\n")


def _iter_interval_sort_chunk(path: Path):
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            chrom_rank, start, end, line = raw.rstrip("\n").split("\t", 3)
            yield (int(chrom_rank), int(start), int(end)), line


def _external_sort_interval_file(
    input_path: Path,
    output_path: Path,
    chrom_sizes: Sequence[tuple[str, int]],
    *,
    records_per_chunk: int = 250_000,
) -> None:
    """Coordinate-sort a BED-like file with bounded memory.

    This is a defensive fallback. Normal combined output remains streaming and
    incurs no sort when its rows are already ordered.
    """

    canonical = [name for name, _length in chrom_sizes]
    rank = {name: index for index, name in enumerate(canonical)}
    headers: list[str] = []
    records: list[tuple[tuple[int, int, int], str]] = []
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=".nucleosuite_interval_sort_", dir=output_path.parent
    ) as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        chunks: list[Path] = []
        with _open_text(input_path) as source:
            for line_number, raw in enumerate(source, 1):
                text = raw.rstrip("\n")
                if not text or text.startswith(("#", "track", "browser")):
                    headers.append(raw if raw.endswith("\n") else raw + "\n")
                    continue
                fields, key = _interval_fields_and_key(
                    text,
                    path=input_path,
                    line_number=line_number,
                    canonical=canonical,
                    rank=rank,
                )
                records.append((key, "\t".join(fields)))
                if len(records) >= records_per_chunk:
                    chunk = temp_dir / f"chunk_{len(chunks):06d}.tsv"
                    _write_interval_sort_chunk(records, chunk)
                    chunks.append(chunk)
                    records = []
        if records:
            chunk = temp_dir / f"chunk_{len(chunks):06d}.tsv"
            _write_interval_sort_chunk(records, chunk)
            chunks.append(chunk)

        with _open_text(output_path, "wt") as destination:
            destination.writelines(headers)
            iterators = [_iter_interval_sort_chunk(chunk) for chunk in chunks]
            for _key, line in heapq.merge(*iterators, key=lambda item: item[0]):
                destination.write(line + "\n")


def _interval_file_has_inversion(
    path: Path,
    chrom_sizes: Sequence[tuple[str, int]],
) -> bool:
    canonical = [name for name, _length in chrom_sizes]
    rank = {name: index for index, name in enumerate(canonical)}
    previous_key: tuple[int, int, int] | None = None
    with _open_text(path) as source:
        for line_number, raw in enumerate(source, 1):
            text = raw.rstrip("\n")
            if not text or text.startswith(("#", "track", "browser")):
                continue
            _fields, key = _interval_fields_and_key(
                text,
                path=path,
                line_number=line_number,
                canonical=canonical,
                rank=rank,
            )
            if previous_key is not None and key < previous_key:
                return True
            previous_key = key
    return False


def _replace_with_sorted_interval_file(
    path: Path,
    chrom_sizes: Sequence[tuple[str, int]],
) -> None:
    suffix = "".join(path.suffixes) or ".bed"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.stem}.sorted.", suffix=suffix, dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        _external_sort_interval_file(path, temporary, chrom_sizes)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _concatenate_intervals(
    inputs: Sequence[Path],
    output: Path,
    chrom_sizes: Sequence[tuple[str, int]],
) -> None:
    """Concatenate BED-like files and sort only if an inversion is detected."""
    canonical = [name for name, _length in chrom_sizes]
    rank = {name: index for index, name in enumerate(canonical)}
    output.parent.mkdir(parents=True, exist_ok=True)
    previous_key: tuple[int, int, int] | None = None
    inversion_found = False
    with _open_text(output, "wt") as destination:
        for path in inputs:
            with _open_text(path) as source:
                for line_number, raw in enumerate(source, 1):
                    text = raw.rstrip("\n")
                    if not text or text.startswith(("#", "track", "browser")):
                        destination.write(raw)
                        continue
                    fields, key = _interval_fields_and_key(
                        text,
                        path=path,
                        line_number=line_number,
                        canonical=canonical,
                        rank=rank,
                    )
                    if previous_key is not None and key < previous_key:
                        inversion_found = True
                    previous_key = key
                    destination.write("\t".join(fields) + "\n")
    if inversion_found:
        _replace_with_sorted_interval_file(output, chrom_sizes)


def _combine_fragment_length_counts(inputs: Sequence[Path], output: Path) -> None:
    counts: Counter[int] = Counter()
    for path in inputs:
        _, rows = _read_tsv(path)
        for row in rows:
            counts[int(row["fragment_length"])] += int(float(row["count"]))
    _write_tsv(
        output,
        ["fragment_length", "count"],
        ({"fragment_length": length, "count": counts[length]} for length in sorted(counts)),
    )



def _combine_metric_value_summary(
    inputs: Sequence[Path],
    output: Path,
    *,
    additive_metrics: set[str],
) -> None:
    sums: Counter[str] = Counter()
    fixed: dict[str, str] = {}
    order: list[str] = []
    for path in inputs:
        _, rows = _read_tsv(path)
        for row in rows:
            metric = row.get("metric", "")
            value = row.get("value", "")
            if metric not in order:
                order.append(metric)
            if metric in additive_metrics:
                sums[metric] += float(value or 0)
            elif metric not in fixed:
                fixed[metric] = value
            elif fixed[metric] != value:
                fixed[metric] = "mixed"
    rows_out = []
    for metric in order:
        if metric in additive_metrics:
            value = sums[metric]
            value = int(value) if value.is_integer() else f"{value:.12g}"
        else:
            value = fixed.get(metric, "")
        rows_out.append({"metric": metric, "value": value})
    _write_tsv(output, ["metric", "value"], rows_out)




def _combine_randomization_qc(inputs: Sequence[Path], output: Path) -> None:
    additive_names = {
        "input", "matched", "uniform", "fallback", "skipped",
        "fragments_excluded_by_blacklist",
        "random_candidates_rejected_by_blacklist",
        "random_candidates_rejected_non_acgt",
        "anchor_start_selected", "anchor_end_selected",
        "anchor_start_matched", "anchor_end_matched",
        "anchor_start", "anchor_end",
        "unique_randomized_coordinates", "duplicate_randomized_fragments",
        "blacklist_intervals", "blacklisted_bases",
    }
    sums: Counter[str] = Counter()
    fixed: dict[str, str] = {}
    order: list[str] = []
    maximum_multiplicity = 0
    for path in inputs:
        _fields, rows = _read_tsv(path)
        for row in rows:
            metric = row.get("metric", "")
            value = row.get("value", "")
            if metric not in order:
                order.append(metric)
            if metric in additive_names or metric.startswith("reason_"):
                sums[metric] += float(value or 0)
            elif metric == "maximum_randomized_multiplicity":
                maximum_multiplicity = max(maximum_multiplicity, int(float(value or 0)))
            elif metric == "collision_fraction":
                continue
            elif metric not in fixed:
                fixed[metric] = value
            elif fixed[metric] != value:
                fixed[metric] = "mixed"
    for required in ("maximum_randomized_multiplicity", "collision_fraction"):
        if required not in order:
            order.append(required)
    randomized_total = sums["matched"] + sums["uniform"] + sums["fallback"]
    duplicate_total = sums["duplicate_randomized_fragments"]
    rows_out = []
    for metric in order:
        if metric in additive_names or metric.startswith("reason_"):
            value = sums[metric]
            value = int(value) if float(value).is_integer() else f"{value:.12g}"
        elif metric == "maximum_randomized_multiplicity":
            value = maximum_multiplicity
        elif metric == "collision_fraction":
            value = f"{(duplicate_total / randomized_total if randomized_total else 0.0):.12g}"
        else:
            value = fixed.get(metric, "")
        rows_out.append({"metric": metric, "value": value})
    _write_tsv(output, ["metric", "value"], rows_out)


def _combine_fragment_summary(inputs: Sequence[Path], output: Path) -> None:
    additive = {
        "total_fragments_filtered_all",
        "total_fragments_used_in_range",
        "unique_bases_covered_by_used_fragments",
        "fragments_written",
        "contigs_written",
    }
    sums: Counter[str] = Counter()
    fixed: dict[str, str] = {}
    order: list[str] = []
    for path in inputs:
        _, rows = _read_tsv(path)
        for row in rows:
            metric, value = row["metric"], row["value"]
            if metric not in order:
                order.append(metric)
            if metric in additive:
                sums[metric] += int(float(value))
            elif metric not in fixed:
                fixed[metric] = value
            elif fixed[metric] != value:
                fixed[metric] = "mixed"
    rows = []
    for metric in order:
        rows.append({"metric": metric, "value": sums[metric] if metric in additive else fixed.get(metric, "")})
    _write_tsv(output, ["metric", "value"], rows)


def _combine_ww_summary(inputs: Sequence[Path], output: Path) -> None:
    counts: Counter[str] = Counter()
    for path in inputs:
        _, rows = _read_tsv(path)
        for row in rows:
            if row["type"] == "all":
                continue
            counts[row["type"]] += int(float(row["fragment_count"]))
    total = sum(counts.values())
    prefix = str(output)
    suffix = "_ww_type_summary.tsv"
    if prefix.endswith(suffix):
        prefix = prefix[: -len(suffix)]
    write_ww_summary(prefix, counts, total)
    png = Path(prefix + "_ww_type_counts.png")
    try:
        with _PLOT_LOCK:
            plot_category_counts(str(output), str(png), title="WW/SS fragment classes")
    except Exception:
        pass


def _combine_ww_type_by_length(inputs: Sequence[Path], output: Path) -> None:
    counts_by_length: dict[int, Counter[str]] = defaultdict(Counter)
    for path in inputs:
        _, rows = _read_tsv(path)
        for row in rows:
            fragment_length = int(float(row["fragment_length"]))
            for group in WW_TYPE_GROUPS:
                counts_by_length[fragment_length][group] += int(
                    float(row.get(f"{group}_count", 0) or 0)
                )
            counts_by_length[fragment_length]["unclassified"] += int(
                float(row.get("unclassified_count", 0) or 0)
            )
    write_length_summary_table(output, counts_by_length)
    plot_ww_type_length_stacked(
        output,
        output.with_name(output.stem + "_stacked.png"),
        title="WW/SS type frequencies by fragment length",
    )


def _combine_dinuc_counts(inputs: Sequence[Path], output: Path) -> None:
    n_valid: Counter[int] = Counter()
    counts: dict[int, Counter[str]] = defaultdict(Counter)
    fragments_used = 0
    fragments_skipped = 0
    for path in inputs:
        fields, rows = _read_tsv(path)
        dinuc_columns = [name for name in fields if name.endswith("_count") and name not in {"fragments_used_count", "fragments_skipped_count"}]
        for row in rows:
            position = int(row["position"])
            n_valid[position] += int(float(row["n_valid"]))
            for column in dinuc_columns:
                counts[position][column[:-6]] += int(float(row[column]))
            if "fragments_used" in row and row["fragments_used"]:
                fragments_used = max(fragments_used, 0) + int(float(row["fragments_used"]))
            if "fragments_skipped" in row and row["fragments_skipped"]:
                fragments_skipped = max(fragments_skipped, 0) + int(float(row["fragments_skipped"]))
    positions = sorted(n_valid)
    fieldnames = ["position", "n_valid", *[f"{name}_count" for name in dinucleotide.DINUCS]]
    if fragments_used or fragments_skipped:
        fieldnames.extend(["fragments_used", "fragments_skipped"])
    rows = []
    for index, position in enumerate(positions):
        row: dict[str, object] = {"position": position, "n_valid": n_valid[position]}
        for name in dinucleotide.DINUCS:
            row[f"{name}_count"] = counts[position][name]
        if "fragments_used" in fieldnames:
            row["fragments_used"] = fragments_used if index == 0 else 0
            row["fragments_skipped"] = fragments_skipped if index == 0 else 0
        rows.append(row)
    _write_tsv(output, fieldnames, rows)


def _profile_from_counts(counts_path: Path, profile_path: Path, fraction: bool) -> None:
    _, rows = _read_tsv(counts_path)
    accumulator = dinucleotide.new_accumulator()
    positions: list[int] = []
    for row in rows:
        position = int(row["position"])
        positions.append(position)
        n_valid = int(float(row["n_valid"]))
        accumulator["n_valid"][position] = n_valid
        for name in dinucleotide.DINUCS:
            accumulator["counts"][position][name] = int(float(row[f"{name}_count"]))
        if row.get("fragments_used"):
            accumulator["fragments_used"] += int(float(row["fragments_used"]))
        if row.get("fragments_skipped"):
            accumulator["fragments_skipped"] += int(float(row["fragments_skipped"]))
    dinucleotide.write_profile(profile_path, accumulator, positions, fraction=fraction, write_count_table=False)
    try:
        with _PLOT_LOCK:
            plot_dinucleotide_profile(str(profile_path), str(profile_path.with_suffix(".png")), title="Dinucleotide profile")
            plot_ww_ss_profile(str(profile_path), str(profile_path.with_name(profile_path.stem + "_ww_ss.png")), title="WW/SS dinucleotide profile")
    except Exception:
        pass


def _find_summary_values(tsv_path: Path, *, kind: str) -> tuple[float, ...]:
    parent = tsv_path.parent
    patterns = ("*_DAC_summary.tsv", "*_DAC_run_summary.tsv") if kind == "dac" else ("*_DCC_summary.tsv", "*_DCC_run_summary.tsv")
    basename = tsv_path.name
    for pattern in patterns:
        for summary in parent.glob(pattern):
            _, rows = _read_tsv(summary)
            for row in rows:
                if Path(row.get("Output", "")).name != basename:
                    continue
                if kind == "dac":
                    return (float(row.get("Total signal", 0.0)),)
                return (
                    float(row.get("Total signal A", 0.0)),
                    float(row.get("Total signal B", 0.0)),
                )
    return (0.0,) if kind == "dac" else (0.0, 0.0)


def _infer_scale(rows: Sequence[dict[str, str]], raw_column: str, scaled_column: str, denominator: float) -> float:
    if denominator <= 0:
        return 1_000_000.0
    for row in rows:
        raw = float(row[raw_column])
        scaled = float(row[scaled_column])
        if raw != 0 and scaled != 0:
            return scaled * denominator / raw
    return 1_000_000.0


def _combine_dac(inputs: Sequence[Path], output: Path) -> None:
    raw: dict[int, float] = defaultdict(float)
    opportunities: dict[int, float] = defaultdict(float)
    total_signal = 0.0
    scale_values: list[float] = []
    for path in inputs:
        _, rows = _read_tsv(path)
        (signal,) = _find_summary_values(path, kind="dac")
        total_signal += signal
        scale_values.append(_infer_scale(rows, "Raw DAC Value", "DAC per million signal-pairs", signal * signal))
        for row in rows:
            distance = int(float(row["Distance"]))
            raw[distance] += float(row["Raw DAC Value"])
            opportunities[distance] += float(row["Opportunities"])
    maximum = max(raw, default=0)
    raw_array = np.zeros(maximum + 1, dtype=float)
    opportunity_array = np.zeros_like(raw_array)
    for distance, value in raw.items():
        raw_array[distance] = value
        opportunity_array[distance] = opportunities[distance]
    normalize = "opportunity_normalized" in output.name
    write_dac_tsv(output, raw_array, opportunity_array, normalize, total_signal, float(np.median(scale_values) if scale_values else 1_000_000.0))
    with _PLOT_LOCK:
        plot_dac_tsv(output, output.with_suffix(".png"), title=output.stem)


def _combine_dcc(inputs: Sequence[Path], output: Path) -> None:
    first_fields, first_rows = _read_tsv(inputs[0])
    x_label = first_fields[0]
    raw: dict[int, float] = defaultdict(float)
    opportunities: dict[int, float] = defaultdict(float)
    total_a = 0.0
    total_b = 0.0
    scale_values: list[float] = []
    for path in inputs:
        _, rows = _read_tsv(path)
        signal_a, signal_b = _find_summary_values(path, kind="dcc")
        total_a += signal_a
        total_b += signal_b
        scale_values.append(_infer_scale(rows, "Raw DCC Value", "DCC per million signal-pairs", signal_a * signal_b))
        for row in rows:
            x = int(float(row[x_label]))
            raw[x] += float(row["Raw DCC Value"])
            opportunities[x] += float(row["Opportunities"])
    signed = x_label == "Lag"
    if signed:
        dmax = max((abs(value) for value in raw), default=0)
        x_values = list(range(-dmax, dmax + 1))
    else:
        dmax = max(raw, default=0)
        x_values = list(range(0, dmax + 1))
    raw_array = np.array([raw[x] for x in x_values], dtype=float)
    opportunity_array = np.array([opportunities[x] for x in x_values], dtype=float)
    normalize = "opportunity_normalized" in output.name
    normalize_totals = "signal_total_normalized" in output.name
    scale = float(np.median(scale_values) if scale_values else 1_000_000.0)
    write_dcc_tsv(output, raw_array, opportunity_array, dmax, signed, normalize, normalize_totals, total_a, total_b, scale)
    with _PLOT_LOCK:
        plot_dcc_tsv(output, output.with_suffix(".png"), title=output.stem)
    try:
        from nucleosuite.dcc import build_reported_dcc
        reported = build_reported_dcc(raw_array, opportunity_array, normalize, normalize_totals, total_a, total_b)
        write_shift_summary(output.with_name(output.stem + "_shift_summary.tsv"), reported, dmax, signed, 25)
    except Exception:
        pass




def _combine_positive_summary(inputs: Sequence[Path], output: Path, counts_path: Path) -> None:
    from nucleosuite.positive_runs import weighted_quantile

    additive = {
        "selected_regions", "scanned_bases", "positive_bases", "nonpositive_bases",
        "missing_bases", "total_runs_observed", "total_runs_retained",
        "retained_positive_bases",
    }
    order: list[str] = []
    sums: Counter[str] = Counter()
    fixed: dict[str, str] = {}
    for path in inputs:
        _, rows = _read_tsv(path)
        for row in rows:
            metric, value = row["metric"], row["value"]
            if metric not in order:
                order.append(metric)
            if metric in additive:
                sums[metric] += float(value or 0)
            elif metric not in fixed:
                fixed[metric] = value
            elif fixed[metric] != value and metric != "bigwig":
                fixed[metric] = "mixed"

    counts: Counter[int] = Counter()
    if counts_path.exists():
        _, rows = _read_tsv(counts_path)
        for row in rows:
            counts[int(float(row["run_length_bp"]))] += int(float(row["count"]))
    total = sum(counts.values())
    mean = sum(length * count for length, count in counts.items()) / total if total else math.nan
    mode = min((length for length, count in counts.items() if count == max(counts.values())), default=math.nan)
    derived = {
        "minimum_run_length_bp": min(counts) if counts else math.nan,
        "q1_run_length_bp": weighted_quantile(counts, 0.25),
        "median_run_length_bp": weighted_quantile(counts, 0.5),
        "mean_run_length_bp": mean,
        "q3_run_length_bp": weighted_quantile(counts, 0.75),
        "maximum_run_length_bp": max(counts) if counts else math.nan,
        "mode_run_length_bp": mode,
    }
    rows_out = []
    for metric in order:
        if metric in additive:
            value = sums[metric]
            value = int(value) if value.is_integer() else f"{value:.12g}"
        elif metric in derived:
            value = derived[metric]
        elif metric == "bigwig":
            value = "combined"
        else:
            value = fixed.get(metric, "")
        rows_out.append({"metric": metric, "value": value})
    _write_tsv(output, ["metric", "value"], rows_out)


def _combine_run_summary(
    inputs: Sequence[Path],
    output: Path,
    *,
    kind: str,
    combined_prefix: Path,
) -> None:
    """Combine DAC/DCC run summaries while retaining one row per analysis state."""
    additive = {
        "Regions",
        "Region-track pairs used",
        "Missing-chromosome skips",
        "Short-or-empty skips",
        "Clipped regions",
        "Signal positions",
        "Non-zero signal positions",
        "Used regions",
        "Regions with A",
        "Regions with B",
        "Missing chromosome A",
        "Missing chromosome B",
        "Short or empty",
        "Clipped regions A",
        "Clipped regions B",
        "Signal positions A",
        "Signal positions B",
        "Non-zero positions A",
        "Non-zero positions B",
        "Total signal",
        "Total signal A",
        "Total signal B",
        "Sparse calculations",
        "FFT calculations",
    }
    maximum = {"BigWig files", "A files", "B files"}
    fields, _ = _read_tsv(inputs[0])
    groups: dict[tuple[str, ...], dict[str, object]] = {}
    sums: dict[tuple[str, ...], Counter[str]] = defaultdict(Counter)
    maxima: dict[tuple[str, ...], dict[str, float]] = defaultdict(dict)
    key_fields = [field for field in ("State", "Mode", "Signal") if field in fields]
    if not key_fields:
        key_fields = [fields[0]] if fields else []

    for summary_path in inputs:
        current_fields, rows = _read_tsv(summary_path)
        if current_fields != fields:
            raise ValueError(f"Header mismatch while combining {summary_path}")
        summary_suffix = f"_{kind.upper()}_summary.tsv"
        summary_prefix = summary_path.name[:-len(summary_suffix)] if summary_path.name.endswith(summary_suffix) else summary_path.stem
        for row in rows:
            key = tuple(row.get(field, "") for field in key_fields)
            target = groups.setdefault(key, dict(row))
            for field in additive:
                if field in row and row[field] not in ("", None):
                    sums[key][field] += float(row[field])
            for field in maximum:
                if field in row and row[field] not in ("", None):
                    maxima[key][field] = max(maxima[key].get(field, 0.0), float(row[field]))
            source_output = Path(row.get("Output", "")).name
            if source_output.startswith(summary_prefix):
                suffix = source_output[len(summary_prefix):]
                target["Output"] = str(combined_prefix.parent / f"{combined_prefix.name}{suffix}")
                if "Shift summary" in fields:
                    target["Shift summary"] = str(
                        combined_prefix.parent / f"{combined_prefix.name}{Path(suffix).stem}_shift_summary.tsv"
                    )

    rows_out: list[dict[str, object]] = []
    for key in groups:
        row = groups[key]
        for field, value in sums[key].items():
            row[field] = int(value) if field not in {"Total signal", "Total signal A", "Total signal B"} and value.is_integer() else f"{value:.12g}"
        for field, value in maxima[key].items():
            row[field] = int(value) if value.is_integer() else f"{value:.12g}"
        rows_out.append(row)
    _write_tsv(output, fields, rows_out)

def _worker_sample_prefix(sample_name: str, worker_name: str) -> str:
    """Return the worker prefix while keeping randomized-control terminal."""
    marker = "_randomized_control"
    if sample_name.endswith(marker):
        return f"{sample_name[:-len(marker)]}_{worker_name}{marker}"
    return f"{sample_name}_{worker_name}"


def _normalise_worker_relative_path(
    path: Path,
    *,
    root: Path,
    sample_name: str | None,
) -> Path:
    """Return a per-worker path relative to its workflow root with sample suffix removed."""
    relative = path.resolve(strict=False).relative_to(root.resolve())
    if sample_name:
        worker_prefix = _worker_sample_prefix(sample_name, root.name)
        if relative.name.startswith(worker_prefix):
            relative = relative.with_name(
                sample_name + relative.name[len(worker_prefix):]
            )
    return relative


def _combine_track_metadata_table(
    inputs: Sequence[Path],
    output: Path,
    *,
    roots: Sequence[Path],
    destination_root: Path,
    sample_name: str | None,
) -> None:
    """Regenerate a track manifest/report with combined output-prefix paths.

    Per-contig manifests and completion reports contain worker-local absolute
    output prefixes. A simple concatenation produces misleading rows. This
    function rewrites each prefix into the combined directory, verifies that
    every worker describes the same requested track set, and writes one row per
    combined specification.
    """
    expected_fields, _ = _read_tsv(inputs[0])
    if "output_prefix" not in expected_fields:
        raise ValueError(f"{inputs[0]} has no output_prefix column")
    canonical_by_input: list[set[tuple[str, ...]]] = []
    canonical_rows: dict[tuple[str, ...], dict[str, object]] = {}
    for path in inputs:
        fields, rows = _read_tsv(path)
        if fields != expected_fields:
            raise ValueError(f"Header mismatch while combining {path}")
        root = next(
            (candidate for candidate in roots if path.is_relative_to(candidate)),
            None,
        )
        if root is None:
            raise ValueError(f"Could not associate metadata table with a worker root: {path}")
        keys: set[tuple[str, ...]] = set()
        for row in rows:
            prefix = Path(row["output_prefix"])
            relative = _normalise_worker_relative_path(
                prefix, root=root, sample_name=sample_name
            )
            transformed: dict[str, object] = dict(row)
            transformed["output_prefix"] = str(destination_root / relative)
            key = tuple(str(transformed.get(field, "")) for field in expected_fields)
            keys.add(key)
            canonical_rows[key] = transformed
        canonical_by_input.append(keys)
    reference = canonical_by_input[0]
    for path, keys in zip(inputs[1:], canonical_by_input[1:]):
        if keys != reference:
            raise ValueError(
                f"Track specification mismatch between {inputs[0]} and {path}"
            )
    rows_out = [canonical_rows[key] for key in sorted(reference)]
    _write_tsv(output, expected_fields, rows_out)


def _supported_tree_tsv_name(name: str) -> bool:
    """Return whether directory-tree combination has an explicit TSV strategy."""
    return (
        _is_track_metadata_name(name)
        or name.endswith(("_fragment_length_counts.tsv", ".fragment_length_counts.tsv"))
        or name.endswith(".randomization_qc.tsv")
        or name.endswith(".relocation_distances.tsv")
        or name.endswith(".merge_summary.tsv")
        or name.endswith(("_fragment_summary.tsv", ".fragments.summary.tsv"))
        or name.endswith("_ww_type_summary.tsv")
        or name.endswith("_ww_type_by_length.tsv")
        or name.endswith("_dinuc_profile_counts.tsv")
        or name.endswith("_dinuc_profile.tsv")
        or name.endswith("_summary.tsv")
        or ("_DAC_" in name and name.endswith(".tsv"))
        or ("_DCC" in name and name.endswith(".tsv"))
    )


def _is_track_metadata_name(name: str) -> bool:
    """Recognise every supported track-metadata table name."""
    return (
        name in {"manifest.tsv", "completion_report.tsv"}
        or name.endswith("_manifest.tsv")
        or name.endswith("_completion_report.tsv")
    )


def _combine_generic_tsv(inputs: Sequence[Path], output: Path) -> None:
    fields, _ = _read_tsv(inputs[0])
    if "count" not in fields:
        _concatenate_text(inputs, output, keep_one_header=True)
        return

    derived = {"fraction", "percent", "density"}
    key_fields = [field for field in fields if field != "count" and field not in derived]
    totals: dict[tuple[str, ...], float] = defaultdict(float)
    for path in inputs:
        current_fields, rows = _read_tsv(path)
        if current_fields != fields:
            raise ValueError(f"Header mismatch while combining {path}")
        for row in rows:
            key = tuple(row[field] for field in key_fields)
            totals[key] += float(row["count"])

    # Fractions and percentages are recalculated from the combined counts. If a
    # label column is present, each label is normalized independently.
    label_index = key_fields.index("label") if "label" in key_fields else None
    denominators: Counter[str] = Counter()
    for key, value in totals.items():
        group = key[label_index] if label_index is not None else "all"
        denominators[group] += value

    rows: list[dict[str, object]] = []
    for key in sorted(totals):
        count = totals[key]
        row: dict[str, object] = dict(zip(key_fields, key))
        row["count"] = int(count) if float(count).is_integer() else f"{count:.12g}"
        group = key[label_index] if label_index is not None else "all"
        denominator = denominators[group]
        if "fraction" in fields:
            row["fraction"] = f"{count / denominator:.12g}" if denominator else "NaN"
        if "percent" in fields:
            row["percent"] = f"{100.0 * count / denominator:.12g}" if denominator else "NaN"
        if "density" in fields:
            row["density"] = f"{count / denominator:.12g}" if denominator else "NaN"
        rows.append(row)
    _write_tsv(output, fields, rows)


def _bigwig_to_bedgraph(
    input_path: Path,
    output_path: Path,
    chrom_order: Sequence[str],
    *,
    chrom_lengths: Mapping[str, int] | None = None,
    chunk_size: int = BIGWIG_COMBINE_CHUNK_SIZE,
    logger: _CombineLogger | None = None,
) -> tuple[int, int]:
    """Stream a BigWig to bedGraph without materialising a whole chromosome.

    ``pyBigWig.intervals(chrom)`` returns all intervals for the chromosome as a
    Python list. Dense PNS, WPS and coverage tracks can make that list several
    gigabytes. Querying fixed genomic chunks bounds memory while preserving the
    original run intervals. Intervals overlapping a chunk boundary are emitted
    only by the chunk that owns their start coordinate.
    """
    if pyBigWig is None:
        raise RuntimeError("BigWig combination requires pyBigWig")
    if chunk_size < 1:
        raise ValueError("BigWig combine chunk size must be positive")
    handle = pyBigWig.open(str(input_path))
    if handle is None:
        raise RuntimeError(f"Could not open BigWig: {input_path}")
    interval_count = 0
    queried_chunks = 0
    try:
        available = handle.chroms()
        available_names = list(available)
        with output_path.open("w", encoding="utf-8") as output:
            for chrom in chrom_order:
                try:
                    source_chrom = resolve_contig_name(
                        chrom, available_names, source_label=f"BigWig {input_path}"
                    )
                except KeyError:
                    continue
                source_length = int(available[source_chrom])
                requested_length = (chrom_lengths or {}).get(chrom, source_length)
                length = min(source_length, int(requested_length))
                chrom_intervals = 0
                chrom_chunks = max(1, math.ceil(length / chunk_size)) if length else 0
                if logger is not None:
                    logger.log(
                        f"Streaming {input_path.name}: {chrom} ({length:,} bp; "
                        f"{chrom_chunks:,} chunks of up to {chunk_size:,} bp)"
                    )
                for chunk_index, chunk_start in enumerate(range(0, length, chunk_size), start=1):
                    chunk_end = min(length, chunk_start + chunk_size)
                    intervals = handle.intervals(source_chrom, chunk_start, chunk_end) or ()
                    queried_chunks += 1
                    for start, end, value in intervals:
                        start_i = int(start)
                        if start_i < chunk_start or start_i >= chunk_end:
                            continue
                        output.write(
                            f"{chrom}\t{start_i}\t{int(end)}\t{float(value):.12g}\n"
                        )
                        interval_count += 1
                        chrom_intervals += 1
                    if logger is not None and (
                        chunk_index == chrom_chunks or chunk_index % 100 == 0
                    ):
                        rss = _rss_mb()
                        memory = "" if rss is None else f"; RSS {rss:.1f} MiB"
                        logger.log(
                            f"  {chrom}: {chunk_index:,}/{chrom_chunks:,} chunks; "
                            f"{chrom_intervals:,} intervals written{memory}"
                        )
    finally:
        handle.close()
    return interval_count, queried_chunks


def _verify_bigwig(path: Path, expected_chroms: Sequence[str]) -> None:
    if pyBigWig is None:
        return
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise RuntimeError(f"Combined BigWig could not be opened: {path}")
    try:
        observed = handle.chroms()
        missing = [chrom for chrom in expected_chroms if chrom not in observed]
        if missing:
            raise RuntimeError(f"Combined BigWig is missing contigs: {', '.join(missing)}")
    finally:
        handle.close()


def _combined_bigwig_marker_path(output: Path) -> Path:
    return output.with_name(output.name + ".combine.complete.json")


def _bigwig_input_signature(inputs: Sequence[Path]) -> list[dict[str, object]]:
    signature: list[dict[str, object]] = []
    for path in inputs:
        if path.exists():
            stat = path.stat()
            size_bytes = int(stat.st_size)
            mtime_ns = int(stat.st_mtime_ns)
        else:
            size_bytes = -1
            mtime_ns = -1
        signature.append(
            {
                "path": str(path.resolve()),
                "size_bytes": size_bytes,
                "mtime_ns": mtime_ns,
            }
        )
    return signature


def _write_combined_bigwig_marker(
    output: Path,
    *,
    inputs: Sequence[Path],
    chrom_sizes: Sequence[tuple[str, int]],
    method: str,
) -> None:
    marker = _combined_bigwig_marker_path(output)
    payload = {
        "schema_version": 1,
        "complete": True,
        "method": str(method),
        "output": str(output.resolve()),
        "output_size_bytes": output.stat().st_size,
        "chrom_sizes": [[chrom, int(size)] for chrom, size in chrom_sizes],
        "inputs": _bigwig_input_signature(inputs),
    }
    temporary = marker.with_name(marker.name + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, marker)


def _combined_bigwig_is_reusable(
    output: Path,
    *,
    inputs: Sequence[Path],
    chrom_sizes: Sequence[tuple[str, int]],
    method: str,
) -> bool:
    marker = _combined_bigwig_marker_path(output)
    if not output.is_file() or not marker.is_file():
        return False
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
        if not payload.get("complete") or payload.get("method") != str(method):
            return False
        if payload.get("chrom_sizes") != [
            [chrom, int(size)] for chrom, size in chrom_sizes
        ]:
            return False
        if payload.get("inputs") != _bigwig_input_signature(inputs):
            return False
        if int(payload.get("output_size_bytes", -1)) != output.stat().st_size:
            return False
        _verify_bigwig(output, [chrom for chrom, _size in chrom_sizes])
        return True
    except Exception:
        return False


def _validate_staged_bedgraph(
    bedgraph: Path,
    source_bigwig: Path,
) -> dict[str, object]:
    metadata_path = bedgraph.with_name(bedgraph.name + ".complete.json")
    if not bedgraph.is_file():
        raise FileNotFoundError(f"Staged bedGraph not found: {bedgraph}")
    if not metadata_path.is_file():
        raise FileNotFoundError(
            f"Completion metadata not found for staged bedGraph: {metadata_path}"
        )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if not metadata.get("complete"):
        raise ValueError(f"Staged bedGraph is not marked complete: {bedgraph}")
    if not metadata.get("sorted") or not metadata.get("nonoverlapping"):
        raise ValueError(f"Staged bedGraph is not validated as sorted: {bedgraph}")
    recorded_source = Path(str(metadata.get("source_bigwig", ""))).resolve()
    if recorded_source != source_bigwig.resolve():
        raise ValueError(
            f"Staged bedGraph source mismatch for {bedgraph}: "
            f"expected {source_bigwig.resolve()}, observed {recorded_source}"
        )
    recorded_size = int(metadata.get("size_bytes", -1))
    actual_size = bedgraph.stat().st_size
    if recorded_size != actual_size:
        raise ValueError(
            f"Staged bedGraph size changed after validation: {bedgraph} "
            f"({recorded_size} recorded, {actual_size} observed)"
        )
    return metadata


def _combine_bigwig_via_bedgraph(
    inputs: Sequence[Path],
    output: Path,
    chrom_sizes: Sequence[tuple[str, int]],
    temp_dir: Path,
    *,
    logger: _CombineLogger | None = None,
    chunk_size: int = BIGWIG_COMBINE_CHUNK_SIZE,
    staged_bedgraphs: Sequence[Path] | None = None,
) -> None:
    executable = shutil.which("bedGraphToBigWig")
    if executable is None:
        raise RuntimeError(
            "bedGraphToBigWig is required for bedGraph combination mode. "
            "Install ucsc-bedgraphtobigwig. Staged bedGraphs have been retained."
        )
    chrom_order = [chrom for chrom, _ in chrom_sizes]
    chrom_lengths = dict(chrom_sizes)
    temp_dir.mkdir(parents=True, exist_ok=True)
    individual: list[Path] = []
    staged_metadata: list[Path] = []
    started = time.monotonic()
    if staged_bedgraphs is not None and len(staged_bedgraphs) != len(inputs):
        raise ValueError("Each BigWig input requires one staged bedGraph")
    if logger is not None:
        total_input_bytes = sum(path.stat().st_size for path in inputs if path.exists())
        source_kind = "staged bedGraphs" if staged_bedgraphs is not None else "BigWig conversion"
        logger.log(
            f"Preparing BigWig {output}: {len(inputs)} input file(s), "
            f"{_format_bytes(total_input_bytes)} total; method: {source_kind}; "
            f"temporary files: {temp_dir}"
        )
    for index, path in enumerate(inputs):
        if logger is not None:
            size = path.stat().st_size if path.exists() else 0
            logger.log(
                f"Input {index + 1}/{len(inputs)}: {path} ({_format_bytes(size)})"
            )
        if staged_bedgraphs is not None:
            bedgraph = Path(staged_bedgraphs[index])
            metadata = _validate_staged_bedgraph(bedgraph, path)
            individual.append(bedgraph)
            staged_metadata.append(bedgraph.with_name(bedgraph.name + ".complete.json"))
            if logger is not None:
                logger.log(
                    f"Using validated staged bedGraph {bedgraph} "
                    f"({_format_bytes(bedgraph.stat().st_size)}; "
                    f"{int(metadata.get('records', 0)):,} records)"
                )
            continue
        bedgraph = temp_dir / f"{output.stem}.{index:04d}.{path.stem}.bedGraph"
        interval_count, queried_chunks = _bigwig_to_bedgraph(
            path,
            bedgraph,
            chrom_order,
            chrom_lengths=chrom_lengths,
            chunk_size=chunk_size,
            logger=logger,
        )
        individual.append(bedgraph)
        if logger is not None:
            logger.log(
                f"Wrote temporary bedGraph {bedgraph} "
                f"({_format_bytes(bedgraph.stat().st_size)}; {interval_count:,} intervals; "
                f"{queried_chunks:,} queried chunks)"
            )
    combined_bedgraph = temp_dir / f"{output.stem}.combined.bedGraph"
    with combined_bedgraph.open("wb") as destination:
        for path in individual:
            with path.open("rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
    chrom_sizes_path = temp_dir / f"{output.stem}.chrom.sizes"
    _write_chrom_sizes(chrom_sizes_path, chrom_sizes)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial_output = output.with_name(output.name + ".partial")
    partial_output.unlink(missing_ok=True)
    if logger is not None:
        logger.log(
            f"Running bedGraphToBigWig for {output}; combined bedGraph size "
            f"{_format_bytes(combined_bedgraph.stat().st_size)}"
        )
    _run_logged_subprocess(
        [executable, str(combined_bedgraph), str(chrom_sizes_path), str(partial_output)],
        logger,
    )
    _verify_bigwig(partial_output, chrom_order)
    os.replace(partial_output, output)
    _write_combined_bigwig_marker(
        output,
        inputs=inputs,
        chrom_sizes=chrom_sizes,
        method="bedgraph",
    )
    if logger is not None:
        logger.log(
            f"Verified {output} ({_format_bytes(output.stat().st_size)}; "
            f"elapsed {time.monotonic() - started:.1f} s)"
        )
    # Cleanup occurs only after successful conversion and verification.
    generated_individual = [] if staged_bedgraphs is not None else individual
    for path in [*generated_individual, combined_bedgraph, chrom_sizes_path]:
        path.unlink(missing_ok=True)
    if staged_bedgraphs is not None:
        for path in [*individual, *staged_metadata]:
            path.unlink(missing_ok=True)
    try:
        temp_dir.rmdir()
    except OSError:
        pass


def _source_bigwig_map(
    inputs: Sequence[Path],
    chrom_sizes: Sequence[tuple[str, int]],
) -> dict[str, tuple[Path, str, int]]:
    if pyBigWig is None:
        raise RuntimeError("Direct BigWig combination requires pyBigWig")
    expected = [chrom for chrom, _size in chrom_sizes]
    mapping: dict[str, tuple[Path, str, int]] = {}
    for path in inputs:
        handle = pyBigWig.open(str(path))
        if handle is None:
            raise RuntimeError(f"Could not open BigWig: {path}")
        try:
            available = handle.chroms()
            names = list(available)
            for chrom in expected:
                try:
                    source_chrom = resolve_contig_name(
                        chrom, names, source_label=f"BigWig {path}"
                    )
                except KeyError:
                    continue
                if chrom in mapping:
                    raise ValueError(
                        f"More than one per-contig BigWig supplies chromosome {chrom}: "
                        f"{mapping[chrom][0]} and {path}"
                    )
                mapping[chrom] = (path, source_chrom, int(available[source_chrom]))
        finally:
            handle.close()
    missing = [chrom for chrom in expected if chrom not in mapping]
    if missing:
        raise ValueError(
            "No per-contig BigWig supplied chromosome(s): " + ", ".join(missing)
        )
    return mapping


def _combine_bigwig_direct(
    inputs: Sequence[Path],
    output: Path,
    chrom_sizes: Sequence[tuple[str, int]],
    *,
    logger: _CombineLogger | None = None,
    chunk_size: int = BIGWIG_COMBINE_CHUNK_SIZE,
) -> None:
    if pyBigWig is None:
        raise RuntimeError("Direct BigWig combination requires pyBigWig")
    if chunk_size < 1:
        raise ValueError("BigWig combine chunk size must be positive")
    started = time.monotonic()
    mapping = _source_bigwig_map(inputs, chrom_sizes)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial_output = output.with_name(output.name + ".partial")
    partial_output.unlink(missing_ok=True)
    writer = pyBigWig.open(str(partial_output), "w")
    if writer is None:
        raise RuntimeError(f"Could not create combined BigWig: {partial_output}")
    total_intervals = 0
    try:
        writer.addHeader([(chrom, int(size)) for chrom, size in chrom_sizes])
        for chrom_index, (chrom, requested_length) in enumerate(chrom_sizes, start=1):
            source_path, source_chrom, source_length = mapping[chrom]
            length = min(int(requested_length), int(source_length))
            source = pyBigWig.open(str(source_path))
            if source is None:
                raise RuntimeError(f"Could not open BigWig: {source_path}")
            chrom_intervals = 0
            chrom_chunks = math.ceil(length / chunk_size) if length else 0
            if logger is not None:
                logger.log(
                    f"Direct stream {chrom_index}/{len(chrom_sizes)}: {source_path} -> "
                    f"{chrom} ({length:,} bp; {chrom_chunks:,} chunks)"
                )
            try:
                for chunk_index, chunk_start in enumerate(
                    range(0, length, chunk_size), start=1
                ):
                    chunk_end = min(length, chunk_start + chunk_size)
                    intervals = source.intervals(
                        source_chrom, chunk_start, chunk_end
                    ) or ()
                    owned = [
                        (int(start), int(end), float(value))
                        for start, end, value in intervals
                        if chunk_start <= int(start) < chunk_end
                    ]
                    if owned:
                        starts, ends, values = zip(*owned)
                        writer.addEntries(
                            [chrom] * len(owned),
                            list(starts),
                            ends=list(ends),
                            values=list(values),
                        )
                        chrom_intervals += len(owned)
                        total_intervals += len(owned)
                    if logger is not None and (
                        chunk_index == chrom_chunks or chunk_index % 100 == 0
                    ):
                        rss = _rss_mb()
                        memory = "" if rss is None else f"; RSS {rss:.1f} MiB"
                        logger.log(
                            f"  {chrom}: {chunk_index:,}/{chrom_chunks:,} chunks; "
                            f"{chrom_intervals:,} intervals copied{memory}"
                        )
            finally:
                source.close()
    finally:
        writer.close()
    _verify_bigwig(partial_output, [chrom for chrom, _size in chrom_sizes])
    os.replace(partial_output, output)
    _write_combined_bigwig_marker(
        output,
        inputs=inputs,
        chrom_sizes=chrom_sizes,
        method="direct",
    )
    if logger is not None:
        logger.log(
            f"Verified direct combined BigWig {output} "
            f"({_format_bytes(output.stat().st_size)}; {total_intervals:,} intervals; "
            f"elapsed {time.monotonic() - started:.1f} s)"
        )


def _combine_bigwig(
    inputs: Sequence[Path],
    output: Path,
    chrom_sizes: Sequence[tuple[str, int]],
    temp_dir: Path,
    *,
    logger: _CombineLogger | None = None,
    chunk_size: int = BIGWIG_COMBINE_CHUNK_SIZE,
    method: str = "bedgraph",
    staged_bedgraphs: Sequence[Path] | None = None,
) -> None:
    selected = str(method).lower()
    if selected == "bedgraphs":
        selected = "bedgraph"
    if selected == "direct":
        _combine_bigwig_direct(
            inputs,
            output,
            chrom_sizes,
            logger=logger,
            chunk_size=chunk_size,
        )
        return
    if selected == "bedgraph":
        _combine_bigwig_via_bedgraph(
            inputs,
            output,
            chrom_sizes,
            temp_dir,
            logger=logger,
            chunk_size=chunk_size,
            staged_bedgraphs=staged_bedgraphs,
        )
        return
    raise ValueError(f"Unknown BigWig combine method: {method}")



def _combine_bams(inputs: Sequence[Path], output: Path) -> None:
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("BAM combination requires pysam") from exc
    output.parent.mkdir(parents=True, exist_ok=True)
    pysam.merge("-f", "-o", str(output), *[str(path) for path in inputs])
    pysam.index(str(output))


def _combine_bigbed(
    inputs: Sequence[Path],
    output_bed: Path,
    chrom_sizes: Sequence[tuple[str, int]] | None = None,
) -> None:
    if pyBigWig is None:
        raise RuntimeError("bigBed combination requires pyBigWig or retained BED output")
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    with output_bed.open("w", encoding="utf-8") as output:
        for path in inputs:
            handle = pyBigWig.open(str(path))
            if handle is None:
                raise RuntimeError(f"Could not open bigBed: {path}")
            try:
                canonical_names = [name for name, _length in chrom_sizes] if chrom_sizes else []
                for source_chrom, size in handle.chroms().items():
                    chrom = source_chrom
                    if canonical_names:
                        try:
                            chrom = resolve_contig_name(
                                source_chrom,
                                canonical_names,
                                source_label="combined chromosome sizes",
                            )
                        except KeyError as exc:
                            raise ValueError(str(exc)) from exc
                    entries = handle.entries(source_chrom, 0, int(size)) or []
                    for start, end, rest in entries:
                        fields = [chrom, str(start), str(end)]
                        if rest:
                            fields.extend(str(rest).split("\t"))
                        output.write("\t".join(fields) + "\n")
            finally:
                handle.close()
    if chrom_sizes and _interval_file_has_inversion(output_bed, chrom_sizes):
        _replace_with_sorted_interval_file(output_bed, chrom_sizes)


def combine_run(
    input_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    chrom_sizes_path: str | Path | None = None,
    combine_tracks: bool = True,
    bigwig_method: str = "direct",
    cores: int = 1,
    streaming_cores: int | None = None,
    indexed_cores: int = 1,
    bigwig_chunk_size: int = BIGWIG_COMBINE_CHUNK_SIZE,
    force: bool = False,
) -> dict[str, object]:
    """Combine one manifest-backed multicontig run."""
    cores = int(cores or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    streaming_cores = int(streaming_cores or cores)
    indexed_cores = int(indexed_cores or 1)
    bigwig_chunk_size = int(bigwig_chunk_size)
    if streaming_cores < 1:
        raise ValueError("--streaming-combine-cores must be at least 1")
    if indexed_cores < 1:
        raise ValueError("--indexed-combine-cores must be at least 1")
    if bigwig_chunk_size < 1:
        raise ValueError("--combine-chunk-bp must be at least 1")
    root = Path(input_dir).resolve()
    manifest = _read_manifest(root)
    if manifest.get("combine_strategy") == "aggregate_matrix":
        from nucleosuite.aggregate_parallel import combine_aggregate_manifest
        resolved_output = Path(output_dir).resolve() if output_dir else None
        return combine_aggregate_manifest(root, manifest, output_dir=resolved_output)
    if manifest.get("combine_strategy") == "rerun_namespace":
        from nucleosuite.partitioned import rerun_manifest
        resolved_output = Path(output_dir).resolve() if output_dir else None
        return rerun_manifest(root, manifest, output_dir=resolved_output)
    if manifest.get("combine_strategy") == "directory_tree":
        settings = dict(manifest.get("directory_tree") or {})
        resolved_output = Path(output_dir).resolve() if output_dir else Path(
            str(manifest.get("combined_dir", root / "combined"))
        ).resolve()
        tree_cores = int(settings.get("cores", cores) or cores)
        return combine_directory_trees(
            settings.get("per_contig_dirs")
            or [entry["directory"] for entry in manifest["per_contig"]],
            resolved_output,
            chrom_sizes=_chrom_sizes_from_manifest(manifest),
            include_roots=settings.get("include_roots"),
            exclude_parts=tuple(settings.get("exclude_parts") or ()),
            sample_name=settings.get("sample_name"),
            bigwig_method=str(settings.get("bigwig_method", bigwig_method)),
            strict_complete=bool(settings.get("strict_complete", True)),
            combine_tracks=bool(settings.get("combine_tracks", combine_tracks)),
            cores=streaming_cores,
            streaming_cores=streaming_cores,
            indexed_cores=indexed_cores,
            bigwig_chunk_size=bigwig_chunk_size,
            force=force,
        )
    output_root = Path(output_dir).resolve() if output_dir else Path(str(manifest.get("combined_dir", root / "combined"))).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    logger = _CombineLogger(output_root / "combine_steps.log")
    logger.log(
        f"Starting manifest combination from {root}; streaming workers="
        f"{streaming_cores}; indexed workers={indexed_cores}; BigWig chunk="
        f"{bigwig_chunk_size:,} bp; force={bool(force)}"
    )
    combined_name = str(manifest.get("combined_name") or manifest.get("base_name") or "combined")
    combined_prefix = output_root / combined_name
    groups = _output_groups(manifest)

    if chrom_sizes_path:
        chrom_sizes = list(read_chrom_sizes_source(chrom_sizes_path))
    else:
        chrom_sizes = _chrom_sizes_from_manifest(manifest)
    if not chrom_sizes:
        raise ValueError("Chromosome sizes are required to combine BigWig outputs")

    written: list[str] = []
    warnings: list[str] = []
    bigwig_groups: list[tuple[str, list[Path]]] = []
    bigbed_groups: list[tuple[str, list[Path]]] = []

    # Sufficient-statistic count files must be written before derived profiles.
    def _suffix_priority(value: str) -> tuple[int, str]:
        if value.endswith("_dinuc_profile_counts.tsv") or value.endswith("_length_counts.tsv") or value.endswith("_ww_type_by_length.tsv"):
            return (0, value)
        if value.endswith("_summary.tsv"):
            return (2, value)
        return (1, value)
    suffixes = sorted(groups, key=_suffix_priority)

    def combine_one_suffix(suffix: str) -> tuple[list[str], list[str], tuple[str, list[Path]] | None, tuple[str, list[Path]] | None]:
        inputs = groups[suffix]
        output = _combined_path(combined_prefix, suffix)
        local_written: list[str] = []
        logger.log(f"START table/interval: {suffix}; inputs={len(inputs)}")
        try:
            if suffix.endswith((".bw", ".bigWig")):
                logger.log(f"QUEUE BigWig: {suffix}")
                return local_written, [], (suffix, inputs), None
            if suffix.endswith(".bam"):
                _combine_bams(inputs, output)
                local_written.extend([str(output), str(output) + ".bai"])
                return local_written, [], None, None
            if suffix.endswith(".bam.bai"):
                return local_written, [], None, None
            if suffix.endswith((".bb", ".bigBed")):
                return local_written, [], None, (suffix, inputs)
            if suffix.endswith((".png", ".svg")):
                return local_written, [], None, None
            if suffix.endswith((".bed", ".bed.gz")):
                _concatenate_intervals(inputs, output, chrom_sizes)
            elif suffix.endswith(".wig.gz"):
                _concatenate_text(inputs, output)
            elif suffix.endswith("_fragment_length_counts.tsv") or suffix.endswith(".fragment_length_counts.tsv"):
                _combine_fragment_length_counts(inputs, output)
            elif suffix.endswith(".randomization_qc.tsv"):
                _combine_randomization_qc(inputs, output)
            elif suffix.endswith(".merge_summary.tsv"):
                _combine_metric_value_summary(
                    inputs, output,
                    additive_metrics={"read_pairs_written", "reads_written", "output_bams"},
                )
            elif suffix.endswith("_fragment_summary.tsv") or suffix.endswith(".fragments.summary.tsv"):
                _combine_fragment_summary(inputs, output)
            elif suffix.endswith("_ww_type_summary.tsv"):
                _combine_ww_summary(inputs, output)
            elif suffix.endswith("_ww_type_by_length.tsv"):
                _combine_ww_type_by_length(inputs, output)
            elif suffix.endswith("_dinuc_profile_counts.tsv"):
                _combine_dinuc_counts(inputs, output)
            elif suffix.endswith("_runs.tsv.gz"):
                _concatenate_text(inputs, output, keep_one_header=True)
            elif suffix.endswith("_dinuc_profile.tsv"):
                counts_suffix = suffix[:-4] + "_counts.tsv"
                counts_path = _combined_path(combined_prefix, counts_suffix)
                if counts_path.exists():
                    fields, _ = _read_tsv(inputs[0])
                    _profile_from_counts(
                        counts_path,
                        output,
                        any(name.endswith("_frac") for name in fields),
                    )
                else:
                    _combine_generic_tsv(inputs, output)
            elif suffix.endswith("_summary.tsv") and any(
                path.name.endswith("_summary.tsv") for path in inputs
            ) and any(
                (path.parent / path.name.replace("_summary.tsv", "_length_counts.tsv")).exists()
                for path in inputs
            ):
                counts_path = output.with_name(
                    output.name.replace("_summary.tsv", "_length_counts.tsv")
                )
                _combine_positive_summary(inputs, output, counts_path)
            elif suffix.endswith(("_DAC_summary.tsv", "_DAC_run_summary.tsv")):
                _combine_run_summary(inputs, output, kind="dac", combined_prefix=combined_prefix)
            elif suffix.endswith(("_DCC_summary.tsv", "_DCC_run_summary.tsv")):
                _combine_run_summary(inputs, output, kind="dcc", combined_prefix=combined_prefix)
            elif "_DAC_" in suffix and suffix.endswith(".tsv"):
                _combine_dac(inputs, output)
            elif "_DCC" in suffix and suffix.endswith(".tsv") and "shift_summary" not in suffix:
                _combine_dcc(inputs, output)
            elif suffix.endswith(".tsv"):
                _combine_generic_tsv(inputs, output)
            else:
                return local_written, [], None, None
            local_written.append(str(output))
            logger.log(
                f"PASS table/interval: {suffix}; output={output}; "
                f"size={_format_bytes(output.stat().st_size)}"
            )
            return local_written, [], None, None
        except Exception as error:
            logger.log(f"FAIL table/interval: {suffix}: {error}")
            return [], [f"{suffix}: {error}"], None, None

    for stage_priority in sorted({_suffix_priority(value)[0] for value in suffixes}):
        stage = [value for value in suffixes if _suffix_priority(value)[0] == stage_priority]
        if not stage:
            continue
        with ThreadPoolExecutor(max_workers=min(streaming_cores, len(stage))) as executor:
            futures = {executor.submit(combine_one_suffix, suffix): suffix for suffix in stage}
            for future in as_completed(futures):
                local_written, local_warnings, bigwig, bigbed = future.result()
                written.extend(local_written)
                warnings.extend(local_warnings)
                if bigwig is not None:
                    bigwig_groups.append(bigwig)
                if bigbed is not None:
                    bigbed_groups.append(bigbed)
    bigwig_groups.sort(key=lambda item: item[0])
    bigbed_groups.sort(key=lambda item: item[0])

    # Recreate the fragment-length plot after counts have been summed.
    for suffix in groups:
        if suffix.endswith("_fragment_length_counts.tsv") or suffix.endswith(".fragment_length_counts.tsv"):
            counts_path = _combined_path(combined_prefix, suffix)
            if counts_path.exists():
                try:
                    from nucleosuite.profile_plots import plot_count_profile
                    png_name = counts_path.name.replace("_fragment_length_counts.tsv", "_fragment_length_distribution.png").replace(".fragment_length_counts.tsv", ".fragment_length_distribution.png")
                    png = counts_path.with_name(png_name)
                    png = plot_count_profile(counts_path, png, x_column="fragment_length", y_column="count", xlabel="Fragment length (bp)", ylabel="Fragment count", title="Fragment-length distribution")
                    written.append(str(png))
                except Exception as error:
                    warnings.append(f"fragment-length plot: {error}")

    for suffix in groups:
        if suffix.endswith("_length_counts.tsv") and any(path.name.endswith("_summary.tsv") for path in groups.get(suffix.replace("_length_counts.tsv", "_summary.tsv"), [])):
            table = _combined_path(combined_prefix, suffix)
            if table.exists():
                try:
                    from nucleosuite.positive_runs import plot_distribution
                    counts = Counter()
                    _, rows = _read_tsv(table)
                    for row in rows:
                        counts[int(float(row["run_length_bp"]))] += int(float(row["count"]))
                    png = table.with_name(table.name.replace("_length_counts.tsv", "_run_length_distribution.png"))
                    png = plot_distribution(png, counts, normalization="count", plot_x_max=550, title="Positive run lengths")
                    written.append(str(png))
                except Exception as error:
                    warnings.append(f"positive-run plot: {error}")

    for suffix in groups:
        if suffix.endswith(".relocation_distances.tsv"):
            table = _combined_path(combined_prefix, suffix)
            if table.exists():
                try:
                    from nucleosuite.profile_plots import plot_count_profile
                    png = table.with_suffix(".png")
                    png = plot_count_profile(
                        str(table), str(png), x_column="relocation_bp", y_column="count",
                        xlabel="Fragment relocation (bp)", ylabel="Fragment count",
                        title="Randomized fragment relocation distances", vertical_zero=True,
                    )
                    written.append(str(png))
                except Exception as error:
                    warnings.append(f"relocation plot: {error}")

    # bigBed conversion is handled before BigWigs, but after all text outputs.
    def combine_one_bigbed(item: tuple[str, list[Path]]) -> tuple[list[str], list[str]]:
        suffix, inputs = item
        bed_suffix = suffix.rsplit(".", 1)[0] + ".bed"
        bed_path = _combined_path(combined_prefix, bed_suffix)
        output_bigbed = _combined_path(combined_prefix, suffix)
        empty_marker = Path(str(output_bigbed) + ".empty")
        logger.log(f"START bigBed: {suffix}; inputs={len(inputs)}")
        try:
            if not inputs:
                bed_path.parent.mkdir(parents=True, exist_ok=True)
                bed_path.write_text("", encoding="utf-8")
                output_bigbed.unlink(missing_ok=True)
                empty_marker.write_text(
                    "status\tempty\nrecord_count\t0\n", encoding="utf-8"
                )
                return [str(bed_path), str(empty_marker)], []
            empty_marker.unlink(missing_ok=True)
            if not bed_path.exists():
                _combine_bigbed(inputs, bed_path, chrom_sizes)
            logger.log(f"PASS bigBed source merge: {suffix}; output={bed_path}")
            return [str(bed_path)], []
        except Exception as error:
            logger.log(f"FAIL bigBed: {suffix}: {error}")
            return [], [f"{suffix}: {error}"]

    if bigbed_groups:
        with ThreadPoolExecutor(max_workers=min(indexed_cores, len(bigbed_groups))) as executor:
            for future in as_completed(
                [executor.submit(combine_one_bigbed, item) for item in bigbed_groups]
            ):
                local_written, local_warnings = future.result()
                written.extend(local_written)
                warnings.extend(local_warnings)

    # BigWig combination is intentionally the final stage.
    if combine_tracks:
        temp_dir = output_root / "temporary_bedgraph_combine"

        def combine_one_bigwig(item: tuple[str, list[Path]]) -> tuple[list[str], list[str]]:
            suffix, inputs = item
            output = _combined_path(combined_prefix, suffix)
            try:
                selected_method = "bedgraph" if bigwig_method == "bedgraphs" else bigwig_method
                logger.log(f"START BigWig: {suffix}; inputs={len(inputs)}; method={selected_method}")
                if not force and _combined_bigwig_is_reusable(
                    output,
                    inputs=inputs,
                    chrom_sizes=chrom_sizes,
                    method=selected_method,
                ):
                    logger.log(f"REUSE BigWig: {suffix}; validated completion marker")
                    return [str(output)], []
                _combine_bigwig(
                    inputs,
                    output,
                    chrom_sizes,
                    temp_dir / output.stem,
                    method=selected_method,
                    chunk_size=bigwig_chunk_size,
                    logger=logger,
                )
                logger.log(f"PASS BigWig: {suffix}; output={output}")
                return [str(output)], []
            except Exception as error:
                logger.log(f"FAIL BigWig: {suffix}: {error}")
                return [], [f"{suffix}: {error}"]

        if bigwig_groups:
            with ThreadPoolExecutor(max_workers=min(indexed_cores, len(bigwig_groups))) as executor:
                for future in as_completed(
                    [executor.submit(combine_one_bigwig, item) for item in bigwig_groups]
                ):
                    local_written, local_warnings = future.result()
                    written.extend(local_written)
                    warnings.extend(local_warnings)
    elif bigwig_groups:
        warnings.append("BigWig combination was skipped by --skip-tracks")

    result = {
        "input_dir": str(root),
        "output_dir": str(output_root),
        "combined_prefix": str(combined_prefix),
        "written": written,
        "warnings": warnings,
        "bigwig_stage_last": True,
        "cores": cores,
        "streaming_combine_cores": streaming_cores,
        "indexed_combine_cores": indexed_cores,
        "bigwig_combine_chunk_bp": bigwig_chunk_size,
        "combine_log": str(logger.path),
        "force": bool(force),
    }
    with (output_root / "combine_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    logger.log(
        f"Combination finished: outputs={len(written):,}; warnings={len(warnings):,}"
    )
    logger.close()
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite combine",
        description=(
            "Combine per-contig NucleoSuite outputs. Raw counts and opportunity "
            "denominators are summed before normalized values and percentages are recalculated."
        ),
        formatter_class=NucleoSuiteHelpFormatter,
    )
    parser.add_argument("--input-dir", required=True, help=f"Multicontig run directory containing {MANIFEST_NAME}.")
    parser.add_argument("--output-dir", help="Combined output directory; default is recorded in the manifest.")
    parser.add_argument("--chrom-sizes", help="Optional chromosome-size table, BAM or CRAM override for BigWig creation.")
    parser.add_argument("--skip-tracks", action="store_true", help="Combine tabular and interval outputs but defer BigWig creation.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute outputs even when compatible completion checkpoints exist.",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help=(
            "Default worker budget for memory-light streaming combines. It does "
            "not increase indexed BigWig/BigBed concurrency (default: 1)."
        ),
    )
    parser.add_argument(
        "--streaming-combine-cores",
        type=int,
        default=None,
        help=(
            "Concurrent memory-light table and interval combines. Default: --cores."
        ),
    )
    parser.add_argument(
        "--indexed-combine-cores",
        type=int,
        default=1,
        help=(
            "Concurrent BigWig/BigBed combines. This finalization-sensitive budget "
            "is independent of --cores (default: 1)."
        ),
    )
    parser.add_argument(
        "--combine-chunk-bp",
        type=int,
        default=BIGWIG_COMBINE_CHUNK_SIZE,
        help="Genomic query chunk used while combining BigWigs (default: 100000 bp).",
    )
    parser.add_argument(
        "--bigwig-method",
        choices=("direct", "bedgraph", "bedgraphs"),
        default="direct",
        help=(
            "Combine BigWigs directly with pyBigWig or through temporary bedGraphs. "
            "Suite bedGraph mode normally supplies pre-staged bedGraphs."
        ),
    )
    from nucleosuite.plotting import add_plotting_arguments
    add_plotting_arguments(parser)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = combine_run(
        args.input_dir,
        output_dir=args.output_dir,
        chrom_sizes_path=args.chrom_sizes,
        combine_tracks=not args.skip_tracks,
        bigwig_method=args.bigwig_method,
        cores=args.cores,
        streaming_cores=args.streaming_combine_cores,
        indexed_cores=args.indexed_combine_cores,
        bigwig_chunk_size=args.combine_chunk_bp,
        force=args.force,
    )
    print(f"Combined outputs: {result['output_dir']}")
    for warning in result["warnings"]:
        print(f"WARNING: {warning}")
    return 0 if not result["warnings"] else 1


if __name__ == "__main__":
    raise SystemExit(main())


def combine_directory_trees(
    per_contig_dirs: Sequence[str | Path],
    output_dir: str | Path,
    *,
    chrom_sizes: Sequence[tuple[str, int]],
    include_roots: Sequence[str] | None = None,
    exclude_parts: Sequence[str] = (),
    sample_name: str | None = None,
    bigwig_method: str = "direct",
    strict_complete: bool = False,
    combine_tracks: bool = True,
    cores: int = 1,
    streaming_cores: int | None = None,
    indexed_cores: int = 1,
    bigwig_chunk_size: int = BIGWIG_COMBINE_CHUNK_SIZE,
    force: bool = False,
) -> dict[str, object]:
    """Combine matching files from complete per-contig workflow directories.

    This is used by ``mnase-suite``. Only directories listed in ``include_roots``
    are considered. BigWig files are queued until all tabular and interval files
    have been written.
    """
    cores = int(cores or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    streaming_cores = int(streaming_cores or cores)
    indexed_cores = int(indexed_cores or 1)
    bigwig_chunk_size = int(bigwig_chunk_size)
    if streaming_cores < 1:
        raise ValueError("--streaming-combine-cores must be at least 1")
    if indexed_cores < 1:
        raise ValueError("--indexed-combine-cores must be at least 1")
    if bigwig_chunk_size < 1:
        raise ValueError("--combine-chunk-bp must be at least 1")
    roots = [Path(path).resolve() for path in per_contig_dirs]
    missing_roots = [root for root in roots if not root.is_dir()]
    if missing_roots:
        raise FileNotFoundError(
            "Missing per-contig workflow directories: "
            + ", ".join(str(path) for path in missing_roots)
        )
    destination_root = Path(output_dir).resolve()
    destination_root.mkdir(parents=True, exist_ok=True)
    combine_log_name = (
        f"{sample_name}_combine_steps.log"
        if sample_name and sample_name.endswith("_randomized_control")
        else "combine_steps.log"
    )
    logger = _CombineLogger(destination_root / combine_log_name)
    logger.log("=" * 72)
    logger.log(
        f"Starting directory-tree combination from {len(roots)} per-contig directories "
        f"into {destination_root}"
    )
    logger.log(
        "Chromosome order: " + ", ".join(f"{chrom} ({size:,} bp)" for chrom, size in chrom_sizes)
    )
    logger.log(
        f"Streaming combine workers: {streaming_cores}; indexed combine workers: "
        f"{indexed_cores}; BigWig chunk: {bigwig_chunk_size:,} bp"
    )
    include = set(include_roots or ())
    groups: dict[Path, list[Path]] = defaultdict(list)
    group_sources: dict[Path, set[Path]] = defaultdict(set)
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            relative = path.relative_to(root)
            if sample_name:
                worker_prefix = _worker_sample_prefix(sample_name, root.name)
                if relative.name.startswith(worker_prefix):
                    relative = relative.with_name(sample_name + relative.name[len(worker_prefix):])
            if include and relative.parts and relative.parts[0] not in include:
                continue
            if any(part in exclude_parts for part in relative.parts):
                continue
            if relative.parts and relative.parts[0] in {"logs", ".done"}:
                continue
            if relative.name.endswith((".bb.empty", ".bigBed.empty")):
                virtual = relative.with_name(relative.name[:-6])
                groups.setdefault(virtual, [])
                group_sources[virtual].add(root)
                continue
            groups[relative].append(path)
            group_sources[relative].add(root)

    incomplete_groups = {
        relative: sorted(set(roots) - sources)
        for relative, sources in group_sources.items()
        if len(sources) != len(roots)
    }
    if incomplete_groups:
        details = "; ".join(
            f"{relative}: missing {','.join(path.name for path in missing)}"
            for relative, missing in sorted(incomplete_groups.items(), key=lambda item: str(item[0]))
        )
        logger.log("Incomplete output groups detected: " + details)
        if strict_complete:
            logger.close()
            raise ValueError(
                "Refusing to create partial combined outputs because one or more "
                "per-contig contributions are missing: " + details
            )

    unsupported_groups = [
        relative for relative in groups
        if relative.name.endswith(".tsv")
        and not _supported_tree_tsv_name(relative.name)
    ]
    if unsupported_groups and strict_complete:
        logger.close()
        raise ValueError(
            "No explicit combination strategy is registered for: "
            + ", ".join(str(path) for path in sorted(unsupported_groups))
        )

    logger.log(
        f"Discovered {len(groups):,} matching output group(s) across "
        f"{sum(len(paths) for paths in groups.values()):,} files"
    )
    written: list[str] = []
    warnings: list[str] = []
    bigwigs: list[tuple[Path, list[Path]]] = []
    bigbeds: list[tuple[Path, list[Path]]] = []

    def priority(relative: Path) -> tuple[int, str]:
        name = relative.name
        if name.endswith("_dinuc_profile_counts.tsv") or name.endswith("_length_counts.tsv") or name.endswith("fragment_length_counts.tsv") or name.endswith("_ww_type_by_length.tsv"):
            return (0, str(relative))
        if name.endswith("_summary.tsv"):
            return (2, str(relative))
        return (1, str(relative))

    ordered_groups = sorted(groups, key=priority)

    def process_tree_group(relative: Path) -> tuple[list[str], list[str], tuple[Path, list[Path]] | None, tuple[Path, list[Path]] | None]:
        inputs = groups[relative]
        output = destination_root / relative
        if relative in incomplete_groups:
            return [], [f"{relative}: skipped because not every selected chromosome contributed"], None, None
        name = relative.name
        if name.endswith((".bw", ".bigWig")):
            return [], [], (relative, inputs), None
        if name.endswith((".bb", ".bigBed")):
            return [], [], None, (relative, inputs)
        logger.log(f"Combining table/interval: {relative} ({len(inputs)} input file(s))")
        try:
            if name.endswith((".png", ".svg", ".xlsx", ".bai")):
                return [], [], None, None
            if name.endswith(".bam"):
                _combine_bams(inputs, output)
                local_written = [str(output), str(output) + ".bai"]
            elif name.endswith((".bed", ".bed.gz")):
                _concatenate_intervals(inputs, output, chrom_sizes)
                local_written = [str(output)]
            elif name.endswith((".wig", ".wig.gz")):
                _concatenate_text(inputs, output)
                local_written = [str(output)]
            elif _is_track_metadata_name(name) or name.endswith("report.tsv"):
                _combine_track_metadata_table(
                    inputs,
                    output,
                    roots=roots,
                    destination_root=destination_root,
                    sample_name=sample_name,
                )
                local_written = [str(output)]
            elif name.endswith("_fragment_length_counts.tsv") or name.endswith(".fragment_length_counts.tsv"):
                _combine_fragment_length_counts(inputs, output)
                local_written = [str(output)]
            elif name.endswith(".randomization_qc.tsv"):
                _combine_randomization_qc(inputs, output)
                local_written = [str(output)]
            elif name.endswith(".relocation_distances.tsv"):
                _combine_generic_tsv(inputs, output)
                local_written = [str(output)]
            elif name.endswith(".merge_summary.tsv"):
                _combine_metric_value_summary(
                    inputs,
                    output,
                    additive_metrics={"read_pairs_written", "reads_written", "output_bams"},
                )
                local_written = [str(output)]
            elif name.endswith("_fragment_summary.tsv") or name.endswith(".fragments.summary.tsv"):
                _combine_fragment_summary(inputs, output)
                local_written = [str(output)]
            elif name.endswith("_ww_type_summary.tsv"):
                _combine_ww_summary(inputs, output)
                local_written = [str(output)]
            elif name.endswith("_ww_type_by_length.tsv"):
                _combine_ww_type_by_length(inputs, output)
                local_written = [str(output)]
            elif name.endswith("_dinuc_profile_counts.tsv"):
                _combine_dinuc_counts(inputs, output)
                local_written = [str(output)]
            elif name.endswith("_dinuc_profile.tsv"):
                counts_path = output.with_name(output.name[:-4] + "_counts.tsv")
                if counts_path.exists():
                    fields, _ = _read_tsv(inputs[0])
                    _profile_from_counts(
                        counts_path,
                        output,
                        any(field.endswith("_frac") for field in fields),
                    )
                else:
                    _combine_generic_tsv(inputs, output)
                local_written = [str(output)]
            elif name.endswith("_runs.tsv.gz"):
                _concatenate_text(inputs, output, keep_one_header=True)
                local_written = [str(output)]
            elif name.endswith("_summary.tsv") and output.with_name(
                name.replace("_summary.tsv", "_length_counts.tsv")
            ).exists():
                _combine_positive_summary(
                    inputs,
                    output,
                    output.with_name(name.replace("_summary.tsv", "_length_counts.tsv")),
                )
                local_written = [str(output)]
            elif "_DAC_" in name and name.endswith(".tsv") and not name.endswith("_DAC_summary.tsv"):
                _combine_dac(inputs, output)
                local_written = [str(output)]
            elif "_DCC" in name and name.endswith(".tsv") and "shift_summary" not in name and not name.endswith("_DCC_summary.tsv"):
                _combine_dcc(inputs, output)
                local_written = [str(output)]
            elif name.endswith(".tsv"):
                message = (
                    f"{relative}: no explicit combination strategy; retained only "
                    "in per-contig outputs"
                )
                if strict_complete:
                    raise ValueError(message)
                logger.log("SKIP table: " + message)
                return [], [message], None, None
            else:
                return [], [], None, None
            if output.exists():
                logger.log(
                    f"Completed table/interval: {relative} "
                    f"({_format_bytes(output.stat().st_size)})"
                )
            return local_written, [], None, None
        except Exception as error:
            message = f"{relative}: {error}"
            logger.log(f"FAILED table/interval: {message}")
            return [], [message], None, None

    for stage_priority in sorted({priority(relative)[0] for relative in ordered_groups}):
        stage = [
            relative for relative in ordered_groups
            if priority(relative)[0] == stage_priority
        ]
        if not stage:
            continue
        logger.log(
            f"Starting table/interval priority stage {stage_priority} with "
            f"{len(stage):,} job(s) and up to {min(streaming_cores, len(stage))} worker(s)"
        )
        with ThreadPoolExecutor(max_workers=min(streaming_cores, len(stage))) as executor:
            futures = {
                executor.submit(process_tree_group, relative): relative
                for relative in stage
            }
            for future in as_completed(futures):
                local_written, local_warnings, bigwig, bigbed = future.result()
                written.extend(local_written)
                warnings.extend(local_warnings)
                if bigwig is not None:
                    bigwigs.append(bigwig)
                if bigbed is not None:
                    bigbeds.append(bigbed)
    bigwigs.sort(key=lambda item: str(item[0]))
    bigbeds.sort(key=lambda item: str(item[0]))

    # Recreate plots supported directly by sufficient-statistic tables.
    for relative in groups:
        output = destination_root / relative
        name = relative.name
        try:
            if name.endswith("_fragment_length_counts.tsv") and output.exists():
                from nucleosuite.profile_plots import plot_count_profile
                png = output.with_name(name.replace("_fragment_length_counts.tsv", "_fragment_length_distribution.png"))
                png = plot_count_profile(output, png, x_column="fragment_length", y_column="count", xlabel="Fragment length (bp)", ylabel="Fragment count", title="Fragment-length distribution")
                written.append(str(png))
            elif name.endswith(".relocation_distances.tsv") and output.exists():
                from nucleosuite.profile_plots import plot_count_profile
                png = output.with_suffix(".png")
                png = plot_count_profile(output, png, x_column="relocation_bp", y_column="count", xlabel="Fragment relocation (bp)", ylabel="Fragment count", title="Randomized fragment relocation distances", vertical_zero=True)
                written.append(str(png))
        except Exception as error:
            warnings.append(f"plot {relative}: {error}")

    # bigBed files are converted to BED first. If UCSC bedToBigBed is available,
    # recreate the requested combined bigBed from the concatenated BED.
    logger.log(
        f"Beginning bigBed stage: {len(bigbeds):,} file(s), up to "
        f"{min(indexed_cores, max(1, len(bigbeds)))} worker(s)"
    )

    def process_tree_bigbed(index: int, item: tuple[Path, list[Path]]) -> tuple[list[str], list[str]]:
        relative, inputs = item
        bed_relative = relative.with_suffix(".bed")
        bed_output = destination_root / bed_relative
        try:
            logger.log(f"START bigBed {index}/{len(bigbeds)}: {relative}")
            output = destination_root / relative
            empty_marker = Path(str(output) + ".empty")
            if not inputs:
                bed_output.parent.mkdir(parents=True, exist_ok=True)
                bed_output.write_text("", encoding="utf-8")
                output.unlink(missing_ok=True)
                empty_marker.write_text(
                    "status\tempty\nrecord_count\t0\n", encoding="utf-8"
                )
                logger.log(f"PASS empty bigBed {index}/{len(bigbeds)}: {relative}")
                return [str(bed_output), str(empty_marker)], []
            empty_marker.unlink(missing_ok=True)
            if not force and _combined_bigwig_is_reusable(
                output,
                inputs=inputs,
                chrom_sizes=chrom_sizes,
                method="bigbed",
            ):
                logger.log(
                    f"REUSE bigBed {index}/{len(bigbeds)}: {relative}; "
                    "validated completion marker"
                )
                reused = [str(output)]
                if bed_output.is_file():
                    reused.append(str(bed_output))
                return reused, []
            local_written: list[str] = []
            if not bed_output.exists():
                _combine_bigbed(inputs, bed_output, chrom_sizes)
                local_written.append(str(bed_output))
            executable = shutil.which("bedToBigBed")
            if executable:
                safe_name = str(relative).replace("/", "_").replace("\\", "_")
                sizes_path = destination_root / f"temporary_bigbed_{safe_name}.chrom.sizes"
                _write_chrom_sizes(sizes_path, chrom_sizes)
                partial_output = output.with_name(output.name + ".partial")
                partial_output.unlink(missing_ok=True)
                try:
                    _run_logged_subprocess(
                        [
                            executable,
                            str(bed_output),
                            str(sizes_path),
                            str(partial_output),
                        ],
                        logger,
                    )
                    _verify_bigwig(
                        partial_output, [chrom for chrom, _size in chrom_sizes]
                    )
                    os.replace(partial_output, output)
                    _write_combined_bigwig_marker(
                        output,
                        inputs=inputs,
                        chrom_sizes=chrom_sizes,
                        method="bigbed",
                    )
                finally:
                    sizes_path.unlink(missing_ok=True)
                local_written.append(str(output))
                logger.log(
                    f"PASS bigBed {index}/{len(bigbeds)}: {relative} "
                    f"({_format_bytes(output.stat().st_size)})"
                )
            else:
                logger.log(f"SKIP bigBed conversion: bedToBigBed not found for {relative}")
            return local_written, []
        except Exception as error:
            message = f"{relative}: {error}"
            logger.log(f"FAIL bigBed {index}/{len(bigbeds)}: {message}")
            return [], [message]

    if bigbeds:
        with ThreadPoolExecutor(max_workers=min(indexed_cores, len(bigbeds))) as executor:
            futures = [
                executor.submit(process_tree_bigbed, index, item)
                for index, item in enumerate(bigbeds, start=1)
            ]
            for future in as_completed(futures):
                local_written, local_warnings = future.result()
                written.extend(local_written)
                warnings.extend(local_warnings)

    # Dense BigWig tracks are deliberately combined after all other outputs.
    selected_bigwig_method = str(bigwig_method).lower()
    if selected_bigwig_method == "bedgraphs":
        selected_bigwig_method = "bedgraph"
    if selected_bigwig_method not in {"direct", "bedgraph"}:
        raise ValueError("BigWig combine method must be direct or bedgraph")
    temp_root = destination_root / "temporary_bedgraph_combine"
    logger.log(
        f"Beginning BigWig stage: {len(bigwigs):,} track(s). "
        f"Method: {selected_bigwig_method}. Chunk size: "
        f"{bigwig_chunk_size:,} bp. Temporary root: {temp_root}. "
        f"Workers: {min(indexed_cores, max(1, len(bigwigs)))}"
    )

    def staged_paths_for(inputs: Sequence[Path]) -> list[Path]:
        staged: list[Path] = []
        for input_path in inputs:
            matching_root = None
            for root in roots:
                try:
                    input_path.relative_to(root)
                except ValueError:
                    continue
                matching_root = root
                break
            if matching_root is None:
                raise ValueError(
                    f"Could not associate BigWig with a per-contig root: {input_path}"
                )
            relative_source = input_path.relative_to(matching_root).with_suffix(
                ".bedGraph"
            )
            staged.append(
                temp_root / "per_contig" / matching_root.name / relative_source
            )
        return staged

    def process_tree_bigwig(index: int, item: tuple[Path, list[Path]]) -> tuple[list[str], list[str]]:
        relative, inputs = item
        try:
            output = destination_root / relative
            logger.log(
                f"START BigWig {index}/{len(bigwigs)}: {relative}; "
                f"method={selected_bigwig_method}"
            )
            if not force and _combined_bigwig_is_reusable(
                output,
                inputs=inputs,
                chrom_sizes=chrom_sizes,
                method=selected_bigwig_method,
            ):
                logger.log(
                    f"REUSE BigWig {index}/{len(bigwigs)}: {relative}; "
                    "completion marker and input signatures match"
                )
                return [str(output)], []
            staged = (
                staged_paths_for(inputs)
                if selected_bigwig_method == "bedgraph"
                else None
            )
            _combine_bigwig(
                inputs,
                output,
                chrom_sizes,
                temp_root / "combined" / relative.parent / relative.stem,
                logger=logger,
                method=selected_bigwig_method,
                chunk_size=bigwig_chunk_size,
                staged_bedgraphs=staged,
            )
            logger.log(f"PASS BigWig {index}/{len(bigwigs)}: {relative}")
            return [str(output)], []
        except Exception as error:
            message = f"{relative}: {error}"
            logger.log(
                f"FAIL BigWig {index}/{len(bigwigs)}: {message}. "
                f"Incomplete or staged files retained under {temp_root}"
            )
            return [], [message]

    if combine_tracks and bigwigs:
        with ThreadPoolExecutor(max_workers=min(indexed_cores, len(bigwigs))) as executor:
            futures = [
                executor.submit(process_tree_bigwig, index, item)
                for index, item in enumerate(bigwigs, start=1)
            ]
            for future in as_completed(futures):
                local_written, local_warnings = future.result()
                written.extend(local_written)
                warnings.extend(local_warnings)
    elif bigwigs:
        warnings.append("BigWig combination was skipped by request")
        logger.log("SKIP BigWig stage: combination disabled")
    # Remove empty staging directories only after successful per-track cleanup.
    if temp_root.exists():
        for directory in sorted(
            (path for path in temp_root.rglob("*") if path.is_dir()),
            key=lambda path: len(path.parts),
            reverse=True,
        ):
            try:
                directory.rmdir()
            except OSError:
                pass
        try:
            temp_root.rmdir()
        except OSError:
            pass

    logger.log(
        f"Combination finished: {len(written):,} output(s), {len(warnings):,} warning(s). "
        f"Detailed log: {logger.path}"
    )
    logger.close()
    return {
        "output_dir": str(destination_root),
        "written": written,
        "warnings": warnings,
        "bigwig_stage_last": True,
        "combine_log": str(destination_root / combine_log_name),
        "temporary_bedgraph_root": str(temp_root),
        "bigwig_method": selected_bigwig_method,
        "cores": cores,
        "streaming_combine_cores": streaming_cores,
        "indexed_combine_cores": indexed_cores,
        "bigwig_combine_chunk_bp": bigwig_chunk_size,
        "force": bool(force),
        "incomplete_groups": [str(path) for path in sorted(incomplete_groups)],
        "unsupported_groups": [str(path) for path in sorted(unsupported_groups)],
        "combined_chromosomes": [chrom for chrom, _size in chrom_sizes],
    }
