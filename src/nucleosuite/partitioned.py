"""Generic per-contig orchestration for BED-driven analysis commands."""

from __future__ import annotations

import argparse
import importlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Sequence

from nucleosuite.io import open_text as open_interval_text
from nucleosuite.parallel import MANIFEST_NAME, _safe_contig, _worker_initializer
from nucleosuite.core.regions import resolve_contig_name


def _open(path: str | Path, mode: str):
    if mode != "rt":
        raise ValueError("Partitioned interval inputs are opened for reading only")
    return open_interval_text(path)


def _is_indexed_interval(path: str | Path) -> bool:
    value = Path(path)
    lower = value.name.lower()
    if lower.endswith((".bb", ".bigbed")):
        return True
    return lower.endswith((".bed.gz", ".bed.bgz", ".tsv.gz", ".tsv.bgz")) and (
        Path(str(value) + ".tbi").exists() or Path(str(value) + ".csi").exists()
    )


def _indexed_contigs(path: str | Path) -> list[str]:
    value = Path(path)
    lower = value.name.lower()
    if lower.endswith((".bb", ".bigbed")):
        try:
            import pyBigWig
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("pyBigWig is required for bigBed parallel input") from exc
        handle = pyBigWig.open(str(value))
        try:
            return [str(name) for name in handle.chroms()]
        finally:
            handle.close()
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pysam is required for tabix parallel input") from exc
    handle = pysam.TabixFile(str(value))
    try:
        return [str(name) for name in handle.contigs]
    finally:
        handle.close()


def _contigs(path: str | Path) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()
    with _open(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            chrom = line.split()[0]
            if chrom not in seen:
                seen.add(chrom)
                output.append(chrom)
    return output


def _filter(path: str | Path, output: Path, contig: str) -> None:
    """Materialise one contig from an indexed interval source."""
    output.parent.mkdir(parents=True, exist_ok=True)
    value = Path(path)
    lower = value.name.lower()
    if lower.endswith((".bb", ".bigbed")):
        try:
            import pyBigWig
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("pyBigWig is required for bigBed parallel input") from exc
        handle = pyBigWig.open(str(value))
        try:
            source_contig = resolve_contig_name(
                contig, list(handle.chroms()), source_label="bigBed index"
            )
            length = int(handle.chroms(source_contig))
            entries = handle.entries(source_contig, 0, length) or []
            with output.open("w", encoding="utf-8") as destination:
                for start, end, rest in entries:
                    suffix = f"\t{rest}" if rest else ""
                    destination.write(f"{contig}\t{int(start)}\t{int(end)}{suffix}\n")
        finally:
            handle.close()
        return
    if _is_indexed_interval(value):
        try:
            import pysam
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("pysam is required for tabix parallel input") from exc
        handle = pysam.TabixFile(str(value))
        try:
            source_contig = resolve_contig_name(
                contig, list(handle.contigs), source_label="tabix index"
            )
            with output.open("w", encoding="utf-8") as destination:
                for raw in handle.fetch(source_contig):
                    fields = raw.rstrip("\r\n").split("\t")
                    if fields:
                        fields[0] = contig
                    destination.write("\t".join(fields) + "\n")
        finally:
            handle.close()
        return
    # This branch is used only by direct tests/helpers.  Parallel orchestration
    # rejects unindexed inputs before creating workers.
    with _open(value, "rt") as source, output.open("w", encoding="utf-8") as destination:
        for raw in source:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(("#", "track", "browser")):
                destination.write(raw if raw.endswith("\n") else raw + "\n")
                continue
            source_contig = line.split()[0]
            try:
                matched = resolve_contig_name(
                    source_contig, [contig], source_label="partition contig"
                )
            except KeyError:
                continue
            if matched == contig:
                destination.write(raw if raw.endswith("\n") else raw + "\n")


def _split_named(value: str) -> tuple[str | None, str]:
    if "=" in value:
        name, path = value.split("=", 1)
        return name, path
    return None, value


def _worker(module_name: str, function_name: str, namespace_data: dict) -> int:
    module = importlib.import_module(module_name)
    function = getattr(module, function_name)
    return int(function(argparse.Namespace(**namespace_data)) or 0)


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def run_partitioned_command(
    command: str,
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
    *,
    runner_module: str,
    runner_function: str,
    primary_attr: str | None,
    primary_list_attr: str | None = None,
    output_prefix_attr: str,
    path_attrs: Sequence[str] = (),
    named_path_list_attrs: Sequence[str] = (),
    clear_output_attrs: Sequence[str] = (),
    output_dir_attr: str | None = None,
    base_name: str | None = None,
) -> int:
    """Run a BED-driven command per contig, then rerun the exact combined analysis.

    The per-contig results are retained. The combined analysis uses the original
    unsplit inputs, which preserves global thresholds, correlations and summary
    statistics for commands that do not expose sufficient-statistic files.
    """
    cores = int(getattr(args, "cores", 1) or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    if getattr(args, "_per_contig_worker", False) or cores == 1:
        return int(serial_runner(args) or 0)
    if primary_list_attr:
        values = getattr(args, primary_list_attr, None) or []
        if not values:
            return int(serial_runner(args) or 0)
        _name, primary = _split_named(str(values[0]))
    else:
        if primary_attr is None:
            raise ValueError("A primary interval input is required for per-contig processing")
        primary = getattr(args, primary_attr)
        if not primary:
            return int(serial_runner(args) or 0)
    interval_paths = [str(primary)]
    for attr in path_attrs:
        value = getattr(args, attr, None)
        if value:
            interval_paths.append(str(value))
    for attr in named_path_list_attrs:
        for value in getattr(args, attr, None) or []:
            _name, path = _split_named(str(value))
            interval_paths.append(path)
    unindexed = [path for path in interval_paths if not _is_indexed_interval(path)]
    if unindexed:
        print(
            "Parallel contig processing was not enabled because these interval "
            "inputs lack bigBed or bgzip/tabix random-access indexes: "
            + ", ".join(unindexed)
            + ". Running serially."
        )
        return int(serial_runner(args) or 0)
    contigs = _indexed_contigs(primary)
    if len(contigs) <= 1:
        return int(serial_runner(args) or 0)

    requested = getattr(args, output_prefix_attr, None)
    base = base_name or (Path(str(requested)).name if requested else command.replace("-", "_"))
    requested_parent = Path(str(requested)).parent if requested else Path(".")
    if output_dir_attr and getattr(args, output_dir_attr, None):
        requested_parent = Path(str(getattr(args, output_dir_attr)))
    root = (
        Path(args.parallel_dir)
        if getattr(args, "parallel_dir", None)
        else requested_parent / f"{base}_multicontig"
    ).resolve()
    per_root = root / "per_contig"
    combined_root = root / "combined"
    per_root.mkdir(parents=True, exist_ok=True)
    combined_root.mkdir(parents=True, exist_ok=True)

    original_data = {
        key: _jsonable(value)
        for key, value in vars(args).items()
        if key not in {"command", "command_function", "command_runner", "func"}
    }
    entries: list[dict[str, object]] = []
    failures: list[str] = []
    with ProcessPoolExecutor(
        max_workers=min(cores, len(contigs)), initializer=_worker_initializer
    ) as executor:
        futures = {}
        for contig in contigs:
            safe = _safe_contig(contig)
            directory = per_root / safe
            directory.mkdir(parents=True, exist_ok=True)
            data = dict(original_data)
            if primary_attr is not None:
                filtered_primary = directory / f"{primary_attr}.bed"
                _filter(primary, filtered_primary, contig)
                data[primary_attr] = str(filtered_primary)
            for attr in path_attrs:
                value = getattr(args, attr, None)
                if value:
                    target = directory / f"{attr}.bed"
                    _filter(value, target, contig)
                    data[attr] = str(target)
            for attr in named_path_list_attrs:
                values = getattr(args, attr, None) or []
                filtered_values: list[str] = []
                for index, value in enumerate(values):
                    name, path = _split_named(str(value))
                    target = directory / f"{attr}_{index}.bed"
                    _filter(path, target, contig)
                    filtered_values.append(f"{name}={target}" if name is not None else str(target))
                data[attr] = filtered_values
            data[output_prefix_attr] = str(directory / f"{base}_{safe}")
            if output_dir_attr:
                data[output_dir_attr] = str(directory)
                # Prefixes are usually basenames when a separate output directory is used.
                data[output_prefix_attr] = f"{base}_{safe}"
            for attr in clear_output_attrs:
                data[attr] = None
            data.update(
                {
                    "cores": 1,
                    "parallel_dir": None,
                    "skip_combine": True,
                    "skip_combined_tracks": True,
                    "_per_contig_worker": True,
                }
            )
            future = executor.submit(_worker, runner_module, runner_function, data)
            futures[future] = (contig, data[output_prefix_attr], directory)
        for future in as_completed(futures):
            contig, prefix, directory = futures[future]
            try:
                code = future.result()
                if code:
                    failures.append(f"{contig}: exit code {code}")
                entries.append({"contig": contig, "prefix": str(directory / str(prefix)) if output_dir_attr else str(prefix), "exit_code": code})
                if code == 0:
                    print(f"Completed {contig}")
            except Exception as error:
                failures.append(f"{contig}: {error}")
                entries.append({"contig": contig, "prefix": str(prefix), "exit_code": 2, "error": str(error)})

    order = {contig: index for index, contig in enumerate(contigs)}
    entries.sort(key=lambda item: order[str(item["contig"])])
    combined_data = dict(original_data)
    if output_dir_attr:
        combined_data[output_dir_attr] = str(combined_root)
        combined_data[output_prefix_attr] = base
    else:
        combined_data[output_prefix_attr] = str(combined_root / base)
    for attr in clear_output_attrs:
        combined_data[attr] = None
    combined_data.update(
        {
            "cores": 1,
            "parallel_dir": None,
            "skip_combine": True,
            "skip_combined_tracks": True,
            "_per_contig_worker": True,
        }
    )
    manifest = {
        "schema_version": 1,
        "command": command,
        "combine_strategy": "rerun_namespace",
        "base_name": base,
        "combined_name": base,
        "root_dir": str(root),
        "per_contig_dir": str(per_root),
        "combined_dir": str(combined_root),
        "per_contig": entries,
        "rerun_module": runner_module,
        "rerun_function": runner_function,
        "rerun_namespace": combined_data,
        "options": {"cores": cores},
    }
    with (root / MANIFEST_NAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, default=str)
        handle.write("\n")
    if failures:
        raise RuntimeError("Per-contig jobs failed: " + "; ".join(failures))
    if not getattr(args, "skip_combine", False):
        result = rerun_manifest(root, manifest)
        print(f"Combined outputs: {result['output_dir']}")
    else:
        print(f"Per-contig outputs: {per_root}")
        print(f"Combine later with: nucleosuite combine --input-dir {root}")
    return 0


def rerun_manifest(root: Path, manifest: dict, *, output_dir: Path | None = None) -> dict[str, object]:
    data = dict(manifest["rerun_namespace"])
    if output_dir is not None:
        old_root = Path(str(manifest.get("combined_dir", root / "combined")))
        for key, value in list(data.items()):
            if isinstance(value, str) and value.startswith(str(old_root)):
                data[key] = str(output_dir / Path(value).relative_to(old_root))
        for key in ("output_dir",):
            if key in data and data[key] is not None:
                data[key] = str(output_dir)
    target = Path(output_dir or manifest.get("combined_dir", root / "combined"))
    target.mkdir(parents=True, exist_ok=True)
    before = {path for path in target.rglob("*") if path.is_file()}
    code = _worker(str(manifest["rerun_module"]), str(manifest["rerun_function"]), data)
    if code:
        raise RuntimeError(f"Combined {manifest.get('command')} analysis exited with code {code}")
    after = {path for path in target.rglob("*") if path.is_file()}
    return {
        "input_dir": str(root),
        "output_dir": str(target),
        "combined_prefix": str(target / str(manifest.get("combined_name", "combined"))),
        "written": [str(path) for path in sorted(after - before)],
        "warnings": [],
        "bigwig_stage_last": True,
    }
