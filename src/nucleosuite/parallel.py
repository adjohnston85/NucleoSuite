"""Per-contig process orchestration for fragment-derived commands."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Mapping, Sequence

from nucleosuite.core.chrom_sizes import read_chrom_sizes_source

MANIFEST_NAME = "nucleosuite_multicontig_manifest.json"

_THREAD_LIMIT_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}

_CHECKPOINT_SCHEMA_VERSION = 1
_RUNTIME_ARGUMENTS = {
    "cores",
    "parallel_dir",
    "skip_combine",
    "skip_combined_tracks",
    "streaming_combine_cores",
    "indexed_combine_cores",
    "combine_chunk_bp",
    "force",
    "_per_contig_worker",
}


def _checkpoint_value(value):
    """Return a deterministic, input-aware representation for restart checks."""
    if isinstance(value, Path):
        value = str(value)
    if isinstance(value, (list, tuple)):
        return [_checkpoint_value(item) for item in value]
    if isinstance(value, dict):
        return {
            str(key): _checkpoint_value(item)
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, str):
        path = Path(value)
        if path.is_file():
            stat = path.stat()
            record: dict[str, object] = {
                "path": str(path.resolve()),
                "size_bytes": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
            }
            if stat.st_size <= 10_000_000:
                record["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            return record
        return value
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _worker_signature(
    command: str,
    namespace_data: Mapping[str, object],
    *,
    contig_specs: Sequence[str],
) -> str:
    parameters = {
        key: _checkpoint_value(value)
        for key, value in sorted(namespace_data.items())
        if key not in _RUNTIME_ARGUMENTS
    }
    payload = {
        "schema_version": _CHECKPOINT_SCHEMA_VERSION,
        "command": command,
        "contig_specs": list(contig_specs),
        "parameters": parameters,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _checkpoint_path(directory: Path, command: str) -> Path:
    return directory / f".{command.replace('-', '_')}.complete.json"


def _output_signature(directory: Path, checkpoint: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(directory.rglob("*")):
        if not path.is_file() or path == checkpoint:
            continue
        if path.name.endswith((".partial", ".status")) or path.name == "track_spec.tsv":
            continue
        stat = path.stat()
        rows.append(
            {
                "path": str(path.relative_to(directory)),
                "size_bytes": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
            }
        )
    return rows


def _checkpoint_reusable(checkpoint: Path, signature: str) -> bool:
    if not checkpoint.is_file():
        return False
    try:
        payload = json.loads(checkpoint.read_text(encoding="utf-8"))
        if (
            payload.get("schema_version") != _CHECKPOINT_SCHEMA_VERSION
            or not payload.get("complete")
            or payload.get("signature") != signature
        ):
            return False
        root = checkpoint.parent
        for item in payload.get("outputs", []):
            path = root / str(item["path"])
            if not path.is_file():
                return False
            stat = path.stat()
            if (
                stat.st_size != int(item["size_bytes"])
                or stat.st_mtime_ns != int(item.get("mtime_ns", -1))
            ):
                return False
        return bool(payload.get("outputs"))
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return False


def _write_checkpoint(
    checkpoint: Path,
    *,
    command: str,
    contig: str,
    signature: str,
) -> None:
    outputs = _output_signature(checkpoint.parent, checkpoint)
    if not outputs:
        raise RuntimeError(f"No completed outputs were found for {command} {contig}")
    payload = {
        "schema_version": _CHECKPOINT_SCHEMA_VERSION,
        "complete": True,
        "command": command,
        "contig": contig,
        "signature": signature,
        "completed_unix_time": time.time(),
        "outputs": outputs,
    }
    temporary = checkpoint.with_name(checkpoint.name + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, checkpoint)


def _prefix_checkpoint_path(prefix: Path, command: str) -> Path:
    return prefix.with_name(f"{prefix.name}.{command.replace('-', '_')}.complete.json")


def _write_prefix_checkpoint(
    checkpoint: Path,
    *,
    prefix: Path,
    command: str,
    signature: str,
) -> None:
    outputs: list[dict[str, object]] = []
    for path in sorted(prefix.parent.glob(prefix.name + "*")):
        if not path.is_file() or path == checkpoint or ".partial" in path.name:
            continue
        stat = path.stat()
        outputs.append(
            {
                "path": path.name,
                "size_bytes": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
            }
        )
    if not outputs:
        raise RuntimeError(f"No completed outputs were found for {command}")
    payload = {
        "schema_version": _CHECKPOINT_SCHEMA_VERSION,
        "complete": True,
        "command": command,
        "contig": "serial",
        "signature": signature,
        "completed_unix_time": time.time(),
        "outputs": outputs,
    }
    temporary = checkpoint.with_name(checkpoint.name + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, checkpoint)


def _run_resumable_serial(
    command: str,
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
) -> int:
    namespace_data = {
        key: value
        for key, value in vars(args).items()
        if key not in {"command_runner", "command_function"}
    }
    requested_prefix = _base_output_prefix(args)
    actual_prefix = Path(
        _resolved_output_prefix(command, requested_prefix, namespace_data)
    ).resolve()
    checkpoint = _prefix_checkpoint_path(actual_prefix, command)
    contig_specs = list(getattr(args, "contigs", None) or ["all"])
    signature = _worker_signature(
        command, namespace_data, contig_specs=contig_specs
    )
    if (
        not bool(getattr(args, "force", False))
        and _checkpoint_reusable(checkpoint, signature)
    ):
        print(
            f"Reused completed {command} outputs for {', '.join(contig_specs)}",
            flush=True,
        )
        return 0
    checkpoint.unlink(missing_ok=True)
    print(
        f"Starting {command}: contigs {', '.join(contig_specs)}; "
        f"output prefix {actual_prefix}",
        flush=True,
    )
    exit_code = int(serial_runner(args) or 0)
    if exit_code == 0:
        _write_prefix_checkpoint(
            checkpoint,
            prefix=actual_prefix,
            command=command,
            signature=signature,
        )
        print(f"Completed {command}: {actual_prefix}", flush=True)
    return exit_code


def _safe_contig(value: str) -> str:
    """Return a filesystem-safe contig label."""
    return value.replace("/", "_").replace("\\", "_").replace(":", "_")


def _worker_initializer() -> None:
    """Prevent nested BLAS/OpenMP pools inside one-contig worker processes."""
    for name, value in _THREAD_LIMIT_ENV.items():
        os.environ[name] = value


def _has_tabix_index(path: Path) -> bool:
    return Path(str(path) + ".tbi").exists() or Path(str(path) + ".csi").exists()


def _fragment_inputs_support_parallel(args: argparse.Namespace) -> tuple[bool, str]:
    """Return whether fragment inputs permit independent contig access.

    BAM inputs are indexed on demand. Fragment interval inputs are parallelised
    only when they are bigBed files or bgzip/tabix-indexed text files. Plain BED
    and ordinary gzip files deliberately fall back to one serial pass instead
    of being rescanned once per worker.
    """
    bam_paths = [Path(value) for value in (getattr(args, "bamfiles", None) or [])]
    if bam_paths:
        try:
            import pysam
        except ImportError as exc:  # pragma: no cover
            return False, f"pysam is unavailable: {exc}"
        for path in bam_paths:
            if not path.exists():
                return False, f"input BAM does not exist: {path}"
            index_candidates = (
                Path(str(path) + ".bai"),
                path.with_suffix(".bai"),
                Path(str(path) + ".csi"),
            )
            if not any(candidate.exists() for candidate in index_candidates):
                try:
                    pysam.index(str(path))
                except (OSError, ValueError, pysam.SamtoolsError) as exc:
                    return False, f"BAM could not be indexed: {path}: {exc}"
        return True, "indexed BAM input"

    fragment_paths = [Path(value) for value in (getattr(args, "fragment_files", None) or [])]
    if not fragment_paths:
        return False, "no fragment input was supplied"
    for path in fragment_paths:
        lower = path.name.lower()
        if lower.endswith((".bb", ".bigbed")):
            continue
        if lower.endswith((".bed.gz", ".bed.bgz", ".tsv.gz", ".tsv.bgz")) and _has_tabix_index(path):
            continue
        return False, (
            f"{path} is not an indexed BAM, bigBed, or bgzip/tabix-indexed interval file"
        )
    return True, "indexed fragment interval input"

_NATIVE_MODULES = {
    "pns": "nucleosuite.cli.pns",
    "wps": "nucleosuite.cli.wps",
    "coverage": "nucleosuite.cli.coverage",
    "dyads": "nucleosuite.cli.dyads",
    "dyad": "nucleosuite.cli.dyads",
    "fragment-ends": "nucleosuite.cli.fragment_ends",
    "dinuc-profile": "nucleosuite.cli.dinuc_profile",
    "ww-types": "nucleosuite.cli.ww_types",
    "fragments": "nucleosuite.fragments_command",
    "fragment-lengths": "nucleosuite.fragment_lengths",
    "randomize-fragments": "nucleosuite.randomize_fragments_command",
}


def add_parallel_arguments(
    parser: argparse.ArgumentParser,
    *,
    combine_resources: bool = False,
    resumable: bool = False,
    cores_option: str = "--cores",
    cores_help: str | None = None,
) -> None:
    core_options = (
        (cores_option, "--cores") if cores_option != "--cores" else ("--cores",)
    )
    parser.add_argument(
        *core_options,
        dest="cores",
        type=int,
        default=1,
        help=cores_help
        or (
            "Number of contigs processed concurrently. Values greater than 1 create "
            "per-contig directories and a combined output directory (default: 1)."
        ),
    )
    parser.add_argument(
        "--parallel-dir",
        default=None,
        help=(
            "Root directory for per-contig and combined outputs. Default: "
            "<output-prefix>_multicontig."
        ),
    )
    parser.add_argument(
        "--skip-combine",
        action="store_true",
        help="Retain per-contig outputs without running the final combine stage.",
    )
    parser.add_argument(
        "--skip-combined-tracks",
        action="store_true",
        help="Combine tabular and interval outputs but defer combined BigWig creation.",
    )
    if combine_resources:
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
                "Concurrent BigWig/BigBed combines. This finalization-sensitive "
                "budget is independent of --cores (default: 1)."
            ),
        )
        parser.add_argument(
            "--combine-chunk-bp",
            type=int,
            default=100_000,
            help="Genomic query chunk used while combining BigWigs (default: 100000 bp).",
        )
    if resumable:
        parser.add_argument(
            "--force",
            action="store_true",
            help="Recompute compatible completed contigs and combined tracks.",
        )
    parser.add_argument("--_per-contig-worker", action="store_true", help=argparse.SUPPRESS)


def _combine_resource_values(args: argparse.Namespace, cores: int) -> tuple[int, int, int]:
    streaming = int(getattr(args, "streaming_combine_cores", None) or cores)
    indexed = int(getattr(args, "indexed_combine_cores", 1) or 1)
    chunk = int(getattr(args, "combine_chunk_bp", 100_000) or 100_000)
    if streaming < 1:
        raise ValueError("--streaming-combine-cores must be at least 1")
    if indexed < 1:
        raise ValueError("--indexed-combine-cores must be at least 1")
    if chunk < 1:
        raise ValueError("--combine-chunk-bp must be at least 1")
    return streaming, indexed, chunk


def _stable_seed(seed: int | None, contig: str) -> int | None:
    if seed is None:
        return None
    digest = hashlib.blake2b(contig.encode("utf-8"), digest_size=4).digest()
    return (int(seed) + int.from_bytes(digest, "big")) % (2**32)


def _resolved_output_prefix(command: str, base: str, args: Mapping[str, object]) -> str:
    lower = int(args.get("frag_lower", 0))
    upper = int(args.get("frag_upper", 0))
    if command == "pns":
        return f"{base}_mode{int(args['mode_length'])}_lower{lower}_upper{upper}"
    if command == "wps":
        value = (
            f"{base}_prot{int(args['protection'])}_lower{lower}_upper{upper}"
        )
        if args.get("randomize_mode") != "none":
            value += f"_rand{args.get('randomize_mode')}"
        return value
    if command == "coverage":
        return f"{base}_coverage_lower{lower}_upper{upper}"
    if command in {"dyads", "dyad"}:
        return f"{base}_dyads_lower{lower}_upper{upper}"
    if command == "fragment-ends":
        return f"{base}_fragment_ends_lower{lower}_upper{upper}"
    if command == "dinuc-profile":
        return f"{base}_dinuc_lower{lower}_upper{upper}"
    if command == "ww-types":
        return f"{base}_wwtypes_lower{lower}_upper{upper}"
    return base


def _worker(command: str, namespace_data: dict, contig_specs: list[str], requested_prefix: str, seed: int | None) -> tuple[int, str]:
    import importlib

    module = importlib.import_module(_NATIVE_MODULES[command])
    data = dict(namespace_data)
    data["contigs"] = [spec.split(":", 1)[0] for spec in contig_specs] if command == "fragment-lengths" else contig_specs
    if command in {"fragments", "randomize-fragments"}:
        data["output_prefix"] = requested_prefix
    elif command == "fragment-lengths":
        data["output"] = requested_prefix + ".tsv"
        data["output_dir"] = str(Path(requested_prefix).parent)
        if data.get("plot"):
            plot_suffix = Path(str(data["plot"])).suffix or ".png"
            data["plot"] = requested_prefix + plot_suffix
    else:
        data["out_prefix"] = requested_prefix
    data["cores"] = 1
    data["parallel_dir"] = None
    data["skip_combine"] = True
    data["skip_combined_tracks"] = True
    data["_per_contig_worker"] = True
    data["seed"] = seed
    namespace = argparse.Namespace(**data)
    exit_code = int(module.run(namespace) or 0)
    actual_prefix = _resolved_output_prefix(command, requested_prefix, data)
    return exit_code, actual_prefix


def _resolve_contig_jobs(args: argparse.Namespace) -> tuple[list[tuple[str, list[str], int]], list[tuple[str, int]]]:
    fasta = None
    source = None
    from nucleosuite.core.fragment_inputs import open_fragment_source
    from nucleosuite.core.regions import expand_contig_tokens

    try:
        if getattr(args, "fasta", None):
            import pysam
            fasta = pysam.FastaFile(args.fasta)
        source = open_fragment_source(
            bam_paths=getattr(args, "bamfiles", None),
            fragment_paths=getattr(args, "fragment_files", None),
            chrom_sizes=getattr(args, "chrom_sizes", None),
            fasta=fasta,
        )
        references = list(source.references)
        lengths = [int(value) for value in source.lengths]
        raw_contigs = getattr(args, "contigs", None)
        if isinstance(raw_contigs, str):
            raw_contigs = [raw_contigs]
        specs = expand_contig_tokens(raw_contigs, references)
        specs_by_contig: dict[str, list[str]] = {}
        for spec in specs:
            contig = spec.split(":", 1)[0]
            specs_by_contig.setdefault(contig, []).append(spec)
        order = {contig: index for index, contig in enumerate(references)}
        length_by_contig = dict(zip(references, lengths))
        jobs = [
            (contig, specs_by_contig[contig], length_by_contig[contig])
            for contig in sorted(specs_by_contig, key=lambda value: order[value])
        ]
        chrom_sizes = [(contig, length_by_contig[contig]) for contig, _, _ in jobs]
        return jobs, chrom_sizes
    finally:
        if source is not None:
            source.close()
        if fasta is not None:
            fasta.close()


def _base_output_prefix(args: argparse.Namespace) -> str:
    from nucleosuite.workflows.common import default_output_prefix, input_paths_from_args
    if getattr(args, "out_prefix", None):
        return str(args.out_prefix)
    if getattr(args, "output_prefix", None):
        return str(args.output_prefix)
    if getattr(args, "output", None):
        output = str(args.output)
        return output[:-4] if output.lower().endswith(".tsv") else output
    return default_output_prefix(input_paths_from_args(args), getattr(args, "contigs", None))


def run_native_per_contig(
    command: str,
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
) -> int:
    """Run a native fragment command serially or once per contig."""
    cores = int(getattr(args, "cores", 1) or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    streaming_cores, indexed_cores, combine_chunk_bp = _combine_resource_values(
        args, cores
    )
    if getattr(args, "_per_contig_worker", False):
        return int(serial_runner(args) or 0)
    if cores == 1:
        return _run_resumable_serial(command, args, serial_runner)

    supported, reason = _fragment_inputs_support_parallel(args)
    if not supported:
        print(
            "Parallel contig processing was not enabled because the input does not "
            f"support efficient independent contig access ({reason}). Running serially."
        )
        return _run_resumable_serial(command, args, serial_runner)

    jobs, chrom_sizes = _resolve_contig_jobs(args)
    if len(jobs) <= 1:
        return _run_resumable_serial(command, args, serial_runner)

    requested_base = _base_output_prefix(args)
    requested_path = Path(requested_base)
    root = (
        Path(args.parallel_dir)
        if getattr(args, "parallel_dir", None)
        else requested_path.parent / f"{requested_path.name}_multicontig"
    ).resolve()
    per_contig_root = root / "per_contig"
    combined_dir = root / "combined"
    per_contig_root.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    base_name = requested_path.name

    namespace_data = {
        key: value
        for key, value in vars(args).items()
        if key not in {"command_runner", "command_function"}
    }
    manifest_entries: list[dict[str, object]] = []
    failures: list[str] = []

    max_workers = min(cores, len(jobs))
    print(f"Processing {len(jobs)} contigs with {max_workers} worker process(es).")
    with ProcessPoolExecutor(
        max_workers=max_workers, initializer=_worker_initializer
    ) as executor:
        future_map = {}
        for contig, specs, _length in jobs:
            safe = contig.replace("/", "_").replace("\\", "_").replace(":", "_")
            contig_dir = per_contig_root / safe
            contig_dir.mkdir(parents=True, exist_ok=True)
            requested_prefix = str(contig_dir / f"{base_name}_{safe}")
            checkpoint = _checkpoint_path(contig_dir, command)
            signature = _worker_signature(
                command, namespace_data, contig_specs=specs
            )
            if (
                not bool(getattr(args, "force", False))
                and _checkpoint_reusable(checkpoint, signature)
            ):
                actual_prefix = _resolved_output_prefix(
                    command, requested_prefix, namespace_data
                )
                print(f"Reused completed {contig}", flush=True)
                manifest_entries.append(
                    {
                        "contig": contig,
                        "contig_specs": specs,
                        "requested_prefix": requested_prefix,
                        "prefix": actual_prefix,
                        "exit_code": 0,
                        "status": "reused",
                        "checkpoint": str(checkpoint),
                    }
                )
                continue
            checkpoint.unlink(missing_ok=True)
            print(
                f"Starting {command}: {contig}; output prefix {requested_prefix}",
                flush=True,
            )
            future = executor.submit(
                _worker,
                command,
                namespace_data,
                specs,
                requested_prefix,
                _stable_seed(getattr(args, "seed", None), contig),
            )
            future_map[future] = (
                contig,
                specs,
                requested_prefix,
                checkpoint,
                signature,
            )

        for future in as_completed(future_map):
            contig, specs, requested_prefix, checkpoint, signature = future_map[future]
            try:
                exit_code, actual_prefix = future.result()
                if exit_code != 0:
                    failures.append(f"{contig}: exit code {exit_code}")
                else:
                    _write_checkpoint(
                        checkpoint,
                        command=command,
                        contig=contig,
                        signature=signature,
                    )
                    print(f"Completed {contig}", flush=True)
                manifest_entries.append(
                    {
                        "contig": contig,
                        "contig_specs": specs,
                        "requested_prefix": requested_prefix,
                        "prefix": actual_prefix,
                        "exit_code": exit_code,
                        "status": "completed" if exit_code == 0 else "failed",
                        "checkpoint": str(checkpoint),
                    }
                )
            except Exception as error:
                failures.append(f"{contig}: {error}")
                manifest_entries.append(
                    {
                        "contig": contig,
                        "contig_specs": specs,
                        "requested_prefix": requested_prefix,
                        "prefix": _resolved_output_prefix(command, requested_prefix, namespace_data),
                        "exit_code": 2,
                        "error": str(error),
                    }
                )

    order = {contig: index for index, (contig, _specs, _length) in enumerate(jobs)}
    manifest_entries.sort(key=lambda item: order[str(item["contig"])])
    combined_requested = str(combined_dir / base_name)
    combined_name = Path(_resolved_output_prefix(command, combined_requested, namespace_data)).name
    manifest = {
        "schema_version": 1,
        "command": command,
        "base_name": base_name,
        "root_dir": str(root),
        "per_contig_dir": str(per_contig_root),
        "combined_dir": str(combined_dir),
        "combined_name": combined_name,
        "chrom_sizes": [{"chrom": chrom, "size": size} for chrom, size in chrom_sizes],
        "per_contig": manifest_entries,
        "options": {
            "cores": cores,
            "seed": getattr(args, "seed", None),
            "skip_combine": bool(getattr(args, "skip_combine", False)),
            "skip_combined_tracks": bool(getattr(args, "skip_combined_tracks", False)),
        },
    }
    with (root / MANIFEST_NAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")

    if failures:
        raise RuntimeError("Per-contig jobs failed: " + "; ".join(failures))

    if not getattr(args, "skip_combine", False):
        from nucleosuite.combine import combine_run
        result = combine_run(
            root,
            combine_tracks=not getattr(args, "skip_combined_tracks", False),
            cores=cores,
            streaming_cores=streaming_cores,
            indexed_cores=indexed_cores,
            bigwig_chunk_size=combine_chunk_bp,
            force=bool(getattr(args, "force", False)),
        )
        print(f"Combined outputs: {result['output_dir']}")
        for warning in result["warnings"]:
            print(f"WARNING: {warning}")
    else:
        print(f"Per-contig outputs: {per_contig_root}")
        print(f"Combine later with: nucleosuite combine --input-dir {root}")
    return 0


def _tracks_worker(
    namespace_data: dict,
    contig_specs: list[str],
    output_dir: str,
    spec_file: str,
    report_path: str | None,
) -> int:
    from nucleosuite.cli.tracks import run as run_tracks

    data = dict(namespace_data)
    data.update(
        {
            "contigs": contig_specs,
            "fragment_range": [],
            "spec_file": [spec_file],
            "output_dir": output_dir,
            "report": report_path,
            "cores": 1,
            "parallel_dir": None,
            "skip_combine": True,
            "skip_combined_tracks": True,
            "_per_contig_worker": True,
        }
    )
    return int(run_tracks(argparse.Namespace(**data)) or 0)


def _tracks_relative_specs(
    args: argparse.Namespace,
) -> tuple[list[object], Path, str]:
    """Resolve track specifications and ensure outputs share one tree root."""
    from nucleosuite.workflows.tracks import load_specs

    copied = copy.deepcopy(args)
    specs = load_specs(copied)
    output_root = Path(str(args.output_dir)).resolve()
    base_name = str(copied.output_prefix or Path(output_root).name or "tracks")
    for spec in specs:
        path = Path(str(spec.output_prefix))
        absolute = path.resolve() if path.is_absolute() else (Path.cwd() / path).resolve()
        try:
            absolute.relative_to(output_root)
        except ValueError as exc:
            raise ValueError(
                "Parallel tracks requires every specification output_prefix to be "
                f"inside --output-dir ({output_root}); observed {absolute}"
            ) from exc
    return specs, output_root, base_name


def _run_resumable_tracks_serial(
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
) -> int:
    _specs, output_root, _base_name = _tracks_relative_specs(args)
    output_root.mkdir(parents=True, exist_ok=True)
    namespace_data = {
        key: value
        for key, value in vars(args).items()
        if key not in {"command_runner", "command_function"}
    }
    contig_specs = list(getattr(args, "contigs", None) or ["all"])
    signature = _worker_signature(
        "tracks", namespace_data, contig_specs=contig_specs
    )
    checkpoint = _checkpoint_path(output_root, "tracks")
    if (
        not bool(getattr(args, "force", False))
        and _checkpoint_reusable(checkpoint, signature)
    ):
        print(
            f"Reused completed tracks outputs for {', '.join(contig_specs)}",
            flush=True,
        )
        return 0
    checkpoint.unlink(missing_ok=True)
    print(
        f"Starting tracks: contigs {', '.join(contig_specs)}; "
        f"output directory {output_root}",
        flush=True,
    )
    exit_code = int(serial_runner(args) or 0)
    if exit_code == 0:
        _write_checkpoint(
            checkpoint,
            command="tracks",
            contig="serial",
            signature=signature,
        )
        print(f"Completed tracks: {output_root}", flush=True)
    return exit_code


def _write_tracks_worker_spec(
    specs: Sequence[object],
    original_root: Path,
    worker_root: Path,
    destination: Path,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("fragment_range", "output_prefix", "tracks", "basic_scope"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for spec in specs:
            source = Path(str(spec.output_prefix))
            absolute = source.resolve() if source.is_absolute() else (Path.cwd() / source).resolve()
            relative = absolute.relative_to(original_root)
            writer.writerow(
                {
                    "fragment_range": spec.fragment_range.label,
                    "output_prefix": str(worker_root / relative),
                    "tracks": ",".join(spec.tracks),
                    "basic_scope": spec.basic_scope,
                }
            )


def run_tracks_per_contig(
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
) -> int:
    """Run the shared multi-track engine once per indexed contig and combine."""
    cores = int(getattr(args, "cores", 1) or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    streaming_cores, indexed_cores, combine_chunk_bp = _combine_resource_values(
        args, cores
    )
    if getattr(args, "_per_contig_worker", False):
        return int(serial_runner(args) or 0)
    if cores == 1:
        return _run_resumable_tracks_serial(args, serial_runner)
    supported, reason = _fragment_inputs_support_parallel(args)
    if not supported:
        print(
            "Parallel contig processing was not enabled because the input does not "
            f"support efficient independent contig access ({reason}). Running serially."
        )
        return _run_resumable_tracks_serial(args, serial_runner)

    jobs, chrom_sizes = _resolve_contig_jobs(args)
    if len(jobs) <= 1:
        return _run_resumable_tracks_serial(args, serial_runner)
    specs, output_root, base_name = _tracks_relative_specs(args)
    work_root = (
        Path(args.parallel_dir).resolve()
        if getattr(args, "parallel_dir", None)
        else output_root.parent / f"{base_name}_multicontig"
    )
    per_root = work_root / "per_contig"
    per_root.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)
    namespace_data = {
        key: value
        for key, value in vars(args).items()
        if key not in {"command_runner", "command_function"}
    }
    entries: list[dict[str, object]] = []
    failures: list[str] = []
    order = {contig: index for index, (contig, _specs, _length) in enumerate(jobs)}
    print(
        f"Processing {len(jobs)} contigs with {min(cores, len(jobs))} "
        "shared-track worker process(es)."
    )
    with ProcessPoolExecutor(
        max_workers=min(cores, len(jobs)), initializer=_worker_initializer
    ) as executor:
        futures = {}
        for contig, contig_specs, _length in jobs:
            safe = _safe_contig(contig)
            worker_root = per_root / safe
            worker_root.mkdir(parents=True, exist_ok=True)
            worker_spec = worker_root / "track_spec.tsv"
            _write_tracks_worker_spec(specs, output_root, worker_root, worker_spec)
            checkpoint = _checkpoint_path(worker_root, "tracks")
            signature = _worker_signature(
                "tracks", namespace_data, contig_specs=contig_specs
            )
            report_path = None
            if getattr(args, "report", None):
                report_source = Path(str(args.report))
                report_absolute = (
                    report_source.resolve()
                    if report_source.is_absolute()
                    else (Path.cwd() / report_source).resolve()
                )
                try:
                    report_relative = report_absolute.relative_to(output_root)
                except ValueError:
                    report_relative = Path(report_source.name)
                report_path = str(worker_root / report_relative)
            if (
                not bool(getattr(args, "force", False))
                and _checkpoint_reusable(checkpoint, signature)
            ):
                print(f"Reused completed {contig}", flush=True)
                entries.append(
                    {
                        "contig": contig,
                        "directory": str(worker_root),
                        "exit_code": 0,
                        "status": "reused",
                        "checkpoint": str(checkpoint),
                    }
                )
                continue
            checkpoint.unlink(missing_ok=True)
            print(
                f"Starting tracks: {contig}; output directory {worker_root}",
                flush=True,
            )
            future = executor.submit(
                _tracks_worker,
                namespace_data,
                contig_specs,
                str(worker_root),
                str(worker_spec),
                report_path,
            )
            futures[future] = (contig, worker_root, checkpoint, signature)
        for future in as_completed(futures):
            contig, worker_root, checkpoint, signature = futures[future]
            try:
                code = int(future.result())
                if code:
                    failures.append(f"{contig}: exit code {code}")
                else:
                    _write_checkpoint(
                        checkpoint,
                        command="tracks",
                        contig=contig,
                        signature=signature,
                    )
                    print(f"Completed {contig}", flush=True)
                entries.append(
                    {
                        "contig": contig,
                        "directory": str(worker_root),
                        "exit_code": code,
                        "status": "completed" if code == 0 else "failed",
                        "checkpoint": str(checkpoint),
                    }
                )
            except Exception as error:
                failures.append(f"{contig}: {error}")
                entries.append(
                    {
                        "contig": contig,
                        "directory": str(worker_root),
                        "exit_code": 2,
                        "error": str(error),
                    }
                )
    entries.sort(key=lambda item: order[str(item["contig"])])
    manifest = {
        "schema_version": 1,
        "command": "tracks",
        "combine_strategy": "directory_tree",
        "base_name": base_name,
        "root_dir": str(work_root),
        "per_contig_dir": str(per_root),
        "combined_dir": str(output_root),
        "chrom_sizes": [{"chrom": chrom, "size": size} for chrom, size in chrom_sizes],
        "per_contig": entries,
        "directory_tree": {
            "per_contig_dirs": [str(item["directory"]) for item in entries],
            "include_roots": None,
            "exclude_parts": ["logs", ".done"],
            "sample_name": None,
            "bigwig_method": "direct",
            "strict_complete": True,
            "combine_tracks": not bool(getattr(args, "skip_combined_tracks", False)),
            "cores": cores,
            "streaming_cores": streaming_cores,
            "indexed_cores": indexed_cores,
            "bigwig_chunk_size": combine_chunk_bp,
        },
        "options": {
            "cores": cores,
            "streaming_combine_cores": streaming_cores,
            "indexed_combine_cores": indexed_cores,
            "combine_chunk_bp": combine_chunk_bp,
            "skip_combine": bool(getattr(args, "skip_combine", False)),
            "skip_combined_tracks": bool(getattr(args, "skip_combined_tracks", False)),
        },
    }
    with (work_root / MANIFEST_NAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    if failures:
        raise RuntimeError("Per-contig tracks jobs failed: " + "; ".join(failures))
    if getattr(args, "skip_combine", False):
        print(f"Per-contig outputs: {per_root}")
        print(f"Combine later with: nucleosuite combine --input-dir {work_root} --cores {cores}")
        return 0

    from nucleosuite.combine import combine_directory_trees

    result = combine_directory_trees(
        [item["directory"] for item in entries],
        output_root,
        chrom_sizes=chrom_sizes,
        bigwig_method="direct",
        strict_complete=True,
        combine_tracks=not bool(getattr(args, "skip_combined_tracks", False)),
        cores=cores,
        streaming_cores=streaming_cores,
        indexed_cores=indexed_cores,
        bigwig_chunk_size=combine_chunk_bp,
        force=bool(getattr(args, "force", False)),
    )
    print(f"Combined outputs: {result['output_dir']}")
    for warning in result["warnings"]:
        print(f"WARNING: {warning}")
    return 0 if not result["warnings"] else 1

_BIGWIG_MODULES = {
    "call-peaks": "nucleosuite.cli.call_peaks",
    "positive-runs": "nucleosuite.positive_runs",
}


def _bigwig_worker(
    command: str,
    namespace_data: dict,
    contig: str,
    requested_prefix: str,
) -> tuple[int, str]:
    import importlib

    module = importlib.import_module(_BIGWIG_MODULES[command])
    data = dict(namespace_data)
    if command == "call-peaks":
        data["regions"] = [contig]
        data["out_prefix"] = requested_prefix
    else:
        data["contigs"] = [contig]
        data["output_prefix"] = requested_prefix
    data["cores"] = 1
    data["parallel_dir"] = None
    data["skip_combine"] = True
    data["skip_combined_tracks"] = True
    data["_per_contig_worker"] = True
    namespace = argparse.Namespace(**data)
    return int(module.run(namespace) or 0), requested_prefix


def _resolve_bigwig_jobs(path: str, tokens: Sequence[str] | None) -> tuple[list[tuple[str, list[str], int]], list[tuple[str, int]]]:
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError("BigWig multicontig processing requires pyBigWig") from exc
    from nucleosuite.core.regions import expand_contig_tokens

    handle = pyBigWig.open(str(path))
    if handle is None:
        raise OSError(f"Could not open BigWig: {path}")
    try:
        sizes = {str(name): int(length) for name, length in handle.chroms().items()}
    finally:
        handle.close()
    specs = expand_contig_tokens(tokens, list(sizes))
    by_contig: dict[str, list[str]] = {}
    for spec in specs:
        by_contig.setdefault(spec.split(":", 1)[0], []).append(spec)
    jobs = [(contig, by_contig[contig], sizes[contig]) for contig in sizes if contig in by_contig]
    return jobs, [(contig, size) for contig, _specs, size in jobs]


def run_bigwig_per_contig(
    command: str,
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
    *,
    bigwig_attr: str,
    selector_attr: str,
    prefix_attr: str,
) -> int:
    cores = int(getattr(args, "cores", 1) or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    if getattr(args, "_per_contig_worker", False) or cores == 1:
        return int(serial_runner(args) or 0)

    jobs, chrom_sizes = _resolve_bigwig_jobs(
        str(getattr(args, bigwig_attr)), getattr(args, selector_attr, None)
    )
    if len(jobs) <= 1:
        return int(serial_runner(args) or 0)

    explicit_prefix = getattr(args, prefix_attr, None)
    if explicit_prefix:
        base_path = Path(str(explicit_prefix))
    else:
        name = Path(str(getattr(args, bigwig_attr))).name
        for suffix in (".bigWig", ".bw"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        base_path = Path(name)
    root = (
        Path(args.parallel_dir)
        if getattr(args, "parallel_dir", None)
        else base_path.parent / f"{base_path.name}_multicontig"
    ).resolve()
    per_contig_root = root / "per_contig"
    combined_dir = root / "combined"
    per_contig_root.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)

    namespace_data = {
        key: value
        for key, value in vars(args).items()
        if key not in {"command_runner", "command_function"}
    }
    entries: list[dict[str, object]] = []
    failures: list[str] = []
    with ProcessPoolExecutor(
        max_workers=min(cores, len(jobs)), initializer=_worker_initializer
    ) as executor:
        futures = {}
        for contig, _specs, _length in jobs:
            safe = _safe_contig(contig)
            directory = per_contig_root / safe
            directory.mkdir(parents=True, exist_ok=True)
            prefix = str(directory / f"{base_path.name}_{safe}")
            future = executor.submit(_bigwig_worker, command, namespace_data, contig, prefix)
            futures[future] = (contig, prefix)
        for future in as_completed(futures):
            contig, prefix = futures[future]
            try:
                exit_code, actual_prefix = future.result()
                if exit_code:
                    failures.append(f"{contig}: exit code {exit_code}")
                entries.append({"contig": contig, "prefix": actual_prefix, "exit_code": exit_code})
                if exit_code == 0:
                    print(f"Completed {contig}")
            except Exception as error:
                failures.append(f"{contig}: {error}")
                entries.append({"contig": contig, "prefix": prefix, "exit_code": 2, "error": str(error)})

    order = {contig: index for index, (contig, _specs, _length) in enumerate(jobs)}
    entries.sort(key=lambda item: order[str(item["contig"])])
    manifest = {
        "schema_version": 1,
        "command": command,
        "base_name": base_path.name,
        "root_dir": str(root),
        "per_contig_dir": str(per_contig_root),
        "combined_dir": str(combined_dir),
        "combined_name": base_path.name,
        "chrom_sizes": [{"chrom": chrom, "size": size} for chrom, size in chrom_sizes],
        "per_contig": entries,
        "options": {
            "cores": cores,
            "skip_combine": bool(getattr(args, "skip_combine", False)),
            "skip_combined_tracks": bool(getattr(args, "skip_combined_tracks", False)),
        },
    }
    with (root / MANIFEST_NAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    if failures:
        raise RuntimeError("Per-contig jobs failed: " + "; ".join(failures))
    if not getattr(args, "skip_combine", False):
        from nucleosuite.combine import combine_run
        result = combine_run(
            root,
            combine_tracks=not getattr(args, "skip_combined_tracks", False),
            cores=cores,
        )
        print(f"Combined outputs: {result['output_dir']}")
        for warning in result["warnings"]:
            print(f"WARNING: {warning}")
    return 0

_REGION_MODULES = {
    "dac": "nucleosuite.dac",
    "dcc": "nucleosuite.dcc",
}


def _read_contigs_from_table(path: str) -> list[str]:
    import gzip

    output: list[str] = []
    seen: set[str] = set()
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            contig = line.split()[0]
            if contig not in seen:
                seen.add(contig)
                output.append(contig)
    return output


def _resolve_region_contigs(args: argparse.Namespace) -> tuple[list[tuple[str, list[str], int]], list[tuple[str, int]]]:
    sizes: dict[str, int] = {}
    if getattr(args, "chrom_sizes", None):
        sizes = dict(read_chrom_sizes_source(args.chrom_sizes))
        available = list(sizes)
    else:
        path = getattr(args, "regions_bed", None) or getattr(args, "genes_bed", None)
        if not path:
            raise ValueError("A chromosome sizes or BED region source is required")
        available = _read_contigs_from_table(path)
        # Sizes are only needed for BigWig combination. Resolve from input BigWigs.
        bigwig_values = (
            getattr(args, "bigwig", None)
            or getattr(args, "bigwig_a", None)
            or []
        )
        try:
            import glob
            import pyBigWig
            expanded: list[str] = []
            for value in bigwig_values:
                expanded.extend(sorted(glob.glob(value)) or [value])
            for value in expanded:
                handle = pyBigWig.open(str(value))
                if handle is None:
                    continue
                try:
                    for chrom, size in handle.chroms().items():
                        sizes.setdefault(str(chrom), int(size))
                finally:
                    handle.close()
        except ImportError:
            pass
    selected = getattr(args, "chromosome", None)
    if selected:
        requested: list[str] = []
        for raw in selected:
            requested.extend(item.strip() for item in str(raw).split(",") if item.strip())
        available = [chrom for chrom in available if chrom in set(requested)]
    jobs = [(chrom, [chrom], int(sizes.get(chrom, 1))) for chrom in available]
    return jobs, [(chrom, int(sizes.get(chrom, 1))) for chrom in available]


def _region_worker(command: str, namespace_data: dict, contig: str, output_dir: str, prefix: str) -> tuple[int, str]:
    import importlib

    module = importlib.import_module(_REGION_MODULES[command])
    data = dict(namespace_data)
    data["chromosome"] = [contig]
    data["scope"] = "chromosome"
    data["output_dir"] = output_dir
    data["out_prefix"] = prefix
    data["cores"] = 1
    data["parallel_dir"] = None
    data["skip_combine"] = True
    data["skip_combined_tracks"] = True
    data["_per_contig_worker"] = True
    namespace = argparse.Namespace(**data)
    return int(module.run(namespace) or 0), str(Path(output_dir) / prefix)


def run_region_per_contig(
    command: str,
    args: argparse.Namespace,
    serial_runner: Callable[[argparse.Namespace], int],
) -> int:
    cores = int(getattr(args, "cores", 1) or 1)
    if cores < 1:
        raise ValueError("--cores must be at least 1")
    if getattr(args, "_per_contig_worker", False) or cores == 1:
        return int(serial_runner(args) or 0)
    jobs, chrom_sizes = _resolve_region_contigs(args)
    if len(jobs) <= 1:
        return int(serial_runner(args) or 0)

    base_name = str(getattr(args, "out_prefix", None) or command)
    base_path = Path(base_name)
    root = (
        Path(args.parallel_dir)
        if getattr(args, "parallel_dir", None)
        else Path(getattr(args, "output_dir", ".")) / f"{base_path.name}_multicontig"
    ).resolve()
    per_contig_root = root / "per_contig"
    combined_dir = root / "combined"
    per_contig_root.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    namespace_data = {
        key: value for key, value in vars(args).items()
        if key not in {"command_runner", "command_function"}
    }
    entries: list[dict[str, object]] = []
    failures: list[str] = []
    with ProcessPoolExecutor(
        max_workers=min(cores, len(jobs)), initializer=_worker_initializer
    ) as executor:
        futures = {}
        for contig, _specs, _size in jobs:
            safe = _safe_contig(contig)
            directory = per_contig_root / safe
            directory.mkdir(parents=True, exist_ok=True)
            prefix = f"{base_path.name}_{safe}"
            future = executor.submit(_region_worker, command, namespace_data, contig, str(directory), prefix)
            futures[future] = (contig, directory / prefix)
        for future in as_completed(futures):
            contig, expected_prefix = futures[future]
            try:
                exit_code, actual_prefix = future.result()
                if exit_code:
                    failures.append(f"{contig}: exit code {exit_code}")
                entries.append({"contig": contig, "prefix": actual_prefix, "exit_code": exit_code})
                if exit_code == 0:
                    print(f"Completed {contig}")
            except Exception as error:
                failures.append(f"{contig}: {error}")
                entries.append({"contig": contig, "prefix": str(expected_prefix), "exit_code": 2, "error": str(error)})
    order = {contig: index for index, (contig, _specs, _size) in enumerate(jobs)}
    entries.sort(key=lambda item: order[str(item["contig"])])
    manifest = {
        "schema_version": 1,
        "command": command,
        "base_name": base_path.name,
        "root_dir": str(root),
        "per_contig_dir": str(per_contig_root),
        "combined_dir": str(combined_dir),
        "combined_name": base_path.name,
        "chrom_sizes": [{"chrom": chrom, "size": size} for chrom, size in chrom_sizes],
        "per_contig": entries,
        "options": {
            "cores": cores,
            "skip_combine": bool(getattr(args, "skip_combine", False)),
            "skip_combined_tracks": bool(getattr(args, "skip_combined_tracks", False)),
        },
    }
    with (root / MANIFEST_NAME).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    if failures:
        raise RuntimeError("Per-contig jobs failed: " + "; ".join(failures))
    if not getattr(args, "skip_combine", False):
        from nucleosuite.combine import combine_run
        result = combine_run(
            root,
            combine_tracks=not getattr(args, "skip_combined_tracks", False),
            cores=cores,
        )
        print(f"Combined outputs: {result['output_dir']}")
        for warning in result["warnings"]:
            print(f"WARNING: {warning}")
    return 0
