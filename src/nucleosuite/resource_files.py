"""Resolve, validate, and copy reference files bundled with NucleoSuite."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections.abc import Mapping, Sequence
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any

from nucleosuite.cli.formatting import NucleoSuiteHelpFormatter


_MANIFEST_RELATIVE_PATH = "resources/manifest.json"

# Resource-name mapping available to internal callers.
# RESOURCE_FILES directly.  The manifest is now the source of truth.
RESOURCE_FILES: dict[str, str] = {}


def load_manifest() -> dict[str, Any]:
    """Load and minimally validate the bundled resource manifest."""
    resource = files("nucleosuite").joinpath(_MANIFEST_RELATIVE_PATH)
    with resource.open("rt", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict):
        raise ValueError("The bundled resource manifest must be a JSON object")
    resources = manifest.get("resources")
    resource_sets = manifest.get("resource_sets")
    if not isinstance(resources, dict) or not isinstance(resource_sets, dict):
        raise ValueError("The resource manifest requires 'resources' and 'resource_sets' objects")
    return manifest


def resource_entries() -> Mapping[str, Mapping[str, str]]:
    """Return named resource entries from the manifest."""
    return load_manifest()["resources"]


# Populate the resource mapping once at import. Values remain package-
# relative paths, matching the former public constant.
RESOURCE_FILES.update(
    {name: str(entry.get("path", "")) for name, entry in resource_entries().items()}
)


def resource_set_entries() -> Mapping[str, Mapping[str, str]]:
    """Return named resource-set entries from the manifest."""
    return load_manifest()["resource_sets"]


def resolve_resource_name(name: str) -> str:
    """Return the package-relative path for a named resource."""
    entries = resource_entries()
    if name not in entries:
        available = ", ".join(sorted(entries))
        raise KeyError(f"Unknown resource {name!r}. Available resources: {available}")
    path = entries[name].get("path")
    if not path:
        raise ValueError(f"Resource {name!r} has no path in the manifest")
    return str(path)


def resolve_set_resource_name(resource_set: str, logical_name: str) -> str:
    """Resolve a logical name such as ``genes`` or ``ctcf`` to a resource name."""
    sets = resource_set_entries()
    if resource_set not in sets:
        available = ", ".join(sorted(sets))
        raise KeyError(f"Unknown resource set {resource_set!r}. Available sets: {available}")
    entry = sets[resource_set]
    if logical_name not in entry or logical_name == "description":
        available = ", ".join(sorted(key for key in entry if key != "description"))
        raise KeyError(
            f"Resource set {resource_set!r} has no logical resource {logical_name!r}. "
            f"Available logical resources: {available}"
        )
    return str(entry[logical_name])


def resource_traversable(name: str):
    """Return the importlib Traversable for a named resource."""
    return files("nucleosuite").joinpath(resolve_resource_name(name))


def materialized_resource_path(name: str):
    """Return a context manager yielding a filesystem path for a resource."""
    return as_file(resource_traversable(name))


def validate_manifest() -> list[tuple[str, str, bool, str]]:
    """Validate all manifest resources and resource-set references."""
    manifest = load_manifest()
    rows: list[tuple[str, str, bool, str]] = []
    resources = manifest["resources"]
    for name in sorted(resources):
        relative = str(resources[name].get("path", ""))
        traversable = files("nucleosuite").joinpath(relative) if relative else None
        exists = bool(traversable is not None and traversable.is_file())
        target = relative
        expected_sha256 = str(resources[name].get("sha256", "")).lower()
        if exists and expected_sha256:
            digest = hashlib.sha256()
            with traversable.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            observed_sha256 = digest.hexdigest()
            exists = observed_sha256 == expected_sha256
            target = (
                f"{relative} sha256={observed_sha256}"
                if exists
                else f"{relative} expected_sha256={expected_sha256} observed_sha256={observed_sha256}"
            )
        rows.append(("resource", name, exists, target))

    for set_name, entry in sorted(manifest["resource_sets"].items()):
        for logical, resource_name in sorted(entry.items()):
            if logical == "description":
                continue
            exists = resource_name in resources
            rows.append((f"set:{set_name}", logical, exists, str(resource_name)))
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nucleosuite resources",
        description="List, resolve, validate, or copy files bundled with NucleoSuite.",
        formatter_class=NucleoSuiteHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="action", required=True)

    list_parser = subparsers.add_parser("list", help="List resources and resource sets.")
    list_parser.add_argument("--resource-set", help="Limit output to one resource set.")

    path_parser = subparsers.add_parser("path", help="Print the installed path of one named resource.")
    path_parser.add_argument("name", help="Installed resource name shown by 'resources list'.")

    show_parser = subparsers.add_parser(
        "show", help="Print the installed path of one logical resource from a resource set."
    )
    show_parser.add_argument("name", help="Logical resource name, for example genes, states, or ctcf.")
    show_parser.add_argument("--resource-set", required=True, help="Resource set containing the requested logical name.")

    validate_parser = subparsers.add_parser("validate", help="Validate the installed manifest and files.")
    validate_parser.add_argument("--resource-set", help="Also require all entries in this resource set.")

    copy_parser = subparsers.add_parser("copy", help="Copy one or more resources to a directory.")
    copy_parser.add_argument("--output-dir", required=True, help="Destination directory for copied resource files.")
    copy_parser.add_argument("--name", action="append", default=[], help="Installed resource name to copy; repeat for multiple resources.")
    copy_parser.add_argument("--resource-set", help="Copy every resource referenced by this resource set.")
    return parser


def _print_list(resource_set: str | None) -> int:
    manifest = load_manifest()
    resources = manifest["resources"]
    sets = manifest["resource_sets"]
    if resource_set:
        if resource_set not in sets:
            raise KeyError(f"Unknown resource set {resource_set!r}")
        print("resource_set\tlogical_name\tresource_name\tdescription")
        for logical, resource_name in sets[resource_set].items():
            if logical == "description":
                continue
            description = resources[resource_name].get("description", "")
            print(f"{resource_set}\t{logical}\t{resource_name}\t{description}")
        return 0

    print("type\tname\tdescription")
    for name, entry in resources.items():
        print(f"resource\t{name}\t{entry.get('description', '')}")
    for name, entry in sets.items():
        print(f"resource_set\t{name}\t{entry.get('description', '')}")
    return 0


def _print_path(resource_name: str) -> int:
    resource = resource_traversable(resource_name)
    if not resource.is_file():
        raise FileNotFoundError(
            f"Bundled resource {resource_name!r} is missing: {resolve_resource_name(resource_name)}"
        )
    with as_file(resource) as path:
        print(path)
    return 0


def _copy_resources(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    selected = list(args.name)
    if args.resource_set:
        entry = resource_set_entries().get(args.resource_set)
        if entry is None:
            raise KeyError(f"Unknown resource set {args.resource_set!r}")
        selected.extend(str(value) for key, value in entry.items() if key != "description")
    if not selected:
        selected = list(resource_entries())

    for name in dict.fromkeys(selected):
        resource = resource_traversable(name)
        if not resource.is_file():
            raise FileNotFoundError(f"Bundled resource {name!r} is missing")
        with as_file(resource) as path:
            destination = output_dir / Path(resolve_resource_name(name)).name
            shutil.copy2(path, destination)
            print(f"{name}\t{destination}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.action == "list":
        return _print_list(args.resource_set)
    if args.action == "path":
        return _print_path(args.name)
    if args.action == "show":
        return _print_path(resolve_set_resource_name(args.resource_set, args.name))
    if args.action == "copy":
        return _copy_resources(args)

    rows = validate_manifest()
    failed = [row for row in rows if not row[2]]
    if args.resource_set:
        set_entry = resource_set_entries().get(args.resource_set)
        if set_entry is None:
            raise KeyError(f"Unknown resource set {args.resource_set!r}")
        required = {str(value) for key, value in set_entry.items() if key != "description"}
        resources = resource_entries()
        for name in required:
            if name not in resources:
                failed.append((f"set:{args.resource_set}", name, False, "missing manifest resource"))
    print("scope\tname\tstatus\ttarget")
    for scope, name, exists, target in rows:
        print(f"{scope}\t{name}\t{'OK' if exists else 'MISSING'}\t{target}")
    if failed:
        raise FileNotFoundError(
            "Resource validation failed for: "
            + ", ".join(f"{scope}/{name}" for scope, name, _, _ in failed)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
