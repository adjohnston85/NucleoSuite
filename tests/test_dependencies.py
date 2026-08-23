from __future__ import annotations

import ast
import sys
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _normalise_conda_name(value: str) -> str:
    import re

    value = value.strip()
    if not value or value.startswith("#"):
        return ""
    return re.split(r"[<>=!~ ]", value, maxsplit=1)[0].lower()


def _environment_dependencies() -> set[str]:
    dependencies: set[str] = set()
    in_dependencies = False
    for raw in (ROOT / "environment.yml").read_text(encoding="utf-8").splitlines():
        line = raw.rstrip()
        if line == "dependencies:":
            in_dependencies = True
            continue
        if not in_dependencies:
            continue
        stripped = line.strip()
        if stripped.startswith("-"):
            name = _normalise_conda_name(stripped[1:])
            if name:
                dependencies.add(name)
    return dependencies


def _recipe_run_dependencies() -> set[str]:
    dependencies: set[str] = set()
    in_run = False
    run_indent = None
    for raw in (ROOT / "recipe" / "meta.yaml").read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if stripped == "run:":
            in_run = True
            run_indent = len(raw) - len(raw.lstrip())
            continue
        if in_run:
            indent = len(raw) - len(raw.lstrip())
            if stripped and indent <= int(run_indent or 0):
                break
            if stripped.startswith("-"):
                name = _normalise_conda_name(stripped[1:])
                if name:
                    dependencies.add(name)
    return dependencies


def _project_dependencies() -> set[str]:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    values = data["project"]["dependencies"]
    return {value.split("[")[0].split(">=")[0].split("==")[0].lower() for value in values}


def _third_party_imports() -> set[str]:
    imports: set[str] = set()
    local_roots = {"nucleosuite", "core", "intervals"}
    for path in (ROOT / "src" / "nucleosuite").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module.split(".")[0])
    return {
        name
        for name in imports
        if name not in sys.stdlib_module_names and name not in local_roots
    }


def test_environment_name_is_nucleosuite():
    first_nonempty = next(
        line.strip()
        for line in (ROOT / "environment.yml").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )
    assert first_nonempty == "name: nucleosuite"


def test_all_python_imports_are_declared_in_package_and_environment():
    import_to_project = {
        "numpy": "numpy",
        "scipy": "scipy",
        "pysam": "pysam",
        "pyBigWig": "pybigwig",
        "matplotlib": "matplotlib",
        "matplotlib_venn": "matplotlib-venn",
        "openpyxl": "openpyxl",
    }
    import_to_conda = {
        **import_to_project,
        "matplotlib": "matplotlib-base",
    }
    observed = _third_party_imports()
    assert observed <= set(import_to_project), f"Unmapped third-party imports: {sorted(observed - set(import_to_project))}"

    project = _project_dependencies()
    environment = _environment_dependencies()
    recipe = _recipe_run_dependencies()
    for module in observed:
        assert import_to_project[module] in project
        assert import_to_conda[module] in environment
        assert import_to_conda[module] in recipe


def test_required_external_tools_are_declared_in_environment_and_recipe():
    required = {
        "bash",
        "samtools",
        "ucsc-bedtobigbed",
        "ucsc-bigbedtobed",
        "ucsc-bedgraphtobigwig",
    }
    assert required <= _environment_dependencies()
    assert required <= _recipe_run_dependencies()


def test_conda_recipe_version_matches_package_version():
    import re
    import tomllib

    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    recipe = (ROOT / "recipe" / "meta.yaml").read_text(encoding="utf-8")
    match = re.search(r'\{% set version = "([^"]+)" %\}', recipe)
    assert match is not None
    assert match.group(1) == project["project"]["version"]
