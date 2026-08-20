from __future__ import annotations

import argparse
import importlib
import re
from pathlib import Path

from nucleosuite.cli.main import DELEGATED_MODULES, build_parser


ROOT = Path(__file__).resolve().parents[1]


def test_required_user_guides_exist() -> None:
    required = [
        ROOT / "README.md",
        ROOT / "docs" / "README.md",
        ROOT / "docs" / "QUICKSTART.md",
        ROOT / "docs" / "CHOOSING_A_COMMAND.md",
        ROOT / "docs" / "WORKFLOWS.md",
        ROOT / "docs" / "COMMAND_REFERENCE.md",
        ROOT / "docs" / "GLOSSARY.md",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
    assert not missing, f"Missing documentation files: {missing}"


def test_every_command_page_explains_what_and_why() -> None:
    command_pages = sorted((ROOT / "docs" / "commands").glob("*.md"))
    assert command_pages
    missing_what = [
        path.name for path in command_pages
        if "## What this command does" not in path.read_text()
    ]
    missing_why = [
        path.name for path in command_pages
        if "## Why use" not in path.read_text()
    ]
    assert not missing_what, f"Command pages missing 'What this command does': {missing_what}"
    assert not missing_why, f"Command pages missing 'Why use it': {missing_why}"


def test_documentation_demonstrates_composable_resource_paths() -> None:
    expected = {
        ROOT / "docs" / "commands" / "resources.md": [
            "nucleosuite resources path gm12878-hg19-states",
            "nucleosuite resources path gm12878-hg19-ctcf",
        ],
        ROOT / "docs" / "commands" / "distances.md": [
            '"$(nucleosuite resources path gm12878-hg19-states)"',
        ],
        ROOT / "docs" / "commands" / "aggregate.md": [
            '"$(nucleosuite resources path gm12878-hg19-ctcf)"',
        ],
        ROOT / "docs" / "QUICKSTART.md": [
            "nucleosuite resources path gm12878-hg19-states",
        ],
    }
    missing: list[str] = []
    for path, snippets in expected.items():
        content = path.read_text()
        for snippet in snippets:
            if snippet not in content:
                missing.append(f"{path.relative_to(ROOT)}: {snippet}")
    assert not missing, "Missing resource-path usage examples:\n" + "\n".join(missing)


def test_user_documentation_does_not_announce_plain_language_style() -> None:
    markdown_files = [ROOT / "README.md", *(ROOT / "docs").rglob("*.md")]
    patterns = [
        re.compile(r"\bin plain language\b", re.I),
        re.compile(r"\bput simply\b", re.I),
        re.compile(r"\bin simple terms\b", re.I),
        re.compile(r"\bsimply stated\b", re.I),
    ]
    failures: list[str] = []
    for markdown in markdown_files:
        text = markdown.read_text()
        for pattern in patterns:
            if pattern.search(text):
                failures.append(str(markdown.relative_to(ROOT)))
                break
    assert not failures, "Documentation contains meta-language labels: " + ", ".join(failures)


def test_workflow_diagrams_are_present() -> None:
    workflows = (ROOT / "docs" / "WORKFLOWS.md").read_text()
    assert workflows.count("```mermaid") >= 10
    for page in (
        ROOT / "docs" / "commands" / "cfdna-suite.md",
        ROOT / "docs" / "commands" / "mnase-suite.md",
        ROOT / "docs" / "commands" / "tracks.md",
    ):
        assert "```mermaid" in page.read_text(), f"Missing workflow diagram: {page.relative_to(ROOT)}"


def test_algorithm_figure_assets_exist_and_are_linked() -> None:
    algorithms = (ROOT / "docs" / "ALGORITHMS.md").read_text()
    figures = [
        "pns_kernels_120_167_180_multipanel_single_legend.png",
        "bns_kernels_120_167_180_mode167.png",
        "tns_kernels_120_167_180_mode167.png",
        "wps_kernels_120_167_180_multiplot.png",
        "dac_periodicity_example.png",
    ]
    missing: list[str] = []
    for filename in figures:
        path = ROOT / "docs" / "images" / filename
        if not path.is_file():
            missing.append(f"missing file: {path.relative_to(ROOT)}")
        if f"images/{filename}" not in algorithms:
            missing.append(f"not linked from ALGORITHMS.md: {filename}")
    assert not missing, "Algorithm figure problems:\n" + "\n".join(missing)


def test_every_primary_command_has_one_documentation_page() -> None:
    parser = build_parser()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    aliases = {"dyad", "peak-call", "peakcall"}
    commands = set(subparsers.choices) - aliases
    pages = {path.stem for path in (ROOT / "docs" / "commands").glob("*.md")}
    assert pages == commands, (
        f"Undocumented commands: {sorted(commands - pages)}; "
        f"orphan command pages: {sorted(pages - commands)}"
    )


def test_command_reference_links_every_command_page() -> None:
    reference = (ROOT / "docs" / "COMMAND_REFERENCE.md").read_text()
    missing = [
        path.name for path in sorted((ROOT / "docs" / "commands").glob("*.md"))
        if f"(commands/{path.name})" not in reference
    ]
    assert not missing, f"Command pages missing from COMMAND_REFERENCE.md: {missing}"


def _primary_command_names() -> set[str]:
    parser = build_parser()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    aliases = {"dyad", "peak-call", "peakcall"}
    return set(subparsers.choices) - aliases


def test_main_readme_links_every_primary_command() -> None:
    readme = (ROOT / "README.md").read_text()
    commands = _primary_command_names()
    missing = [
        command for command in sorted(commands)
        if f"(docs/commands/{command}.md)" not in readme
    ]
    assert not missing, f"Primary commands missing from README.md: {missing}"


def test_command_pages_have_matching_heading_and_back_link() -> None:
    failures: list[str] = []
    for path in sorted((ROOT / "docs" / "commands").glob("*.md")):
        text = path.read_text()
        expected_heading = f"# `nucleosuite {path.stem}`"
        first_line = text.splitlines()[0] if text.splitlines() else ""
        if first_line != expected_heading:
            failures.append(f"{path.name}: heading is {first_line!r}")
        if "../COMMAND_REFERENCE.md" not in text:
            failures.append(f"{path.name}: missing command-reference back link")
    assert not failures, "Command-page navigation problems:\n" + "\n".join(failures)


def test_flank_spacing_is_discoverable_in_relevant_guides() -> None:
    required = [
        ROOT / "README.md",
        ROOT / "docs" / "README.md",
        ROOT / "docs" / "CHOOSING_A_COMMAND.md",
        ROOT / "docs" / "COMMAND_REFERENCE.md",
        ROOT / "docs" / "QUICKSTART.md",
        ROOT / "docs" / "WORKFLOWS.md",
        ROOT / "docs" / "ALGORITHMS.md",
        ROOT / "docs" / "FILE_FORMATS.md",
        ROOT / "docs" / "GLOSSARY.md",
        ROOT / "docs" / "PLOTTING.md",
    ]
    missing = [
        str(path.relative_to(ROOT)) for path in required
        if "flank-spacing" not in path.read_text()
    ]
    assert not missing, f"flank-spacing missing from relevant guides: {missing}"


def _all_parser_options(parser: argparse.ArgumentParser) -> set[str]:
    options: set[str] = set()
    seen: set[int] = set()

    def walk(current: argparse.ArgumentParser) -> None:
        if id(current) in seen:
            return
        seen.add(id(current))
        for action in current._actions:
            options.update(action.option_strings)
            if isinstance(action, argparse._SubParsersAction):
                for child in action.choices.values():
                    walk(child)

    walk(parser)
    return options


def _documented_command_option_map() -> dict[str, set[str]]:
    root_parser = build_parser()
    root_subparsers = next(
        action for action in root_parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    option_map = {
        command: _all_parser_options(parser)
        for command, parser in root_subparsers.choices.items()
    }

    for command, module_name in DELEGATED_MODULES.items():
        if command in {"mnase-suite", "cfdna-suite"}:
            continue
        module = importlib.import_module(module_name)
        builder = getattr(module, "build_parser", None)
        if builder is not None:
            option_map[command] = _all_parser_options(builder())

    # `nucleosuite plot` uses a two-stage parser: plot-family options are
    # exposed only after the input/metadata identifies the relevant renderer.
    from nucleosuite import replot as replot_module
    dynamic_plot_options = set(option_map.get("plot", set()))
    for plot_type in replot_module.PLOT_TYPES:
        if plot_type == "auto":
            continue
        dynamic_plot_options.update(_all_parser_options(replot_module.build_parser(plot_type)))
    option_map["plot"] = dynamic_plot_options

    from nucleosuite.cli.aggregate import add_aggregate_parser

    aggregate_root = argparse.ArgumentParser(prog="nucleosuite")
    aggregate_subparsers = aggregate_root.add_subparsers(dest="command", required=True)
    aggregate_parser = add_aggregate_parser(aggregate_subparsers)
    option_map["aggregate"] = _all_parser_options(aggregate_parser)

    for alias, canonical in {
        "dyad": "dyads",
        "peak-call": "call-peaks",
        "peakcall": "call-peaks",
        "region-peak-extractor": "region-extract",
        "pns-region-extractor": "region-extract",
    }.items():
        option_map[alias] = option_map[canonical]
    return option_map


def test_fenced_nucleosuite_examples_use_current_cli_options() -> None:
    option_map = _documented_command_option_map()
    markdown_files = [ROOT / "README.md", *(ROOT / "docs").rglob("*.md")]
    failures: list[str] = []
    for markdown in markdown_files:
        text = markdown.read_text()
        for block in re.findall(r"(?ms)^```(?:bash|sh)\s*\n(.*?)\n```", text):
            normalized = re.sub(r"\\\s*\n", " ", block)
            for match in re.finditer(
                r"(?<![$(\w-])nucleosuite\s+([a-z0-9-]+)\s+([^\n;]*)",
                normalized,
            ):
                command = match.group(1)
                if command in {"mnase-suite", "cfdna-suite"} or command not in option_map:
                    continue
                used = set(re.findall(r"(?<![\w-])--[a-z0-9-]+", match.group(2)))
                unknown = sorted(used - option_map[command])
                if unknown:
                    failures.append(
                        f"{markdown.relative_to(ROOT)}: nucleosuite {command}: {unknown}"
                    )
    assert not failures, "Documentation examples use unsupported CLI options:\n" + "\n".join(failures)


def test_markdown_has_no_literal_newline_escape_artifacts() -> None:
    failures: list[str] = []
    for markdown in [ROOT / "README.md", *(ROOT / "docs").rglob("*.md")]:
        text = markdown.read_text()
        if re.search(r"\\n(?:#|```|[-*] |[A-Za-z])", text):
            failures.append(str(markdown.relative_to(ROOT)))
    assert not failures, "Markdown contains literal \\n escape artifacts: " + ", ".join(failures)


def test_local_markdown_links_resolve() -> None:
    markdown_files = [ROOT / "README.md"]
    markdown_files.extend((ROOT / "docs").rglob("*.md"))
    pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
    broken: list[str] = []
    for markdown in markdown_files:
        for target in pattern.findall(markdown.read_text()):
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            target_path = target.split("#", 1)[0]
            if not target_path:
                continue
            resolved = (markdown.parent / target_path).resolve()
            if not resolved.exists():
                broken.append(f"{markdown.relative_to(ROOT)} -> {target}")
    assert not broken, "Broken local Markdown links:\n" + "\n".join(broken)


def _walk_parsers(
    name: str,
    parser: argparse.ArgumentParser,
    seen: set[int],
) -> list[str]:
    if id(parser) in seen:
        return []
    seen.add(id(parser))
    missing: list[str] = []
    for action in parser._actions:
        if action.option_strings and action.help is None:
            missing.append(f"{name}: {', '.join(action.option_strings)}")
        if isinstance(action, argparse._SubParsersAction):
            for child_name, child in action.choices.items():
                missing.extend(_walk_parsers(f"{name}/{child_name}", child, seen))
    return missing


def test_every_argparse_option_has_help_text() -> None:
    parsers: list[tuple[str, argparse.ArgumentParser]] = [("nucleosuite", build_parser())]
    for command, module_name in DELEGATED_MODULES.items():
        if command in {"mnase-suite", "cfdna-suite"}:
            continue
        module = importlib.import_module(module_name)
        builder = getattr(module, "build_parser", None)
        if builder is not None:
            parsers.append((command, builder()))

    from nucleosuite.cli.aggregate import add_aggregate_parser

    aggregate_root = argparse.ArgumentParser(prog="nucleosuite")
    aggregate_subparsers = aggregate_root.add_subparsers(dest="command", required=True)
    add_aggregate_parser(aggregate_subparsers)
    parsers.append(("aggregate", aggregate_root))

    missing: list[str] = []
    for name, parser in parsers:
        missing.extend(_walk_parsers(name, parser, set()))
    assert not missing, "CLI options without help text:\n" + "\n".join(missing)


def test_suite_help_covers_every_public_shell_option() -> None:
    resource_dir = ROOT / "src" / "nucleosuite" / "resources"
    internal_options = {
        "--analysis-chrom-sizes-source",
        "--randomized-control-input",
        "--provenance-bam",
        "--provenance-fragment",
        "--trusted-combined-prerequisites",
        "--combine-prerequisites-only",
        "--validate-only",
        "--broad-frag-lower",
        "--broad-frag-upper",
    }
    failures: list[str] = []
    for filename in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = (resource_dir / filename).read_text()
        usage_body = text.split("\nusage() {", 1)[1]
        usage = usage_body.split("cat <<'EOF'", 1)[1].split("\nEOF", 1)[0]
        parser_body = text.split("while [[ $# -gt 0 ]]", 1)[1]
        parser_body = parser_body.split("\n    esac", 1)[0].split("\n  esac", 1)[0]
        accepted = set(re.findall(r"(?m)^\s*(--[a-z0-9-]+)(?=[)=])", parser_body))
        documented = set(re.findall(r"(?<![\w-])--[a-z][a-z0-9-]*", usage))
        missing = sorted(accepted - documented - internal_options)
        failures.extend(f"{filename}: {option}" for option in missing)
    assert not failures, "Public suite options missing from --help:\n" + "\n".join(failures)


def test_user_documentation_contains_only_current_method_language() -> None:
    markdown_files = [ROOT / "README.md", ROOT / "src" / "nucleosuite" / "resources" / "README.md"]
    markdown_files.extend((ROOT / "docs").rglob("*.md"))
    banned = {
        "stale Alignment glossary entry": re.compile(r"(?m)^## Alignment\s*$"),
        "misleading positive-only PNS": re.compile(r"positive-only PNS", re.I),
        "misleading WPS local prominence": re.compile(r"local prominence", re.I),
        "shift metaphor for DCC": re.compile(r"one signal is shifted relative", re.I),
        "hidden compatibility option": re.compile(r"max-anchor-tries", re.I),
        "old leftover-set semantics": re.compile(
            r"(?:overlaps involving repression.*(?:become|fall into).*leftover|"
            r"leftover-set-name.*not assigned to another final set)",
            re.I,
        ),
        "removed suite behavior": re.compile(r"intentionally (?:omits|does not calculate)", re.I),
        "breakpoint-like terminology": re.compile(r"breakpoint-like", re.I),
        "PNS reference-boundary caveat": re.compile(
            r"(?:reference-boundary clipping|clipping at reference boundaries)", re.I
        ),
        "PNS unclipped-fragment caveat": re.compile(r"complete, unclipped", re.I),
    }
    failures: list[str] = []
    for markdown in markdown_files:
        text = markdown.read_text()
        for label, pattern in banned.items():
            if pattern.search(text):
                failures.append(f"{markdown.relative_to(ROOT)}: {label}")
    assert not failures, "Stale or misleading documentation language:\n" + "\n".join(failures)


def test_documentation_math_uses_supported_fenced_syntax() -> None:
    markdown_files = [ROOT / "README.md", *(ROOT / "docs").rglob("*.md")]
    incompatible_delimiters: list[str] = []
    malformed: list[str] = []
    unsupported: list[str] = []
    for markdown in markdown_files:
        text = markdown.read_text()
        if "$$" in text or re.search(r"\\[()\[\]]", text):
            incompatible_delimiters.append(str(markdown.relative_to(ROOT)))

        open_fence: tuple[int, str] | None = None
        for line_number, line in enumerate(text.splitlines(), start=1):
            if not line.startswith("```"):
                continue
            fence_info = line[3:].strip()
            if open_fence is None:
                open_fence = (line_number, fence_info)
            else:
                if fence_info:
                    malformed.append(
                        f"{markdown.relative_to(ROOT)}:{line_number}: "
                        "closing fence must not contain an info string"
                    )
                open_fence = None
        if open_fence is not None:
            malformed.append(
                f"{markdown.relative_to(ROOT)}:{open_fence[0]}: unclosed code fence"
            )

        for match in re.finditer(r"\\(operatorname|begin|end)\b", text):
            line_number = text.count("\n", 0, match.start()) + 1
            unsupported.append(
                f"{markdown.relative_to(ROOT)}:{line_number}: \\{match.group(1)}"
            )

        blocks = list(re.finditer(r"(?ms)^```math[ \t]*\n(.*?)\n```[ \t]*$", text))
        opening_fences = len(re.findall(r"(?m)^```math[ \t]*$", text))
        if opening_fences != len(blocks):
            malformed.append(f"{markdown.relative_to(ROOT)}: unclosed math fence")

        for block in blocks:
            formula = block.group(1)
            environment_stack: list[str] = []
            for action, environment in re.findall(r"\\(begin|end)\{([^{}]+)\}", formula):
                if action == "begin":
                    environment_stack.append(environment)
                elif not environment_stack or environment_stack.pop() != environment:
                    malformed.append(
                        f"{markdown.relative_to(ROOT)}: unmatched \\end{{{environment}}}"
                    )
            if environment_stack:
                malformed.append(
                    f"{markdown.relative_to(ROOT)}: unclosed \\begin{{{environment_stack[-1]}}}"
                )

            brace_depth = 0
            for index, character in enumerate(formula):
                escaped = index > 0 and formula[index - 1] == "\\"
                if character == "{" and not escaped:
                    brace_depth += 1
                elif character == "}" and not escaped:
                    brace_depth -= 1
                    if brace_depth < 0:
                        break
            if brace_depth != 0:
                malformed.append(
                    f"{markdown.relative_to(ROOT)}: unbalanced braces in {formula.strip()[:60]}"
                )
    assert not incompatible_delimiters, "GitHub-incompatible math delimiters: " + ", ".join(incompatible_delimiters)
    assert not malformed, "Malformed fenced mathematics:\n" + "\n".join(malformed)
    assert not unsupported, (
        "Renderer-incompatible mathematical constructs:\n" + "\n".join(unsupported)
    )


def test_documentation_has_no_embedded_control_characters() -> None:
    failures: list[str] = []
    for markdown in [ROOT / "README.md", *(ROOT / "docs").rglob("*.md")]:
        text = markdown.read_text()
        controls = [character for character in text if ord(character) < 32 and character not in "\n\t"]
        if controls:
            failures.append(str(markdown.relative_to(ROOT)))
    assert not failures, "Documentation contains embedded control characters: " + ", ".join(failures)


def test_tns_documentation_uses_renderer_safe_latex() -> None:
    text = (ROOT / "docs" / "ALGORITHMS.md").read_text()
    assert r"q_n(j)=\min\left(j,n-1-j\right)." in text
    assert r"u^{TNS}_{m,L}(j)=\frac" in text
    assert r"\qquad(L\lt m)" in text
