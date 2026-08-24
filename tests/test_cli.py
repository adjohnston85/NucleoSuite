"""Smoke tests for the NucleoSuite command surface."""

from __future__ import annotations

from pathlib import Path

import argparse
import re

import pytest

from nucleosuite import __version__
from nucleosuite.cli.main import DELEGATED_COMMANDS, build_parser


def subparser_choices(parser: argparse.ArgumentParser) -> set[str]:
    action = next(
        item for item in parser._actions
        if isinstance(item, argparse._SubParsersAction)
    )
    return set(action.choices)


def test_version_is_current_release():
    assert __version__ == "0.10.6"


def test_all_primary_commands_are_registered():
    choices = subparser_choices(build_parser())
    assert {
        "tracks", "pns", "wps", "coverage", "dyads", "dyad", "fragment-ends", "mean-scale",
        "dinuc-profile", "ww-types", "call-peaks", "peak-call",
        "aggregate", "compare-positions", "dac", "dcc", "distances", "fragment-lengths",
        "fragment-heatmap", "filter-peaks", "gene-sets", "gene-expression", "tss-expression-quintiles", "mnase-suite", "cfdna-suite", "chip-suite", "chip-compare", "chrom-sizes", "nrl", "plot", "positive-runs", "peak-score-frequency", "pns-peak-fdr", "peak-states", "resources", "region-extract", "validate-inputs",
    } <= choices


def _option_help_blocks(help_text: str) -> list[str]:
    """Return argparse option/positional blocks including wrapped lines."""
    starts = [match.start() for match in re.finditer(r"(?m)^  \S", help_text)]
    return [
        help_text[start : starts[index + 1] if index + 1 < len(starts) else len(help_text)]
        for index, start in enumerate(starts)
    ]


def test_every_command_help_shows_each_default_at_most_once(capsys):
    from nucleosuite.cli import main as cli_main

    commands = sorted(subparser_choices(build_parser()))
    invocations = [[command, "--help"] for command in commands]
    invocations.extend(
        [
            ["dcc", "bigwig", "--help"],
            ["dcc", "bam", "--help"],
            ["region-peak-extractor", "--help"],
            ["pns-region-extractor", "--help"],
        ]
    )
    for invocation in invocations:
        try:
            exit_code = cli_main.main(invocation)
        except SystemExit as error:
            exit_code = error.code
        assert exit_code == 0
        captured = capsys.readouterr()
        for block in _option_help_blocks(captured.out):
            annotations = re.findall(r"\(default:[^)]*\)", block, flags=re.IGNORECASE)
            assert len(annotations) <= 1, (
                f"duplicated defaults in {' '.join(invocation)}:\n{block}"
            )



def test_help_all_is_available_for_every_primary_command(capfd):
    from nucleosuite.cli import main as cli_main

    for command in sorted(subparser_choices(build_parser())):
        assert cli_main.main([command, "--help-all"]) == 0
        capfd.readouterr()


def test_core_help_hides_curated_advanced_options(capsys):
    from nucleosuite.cli import main as cli_main

    examples = (
        ("distances", "--pct-bin-seed"),
        ("pns", "--chunk-bp"),
        ("compare-positions", "--stats-test"),
        ("plot", "--mpl-rc"),
    )
    for command, advanced_option in examples:
        assert cli_main.main([command, "--help"]) == 0
        core = capsys.readouterr().out
        assert "--help-all" in core
        assert advanced_option not in core

        assert cli_main.main([command, "--help-all"]) == 0
        extended = capsys.readouterr().out
        assert advanced_option in extended
        assert len(extended) > len(core)


def test_analysis_help_is_layered_and_plotting_stays_separate(capsys):
    from nucleosuite.cli import main as cli_main

    for command in ("tracks", "distances", "nrl", "dac", "positive-runs", "peak-states"):
        assert cli_main.main([command, "--help"]) == 0
        normal = capsys.readouterr().out
        assert "--help-all" in normal
        assert "--help-plotting" in normal
        assert "--plot-format" not in normal

        assert cli_main.main([command, "--help-all"]) == 0
        extended = capsys.readouterr().out
        assert len(extended) > len(normal)
        assert "Full command-specific help shown" in extended
        assert "--plot-format" not in extended

        with pytest.raises(SystemExit) as error:
            cli_main.main([command, "--help-plotting"])
        assert error.value.code == 0
        plotting = capsys.readouterr().out
        assert "--plot-format" in plotting
        assert "--plot-point-label-offset" in plotting


def test_suite_help_hides_plot_block_and_expands_on_request(capfd):
    from nucleosuite.cli import main as cli_main

    for command in ("mnase-suite", "cfdna-suite"):
        assert cli_main.main([command, "--help"]) == 0
        normal = capfd.readouterr().out
        assert "--help-plotting" in normal
        assert "--plot-format" not in normal

        assert "--help-all" in normal

        assert cli_main.main([command, "--help-all"]) == 0
        full = capfd.readouterr().out
        assert len(full) > len(normal)
        assert "--analysis-cores" in full

        assert cli_main.main([command, "--help-plotting"]) == 0
        expanded = capfd.readouterr().out
        assert "--plot-format" in expanded
        assert "--plot-point-label-offset" in expanded


def test_plot_command_uses_core_and_extended_help(capsys):
    from nucleosuite.cli import main as cli_main

    assert cli_main.main(["plot", "--help"]) == 0
    core = capsys.readouterr().out
    assert "--help-all" in core
    assert "--x-min" in core
    assert "--mpl-rc" not in core

    assert cli_main.main(["plot", "--help-all"]) == 0
    extended = capsys.readouterr().out
    assert "--x-major-tick" in extended
    assert "--mpl-rc" in extended
    assert len(extended) > len(core)

def test_region_extractor_aliases_are_delegated():
    assert DELEGATED_COMMANDS["region-extract"][0] is DELEGATED_COMMANDS["region-peak-extractor"][0]
    assert DELEGATED_COMMANDS["region-extract"][0] is DELEGATED_COMMANDS["pns-region-extractor"][0]


def test_peak_commands_expose_only_bed_output():
    parser = build_parser()
    action = next(
        item for item in parser._actions
        if isinstance(item, argparse._SubParsersAction)
    )
    for command in ("pns", "wps"):
        help_text = action.choices[command].format_help()
        assert "--peak-format" not in help_text
        assert "rich" not in help_text.lower()


def test_orc_command_console_messages_are_paired(monkeypatch, capsys, tmp_path):
    from nucleosuite.cli import main as cli_main
    from nucleosuite.console import CommandMessagePair

    monkeypatch.setattr(
        cli_main,
        "message_pair",
        lambda: CommandMessagePair("Zug zug.", "Work complete."),
    )
    monkeypatch.setitem(
        cli_main.DELEGATED_COMMANDS,
        "test-success",
        (lambda argv=None: 0, "Test successful command."),
    )
    monkeypatch.chdir(tmp_path)

    assert cli_main.main(["test-success"]) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].startswith("[test-success] Command log: ")
    assert lines[1:] == ["Zug zug.", "Work complete."]
    logs = list((tmp_path / "nucleosuite_logs" / "logs" / "commands").glob("*.log"))
    assert len(logs) == 1
    log_text = logs[0].read_text()
    assert "[COMMAND] nucleosuite test-success" in log_text
    assert "exit_code=0" in log_text


def test_peasant_command_console_messages_are_paired(monkeypatch, capsys, tmp_path):
    from nucleosuite.cli import main as cli_main
    from nucleosuite.console import CommandMessagePair

    monkeypatch.setattr(
        cli_main,
        "message_pair",
        lambda: CommandMessagePair("Yes, milord.", "Job's done."),
    )
    monkeypatch.setitem(
        cli_main.DELEGATED_COMMANDS,
        "test-success",
        (lambda argv=None: 0, "Test successful command."),
    )
    monkeypatch.chdir(tmp_path)

    assert cli_main.main(["test-success"]) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].startswith("[test-success] Command log: ")
    assert lines[1:] == ["Yes, milord.", "Job's done."]


def test_failed_command_does_not_print_completion(monkeypatch, capsys, tmp_path):
    from nucleosuite.cli import main as cli_main
    from nucleosuite.console import CommandMessagePair

    monkeypatch.setattr(
        cli_main,
        "message_pair",
        lambda: CommandMessagePair("Okay.", "Job's done."),
    )
    monkeypatch.setitem(
        cli_main.DELEGATED_COMMANDS,
        "test-failure",
        (lambda argv=None: 2, "Test failed command."),
    )
    monkeypatch.chdir(tmp_path)

    assert cli_main.main(["test-failure"]) == 2
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].startswith("[test-failure] Command log: ")
    assert lines[1:] == ["Okay."]
    log = next((tmp_path / "nucleosuite_logs" / "logs" / "commands").glob("*.log"))
    assert "exit_code=2" in log.read_text()


def test_executing_commands_have_generic_progress_bookends(monkeypatch, capsys, tmp_path):
    from nucleosuite.cli import main as cli_main
    from nucleosuite.console import CommandMessagePair

    monkeypatch.setattr(
        cli_main,
        "message_pair",
        lambda: CommandMessagePair("Okay.", "Job's done."),
    )
    monkeypatch.setitem(
        cli_main.DELEGATED_COMMANDS,
        "test-progress",
        (lambda argv=None: 0, "Test command progress."),
    )
    monkeypatch.chdir(tmp_path)

    assert cli_main.main(["test-progress"]) == 0
    captured = capsys.readouterr()
    assert "[test-progress] Starting analysis" in captured.err
    assert "[test-progress] Analysis completed; elapsed" in captured.err


def test_resource_commands_are_machine_readable_without_console_messages(capsys):
    from nucleosuite.cli import main as cli_main

    assert cli_main.main(["resources", "show", "genes", "--resource-set", "hg19-gm12878"]) == 0
    output = capsys.readouterr().out.strip().splitlines()
    assert len(output) == 1
    assert output[0].endswith("hg19_ensembl_genes.bed")
    assert output[0] not in {
        "Zug zug.", "Dabu.", "Swobu.", "Lok'tar.",
        "Okay.", "Right-o.", "Alright.", "Yes, milord.",
        "Work complete.", "Job's done.",
    }


def test_all_console_message_pairs_use_the_matching_completion_family():
    from nucleosuite.console import (
        MESSAGE_PAIRS,
        ORC_COMPLETION_MESSAGE,
        ORC_STARTUP_MESSAGES,
        PEASANT_COMPLETION_MESSAGE,
        PEASANT_STARTUP_MESSAGES,
    )

    pairs = {pair.startup: pair.completion for pair in MESSAGE_PAIRS}
    assert set(pairs) == set(ORC_STARTUP_MESSAGES) | set(PEASANT_STARTUP_MESSAGES)
    assert all(pairs[startup] == ORC_COMPLETION_MESSAGE for startup in ORC_STARTUP_MESSAGES)
    assert all(pairs[startup] == PEASANT_COMPLETION_MESSAGE for startup in PEASANT_STARTUP_MESSAGES)


def test_logo_is_the_embedded_fixed_artwork():
    from nucleosuite.logo import LOGO, LOGO_LINES, render_logo

    assert render_logo() == LOGO
    assert LOGO.splitlines() == list(LOGO_LINES)
    assert len(LOGO_LINES) == 36
    assert LOGO_LINES[0] == "                       --------O......O------"
    assert LOGO_LINES[-1] == "                       -------o....o---------"
    assert "▄   ▄ ▄   ▄" in LOGO


def test_top_level_help_is_prefixed_by_logo():
    from nucleosuite.logo import render_logo

    help_text = build_parser().format_help()
    logo = render_logo()
    assert help_text.startswith(logo + "\n\nusage: nucleosuite")
    assert "Turn paired-end sequencing fragments into nucleosome signal tracks" in help_text


def test_registered_subcommand_help_does_not_include_logo():
    from nucleosuite.logo import render_logo

    parser = build_parser()
    action = next(
        item for item in parser._actions
        if isinstance(item, argparse._SubParsersAction)
    )
    logo_first_line = render_logo().splitlines()[0]
    assert logo_first_line not in action.choices["pns"].format_help()
    assert logo_first_line not in action.choices["wps"].format_help()


def test_top_level_help_prints_without_warcraft_messages(capsys):
    from nucleosuite.cli import main as cli_main

    try:
        cli_main.main(["--help"])
    except SystemExit as error:
        assert error.code == 0

    from nucleosuite.logo import render_logo

    output = capsys.readouterr().out
    assert output.startswith(render_logo() + "\n\nusage: nucleosuite")
    assert "usage: nucleosuite" in output
    for message in (
        "Zug zug.", "Dabu.", "Swobu.", "Lok'tar.",
        "Okay.", "Right-o.", "Alright.", "Yes, milord.",
        "Work complete.", "Job's done.",
    ):
        assert message not in output


def test_top_level_version_is_prefixed_by_logo(capsys):
    from nucleosuite.cli import main as cli_main
    from nucleosuite.logo import render_logo

    try:
        cli_main.main(["--version"])
    except SystemExit as error:
        assert error.code == 0

    output = capsys.readouterr().out
    assert output == f"{render_logo()}\n\nnucleosuite {__version__}\n"
    for message in (
        "Zug zug.", "Dabu.", "Swobu.", "Lok'tar.",
        "Okay.", "Right-o.", "Alright.", "Yes, milord.",
        "Work complete.", "Job's done.",
    ):
        assert message not in output


def test_module_entrypoint_help_and_version_include_logo():
    """Exercise the same module entry point used by the console script."""
    import os
    import subprocess
    import sys

    from nucleosuite import __version__
    from nucleosuite.logo import render_logo

    env = os.environ.copy()
    source_root = str(Path(__file__).resolve().parents[1] / "src")
    env["PYTHONPATH"] = source_root + os.pathsep + env.get("PYTHONPATH", "")

    for option, suffix in (
        ("--help", "usage: nucleosuite"),
        ("--version", f"nucleosuite {__version__}"),
    ):
        result = subprocess.run(
            [sys.executable, "-m", "nucleosuite", option],
            check=False,
            capture_output=True,
            text=True,
            env=env,
        )
        assert result.returncode == 0
        assert result.stdout.startswith(render_logo() + "\n\n")
        assert suffix in result.stdout


def test_delegated_help_does_not_print_job_messages(capsys):
    from nucleosuite.cli import main as cli_main

    assert cli_main.main(["combine", "-h"]) == 0
    output = capsys.readouterr().out
    assert "usage: nucleosuite combine" in output
    for message in (
        "Zug zug.", "Dabu.", "Swobu.", "Lok'tar.",
        "Okay.", "Right-o.", "Alright.", "Yes, milord.",
        "Work complete.", "Job's done.",
    ):
        assert message not in output


def test_unknown_command_does_not_print_job_messages(capsys):
    from nucleosuite.cli import main as cli_main

    with pytest.raises(SystemExit) as error:
        cli_main.main(["not-a-command"])
    assert error.value.code == 2
    captured = capsys.readouterr()
    assert "Zug zug." not in captured.out
    assert "Right-o." not in captured.out


def test_suite_argument_error_does_not_print_job_messages(capsys):
    from nucleosuite.cli import main as cli_main

    with pytest.raises(SystemExit):
        cli_main.main(["cfdna-suite"])
    captured = capsys.readouterr()
    combined = captured.out + captured.err
    assert "Zug zug." not in combined
    assert "Lok'tar." not in combined
    assert "Right-o." not in combined
    assert "Alright." not in combined
