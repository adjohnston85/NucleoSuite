"""Small end-to-end regression for the suite-to-tracks directory contract."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _write_tiny_inputs(directory: Path) -> tuple[Path, Path, Path]:
    pysam = pytest.importorskip("pysam")
    pytest.importorskip("pyBigWig")

    fasta = directory / "tiny.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 3000 + "\n", encoding="utf-8")
    pysam.faidx(str(fasta))

    fragments = directory / "tiny_fragments.bed"
    rows: list[str] = []
    lengths = (145, 147, 161, 167)
    for group, centre in enumerate(range(700, 10_000, 700)):
        for offset, length in enumerate(lengths):
            start = centre - length // 2 + offset - 2
            rows.append(f"chr1\t{start}\t{start + length}\tfragment_{group}_{length}\n")
    fragments.write_text("".join(rows), encoding="utf-8")

    ctcf = directory / "ctcf.bed"
    ctcf.write_text("chr1\t5000\t5001\tctcf\t0\t+\n", encoding="utf-8")
    return fasta, fragments, ctcf


@pytest.mark.parametrize("script_name", ["mnase_full_suite.sh", "cfdna_full_suite.sh"])
def test_suite_executes_combined_tracks_inside_declared_output_dir(
    tmp_path: Path,
    script_name: str,
) -> None:
    fasta, fragments, ctcf = _write_tiny_inputs(tmp_path)
    outdir = tmp_path / script_name.removesuffix("_full_suite.sh")
    executable = tmp_path / "nucleosuite-source"
    executable.write_text(
        f"#!{sys.executable}\n"
        "from nucleosuite.cli.main import main\n"
        "raise SystemExit(main())\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)

    command = [
        "bash",
        str(ROOT / "src" / "nucleosuite" / "resources" / script_name),
        "--fragments", str(fragments),
        "--fasta", str(fasta),
        "--ctcf-bed", str(ctcf),
        "--sample-name", "integration",
        "--outdir", str(outdir),
        "--contigs", "chr1",
        "--interval-format", "bed",
        "--no-blacklist",
        "--combine-prerequisites-only",
    ]
    environment = {
        **os.environ,
        "NUCLEOSUITE_BIN": str(executable),
        "PYTHON_BIN": sys.executable,
        "PYTHONPATH": str(ROOT / "src"),
        "MPLCONFIGDIR": str(tmp_path / "mplconfig"),
    }
    completed = subprocess.run(
        command,
        check=False,
        text=True,
        capture_output=True,
        env=environment,
        timeout=180,
    )

    stage_log = outdir / "logs" / "01_combined_tracks.log"
    report = outdir / "01_combined_tracks" / "completion_report.tsv"
    combined_output_dir = str((outdir / "01_combined_tracks").resolve())
    diagnostics = "\n".join(
        [completed.stdout, completed.stderr, stage_log.read_text() if stage_log.exists() else ""]
    )

    assert stage_log.is_file(), diagnostics
    assert report.is_file(), diagnostics
    assert "--output-dir" in stage_log.read_text()
    assert combined_output_dir in stage_log.read_text()
    assert "Parallel tracks requires every specification output_prefix" not in diagnostics
