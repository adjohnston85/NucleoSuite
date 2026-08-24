from pathlib import Path

import pytest

from nucleosuite.cli.suite_paired import (
    annotate_suite_combined_peaks,
    extract_paired_options,
)


def _write_peaks(path: Path, scores):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(
            f"chr1\t{i * 100}\t{i * 100 + 80}\tp{i}\t{score}\t.\t"
            f"{i * 100 + 40}\t{i * 100 + 41}\n"
            for i, score in enumerate(scores)
        ),
        encoding="utf-8",
    )


def test_extract_paired_options_consumes_wrapper_arguments():
    paired, fdr, remaining = extract_paired_options(
        ["--bam", "sample.bam", "--with-randomized-control", "--fdr", "0.05"]
    )
    assert paired is True
    assert fdr == 0.05
    assert remaining == ["--bam", "sample.bam"]


def test_suite_fdr_requires_paired_execution():
    with pytest.raises(ValueError, match="requires --with-randomized-control"):
        extract_paired_options(["--fdr", "0.05"])


def test_paired_suite_annotates_combined_nucleosome_and_breakpoint_beds(tmp_path: Path):
    scaled = tmp_path / "combined" / "01_combined_tracks" / "scaled"
    for suffix in ("nucleosome_regions", "breakpoint_peaks"):
        _write_peaks(
            scaled / f"sample_PNS_{suffix}_mean_scaled.bed", [30, 20, 10]
        )
        _write_peaks(
            scaled / f"sample_randomized_control_PNS_{suffix}_mean_scaled.bed", [15, 5]
        )

    outputs = annotate_suite_combined_peaks(
        tmp_path, suite_name="cfdna-suite", fdr_threshold=None
    )

    assert set(outputs) == {"nucleosome", "breakpoint"}
    for result in outputs.values():
        rows = result.annotated_path.read_text().splitlines()
        assert len(rows) == 3
        assert all(len(row.split("\t")) == 9 for row in rows)
