from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from nucleosuite.category_aggregate import split_region_bed_by_category


def test_category_split_is_one_based_and_deterministic(tmp_path: Path):
    bed = tmp_path / "tss.bed"
    bed.write_text(
        "chr1\t10\t11\tweak_genes\t0\t+\n"
        "chr1\t20\t21\tactive_genes\t0\t-\n"
        "chr1\t30\t31\tweak_genes\t0\t+\n"
    )
    categories, paths, counts = split_region_bed_by_category(
        bed, category_col=4, skip_header=False, destination=tmp_path / "split"
    )
    assert categories == ["active_genes", "weak_genes"]
    assert counts == {"active_genes": 1, "weak_genes": 2}
    assert paths["active_genes"].read_text().startswith("chr1\t20\t21")


def test_category_split_rejects_missing_column(tmp_path: Path):
    bed = tmp_path / "bad.bed"
    bed.write_text("chr1\t10\t11\n")
    with pytest.raises(ValueError, match="--category-col 4"):
        split_region_bed_by_category(
            bed, category_col=4, skip_header=False, destination=tmp_path / "split"
        )


def test_aggregate_help_exposes_category_mode():
    import argparse
    from nucleosuite.cli.aggregate import add_aggregate_parser

    root = argparse.ArgumentParser()
    parser = add_aggregate_parser(root.add_subparsers(dest="command", required=True))
    help_text = parser.format_help()
    assert "--category-col" in help_text
    assert "--category-profile-output" in help_text
    assert "--category-plot-output" in help_text
    assert "--category-nrl-summary-output" in help_text
    assert "State labels are not used to group" in help_text


def test_category_aggregate_end_to_end(tmp_path: Path):
    pyBigWig = pytest.importorskip("pyBigWig")
    from nucleosuite.cli.aggregate import main as aggregate_main

    chrom_size = 12000
    x = np.arange(chrom_size, dtype=float)
    values = 10.0 + 2.0 * np.cos(2.0 * np.pi * x / 185.0)
    bigwig = tmp_path / "synthetic_PNS.bw"
    bw = pyBigWig.open(str(bigwig), "w")
    bw.addHeader([("chr1", chrom_size)])
    bw.addEntries(
        ["chr1"] * chrom_size,
        list(range(chrom_size)),
        ends=list(range(1, chrom_size + 1)),
        values=values.tolist(),
    )
    bw.close()

    bed = tmp_path / "gene_sets_final_tss.bed"
    with bed.open("w") as out:
        for center, category, strand in [
            (3000, "active_genes", "+"),
            (4000, "repressed_genes", "-"),
            (5000, "active_genes", "+"),
            (6000, "repressed_genes", "-"),
            (7000, "active_genes", "+"),
            (8000, "repressed_genes", "-"),
        ]:
            out.write(f"chr1\t{center}\t{center + 1}\t{category}\t0\t{strand}\n")

    outdir = tmp_path / "out"
    assert aggregate_main(
        [
            "--bigwig", str(bigwig),
            "--region-bed", str(bed),
            "--category-col", "4",
            "--strand-col", "6",
            "--missing-strand", "error",
            "--window-half", "700",
            "--zero-thresh", "0",
            "--max-score", "inf",
            "--nrl-peak-resolution", "160",
            "--nrl-regression-max", "650",
            "--output-dir", str(outdir),
            "--output-prefix", "synthetic_TSS",
        ]
    ) == 0

    profiles = list(outdir.glob("*_category_profiles.tsv"))
    plots = list(outdir.glob("*_category_profiles.png"))
    summaries = list(outdir.glob("*_category_nrl_summary.tsv"))
    assert len(profiles) == len(plots) == len(summaries) == 1
    assert plots[0].stat().st_size > 0

    with profiles[0].open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == ["relative_position", "active_genes", "repressed_genes"]
        assert sum(1 for _ in reader) == 1401

    with summaries[0].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {(row["category"], row["direction"]) for row in rows} == {
        ("active_genes", "positive"),
        ("active_genes", "negative"),
        ("repressed_genes", "positive"),
        ("repressed_genes", "negative"),
    }
    assert {row["valid_region_count"] for row in rows} == {"3"}
    assert list((outdir / "active_genes").glob("*_aggregate_nrl_positive_regression.tsv"))
    assert list((outdir / "repressed_genes").glob("*_aggregate_nrl_negative_regression.tsv"))
