from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from nucleosuite import tss_expression_quintiles as teq


class FakeBigWig:
    def __init__(self, values: np.ndarray):
        self._values = np.asarray(values, dtype=float)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def chroms(self):
        return {"chr1": len(self._values)}

    def values(self, chrom, start, end, numpy=False):
        assert chrom == "chr1"
        result = self._values[start:end].copy()
        return result if numpy else result.tolist()


def test_tss_expression_quintiles_outputs_five_equal_groups(tmp_path, monkeypatch):
    genes = tmp_path / "genes.bed"
    gene_rows = []
    for index in range(10):
        start = 100 + index * 40
        strand = "+" if index % 2 == 0 else "-"
        gene_rows.append(
            f"chr1\t{start}\t{start + 20}\tENSG{index + 1:011d}\tG{index + 1}\t{strand}\n"
        )
    genes.write_text("".join(gene_rows))

    expression = tmp_path / "expression.tsv.gz"
    with gzip.open(expression, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Gene", "Gene name", "Tissue", "nTPM"])
        for index in range(10):
            writer.writerow([f"ENSG{index + 1:011d}", f"G{index + 1}", "bone marrow", index])

    signal_path = tmp_path / "signal.bw"
    signal_path.touch()
    values = np.arange(1000, dtype=float)
    monkeypatch.setitem(sys.modules, "pyBigWig", SimpleNamespace(open=lambda _: FakeBigWig(values)))

    prefix = tmp_path / "out" / "sample"
    code = teq.main(
        [
            "--signal", str(signal_path),
            "--sample", "S1",
            "--signal-label", "PNS",
            "--expression", str(expression),
            "--tissue", "bone_marrow",
            "--genes-bed", str(genes),
            "--window", "2",
            "--output-prefix", str(prefix),
        ]
    )
    assert code == 0

    summary = list(csv.DictReader(Path(f"{prefix}_tss_expression_quintile_summary.tsv").open(), delimiter="\t"))
    assert len(summary) == 5
    assert [int(row["assigned_gene_count"]) for row in summary] == [2, 2, 2, 2, 2]
    assert summary[0]["minimum_nTPM"] == "0.0"
    assert summary[-1]["maximum_nTPM"] == "9.0"

    profiles = list(csv.DictReader(Path(f"{prefix}_tss_expression_quintiles.tsv").open(), delimiter="\t"))
    assert len(profiles) == 5 * 5
    assert {row["quintile"] for row in profiles} == {
        "Q1_lowest", "Q2_20_40_percent", "Q3_middle", "Q4_60_80_percent", "Q5_highest"
    }
    assert Path(f"{prefix}_tss_expression_quintiles.png").stat().st_size > 0


def test_profile_selector_accepts_underscores(tmp_path):
    expression = tmp_path / "expression.tsv"
    expression.write_text(
        "Gene\tGene name\tTissue\tnTPM\n"
        "ENSG00000000001\tG1\tbone marrow\t2.0\n"
    )
    tissue, records, duplicates, invalid = teq.read_profile_expression(
        expression,
        profile="bone_marrow",
        gene_column="Gene",
        name_column="Gene name",
        profile_column="Tissue",
        value_column="nTPM",
    )
    assert tissue == "bone marrow"
    assert records["ENSG00000000001"].value == 2.0
    assert duplicates == 0
    assert invalid == 0
