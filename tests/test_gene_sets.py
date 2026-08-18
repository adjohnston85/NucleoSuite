"""Tests for rule-based gene-set generation."""

from __future__ import annotations

import csv
import os
from pathlib import Path

from nucleosuite.gene_sets import (
    build_gene_sets,
    filter_blacklisted_gene_anchors,
    load_rules,
    read_genes,
    read_states,
)
from nucleosuite.core.blacklist import load_blacklist_unbounded


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_default_style_rules_remove_candidate_overlap(tmp_path: Path):
    genes_path = write(
        tmp_path / "genes.bed",
        "chr1\t0\t100\tactive_only\n"
        "chr1\t100\t200\tshared\n"
        "chr1\t200\t300\trepressed_only\n"
        "chr1\t300\t400\tunassigned\n",
    )
    states_path = write(
        tmp_path / "states.bed",
        "chr1\t0\t200\t1_Active_Promoter\n"
        "chr1\t0\t200\t10_Txn_Elongation\n"
        "chr1\t100\t300\t12_Repressed\n",
    )
    config_path = write(
        tmp_path / "sets.tsv",
        "set_name\tinclude_rule\n"
        "active\t1_Active_Promoter & (9_Txn_Transition | 10_Txn_Elongation)\n"
        "repressed\t12_Repressed\n",
    )

    rules = load_rules(config_path, [])
    outputs = build_gene_sets(
        read_genes(genes_path),
        read_states(states_path),
        rules,
        tmp_path / "out",
        venn_sets=["active", "repressed"],
    )

    active_final = (tmp_path / "out/final_sets/active.bed").read_text().splitlines()
    repressed_final = (tmp_path / "out/final_sets/repressed.bed").read_text().splitlines()
    overlap = (tmp_path / "out/overlaps/shared_by_multiple_sets.bed").read_text().splitlines()
    assert any("active_only" in row for row in active_final)
    assert not any("shared" in row for row in active_final)
    assert any("repressed_only" in row for row in repressed_final)
    assert not any("shared" in row for row in repressed_final)
    assert len(overlap) == 1 and "shared" in overlap[0]
    assert outputs["venn"].is_file()

    with outputs["summary"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    summary = {row["set_name"]: row for row in rows}
    assert summary["active"]["candidate_gene_count"] == "2"
    assert summary["active"]["overlap_removed_count"] == "1"
    assert summary["active"]["final_gene_count"] == "1"
    assert summary["repressed"]["candidate_gene_count"] == "2"
    assert summary["repressed"]["final_gene_count"] == "1"


def test_exclusion_operator_is_rejected(tmp_path: Path):
    config_path = write(
        tmp_path / "sets.tsv",
        "set_name\tinclude_rule\nactive\t1_Active_Promoter & !12_Repressed\nrepressed\t12_Repressed\n",
    )
    try:
        load_rules(config_path, [])
    except ValueError as error:
        assert "inclusion only" in str(error)
    else:
        raise AssertionError("Expected exclusion rule to be rejected")


def test_directed_expression_categories_and_leftover_state_bed(tmp_path: Path):
    gene_names = [
        "active_only",
        "active_weak",
        "weak_only",
        "repressed_only",
        "active_repressed",
        "weak_repressed",
        "triple_overlap",
        "none",
    ]
    genes_path = tmp_path / "genes.bed"
    with genes_path.open("w") as handle:
        for index, name in enumerate(gene_names):
            start = index * 100
            handle.write(
                f"chr1\t{start}\t{start + 100}\tENSG{index:011d}\t{name}\t+\n"
            )

    state_map = {
        "active_only": ["1_Active_Promoter", "10_Txn_Elongation"],
        "active_weak": [
            "1_Active_Promoter",
            "2_Weak_Promoter",
            "10_Txn_Elongation",
        ],
        "weak_only": ["2_Weak_Promoter", "11_Weak_Txn"],
        "repressed_only": ["12_Repressed"],
        "active_repressed": [
            "1_Active_Promoter",
            "10_Txn_Elongation",
            "12_Repressed",
        ],
        "weak_repressed": ["2_Weak_Promoter", "11_Weak_Txn", "12_Repressed"],
        "triple_overlap": [
            "1_Active_Promoter",
            "2_Weak_Promoter",
            "10_Txn_Elongation",
            "12_Repressed",
        ],
        "none": [],
    }
    states_path = tmp_path / "states.bed"
    with states_path.open("w") as handle:
        for index, name in enumerate(gene_names):
            start = index * 100
            for state in state_map[name]:
                handle.write(f"chr1\t{start}\t{start + 100}\t{state}\n")

    config_path = write(
        tmp_path / "sets.tsv",
        "set_name\tinclude_rule\texclude_if_candidate\n"
        "active_genes\t1_Active_Promoter & (9_Txn_Transition | 10_Txn_Elongation)\trepressed_genes\n"
        "weak_genes\t2_Weak_Promoter & (9_Txn_Transition | 10_Txn_Elongation | 11_Weak_Txn)\tactive_genes,repressed_genes\n"
        "repressed_genes\t12_Repressed\tactive_genes,weak_genes\n",
    )

    outputs = build_gene_sets(
        read_genes(genes_path),
        read_states(states_path),
        load_rules(config_path, []),
        tmp_path / "out",
        venn_sets=["active_genes", "weak_genes", "repressed_genes"],
        leftover_set_name="leftover_genes",
    )

    id_by_name = {name: f"ENSG{index:011d}" for index, name in enumerate(gene_names)}

    def identifiers(path: Path) -> set[str]:
        rows = [line.split("\t") for line in path.read_text().splitlines() if line]
        assert all(len(row) == 6 for row in rows)
        assert all(row[4] == "0" and row[5] in {"+", "-", "."} for row in rows)
        return {row[3] for row in rows}

    assert identifiers(tmp_path / "out/final_sets/active_genes.bed") == {
        id_by_name["active_only"],
        id_by_name["active_weak"],
    }
    assert identifiers(tmp_path / "out/final_sets/weak_genes.bed") == {id_by_name["weak_only"]}
    assert identifiers(tmp_path / "out/final_sets/repressed_genes.bed") == {
        id_by_name["repressed_only"],
    }
    assert identifiers(tmp_path / "out/final_sets/leftover_genes.bed") == {
        id_by_name["none"],
    }

    state_rows = [line.split("\t") for line in outputs["final_state_interval"].read_text().splitlines()]
    assert all(len(row) == 6 for row in state_rows)
    assert all(row[4] == "0" and row[5] == "+" for row in state_rows)
    assert {row[3] for row in state_rows} == {
        "active_genes",
        "weak_genes",
        "repressed_genes",
        "leftover_genes",
    }

    with outputs["assignments"].open() as handle:
        assignments = {row["gene_name"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert assignments["active_weak"]["final_set"] == "active_genes"
    assert assignments["none"]["final_set"] == "leftover_genes"
    for name in ("active_repressed", "weak_repressed", "triple_overlap"):
        assert assignments[name]["candidate_sets"]
        assert assignments[name]["final_set"] == ""
    assert len(state_rows) == 5


def test_six_column_gene_input_converts_to_valid_bigbed(tmp_path: Path, monkeypatch):
    tools = tmp_path / "tools"
    tools.mkdir()
    converter = tools / "bedToBigBed"
    converter.write_text(
        "#!/usr/bin/env python3\n"
        "import shutil, sys\n"
        "shutil.copyfile(sys.argv[-3], sys.argv[-1])\n"
    )
    converter.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ.get('PATH', '')}")

    genes_path = write(
        tmp_path / "genes.bed",
        "chr1\t0\t100\tENSG000001\tGENE_A\t+\n"
        "chr1\t100\t200\tENSG000002\tGENE_B\t-\n",
    )
    states_path = write(
        tmp_path / "states.bed",
        "chr1\t0\t100\t1_Active_Promoter\n"
        "chr1\t0\t100\t10_Txn_Elongation\n"
        "chr1\t100\t200\t12_Repressed\n",
    )
    config_path = write(
        tmp_path / "sets.tsv",
        "set_name\tinclude_rule\texclude_if_candidate\n"
        "active_genes\t1_Active_Promoter & 10_Txn_Elongation\trepressed_genes\n"
        "repressed_genes\t12_Repressed\tactive_genes\n",
    )

    outputs = build_gene_sets(
        read_genes(genes_path),
        read_states(states_path),
        load_rules(config_path, []),
        tmp_path / "out",
        interval_format="both",
        chrom_sizes={"chr1": 1000},
        leftover_set_name="leftover_genes",
    )

    candidate_bed = tmp_path / "out/candidate_sets/active_genes.bed"
    candidate_bb = candidate_bed.with_suffix(".bb")
    assert candidate_bed.exists() and candidate_bb.exists()
    for row in candidate_bed.read_text().splitlines():
        fields = row.split("\t")
        assert len(fields) == 6
        assert fields[3].startswith("ENSG")
        assert fields[4] == "0"
        assert fields[5] in {"+", "-", "."}

    state_bed = outputs["final_state_interval"]
    assert state_bed.suffix == ".bed"
    assert state_bed.with_suffix(".bb").exists()
    assert all(len(row.split("\t")) == 6 for row in state_bed.read_text().splitlines())


def test_gene_and_state_contigs_resolve_against_chrom_sizes(tmp_path: Path):
    genes_path = write(
        tmp_path / "genes.bed",
        "chr20\t0\t100\tENSG000001\tGENE_A\t+\n",
    )
    states_path = write(
        tmp_path / "states.bed",
        "chr20\t0\t100\t12_Repressed\n",
    )
    genes = read_genes(genes_path, chrom_sizes={"20": 1000})
    states = read_states(states_path, chrom_sizes={"20": 1000})
    assert genes[0].chrom == "20"
    assert genes[0].fields[0] == "20"
    assert states[0].chrom == "20"


def test_gene_set_blacklist_filter_uses_tss_anchor_not_complete_gene(tmp_path: Path):
    genes_path = write(
        tmp_path / "genes.bed",
        "chr1\t10\t100\tPLUS\tPlus\t+\n"
        "chr1\t120\t200\tMINUS\tMinus\t-\n",
    )
    blacklist_path = write(
        tmp_path / "blacklist.bed",
        "chr1\t50\t60\nchr1\t199\t200\n",
    )
    retained, excluded = filter_blacklisted_gene_anchors(
        read_genes(genes_path), load_blacklist_unbounded(blacklist_path)
    )
    assert [gene.gene_id for gene in retained] == ["PLUS"]
    assert excluded == 1


def test_member_filename_prefix_is_opt_in_for_randomized_controls(tmp_path: Path):
    genes_path = write(
        tmp_path / "genes.bed",
        "chr1\t10\t100\tGENE\tGene\t+\n",
    )
    states_path = write(tmp_path / "states.bed", "chr1\t10\t100\tACTIVE\n")
    config_path = write(
        tmp_path / "sets.tsv",
        "set_name\tinclude_rule\nactive_genes\tACTIVE\nother_genes\tOTHER\n",
    )
    prefix = "sample_randomized_control_gene_sets"
    build_gene_sets(
        read_genes(genes_path),
        read_states(states_path),
        load_rules(config_path, []),
        tmp_path / "out",
        output_prefix=prefix,
        venn_sets=["active_genes", "other_genes"],
        prefix_member_files=True,
    )
    assert (tmp_path / "out/final_sets" / f"{prefix}_active_genes.bed").is_file()
    assert (tmp_path / "out/overlaps" / f"{prefix}_shared_by_multiple_sets.bed").is_file()


def test_final_gene_set_tss_outputs_use_strand_aware_coordinates(tmp_path: Path):
    genes_path = write(
        tmp_path / "genes.bed",
        "chr1\t100\t200\tGENE_PLUS\tPlus\t+\n"
        "chr1\t300\t450\tGENE_MINUS\tMinus\t-\n"
        "chr1\t500\t600\tGENE_OTHER\tOther\t+\n",
    )
    states_path = write(
        tmp_path / "states.bed",
        "chr1\t100\t200\tACTIVE\n"
        "chr1\t300\t450\tACTIVE\n"
        "chr1\t500\t600\tOTHER\n",
    )
    config_path = write(
        tmp_path / "sets.tsv",
        "set_name\tinclude_rule\n"
        "active_genes\tACTIVE\n"
        "other_genes\tOTHER\n",
    )
    outputs = build_gene_sets(
        read_genes(genes_path),
        read_states(states_path),
        load_rules(config_path, []),
        tmp_path / "out",
        venn_sets=["active_genes", "other_genes"],
    )

    tss_rows = [
        line.split("\t")
        for line in (tmp_path / "out/final_tss/active_genes.bed").read_text().splitlines()
    ]
    assert tss_rows == [
        ["chr1", "100", "101", "GENE_PLUS", "0", "+"],
        ["chr1", "449", "450", "GENE_MINUS", "0", "-"],
    ]
    labelled_rows = [
        line.split("\t") for line in outputs["final_tss_interval"].read_text().splitlines()
    ]
    assert labelled_rows == [
        ["chr1", "100", "101", "active_genes", "0", "+"],
        ["chr1", "449", "450", "active_genes", "0", "-"],
        ["chr1", "500", "501", "other_genes", "0", "+"],
    ]
    with outputs["summary"].open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["final_tss_interval"].endswith("final_tss/active_genes.bed")
