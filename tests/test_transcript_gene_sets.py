"""Transcript-level promoter logic and Ensembl annotation conversion."""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from nucleosuite.gene_sets import build_gene_sets, load_rules, read_genes, read_states
from nucleosuite.transcript_tss import extract_transcript_tss, read_transcript_tss
from nucleosuite.resource_files import materialized_resource_path


def _write(path, data):
    path.write_text(data, encoding="utf-8")
    return path


def test_ensembl_transcript_gtf_conversion_and_validation(tmp_path):
    genes = read_genes(_write(tmp_path / "genes.bed",
        "chr1\t100\t300\tENSG000001\tA\t+\n"
        "chr1\t400\t800\tENSG000002\tB\t-\n"))
    gtf = _write(tmp_path / "annotation.gtf",
        '1\tEnsembl\tgene\t101\t300\t.\t+\t.\tgene_id "ENSG000001";\n'
        '1\tEnsembl\ttranscript\t101\t300\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000001";\n'
        '1\tEnsembl\ttranscript\t201\t300\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000002";\n'
        '1\tEnsembl\ttranscript\t401\t750\t.\t-\t.\tgene_id "ENSG000002"; transcript_id "ENST000003";\n'
        '1\tEnsembl\ttranscript\t401\t800\t.\t-\t.\tgene_id "ENSG000002"; transcript_id "ENST000004";\n'
        '1\tEnsembl\ttranscript\t401\t800\t.\t+\t.\tgene_id "ENSG000002"; transcript_id "WRONG_STRAND";\n')
    result = tmp_path / "tss.tsv.gz"
    assert extract_transcript_tss(gtf, genes, result) == (4, 2)
    mapped = read_transcript_tss(result, genes)
    assert [(x.transcript_id, x.start) for x in mapped["ENSG000001"]] == [
        ("ENST000001", 100), ("ENST000002", 200)]
    assert [(x.transcript_id, x.start) for x in mapped["ENSG000002"]] == [
        ("ENST000003", 749), ("ENST000004", 799)]
    bad = _write(tmp_path / "bad.tsv", "\t".join(("gene_id", "transcript_id", "chrom", "tss_start", "tss_end", "strand"))
        + "\nENSG000001\tENSTX\tchr1\t301\t302\t+\n")
    with pytest.raises(ValueError, match="does not match"):
        read_transcript_tss(bad, genes)


def test_all_transcript_promoters_and_trimmed_body_plus_minus(tmp_path):
    genes_path = _write(tmp_path / "genes.bed",
        "chr1\t100\t500\tACTIVE_LATE\tActiveLate\t+\n"
        "chr1\t600\t1000\tACTIVE_MINUS\tActiveMinus\t-\n"
        "chr1\t1100\t1500\tWEAK_ONLY\tWeakOnly\t+\n"
        "chr1\t1600\t2000\tREP_MULTI\tRepMulti\t+\n"
        "chr1\t2100\t2500\tREP_PURE\tRepPure\t+\n"
        "chr1\t2600\t3000\tACTIVE_WEAK\tActiveWeak\t+\n"
        "chr1\t3100\t3500\tEARLY_BAD\tEarlyBad\t+\n")
    states_path = _write(tmp_path / "states.bed",
        # Earlier nonqualifying promoter + repressed state must fall outside trimmed gene.
        "chr1\t100\t110\t12_Repressed\n"
        "chr1\t200\t201\t1_Active_Promoter\n"
        "chr1\t300\t301\t1_Active_Promoter\n"
        "chr1\t350\t360\t10_Txn_Elongation\n"
        # Minus: choose highest transcript TSS, trimming right boundary.
        "chr1\t920\t930\t12_Repressed\n"
        "chr1\t699\t700\t1_Active_Promoter\n"
        "chr1\t799\t800\t1_Active_Promoter\n"
        "chr1\t650\t660\t10_Txn_Elongation\n"
        # Weak TSS is internal.
        "chr1\t1200\t1201\t2_Weak_Promoter\n"
        "chr1\t1300\t1320\t11_Weak_Txn\n"
        # Repressed has a promoter at the SECOND transcript TSS -> exclude.
        "chr1\t1700\t1710\t12_Repressed\n"
        "chr1\t1800\t1801\t1_Active_Promoter\n"
        "chr1\t2200\t2210\t12_Repressed\n"
        # Both types of promoter for one gene, active takes precedence.
        "chr1\t2600\t2601\t1_Active_Promoter\n"
        "chr1\t2700\t2701\t2_Weak_Promoter\n"
        "chr1\t2800\t2810\t10_Txn_Elongation\n"
        # First active promoter is in a repressed interval; later promoter is clean.
        "chr1\t3100\t3150\t12_Repressed\n"
        "chr1\t3110\t3111\t1_Active_Promoter\n"
        "chr1\t3200\t3201\t1_Active_Promoter\n"
        "chr1\t3350\t3360\t10_Txn_Elongation\n")
    genes = read_genes(genes_path)
    rows = [(g, t, pos, strand) for g, t, pos, strand in [
        ("ACTIVE_LATE", "tx1", 100, "+"), ("ACTIVE_LATE", "tx2", 200, "+"),
        ("ACTIVE_LATE", "tx3", 300, "+"),
        ("ACTIVE_MINUS", "tx4", 929, "-"), ("ACTIVE_MINUS", "tx5", 799, "-"),
        ("ACTIVE_MINUS", "tx6", 699, "-"),
        ("WEAK_ONLY", "tx7", 1100, "+"), ("WEAK_ONLY", "tx8", 1200, "+"),
        ("REP_MULTI", "tx9", 1600, "+"), ("REP_MULTI", "tx10", 1800, "+"),
        ("REP_PURE", "tx11", 2100, "+"), ("REP_PURE", "tx12", 2200, "+"),
        ("ACTIVE_WEAK", "tx13", 2600, "+"), ("ACTIVE_WEAK", "tx14", 2700, "+"),
        ("EARLY_BAD", "tx15", 3110, "+"), ("EARLY_BAD", "tx16", 3200, "+"),
    ]]
    path = tmp_path / "tss.tsv.gz"
    with gzip.open(path, "wt") as stream:
        stream.write("gene_id\ttranscript_id\tchrom\ttss_start\ttss_end\tstrand\n")
        for gid, tid, start, strand in rows:
            stream.write(f"{gid}\t{tid}\tchr1\t{start}\t{start+1}\t{strand}\n")
    with materialized_resource_path("default-gene-sets") as cfg:
        rules = load_rules(cfg, [])
    outputs = build_gene_sets(genes, read_states(states_path), rules, tmp_path / "result",
                               leftover_set_name="leftover_genes",
                               transcript_tss=read_transcript_tss(path, genes))
    def intervals(name):
        return {row[3]: (int(row[1]), int(row[2])) for row in (
            line.split("\t") for line in (tmp_path / "result/final_sets" / f"{name}.bed").read_text().splitlines())}
    assert intervals("active_genes") == {
        "ACTIVE_LATE": (200, 500),
        "ACTIVE_MINUS": (600, 800),
        "ACTIVE_WEAK": (2600, 3000),
        "EARLY_BAD": (3200, 3500),
    }
    assert intervals("weak_genes") == {"WEAK_ONLY": (1200, 1500)}
    assert intervals("repressed_genes") == {"REP_PURE": (2100, 2500)}
    assert "REP_MULTI" in intervals("leftover_genes")
    with outputs["assignments"].open() as handle:
        assignments = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert assignments["REP_MULTI"]["final_set"] == "leftover_genes"
    assert assignments["REP_MULTI"]["tss_intersecting_states"].find("1_Active_Promoter") >= 0
    assert assignments["ACTIVE_MINUS"]["selected_transcript_id"] == "tx5"
    assert assignments["ACTIVE_LATE"]["original_start"] == "100"
    assert assignments["ACTIVE_LATE"]["start"] == "200"
    assert "12_Repressed" not in assignments["ACTIVE_LATE"]["intersecting_states"]
    assert "12_Repressed" in assignments["ACTIVE_LATE"]["original_intersecting_states"]
    assert assignments["ACTIVE_LATE"]["transcript_tss_count"] == "3"
    selected_rows = list(csv.DictReader(outputs["selected_transcript_tss"].open(), delimiter="\t"))
    assert len({row["gene_id"] for row in selected_rows}) == len(selected_rows)
    assert next(row for row in selected_rows if row["gene_id"] == "ACTIVE_LATE")["transcript_id"] == "tx2"
    final_tss = (tmp_path / "result/final_tss/active_genes.bed").read_text()
    assert "chr1\t799\t800\tACTIVE_MINUS\t0\t-" in final_tss


def test_command_accepts_gtf_and_emits_auditable_outputs(tmp_path):
    from nucleosuite.gene_sets import main

    genes = _write(tmp_path / "genes.bed", "chr1\t100\t400\tENSG001\tA\t+\n")
    states = _write(tmp_path / "states.bed",
                    "chr1\t200\t201\t1_Active_Promoter\n"
                    "chr1\t300\t320\t10_Txn_Elongation\n"
                    "chr1\t100\t110\t12_Repressed\n")
    gtf = _write(tmp_path / "annotation.gtf",
                 '1\tEnsembl\ttranscript\t101\t400\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001";\n'
                 '1\tEnsembl\ttranscript\t201\t400\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST002";\n')
    config = _write(tmp_path / "rules.tsv", "set_name\tinclude_rule\trequired_tss_state\tforbidden_tss_states\tforbidden_gene_states\n"
                    "active\t10_Txn_Elongation\t1_Active_Promoter\t\t12_Repressed\n"
                    "repressed\t12_Repressed\t\t1_Active_Promoter\t\n")
    assert main(["--genes-bed", str(genes), "--states-bed", str(states),
                 "--config", str(config), "--transcript-gtf", str(gtf),
                 "--output-dir", str(tmp_path / "out")]) == 0
    assert (tmp_path / "out/final_sets/active.bed").read_text().strip() == "chr1\t200\t400\tENSG001\t0\t+"
    selected = tmp_path / "out/genes_states_gene_sets_selected_transcript_tss.tsv"
    # The default prefix derives from the two input basenames.
    selected = next((tmp_path / "out").glob("*_selected_transcript_tss.tsv"))
    assert "ENST002" in selected.read_text()
    assert next((tmp_path / "out").glob("*_ensembl87_transcript_tss.tsv.gz")).stat().st_size > 50


def test_default_config_reports_download_failure_clearly(monkeypatch, tmp_path):
    from nucleosuite import transcript_tss as module
    genes = read_genes(_write(tmp_path / "genes.bed", "chr1\t100\t300\tENSG001\tA\t+\n"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    def fail(*args, **kwargs):
        raise OSError("outbound network unavailable")
    monkeypatch.setattr(module, "urlopen", fail)
    with pytest.raises(RuntimeError, match="Could not obtain the Ensembl"):
        module.ensure_ensembl_tss(genes)
    assert not module.ensembl_cache_path().exists()


def test_gene_missing_transcripts_is_unclassified_not_leftover(tmp_path):
    genes = read_genes(_write(tmp_path / "genes.bed",
        "chr1\t100\t400\tMATCHED\tA\t+\n"
        "chr1\t500\t800\tNOT_ANNOTATED\tB\t+\n"))
    states = read_states(_write(tmp_path / "states.bed",
        "chr1\t100\t400\t12_Repressed\n"
        "chr1\t500\t800\t12_Repressed\n"
        "chr1\t350\t351\t1_Active_Promoter\n"))
    annotation = _write(tmp_path / "transcripts.tsv",
        "gene_id\ttranscript_id\tchrom\ttss_start\ttss_end\tstrand\n"
        "MATCHED\tTX1\tchr1\t100\t101\t+\n")
    cfg = _write(tmp_path / "rules.tsv",
        "set_name\tinclude_rule\tforbidden_tss_states\texclude_if_candidate\n"
        "repressed\t12_Repressed\t1_Active_Promoter\t\n"
        "other\t1_Active_Promoter\t\trepressed\n")
    outputs = build_gene_sets(genes, states, load_rules(cfg, []), tmp_path / "out",
                              transcript_tss=read_transcript_tss(annotation, genes),
                              leftover_set_name="leftover_genes")
    assert (tmp_path / "out/final_sets/repressed.bed").read_text().count("MATCHED") == 1
    assert (tmp_path / "out/final_sets/leftover_genes.bed").read_text() == ""
    with outputs["assignments"].open() as stream:
        rows = {r["gene_id"]: r for r in csv.DictReader(stream, delimiter="\t")}
    assert rows["NOT_ANNOTATED"]["final_set"] == ""
    assert rows["NOT_ANNOTATED"]["tss_annotation_status"] == "missing"


def test_bundled_ensembl_annotations_are_present_and_gene_matched():
    """Release-87 transcript TSSs must be included in the installed resources."""
    from nucleosuite.resource_files import validate_manifest
    import gzip

    with materialized_resource_path("hg19-ensembl87-gtf") as gtf:
        assert gtf.stat().st_size > 40_000_000
        with gzip.open(gtf, "rt") as handle:
            assert "GRCh37" in handle.readline()
    with materialized_resource_path("hg19-genes") as bed:
        genes = read_genes(bed)
    with materialized_resource_path("hg19-ensembl87-transcript-tss") as tss:
        mapped = read_transcript_tss(tss, genes)
        assert len(mapped) == 19_383
        assert sum(len(rows) for rows in mapped.values()) == 143_066
        assert "ENSG00000186092" in mapped
    assert all(ok for _, _, ok, _ in validate_manifest())


def test_default_gene_sets_use_bundled_tss_without_network(monkeypatch, tmp_path):
    """The standard command loads shipped transcript TSS positions offline."""
    from nucleosuite.gene_sets import main
    from nucleosuite import transcript_tss

    def no_download(*_args, **_kwargs):
        raise AssertionError("The bundled annotation should be used without a download")

    monkeypatch.setattr(transcript_tss, "urlopen", no_download)
    monkeypatch.setattr(transcript_tss, "ensure_ensembl_tss", no_download)
    genes = _write(tmp_path / "genes.bed", "chr1\t69090\t70008\tENSG00000186092\tOR4F5\t+\n")
    states = _write(tmp_path / "states.bed",
                    "chr1\t69090\t69091\t1_Active_Promoter\n"
                    "chr1\t69100\t69200\t10_Txn_Elongation\n")
    with materialized_resource_path("default-gene-sets") as rules:
        assert main(["--genes-bed", str(genes), "--states-bed", str(states),
                     "--config", str(rules), "--output-dir", str(tmp_path / "result")]) == 0
    final = (tmp_path / "result/final_sets/active_genes.bed").read_text()
    assert "ENSG00000186092" in final
    selected = next((tmp_path / "result").glob("*_selected_transcript_tss.tsv"))
    assert "ENST00000335137" in selected.read_text()
