"""Tests for pooled and individual gene-region DAC grouping."""

from pathlib import Path

from nucleosuite.dac import (
    expand_gene_queries,
    group_regions_by_state,
    read_selected_gene_regions,
)


def test_selected_genes_support_ensembl_and_gene_names(tmp_path: Path):
    genes = tmp_path / "genes.bed"
    genes.write_text(
        "chr1\t0\t100\tENSG00000141510.7\tTP53\t+\n"
        "chr1\t200\t350\tENSG00000146648\tEGFR\t-\n"
        "chr2\t100\t250\tENSG00000136997\tMYC\t+\n"
    )
    gene_list = tmp_path / "genes.txt"
    gene_list.write_text("EGFR\n")
    queries = expand_gene_queries(["ENSG00000141510,MYC"], [str(gene_list)])

    regions, rows = read_selected_gene_regions(
        path=str(genes),
        queries=queries,
        gene_id_column=4,
        gene_name_column=5,
        mode="both",
        pool_name="selected_genes",
        label_mode="name-id",
        selected_chromosomes=None,
        min_region_length=2,
    )
    grouped = group_regions_by_state(regions)
    assert len(rows) == 3
    assert len(grouped["selected_genes"]) == 3
    assert "TP53__ENSG00000141510.7" in grouped
    assert "EGFR__ENSG00000146648" in grouped
    assert "MYC__ENSG00000136997" in grouped


def test_selected_genes_pooled_mode_creates_one_group(tmp_path: Path):
    genes = tmp_path / "genes.bed"
    genes.write_text(
        "chr1\t0\t100\tENSG1\tGENE1\t+\n"
        "chr1\t200\t300\tENSG2\tGENE2\t+\n"
    )
    regions, _rows = read_selected_gene_regions(
        path=str(genes),
        queries=["GENE1", "ENSG2"],
        gene_id_column=4,
        gene_name_column=5,
        mode="pooled",
        pool_name="my_pool",
        label_mode="name-id",
        selected_chromosomes=None,
        min_region_length=2,
    )
    assert set(group_regions_by_state(regions)) == {"my_pool"}
