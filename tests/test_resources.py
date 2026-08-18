"""Tests for bundled package resources."""

from __future__ import annotations

from importlib.resources import as_file, files

from nucleosuite.resource_files import RESOURCE_FILES


def test_all_declared_resources_are_packaged():
    for relative in RESOURCE_FILES.values():
        resource = files("nucleosuite").joinpath(relative)
        with as_file(resource) as path:
            assert path.is_file()
            assert path.stat().st_size > 0


def test_hg19_gene_resource_has_expression_matching_columns():
    from nucleosuite.resource_files import materialized_resource_path

    with materialized_resource_path("hg19-genes") as path:
        seen = set()
        count = 0
        with path.open() as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                assert len(fields) >= 6
                assert fields[3].startswith("ENSG")
                assert fields[5] in {"+", "-"}
                assert fields[3] not in seen
                seen.add(fields[3])
                count += 1
        assert count == 19396


def test_hg19_ctcf_resource_is_stranded_bed6():
    from nucleosuite.resource_files import materialized_resource_path

    with materialized_resource_path("gm12878-hg19-ctcf") as path:
        names = set()
        count = 0
        with path.open() as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                assert len(fields) == 6
                assert fields[0].startswith("chr")
                assert int(fields[1]) < int(fields[2])
                assert fields[3].startswith("CTCF_MA0139.1_")
                assert 0 <= int(fields[4]) <= 1000
                assert fields[5] in {"+", "-"}
                assert fields[3] not in names
                names.add(fields[3])
                count += 1
        assert count == 21536


def test_bundled_expression_resources_are_available():
    import csv
    from nucleosuite.io import open_text
    from nucleosuite.resource_files import materialized_resource_path

    with materialized_resource_path("hpa-tissue-expression") as path:
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            assert reader.fieldnames == ["Gene", "Gene name", "Tissue", "nTPM"]
            assert next(reader)["Tissue"]

    with materialized_resource_path("hpa-cell-line-metadata") as path:
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            assert "Cellosaurus ID" in (reader.fieldnames or [])
            assert next(reader)["Cell line"]


def test_official_hg19_blacklist_is_packaged_and_verified():
    import gzip
    import hashlib
    from nucleosuite.resource_files import load_manifest, materialized_resource_path

    expected = load_manifest()["resources"]["hg19-blacklist-v2"]["sha256"]
    with materialized_resource_path("hg19-blacklist-v2") as path:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
        with gzip.open(path, "rt") as handle:
            rows = [line for line in handle if line.strip()]
    assert len(rows) == 834


def test_default_gene_set_resource_encodes_requested_precedence():
    from nucleosuite.resource_files import materialized_resource_path

    with materialized_resource_path("default-gene-sets") as path:
        rows = path.read_text().splitlines()
    assert rows[1].endswith("\trepressed_genes")
    assert rows[2].endswith("\tactive_genes,repressed_genes")
    assert rows[3] == "repressed_genes\t12_Repressed\tactive_genes,weak_genes"


def test_chromatin_state_resource_has_distinct_requested_colours():
    from nucleosuite.resource_files import materialized_resource_path

    colours = {}
    with materialized_resource_path("gm12878-hg19-states") as path:
        with path.open() as handle:
            for raw in handle:
                fields = raw.rstrip("\n").split("\t")
                if fields[3] in colours:
                    assert colours[fields[3]] == fields[8]
                else:
                    colours[fields[3]] = fields[8]

    assert colours["4_Strong_Enhancer"] == "250,202,0"
    assert colours["5_Strong_Enhancer"] == "230,175,0"
    assert colours["6_Weak_Enhancer"] == "255,252,4"
    assert colours["7_Weak_Enhancer"] == "230,220,0"
    assert colours["9_Txn_Transition"] == "0,176,80"
    assert colours["10_Txn_Elongation"] == "0,125,60"
    assert colours["13_Heterochrom/lo"] == "169,169,169"
    assert colours["14_Repetitive/CNV"] == "245,245,245"
    assert colours["15_Repetitive/CNV"] == "220,220,220"
