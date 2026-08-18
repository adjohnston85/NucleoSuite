from nucleosuite.core.regions import expand_contig_tokens


def test_autosomes_are_header_driven():
    references = ["chr1", "chr2", "chrX", "chrM", "scaffold_7"]
    assert expand_contig_tokens(["autosomes"], references) == ["chr1", "chr2"]


def test_numeric_range_resolves_chr_names():
    references = ["chr1", "chr2", "chr3"]
    assert expand_contig_tokens(["1-3"], references) == references
