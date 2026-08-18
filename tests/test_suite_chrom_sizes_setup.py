from pathlib import Path


def _resource(name: str) -> str:
    return (Path(__file__).parents[1] / "src" / "nucleosuite" / "resources" / name).read_text()


def test_cfdna_suite_creates_chrom_sizes_before_prepare_regions():
    text = _resource("cfdna_full_suite.sh")
    chrom = text.index('run_step "00_chrom_sizes"')
    prepare = text.index('run_step "00_prepare_regions"')
    genes = text.index('nucleosuite gene-sets') if 'nucleosuite gene-sets' in text else text.index('$NUCLEOSUITE_BIN" gene-sets')
    assert chrom < prepare < genes
    prepare_block = text[prepare:text.index("\nPY\n", prepare)]
    assert '"$CHROM_SIZES"' in prepare_block
    assert 'Path(sizes_out).write_text' not in prepare_block


def test_mnase_suite_creates_chrom_sizes_before_tracks_and_regions():
    text = _resource("mnase_full_suite.sh")
    chrom = text.index('run_step "00_chrom_sizes"')
    prepare = text.index('run_step "00_prepare_regions"')
    genes = text.index('run_step "00_gene_sets"')
    tracks = text.index('run_step "01_combined_tracks"')
    assert chrom < prepare < genes < tracks
    prepare_block = text[prepare:text.index("\nPY\n", prepare)]
    assert '"$CHROM_SIZES"' in prepare_block
    assert 'pyBigWig' not in prepare_block
    assert 'Path(chrom_sizes).write_text' not in prepare_block


def test_chrom_size_setup_uses_selected_fasta_contigs():
    for name in ("cfdna_full_suite.sh", "mnase_full_suite.sh"):
        text = _resource(name)
        block_start = text.index('run_step "00_chrom_sizes"')
        block = text[block_start:text.index("\nPY\n", block_start)]
        assert 'pysam.FastaFile' in block
        assert 'expand_contig_tokens' in block
        assert 'analysis.chrom.sizes' in text


def test_each_suite_has_one_authoritative_chrom_size_setup_step():
    for name in ("cfdna_full_suite.sh", "mnase_full_suite.sh"):
        text = _resource(name)
        assert text.count('run_step "00_chrom_sizes"') == 1
        assert "selected.chrom.sizes" in text
