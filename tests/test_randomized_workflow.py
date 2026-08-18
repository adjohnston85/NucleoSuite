from pathlib import Path


SUITES = ("mnase_full_suite.sh", "cfdna_full_suite.sh")


def _suite(name: str) -> str:
    return Path("src/nucleosuite/resources", name).read_text()


def test_randomized_mode_materializes_one_control_before_full_tracks():
    for name in SUITES:
        script = _suite(name)
        assert 'RUN_MODE="observed"' in script
        assert '--randomize) RUN_MODE="randomized"' in script
        assert 'SAMPLE" != *_randomized_control*' in script
        assert script.count('run_step "00_randomize_fragments"') == 1
        random_pos = script.index('run_step "00_randomize_fragments"')
        tracks_pos = script.index('run_step "01_combined_tracks"')
        assert random_pos < tracks_pos
        assert 'ANALYSIS_INPUT_ARGS=(--fragments "$RANDOM_BED")' in script
        assert '14_randomized_controls' not in script


def test_randomized_mode_requires_atomic_control_before_downstream_work():
    for name in SUITES:
        script = _suite(name)
        assert '[[ -f "$RANDOM_MARKER" && -s "$RANDOM_BED" ]] || fatal' in script
        assert 'Randomized fragment generation failed; downstream analysis was not started' in script
        assert '${SUPPORT_PREFIX}${PARAM_HASH}_00_randomize_fragments.done' in script


def test_randomized_full_tree_includes_downstream_scientific_stages():
    for name in SUITES:
        script = _suite(name)
        for token in (
            'DAC_DIR="$OUTDIR/02_dac"',
            'DCC_DIR="$OUTDIR/03_dcc"',
            'NRL_DIR="$OUTDIR/04_nrl"',
            'DIST_DIR="$OUTDIR/07_distances"',
            'GENE_EXPRESSION_DIR="$OUTDIR/11_gene_expression"',
            'POSITIVE_RUNS_DIR="$OUTDIR/12_positive_runs"',
            'PEAK_ANALYSIS_DIR="$OUTDIR/13_peak_analysis"',
        ):
            assert token in script


def test_scientific_defaults_and_cfdna_dcc_scope():
    mnase = _suite("mnase_full_suite.sh")
    cfdna = _suite("cfdna_full_suite.sh")
    assert 'PNS_FRAG_LOWER=120' in mnase
    assert 'PNS_FRAG_UPPER=180' in mnase
    assert 'PNS_MODE_LENGTH=147' in mnase
    assert 'PNS_FRAG_LOWER=137' in cfdna
    assert 'PNS_FRAG_UPPER=197' in cfdna
    for script in (mnase, cfdna):
        assert 'WPS_FRAG_LOWER=120' in script
        assert 'WPS_FRAG_UPPER=180' in script
        assert 'WPS_PROTECTION=120' in script
        assert '_maxdup' not in script
        assert '(( FRAG_COUNT_MIN < RANDOM_LOWER ))' in script
        assert '(( FRAG_COUNT_MAX > RANDOM_UPPER ))' in script
        assert '(( HEATMAP_MIN_FRAG < RANDOM_LOWER ))' in script
        assert '(( HEATMAP_MAX_FRAG > RANDOM_UPPER ))' in script
    assert 'PNS_MODE_LENGTH=167' in cfdna
    assert 'RANGE_SPECS=("145:147" "160:162" "166:168")' in cfdna
    assert 'dyad_vs_left' not in cfdna
    assert 'dyad_vs_right' not in cfdna


def test_resume_hashes_include_input_content_identities_and_source_provenance():
    for name in SUITES:
        script = _suite(name)
        assert "append_input_identities" in script
        assert "sha256_full" in script
        assert "sha256_sampled_first_middle_last_1MiB" in script
        assert "INPUT_IDENTITY_" in script
        assert "SOURCE_BAM_" in script
        assert "SOURCE_FRAGMENT_" in script
        assert '"source_inputs"' in script
        assert '"active_inputs"' in script
        assert 'active_bams = [] if randomized_bed else parameter_bams' in script
        assert '[str(Path(randomized_bed).resolve())]' in script

    for wrapper_name in ("mnase_suite.py", "cfdna_suite.py"):
        wrapper = Path("src/nucleosuite/cli", wrapper_name).read_text()
        assert 'final_args.extend(["--provenance-bam", path])' in wrapper
        assert 'final_args.extend(["--provenance-fragment", path])' in wrapper
