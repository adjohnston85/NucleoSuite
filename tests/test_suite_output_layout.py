from pathlib import Path


def _script(name: str) -> str:
    return Path("src/nucleosuite/resources", name).read_text()


def test_both_suites_use_the_same_numbered_tree_without_random_tail():
    required = [
        'SETUP_DIR="$OUTDIR/00_setup"',
        'GENE_SET_DIR="$OUTDIR/00_gene_sets"',
        'COMBINED_TRACK_DIR="$OUTDIR/01_combined_tracks"',
        'DAC_DIR="$OUTDIR/02_dac"',
        'DCC_DIR="$OUTDIR/03_dcc"',
        'NRL_DIR="$OUTDIR/04_nrl"',
        'TSS_AGG_DIR="$OUTDIR/06_tss_aggregation"',
        'DIST_DIR="$OUTDIR/07_distances"',
        'REGION_DIR="$OUTDIR/08_region_extract"',
        'FRAG_DIR="$OUTDIR/09_fragment_lengths"',
        'HEATMAP_DIR="$OUTDIR/10_fragment_heatmaps"',
        'GENE_EXPRESSION_DIR="$OUTDIR/11_gene_expression"',
        'POSITIVE_RUNS_DIR="$OUTDIR/12_positive_runs"',
        'PEAK_ANALYSIS_DIR="$OUTDIR/13_peak_analysis"',
    ]
    for name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = _script(name)
        assert all(token in text for token in required)
        assert "14_randomized_controls" not in text
        assert "/observed" not in text


def test_observed_support_names_are_unprefixed_and_randomized_names_are_marked():
    for name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = _script(name)
        assert 'SUPPORT_PREFIX=""' in text
        assert 'SUPPORT_PREFIX="${SAMPLE}_"' in text
        for suffix in (
            "run_parameters.tsv",
            "failed_steps.tsv",
            "run_manifest.json",
            "completion_report.tsv",
            "manifest.tsv",
        ):
            assert "${SUPPORT_PREFIX}" + suffix in text


def test_schedulers_refill_available_slots():
    for name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = _script(name)
        assert "wait -n -p completed_pid" in text
        assert "wait_one_queued_step" in text
        assert 'ASYNC_CORES="${NUCLEOSUITE_SUITE_CORES:-1}"' in text


def test_combined_wrapper_support_outputs_are_marked_only_for_randomized_runs():
    for name in ("mnase_suite.py", "cfdna_suite.py"):
        text = Path("src/nucleosuite/cli", name).read_text()
        assert 'support_prefix = f"{suite_sample}_" if randomized_mode else ""' in text
        assert 'f"{support_prefix}analysis.chrom.sizes"' in text
        assert 'f"{support_prefix}combined_chromosomes.tsv"' in text
        assert 'f"{support_prefix}combine_steps.log"' in text
