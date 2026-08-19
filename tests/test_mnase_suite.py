"""Tests for the packaged MNase full-suite workflow."""

from __future__ import annotations

import subprocess
from pathlib import Path
from importlib.resources import as_file, files


def test_packaged_mnase_suite_has_valid_bash_syntax():
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    with as_file(resource) as path:
        completed = subprocess.run(["bash", "-n", str(path)], check=False)
    assert completed.returncode == 0


def test_packaged_mnase_suite_help_mentions_region_inputs():
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    with as_file(resource) as path:
        completed = subprocess.run(
            ["bash", str(path), "--help"], check=False, text=True, capture_output=True
        )
    assert completed.returncode == 0
    assert "--states-bed" in completed.stdout
    assert "--genes-bed" in completed.stdout
    assert "--gene-set-config" in completed.stdout
    assert "--resource-set" in completed.stdout
    assert "--expression" in completed.stdout
    assert "--fragments" in completed.stdout
    assert "--randomize" in completed.stdout
    assert "--randomize-seed" in completed.stdout
    assert "--skip-randomized-controls" not in completed.stdout
    assert "--blacklist-bed" in completed.stdout
    assert "--no-blacklist" in completed.stdout
    assert "--combine-cores" in completed.stdout
    assert "--streaming-combine-cores" in completed.stdout
    assert "--indexed-combine-cores" in completed.stdout
    assert "--combine-chunk-bp" in completed.stdout
    assert "--analysis-cores" in completed.stdout
    assert "--memory-intensive-analysis-cores" in completed.stdout
    assert "--resume" in completed.stdout
    assert "--force" in completed.stdout
    assert "--dry-run" in completed.stdout
    assert "Default: nTPM" in completed.stdout
    assert "--state-distance-max" in completed.stdout
    assert "--position-percentile-interval" in completed.stdout
    assert "--gene-fft-window" in completed.stdout
    assert "--sample-name" in completed.stdout
    assert "FILE_OR_GLOB" in completed.stdout
    assert "Default: 2000" in completed.stdout
    assert "--nrl-max-distance" in completed.stdout
    assert "Default: 1200" in completed.stdout
    assert "--nrl-peak-resolution" in completed.stdout
    assert "Default: 160" in completed.stdout
    assert "--distance-x-major-tick" in completed.stdout
    assert "--distance-x-minor-tick" in completed.stdout
    assert "--score-z-limit" in completed.stdout
    assert "--distance-histogram-x-max" in completed.stdout
    assert "--percentile-boxplot-y-max" in completed.stdout
    assert "--positive-runs-plot-x-max" in completed.stdout
    assert "--skip-positive-runs" in completed.stdout
    assert "--peak-score-normalization" in completed.stdout
    assert "--skip-peak-score-frequency" in completed.stdout


def test_packaged_mnase_suite_accepts_multiple_contig_tokens(tmp_path):
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    missing_bam = tmp_path / "missing.bam"
    missing_fasta = tmp_path / "genome.fa"
    missing_fasta.write_text(">chr1\nA\n")
    with as_file(resource) as path:
        completed = subprocess.run(
            [
                "bash",
                str(path),
                "--bam",
                str(missing_bam),
                "--fasta",
                str(missing_fasta),
                "--ctcf-bed",
                str(tmp_path / "ctcf.bed"),
                "--outdir",
                str(tmp_path / "out"),
                "--contigs",
                "chr1-22",
                "chrX",
                "chrY",
            ],
            check=False,
            text=True,
            capture_output=True,
        )
    assert completed.returncode != 0
    assert "unknown option" not in completed.stderr
    assert "BAM input did not match any files" in completed.stderr


def test_packaged_mnase_suite_forwards_contigs_as_array():
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    text = resource.read_text()
    assert 'CONTIGS=("autosomes")' in text
    assert 'POSITION_PERCENTILE_INTERVAL=25' in text
    assert 'PNS_SMOOTH_WINDOW=0' in text
    assert 'PEAK_SMOOTH_WINDOW=0' in text
    assert 'PNS_MAX_NEG_RUN=0' in text
    assert 'PNS_TRACK_LIST="pns,posPNS,coverage,dyad,fragment_ends,fragment_left_ends,fragment_right_ends,pns_peaks"' in text
    assert 'MAX_DUPLICATES=1' in text
    assert 'MAX_PER_COORDINATE=0' in text
    assert 'run_step "01_combined_tracks"' in text
    assert '--spec-file "$OBS_TRACK_SPEC"' in text
    assert '--pns-max-neg-run "$PNS_MAX_NEG_RUN"' in text
    assert 'BIGBED_SCORE_SCALE=1000' in text
    assert '--bigbed-score-scale N' in text
    assert '--bigbed-score-scale "$BIGBED_SCORE_SCALE"' in text
    assert '-c "${CONTIGS[@]}"' in text
    assert '--contigs VALUE [VALUE ...]' in text
    assert 'BAM_INPUTS=()' in text
    assert '-b "${BAMS[@]}"' in text
    assert '--bam-a "${BAMS[@]}"' in text
    assert 'ANALYSIS_INPUT_ARGS=(-b "${BAMS[@]}")' in text
    assert 'ANALYSIS_INPUT_ARGS=(--fragments "${FRAGMENTS[@]}")' in text
    assert (
        '"${ANALYSIS_INPUT_ARGS[@]}" "${BLACKLIST_ARGS[@]}" --contigs'
        in text
    )


def test_packaged_mnase_suite_reports_unmatched_quoted_bam_glob(tmp_path):
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nA\n")
    with as_file(resource) as path:
        completed = subprocess.run(
            [
                "bash", str(path),
                "--bam", str(tmp_path / "sample_chr*.bam"),
                "--fasta", str(fasta),
                "--ctcf-bed", str(tmp_path / "ctcf.bed"),
                "--outdir", str(tmp_path / "out"),
                "--interval-format", "bed",
            ],
            check=False, text=True, capture_output=True,
        )
    assert completed.returncode != 0
    assert "BAM input did not match any files" in completed.stderr


def test_packaged_mnase_suite_accepts_multiple_bam_tokens_before_next_option(tmp_path):
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    bam1 = tmp_path / "sample_chr1.bam"
    bam2 = tmp_path / "sample_chr2.bam"
    bam1.touch(); bam2.touch()
    missing_fasta = tmp_path / "missing.fa"
    with as_file(resource) as path:
        completed = subprocess.run(
            [
                "bash", str(path),
                "--bam", str(bam1), str(bam2),
                "--fasta", str(missing_fasta),
                "--ctcf-bed", str(tmp_path / "ctcf.bed"),
                "--outdir", str(tmp_path / "out"),
            ],
            check=False, text=True, capture_output=True,
        )
    assert completed.returncode != 0
    assert "unknown option" not in completed.stderr
    assert "FASTA not found" in completed.stderr


def test_packaged_mnase_suite_expands_matching_quoted_bam_glob(tmp_path):
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    (tmp_path / "sample_chr1.bam").touch()
    (tmp_path / "sample_chr2.bam").touch()
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nA\n")
    missing_ctcf = tmp_path / "missing_ctcf.bed"
    with as_file(resource) as path:
        completed = subprocess.run(
            [
                "bash", str(path),
                "--bam", str(tmp_path / "sample_chr*.bam"),
                "--fasta", str(fasta),
                "--ctcf-bed", str(missing_ctcf),
                "--outdir", str(tmp_path / "out"),
                "--interval-format", "bed",
            ],
            env={**__import__("os").environ, "NUCLEOSUITE_BIN": "true"},
            check=False, text=True, capture_output=True,
        )
    assert completed.returncode != 0
    assert "BAM input did not match" not in completed.stderr
    assert "CTCF BED not found" in completed.stderr


def test_packaged_mnase_suite_runs_three_nrl_profiles_per_dac():
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    text = resource.read_text()
    assert 'nrl_${NRL_MIN_DISTANCE}_${NRL_MAX_DISTANCE}' in text
    assert 'periodicity_${SHORT_PERIODICITY_MIN}_${SHORT_PERIODICITY_MAX}' in text
    assert 'periodicity_${NUCLEOSOME_PERIODICITY_MIN}_${NUCLEOSOME_PERIODICITY_MAX}' in text
    assert '"$SHORT_PERIODICITY_MIN" "$SHORT_PERIODICITY_MAX" 0' in text
    assert '"$NUCLEOSOME_PERIODICITY_MIN" "$NUCLEOSOME_PERIODICITY_MAX" 0' in text
    assert '--peak-resolution "$peak_resolution"' in text
    assert 'NRL_PEAK_RESOLUTION=160' in text
    assert 'DISTANCE_NRL_TICK_ARGS' in text


def test_packaged_mnase_suite_uses_pns_for_optional_expression_and_state_overlays():
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    text = resource.read_text()
    assert 'set -Eeuo pipefail' in text
    assert 'if [[ -n "$EXPRESSION" && "$SKIP_GENE_EXPRESSION" -eq 0 ]]' in text
    assert '--signal-type "$signal_type" --analysis all' in text
    assert 'GENE_EXPRESSION_SIGNALS="pns"' in text
    assert 'run_gene_expression PNS pns' in text
    assert 'run_gene_expression WPS wps' in text
    assert '--peaks "${SAMPLE}=${peaks}"' in text
    assert '--signal "${SAMPLE}=${signal}"' in text
    assert 'run_gene_expression PNS pns "$PNS_CALL_NUC" "$PNS_ANALYSIS_BW"' in text
    assert 'for kind in PNS WPS' in text
    assert '--state-overlay-plot' in text
    assert '--state-color-column 9' in text
    assert '--state-overlay-smooth-window "$STATE_DISTANCE_SMOOTH_WINDOW"' in text
    assert '13_compare_positions_pns_vs_wps' in text
    assert '--bed-a "$pns" --bed-b "$wps"' in text
    assert '--percentile-interval "$interval"' in text
    assert 'PNS_percentiles_vs_all_WPS' in text
    assert 'WPS_percentiles_vs_all_PNS' in text
    assert 'return "$status"' not in text
    assert 'FAILED_STEPS_TSV="$OUTDIR/${SUPPORT_PREFIX}failed_steps.tsv"' in text

def test_packaged_mnase_suite_pools_gene_categories_and_retains_combined_chromosome_dac():
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert 'VENN_SETS="active_genes,weak_genes,repressed_genes"' in text
    assert '--leftover-set-name leftover_genes' in text
    assert 'GENE_STATE_INTERVAL="$GENE_SET_DIR/${GENE_SET_OUTPUT_PREFIX}_final_states.bed"' in text
    assert 'run_dac_scope "02_dac_${step_label}_combined_chromosomes"' in text
    assert 'run_dac_scope "02_dac_${step_label}_gene_sets"' in text
    assert '$DAC_DIR/dyads/' in text
    assert '$DAC_DIR/type_dyads/' in text
    assert '$DAC_DIR/pns' not in text
    assert '$DAC_DIR/wps' not in text

def test_packaged_mnase_suite_uses_primary_source_equivalent_wps_calls():
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    text = resource.read_text()
    assert 'WPS_CALL_PREFIX="$WPS_PREFIX"' in text
    assert 'WPS_CALL_NUC="${WPS_PREFIX}_nucleosome_regions.${INTERVAL_EXT}"' in text
    assert '--input-bigwig "$WPS_BW" --out-prefix "$WPS_CALL_PREFIX" --method wps' not in text


def test_packaged_mnase_suite_runs_positive_runs_separately_for_pns_and_wps():
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert 'POSITIVE_RUNS_DIR="$OUTDIR/12_positive_runs"' in text
    assert '12_positive_runs_pns' in text
    assert '12_positive_runs_wps' in text
    assert '--bigwig "$PNS_ANALYSIS_BW"' in text
    assert '--bigwig "$WPS_SMOOTH_BW"' in text
    assert 'PNS_POSITIVE_DIR="$POSITIVE_RUNS_DIR/pns"' in text
    assert 'WPS_POSITIVE_DIR="$POSITIVE_RUNS_DIR/wps"' in text
    assert "14_randomized_controls" not in text

def test_mnase_suite_uses_same_signal_type_dcc_and_updated_ranges():
    script = Path("src/nucleosuite/resources/mnase_full_suite.sh").read_text()
    assert 'AGGREGATE_WINDOW_HALF=2500' in script
    assert 'DCC_DMAX=500' in script
    assert 'FRAG_COUNT_MIN=100' in script
    assert 'FRAG_COUNT_MAX=1000' in script
    assert 'HEATMAP_MAX_FRAG=500' in script
    assert 'dyad_vs_${signal_type}' not in script
    assert '03_dcc_${signal_type}_range_self' in script
    assert 'left_end_vs_left_end' in script
    assert 'right_end_vs_right_end' in script
    assert '--sort-mode mean_absolute' in script

def test_mnase_suite_plots_active_mode_peak_score_frequencies():
    script = Path("src/nucleosuite/resources/mnase_full_suite.sh").read_text()
    assert 'peak-score-frequency' in script
    assert '13_peak_scores_pns_nucleosome' in script
    assert '13_peak_scores_pns_breakpoint' in script
    assert '13_peak_scores_wps_nucleosome' in script
    assert '13_peak_scores_wps_breakpoint' in script
    assert '--peaks "${RUN_MODE}=$peaks"' in script
    assert 'RAND_PNS_NUC=' not in script

def test_multicontig_combined_pass_reuses_combined_upstream_outputs():
    wrapper = Path("src/nucleosuite/cli/mnase_suite.py").read_text()
    script = Path("src/nucleosuite/resources/mnase_full_suite.sh").read_text()
    assert 'final_args.append("--resume")' in wrapper
    assert '--resume) REUSE_EXISTING_OUTPUTS=1' in script
    assert '"$REUSE_EXISTING_OUTPUTS" -eq 1' in script
    assert '"$FORCE" -eq 0 && "$REUSE_EXISTING_OUTPUTS" -eq 1' in script


def test_packaged_mnase_suite_accepts_wrapper_cores_forms(tmp_path):
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nA\n")
    missing_bam = tmp_path / "missing.bam"
    common = [
        "--bam", str(missing_bam),
        "--fasta", str(fasta),
        "--ctcf-bed", str(tmp_path / "ctcf.bed"),
        "--outdir", str(tmp_path / "out"),
    ]
    with as_file(resource) as path:
        for cores_args in (["--cores", "4"], ["--cores=4"]):
            completed = subprocess.run(
                ["bash", str(path), *common, *cores_args, "--validate-only"],
                check=False,
                text=True,
                capture_output=True,
            )
            assert completed.returncode != 0
            assert "unknown option" not in completed.stderr
            assert "BAM input did not match any files" in completed.stderr


def test_packaged_mnase_suite_rejects_invalid_wrapper_cores(tmp_path):
    resource = files("nucleosuite").joinpath("resources/mnase_full_suite.sh")
    with as_file(resource) as path:
        completed = subprocess.run(
            ["bash", str(path), "--cores", "0", "--validate-only"],
            check=False,
            text=True,
            capture_output=True,
        )
    assert completed.returncode == 2
    assert "--cores must be a positive integer" in completed.stderr


def test_peak_score_frequency_output_is_derived_after_prefix_assignment():
    for script_name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = files("nucleosuite").joinpath("resources", script_name).read_text()
        assert 'local prefix="$output_dir/${SAMPLE}_${label}" output="${prefix}_score_frequency.tsv"' not in text
        assert 'local prefix output' in text
        assert 'output="${prefix}_bins*_score_frequency.tsv"' in text


def test_cfdna_dac_validation_uses_current_dac_columns():
    text = files("nucleosuite").joinpath("resources", "cfdna_full_suite.sh").read_text()
    for column in ("Distance", "DAC Value", "DAC Value Percent", "Raw DAC Value", "Opportunities"):
        assert column in text
    assert "('distance','dac')" not in text


def test_suites_use_integer_peak_score_bins():
    for script_name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = files("nucleosuite").joinpath("resources", script_name).read_text()
        assert "--integer-bins" in text
        assert "--peak-score-bins" not in text


def test_suites_include_tissue_expression_quintile_tss_analysis():
    for script_name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = files("nucleosuite").joinpath("resources", script_name).read_text()
        assert 'TSS_EXPRESSION_TISSUE="bone_marrow"' in text
        assert 'TSS_EXPRESSION_WINDOW=2000' in text
        assert '--tss-expression-tissue' in text
        assert '--skip-tss-expression-quintiles' in text
        assert 'tss-expression-quintiles' in text
        assert 'run_tss_expression_quintiles PNS' in text
        assert 'run_tss_expression_quintiles WPS' in text


def test_cfdna_tss_quintiles_use_requested_pns_and_wps_defaults():
    text = files("nucleosuite").joinpath("resources/cfdna_full_suite.sh").read_text()
    assert 'PNS_FRAG_LOWER=137' in text
    assert 'PNS_FRAG_UPPER=197' in text
    assert 'PNS_MODE_LENGTH=167' in text
    assert 'WPS_FRAG_LOWER=120' in text
    assert 'WPS_FRAG_UPPER=180' in text
    assert 'WPS_PROTECTION=120' in text
    assert 'run_tss_expression_quintiles PNS "$PNS_BW"' in text
    assert 'run_tss_expression_quintiles WPS "$WPS_BW"' in text


def test_suites_use_pns_coverage_and_range_limited_wps_auxiliary_tracks():
    mnase = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    cfdna = files("nucleosuite").joinpath("resources/cfdna_full_suite.sh").read_text()

    for script in (mnase, cfdna):
        assert 'COVERAGE_DIR="$COMBINED_TRACK_DIR/coverage"' not in script
        assert 'PNS_COVERAGE_BW="${PNS_PREFIX}_coverage.bw"' in script
        assert '--coverage-bw "$PNS_COVERAGE_BW"' in script
        assert "'coverage,sm_mWPS,wps,wps_smoothed,mWPS,dyad,wps_peaks'" in script

    assert "printf '%s-%s\\t%s\\t%s\\trange\\n' \"$WPS_FRAG_LOWER\"" in mnase
    assert "printf '%s-%s\\t%s\\t%s\\tall\\n' \"$WPS_FRAG_LOWER\"" not in mnase
    assert '$COV_PREFIX' not in mnase

    assert "printf '%s\\t%s\\t%s\\trange\\n' \"${WPS_FRAG_LOWER}-${WPS_FRAG_UPPER}\"" in cfdna
    assert "printf '%s\\t%s\\t%s\\tall\\n' \"${WPS_FRAG_LOWER}-${WPS_FRAG_UPPER}\"" not in cfdna
    assert '$COVERAGE_PREFIX' not in cfdna


def test_standalone_wps_limits_auxiliary_tracks_to_its_fragment_range():
    text = Path("src/nucleosuite/workflows/wps.py").read_text()
    loop_start = text.index("            for fragment_start, fragment_end in fragments:")
    loop_end = text.index("            basic_tracks.cap_sparse_arrays", loop_start)
    loop = text[loop_start:loop_end]
    range_guard = "                if fragment_end - fragment_start in fragment_length_set:"
    assert range_guard in loop
    assert loop.index(range_guard) < loop.index("                    basic_tracks.add_fragment(")
    assert loop.index(range_guard) < loop.index("                    wps_scoring.add_fragment(")


def test_multicontig_suites_default_to_combine_first_execution():
    for module_name in ("mnase_suite.py", "cfdna_suite.py"):
        text = Path("src/nucleosuite/cli", module_name).read_text()
        assert 'scope = "combined-only"' in text
        assert 'worker_args.append("--combine-prerequisites-only")' in text
        assert '_replace_multi_option(final_args, "--contigs", active_contigs)' in text
        assert "strict_complete=True" in text

    for script_name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = files("nucleosuite").joinpath("resources", script_name).read_text()
        assert 'COMBINE_PREREQUISITES_ONLY=0' in text
        assert '--analysis-scope VALUE' in text
        assert 'fi  # end analytical stages skipped by combine-prerequisites-only' in text
        assert '--scope all' not in text
        assert '--scope combined_chromosomes' in text
        assert 'whole_genome' not in text
