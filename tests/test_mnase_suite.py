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
        completed = subprocess.run(["bash", str(path), "--help"], check=False, text=True, capture_output=True)
    assert completed.returncode == 0
    for option in ("--states-bed", "--genes-bed", "--gene-set-config", "--resource-set", "--expression", "--fragments", "--blacklist-bed", "--no-blacklist", "--nrl-max-distance", "--nrl-peak-resolution", "--distance-adjacent-max", "--distance-long-max", "--distance-long-max-order"):
        assert option in completed.stdout
    assert "Default: 1500" in completed.stdout
    assert "Default: 160" in completed.stdout
    assert "--wps-" not in completed.stdout.lower()
    assert "--dcc-" not in completed.stdout.lower()
    assert "--exact-size N" in completed.stdout


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
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert 'CONTIGS=("autosomes")' in text
    assert 'PNS_TRACK_LIST="pns,posPNS,coverage,dyad,fragment_ends,fragment_left_ends,fragment_right_ends,pns_peaks"' in text
    assert 'FINE_FRAG_LOWER=146' in text and 'FINE_FRAG_UPPER=148' in text
    assert 'EXACT_SIZE=147' in text
    assert 'DINUC_EXACT_A=145' in text and 'DINUC_EXACT_B=147' in text
    assert '--spec-file "$OBS_TRACK_SPEC"' in text
    assert '-c "${CONTIGS[@]}"' in text
    assert 'ANALYSIS_INPUT_ARGS=(-b "${BAMS[@]}")' in text
    assert 'ANALYSIS_INPUT_ARGS=(--fragments "${FRAGMENTS[@]}")' in text
    assert 'DCC_INPUT_A' not in text


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
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert 'NRL_MIN_DISTANCE=1' in text and 'NRL_MAX_DISTANCE=1500' in text
    assert 'SHORT_PERIODICITY_MAX=144' in text
    assert 'INTERMEDIATE_PERIODICITY_MIN=150' in text and 'INTERMEDIATE_PERIODICITY_MAX=220' in text
    assert 'INTERMEDIATE_PERIODICITY_RESOLUTION=8' in text
    assert '--skip-first-peaks "$skip_first"' in text
    assert '"$NRL_MIN_DISTANCE" "$NRL_MAX_DISTANCE" "$NRL_PEAK_RESOLUTION" 1' in text
    assert '"$SHORT_PERIODICITY_MIN" "$SHORT_PERIODICITY_MAX" 1 0' in text
    assert '"$INTERMEDIATE_PERIODICITY_MIN" "$INTERMEDIATE_PERIODICITY_MAX" "$INTERMEDIATE_PERIODICITY_RESOLUTION" 0' in text


def test_packaged_mnase_suite_uses_pns_for_optional_expression_and_state_overlays():
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert 'run_gene_expression PNS pns "$PNS_CALL_NUC" "$PNS_ANALYSIS_BW"' in text
    assert 'run_gene_expression WPS' not in text
    assert '--state-overlay-plot' in text
    assert '13_compare_positions_pns_vs_wps' not in text
    assert 'WPS_CALL_NUC' not in text
    assert 'DCC_DIR=' not in text


def test_packaged_mnase_suite_pools_gene_categories_and_retains_combined_chromosome_dac():
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert 'VENN_SETS="active_genes,weak_genes,repressed_genes"' in text
    assert '--leftover-set-name leftover_genes' in text
    assert '02_dac_dyad_range_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}_combined_chromosomes' in text
    assert '02_dac_dyad_range_${FINE_FRAG_LOWER}_${FINE_FRAG_UPPER}_gene_sets' in text
    assert '$DAC_DIR/type_dyads/' not in text
    assert 'dyad_exact' not in text[text.index('# 02_dac:'):text.index('# 04_nrl:')]


def test_packaged_mnase_suite_has_no_wps_outputs():
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert 'WPS_CALL_PREFIX' not in text
    assert 'WPS_CALL_NUC' not in text
    assert 'sm_mWPS' not in text
    assert 'wps_peaks' not in text


def test_packaged_mnase_suite_runs_positive_runs_for_pns_only():
    text = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()
    assert '12_positive_runs_pns' in text
    assert '12_positive_runs_wps' not in text
    assert '--bigwig "$PNS_ANALYSIS_BW"' in text


def test_mnase_suite_removes_dcc_and_uses_updated_ranges():
    script = Path("src/nucleosuite/resources/mnase_full_suite.sh").read_text()
    assert 'FINE_FRAG_LOWER=146' in script
    assert 'FINE_FRAG_UPPER=148' in script
    assert 'EXACT_SIZE=147' in script
    assert 'DCC_DIR=' not in script
    assert ' dcc ' not in script
    assert 'FRAG_COUNT_MIN=100' in script and 'FRAG_COUNT_MAX=1000' in script


def test_mnase_suite_plots_pns_peak_score_frequencies_only():
    script = Path("src/nucleosuite/resources/mnase_full_suite.sh").read_text()
    assert '13_peak_scores_pns_nucleosome' in script
    assert '13_peak_scores_pns_breakpoint' in script
    assert '13_peak_scores_wps' not in script
    assert '--peaks "${RUN_MODE}=$peaks"' in script


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


def test_suites_include_pns_tissue_expression_quintile_tss_analysis():
    for script_name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = files("nucleosuite").joinpath("resources", script_name).read_text()
        assert 'TSS_EXPRESSION_TISSUE="bone_marrow"' in text
        assert 'tss-expression-quintiles' in text
        assert 'run_tss_expression_quintiles PNS' in text
        assert 'run_tss_expression_quintiles WPS' not in text


def test_cfdna_tss_quintiles_use_scaled_pns_and_requested_ranges():
    text = files("nucleosuite").joinpath("resources/cfdna_full_suite.sh").read_text()
    assert 'PNS_FRAG_LOWER=137' in text and 'PNS_FRAG_UPPER=197' in text and 'PNS_MODE_LENGTH=167' in text
    assert 'EXACT_LENGTHS=(145 161 167)' in text
    assert 'RANGE_SPECS=("144:146" "160:162" "166:168")' in text
    assert 'run_tss_expression_quintiles PNS "$PNS_ANALYSIS_BW"' in text
    assert 'WPS' not in text


def test_suites_scale_pns_pospns_coverage_and_peak_beds_post_combine():
    for script_name in ("mnase_full_suite.sh", "cfdna_full_suite.sh"):
        text = files("nucleosuite").joinpath("resources", script_name).read_text()
        assert 'PNS_COVERAGE_BW="${PNS_PREFIX}_coverage.bw"' in text
        assert 'PNS_POS_BW="${PNS_PREFIX}_posPNS.bw"' in text
        assert '01_scale_nucleosome_peak_scores' in text
        assert '01_scale_breakpoint_peak_scores' in text
        assert '_nucleosome_regions_mean_scaled.${INTERVAL_EXT}' in text
        assert '_breakpoint_peaks_mean_scaled.${INTERVAL_EXT}' in text
        assert '--score-column 5 --scale 100' in text
        assert 'PNS_ANALYSIS_BW="$PNS_SCALED_BW"' in text
        assert 'WPS' not in text


def test_suites_use_mean_scaled_peak_beds_for_downstream_and_no_histogram_rescaling():
    cfdna = files("nucleosuite").joinpath("resources/cfdna_full_suite.sh").read_text()
    mnase = files("nucleosuite").joinpath("resources/mnase_full_suite.sh").read_text()

    assert 'PNS_NUC="$PNS_NUC_SCALED"' in cfdna
    assert 'PNS_BRK="$PNS_BRK_SCALED"' in cfdna
    assert '--regions "$PNS_NUC_RAW" --score-column 5 --scale 100' in cfdna
    assert '--score-column 5 --score-scale 1 --integer-bins' in cfdna
    assert '--score-scale 1000' not in cfdna

    assert 'PNS_CALL_NUC="$PNS_CALL_NUC_SCALED"' in mnase
    assert 'PNS_CALL_BRK="$PNS_CALL_BRK_SCALED"' in mnase
    assert '--regions "$PNS_CALL_NUC_RAW" --score-column 5 --scale 100' in mnase
    assert '--score-column 5 --score-scale 1 --integer-bins' in mnase
    assert '--score-scale 1000' not in mnase


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
