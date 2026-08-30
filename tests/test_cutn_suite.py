import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

pyBigWig = pytest.importorskip("pyBigWig")

from nucleosuite.cutn_compare import (
    _moderated_interaction_statistics,
    compare_stage1,
)
from nucleosuite.cutn_aggregate import run_cluster_aggregate
from nucleosuite.cutn_peaks import (
    CompetitivePeak,
    ReplicatePeakStatistics,
    analyze_cutn_replicate_peaks,
    cluster_seeded_gate_peaks,
)
from nucleosuite.cutn_suite import (
    _estimate_modes,
    _locate_completed_prefix,
    _resolved_prefix,
    _scoring_fragment_range,
    _validate,
    build_parser,
    run,
)
from nucleosuite.mode_estimation import ModeEstimate
from nucleosuite.progress import ProgressReporter
from nucleosuite.peak_fdr import PeakRow


def _row(start: int, score: float, line: int = 1) -> PeakRow:
    centre = start + 40
    return PeakRow(
        ("chr1", str(start), str(start + 80), f"peak{line}", str(score), ".", str(centre), str(centre + 1)),
        score,
        line,
    )


def _cluster_record(index: int, state: str) -> ReplicatePeakStatistics:
    row = _row(index * 100, 10.0, index + 1)
    gate = state in {"S", "G"}
    pvalue = 0.01 if state == "S" else 0.2
    peak = CompetitivePeak(
        row=row,
        chrom="chr1",
        start=index * 100,
        end=index * 100 + 80,
        summit=index * 100 + 40,
        source="target",
        winner=gate,
        matched_score=50.0,
        signal_score=100.0 if gate else 40.0,
        competition_score=50.0 if gate else 0.0,
        qvalue=min(1.0, pvalue * 2),
    )
    return ReplicatePeakStatistics(
        peak=peak,
        treatment_scores=(100.0, 105.0),
        control_scores=(50.0, 55.0),
        treatment_mean=102.5,
        control_mean=52.5,
        mean_difference=50.0,
        minimum_treatment=100.0,
        maximum_control=55.0,
        conservative_excess=45.0 if gate else 0.0,
        conservative_fold_enrichment=101.0 / 56.0,
        conservative_log2_enrichment=0.85,
        all_controls_gate=gate,
        pvalue=pvalue,
        qvalue=min(1.0, pvalue * 2),
    )


def test_cutn_suite_defaults_to_pns_auto_mode_scoring_flank_and_coverage_range(
    tmp_path: Path,
):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(target),
            "--control1-bam", str(control),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    assert args.scoring_method == "pns"
    assert args.mode == "auto"
    assert args.frag_mode_padding == 30
    assert args.score_frag_lower is None
    assert args.score_frag_upper is None
    assert (args.coverage_frag_lower, args.coverage_frag_upper) == (1, 1000)
    assert args.mode_strategy == "pooled"
    assert args.bam_mode == "replicates"
    assert args.cluster_seed_p_value == 0.05
    assert args.cluster_seed_mode == "auto"
    assert args.cluster_seed_gate_mode == "auto"
    assert args.stage1_gate_mode == "auto"
    assert args.stage1_coverage_statistic == "mean"
    assert args.cluster_max_non_member_gap == 1
    assert args.min_cluster_members == 2
    assert args.cluster_aggregate_nrl_resolution == 130
    assert args.cluster_aggregate_nrl_min_order == 0
    assert args.cluster_aggregate_nrl_max_order == 3
    assert not hasattr(args, "compare_feature_level")
    assert not hasattr(args, "peak_match_distance")
    assert not hasattr(args, "stage1_control_mode")
    assert args.treatment1_bam == [str(target)]
    assert args.control1_bam == [str(control)]


def test_cutn_suite_emits_score_and_coverage_in_one_tracks_pass_per_group(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(target),
            "--control1-bam", str(control),
            "--outdir", str(tmp_path / "out"),
            "--sample-name", "sample",
            "--mode", "167",
            "--contigs", "chr1",
        ]
    )
    tracks_commands: list[list[str]] = []

    def fake_command(command: list[str]) -> None:
        if command[0] == "tracks":
            tracks_commands.append(command)
            spec_path = Path(command[command.index("--spec-file") + 1])
            lines = spec_path.read_text(encoding="utf-8").splitlines()[1:]
            assert len(lines) == 2
            score_fields = lines[0].split("\t")
            coverage_fields = lines[1].split("\t")
            assert score_fields[0] == "137-197"
            assert score_fields[2] == "pns,posPNS"
            assert coverage_fields[0] == "1-1000"
            assert coverage_fields[2] == "coverage"
            score_prefix = Path(score_fields[1])
            coverage_prefix = Path(coverage_fields[1])
            score_prefix.parent.mkdir(parents=True, exist_ok=True)
            Path(f"{score_prefix}_pns.bw").touch()
            Path(f"{score_prefix}_posPNS.bw").touch()
            Path(f"{coverage_prefix}_coverage.bw").touch()
            return
        assert command[0] == "call-peaks"
        assert command[command.index("--call-type") + 1] == "nucleosome"
        peak_prefix = Path(command[command.index("--out-prefix") + 1])
        peak_prefix.parent.mkdir(parents=True, exist_ok=True)
        Path(f"{peak_prefix}_nucleosome_regions.bed").write_text(
            "chr1\t10\t20\tpeak1\t1\t.\t15\t16\n", encoding="utf-8"
        )

    def fake_scale(_score, _reference, output, **_kwargs):
        assert str(_score).endswith('_coverage.bw')
        assert _score == _reference
        assert _kwargs['scale'] == 100.0
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).touch()
        return Path(output), 1.0, 1

    selected_clusters = tmp_path / "selected_clusters.bed"
    selected_clusters.write_text("", encoding="utf-8")
    with (
        patch("nucleosuite.cutn_suite._run_nucleosuite", side_effect=fake_command),
        patch("nucleosuite.cutn_suite.scale_bigwig_by_reference", side_effect=fake_scale),
        patch(
            "nucleosuite.cutn_suite.analyze_cutn_replicate_peaks",
            return_value={"selected_clusters": selected_clusters},
        ) as analyze_mock,
    ):
        assert run(args) == 0

    assert len(tracks_commands) == 2
    assert all(command[0] == "tracks" for command in tracks_commands)
    assert all("--output-dir" in command for command in tracks_commands)
    assert all("--scoring-method" not in command for command in tracks_commands)
    summary = (tmp_path / "out" / "sample_cutn_suite_summary.tsv").read_text()
    assert "target_raw_coverage_track\t" in summary
    assert "control_raw_coverage_track\t" in summary
    assert "condition_mean_treatment_coverage\t" in summary
    assert "condition_mean_control_coverage\t" in summary
    call = analyze_mock.call_args.kwargs
    assert call["member_gate_mode"] == "all-controls"
    assert call["cluster_member_mode"] == "seed-and-gated"
    assert len(call["target_replicate_bigwigs"]) == 1
    assert len(call["control_replicate_bigwigs"]) == 1
    manifest = json.loads((tmp_path / "out" / "cutn_stage1_manifest.json").read_text())
    assert manifest["control_candidate_peaks"] is None
    assert manifest["stage1_gate_mode"] == "all-controls"
    assert manifest["condition_mean_treatment_cluster_aggregate_score"]
    for record in manifest['treatment_replicates'] + manifest['control_replicates']:
        assert record['analysis_score'] == record['score']
        assert record['score_scaling'] == 'native'


def test_explicit_cutn_mode_overrides_automatic_estimation(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(target),
            "--control1-bam", str(control),
            "--outdir", str(tmp_path / "out"),
            "--mode", "167",
        ]
    )
    assert args.mode == 167


def test_cutn_real_bams_keep_native_pns_and_normalize_only_coverage(tmp_path):
    import numpy as np
    import pysam
    import pyBigWig
    from nucleosuite.scoring.pns import precompute_distributions

    intervals = [(500, 652), (1200, 1340), (1900, 2080), (2800, 3100)]
    paths = []
    for role in ('treatment', 'control'):
        bam = tmp_path / f'{role}.bam'
        with pysam.AlignmentFile(str(bam), 'wb', header={
            'HD': {'VN': '1.6', 'SO': 'coordinate'},
            'SQ': [{'SN': 'chr1', 'LN': 4000}],
        }) as handle:
            for i, (start, end) in enumerate(intervals):
                for first in (True, False):
                    read = pysam.AlignedSegment()
                    read.query_name = f'fragment{i}'
                    read.query_sequence = 'A' * 50
                    read.flag = 99 if first else 147
                    read.reference_id = read.next_reference_id = 0
                    read.reference_start = start if first else end - 50
                    read.next_reference_start = end - 50 if first else start
                    read.template_length = (end - start) * (1 if first else -1)
                    read.mapping_quality = 60
                    read.cigar = [(0, 50)]
                    handle.write(read)
        pysam.index(str(bam))
        paths.append(bam)
    out = tmp_path / 'out'
    args = build_parser().parse_args([
        '--treatment1-bam', str(paths[0]), '--control1-bam', str(paths[1]),
        '--outdir', str(out), '--mode', '152', '--contigs', 'chr1',
        '--cores', '1', '--skip-cluster-aggregate',
    ])
    assert run(args) == 0
    manifest = json.loads((out / 'cutn_stage1_manifest.json').read_text())
    kernels, _ = precompute_distributions([140, 152, 180], 152)
    expected = np.zeros(4000)
    for start, end in intervals[:3]:
        kernel = kernels[end-start]
        left = start-max(0, 152-(end-start))
        expected[left:left+len(kernel)] += kernel
    for record in manifest['treatment_replicates'] + manifest['control_replicates']:
        assert record['analysis_score'] == record['score']
        with pyBigWig.open(record['score']) as bw:
            observed = np.nan_to_num(bw.values('chr1', 0, 4000, numpy=True))
        np.testing.assert_allclose(observed, expected, rtol=1e-6, atol=1e-7)
        with pyBigWig.open(record['coverage']) as bw:
            coverage = np.nan_to_num(bw.values('chr1', 0, 4000, numpy=True))
        assert coverage.sum() == sum(end-start for start, end in intervals)
        with pyBigWig.open(record['scaled_coverage']) as bw:
            scaled = np.nan_to_num(bw.values('chr1', 0, 4000, numpy=True))
        assert np.isclose(scaled[scaled != 0].mean(), 100)


def test_cutn_auto_mode_prints_treatment_control_and_pooled_estimates(
    tmp_path: Path, capsys
):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(target),
            "--control1-bam", str(control),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    treatment = ModeEstimate(
        153, 152, 154, 1000, 900, True, 3, (1, 8, 1), 152, 154
    )
    control_estimate = ModeEstimate(
        151, 150, 152, 1000, 900, True, 3, (1, 8, 1), 150, 152
    )
    pooled = ModeEstimate(
        152, 151, 153, 2000, 900, True, 3, (1, 8, 1), 151, 153
    )

    with (
        patch(
            "nucleosuite.cutn_suite.estimate_bam_fragment_mode",
            side_effect=(treatment, control_estimate),
        ),
        patch(
            "nucleosuite.cutn_suite.pooled_mode_estimate",
            return_value=pooled,
        ),
    ):
        _estimate_modes(
            args,
            [(args.treatment1_bam, args.control1_bam)],
            ProgressReporter("cutn-suite"),
        )

    output = capsys.readouterr().err
    assert "Condition 1 treatment fragment mode: 153 bp" in output
    assert "Condition 1 control fragment mode: 151 bp" in output
    assert "Condition 1 pooled fragment mode: 152 bp" in output


def test_explicit_cutn_mode_does_not_validate_unused_auto_search_bounds(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(target),
            "--control1-bam", str(control),
            "--outdir", str(tmp_path / "out"),
            "--mode", "167",
            "--mode-search-lower", "300",
            "--mode-search-upper", "200",
        ]
    )
    _validate(args)


def test_explicit_cutn_score_range_overrides_mode_flank(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(target),
            "--control1-bam", str(control),
            "--outdir", str(tmp_path / "out"),
            "--mode", "167",
            "--score-frag-lower", "125",
            "--score-frag-upper", "205",
        ]
    )
    _validate(args)
    assert _scoring_fragment_range(args, 167) == (125, 205)


def test_cutn_score_padding_and_individual_bound_overrides(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()

    padded = build_parser().parse_args([
        "--treatment1-bam", str(target),
        "--control1-bam", str(control),
        "--outdir", str(tmp_path / "out1"),
        "--mode", "165",
        "--frag-mode-padding", "25",
    ])
    _validate(padded)
    assert _scoring_fragment_range(padded, 165) == (140, 190)

    lower_only = build_parser().parse_args([
        "--treatment1-bam", str(target),
        "--control1-bam", str(control),
        "--outdir", str(tmp_path / "out2"),
        "--mode", "165",
        "--score-frag-lower", "145",
    ])
    _validate(lower_only)
    assert _scoring_fragment_range(lower_only, 165) == (145, 195)

    upper_only = build_parser().parse_args([
        "--treatment1-bam", str(target),
        "--control1-bam", str(control),
        "--outdir", str(tmp_path / "out3"),
        "--mode", "165",
        "--score-frag-upper", "188",
    ])
    _validate(upper_only)
    assert _scoring_fragment_range(upper_only, 165) == (135, 188)


def test_cutn_suite_locates_parameterized_multicontig_score_outputs(tmp_path: Path):
    requested = tmp_path / "tracks" / "cutn_target"
    direct = _resolved_prefix(requested, "pns", 153, 120, 500)
    root = requested.parent / f"{requested.name}_multicontig"
    combined = root / "combined"
    combined.mkdir(parents=True)
    combined_name = direct.name
    (root / "nucleosuite_multicontig_manifest.json").write_text(
        json.dumps(
            {"combined_dir": str(combined), "combined_name": combined_name}
        ),
        encoding="utf-8",
    )
    for suffix in ("_pns.bw", "_posPNS.bw", "_coverage.bw"):
        (combined / f"{combined_name}{suffix}").touch()

    located = _locate_completed_prefix(
        requested, direct, ("_pns.bw", "_posPNS.bw", "_coverage.bw")
    )

    assert located == combined / combined_name


def test_cutn_suite_locates_multicontig_peak_outputs(tmp_path: Path):
    requested = tmp_path / "peaks" / "cutn_target_pns_mean_scaled"
    root = requested.parent / f"{requested.name}_multicontig"
    combined = root / "combined"
    combined.mkdir(parents=True)
    (root / "nucleosuite_multicontig_manifest.json").write_text(
        json.dumps(
            {"combined_dir": str(combined), "combined_name": requested.name}
        ),
        encoding="utf-8",
    )
    (combined / f"{requested.name}_nucleosome_regions.bed").touch()

    located = _locate_completed_prefix(
        requested, requested, ("_nucleosome_regions.bed",)
    )

    assert located == combined / requested.name


def test_cutn_suite_prefers_existing_serial_outputs(tmp_path: Path):
    requested = tmp_path / "tracks" / "cutn_target"
    direct = _resolved_prefix(requested, "pns", 153, 120, 500)
    direct.parent.mkdir(parents=True)
    for suffix in ("_pns.bw", "_posPNS.bw", "_coverage.bw"):
        Path(f"{direct}{suffix}").touch()

    located = _locate_completed_prefix(
        requested, direct, ("_pns.bw", "_posPNS.bw", "_coverage.bw")
    )

    assert located == direct






def _write_bigwig(path: Path, values: list[float]) -> None:
    handle = pyBigWig.open(str(path), "w")
    handle.addHeader([("chr1", len(values))])
    handle.addEntries("chr1", 0, values=values, span=1, step=1)
    handle.close()










def test_stage1_replicate_statistics_report_raw_p_without_control_peak_calls(tmp_path: Path):
    target_bed = tmp_path / "target_candidates.bed"
    target_bed.write_text(
        "chr1\t100\t180\tpeak1\t8\t.\t140\t141\n", encoding="utf-8"
    )
    treatment_paths = [tmp_path / f"treatment{index}.bw" for index in range(2)]
    control_paths = [tmp_path / f"control{index}.bw" for index in range(2)]
    for path, score in zip(treatment_paths, (120.0, 120.0)):
        _write_bigwig(
            path,
            [score if 100 <= position < 180 else 0.0 for position in range(300)],
        )
    for path, score in zip(control_paths, (50.0, 50.0)):
        _write_bigwig(
            path,
            [score if 100 <= position < 180 else 0.0 for position in range(300)],
        )
    mean_path = tmp_path / "treatment_mean.bw"
    _write_bigwig(
        mean_path,
        [120.0 if 100 <= position < 180 else 0.0 for position in range(300)],
    )

    outputs = analyze_cutn_replicate_peaks(
        target_bed,
        output_dir=tmp_path / "results",
        target_replicate_bigwigs=treatment_paths,
        control_replicate_bigwigs=control_paths,
        target_mean_bigwig=mean_path,
        seed_mode="pvalue",
        seed_gate_mode="mean",
        member_gate_mode="all-controls",
        compute_pvalues=True,
        minimum_cluster_members=1,
    )

    annotated = outputs["annotated_peaks"].read_text().strip().split("\t")
    significant = outputs["significant_peaks"].read_text().strip().split("\t")
    statistics = outputs["competition_table"].read_text().splitlines()
    assert annotated[4] == "120"
    assert len(annotated) == 9
    assert annotated[-1] == "0"
    assert significant == annotated
    seed = outputs["seed_peaks"].read_text().strip().split("\t")
    assert seed == annotated
    assert outputs["annotated_peaks"].name == "target_peaks_replicate_statistics.bed"
    assert outputs["seed_peaks"].name == "target_seed_peaks_S-pvalue-mean.bed"
    assert outputs["competition_table"].name == "target_peak_replicate_statistics.tsv"
    assert "treatment_replicate_means" in statistics[0]
    assert "control_replicate_means" in statistics[0]
    header = statistics[0].split("\t")
    row = statistics[1].split("\t")
    assert row[header.index("is_seed")] == "true"
    assert row[header.index("is_gated_member")] == "true"
    assert row[header.index("p_value")] == "0"


def test_stage1_single_replicate_can_use_gate_only_seeding(tmp_path: Path):
    target_bed = tmp_path / "target_candidates.bed"
    target_bed.write_text(
        "chr1\t100\t180\tpeak1\t8\t.\t140\t141\n", encoding="utf-8"
    )
    treatment = tmp_path / "treatment.bw"
    control = tmp_path / "control.bw"
    mean_path = tmp_path / "treatment_mean.bw"
    for path, score in ((treatment, 120.0), (control, 50.0), (mean_path, 120.0)):
        _write_bigwig(
            path,
            [score if 100 <= position < 180 else 0.0 for position in range(300)],
        )

    outputs = analyze_cutn_replicate_peaks(
        target_bed,
        output_dir=tmp_path / "results",
        target_replicate_bigwigs=[treatment],
        control_replicate_bigwigs=[control],
        target_mean_bigwig=mean_path,
        seed_mode="gated",
        seed_gate_mode="all-controls",
        member_gate_mode="all-controls",
        compute_pvalues=False,
        minimum_cluster_members=1,
    )

    assert outputs["annotated_peaks"].read_text().strip().endswith("\t.")
    assert outputs["selected_peaks"].read_text().strip().endswith("\t.")
    assert outputs["seed_peaks"].read_text().strip().endswith("\t.")
    assert outputs["selected_clusters"].read_text().strip()


def test_cluster_aggregate_writes_combined_heatmap_profiles_and_directional_nrl(
    tmp_path: Path,
):
    length = 1200
    positions = np.arange(length)
    signal1 = (
        4.0 * np.cos(2.0 * np.pi * (positions - 400) / 180.0)
        + 5.0 * np.exp(-0.5 * ((positions - 400) / 20.0) ** 2)
    ).tolist()
    signal2 = (
        3.5 * np.cos(2.0 * np.pi * (positions - 800) / 180.0)
        + 4.5 * np.exp(-0.5 * ((positions - 800) / 20.0) ** 2)
    ).tolist()
    mean = ((np.asarray(signal1) + np.asarray(signal2)) / 2.0).tolist()
    paths = [tmp_path / "rep1.bw", tmp_path / "rep2.bw", tmp_path / "mean.bw"]
    for path, values in zip(paths, (signal1, signal2, mean)):
        _write_bigwig(path, values)
    anchors = tmp_path / "anchors.bed"
    anchors.write_text(
        "chr1\t400\t401\tcluster1\t1\t.\t400\t401\n"
        "chr1\t800\t801\tcluster2\t1\t.\t800\t801\n",
        encoding="utf-8",
    )

    outputs = run_cluster_aggregate(
        mean_score=paths[2],
        replicate_scores=paths[:2],
        anchor_bed=anchors,
        output_dir=tmp_path / "aggregate",
        label="treatment",
        window_half=300,
        maximum_heatmap_rows=10,
        bootstrap_replicates=5,
        nrl_peak_resolution=130,
        nrl_min_order=0,
        nrl_max_order=3,
    )

    assert outputs["status"] == "complete"
    combined = outputs["combined"]
    assert Path(combined["heatmap"]).is_file()
    assert Path(combined["heatmap_matrix"]).is_file()
    assert Path(combined["nrl_summary"]).is_file()
    nrl_summary = Path(combined["nrl_summary"]).read_text().splitlines()
    assert "regression_min_peak_order" in nrl_summary[0]
    assert "\t0\t3\t" in nrl_summary[1]
    assert Path(outputs["replicate_overlay_plot"]).is_file()
    assert Path(outputs["bootstrap_profile"]).is_file()


def test_stage1_clusters_use_p_seeds_and_all_gated_members(tmp_path: Path):
    target_bed = tmp_path / "target_candidates.bed"
    target_bed.write_text(
        "chr1\t20\t60\tpeak1\t8\t.\t40\t41\n"
        "chr1\t100\t140\tpeak2\t8\t.\t120\t121\n"
        "chr1\t180\t220\tpeak3\t8\t.\t200\t201\n",
        encoding="utf-8",
    )
    treatment_paths = [tmp_path / f"treatment{index}.bw" for index in range(2)]
    control_paths = [tmp_path / f"control{index}.bw" for index in range(2)]
    treatment_scores = ((120.0, 120.0, 150.0), (122.0, 200.0, 152.0))
    control_scores = ((50.0, 50.0, 60.0), (52.0, 110.0, 62.0))

    def values(scores: tuple[float, float, float]) -> list[float]:
        output = [0.0] * 300
        for (start, end), score in zip(((20, 60), (100, 140), (180, 220)), scores):
            output[start:end] = [score] * (end - start)
        return output

    for path, scores in zip(treatment_paths, treatment_scores):
        _write_bigwig(path, values(scores))
    for path, scores in zip(control_paths, control_scores):
        _write_bigwig(path, values(scores))
    mean_path = tmp_path / "treatment_mean.bw"
    _write_bigwig(mean_path, values((121.0, 160.0, 151.0)))

    outputs = analyze_cutn_replicate_peaks(
        target_bed,
        output_dir=tmp_path / "results",
        target_replicate_bigwigs=treatment_paths,
        control_replicate_bigwigs=control_paths,
        target_mean_bigwig=mean_path,
        cluster_seed_pvalue=0.05,
        minimum_cluster_members=2,
    )

    statistics = outputs["competition_table"].read_text().splitlines()
    header = statistics[0].split("\t")
    rows = [line.split("\t") for line in statistics[1:]]
    assert all(row[header.index("is_gated_member")] == "true" for row in rows)
    assert float(rows[0][header.index("p_value")]) < 0.05
    assert float(rows[1][header.index("p_value")]) > 0.05
    assert float(rows[2][header.index("p_value")]) < 0.05
    lines = outputs["cluster_table"].read_text().splitlines()
    cluster_header = lines[0].split("\t")
    cluster = lines[1].split("\t")
    assert cluster[:4] == ["cutn_cluster_1", "chr1", "20", "220"]
    assert cluster[cluster_header.index("seed_peak_count")] == "2"
    assert cluster[cluster_header.index("member_count")] == "3"
    assert cluster[cluster_header.index("bridged_non_member_peak_count")] == "0"


def test_seeded_clusters_end_at_last_gated_member_and_require_a_seed():
    records = [_cluster_record(index, state) for index, state in enumerate("SGxxGSG")]
    clusters = cluster_seeded_gate_peaks(records)
    assert [[peak.summit for peak in cluster.significant_peaks] for cluster in clusters] == [
        [40, 140],
        [440, 540, 640],
    ]
    assert [(cluster.start, cluster.end) for cluster in clusters] == [(0, 180), (400, 680)]
    assert [cluster.bridged_non_member_peak_count for cluster in clusters] == [0, 0]

    records = [_cluster_record(index, state) for index, state in enumerate("GSGxxGG")]
    clusters = cluster_seeded_gate_peaks(records)
    assert len(clusters) == 1
    assert [peak.summit for peak in clusters[0].significant_peaks] == [40, 140, 240]


def test_seeded_clusters_bridge_one_non_gated_peak_without_counting_it():
    records = [_cluster_record(index, state) for index, state in enumerate("SxG")]
    clusters = cluster_seeded_gate_peaks(records)
    assert len(clusters) == 1
    cluster = clusters[0]
    assert [peak.summit for peak in cluster.significant_peaks] == [40, 240]
    assert cluster.bridged_non_member_peak_count == 1
    assert cluster.seed_peak_count == 1
    assert cluster.score == 90.0

    assert cluster_seeded_gate_peaks([_cluster_record(0, "S")]) == []
    assert cluster_seeded_gate_peaks(
        [_cluster_record(0, "G"), _cluster_record(1, "G")]
    ) == []


def test_seeded_clusters_enforce_1000_bp_adjacent_gated_summit_limit():
    at_limit = cluster_seeded_gate_peaks(
        [_cluster_record(0, "S"), _cluster_record(10, "G")],
        max_cluster_gap=1000,
    )
    assert len(at_limit) == 1

    beyond_limit = cluster_seeded_gate_peaks(
        [_cluster_record(0, "S"), _cluster_record(11, "G")],
        max_cluster_gap=1000,
    )
    assert beyond_limit == []


def test_stage2_empirical_bayes_moderates_region_variances():
    groups = []
    for index in range(30):
        offset = index * 0.01
        groups.append(
            [
                [2.0 + offset, 2.2 + offset],
                [1.0, 1.1],
                [3.0 + offset, 3.3 + offset],
                [1.0, 1.2],
            ]
        )
    statistics, metadata = _moderated_interaction_statistics(groups)
    assert metadata["available"] is True
    assert metadata["moderated"] is True
    assert float(metadata["prior_degrees_freedom"]) > 0
    assert len(statistics) == 30
    assert all(0 <= row["moderated_pvalue"] <= 1 for row in statistics)


def test_replicate_mode_allows_unequal_unpaired_treatment_control_groups(tmp_path: Path):
    paths = [tmp_path / name for name in ("t1.bam", "t2.bam", "c1.bam")]
    for path in paths:
        path.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(paths[0]), str(paths[1]),
            "--control1-bam", str(paths[2]),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    _validate(args)


def test_condition2_must_be_complete(tmp_path: Path):
    paths = [tmp_path / name for name in ("t1.bam", "c1.bam", "t2.bam")]
    for path in paths:
        path.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(paths[0]),
            "--control1-bam", str(paths[1]),
            "--treatment2-bam", str(paths[2]),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    with pytest.raises(ValueError, match="must be supplied together"):
        _validate(args)


def test_four_group_dry_run_plans_stage2(tmp_path: Path, capsys):
    paths = [tmp_path / f"group{index}.bam" for index in range(4)]
    for path in paths:
        path.touch()
    args = build_parser().parse_args(
        [
            "--treatment1-bam", str(paths[0]),
            "--control1-bam", str(paths[1]),
            "--treatment2-bam", str(paths[2]),
            "--control2-bam", str(paths[3]),
            "--outdir", str(tmp_path / "out"),
            "--dry-run",
        ]
    )
    assert run(args) == 0
    output = capsys.readouterr().out
    assert "conditions\t2" in output
    assert "stage2-log2-moderated-four-group-interaction-bh" in output


def _stage1_manifest(
    root: Path,
    name: str,
    treatment_values: list[list[float]],
    control_values: list[list[float]],
) -> Path:
    root.mkdir(parents=True)
    significant_peaks = root / "peaks.bed"
    significant_clusters = root / "clusters.bed"
    significant_peaks.write_text(
        "chr1\t80\t180\tpeak\t10\t.\t130\t131\t0.01\n", encoding="utf-8"
    )
    significant_clusters.write_text(
        "chr1\t60\t220\tcluster\t10\t.\t130\t131\t0.01\n", encoding="utf-8"
    )
    treatment_records = []
    control_records = []
    for index, target in enumerate(treatment_values, 1):
        target_path = root / f"target_{index}.bw"
        _write_bigwig(target_path, target)
        treatment_records.append(
            {
                "replicate": index,
                "analysis_score": str(target_path),
                "scaled_coverage": str(target_path),
            }
        )
    for index, control in enumerate(control_values, 1):
        control_path = root / f"control_{index}.bw"
        _write_bigwig(control_path, control)
        control_records.append(
            {
                "replicate": index,
                "analysis_score": str(control_path),
                "scaled_coverage": str(control_path),
            }
        )
    manifest = {
        "schema": "nucleosuite_cutn_stage1",
        "schema_version": 2,
        "condition_name": name,
        "bam_mode": "replicates",
        "scoring_method": "pns",
        "score_track": "pns",
        "positive_track": "posPNS",
        "peak_discovery_track": "pns",
        "peak_measurement_track": "coverage_divided_by_nonzero_mean_x100",
        "target_mode": 167,
        "control_mode": 167,
        "frag_lower": 120,
        "frag_upper": 500,
        "contigs": ["chr1"],
        "treatment_replicates": treatment_records,
        "control_replicates": control_records,
        "condition_mean_treatment_score": treatment_records[0]["analysis_score"],
        "condition_mean_control_score": control_records[0]["analysis_score"],
        "condition_mean_treatment_coverage": treatment_records[0]["scaled_coverage"],
        "condition_mean_control_coverage": control_records[0]["scaled_coverage"],
        "significant_peaks": str(significant_peaks),
        "significant_clusters": str(significant_clusters),
    }
    path = root / "cutn_stage1_manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def test_cutn_compare_uses_scaled_bigwig_replicates_without_bams(tmp_path: Path):
    length = 300
    first_treatment = [[0.0] * length for _ in range(2)]
    first_control = [[0.0] * length for _ in range(2)]
    second_treatment = [[0.0] * length for _ in range(2)]
    second_control = [[0.0] * length for _ in range(2)]
    first_treatment[0][80:180] = [1.0] * 100
    first_treatment[1][80:180] = [2.0] * 100
    second_treatment[0][80:180] = [10.0] * 100
    second_treatment[1][80:180] = [12.0] * 100
    first = _stage1_manifest(tmp_path / "first", "wild_type", first_treatment, first_control)
    second = _stage1_manifest(tmp_path / "second", "mutant", second_treatment, second_control)

    comparison = compare_stage1(first, second, outdir=tmp_path / "comparison")

    payload = json.loads(comparison.read_text())
    assert payload["feature_level"] == "clusters"
    assert "peaks" not in payload["results"]
    assert payload["results"]["clusters"]["inferential_fdr_available"] is True
    row = (tmp_path / "comparison" / "differential_clusters.tsv").read_text().splitlines()[1]
    fields = row.split("\t")
    assert fields[-1] == "significant_gain"
    header = (tmp_path / "comparison" / "differential_clusters.tsv").read_text().splitlines()[0].split("\t")
    assert fields[header.index("region_origin")] == "overlap_union"
    gain_fields = (
        tmp_path / "comparison" / "differential_clusters_fdr0.05_gains.bed"
    ).read_text().strip().split("\t")
    assert gain_fields[-1] == "overlap_union"
    assert Path(payload["cluster_overlap"]["summary"]).is_file()
    assert Path(payload["cluster_overlap"]["venn_plot"]).is_file()
    assert payload["cluster_aligned_aggregates"]["status"] == "unavailable"


def test_cutn_compare_uses_unpaired_unequal_four_group_interaction(tmp_path: Path):
    length = 300

    def tracks(count: int, score: float) -> list[list[float]]:
        output = [[0.0] * length for _ in range(count)]
        for values in output:
            values[80:180] = [score] * 100
        return output

    first = _stage1_manifest(
        tmp_path / "first",
        "wild_type",
        tracks(3, 10.0),
        tracks(2, 2.0),
    )
    second = _stage1_manifest(
        tmp_path / "second",
        "mutant",
        tracks(2, 30.0),
        tracks(3, 4.0),
    )

    comparison = compare_stage1(first, second, outdir=tmp_path / "comparison")

    payload = json.loads(comparison.read_text())
    assert payload["results"]["clusters"]["inferential_fdr_available"] is True
    lines = (tmp_path / "comparison" / "differential_clusters.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    assert "condition1_treatment_replicate_3" in header
    assert "condition1_control_replicate_3" not in header
    assert "condition2_control_replicate_3" in header
    # The cluster spans 160 bp, of which 100 bp carries coverage.
    expected = ((30.0 - 4.0) - (10.0 - 2.0)) * 100.0 / 160.0
    assert float(row[header.index("raw_interaction_difference")]) == expected
    assert float(row[header.index("log2_interaction_difference")]) > 0
    assert row[header.index("replicate_consistency")] == "robust_gain"
    assert row[-1] == "significant_gain"


def test_cutn_compare_uses_connected_union_for_overlapping_cluster_coordinates(tmp_path: Path):
    length = 300
    values = [[0.0] * length for _ in range(2)]
    first = _stage1_manifest(tmp_path / "first", "first", values, values)
    second = _stage1_manifest(tmp_path / "second", "second", values, values)
    (tmp_path / "first" / "clusters.bed").write_text(
        "chr1\t80\t150\tcluster1\t10\t.\t120\t121\t0.01\n", encoding="utf-8"
    )
    (tmp_path / "second" / "clusters.bed").write_text(
        "chr1\t120\t200\tcluster2\t10\t.\t160\t161\t0.01\n", encoding="utf-8"
    )

    compare_stage1(first, second, outdir=tmp_path / "comparison")

    lines = (tmp_path / "comparison" / "differential_clusters.tsv").read_text().splitlines()
    header, fields = lines[0].split("\t"), lines[1].split("\t")
    assert fields[:3] == ["chr1", "80", "200"]
    assert fields[header.index("condition1_stage1_support")] == "true"
    assert fields[header.index("condition2_stage1_support")] == "true"
    assert fields[header.index("region_origin")] == "overlap_union"
    assert fields[header.index("condition1_cluster_ids")] == "cluster1"
    assert fields[header.index("condition2_cluster_ids")] == "cluster2"


def test_cutn_compare_retains_condition_specific_cluster_regions(tmp_path: Path):
    length = 300
    values = [[0.0] * length for _ in range(2)]
    first = _stage1_manifest(tmp_path / "first", "first", values, values)
    second = _stage1_manifest(tmp_path / "second", "second", values, values)
    (tmp_path / "first" / "clusters.bed").write_text(
        "chr1\t80\t120\tcluster1\t10\t.\t100\t101\t0.01\n", encoding="utf-8"
    )
    (tmp_path / "second" / "clusters.bed").write_text(
        "chr1\t180\t220\tcluster2\t10\t.\t200\t201\t0.01\n", encoding="utf-8"
    )

    compare_stage1(first, second, outdir=tmp_path / "comparison")

    lines = (tmp_path / "comparison" / "differential_clusters.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    coordinates = []
    for line in lines[1:]:
        row = line.split("\t")
        coordinates.append(
            [
                *row[:4],
                row[header.index("condition1_stage1_support")],
                row[header.index("condition2_stage1_support")],
                row[header.index("region_origin")],
            ]
        )
    assert coordinates == [
        ["chr1", "80", "120", "100", "true", "false", "condition1_only"],
        ["chr1", "180", "220", "200", "false", "true", "condition2_only"],
    ]


def test_cutn_compare_does_not_merge_nearby_nonoverlapping_clusters(tmp_path: Path):
    length = 300
    values = [[0.0] * length for _ in range(2)]
    first = _stage1_manifest(tmp_path / "first", "first", values, values)
    second = _stage1_manifest(tmp_path / "second", "second", values, values)
    (tmp_path / "first" / "clusters.bed").write_text(
        "chr1\t80\t120\tcluster1\t10\t.\t100\t101\t0.01\n", encoding="utf-8"
    )
    (tmp_path / "second" / "clusters.bed").write_text(
        "chr1\t130\t170\tcluster2\t10\t.\t140\t141\t0.01\n", encoding="utf-8"
    )

    manifest = compare_stage1(first, second, outdir=tmp_path / "comparison")

    lines = (tmp_path / "comparison" / "differential_clusters.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    origins = [row.split("\t")[header.index("region_origin")] for row in lines[1:]]
    assert origins == ["condition1_only", "condition2_only"]
    payload = json.loads(manifest.read_text())
    assert payload["results"]["clusters"]["region_origin_counts"] == {
        "condition1_only": 1,
        "condition2_only": 1,
    }
    summary = Path(payload["cluster_overlap"]["summary"]).read_text()
    assert "overlapping_cluster_bp\t0" in summary


def test_cutn_compare_preserves_one_to_many_cluster_overlap_and_base_counts(
    tmp_path: Path,
):
    values = [[0.0] * 300 for _ in range(2)]
    first = _stage1_manifest(tmp_path / "first", "first", values, values)
    second = _stage1_manifest(tmp_path / "second", "second", values, values)
    (tmp_path / "first" / "clusters.bed").write_text(
        "chr1\t40\t240\tfirst_1\t10\t.\t100\t101\t0.01\n",
        encoding="utf-8",
    )
    (tmp_path / "second" / "clusters.bed").write_text(
        "chr1\t50\t110\tsecond_1\t8\t.\t80\t81\t0.01\n"
        "chr1\t170\t230\tsecond_2\t9\t.\t200\t201\t0.01\n",
        encoding="utf-8",
    )

    manifest = compare_stage1(first, second, outdir=tmp_path / "comparison")
    payload = json.loads(manifest.read_text())
    mapping = Path(payload["cluster_overlap"]["component_mapping"]).read_text().splitlines()
    header, row = mapping[0].split("\t"), mapping[1].split("\t")
    assert row[header.index("relationship")] == "1_to_many"
    assert row[header.index("condition1_cluster_ids")] == "first_1"
    assert row[header.index("condition2_cluster_ids")] == "second_1;second_2"
    summary = Path(payload["cluster_overlap"]["summary"]).read_text()
    assert "condition1_cluster_bp\t200" in summary
    assert "condition2_cluster_bp\t120" in summary
    assert "overlapping_cluster_bp\t120" in summary
    assert "union_cluster_bp\t200" in summary
