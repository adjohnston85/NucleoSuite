import json
from pathlib import Path
from unittest.mock import patch

import pytest

from nucleosuite import cutn_suite


def _source_manifest(root: Path, *, condition_name: str = "condition1") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    treatment = []
    control = []
    for role, records in (("treatment", treatment), ("control", control)):
        for index in (1, 2):
            bam = root / f"{role}_{index}.bam"
            bam.touch()
            score = root / f"{role}_{index}_scaled_score.bw"
            coverage = root / f"{role}_{index}_scaled_coverage.bw"
            raw_coverage = root / f"{role}_{index}_coverage.bw"
            score.touch(); coverage.touch(); raw_coverage.touch()
            counts = root / f"{role}_{index}_fragment_length_counts.tsv"
            counts.write_text(
                "fragment_length\tcount\n150\t2\n152\t7\n154\t3\n",
                encoding="utf-8",
            )
            summary = root / f"{role}_{index}_fragment_summary.tsv"
            summary.write_text(
                "metric\tvalue\ntotal_fragments_used_in_range\t12\n",
                encoding="utf-8",
            )
            records.append(
                {
                    "replicate": index,
                    "bams": [str(bam)],
                    "score": str(score),
                    "positive_score": str(score),
                    "coverage": str(raw_coverage),
                    "analysis_score": str(score),
                    "scaled_coverage": str(coverage),
                    "fragment_length_counts": str(counts),
                    "fragment_summary": str(summary),
                    "positive_score_mean": 0.5,
                    "coverage_nonzero_mean": 2.0,
                }
            )
    manifest = {
        "schema": "nucleosuite-cutn-stage1",
        "schema_version": 4,
        "nucleosuite_version": "0.10.12",
        "condition_name": condition_name,
        "sample_name": "cutn",
        "bam_mode": "replicates",
        "scoring_method": "pns",
        "score_track": "pns",
        "positive_track": "posPNS",
        "target_mode": 152,
        "control_mode": 152,
        "frag_lower": 122,
        "frag_upper": 182,
        "score_fragment_flank": 30,
        "target_score_frag_lower": 122,
        "target_score_frag_upper": 182,
        "control_score_frag_lower": 122,
        "control_score_frag_upper": 182,
        "coverage_frag_lower": 1,
        "coverage_frag_upper": 1000,
        "contigs": ["chr1"],
        "blacklist_bed": None,
        "max_duplicates": 1,
        "dedup_scope": "all_bams",
        "cluster_seed_p_value": 0.05,
        "cluster_max_non_member_gap": 1,
        "max_cluster_gap": 1000,
        "minimum_cluster_members": 2,
        "cluster_seed_mode_requested": "auto",
        "cluster_seed_gate_mode_requested": "auto",
        "stage1_gate_mode_requested": "auto",
        "stage1_coverage_statistic": "mean",
        "stage1_gate_mode": "all-controls",
        "cluster_member_mode": "seed-and-gated",
        "treatment_replicates": treatment,
        "control_replicates": control,
        "cluster_aggregate_parameters": {
            "window_half": 1000,
            "maximum_heatmap_rows": 5000,
            "bootstrap_replicates": 200,
            "nrl_peak_resolution": 140.0,
            "nrl_min_order": 0,
            "nrl_max_order": 3,
        },
    }
    path = root / "cutn_stage1_manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def test_condition_labels_default_to_condition1_condition2(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = cutn_suite.build_parser().parse_args(
        [
            "--treatment1-bam", str(target),
            "--control1-bam", str(control),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    assert args.condition1_name == "condition1"
    assert args.condition2_name == "condition2"


def test_inspect_run_reports_bams_and_per_sample_fragment_mode(tmp_path: Path, capsys):
    _source_manifest(tmp_path / "run", condition_name="WT")
    assert cutn_suite.inspect_run(tmp_path / "run") == 0
    output = capsys.readouterr().out
    assert "Condition 1: WT" in output
    assert "treatment_1.bam" in output
    assert "control_2.bam" in output
    assert "fragment mode: 152 bp" in output
    assert "fragments used (coverage range): 12" in output
    assert "analysis modes: treatment=152 bp; control=152 bp" in output


def test_rerun_reuses_bigwigs_excludes_sample_and_changes_downstream_parameters(tmp_path: Path):
    root = tmp_path / "run"
    _source_manifest(root)
    commands = []

    def fake_average(paths, output):
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.touch()
        return output

    def fake_command(command):
        commands.append(list(command))
        assert command[0] == "call-peaks"
        prefix = Path(command[command.index("--out-prefix") + 1])
        prefix.parent.mkdir(parents=True, exist_ok=True)
        Path(f"{prefix}_nucleosome_regions.bed").write_text(
            "chr1\t10\t90\tpeak1\t1\t.\t50\t51\n", encoding="utf-8"
        )

    def fake_analyze(*_args, **kwargs):
        output_dir = Path(kwargs["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        clusters = output_dir / "clusters.bed"
        clusters.write_text("", encoding="utf-8")
        return {"selected_clusters": clusters}

    def fake_anchor(_clusters, output):
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("", encoding="utf-8")
        return output

    with (
        patch("nucleosuite.cutn_suite.average_bigwigs", side_effect=fake_average),
        patch("nucleosuite.cutn_suite._run_nucleosuite", side_effect=fake_command),
        patch("nucleosuite.cutn_suite.analyze_cutn_replicate_peaks", side_effect=fake_analyze),
        patch("nucleosuite.cutn_suite.write_cluster_anchor_bed", side_effect=fake_anchor),
    ):
        assert cutn_suite.main(
            [
                "--rerun-from", str(root),
                "--exclude-sample", "treatment_2.bam",
                "--peak-max-neg-run", "2",
                "--cluster-seed-p-value", "0.01",
                "--cluster-member-mode", "significant-only",
                "--min-cluster-members", "1",
                "--skip-cluster-aggregate",
            ]
        ) == 0

    assert commands and all(command[0] == "call-peaks" for command in commands)
    assert not any(command[0] == "tracks" for command in commands)
    assert "--max-neg-run" in commands[0]
    assert commands[0][commands[0].index("--max-neg-run") + 1] == "2"

    rerun = root / "rerun_excluding_treatment_2_01"
    manifest = json.loads((rerun / "cutn_stage1_manifest.json").read_text())
    assert len(manifest["treatment_replicates"]) == 1
    assert len(manifest["control_replicates"]) == 2
    assert manifest["cluster_seed_p_value"] == 0.01
    assert manifest["cluster_member_mode"] == "significant-only"
    assert manifest["minimum_cluster_members"] == 1
    assert manifest["peak_max_neg_run"] == 2
    assert manifest["rerun_excluded_bams"] == [str(root / "treatment_2.bam")]
    run_manifest = json.loads((rerun / "cutn_suite_run_manifest.json").read_text())
    assert run_manifest["source_run"] == str(root.resolve())
    assert run_manifest["changed_parameters"]["peak_max_neg_run"] == 2


def test_rerun_directory_number_increments(tmp_path: Path):
    root = tmp_path / "run"
    root.mkdir()
    (root / "rerun_excluding_sample_01").mkdir()
    (root / "rerun_excluding_sample_02").mkdir()
    assert cutn_suite._next_rerun_directory(root, ["sample.bam"]).name == "rerun_excluding_sample_03"
    assert cutn_suite._next_rerun_directory(root, ["a.bam", "b.bam"]).name == "rerun_excluding_2_samples_01"


def test_rerun_rejects_bigwig_generating_parameter_changes(tmp_path: Path):
    root = tmp_path / "run"
    _source_manifest(root)
    with pytest.raises(ValueError, match="cannot be changed"):
        cutn_suite.main(
            [
                "--rerun-from", str(root),
                "--mode", "170",
            ]
        )
