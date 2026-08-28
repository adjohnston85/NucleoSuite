from pathlib import Path
import pytest

from nucleosuite.cutn_suite import build_parser
from nucleosuite.cutn_peaks import CompetitivePeak, ReplicatePeakStatistics, cluster_seeded_gate_peaks
from nucleosuite.peak_fdr import PeakRow


def _record(i, *, mean_diff, all_gate, p):
    start = i * 100
    row = PeakRow(("chr1", str(start), str(start+80), f"p{i}", "10", ".", str(start+40), str(start+41)), 10.0, i+1)
    peak = CompetitivePeak(row, "chr1", start, start+80, start+40, "target", winner=True, signal_score=100.0)
    return ReplicatePeakStatistics(
        peak=peak, treatment_scores=(10.0, 20.0), control_scores=(12.0, 13.0),
        treatment_mean=15.0, control_mean=12.5, mean_difference=mean_diff,
        minimum_treatment=10.0, maximum_control=13.0,
        conservative_excess=max(10.0-13.0, 0.0), conservative_fold_enrichment=11/14,
        conservative_log2_enrichment=-0.347923, all_controls_gate=all_gate, pvalue=p, qvalue=p,
    )


def test_cutn_suite_public_clustering_defaults(tmp_path: Path):
    t = tmp_path / "t.bam"; c = tmp_path / "c.bam"; t.touch(); c.touch()
    args = build_parser().parse_args(["--treatment1-bam", str(t), "--control1-bam", str(c), "--outdir", str(tmp_path/'out')])
    assert args.cluster_seed_mode == "auto"
    assert args.cluster_seed_gate_mode == "auto"
    assert args.stage1_gate_mode == "auto"
    assert args.stage1_coverage_statistic == "mean"
    assert args.cluster_member_mode == "seed-and-gated"
    assert args.cluster_max_non_member_gap == 1
    assert args.min_cluster_members == 2
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--target-bam", str(t), "--control-bam", str(c), "--outdir", str(tmp_path/'bad')])


def test_mean_gate_can_pass_when_all_controls_gate_fails():
    records = [_record(0, mean_diff=2.5, all_gate=False, p=0.01), _record(1, mean_diff=2.5, all_gate=False, p=0.2)]
    clusters = cluster_seeded_gate_peaks(records, seed_gate_mode="mean", member_gate_mode="mean", member_mode="seed-and-gated")
    assert len(clusters) == 1
    assert clusters[0].score == pytest.approx(5.0)
    assert cluster_seeded_gate_peaks(records, seed_gate_mode="all-controls", member_gate_mode="all-controls") == []


def test_significant_only_and_non_member_gap_split_clusters():
    records = [
        _record(0, mean_diff=2, all_gate=True, p=0.01),  # S
        _record(1, mean_diff=2, all_gate=True, p=0.2),   # G -> non-member
        _record(2, mean_diff=-1, all_gate=False, p=0.2), # x -> second non-member => split
        _record(3, mean_diff=2, all_gate=True, p=0.01),  # S
        _record(4, mean_diff=2, all_gate=True, p=0.01),  # S
    ]
    clusters = cluster_seeded_gate_peaks(
        records, seed_gate_mode="mean", member_gate_mode="mean", member_mode="significant-only",
        maximum_non_member_gap=1, minimum_member_peaks=1,
    )
    assert len(clusters) == 2
    assert [len(c.significant_peaks) for c in clusters] == [1, 2]


def test_cutn_suite_shared_tracks_pass_contains_score_and_broad_coverage(tmp_path: Path):
    from unittest.mock import patch
    from nucleosuite.cutn_suite import _generate_tracks
    from nucleosuite.progress import ProgressReporter

    bam = tmp_path / "sample.bam"
    bam.touch()
    args = build_parser().parse_args([
        "--treatment1-bam", str(bam),
        "--control1-bam", str(bam),
        "--outdir", str(tmp_path / "out"),
        "--sample-name", "cutn",
        "--mode", "152",
        "--contigs", "chr1",
    ])
    tracks_dir = tmp_path / "tracks"
    scaled_dir = tmp_path / "scaled"
    setup_dir = tmp_path / "setup"
    for directory in (tracks_dir, scaled_dir, setup_dir):
        directory.mkdir()

    commands = []

    def fake_run(command):
        commands.append(command)
        assert command[0] == "tracks"
        assert "coverage" not in command[:1]
        spec = Path(command[command.index("--spec-file") + 1])
        rows = [line.split("\t") for line in spec.read_text().splitlines()[1:]]
        assert rows[0][0] == "122-182"
        assert rows[0][2] == "pns,posPNS"
        assert rows[1][0] == "1-1000"
        assert rows[1][2] == "coverage"
        for _range, prefix, track_names, _scope in rows:
            prefix = Path(prefix)
            for track in track_names.split(","):
                Path(f"{prefix}_{track}.bw").touch()

    def fake_scale(_score, _reference, output, **_kwargs):
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.touch()
        return output, 2.0, 10

    with (
        patch("nucleosuite.cutn_suite._run_nucleosuite", side_effect=fake_run),
        patch("nucleosuite.cutn_suite.scale_bigwig_by_reference", side_effect=fake_scale),
    ):
        record = _generate_tracks(
            args=args,
            reporter=ProgressReporter("test"),
            bam_paths=[str(bam)],
            label="target",
            mode=152,
            tracks_dir=tracks_dir,
            analysis_dir=scaled_dir,
            setup_dir=setup_dir,
            score_track="pns",
            positive_track="posPNS",
            seed=123,
        )

    assert len(commands) == 1
    assert record["score_frag_lower"] == 122
    assert record["score_frag_upper"] == 182
    assert record["coverage_frag_lower"] == 1
    assert record["coverage_frag_upper"] == 1000


def test_stage2_manifest_uses_native_pns_score_tracks(tmp_path: Path):
    from nucleosuite.cutn_compare import _manifest_score_tracks

    mean = tmp_path / "mean_pns.bw"
    rep1 = tmp_path / "rep1_pns.bw"
    rep2 = tmp_path / "rep2_pns.bw"
    for path in (mean, rep1, rep2):
        path.touch()
    manifest = {
        "scoring_method": "pns",
        "positive_track": "posPNS",
        "condition_mean_treatment_cluster_aggregate_score": str(mean),
        "treatment_replicates": [
            {"analysis_score": str(rep1)},
            {"analysis_score": str(rep2)},
        ],
    }
    resolved = _manifest_score_tracks(manifest)
    assert resolved is not None
    assert resolved[0] == mean.resolve()
    assert resolved[1] == [rep1.resolve(), rep2.resolve()]
    assert resolved[2:] == ("pns", "posPNS")


def test_stage1_statistics_and_seed_beds_append_raw_p(tmp_path: Path, monkeypatch):
    import nucleosuite.cutn_peaks as cp

    class Handle:
        def __init__(self, score):
            self.score = score
        def close(self):
            pass

    target_bed = tmp_path / "target_candidates.bed"
    target_bed.write_text("chr1\t100\t180\tpeak1\t8\t.\t140\t141\n", encoding="utf-8")
    handles = [Handle(120.0), Handle(122.0), Handle(50.0), Handle(52.0), Handle(121.0)]
    monkeypatch.setattr(cp, "open_bigwigs", lambda _paths: handles)
    monkeypatch.setattr(cp, "interval_mean", lambda handle, *_args: handle.score)

    outputs = cp.analyze_cutn_replicate_peaks(
        target_bed,
        output_dir=tmp_path / "out",
        target_replicate_bigwigs=["t1.bw", "t2.bw"],
        control_replicate_bigwigs=["c1.bw", "c2.bw"],
        target_mean_bigwig="treatment_mean.bw",
        minimum_cluster_members=1,
    )

    assert outputs["annotated_peaks"].name == "target_peaks_replicate_statistics.bed"
    assert outputs["seed_peaks"].name == "target_seed_peaks_S-pvalue-mean.bed"
    assert outputs["competition_table"].name == "target_peak_replicate_statistics.tsv"
    annotated = outputs["annotated_peaks"].read_text().strip().split("\t")
    seed = outputs["seed_peaks"].read_text().strip().split("\t")
    assert len(annotated) == 9
    assert seed == annotated
    assert float(annotated[-1]) < 0.05
