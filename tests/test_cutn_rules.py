from pathlib import Path

from nucleosuite.cutn_compare import ClusterRecord, _consensus_cluster_regions
from nucleosuite.cutn_suite import _resolve_cluster_rules, build_parser


def _args(tmp_path: Path):
    t = tmp_path / "t.bam"; c = tmp_path / "c.bam"
    t.touch(); c.touch()
    args = build_parser().parse_args([
        "--treatment1-bam", str(t), "--control1-bam", str(c),
        "--outdir", str(tmp_path / "out"),
    ])
    args._explicit_options = set()
    return args


def test_automatic_cluster_rules_use_gate_only_when_either_group_has_fewer_than_three(tmp_path, capsys):
    args = _args(tmp_path)
    rules = _resolve_cluster_rules(
        args,
        [Path("t1"), Path("t2")],
        [Path("c1"), Path("c2"), Path("c3")],
    )
    assert rules == {
        "seed_mode": "gated",
        "seed_gate_mode": "all-controls",
        "member_gate_mode": "all-controls",
        "compute_pvalues": False,
    }
    text = capsys.readouterr().err
    assert "One or both groups have fewer than 3 replicates" in text
    assert "Peak p-values are not used because replicate counts are insufficient" in text


def test_automatic_cluster_rules_use_pvalue_mean_seed_at_three_by_three(tmp_path, capsys):
    args = _args(tmp_path)
    rules = _resolve_cluster_rules(
        args,
        [Path("t1"), Path("t2"), Path("t3")],
        [Path("c1"), Path("c2"), Path("c3")],
    )
    assert rules == {
        "seed_mode": "pvalue",
        "seed_gate_mode": "mean",
        "member_gate_mode": "all-controls",
        "compute_pvalues": True,
    }
    text = capsys.readouterr().err
    assert "raw p < 0.05 AND mean treatment > mean control" in text
    assert "all treatment replicates > all control replicates" in text


def test_overlapping_cluster_locus_measures_actual_intersection():
    first = [ClusterRecord("c1", "chr1", 100, 220, 160, 10.0, 1)]
    second = [ClusterRecord("c2", "chr1", 180, 300, 240, 9.0, 2)]
    regions = _consensus_cluster_regions(first, second, chroms={"chr1": 1000})
    assert len(regions) == 1
    assert (regions[0].start, regions[0].end) == (100, 300)
    assert regions[0].measurement_intervals == ((180, 220),)
