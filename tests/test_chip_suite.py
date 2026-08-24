import json
from pathlib import Path
from unittest.mock import patch

import pyBigWig
import pytest

from nucleosuite.chip_compare import compare_stage1
from nucleosuite.chip_peaks import (
    analyze_chip_peaks,
    analyze_chip_replicate_peaks,
    assign_competition_qvalues,
    compete_peaks,
    compete_peaks_all_controls,
    compete_peaks_with_bigwigs,
)
from nucleosuite.chip_suite import (
    _locate_completed_prefix,
    _resolved_prefix,
    _validate,
    build_parser,
    run,
)
from nucleosuite.peak_fdr import PeakRow


def _row(start: int, score: float, line: int = 1) -> PeakRow:
    centre = start + 40
    return PeakRow(
        ("chr1", str(start), str(start + 80), f"peak{line}", str(score), ".", str(centre), str(centre + 1)),
        score,
        line,
    )


def test_chip_suite_defaults_to_tns_auto_mode_and_120_500(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--target-bam", str(target),
            "--control-bam", str(control),
            "--outdir", str(tmp_path / "out"),
        ]
    )
    assert args.scoring_method == "tns"
    assert args.mode == "auto"
    assert (args.frag_lower, args.frag_upper) == (120, 500)
    assert args.mode_strategy == "pooled"
    assert args.bam_mode == "replicates"
    assert args.stage1_control_mode == "all-controls"
    assert args.treatment1_bam == [str(target)]
    assert args.control1_bam == [str(control)]


def test_chip_suite_emits_target_and_control_raw_coverage(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--target-bam", str(target),
            "--control-bam", str(control),
            "--outdir", str(tmp_path / "out"),
            "--sample-name", "sample",
            "--mode", "167",
            "--contigs", "chr1",
        ]
    )
    pns_commands: list[list[str]] = []

    def fake_command(command: list[str]) -> None:
        if command[0] == "pns":
            pns_commands.append(command)
            base = Path(command[command.index("--out-prefix") + 1])
            prefix = _resolved_prefix(base, "tns", 167, 120, 500)
            prefix.parent.mkdir(parents=True, exist_ok=True)
            for suffix in ("_tns.bw", "_posTNS.bw", "_coverage.bw"):
                Path(f"{prefix}{suffix}").touch()
            return
        peak_prefix = Path(command[command.index("--out-prefix") + 1])
        peak_prefix.parent.mkdir(parents=True, exist_ok=True)
        Path(f"{peak_prefix}_nucleosome_regions.bed").write_text(
            "chr1\t10\t20\tpeak1\t1\t.\t15\t16\n", encoding="utf-8"
        )

    def fake_scale(_score, _reference, output, **_kwargs):
        Path(output).touch()
        return Path(output), 1.0, 1

    with (
        patch("nucleosuite.chip_suite._run_nucleosuite", side_effect=fake_command),
        patch("nucleosuite.chip_suite.scale_bigwig_by_reference", side_effect=fake_scale),
            patch(
                "nucleosuite.chip_suite.analyze_chip_replicate_peaks",
                return_value={},
            ) as analyze_mock,
    ):
        assert run(args) == 0

    assert len(pns_commands) == 2
    assert all(
        command[command.index("--other-tracks") + 1] == "coverage"
        and command[command.index("--other-format") + 1] == "bigwig"
        for command in pns_commands
    )
    summary = (tmp_path / "out" / "sample_chip_suite_summary.tsv").read_text()
    assert "target_raw_coverage_track\t" in summary
    assert "control_raw_coverage_track\t" in summary
    assert "condition_mean_treatment_coverage\t" in summary
    assert "condition_mean_control_coverage\t" in summary
    call = analyze_mock.call_args.kwargs
    assert len(call["target_replicate_bigwigs"]) == 1
    assert len(call["control_replicate_bigwigs"]) == 1
    manifest = json.loads(
        (tmp_path / "out" / "chip_stage1_manifest.json").read_text()
    )
    assert manifest["control_candidate_peaks"] is None
    assert manifest["stage1_statistics"] == "one_sided_welch_bh_all_candidates"


def test_explicit_chip_mode_overrides_automatic_estimation(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--target-bam", str(target),
            "--control-bam", str(control),
            "--outdir", str(tmp_path / "out"),
            "--mode", "167",
        ]
    )
    assert args.mode == 167


def test_explicit_chip_mode_does_not_validate_unused_auto_search_bounds(tmp_path: Path):
    target = tmp_path / "target.bam"
    control = tmp_path / "control.bam"
    target.touch(); control.touch()
    args = build_parser().parse_args(
        [
            "--target-bam", str(target),
            "--control-bam", str(control),
            "--outdir", str(tmp_path / "out"),
            "--mode", "167",
            "--frag-upper", "200",
        ]
    )
    _validate(args)


def test_chip_suite_locates_parameterized_multicontig_score_outputs(tmp_path: Path):
    requested = tmp_path / "tracks" / "chip_target"
    direct = _resolved_prefix(requested, "tns", 153, 120, 500)
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
    for suffix in ("_tns.bw", "_posTNS.bw", "_coverage.bw"):
        (combined / f"{combined_name}{suffix}").touch()

    located = _locate_completed_prefix(
        requested, direct, ("_tns.bw", "_posTNS.bw", "_coverage.bw")
    )

    assert located == combined / combined_name


def test_chip_suite_locates_multicontig_peak_outputs(tmp_path: Path):
    requested = tmp_path / "peaks" / "chip_target_tns_mean_scaled"
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


def test_chip_suite_prefers_existing_serial_outputs(tmp_path: Path):
    requested = tmp_path / "tracks" / "chip_target"
    direct = _resolved_prefix(requested, "tns", 153, 120, 500)
    direct.parent.mkdir(parents=True)
    for suffix in ("_tns.bw", "_posTNS.bw", "_coverage.bw"):
        Path(f"{direct}{suffix}").touch()

    located = _locate_completed_prefix(
        requested, direct, ("_tns.bw", "_posTNS.bw", "_coverage.bw")
    )

    assert located == direct


def test_target_control_peak_competition_uses_control_for_ties():
    target, control = compete_peaks(
        [_row(100, 8), _row(300, 10, 2)],
        [_row(102, 8), _row(700, 9, 2)],
        match_distance=84,
    )
    assert target[0].winner is False
    assert control[0].winner is True
    assert target[1].winner is True
    assert control[1].winner is True
    annotated, _ = assign_competition_qvalues(target, control)
    assert annotated[0].qvalue == 1.0
    assert 0 <= annotated[1].qvalue <= 1


def test_chip_peak_outputs_preserve_target_bed_and_append_fdr(tmp_path: Path):
    target = tmp_path / "target.bed"
    control = tmp_path / "control.bed"
    target.write_text(
        "chr1\t100\t180\tt1\t20\t.\t140\t141\n"
        "chr1\t280\t360\tt2\t18\t.\t320\t321\n",
        encoding="utf-8",
    )
    control.write_text(
        "chr1\t500\t580\tc1\t5\t.\t540\t541\n",
        encoding="utf-8",
    )
    outputs = analyze_chip_peaks(
        target, control, output_dir=tmp_path / "out", peak_fdr=1.0, cluster_fdr=1.0
    )
    rows = outputs["annotated_peaks"].read_text().splitlines()
    assert len(rows) == 2
    assert all(len(row.split("\t")) == 9 for row in rows)
    assert outputs["cluster_table"].is_file()


def _write_bigwig(path: Path, values: list[float]) -> None:
    handle = pyBigWig.open(str(path), "w")
    handle.addHeader([("chr1", len(values))])
    handle.addEntries("chr1", 0, values=values, span=1, step=1)
    handle.close()


def test_stage1_uses_control_bigwig_max_inside_target_interval(tmp_path: Path):
    target_bw = tmp_path / "target.bw"
    control_bw = tmp_path / "control.bw"
    target_values = [0.0] * 300
    control_values = [0.0] * 300
    target_values[100:180] = [8.0] * 80
    control_values[150] = 10.0
    _write_bigwig(target_bw, target_values)
    _write_bigwig(control_bw, control_values)
    target, _control = compete_peaks_with_bigwigs(
        [_row(100, 8.0)],
        [],
        target_bigwig=target_bw,
        control_bigwig=control_bw,
    )
    assert target[0].winner is False
    assert target[0].matched_score == 10.0
    assert target[0].competition_score == 0.0


def test_stage1_bed_score_is_scaled_coverage_max(tmp_path: Path):
    target_bed = tmp_path / "target.bed"
    control_bed = tmp_path / "control.bed"
    target_bed.write_text(
        "chr1\t100\t180\tt1\t8\t.\t140\t141\n", encoding="utf-8"
    )
    control_bed.write_text("", encoding="utf-8")
    target_bw = tmp_path / "target_coverage_scaled.bw"
    control_bw = tmp_path / "control_coverage_scaled.bw"
    target_values = [0.0] * 300
    control_values = [0.0] * 300
    target_values[150] = 275.0
    control_values[150] = 25.0
    _write_bigwig(target_bw, target_values)
    _write_bigwig(control_bw, control_values)
    outputs = analyze_chip_peaks(
        target_bed,
        control_bed,
        output_dir=tmp_path / "out",
        peak_fdr=1.0,
        cluster_fdr=1.0,
        target_bigwig=target_bw,
        control_bigwig=control_bw,
    )
    fields = outputs["annotated_peaks"].read_text().strip().split("\t")
    assert fields[4] == "275"
    assert len(fields) == 9


def test_stage1_all_controls_requires_every_treatment_to_exceed_every_control(
    tmp_path: Path,
):
    treatment_paths = [tmp_path / "treatment1.bw", tmp_path / "treatment2.bw"]
    control_paths = [tmp_path / "control1.bw", tmp_path / "control2.bw"]
    treatment_values = [[0.0] * 300 for _ in treatment_paths]
    control_values = [[0.0] * 300 for _ in control_paths]
    treatment_values[0][100:180] = [140.0] * 80
    treatment_values[1][100:180] = [90.0] * 80
    control_values[0][100:180] = [80.0] * 80
    control_values[1][100:180] = [100.0] * 80
    for path, values in zip(treatment_paths, treatment_values):
        _write_bigwig(path, values)
    for path, values in zip(control_paths, control_values):
        _write_bigwig(path, values)
    treatment_mean = tmp_path / "treatment_mean.bw"
    control_mean = tmp_path / "control_mean.bw"
    _write_bigwig(treatment_mean, [115.0 if 100 <= i < 180 else 0.0 for i in range(300)])
    _write_bigwig(control_mean, [90.0 if 100 <= i < 180 else 0.0 for i in range(300)])

    target, _control = compete_peaks_all_controls(
        [_row(100, 8.0)],
        [],
        target_bigwigs=treatment_paths,
        control_bigwigs=control_paths,
        target_mean_bigwig=treatment_mean,
        control_mean_bigwig=control_mean,
    )

    assert target[0].winner is False
    assert target[0].signal_score == 115.0
    assert target[0].matched_score == 100.0
    assert target[0].competition_score == 0.0
    assert target[0].treatment_replicate_scores == (140.0, 90.0)
    assert target[0].control_replicate_scores == (80.0, 100.0)


def test_stage1_all_controls_accepts_peak_when_minimum_treatment_exceeds_maximum_control(
    tmp_path: Path,
):
    paths = [tmp_path / f"track{index}.bw" for index in range(6)]
    scores = [140.0, 110.0, 80.0, 100.0, 125.0, 90.0]
    for path, score in zip(paths, scores):
        _write_bigwig(
            path,
            [score if 100 <= index < 180 else 0.0 for index in range(300)],
        )

    target, _control = compete_peaks_all_controls(
        [_row(100, 8.0)],
        [],
        target_bigwigs=paths[:2],
        control_bigwigs=paths[2:4],
        target_mean_bigwig=paths[4],
        control_mean_bigwig=paths[5],
    )

    assert target[0].winner is True
    assert target[0].competition_score == 10.0


def test_stage1_replicate_statistics_calls_fdr_without_control_peaks(tmp_path: Path):
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

    outputs = analyze_chip_replicate_peaks(
        target_bed,
        output_dir=tmp_path / "results",
        target_replicate_bigwigs=treatment_paths,
        control_replicate_bigwigs=control_paths,
        target_mean_bigwig=mean_path,
        peak_fdr=0.05,
        cluster_fdr=0.05,
        minimum_significant_peaks=1,
    )

    annotated = outputs["annotated_peaks"].read_text().strip().split("\t")
    significant = outputs["significant_peaks"].read_text().strip().split("\t")
    statistics = outputs["competition_table"].read_text().splitlines()
    assert annotated[4] == "120"
    assert annotated[-1] == "0"
    assert significant == annotated
    assert "treatment_replicate_maxima" in statistics[0]
    assert "control_replicate_maxima" in statistics[0]
    assert statistics[1].split("\t")[-3:] == ["true", "0", "0"]


def test_stage1_single_replicate_reports_unavailable_fdr(tmp_path: Path):
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

    outputs = analyze_chip_replicate_peaks(
        target_bed,
        output_dir=tmp_path / "results",
        target_replicate_bigwigs=[treatment],
        control_replicate_bigwigs=[control],
        target_mean_bigwig=mean_path,
        minimum_significant_peaks=1,
    )

    assert outputs["annotated_peaks"].read_text().strip().endswith("\t.")
    assert outputs["significant_peaks"].read_text() == ""


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
    assert "stage2-four-group-interaction-bh" in output


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
                "scaled_score": str(target_path),
                "scaled_coverage": str(target_path),
            }
        )
    for index, control in enumerate(control_values, 1):
        control_path = root / f"control_{index}.bw"
        _write_bigwig(control_path, control)
        control_records.append(
            {
                "replicate": index,
                "scaled_score": str(control_path),
                "scaled_coverage": str(control_path),
            }
        )
    manifest = {
        "schema": "nucleosuite_chip_stage1",
        "schema_version": 2,
        "condition_name": name,
        "bam_mode": "replicates",
        "scoring_method": "tns",
        "score_track": "tns",
        "positive_track": "posTNS",
        "peak_discovery_track": "tns",
        "peak_measurement_track": "coverage_divided_by_nonzero_mean_x100",
        "target_mode": 167,
        "control_mode": 167,
        "frag_lower": 120,
        "frag_upper": 500,
        "contigs": ["chr1"],
        "treatment_replicates": treatment_records,
        "control_replicates": control_records,
        "condition_mean_treatment_score": treatment_records[0]["scaled_score"],
        "condition_mean_control_score": control_records[0]["scaled_score"],
        "condition_mean_treatment_coverage": treatment_records[0]["scaled_coverage"],
        "condition_mean_control_coverage": control_records[0]["scaled_coverage"],
        "significant_peaks": str(significant_peaks),
        "significant_clusters": str(significant_clusters),
    }
    path = root / "chip_stage1_manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def test_chip_compare_uses_scaled_bigwig_replicates_without_bams(tmp_path: Path):
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
    assert payload["results"]["peaks"]["inferential_fdr_available"] is True
    row = (tmp_path / "comparison" / "differential_peaks.tsv").read_text().splitlines()[1]
    fields = row.split("\t")
    assert fields[-1] == "significant_gain"
    header = (tmp_path / "comparison" / "differential_peaks.tsv").read_text().splitlines()[0].split("\t")
    assert fields[header.index("region_origin")] == "overlap_union"
    gain_fields = (
        tmp_path / "comparison" / "differential_peaks_fdr0.05_gains.bed"
    ).read_text().strip().split("\t")
    assert gain_fields[-1] == "overlap_union"
    cluster_lines = (
        tmp_path / "comparison" / "differential_clusters.tsv"
    ).read_text().splitlines()
    cluster_header = cluster_lines[0].split("\t")
    cluster_row = cluster_lines[1].split("\t")
    assert cluster_row[cluster_header.index("region_origin")] == "overlap_union"


def test_chip_compare_uses_unpaired_unequal_four_group_interaction(tmp_path: Path):
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
    assert payload["results"]["peaks"]["inferential_fdr_available"] is True
    lines = (tmp_path / "comparison" / "differential_peaks.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    assert "condition1_treatment_replicate_3" in header
    assert "condition1_control_replicate_3" not in header
    assert "condition2_control_replicate_3" in header
    assert row[header.index("interaction_difference")] == "18"
    assert row[-1] == "significant_gain"


def test_chip_compare_uses_union_for_overlapping_peak_coordinates(tmp_path: Path):
    length = 300
    values = [[0.0] * length for _ in range(2)]
    first = _stage1_manifest(tmp_path / "first", "first", values, values)
    second = _stage1_manifest(tmp_path / "second", "second", values, values)
    (tmp_path / "first" / "peaks.bed").write_text(
        "chr1\t80\t150\tpeak1\t10\t.\t120\t121\t0.01\n", encoding="utf-8"
    )
    (tmp_path / "second" / "peaks.bed").write_text(
        "chr1\t120\t200\tpeak2\t10\t.\t160\t161\t0.01\n", encoding="utf-8"
    )

    compare_stage1(first, second, outdir=tmp_path / "comparison", feature_level="peaks")

    row = (tmp_path / "comparison" / "differential_peaks.tsv").read_text().splitlines()[1]
    fields = row.split("\t")
    assert fields[:3] == ["chr1", "80", "200"]
    assert fields[4:6] == ["true", "true"]
    assert fields[6] == "overlap_union"


def test_chip_compare_retains_condition_specific_peak_regions(tmp_path: Path):
    length = 300
    values = [[0.0] * length for _ in range(2)]
    first = _stage1_manifest(tmp_path / "first", "first", values, values)
    second = _stage1_manifest(tmp_path / "second", "second", values, values)
    (tmp_path / "first" / "peaks.bed").write_text(
        "chr1\t80\t120\tpeak1\t10\t.\t100\t101\t0.01\n", encoding="utf-8"
    )
    (tmp_path / "second" / "peaks.bed").write_text(
        "chr1\t180\t220\tpeak2\t10\t.\t200\t201\t0.01\n", encoding="utf-8"
    )

    compare_stage1(first, second, outdir=tmp_path / "comparison", feature_level="peaks")

    rows = (tmp_path / "comparison" / "differential_peaks.tsv").read_text().splitlines()[1:]
    coordinates = [row.split("\t")[:7] for row in rows]
    assert coordinates == [
        ["chr1", "80", "120", "100", "true", "false", "condition1_only"],
        ["chr1", "180", "220", "200", "false", "true", "condition2_only"],
    ]


def test_chip_compare_labels_proximity_unions(tmp_path: Path):
    length = 300
    values = [[0.0] * length for _ in range(2)]
    first = _stage1_manifest(tmp_path / "first", "first", values, values)
    second = _stage1_manifest(tmp_path / "second", "second", values, values)
    (tmp_path / "first" / "peaks.bed").write_text(
        "chr1\t80\t120\tpeak1\t10\t.\t100\t101\t0.01\n", encoding="utf-8"
    )
    (tmp_path / "second" / "peaks.bed").write_text(
        "chr1\t130\t170\tpeak2\t10\t.\t140\t141\t0.01\n", encoding="utf-8"
    )

    manifest = compare_stage1(
        first,
        second,
        outdir=tmp_path / "comparison",
        feature_level="peaks",
        peak_match_distance=50,
    )

    lines = (tmp_path / "comparison" / "differential_peaks.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    assert row[:3] == ["chr1", "80", "170"]
    assert row[header.index("region_origin")] == "proximity_union"
    payload = json.loads(manifest.read_text())
    assert payload["results"]["peaks"]["region_origin_counts"] == {
        "proximity_union": 1
    }
