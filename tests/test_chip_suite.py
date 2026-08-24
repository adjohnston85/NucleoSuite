from pathlib import Path

from nucleosuite.chip_peaks import (
    analyze_chip_peaks,
    assign_competition_qvalues,
    compete_peaks,
)
from nucleosuite.chip_suite import _validate, build_parser
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
