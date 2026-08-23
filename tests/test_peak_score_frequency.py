from pathlib import Path

import csv

from nucleosuite.peak_score_frequency import main


def _write_bed(path: Path, scores: list[float]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        for index, score in enumerate(scores):
            start = index * 20
            handle.write(f"chr1\t{start}\t{start + 10}\tpeak{index}\t{score}\t.\n")


def test_peak_score_frequency_overlays_observed_and_randomized(tmp_path: Path):
    observed = tmp_path / "observed.bed"
    randomized = tmp_path / "randomized.bed"
    _write_bed(observed, [1, 2, 3, 4, 5])
    _write_bed(randomized, [1, 1, 2, 2, 3])
    prefix = tmp_path / "scores"

    assert main([
        "--peaks", f"observed={observed}",
        "--peaks", f"randomized={randomized}",
        "--output-prefix", str(prefix),
        "--bins", "4",
    ]) == 0
    prefix = Path(f"{prefix}_bins4_scorescale1_scoreminnone")

    with open(f"{prefix}_score_frequency.tsv", "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    totals = {}
    for row in rows:
        totals[row["dataset"]] = totals.get(row["dataset"], 0) + int(row["count"])
    assert totals == {"observed": 5, "randomized": 5}
    assert not Path(f"{prefix}_scores.tsv.gz").exists()
    assert Path(f"{prefix}_score_summary.tsv").is_file()
    assert Path(f"{prefix}_score_frequency.png").is_file()


def test_peak_score_frequency_defaults_to_integer_score_bins(tmp_path: Path):
    observed = tmp_path / "observed_integer.bed"
    _write_bed(observed, [0.2, 0.5, 1.49, 1.5, 3.6])
    prefix = tmp_path / "integer_scores"

    assert main([
        "--peaks", f"observed={observed}",
        "--output-prefix", str(prefix),
    ]) == 0
    prefix = Path(f"{prefix}_binsinteger_scorescale1_scoreminnone")

    with open(f"{prefix}_score_frequency.tsv", "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert list(rows[0]) == [
        "dataset", "score", "count", "fraction", "percent",
        "cumulative_count", "cumulative_fraction", "cumulative_percent",
    ]
    observed_rows = {int(row["score"]): int(row["count"]) for row in rows if int(row["count"]) > 0}
    assert observed_rows == {0: 1, 1: 2, 2: 1, 4: 1}
    assert int(rows[-1]["cumulative_count"]) == 5


def test_peak_score_frequency_individual_scores_are_opt_in(tmp_path: Path):
    observed = tmp_path / "observed_detail.bed"
    _write_bed(observed, [1, 2, 3])
    prefix = tmp_path / "detail_scores"
    assert main([
        "--peaks", f"observed={observed}",
        "--output-prefix", str(prefix),
        "--write-detail-tables",
    ]) == 0
    prefix = Path(f"{prefix}_binsinteger_scorescale1_scoreminnone")
    assert Path(f"{prefix}_scores.tsv.gz").is_file()


def test_peak_score_frequency_score_scale_applies_before_integer_binning(tmp_path: Path):
    observed = tmp_path / "pns_scores.bed"
    _write_bed(observed, [0.004, 0.006, 0.014, 0.016])
    prefix = tmp_path / "scaled_scores"

    assert main([
        "--peaks", f"PNS={observed}",
        "--output-prefix", str(prefix),
        "--score-scale", "100",
    ]) == 0
    prefix = Path(f"{prefix}_binsinteger_scorescale100_scoreminnone")
    with open(f"{prefix}_score_frequency.tsv", "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    counts = {int(row["score"]): int(row["count"]) for row in rows}
    assert counts == {0: 1, 1: 2, 2: 1}

    # Raw-score summaries remain in the units of the supplied BED scores.
    with open(f"{prefix}_score_summary.tsv", "rt", encoding="utf-8") as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert float(summary["maximum"]) == 0.016

    metadata = Path(f"{prefix}_score_frequency_metadata.tsv").read_text()
    assert "score_scale\t100.0" in metadata
    assert "x_label\tPeak score ×100" in metadata


def test_peak_score_frequency_default_scale_is_one_for_all_input_formats(tmp_path: Path):
    import numpy as np
    from nucleosuite.peak_score_frequency import PeakScoreSet, _effective_score_scale
    def ds(path):
        return PeakScoreSet("test", path, np.array([], dtype=float), [])
    assert _effective_score_scale(ds(tmp_path / "peaks.bed"), None) == 1.0
    assert _effective_score_scale(ds(tmp_path / "peaks.bed.gz"), None) == 1.0
    assert _effective_score_scale(ds(tmp_path / "peaks.bb"), None) == 1.0
    assert _effective_score_scale(ds(tmp_path / "peaks.bigBed"), 5.0) == 5.0


def test_peak_score_frequency_default_output_prefix_uses_input_basename(tmp_path: Path, monkeypatch):
    observed = tmp_path / "sample_nucleosome_regions.bed"
    _write_bed(observed, [0.1, 0.2, 0.3])
    monkeypatch.chdir(tmp_path)
    assert main(["--peaks", str(observed)]) == 0
    expected = tmp_path / "sample_nucleosome_regions_peak_score_frequency_binsinteger_scorescale1_scoreminnone_score_frequency.tsv"
    assert expected.is_file()
