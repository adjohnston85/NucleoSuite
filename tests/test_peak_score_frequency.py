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
    prefix = Path(f"{prefix}_bins4_scoreminnone_scoremaxnone_normcount")

    with open(f"{prefix}_score_frequency.tsv", "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    totals = {}
    for row in rows:
        totals[row["dataset"]] = totals.get(row["dataset"], 0) + int(row["count"])
    assert totals == {"observed": 5, "randomized": 5}
    assert Path(f"{prefix}_scores.tsv.gz").is_file()
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
    prefix = Path(f"{prefix}_binsinteger_scoreminnone_scoremaxnone_normcount")

    with open(f"{prefix}_score_frequency.tsv", "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert list(rows[0]) == [
        "dataset", "score", "count", "fraction", "percent",
        "cumulative_count", "cumulative_fraction", "cumulative_percent",
    ]
    observed_rows = {int(row["score"]): int(row["count"]) for row in rows}
    assert observed_rows == {0: 1, 1: 2, 2: 1, 3: 0, 4: 1}
    assert int(rows[-1]["cumulative_count"]) == 5
