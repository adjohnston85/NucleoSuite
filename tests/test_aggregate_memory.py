"""Bounded-memory aggregate and multicontig accumulator tests."""

from __future__ import annotations

import json
import gzip
from pathlib import Path

import numpy as np
import pytest

from nucleosuite.aggregate_parallel import _read_accumulator
from nucleosuite.align import write_aggregate_accumulator


def _read_profile(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter="\t", skiprows=1, usecols=1)


def test_aggregate_accumulator_round_trip_preserves_exact_totals(tmp_path: Path):
    path = tmp_path / "part.tsv.gz"
    totals = np.array([4.5, -2.0, 0.0], dtype=float)
    counts = np.array([3, 2, 0], dtype=np.int64)
    write_aggregate_accumulator(path, totals, counts, valid_total=3)

    positions, observed_totals, observed_counts, valid_total = _read_accumulator(path)
    assert np.array_equal(positions, np.array([-1, 0, 1]))
    assert np.array_equal(observed_totals, totals)
    assert np.array_equal(observed_counts, counts)
    assert valid_total == 3


def test_multicontig_mean_only_uses_compact_accumulators(tmp_path: Path):
    pybigwig = pytest.importorskip("pyBigWig")
    from nucleosuite.cli.aggregate import main as aggregate_main

    bigwig = tmp_path / "signal.bw"
    handle = pybigwig.open(str(bigwig), "w")
    handle.addHeader([("chr1", 100), ("chr2", 100)])
    handle.addEntries(
        ["chr1"] * 100 + ["chr2"] * 100,
        list(range(100)) * 2,
        ends=list(range(1, 101)) * 2,
        values=[float(value) for value in range(100)]
        + [float(value + 100) for value in range(100)],
    )
    handle.close()
    regions = tmp_path / "regions.bed"
    regions.write_text(
        "chr1\t19\t20\ta\t0\t+\n"
        "chr1\t39\t40\tb\t0\t+\n"
        "chr2\t19\t20\tc\t0\t+\n"
        "chr2\t39\t40\td\t0\t+\n"
    )

    serial_dir = tmp_path / "serial"
    parallel_dir = tmp_path / "parallel"
    common = [
        "--bigwig", str(bigwig),
        "--region-bed", str(regions),
        "--window-half", "1",
        "--zero-thresh", "0",
        "--max-score", "inf",
        "--no-nrl",
        "--output-prefix", "sample",
    ]
    assert aggregate_main([*common, "--cores", "1", "--output-dir", str(serial_dir)]) == 0
    assert aggregate_main([*common, "--cores", "2", "--output-dir", str(parallel_dir)]) == 0

    serial_profile = next(serial_dir.glob("*_aggregate_all.tsv"))
    root = next(parallel_dir.glob("*_multicontig"))
    parallel_profile = next((root / "combined").glob("*_aggregate_all.tsv"))
    assert np.array_equal(_read_profile(serial_profile), _read_profile(parallel_profile))

    manifest = json.loads(
        (root / "nucleosuite_multicontig_manifest.json").read_text()
    )
    assert manifest["combine_strategy"] == "aggregate_accumulator"
    assert len(manifest["per_contig"]) == 2
    for entry in manifest["per_contig"]:
        assert Path(entry["accumulator"]).is_file()
        prefix = Path(entry["prefix"])
        assert not (prefix.parent / f"{prefix.name}_heatmap_matrix.tsv.gz").exists()

    detail_dir = tmp_path / "detail"
    assert aggregate_main(
        [
            *common,
            "--cores", "2",
            "--output-dir", str(detail_dir),
            "--write-detail-tables",
            "--max-heatmap-rows", "1",
            "--subsample-mode", "random",
            "--seed", "7",
        ]
    ) == 0
    detail_root = next(detail_dir.glob("*_multicontig"))
    detail_manifest = json.loads(
        (detail_root / "nucleosuite_multicontig_manifest.json").read_text()
    )
    assert detail_manifest["combine_strategy"] == "aggregate_matrix"
    for entry in detail_manifest["per_contig"]:
        prefix = Path(entry["prefix"])
        matrix = prefix.parent / f"{prefix.name}_heatmap_matrix.tsv.gz"
        with gzip.open(matrix, "rt", encoding="utf-8") as handle:
            assert sum(1 for _line in handle) == 2
    combined_matrix = next(
        (detail_root / "combined").glob("*_heatmap_matrix.tsv.gz")
    )
    with gzip.open(combined_matrix, "rt", encoding="utf-8") as handle:
        assert sum(1 for _line in handle) == 2
    detail_profile = next(
        (detail_root / "combined").glob("*_aggregate_all.tsv")
    )
    assert np.array_equal(_read_profile(serial_profile), _read_profile(detail_profile))

    stopped_serial_dir = tmp_path / "stopped_serial"
    stopped_parallel_dir = tmp_path / "stopped_parallel"
    stopped = [*common, "--stop-after-valid", "3"]
    assert aggregate_main(
        [*stopped, "--cores", "1", "--output-dir", str(stopped_serial_dir)]
    ) == 0
    assert aggregate_main(
        [*stopped, "--cores", "2", "--output-dir", str(stopped_parallel_dir)]
    ) == 0
    stopped_serial_profile = next(stopped_serial_dir.glob("*_aggregate_all.tsv"))
    stopped_root = next(stopped_parallel_dir.glob("*_multicontig"))
    stopped_parallel_profile = next(
        (stopped_root / "combined").glob("*_aggregate_all.tsv")
    )
    assert np.array_equal(
        _read_profile(stopped_serial_profile),
        _read_profile(stopped_parallel_profile),
    )
