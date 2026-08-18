from __future__ import annotations

import gzip
import sys
import types
from pathlib import Path

import numpy as np

from nucleosuite.positive_runs import PositiveRunScanner, RunCollector, SelectedRegion, run_analysis, select_regions


def test_scanner_preserves_runs_across_chunks_and_breaks_on_missing(tmp_path: Path):
    output = tmp_path / "runs.tsv.gz"
    with gzip.open(output, "wt") as handle:
        handle.write("header\n")
        collector = RunCollector(handle, min_run_length=1, max_run_length=None)
        scanner = PositiveRunScanner(collector, threshold=0.0)
        scanner.consume("chr1", 0, np.array([-1.0, 1.0, 2.0, 3.0]))
        scanner.consume("chr1", 4, np.array([4.0, np.nan, 5.0, 6.0, 0.0]))
        scanner.finish_region()

    assert collector.counts == {4: 1, 2: 1}
    with gzip.open(output, "rt") as handle:
        rows = handle.read().splitlines()[1:]
    assert rows[0].split("\t")[:4] == ["chr1", "1", "5", "4"]
    assert rows[1].split("\t")[:4] == ["chr1", "6", "8", "2"]
    assert collector.summary.scanned_bases == 9
    assert collector.summary.missing_bases == 1
    assert collector.summary.positive_bases == 6


def test_run_analysis_writes_plot_data_and_summary(tmp_path: Path, monkeypatch):
    input_bw = tmp_path / "signal.bw"
    input_bw.write_bytes(b"placeholder")

    class FakeBigWig:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return None

        def chroms(self):
            return {"chr1": 10}

        def values(self, chrom, start, end, numpy=False):
            data = np.array([-1, 1, 2, 0, 3, 4, 5, -1, 2, 0], dtype=float)
            result = data[start:end]
            return result if numpy else result.tolist()

    fake_module = types.SimpleNamespace(open=lambda path: FakeBigWig())
    monkeypatch.setitem(sys.modules, "pyBigWig", fake_module)

    outputs = run_analysis(
        input_bw,
        tmp_path / "positive",
        chunk_size=3,
        plot_x_max=20,
    )
    assert all(path.exists() for path in outputs.values())
    counts = outputs["counts"].read_text().splitlines()
    assert counts[0] == "run_length_bp\tcount\tfraction\tpercent"
    assert "1\t1\t" in counts[1]
    assert any(line.startswith("2\t1\t") for line in counts[1:])
    assert any(line.startswith("3\t1\t") for line in counts[1:])
    summary = outputs["summary"].read_text()
    assert "total_runs_retained\t3" in summary
    with gzip.open(outputs["runs"], "rt") as handle:
        rows = handle.read().splitlines()
    assert len(rows) == 4


def test_select_regions_resolves_chr_alias():
    regions = select_regions(["20"], {"chr20": 1000})
    assert regions == [SelectedRegion("chr20", 0, 1000)]
