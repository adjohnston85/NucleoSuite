"""Fragment-size NRL calculation and output tests."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from nucleosuite.fragment_lengths import (
    analyse_fragment_size_nrl,
    build_parser,
    effective_plot_maximum,
    fragment_size_nrl_quality,
    plot_distributions,
    write_fragment_size_nrl_outputs,
)
from nucleosuite.replot import _read_table, detect_plot_type


def _multinucleosome_counts() -> Counter[int]:
    counts: Counter[int] = Counter()
    for centre, amplitude in zip(
        (150, 330, 510, 690, 870),
        (1000, 800, 600, 400, 250),
    ):
        for length in range(1, 1001):
            value = int(round(amplitude * np.exp(-0.5 * ((length - centre) / 18.0) ** 2)))
            if value:
                counts[length] += value
    return counts


def test_fragment_length_defaults_cover_longest_count_up_to_1000() -> None:
    args = build_parser().parse_args(["--fragments", "sample.bed"])
    assert args.plot_min == 0
    assert args.plot_max == 1000
    assert args.nrl_peak_resolution == 160
    assert args.nrl_min_length == 100
    assert args.nrl_max_length == 1000
    assert effective_plot_maximum({"all": Counter({180: 2, 720: 1})}) == 720
    assert effective_plot_maximum({"all": Counter({180: 2, 1200: 1})}) == 1000


def test_distribution_plot_stops_at_longest_counted_fragment(tmp_path, monkeypatch) -> None:
    import matplotlib.pyplot as plt

    original_close = plt.close
    monkeypatch.setattr(plt, "close", lambda figure=None: None)
    plot_distributions(
        {"all": Counter({100: 2, 720: 1})},
        tmp_path / "lengths.png",
        minimum=0,
        maximum=1000,
        density=True,
    )
    figure = plt.gcf()
    assert figure.axes[0].lines[0].get_xdata()[0] == 0
    assert figure.axes[0].lines[0].get_xdata()[-1] == 720
    original_close(figure)


def test_fragment_size_nrl_uses_resolution_caller_and_regression() -> None:
    result = analyse_fragment_size_nrl({"all": _multinucleosome_counts()})[0]
    assert result.detection_window == 51
    assert result.local_max_window == 21
    assert [peak.distance for peak in result.peaks] == [150, 330, 510, 690, 870]
    assert result.regression.slope == pytest.approx(180.0)
    assert result.regression.r_squared == pytest.approx(1.0)
    assert fragment_size_nrl_quality(result) == "pass"


def test_fragment_size_nrl_outputs_are_replottable(tmp_path: Path) -> None:
    result = analyse_fragment_size_nrl({"all": _multinucleosome_counts()})[0]
    table = tmp_path / "sample_fragment_lengths.tsv"
    table.write_text("fragment_length\tcount\n", encoding="utf-8")
    outputs = write_fragment_size_nrl_outputs([result], table, dpi=80)
    assert len(outputs) == 6
    assert all(path.is_file() for path in outputs)

    stem = tmp_path / "sample_fragment_lengths_peakres160_min100_max1000"
    profile = Path(f"{stem}_fragment_size_nrl_profile.tsv")
    peaks = Path(f"{stem}_fragment_size_nrl_peaks.tsv")
    profile_headers, _ = _read_table(profile)
    peak_headers, _ = _read_table(peaks)
    assert detect_plot_type(profile, profile_headers) == "fragment-size-nrl-profile"
    assert detect_plot_type(peaks, peak_headers) == "fragment-size-nrl-regression"

    with Path(f"{stem}_fragment_size_nrl_summary.tsv").open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["nrl_method"] == "fragment_size_distribution"
    assert float(row["nrl_bp"]) == pytest.approx(180.0)
    assert row["quality_status"] == "pass"
