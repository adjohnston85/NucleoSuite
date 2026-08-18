"""Regression tests for DAC and DCC plot-only distance filtering."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.figure
import numpy as np

from nucleosuite.dac import plot_dac_tsv, write_dac_tsv
from nucleosuite.dcc import plot_dcc_tsv, write_dcc_tsv


def _capture_first_line(monkeypatch):
    captured: dict[str, object] = {}

    def capture_savefig(figure, *_args, **_kwargs):
        figure.canvas.draw()
        axis = figure.axes[0]
        captured["x"] = np.asarray(axis.lines[0].get_xdata(), dtype=float)
        captured["y"] = np.asarray(axis.lines[0].get_ydata(), dtype=float)
        captured["ticks"] = np.asarray(axis.get_xticks(), dtype=float)

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", capture_savefig)
    return captured


def test_dac_tsv_retains_distance_one_but_png_omits_it(tmp_path: Path, monkeypatch):
    tsv = tmp_path / "dac.tsv"
    write_dac_tsv(
        str(tsv),
        raw_dac=np.asarray([0.0, 1000.0, 4.0, 3.0]),
        opportunities=np.asarray([0.0, 10.0, 10.0, 10.0]),
        normalize_dac=False,
        total_signal=20.0,
        cpm_scale=1_000_000.0,
    )
    rows = tsv.read_text().splitlines()
    assert any(row.startswith("1\t") for row in rows[1:])

    captured = _capture_first_line(monkeypatch)
    plot_dac_tsv(str(tsv), str(tmp_path / "dac.png"))
    assert captured["x"].tolist() == [2.0, 3.0]


def test_absolute_dcc_tsv_retains_distance_one_but_png_omits_it(tmp_path: Path, monkeypatch):
    tsv = tmp_path / "dcc_absolute.tsv"
    write_dcc_tsv(
        str(tsv),
        raw_dcc=np.asarray([2.0, 1000.0, 4.0]),
        opportunities=np.asarray([10.0, 10.0, 10.0]),
        dmax=2,
        signed_lags=False,
        normalize_dcc=False,
        normalize_by_signal_totals=False,
        total_signal_a=10.0,
        total_signal_b=10.0,
        cpm_scale=1_000_000.0,
    )
    rows = tsv.read_text().splitlines()
    assert any(row.startswith("1\t") for row in rows[1:])

    captured = _capture_first_line(monkeypatch)
    plot_dcc_tsv(str(tsv), str(tmp_path / "dcc_absolute.png"))
    assert captured["x"].tolist() == [0.0, 2.0]


def test_signed_dcc_png_omits_both_one_base_pair_lags(tmp_path: Path, monkeypatch):
    tsv = tmp_path / "dcc_signed.tsv"
    write_dcc_tsv(
        str(tsv),
        raw_dcc=np.asarray([3.0, 900.0, 2.0, 1000.0, 4.0]),
        opportunities=np.asarray([10.0, 10.0, 10.0, 10.0, 10.0]),
        dmax=2,
        signed_lags=True,
        normalize_dcc=False,
        normalize_by_signal_totals=False,
        total_signal_a=10.0,
        total_signal_b=10.0,
        cpm_scale=1_000_000.0,
    )
    rows = tsv.read_text().splitlines()
    assert any(row.startswith("-1\t") for row in rows[1:])
    assert any(row.startswith("1\t") for row in rows[1:])

    captured = _capture_first_line(monkeypatch)
    plot_dcc_tsv(str(tsv), str(tmp_path / "dcc_signed.png"))
    assert captured["x"].tolist() == [-2.0, 0.0, 2.0]
