#!/usr/bin/env python3
"""Draw exact DAC/DCC examples as editable SVG and PNG figures.

Run with NucleoSuite installed, or from its source tree with PYTHONPATH=src:
    python plot_dac_dcc_examples.py --output-dir figures

The calculations use NucleoSuite's own sparse, FFT and normalization functions.
No package or documentation files are modified by this script.
"""

from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

from nucleosuite.dac import (
    update_dac_sparse, update_dac_fft, opportunity_vector, build_reported_dac,
)
from nucleosuite.dcc import (
    update_dense_dcc_sparse, update_dense_dcc_fft,
    signed_opportunity_vector, build_reported_dcc,
)


BLUE = "#377bb5"
AMBER = "#e29331"
GREEN = "#2b9b74"
PURPLE = "#9467bd"
INK = "#263544"
GREY = "#aab2ba"


def panel(ax, letter, title):
    ax.set_title(f"{letter}   {title}", loc="left", fontsize=12,
                 fontweight="bold", color=INK, pad=13)


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#aab2ba")
    ax.tick_params(color="#aab2ba", labelsize=9)
    ax.grid(axis="y", color="#e8ebee", lw=.8)
    ax.set_axisbelow(True)


def track(ax, positions, baseline, color, label=None, heights=None):
    positions = np.asarray(positions)
    if heights is None:
        heights = np.full(positions.size, .40)
    ax.axhline(baseline, color="#d9dfe4", lw=.8, zorder=0)
    ax.vlines(positions, baseline, baseline + heights, color=color, lw=2.4)
    ax.scatter(positions, baseline + heights, color=color, s=26, zorder=3)
    if label:
        ax.text(-.018, baseline + .20, label,
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=10, color=color)


def raw_profile(ax, x, y, color, peak_labels=None):
    ax.plot(x, np.zeros_like(y), color=color, lw=1)
    nonzero = y > 1e-12
    ax.vlines(x[nonzero], 0, y[nonzero], color=color, lw=2.5)
    ax.scatter(x[nonzero], y[nonzero], s=36, color=color, zorder=4)
    if peak_labels:
        for at, label in peak_labels.items():
            value = y[np.flatnonzero(x == at)[0]]
            ax.annotate(label, (at, value), xytext=(0, 9), textcoords="offset points",
                        ha="center", va="bottom", fontsize=10, color=color)
    clean(ax)


def canvas(title, subtitle):
    fig = plt.figure(figsize=(12.5, 10.2), facecolor="white")
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.40, 1.45],
                          left=.115, right=.965, top=.86, bottom=.14,
                          hspace=.82, wspace=.34)
    fig.text(.04, .973, title, fontsize=21, weight="bold", color=INK, va="top")
    fig.text(.04, .937, subtitle, fontsize=11, color="#5b6671", va="top")
    return fig, [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, :]),
                 fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])]


def save(fig, destination, stem):
    for extension in ("png", "svg"):
        fig.savefig(destination / f"{stem}.{extension}", dpi=180,
                    facecolor="white", metadata={"Creator": "NucleoSuite figure examples"})
    plt.close(fig)


def dac_example(destination):
    n, dmax, spacing = 925, 800, 185
    positions = np.arange(5) * spacing
    values = np.zeros(n)
    values[positions] = 1
    raw = np.zeros(dmax + 1)
    update_dac_sparse(raw, values, dmax)
    fft = np.zeros_like(raw)
    update_dac_fft(fft, values, dmax)
    np.testing.assert_allclose(raw, fft, atol=1e-12)
    opportunities = opportunity_vector(n, dmax)
    normalized = build_reported_dac(raw, opportunities, True)
    distances = np.arange(1, dmax + 1)
    selected = np.arange(1, 5) * spacing
    np.testing.assert_array_equal(raw[selected], [4, 3, 2, 1])
    np.testing.assert_allclose(normalized[selected], 1 / spacing)

    fig, (a, b, c, d) = canvas(
        "DAC  |  Repeating distances within one signal",
        "Five dyads, 185 bp apart. Each dyad contributes 1 at one genomic base; all other values are 0.",
    )
    panel(a, "A", "Start with a regularly spaced dyad signal")
    a.axhline(0, color=GREY, lw=.8)
    a.vlines(positions, 0, 1, lw=2.4, color=BLUE)
    a.scatter(positions, np.ones(5), color=BLUE, s=33, zorder=3)
    for left, right in zip(positions[:-1], positions[1:]):
        a.annotate("", (left, 1.28), (right, 1.28),
                   arrowprops={"arrowstyle": "<->", "lw": 1, "color": INK})
        a.text((left+right)/2, 1.40, "185 bp", ha="center", fontsize=9, color=INK)
    a.set(xlim=(-25, n), ylim=(-.08, 1.75), yticks=[0, 1],
          xticks=[0, 185, 370, 555, 740, 924], ylabel="Signal", xlabel="Genomic position (bp)")
    clean(a)

    panel(b, "B", "At distance 185 bp, four pairs contribute 1 × 1")
    matching = positions[:-1]
    for x in matching:
        b.axvspan(x-5, x+5, color=GREEN, alpha=.10, lw=0)
    track(b, positions, 2, BLUE, r"$S(x)$")
    track(b, matching, 1, AMBER, r"$S(x+185)$")
    track(b, matching, 0, GREEN, "Product")
    b.text(740, .35, "Sum of products = 4", color=GREEN, fontsize=12,
           ha="center", weight="bold")
    b.set(xlim=(-25, n), ylim=(-.1, 2.7), yticks=[],
          xticks=[0, 185, 370, 555, 740, 924], xlabel="Pair starting position x (bp)")
    b.spines[["top", "right", "left"]].set_visible(False)
    b.spines["bottom"].set_color(GREY)
    b.tick_params(axis="x", labelsize=9, color=GREY)

    panel(c, "C", "Repeat the calculation at every distance")
    raw_profile(c, distances, raw[1:], GREEN,
                {int(x): str(int(raw[x])) for x in selected})
    c.set(xlim=(0, 800), ylim=(-.12, 4.9), yticks=range(5), xticks=[0, *selected],
          xlabel="Distance d (bp)", ylabel="Raw DAC\n(sum of pair products)")

    panel(d, "D", "Divide by the available position pairs")
    raw_profile(d, distances, normalized[1:], PURPLE,
                {int(x): f"{int(raw[x])}/{int(opportunities[x])}" for x in selected})
    d.set(xlim=(0, 800), ylim=(-.0002, .0070), xticks=[0, *selected],
          xlabel="Distance d (bp)", ylabel="Opportunity-normalized DAC")
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, -3))
    d.yaxis.set_major_formatter(formatter)
    fig.text(.115, .035,
             "This 925 bp unmasked region has 925 − d possible position pairs. Zero-signal bases still count.\n"
             "In this constructed example, all four normalized peaks equal 1/185 ≈ 0.00541.",
             fontsize=9.5, color="#5b6671", va="bottom")
    save(fig, destination, "dac_periodicity_example")
    np.savetxt(destination / "DAC-example-data.tsv",
               np.column_stack((distances, raw[1:], opportunities[1:], normalized[1:])),
               delimiter="\t", header="distance_bp\traw_DAC\topportunities\tnormalized_DAC",
               comments="", fmt="%.12g")
    return {"region_length": n, "dyads": positions.tolist(),
            "peak_distances": selected.tolist(), "raw_values": raw[selected].tolist(),
            "opportunities": opportunities[selected].tolist(),
            "normalized_values": normalized[selected].tolist()}


def dcc_example(destination):
    n, dmax, offset = 640, 40, 10
    positions_a = np.asarray([100, 300, 530])
    positions_b = positions_a + offset
    values_a, values_b = np.zeros(n), np.zeros(n)
    values_a[positions_a] = 1
    values_b[positions_b] = 1
    raw = np.zeros(2*dmax + 1)
    update_dense_dcc_sparse(raw, values_a, values_b, dmax)
    fft = np.zeros_like(raw)
    update_dense_dcc_fft(fft, values_a, values_b, dmax)
    np.testing.assert_allclose(raw, fft, atol=1e-12)
    opportunities = signed_opportunity_vector(n, dmax)
    normalized = build_reported_dcc(raw, opportunities, True, False, 3, 3)
    lags = np.arange(-dmax, dmax + 1)
    assert raw[offset+dmax] == 3 and np.count_nonzero(raw) == 1
    assert opportunities[offset+dmax] == 630
    np.testing.assert_allclose(normalized[offset+dmax], 3/630)

    fig, (a, b, c, d) = canvas(
        "DCC  |  Where does signal B occur relative to A?",
        "Three A peaks and three B peaks, each with height 1. B occurs 10 bp downstream of each A peak.",
    )
    panel(a, "A", "Keep the two input signals on separate rows")
    track(a, positions_a, 1, BLUE, "Signal A")
    track(a, positions_b, 0, AMBER, "Signal B")
    for pa, pb in zip(positions_a, positions_b):
        a.text(pa, 1.57, str(pa), ha="center", color=BLUE, fontsize=9)
        a.text(pb, -.20, str(pb), ha="center", color=AMBER, fontsize=9)
        a.annotate("", (pb, .45), (pa, .95),
                   arrowprops={"arrowstyle": "->", "color": INK, "lw": 1})
        a.text(pb+13, .70, "+10 bp", color=INK, fontsize=9, va="center")
    a.set(xlim=(0, n), ylim=(-.36, 1.92), yticks=[],
          xticks=[0, 100, 200, 300, 400, 500, 600], xlabel="Genomic position (bp)")
    a.spines[["top", "right", "left"]].set_visible(False)
    a.spines["bottom"].set_color(GREY)
    a.tick_params(axis="x", labelsize=9, color=GREY)

    panel(b, "B", "At lag +10 bp, compare A(x) with B(x + 10)")
    for x in positions_a:
        b.axvspan(x-5, x+5, color=GREEN, alpha=.10, lw=0)
    track(b, positions_a, 2, BLUE, r"$A(x)$")
    track(b, positions_b-offset, 1, AMBER, r"$B(x+10)$")
    track(b, positions_a, 0, GREEN, "Product")
    b.text(315, .5, "1 × 1", color=GREEN, ha="left", fontsize=9)
    b.text(595, 1.15, "3 matches\n\nRaw DCC = 3", color=GREEN,
           fontsize=10, ha="center", va="center", weight="bold")
    b.set(xlim=(0, n), ylim=(-.1, 2.7), yticks=[],
          xticks=[0, 100, 200, 300, 400, 500, 600], xlabel="A position x (bp)")
    b.spines[["top", "right", "left"]].set_visible(False)
    b.spines["bottom"].set_color(GREY)
    b.tick_params(axis="x", labelsize=9, color=GREY)

    panel(c, "C", "Signed lags retain the direction")
    raw_profile(c, lags, raw, GREEN, {10: "3 matches"})
    c.axvline(0, color=GREY, ls=":", lw=1)
    c.text(-22, 3.58, "B upstream", ha="center", fontsize=9, color="#5b6671")
    c.text(22, 3.58, "B downstream", ha="center", fontsize=9, color="#5b6671")
    c.set(xlim=(-40, 40), ylim=(-.1, 4.0), yticks=range(4),
          xticks=[-40, -20, 0, 10, 20, 40],
          xlabel="Lag ℓ = position B − position A (bp)", ylabel="Raw signed DCC")

    panel(d, "D", "Normalize by available A/B position pairs")
    raw_profile(d, lags, normalized, PURPLE, {10: "3/630 = 0.00476"})
    d.axvline(0, color=GREY, ls=":", lw=1)
    d.set(xlim=(-40, 40), ylim=(-.0002, .0066), xticks=[-40, -20, 0, 10, 20, 40],
          xlabel="Signed lag ℓ (bp)", ylabel="Opportunity-normalized DCC")
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, -3))
    d.yaxis.set_major_formatter(formatter)
    fig.text(.115, .035,
             "Signed lags shown (−40 to +40 bp): use --signed-lags. At +10 bp, a 640 bp region has 630 possible A/B pairs.\n"
             "A peak at +10 bp means B repeatedly occurs 10 bp downstream of A in the active coordinate orientation.",
             fontsize=9.5, color="#5b6671", va="bottom")
    save(fig, destination, "dcc_shift_example")
    np.savetxt(destination / "DCC-example-data.tsv",
               np.column_stack((lags, raw, opportunities, normalized)),
               delimiter="\t", header="lag_bp\traw_DCC\topportunities\tnormalized_DCC",
               comments="", fmt="%.12g")
    return {"region_length": n, "positions_A": positions_a.tolist(),
            "positions_B": positions_b.tolist(), "peak_lag": offset,
            "raw_value": 3, "opportunities": 630, "normalized_value": 3/630}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "savefig.facecolor": "white"})
    results = {"DAC": dac_example(args.output_dir), "DCC": dcc_example(args.output_dir),
               "verification": "Sparse and FFT results agree for both examples."}
    (args.output_dir / "example-checks.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
