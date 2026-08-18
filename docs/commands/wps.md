# `nucleosuite wps`

## What this command does

`wps` calculates window protection score (WPS) tracks from paired-end fragments and can call nucleosome and breakpoint features from a selected WPS-derived signal.

## Why use it

Use WPS to compare complete-window protection with fragment termination at each genomic position. WPS signals and calls can be used for cfDNA footprinting, regional profiles, spacing, and callset comparisons.

## How it works

For each genomic window centre, one fragment contributes:

- **+1** if it spans the complete protection window;
- **−1** if a fragment endpoint lies inside the window; or
- **0** if neither condition applies.

Raw WPS is the sum of those fragment contributions. With the default 120 bp protection window, a 120 bp fragment has one central +1 position; longer fragments have a wider +1 interior.

NucleoSuite can smooth raw WPS and subtract a 1,000 bp running-median baseline. The default peak caller uses the smoothed, median-centred track `sm_mWPS`.

See [Windowed protection score](../ALGORITHMS.md#windowed-protection-score) for the kernel examples, exact preprocessing, and peak-calling steps.

## Typical use

```bash
nucleosuite wps \
  --bam sample.bam \
  --frag-lower 120 \
  --frag-upper 180 \
  --protection 120 \
  --score-format bigwig \
  --out-prefix sample_wps
```

## Defaults

NucleoSuite's WPS implementation was written to reproduce the L-WPS algorithm used by [Snyder et al.](https://doi.org/10.1016/j.cell.2015.11.050). The defaults below match the parameters reported for that method.

| Setting | Default |
|---|---:|
| Fragment range | 120–180 bp |
| Protection window | 120 bp |
| Savitzky–Golay window/order | 21 bp / 2 |
| Running-median baseline | 1,000 bp |
| Peak input | `sm_mWPS` |
| Positive-region range | 50–450 bp |
| Retained long-region subrun size | 50–150 bp |
| Positive-position merge gap | 5 bp |
| Peak maximum cutoff | strictly greater than 5 |

## Peak calling

The WPS caller:

1. finds positive regions of the selected adjusted WPS signal;
2. keeps the stronger, above-median parts of each region;
3. selects one or more appropriately sized peak-like blocks;
4. requires the retained block to exceed the WPS cutoff; and
5. reports the retained block midpoint as the call centre.

See [WPS peak calling](../ALGORITHMS.md#wps-peak-calling) for the exact segmentation and selection rules.

## What it writes

WPS can write raw, smoothed, median-adjusted, and smoothed median-adjusted signals. It can also write coverage and dyad tracks, nucleosome and breakpoint calls, and fragment summaries.

Use `--score-tracks` to choose signal outputs, `--peak-track` to choose the calling signal, and `--peak-caller none` to omit peak calling.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. (2016). Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. *Cell* 164, 57–68. https://doi.org/10.1016/j.cell.2015.11.050
- Savitzky A, Golay MJE. (1964). Smoothing and Differentiation of Data by Simplified Least Squares Procedures. *Analytical Chemistry* 36, 1627–1639. https://doi.org/10.1021/ac60214a047
