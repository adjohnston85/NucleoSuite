# `nucleosuite wps`

## What this command does

`wps` calculates window protection score (WPS) tracks from paired-end fragments and can call nucleosome and breakpoint features from a selected WPS-derived signal.

## Why use it

Use WPS to compare complete-window protection with fragment termination at each genomic position. WPS signals and calls can be used for cfDNA footprinting, regional profiles, spacing, and callset comparisons.

## How it works

Before WPS is calculated, the default `--mode auto` estimates the dominant accepted fragment length from an unsmoothed fragment-length histogram. The resolved mode becomes the centred protection-window width. This adapts the protected window to the library instead of assuming that every assay has the same modal protected length. Use an integer such as `--mode 120` to reproduce a fixed 120 bp L-WPS protection window.

For each genomic window centre, one fragment contributes:

- **+1** if it spans the complete protection window;
- **−1** if a fragment endpoint lies inside the window; or
- **0** if neither condition applies.

Raw WPS is the sum of those fragment contributions. With a 120 bp protection window, a 120 bp fragment has one central +1 position; longer fragments have a wider +1 interior.

NucleoSuite can smooth raw WPS and subtract a 1,000 bp running-median baseline. The default peak caller uses the smoothed, median-centred track `sm_mWPS`.

See [Windowed protection score](../ALGORITHMS.md#windowed-protection-score) for the kernel examples, exact preprocessing, and peak-calling steps.

## Basic usage

```bash
nucleosuite wps \
  --bam sample.bam \
  --frag-lower 120 \
  --frag-upper 180 \
  --score-format bigwig \
  --out-prefix sample_wps
```

## Defaults

NucleoSuite's WPS implementation was written to reproduce the L-WPS algorithm used by [Snyder et al.](https://doi.org/10.1016/j.cell.2015.11.050). Signal preprocessing and peak-calling defaults follow that method, while the protection width now defaults to an automatically estimated fragment mode. Supply `--mode 120` for the original fixed 120 bp protection width.

| Setting | Default |
|---|---:|
| Fragment range | 120–180 bp |
| Protection window | Automatically estimated fragment mode |
| Savitzky–Golay window/order | 21 bp / 2 |
| Running-median baseline | 1,000 bp |
| Peak input | `sm_mWPS` |
| Positive-region range | 50–450 bp |
| Retained long-region subrun size | 50–150 bp |
| Positive-position merge gap | 5 bp |
| Peak maximum cutoff | strictly greater than 5 |

Automatic estimation uses raw integer histogram counts by default. Optional `--mode-histogram-smoothing binomial` applies the normalized `1,4,6,4,1` kernel. This option affects only fragment-mode estimation. It does not change the 21 bp Savitzky–Golay smoothing applied to the genomic WPS signal; use `--sg-window 0` to disable that separate signal-processing step.

The resolved mode is printed and written to `*_fragment_mode_estimation.tsv` with its bootstrap interval, sampling counts, search range, convergence result, seed, smoothing setting, and histogram. Output filenames contain the resolved numeric protection length.

## Peak calling

The WPS caller:

1. finds positive regions of the selected adjusted WPS signal;
2. keeps the stronger, above-median parts of each region;
3. selects one or more appropriately sized peak-like blocks;
4. requires the retained block to exceed the WPS cutoff; and
5. reports the retained block midpoint as the call centre.

See [WPS peak calling](../ALGORITHMS.md#wps-peak-calling) for the exact segmentation and selection rules.

## Outputs

WPS can write raw, smoothed, median-adjusted, and smoothed median-adjusted signals. It can also write coverage and dyad tracks, nucleosome and breakpoint calls, and fragment summaries.

Use `--score-tracks` to choose signal outputs, `--peak-track` to choose the calling signal, and `--peak-caller none` to omit peak calling.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. (2016). Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. *Cell* 164, 57–68. https://doi.org/10.1016/j.cell.2015.11.050
- Savitzky A, Golay MJE. (1964). Smoothing and Differentiation of Data by Simplified Least Squares Procedures. *Analytical Chemistry* 36, 1627–1639. https://doi.org/10.1021/ac60214a047
