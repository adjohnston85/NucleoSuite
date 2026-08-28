# `nucleosuite wps`

## What this command does

`wps` calculates window protection score (WPS) tracks from paired-end fragments and can call nucleosome and breakpoint features from a selected WPS-derived signal.

## Why use it

Use WPS to compare complete-window protection with fragment termination at each genomic position. WPS signals and calls can be used for cfDNA footprinting, regional profiles, spacing, and callset comparisons.

## How it works

WPS uses a fixed centred protection window. The default protection width is **120 bp**, and the default accepted fragment range is **120–180 bp**. These values are independent of fragment-mode estimation.

For each genomic window centre, one fragment contributes:

- **+1** if it spans the complete protection window;
- **−1** if a fragment endpoint lies inside the window; or
- **0** if neither condition applies.

Raw WPS is the sum of those fragment contributions. NucleoSuite can smooth raw WPS and subtract a 1,000 bp running-median baseline. The default peak caller uses the smoothed, median-centred `sm_mWPS` track.

See [Windowed protection score](../ALGORITHMS.md#windowed-protection-score) for the exact kernel geometry and peak-calling steps.

## Basic usage

```bash
nucleosuite wps \
  --bam sample.bam \
  --out-prefix sample_wps
```

The defaults used by this command are equivalent to:

```bash
nucleosuite wps \
  --bam sample.bam \
  --frag-lower 120 \
  --frag-upper 180 \
  --protection 120 \
  --out-prefix sample_wps
```

## Core options

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

Use `--protection` to change the protection-window width. `--frag-lower` and `--frag-upper` change the accepted fragment range. `--sg-window 0` disables Savitzky–Golay smoothing.

## Outputs

WPS can write raw, smoothed, median-adjusted, and smoothed median-adjusted signals. It can also write coverage and dyad tracks, nucleosome and breakpoint calls, and fragment summaries.

Use `--score-tracks` to choose signal outputs, `--peak-track` to choose the calling signal, and `--peak-caller none` to omit peak calling.

## Related commands

- [`pns`](pns.md) — generate the length-adaptive PNS signal and peak calls.
- [`call-peaks`](call-peaks.md) — call peaks from an existing score track.
- [`tracks`](tracks.md) — generate several fragment-derived tracks in one pass.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. (2016). Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin. *Cell* 164, 57–68. https://doi.org/10.1016/j.cell.2015.11.050
- Savitzky A, Golay MJE. (1964). Smoothing and Differentiation of Data by Simplified Least Squares Procedures. *Analytical Chemistry* 36, 1627–1639. https://doi.org/10.1021/ac60214a047
