# `nucleosuite peak-score-frequency`

## What this command does

`peak-score-frequency` compares the score distributions of one or more peak callsets using shared histogram bins.

## Why use it

Use it to compare proportions of weak and strong peaks between:

- observed versus randomized controls;
- PNS versus WPS score distributions;
- two samples; or
- two parameter settings of the same caller.

[`compare-positions`](compare-positions.md) measures whether the callsets contain the same genomic positions.

## Typical use

```bash
nucleosuite peak-score-frequency \
  --peaks observed=sample_nucleosome_regions.bed \
  --peaks randomized=sample_randomized_nucleosome_regions.bed \
  --score-column 5 \
  --integer-bins \
  --normalization count \
  --output-prefix sample_nucleosome_scores
```

## Choose the histogram representation

Integer-bin mode rounds scores to the nearest integer and writes every integer from the selected minimum to maximum, including zero-count values between observed scores.

Use `--bins N` for a fixed number of shared continuous bins or `--bin-width N` for a chosen continuous bin width.

Shared boundaries are used across all labelled inputs so their distributions are directly comparable.

## Choose what the y-axis means

`--normalization` can report:

- `count` — number of peaks in each bin;
- `fraction` — fraction of peaks;
- `percent` — percentage of peaks; or
- `density` — frequency adjusted for continuous bin width.

The underlying frequency table retains the corresponding count and normalized columns.

## What it writes

Outputs include:

- the finite individual scores used;
- the shared score-frequency table;
- per-input score summaries; and
- a comparison figure.

Plot x/y limits do not remove values from the TSV or summary.

## Blacklist handling

`--blacklist-bed` removes complete overlapping peak intervals before the score distribution is calculated.

## Plot customization

Figures use the shared options in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
