# `nucleosuite pns-peak-fdr`

## What this command does

`pns-peak-fdr` assigns an empirical false-discovery rate to every observed PNS peak by comparing its score with peaks called after the same sample fragments were coordinate-randomized.

The observed and randomized peaks must have been produced with identical PNS mode, fragment limits, contigs, blacklist, score scaling and peak-calling parameters. Randomized peaks are a positional null, so their genomic coordinates are not matched to observed peaks.

## Why use it

Use it to quantify how often peaks at least as strong as each observed peak arise after fragment coordinates are randomized, while retaining a complete annotated peak file for threshold exploration and reproducible filtering.

## Typical run

```bash
nucleosuite pns-peak-fdr \
  sample_nucleosome_regions.bed \
  sample_randomized_control_nucleosome_regions.bed
```

Without `--fdr`, every sample peak is retained. The complete input record is preserved and one final `empirical_fdr` column is appended.

To also write a filtered BED:

```bash
nucleosuite pns-peak-fdr \
  sample_nucleosome_regions.bed \
  sample_randomized_control_nucleosome_regions.bed \
  --fdr 0.05
```

## Multiple randomizations

One or more independently randomized callsets may be supplied:

```bash
nucleosuite pns-peak-fdr sample_peaks.bed \
  random_01_peaks.bed random_02_peaks.bed random_03_peaks.bed
```

At score threshold $s$, with $B$ randomized callsets, the estimated FDR is:

```math
\widehat{\mathrm{FDR}}(s)=\min\left(1,
\frac{1+\sum_{b=1}^{B}R_b(s)}{B\,\max(1,S(s))}\right)
```

where $S(s)$ is the number of observed peaks at or above the threshold and $R_b(s)$ is the corresponding count in randomized callset $b$. The pseudocount prevents a zero estimate when no randomized peak reaches a high score. Monotonic q-values are calculated from these threshold estimates and written as `empirical_fdr`.

## Inputs and outputs

Column 5 is the score by default; select another 1-based column with `--score-column`.

The complete annotated BED is always written. When `--fdr` is supplied, an additional BED contains only peaks whose appended empirical FDR is at or below the requested threshold. A summary TSV records inputs, counts, pseudocount and output paths.

[Back to the command reference](../COMMAND_REFERENCE.md)
