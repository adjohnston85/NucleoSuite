# `nucleosuite empirical-peak-fdr`

## What this command does

`empirical-peak-fdr` compares an observed peak callset with one or more peak callsets produced from fragment-randomized controls.

## Why use it

Use this command when observed and randomized peak files already exist and you want a positional-null significance estimate without matching randomized coordinates to observed coordinates.

## Basic usage

```bash
nucleosuite empirical-peak-fdr \
  sample_peaks.bed \
  randomized_peaks.bed \
  --output-prefix sample_randomized_comparison
```

Multiple randomized callsets can be supplied:

```bash
nucleosuite empirical-peak-fdr \
  sample_peaks.bed \
  random_1.bed random_2.bed random_3.bed \
  --fdr 0.05
```

## Method

For every observed peak, the command reports a pooled empirical upper-tail p-value from randomized peak scores and a monotonic empirical FDR/q-value based on observed-versus-randomized score-threshold counts.

The original observed BED fields are preserved. Two fields are appended:

1. `empirical_p_value`
2. `empirical_fdr`

`--score-column` selects the 1-based score column used in every input. `--fdr` optionally writes an additional FDR-filtered BED while retaining the complete annotated output.

See [Empirical randomized-peak FDR](../ALGORITHMS.md#empirical-randomized-peak-fdr) for the equations.

## Outputs

- complete observed peak BED with empirical p-value and FDR appended;
- optional FDR-filtered BED;
- summary TSV describing the observed and randomized inputs and counts.

## Related commands

- [`randomize-fragments`](randomize-fragments.md) — generate coordinate-randomized fragment controls.
- [`nuc-score`](nuc-score.md) — generate SNS/PNS/BNS/TNS score tracks and peaks.

[Back to the command reference](../COMMAND_REFERENCE.md)
