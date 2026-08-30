# `nucleosuite positive-runs`

## What this command does

`positive-runs` finds continuous stretches of a BigWig signal that remain above 0 and reports their lengths. By default, runs are defined by values strictly greater than 0, but an alternative threshold can be supplied with `--threshold`.

## Why use it

Use it to measure the lengths of positive-signal regions and compare those run-length distributions between tracks. Raise `--threshold` when a cutoff above zero is more useful for the analysis.

## How it works

A run begins when the signal becomes greater than `--threshold` and continues until the signal falls to or below that threshold or another boundary ends the run.

Missing/non-finite positions, genomic gaps, selected-region boundaries, contig boundaries, and blacklist masks also terminate a run.

The exact interval and signal-area definitions are in [Positive-signal runs](../ALGORITHMS.md#positive-signal-runs).

## Basic usage

```bash
nucleosuite positive-runs \
  --bigwig sample_pns.bw \
  --threshold 0 \
  --output-prefix sample_pns_positive_runs
```

Use the same threshold when comparing samples if the run-length distributions are meant to be directly comparable.

## Outputs

For each retained run, NucleoSuite reports:

- interval coordinates;
- run length;
- maximum signal;
- mean signal; and
- summed signal area.

Count/distribution tables and plots summarize those runs.

For a three-base run with values `[2,4,6]`, the length is 3 bp, maximum is 6, mean is 4, and summed signal area is 12 score·bp. The sum grows with run length; the mean describes average signal per base.

## Relation to peak calling

A positive run uses only the threshold and interval boundaries. PNS and WPS peak callers apply their signal-specific region, median, and cutoff rules.

## Plot customization

Run-length figures use the shared plotting options in [Plot customization](../PLOTTING.md). The default `--plot-x-max 550` is a ceiling: if the longest observed run is shorter than 550 bp, the x-axis stops at the observed maximum instead of extending to 550 bp.

[Back to the command reference](../COMMAND_REFERENCE.md)
