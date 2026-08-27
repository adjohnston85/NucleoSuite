# `nucleosuite mean-scale`

## What this command does

`mean-scale` expresses BigWig signal values or BED-family interval scores relative to a reference mean:

```math
x_{\mathrm{scaled}} = \frac{x}{\mu_{\mathrm{reference}}} \times S
```

where `S` is `--scale` (default **100**). At the default scale, 100 equals the reference mean, 50 is half the reference mean, and 200 is twice the reference mean.

The command accepts:

- BigWig (`.bw`, `.bigWig`);
- BED;
- BED.gz; and
- bigBed (`.bb`, `.bigBed`).

## Why use it

Use `mean-scale` when signals or interval scores need a common, interpretable scale relative to a measured reference mean, while preserving the original genomic coordinates and file structure.

## BigWig input

For BigWig input, the default reference mean is calculated across finite, non-zero BigWig bases:

```bash
nucleosuite mean-scale coverage.bw
```

Zero-valued bases and missing/non-finite values do not contribute to this calculated mean. Output remains BigWig.

A BED, BED.gz or bigBed can instead define the reference mean through its score column:

```bash
nucleosuite mean-scale PNS.bw \
  --regions nucleosome_regions.bed \
  --score-column 5 \
  --scale 100
```

Or supply a known reference directly:

```bash
nucleosuite mean-scale PNS.bw \
  --reference-mean 16.7644 \
  --scale 100
```

`--normalization-mean` is an alias for `--reference-mean`. `--regions` and `--reference-mean` are mutually exclusive.

A second BigWig can define the reference mean. This is used by `cutn-suite` to divide a centred score track by its matching positive-score mean:

```bash
nucleosuite mean-scale target_tns.bw \
  --reference-bigwig target_posTNS.bw \
  --scale 1
```

`--reference-bigwig`, `--regions`, and `--reference-mean` are mutually exclusive.

## BED, BED.gz and bigBed input

For a BED-family primary input, the default reference is the mean of the finite values in `--score-column` (column 5 by default):

```bash
nucleosuite mean-scale nucleosome_regions.bed
```

If the input scores are `1`, `2`, and `3`, their mean is `2`. With the default `--scale 100`, the output scores are `50`, `100`, and `150`.

All BED columns are preserved; only the selected score column is transformed.

An alternate BED-family reference can be supplied with `--regions`, or `--reference-mean` can be used exactly as for BigWig input.

## Integer scores and clamping

BED and BED.gz output retains floating-point scaled scores by default. Use `--integer-scores` to round them to integers:

```bash
nucleosuite mean-scale peaks.bed \
  --scale 1000 \
  --integer-scores
```

Scaled scores can be clamped explicitly:

```bash
nucleosuite mean-scale peaks.bed \
  --scale 1000 \
  --integer-scores \
  --clamp-min 0 \
  --clamp-max 1000
```

For **bigBed output**, integer conversion and clamping to the standard **0–1000** BED score range are automatic. Explicit `--clamp-min` or `--clamp-max` can narrow that range, but bigBed output can never exceed 0–1000.

## Output format

BED-family input is written in the same format by default:

```text
BED     -> BED
BED.gz  -> BED.gz
bigBed  -> bigBed
```

Override this with:

```text
--output-format bed
--output-format bed.gz
--output-format bigbed
```

For example:

```bash
nucleosuite mean-scale peaks.bed \
  --scale 1000 \
  --output-format bigbed \
  --chrom-sizes genome.chrom.sizes
```

When converting BED or BED.gz to bigBed, chromosome sizes are required. A bigBed primary input can provide its embedded chromosome sizes automatically.

BigWig input is written as BigWig and is not converted to an interval format.

## Automatic output names

If `--output` is omitted, the primary input basename and analysis-defining settings are used automatically. Examples include:

```text
coverage_meanscale_bwnonzero_x100.bw
PNS_meanscale_regions-nucleosome-regions-col5_x100.bw
nucleosome_regions_meanscale_scores-col5_x100.bed
```

Use `--output`/`-o` to override the filename.

A companion `*_mean_scale_summary.tsv` records the input, output format, reference mode, exact reference mean, scale, multiplier, score column where applicable, and interval score conversion/clamping settings.

## Suite use

The cfDNA and MNase suites perform their mean-scaling stage after chromosome combination. Combined coverage and posSNS are mean-scaled to 100. SNS uses the mean score of the raw combined nucleosome calls as its reference, and the combined nucleosome-region and breakpoint-peak BED scores are each mean-scaled to 100. Those mean-scaled peak BEDs are then used for downstream peak-based suite analyses.

[Back to the command reference](../COMMAND_REFERENCE.md)
