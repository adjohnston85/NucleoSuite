# `nucleosuite mean-scale`

## What this command does

`mean-scale` writes a BigWig whose values are expressed relative to a reference mean:

```math
x_{\mathrm{scaled}} = \frac{x}{\mu_{\mathrm{reference}}} \times S
```

where `S` is `--scale` (default **100**). With the default scale, a value of 100 is equal to the reference mean, 50 is half the reference mean, and 200 is twice the reference mean.

This is **mean scaling** (mean-ratio scaling), rather than centring or z-score normalization.

## Why use it

Use mean scaling when tracks from one dataset should be expressed relative to a biologically or technically meaningful reference mean while preserving their proportional signal structure. Typical uses include expressing coverage relative to mean non-zero coverage, expressing PNS signal relative to the mean score of nucleosome-protection peaks, or applying one previously calculated reference mean consistently across several related BigWigs.

## Reference-mean modes

The command chooses the reference mean in one of three ways.

### 1. BigWig-derived mean (default)

If neither `--regions` nor `--reference-mean` is supplied, the reference mean is calculated across all **finite, non-zero BigWig base values**. Zero-valued bases and missing/non-finite values do not contribute to the mean.

```bash
nucleosuite mean-scale coverage.bw
```

This is suitable for a coverage track when the desired reference is the mean non-zero coverage.

### 2. Region-score mean

Supply a BED, BED.gz or bigBed file with `--regions`. The reference mean is the unweighted arithmetic mean of the finite region scores. Scores are read from **column 5 by default** and can be changed with `--score-column`.

```bash
nucleosuite mean-scale PNS.bw \
  --regions nucleosome_protection_peaks.bb
```

For example, if the mean nucleosome-protection peak score is 16.7644, every PNS BigWig value is divided by 16.7644 and multiplied by 100.

The region file is used only to determine the reference mean; this command writes a scaled BigWig and does not modify the supplied BED/bigBed.

### 3. Supplied reference mean

Use `--reference-mean` when the desired mean is already known:

```bash
nucleosuite mean-scale PNS.bw \
  --reference-mean 16.7644
```

`--normalization-mean` is accepted as an alias for `--reference-mean`.

`--regions` and `--reference-mean` are mutually exclusive.

## Scale

The default is:

```text
--scale 100
```

A different multiplier can be supplied explicitly:

```bash
nucleosuite mean-scale PNS.bw \
  --reference-mean 16.7644 \
  --scale 1000
```

The reference mean must be finite and non-zero, and `--scale` must be finite and greater than zero.

## Outputs

The default output filename records the reference mode and scale. Examples include:

```text
coverage_meanscale_bwnonzero_x100.bw
PNS_meanscale_regions-nucleosome-protection-peaks-col5_x100.bw
PNS_meanscale_mean-16p7644_x100.bw
```

Use `--output FILE.bw` to choose the BigWig filename explicitly.

A companion `*_mean_scale_summary.tsv` records:

- input and output BigWigs;
- reference-mean mode;
- exact reference mean;
- scale and resulting multiplier;
- region file and score column when applicable; and
- the number of finite values that contributed to a calculated mean.

## Example matching a PNS resource workflow

Coverage scaled relative to mean non-zero coverage:

```bash
nucleosuite mean-scale CH01_coverage.bw \
  --scale 100
```

PNS signal scaled relative to the mean score of nucleosome-protection peaks:

```bash
nucleosuite mean-scale CH01_PNS.bw \
  --regions CH01_nucleosome_protection_peaks.bb \
  --score-column 5 \
  --scale 100
```

A second signal track can be placed on exactly the same reference scale by supplying the previously calculated mean directly:

```bash
nucleosuite mean-scale CH01_breakpoint_signal.bw \
  --reference-mean 16.7644 \
  --scale 100
```

[Back to the command reference](../COMMAND_REFERENCE.md)
