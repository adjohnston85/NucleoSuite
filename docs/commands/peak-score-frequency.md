# `nucleosuite peak-score-frequency`

## What this command does

`peak-score-frequency` compares the score distributions of one or more peak callsets using shared histogram bins.

## Why use it

Use `peak-score-frequency` to inspect or compare the distribution of peak strengths across one or more callsets on shared bins, with optional display scaling when the stored score units are inconvenient for plotting.

## Typical use

A single input needs no explicit output prefix:

```bash
nucleosuite peak-score-frequency \
  --peaks sample_nucleosome_regions.bed
```

The primary input basename is used for the automatic output name.

Multiple labelled inputs can be overlaid:

```bash
nucleosuite peak-score-frequency \
  --peaks observed=sample_nucleosome_regions.bed \
  --peaks control=control_nucleosome_regions.bed \
  --score-column 5 \
  --integer-bins \
  --normalization count
```

## Score scaling

`--score-scale` multiplies scores **for histogram binning and plotting**. Its standalone default is:

```text
--score-scale 1
```

So BED, BED.gz and bigBed scores are all plotted as stored unless scaling is explicitly requested.

For decimal-valued scores that should be displayed after multiplication by 1000:

```bash
nucleosuite peak-score-frequency \
  --peaks PNS=sample_nucleosome_regions.bed \
  --score-scale 1000
```

The x-axis then reads:

```text
Peak score ×1000
```

The compact frequency table uses the scaled score axis. Per-input summary tables and optional individual-score detail tables remain in the original input-score units.

The cfDNA and MNase suites first mean-scale their combined SNS-derived nucleosome-region and breakpoint-peak BED scores to 100, then pass those normalized BEDs to `peak-score-frequency` with `--score-scale 1`. No additional suite-level histogram scaling is applied.

## Histogram representation

Integer-bin mode rounds display scores to the nearest integer and writes every integer from the selected minimum to maximum, including zero-count values between observed scores. It is the default when neither `--bins` nor `--bin-width` is supplied.

Use `--bins N` for a fixed number of shared continuous bins or `--bin-width N` for a chosen continuous bin width. Shared boundaries are used across all labelled inputs.

## Normalization

`--normalization` can report:

- `count` — number of peaks in each bin;
- `fraction` — fraction of peaks;
- `percent` — percentage of peaks; or
- `density` — frequency adjusted for continuous bin width.

## Outputs

Default outputs include the shared score-frequency table, per-input score summary, and figure. `--write-detail-tables` additionally writes the compressed table of individual finite peak scores.

Automatic output names begin with the first peak-input basename and include central histogram parameters. `--output-prefix`/`-o` overrides the automatic base.

## Blacklist handling

`--blacklist-bed` excludes complete overlapping peak intervals before the score distribution is calculated.

## Plot customization

Figures use the shared options in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
