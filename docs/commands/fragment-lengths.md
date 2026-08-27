# `nucleosuite fragment-lengths`

## What this command does

`fragment-lengths` counts how many accepted fragments occur at each integer fragment length. It can calculate one whole-sample distribution or separate distributions for labelled genomic regions such as chromatin states. By default it also calls mono-, di-, tri- and higher multinucleosome fragment-size peaks and regresses their summit positions to estimate fragment-size NRL.

## Why use it

Use it to:

- compare fragment-size distributions between samples;
- compare fragment lengths between chromatin states or other region classes;
- prepare count tables for [`fragment-heatmap`](fragment-heatmap.md); or
- check which fragment lengths dominate the exact fragment population used by another analysis; or
- estimate NRL from a multinucleosome fragment-size ladder, including low-coverage data for which a dyad autocorrelation may not be practical.

## Basic usage

```bash
nucleosuite fragment-lengths \
  --bam sample.bam \
  --min-length 1 \
  --max-length 1000 \
  --output sample_lengths.tsv \
  --plot sample_lengths.png
```

## Count fragment lengths by chromatin state

The bundled GM12878 state BED can be used directly:

```bash
nucleosuite fragment-lengths \
  --bam sample.bam \
  --bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --bed-label-column 4 \
  --output sample_lengths_by_state.tsv
```

A fragment is assigned to a region using its midpoint. The exact midpoint definition is given in [Fragment-length counts](../ALGORITHMS.md#fragment-length-counts-and-heatmap-transformations).

## What the counts mean

The primary TSV contains raw counts for each fragment length. `fragment-heatmap` derives percentages or z-scores from this table for plotting.

## Fragment-size NRL

For each fragment-length profile, the command analyses lengths from 100 bp through the longest counted fragment, capped at 1000 bp by default. The broad peak detector uses the same resolution-driven method as [`nrl`](nrl.md): the default 160 bp resolution gives 61 bp detection smoothing and 21 bp local-maximum refinement. The called peak summits are ordered as the mono-, di-, tri- and higher multinucleosome series. Fragment length is regressed against peak number, and the slope is reported as the fragment-size NRL.

The regression summary includes the number of retained peaks, NRL, intercept, R², slope standard error, mean adjacent spacing and a quality status. Fewer than three peaks are marked `insufficient_peaks`; a fit with R² below 0.9 is marked `low_r_squared`. The numerical result is still written so the caller is transparent, but these statuses should not be treated as a confident NRL estimate.

Change the analysis range or caller resolution with:

```bash
nucleosuite fragment-lengths \
  --bam sample.bam \
  --nrl-min-length 120 \
  --nrl-max-length 900 \
  --nrl-peak-resolution 160 \
  --output sample_lengths.tsv
```

Use `--no-fragment-size-nrl` when only the count table is wanted. This fragment-size method is kept distinct from dyad-distance DAC/DCC NRL because NRL values should be compared using a consistent method. The method follows the fragment-size strategy described by [Bikova, Clarkson and Teif (2026)](https://academic.oup.com/nar/article/54/5/gkag074/8506906).

## Inputs

Use indexed paired-end BAM input or materialized fragment BED/BED.gz/bigBed input. Fragment filters and duplicate-coordinate settings determine which fragments enter the counts.

## Outputs

Outputs include the raw fragment-length count table and, when `--plot` is supplied, a fragment-length distribution figure. The plot starts at 0 bp and stops at the longest counted fragment or 1000 bp, whichever is shorter. `--plot-max` can impose a lower upper limit. Region-aware runs include separate rows/profiles for the requested region labels.

Unless `--no-fragment-size-nrl` is used, each profile also writes:

- `_fragment_size_nrl_profile.tsv` and `.png`, containing the unsmoothed density, 21 bp local refinement curve, 61 bp detection curve and called peaks;
- `_fragment_size_nrl_peaks.tsv`, containing peak number, observed and fitted fragment length, residual and peak signals;
- `_fragment_size_nrl_regression.tsv` and `.png`, containing the fit summary and the square open-circle/dotted-line regression figure; and
- `_fragment_size_nrl_summary.tsv`, collecting all label-level fits.

For labelled region runs, the label is inserted before `_fragment_size_nrl` in the per-profile filenames.

## Chromosome-wise processing

With multiple contigs, NucleoSuite sums per-contig raw counts and regenerates percentages and plots from the combined table.

## Plot customization

Distribution figures use the shared plotting options in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
