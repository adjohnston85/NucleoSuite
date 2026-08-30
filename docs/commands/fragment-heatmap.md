# `nucleosuite fragment-heatmap`

## What this command does

`fragment-heatmap` turns one or more fragment-length count tables into a matrix and heatmap. It is designed for comparing either **different samples** or **different region classes from the same sample** across fragment length.

## Why use it

Use it when a line plot becomes difficult to compare across many samples/states, or when you want to highlight fragment lengths that are relatively enriched or depleted in one profile.

## Start from raw counts

The input should be fragment-length counts, typically produced by [`fragment-lengths`](fragment-lengths.md). The command then applies the selected normalization for visualization.

## Choose the normalization for the comparison

### `profile-percent`

Each row sums to 100%. Use this when you want to compare the **shape of each sample/profile's own fragment-length distribution**.

### `fragment-percent`

Each fragment-length column sums to 100% across profiles, showing **which profiles contribute most strongly to each fragment length**.

### `fragment-zscore` — default

The command first converts each row to profile percentages, then calculates a z-score across profiles separately for each fragment length. This highlights lengths that are higher or lower than expected for that length across the compared profiles.

With the default blue-white-orange palette:

- blue = below the across-profile mean for that fragment length;
- white = near the mean;
- orange = above the mean.

### Min-max scaling

Min-max modes rescale a row or column to 0–1. Use these when relative shape matters more than absolute differences.

The exact formulas for all transformations are in [Fragment-length counts and heatmap transformations](../ALGORITHMS.md#fragment-length-counts-and-heatmap-transformations).

### Worked example

Suppose only two fragment lengths are selected. Profile A has counts `[40,60]` at lengths 145 and 167 bp; profile B has `[80,20]`.

| Representation | A: 145 bp | A: 167 bp | B: 145 bp | B: 167 bp |
|---|---:|---:|---:|---:|
| Raw counts | 40 | 60 | 80 | 20 |
| `profile-percent` | 40% | 60% | 80% | 20% |
| `fragment-percent` | 33.33% | 75% | 66.67% | 25% |
| `fragment-zscore` | −1 | +1 | +1 | −1 |

Profile percentages compare the distribution within each profile. Fragment percentages compare profiles at one length: A supplies `40 / (40 + 80) = 33.33%` of the 145 bp fragments. Z-scores compare profile percentages with their across-profile mean and population standard deviation; they are neither percentages nor significance tests. All these calculations use the selected fragment-length range. An empty row/column, or a z-score column with no variation, is reported as zero.

## Optional downsampling

Downsampling is applied to raw integer counts before normalization. Use it when profiles have very different total fragment counts and you want to compare equal-sized random samples.

`--downsample-to min` uses the smallest positive profile total as the target. Profiles already at or below the target are retained at their original counts.

## Basic usage

```bash
nucleosuite fragment-heatmap \
  --input sample1_fragment_length_counts.tsv \
  --input sample2_fragment_length_counts.tsv \
  --normalization fragment-zscore \
  --out-prefix fragment_length_comparison
```

## Outputs

The command writes:

- `_heatmap.png`;
- `_normalised_matrix.tsv`, containing the displayed values in clustered row order;
- `_clustered_profiles.tsv` and `_clustered_fragment_stats.tsv`;
- `_heatmap_plot_metadata.tsv` and `_heatmap_linkage.tsv`, which let `nucleosuite plot` faithfully restore the complete figure.

Add `--write-detail-tables` when the additional Excel workbook is required.

Recreate the original figure with:

```bash
nucleosuite plot fragment_length_comparison_normalised_matrix.tsv
```

The normalized matrix remains a default output because it is the compact numerical source required to reproduce the heatmap. The Excel workbook is a larger supporting output and is written only with `--write-detail-tables`.

## Plot customization

Heatmaps retain their command-specific palette/normalization options and also accept the shared figure-size, format, title, label, tick, and transparency options in [Plot customization](../PLOTTING.md).

## Automatic output naming

If `--out-prefix` is omitted, NucleoSuite derives the prefix from the first input-table basename and appends `_fragment_heatmap`. An explicit prefix overrides the automatic name.

[Back to the command reference](../COMMAND_REFERENCE.md)
