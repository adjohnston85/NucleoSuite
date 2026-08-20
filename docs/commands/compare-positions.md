# `nucleosuite compare-positions`

## What this command does

`compare-positions` compares one main nucleosome callset with one or more comparison callsets using one-to-one positional matching. For each comparison, matched pairs are ranked by the **main BED score** and divided into percentile groups; quartiles are used by default.

## Why use it

Use this command when several nucleosome callsets should be compared against one common reference while retaining the main callset score for stratification and statistical analysis.

## Typical use

```bash
nucleosuite compare-positions \
  --main-bed PNS=PNS_nucleosomes.bed \
  --compare-bed iNPS=iNPS_nucleosomes.bed \
  --compare-bed DANPOS=DANPOS_nucleosomes.bed \
  --compare-bed WPS=WPS_nucleosomes.bed \
  --stats \
  --output-prefix Gaffney32_position_compare
```

Both the main BED and comparison BEDs accept `LABEL=path.bed`. The label is used consistently in plot labels, legends and output tables. If no label is supplied, the BED filename is used.


## Matching

Each main-versus-comparison pair is searched once. The callset with fewer positions is used as the query and the larger callset as the target. Matching is one-to-one, so a position can occur in at most one accepted pair.

Regardless of search direction, every accepted pair is represented as a main call and a comparison call. Signed distance is:

```text
comparison summit - main summit
```

Use `--max-distance` to reject pairs beyond a specified absolute summit distance.

By default, the position is the integer midpoint of BED start and end. Explicit summit columns can be supplied with `--main-summit-column` and `--compare-summit-column`. Scores default to BED column 5 and can be changed with `--main-score-column` and `--compare-score-column`.

## Main-score percentile groups

Matched pairs are sorted by the main BED score and divided into equal-frequency percentile groups. The default:

```text
--percentile-interval 25
```

produces `0-25`, `25-50`, `50-75`, and `75-100` groups. Percentiles are assigned independently for each comparison because different callsets can match different subsets of the main BED.

## Plots

The command produces:

- a combined **signed matched-position distance distribution** for all comparisons. The default display range is **-250 to 250 bp**;
- a combined **score-correlation-by-distance-bin** plot. These bins use absolute summit distance;
- a grouped **main-callset score percentile distance boxplot**, with comparison callsets side-by-side within each percentile group. The default displayed y-axis range is **0-200 bp**. Outliers beyond the 1.5 × IQR whiskers are shown by default; use `--hide-boxplot-outliers` to hide them;
- a separate **main-versus-comparison score agreement** plot for each comparison, coloured by absolute summit distance. The default colour scale is capped at **100 bp**;
- a separate **main-callset score versus matched distance** plot for each comparison. Absolute distance and Spearman correlation are the defaults. The distance axis uses the data range unless the user sets a plotting limit.

Display limits do not remove matched pairs from the underlying tables or statistical calculations. Use `--max-distance` when an actual matching cutoff is required.

## Statistical tests within percentile groups

Add `--stats` to compare the comparison callsets **within each percentile group**. With three comparison callsets, all three pairwise comparisons are tested separately within every group.

The default non-parametric analysis uses a paired Wilcoxon signed-rank test when both comparison distributions contain the same main nucleosome calls. Otherwise it uses a Mann-Whitney U test. `--stats-test parametric` uses a paired t-test or Welch's t-test instead.

Holm multiple-testing correction is applied separately within each percentile group by default. Use:

```text
--p-adjust none
```

for unadjusted p-values.

Plot annotations can show adjusted/raw p-values with:

```text
--p-display value
```

or significance stars with:

```text
--p-display stars
```

The statistics TSV records the test, pairing status, sample sizes, test statistic, raw p-value, adjusted p-value, and significance class.

## Main score versus matched distance

The default relationship is the labelled main-callset peak score versus absolute matched distance:

```text
--score-distance-type absolute
--score-distance-correlation spearman
--score-distance-plot hexbin
```

Use `--score-distance-type signed` to retain upstream/downstream direction, or select Pearson/both correlations with `--score-distance-correlation`.

The corresponding statistics table reports Spearman and Pearson correlations, linear-regression slope/intercept/R-squared, and mean/median absolute distance. `--plot-max-points` limits only the points rendered in scatter/hexbin figures; statistics use all matched pairs.

## Outputs

Default outputs retain the compact tables needed for summaries, statistics, and faithful replotting:

```text
<prefix>_summary.tsv
<prefix>_distance_histogram.tsv
<prefix>_correlation_by_distance.tsv
<prefix>_percentile_summary.tsv
<prefix>_percentile_boxplot.tsv
<prefix>_main_score_vs_distance_statistics.tsv
<prefix>_<comparison>_score_agreement.tsv.gz
<prefix>_<comparison>_main_score_vs_distance.tsv.gz
<prefix>_distance_histogram.png
<prefix>_correlation_by_distance.png
<prefix>_percentile_distance_boxplot.png
<prefix>_<comparison>_score_agreement.png
<prefix>_<comparison>_main_score_vs_distance.png
```

With `--stats` and at least two comparison BEDs, `<prefix>_percentile_statistics.tsv` is also written.

The two compressed per-comparison plot-source tables contain only the sampled points needed to recreate their figures plus the full-data statistics used for annotations. The compact percentile-boxplot source stores box/whisker statistics, actual outliers, comparison labels, and optional statistical annotations rather than every matched pair. All default plot-source tables can be passed directly to [`nucleosuite plot`](plot.md).

Large matched-pair tables are opt-in. Add:

```text
--write-detail-tables
```

to additionally write `<prefix>_percentile_distances.tsv` and each `<prefix>_<comparison>_pairs.tsv`. These files contain one row per matched pair and can be very large for whole-genome callsets.

## Blacklist handling

`--blacklist-bed` removes complete overlapping BED records before matching. The same blacklist is applied to the main BED and every comparison BED.

## Plot customization

Comparison figures use the shared plotting interface described in [Plot customization](../PLOTTING.md). Generated TSV files can also be replotted with [`nucleosuite plot`](plot.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
