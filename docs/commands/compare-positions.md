# `nucleosuite compare-positions`

## What this command does

`compare-positions` compares one main nucleosome callset with one or more comparison callsets using one-to-one positional matching. It can also sample one or more BigWig tracks at the main-callset summit positions and use those values as score-only comparators. Matched BED pairs are ranked by the **main-callset score** and divided into percentile groups; quartiles are used by default.

## Why use it

Use it to compare several nucleosome callsets against one common reference while retaining the main-callset score for stratification, grouped summaries and optional within-percentile tests. `--score-bigwig` is useful when you also want to ask how a signal such as coverage, accessibility, or another continuous track correlates with the main peak scores without treating that track as a positional callset.

## Typical use

```bash
nucleosuite compare-positions \
  --main-bed PNS=PNS_nucleosomes.bed \
  --compare-bed iNPS=iNPS_nucleosomes.bed \
  --compare-bed DANPOS=DANPOS_nucleosomes.bed \
  --compare-bed WPS=WPS_nucleosomes.bed \
  --score-bigwig Coverage=coverage.bw \
  --output-prefix Gaffney32_position_compare
```

Both `--main-bed` and `--compare-bed` accept `LABEL=path.bed`. `--score-bigwig` similarly accepts `LABEL=track.bw` and can be repeated. The supplied labels are used in plot axes, legends and output tables. If no label is supplied, the filename is used.

## BigWig score comparators

Each `--score-bigwig` track is sampled at the **retained main-callset summit coordinate**. The summit is taken from `--main-summit-column` when supplied; otherwise the BED interval midpoint is used. The BigWig value at `[summit, summit + 1)` is paired directly with the main BED score at that same position.

BigWig comparators are **score-only**. They are included in the per-comparator score-agreement plots and score-correlation summary statistics, but are deliberately excluded from:

- positional matching;
- signed/absolute distance distributions;
- score-correlation-by-distance bins;
- percentile distance boxplots;
- percentile distance trends; and
- within-percentile distance tests.

Missing or non-finite BigWig values, absent BigWig contigs, and out-of-range summit coordinates are represented as zero for score comparison. The summary table reports how many sampled values were converted to zero by this rule.

Example:

```bash
nucleosuite compare-positions \
  --main-bed PNS=PNS_nucleosomes.bed \
  --compare-bed iNPS=iNPS_nucleosomes.bed \
  --score-bigwig Coverage=coverage.bw \
  --score-bigwig PNS_signal=PNS.bw
```

This produces ordinary positional outputs for `iNPS`, plus independent score-agreement plots for `Coverage` and `PNS_signal` against the PNS peak scores.

## Matching

Each main-versus-comparison pair is searched once. The callset with fewer positions is used as the query and the larger callset as the target. Matching is one-to-one, so a position can occur in at most one accepted pair.

Regardless of search direction, every accepted pair is represented as a main call and a comparison call. Signed distance is:

```text
comparison summit - main summit
```

Use `--max-distance` to reject pairs beyond a specified absolute summit distance.

By default, position is the integer midpoint of BED start and end. Explicit summit columns can be supplied with `--main-summit-column` and `--compare-summit-column`. Scores default to BED column 5 and can be changed with `--main-score-column` and `--compare-score-column`.

## Main-score percentile groups

Matched pairs are sorted by the main-callset score and divided into equal-frequency percentile groups. The default:

```text
--percentile-interval 25
```

produces `0-25`, `25-50`, `50-75`, and `75-100` groups. Percentiles are assigned independently for each comparison because different comparison callsets can match different subsets of the main BED.

## Plots

The default figures are:

- a combined **signed matched-position distance distribution** for all comparisons, displayed from **-250 to 250 bp** by default;
- a combined **score-correlation-by-distance-bin** plot using absolute summit distance bins;
- a grouped **main-score percentile distance boxplot**, with comparison callsets side-by-side within each percentile group and a default display range of **0-200 bp**;
- a combined **1%-percentile distance trend**, showing median absolute matched distance for each comparison with an interquartile-range ribbon;
- a separate **main-versus-comparison score-agreement** plot for each comparison BED, coloured by absolute summit distance with a default colour scale of **0-100 bp**; and
- a separate **main-score-versus-BigWig-value agreement** plot for each `--score-bigwig`, with no distance colouring because no positional matching is involved.

Boxplot outliers and statistical annotations are hidden by default. Display limits affect figures only; use `--max-distance` when an actual matching cutoff is required.

## Statistical tests within percentile groups

Add `--stats` to compare comparison callsets **within each percentile group**. With three comparison callsets, all three pairwise comparisons are tested separately within every group.

The default non-parametric analysis uses a paired Wilcoxon signed-rank test when both comparison distributions contain the same main nucleosome calls. Otherwise it uses a Mann-Whitney U test. `--stats-test parametric` uses a paired t-test or Welch's t-test instead.

Holm correction is applied separately within each percentile group by default. `--p-adjust none` disables the correction. Plot annotations can show p-values or significance stars with `--p-display value` or `--p-display stars`.

## Outputs

Compact summary and plot-source tables are written by default:

```text
<prefix>_summary.tsv
<prefix>_distance_histogram.tsv
<prefix>_correlation_by_distance.tsv
<prefix>_percentile_summary.tsv
<prefix>_percentile_boxplot.tsv
<prefix>_percentile_distance_trend.tsv
<prefix>_<comparison>_score_agreement.tsv.gz
```

The score-agreement source is produced for both BED and BigWig comparators. BigWig score sources store only the deterministic plotted sample plus full-data correlation statistics, so they remain compact even when millions of main positions are sampled.

Corresponding figures include:

```text
<prefix>_distance_histogram.png
<prefix>_correlation_by_distance.png
<prefix>_percentile_distance_boxplot.png
<prefix>_percentile_distance_trend.png
<prefix>_<comparison>_score_agreement.png
```

With `--stats` and at least two comparison BEDs, `<prefix>_percentile_statistics.tsv` is also written.

The compact boxplot source stores box/whisker statistics and the values needed to reproduce optional outliers/statistical annotations. The 1%-trend source stores comparison, percentile, observation count, median absolute distance, and the 25th/75th percentiles. These sources can be passed directly to [`nucleosuite plot`](plot.md).

Large one-row-per-match tables are opt-in. Add:

```text
--write-detail-tables
```

to additionally write `<prefix>_percentile_distances.tsv`, each `<prefix>_<comparison>_pairs.tsv` for BED comparisons, and `<prefix>_<label>_score_bigwig_values.tsv` for each BigWig comparator. The latter contains one row per retained main position with chromosome, summit, main score and sampled BigWig value.

## Blacklist handling

`--blacklist-bed` removes complete overlapping BED records before matching. The same blacklist is applied to the main BED and every comparison BED. BigWig score sampling is then performed only at the retained main positions, so main records removed by the blacklist are not sampled.

## Plot customization

Comparison figures use the shared plotting interface described in [Plot customization](../PLOTTING.md). Generated compact plot-source TSV files can be replotted with [`nucleosuite plot`](plot.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
