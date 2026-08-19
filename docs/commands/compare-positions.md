# `nucleosuite compare-positions`

## What this command does

`compare-positions` compares one main nucleosome callset with one or more comparison callsets. Each comparison is analysed independently against the same main BED. The command matches positions one-to-one, measures the positional distance between matched calls, ranks matched pairs by the main BED score, and divides them into main-score percentile groups.

Quartiles are used by default. When multiple comparison BEDs are supplied, plots are combined where this makes the comparisons easier to interpret, including grouped percentile-distance boxplots and overlaid distance/correlation profiles.

## Why use it

Use it when you want to compare several nucleosome callers, datasets, or parameter settings against one common reference callset while keeping the main callset as the basis for score stratification.

## Typical use

```bash
nucleosuite compare-positions \
  --main-bed PNS_nucleosomes.bed \
  --compare-bed iNPS=iNPS_nucleosomes.bed \
  --compare-bed DANPOS=DANPOS_nucleosomes.bed \
  --main-score-column 5 \
  --compare-score-column 5 \
  --output-prefix Gaffney32_position_compare
```

The label before `=` is optional. Without it, the comparison BED basename is used.

## One search per comparison

For each comparison BED, the command performs one nearest-position search:

1. the main BED and that comparison BED are counted after filtering;
2. the BED with fewer positions is used as the query set;
3. the BED with more positions is used as the target set;
4. one-to-one unique matching is applied, so a position can participate in at most one matched pair;
5. matched pairs always retain the main BED call and main BED score regardless of which input was the query set.

This means a comparison with fewer positions than the main BED searches comparison → main, while a larger comparison callset searches main → comparison. There is no reciprocal second search.

## Summit and score columns

By default, BED interval midpoints are used as positions and column 5 is used as the score.

Use:

```text
--main-summit-column
--compare-summit-column
--main-score-column
--compare-score-column
```

when explicit summit or score columns are required. The comparison summit and score column settings are applied to all `--compare-bed` inputs.

## Main-score percentile groups

After one-to-one matching is complete for a comparison, the matched pairs are sorted by the **main BED score**. They are then divided into equal-frequency percentile groups.

The default:

```text
--percentile-interval 25
```

produces:

```text
0-25
25-50
50-75
75-100
```

The percentile assignment is calculated independently for each main-versus-comparison match set because different comparison BEDs can match different subsets of the main callset.

The grouped percentile-distance boxplot places all comparison callsets side-by-side within each main-score percentile group.

## Statistical comparisons within percentile groups

Statistics are optional:

```bash
nucleosuite compare-positions \
  --main-bed PNS.bed \
  --compare-bed iNPS=iNPS.bed \
  --compare-bed DANPOS=DANPOS.bed \
  --compare-bed NucPos=NucPos.bed \
  --stats \
  --p-display value
```

Pairwise tests are performed **separately within each percentile group**. For three comparison callsets, each percentile group therefore tests:

```text
iNPS vs DANPOS
iNPS vs NucPos
DANPOS vs NucPos
```

The default statistical family is non-parametric:

```text
--stats-test nonparametric
```

When all observations in two comparison distributions within a percentile group can be aligned to the same main nucleosome calls, the observations are paired and a two-sided Wilcoxon signed-rank test is used. If complete main-call pairing is not available, a two-sided Mann-Whitney U test uses the full two distributions.

With:

```text
--stats-test parametric
```

paired observations use a paired t-test and unpaired observations use Welch's t-test.

Multiple-testing correction is applied independently within each percentile group using Holm correction by default:

```text
--p-adjust holm
```

Use `--p-adjust none` for raw p-values.

Plot annotations can show p-values:

```text
--p-display value
```

or significance stars:

```text
--p-display stars
```

using:

```text
ns      p >= 0.05
*       p < 0.05
**      p < 0.01
***     p < 0.001
****    p < 0.0001
```

When Holm correction is enabled, the adjusted p-value is used for the displayed value or significance class.

The statistics TSV reports the percentile group, comparison pair, test used, whether pairing was possible, sample counts, test statistic, raw p-value, adjusted p-value, and significance class.

## Main peak score versus matched distance

Each comparison also receives a separate plot of main BED peak score against the distance to the matched comparison call.

The default uses absolute distance:

```text
--score-distance-type absolute
```

Use signed distance to retain upstream/downstream direction:

```text
--score-distance-type signed
```

The default correlation is Spearman:

```text
--score-distance-correlation spearman
```

Pearson or both can be displayed with `pearson` or `both`.

The statistics table for these plots reports, for each comparison:

- matched-pair count;
- Spearman rho and p-value;
- Pearson r and p-value;
- linear-regression slope and intercept;
- linear R-squared and slope p-value;
- median and mean absolute matched distance.

The default plot rendering is a hexbin density representation so large callsets remain readable:

```text
--score-distance-plot hexbin
```

Use `--score-distance-plot scatter` for individual points. `--plot-max-points` limits only the number of points drawn; correlation and summary statistics are calculated from all matched pairs.

## Other comparison plots

The command retains the existing position/score comparison concepts while combining multiple comparisons where practical:

- **distance distribution** — all comparison callsets are overlaid in one plot;
- **score correlation by distance bin** — comparison trajectories are overlaid;
- **main-score percentile distance boxplot** — comparison boxes are side-by-side within each percentile group;
- **main-versus-comparison score agreement** — written separately for each comparison because overlaid scatter plots would obscure the relationships;
- **main score versus matched distance** — written separately for each comparison.

## Outputs

The combined outputs include:

```text
<prefix>_summary.tsv
<prefix>_distance_histogram.tsv
<prefix>_correlation_by_distance.tsv
<prefix>_percentile_distances.tsv
<prefix>_percentile_summary.tsv
<prefix>_main_score_vs_distance_statistics.tsv
<prefix>_distance_histogram.png
<prefix>_correlation_by_distance.png
<prefix>_percentile_distance_boxplot.png
```

With `--stats` and at least two comparison BEDs:

```text
<prefix>_percentile_statistics.tsv
```

Each comparison also receives:

```text
<prefix>_<comparison>_pairs.tsv
<prefix>_<comparison>_score_agreement.png
<prefix>_<comparison>_main_score_vs_distance.png
```

Use `--skip-pairs-tsv` to suppress the large detailed matched-pair tables.

Plot metadata sidecars record the complete command invocation and parameters.

## Blacklist handling

`--blacklist-bed` removes complete overlapping BED records before matching. The same blacklist is applied to the main BED and all comparison BEDs.

## Plot customization

Comparison figures use the shared plotting interface described in [Plot customization](../PLOTTING.md). The generated TSV files can also be replotted with [`nucleosuite plot`](plot.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
