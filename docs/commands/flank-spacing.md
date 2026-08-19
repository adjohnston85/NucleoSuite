# `nucleosuite flank-spacing`

## What this command does

`flank-spacing` measures the distance between the nearest nucleosome centre strictly upstream of each BED reference site and the nearest nucleosome centre strictly downstream. Reference sites are grouped by a category column, column 4 by default, and one spacing distribution is calculated per category.

By default, the command plots kernel-density curves. Raw-count curves are available with `--distribution count`.

Each category is ranked by the ratio of the distribution height at two user-selected spacing positions. The defaults are 190 bp and 260 bp, so the default ranking statistic is:

```text
height at 190 bp / height at 260 bp
```

The lowest ratio is rank 1. This prioritises categories with a relatively stronger widened-spacing component at 260 bp.

## Why use it

Use this command when a set of reference sites is divided into biological or technical categories and you want to compare how the nucleosomes immediately flanking those sites are spaced.

The categories can represent cell lines, chromatin classes, transcription-factor site classes, experimental conditions, annotations, or any other labels stored in the reference BED. DNase-hypersensitive-site callsets are one example, analogous to the comparison in Snyder et al. (2016) Figure 5A, but the command is not specific to DHS data.

## Typical command

```bash
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosomes.bed \
  --region-bed categorized_sites.bed \
  --category-col 4 \
  --output-prefix sample_flank_spacing
```

For each reference site, NucleoSuite uses the interval midpoint unless `--point-col` supplies an exact coordinate. Nucleosome centres are likewise taken from interval midpoints unless `--nucleosome-center-col` is supplied. A standalone example is provided in [`examples/flank_spacing.sh`](../../examples/flank_spacing.sh).

## Density versus raw counts

Density is the default:

```bash
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosomes.bed \
  --region-bed categorized_sites.bed
```

To plot raw spacing counts instead:

```bash
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosomes.bed \
  --region-bed categorized_sites.bed \
  --distribution count
```

Density curves are estimated using all valid spacing observations in the category. `--x-max` controls the displayed/output grid range; it does not discard wider spacing observations before density estimation.

## Ranking categories

The two evaluation positions can be changed independently:

```bash
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosomes.bed \
  --region-bed categorized_sites.bed \
  --ratio-x1 185 \
  --ratio-x2 270
```

The statistic is always `y(x1) / y(x2)`, and categories are ranked from the lowest ratio to the highest. If the denominator is zero and the numerator is positive, the ratio is infinite and sorts after finite ratios. Categories for which both evaluated values are zero sort last.

## Highlighting the top categories

By default, the top seven ranked categories receive distinct colours and appear in the legend. All remaining categories are grey.

```bash
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosomes.bed \
  --region-bed categorized_sites.bed \
  --top-categories 10
```

Grey background curves are drawn first. Highlighted curves are layered by rank so rank 1 is drawn last and therefore remains visible above every other curve. The legend is ordered rank 1, rank 2, and so on.

## Plot range

The default displayed range is 0-500 bp. Change it with:

```bash
--x-min 0 --x-max 600
```

Both ratio positions must lie inside this range.

## Outputs

The command writes:

- `*_sites.tsv`: one row per input reference site, including the selected upstream/downstream nucleosome centres and resulting spacing;
- `*_distributions.tsv`: the plotted density or count value at each spacing position for every category;
- `*_ranking.tsv`: category rank, site counts, values at the two ratio positions, and the ratio;
- a PNG or SVG distribution plot; and
- a plot metadata TSV containing the complete command parameter set.

Automatic output names contain no more than three central analysis-parameter tokens. The complete parameter set is retained in the plot metadata sidecar.

## Plot customization

The plot accepts the shared options described in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
