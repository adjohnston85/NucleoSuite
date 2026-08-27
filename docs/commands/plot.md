# `nucleosuite plot`

## What this command does

`plot` recreates a NucleoSuite figure from its compact plot-source table and metadata sidecar without rerunning the underlying genomic analysis.

## Why use it

Use it to change figure presentation or plot-specific construction settings after an analysis has finished, without recalculating the genomic result.

## Metadata-driven replots

The simplest form is:

```bash
nucleosuite plot sample_output.tsv
```

For NucleoSuite-generated plot sources, the associated metadata records the plot family and the parameters used to construct the original figure. `plot` reads those values first and uses them as the replot defaults. A single source table may be associated with more than one native figure; when that happens, `plot` recreates all applicable figures by default.

If the original metadata file is edited, the next replot uses the edited values. Command-line arguments supplied to `nucleosuite plot` override matching metadata values for that run. Replotting never rewrites the original analysis metadata. Each recreated figure instead receives its own `_replot_metadata.tsv` sidecar, which can be replaced by later replots.

The parser is two-stage: universal figure controls are always available, while plot-specific controls are added only after the source table and metadata identify the plot family. This keeps analysis-specific options out of unrelated replots.

For example, a percentile-boxplot source may expose an outlier control, while a DAC source may expose peak-detection controls. Those options are available only for the corresponding plot type.

## Automatic detection

`plot` recognises the major NucleoSuite plot families, including:

- DAC and DCC profiles;
- NRL profiles and regression-point tables;
- distance distributions, state overlays, percentile curves/tables, and NRL regressions;
- aggregate mean profiles and heatmap matrices;
- fragment-length profiles and normalized fragment heatmap matrices;
- positive-run distributions;
- peak-score frequency distributions;
- flank-spacing category distributions;
- peak-state stacked compositions;
- compare-positions signed distance distributions, BED score-agreement plots, BigWig score-only agreement plots, distance-bin correlations, grouped percentile boxplots, and 1%-percentile median/IQR trends;
- gene-expression spacing/FFT outputs;
- TSS-expression-quintile profiles;
- dinucleotide and WW/SS profiles;
- gene-set candidate-overlap plots;
- multi-profile overlays; and
- fragment-relocation/count profiles.

### Multiple plots from one source table

Some compact outputs contain the data for multiple figures. For example, a dinucleotide-profile TSV contains both the 16 individual dinucleotides and the derived WW/SS profile. Therefore:

```bash
nucleosuite plot sample_dinuc_profile.tsv
```

recreates both the dinucleotide figure and the WW/SS figure.

List the figures available from an input without rendering them:

```bash
nucleosuite plot sample_dinuc_profile.tsv --list-plots
```

Render only a subset with `--plots`. The option accepts a comma-separated list and can also be repeated:

```bash
nucleosuite plot sample_dinuc_profile.tsv --plots dinuc
nucleosuite plot sample_dinuc_profile.tsv --plots ww-ss
nucleosuite plot sample_dinuc_profile.tsv --plots dinuc,ww-ss
```

When several plots are recreated, each receives a distinct automatic output name. If `--output` is supplied, it applies to the primary plot and additional plots retain their automatic names. Use `--plots` to select a single figure when one explicit output name is required.

When a table is unusual or has been renamed, specify the source command or exact plot family:

```bash
nucleosuite plot output.tsv --from-command distances
```

```bash
nucleosuite plot output.tsv --plot-type generic-line \
  --x-column distance_bp \
  --y-column count
```

## Reconstruct zero-count distance positions

Distance tables produced by older NucleoSuite versions may omit distances with a raw count of zero. For a distance source, `plot` can reconstruct those missing integer-bp positions in memory before drawing the raw profile:

```bash
nucleosuite plot sample_distances.tsv --include-zero-distances
```

The source TSV and its original metadata are not modified. The reconstruction is limited to the displayed distance range (`--x-min` / `--x-max`, normally restored from the source metadata). Use `--no-include-zero-distances` when the sparse observed-only representation is desired. New `distances` outputs include zero-count positions by default.

## Major and minor ticks

Base-pair plots use readable multiples of 10 where practical. Set exact tick spacing with:

```text
--x-major-tick FLOAT
--x-minor-tick FLOAT
--y-major-tick FLOAT
--y-minor-tick FLOAT
```

Minor tick labels are hidden. For example, a DCC plot with labelled 10 bp ticks and unlabelled 5 bp ticks can be written with:

```bash
nucleosuite plot sample_dcc.tsv \
  --x-major-tick 10 \
  --x-minor-tick 5
```

## Major and minor grid lines

Each tick family can control its own grid lines:

```text
--x-major-grid / --no-x-major-grid
--x-minor-grid / --no-x-minor-grid
--y-major-grid / --no-y-major-grid
--y-minor-grid / --no-y-minor-grid
```

Their styles can be tuned independently:

```text
--major-grid-color COLOR
--minor-grid-color COLOR
--major-grid-alpha FLOAT
--minor-grid-alpha FLOAT
--major-grid-width FLOAT
--minor-grid-width FLOAT
--major-grid-style STYLE
--minor-grid-style STYLE
```

The compact signed-lag DCC default uses solid 10 bp major guides and dashed 5 bp minor guides. Other plot families use clean axes unless grid lines are requested.

## DAC publication-style replots

DAC replots show the raw DAC profile only by default. No smoothing, peak detection, peak markers, or peak labels are added unless peak detection is explicitly requested.

Enable peak detection with `--detect-peaks`. DAC then uses the same resolution-driven smoothing and peak-calling method as the standalone `nrl` command. `--peak-resolution` controls the minimum peak spacing and both smoothing scales; the default 160 bp resolution gives 61 bp detection smoothing and 21 bp local-maximum refinement. The raw profile is drawn in grey behind those two smoothed profiles.

```bash
nucleosuite plot sample_dac.tsv \
  --detect-peaks \
  --peak-resolution 160 \
  --label-peaks peaks \
  --nrl-inset on
```

Peak labels are centred directly over the called peak position and placed above the marker. Use `--label-peaks none` to keep the NRL-style smoothing and called-peak markers without text labels. NRL profile replots label their called peaks by default.

Use `--nrl-inset off` to remove the inset or adjust its position with:

```text
--inset-bounds X Y WIDTH HEIGHT
```

## Heatmap saturation and ticks

Heatmap matrices can be rescaled without modifying their values:

```bash
nucleosuite plot sample_heatmap_matrix.tsv.gz \
  --vmin -10 \
  --vmax 10 \
  --x-major-tick 100 \
  --x-minor-tick 10
```

Values below `--vmin` or above `--vmax` saturate at the ends of the colour scale. This is useful when a small number of extreme values would otherwise compress most of the heatmap into a narrow colour range.

Aggregate heatmap matrices are written when `aggregate --write-detail-tables` is requested and can then be replotted directly. Fragment-heatmap `_normalised_matrix.tsv` outputs remain default because the matrix is the compact source needed to reproduce that heatmap. Fragment-heatmap runs also write `_heatmap_plot_metadata.tsv` and `_heatmap_linkage.tsv`; when these sit beside the matrix, `plot` restores the original dendrogram, category strip and legend, palette, colour limits, labels, and layout.

## Compact plot-source tables

Large per-region or per-match detail tables are disabled by default in commands that can generate very large supporting outputs. Plot-producing analyses retain a smaller source table containing the information required to reproduce each default figure. For example, `compare-positions` writes `_percentile_boxplot.tsv` plus compressed per-comparison score plot sources even when the full matched-pair tables are omitted.

These source tables can be passed directly to `nucleosuite plot`. Plot metadata sidecars record the associated source table and detected plot family where applicable.

For compare-position percentile boxplots, outliers beyond the 1.5 × IQR whiskers are hidden by default. After the source and metadata identify a boxplot, they can be enabled or disabled with:

```text
--show-boxplot-outliers
--hide-boxplot-outliers
```

This is a plot-specific control; `showfliers` is not a Matplotlib `rcParam`.

## General figure controls

Distance-distribution replots reproduce the analysis output: each neighbour order is a differently coloured solid line and the legend identifies the order. Distance-command and standalone NRL regression replots use open circles with a dotted fit and default to a square 6.5 × 6.5 inch figure. Use `--mpl-kw` to change the line or marker appearance. Explicit `--width` or `--height` values override the corresponding default dimension.

Fragment-length tables are source-aware. Auto mode recreates the within-window density used by `fragment-lengths --plot`, including its 0 bp start and observed-length upper bound capped at 1000 bp. `.fragment_length_counts.tsv` files from `fragments` recreate that command's raw-count figure. Override normalization only when needed with `--normalization count` or `--normalization density`.

The fragment-size NRL profile and peak tables are directly re-plottable. The profile table recreates the raw, local-refinement and broad-detection layers with the called summits. The peaks table recreates the square regression with open circles and a dotted fit:

```bash
nucleosuite plot sample_fragment_lengths_fragment_size_nrl_profile.tsv
nucleosuite plot sample_fragment_lengths_fragment_size_nrl_peaks.tsv
```

Aggregate NRL outputs are also source-aware. The aggregate profile TSV recreates the one unified negative-to-positive peak-calling plot, including smoothing across position 0 and any shaded regression exclusion interval. Each directional regression TSV recreates its own square open-circle/dotted-line fit, including an order-0 shared central peak and retained peak-number gaps when present:

```bash
nucleosuite plot sample_aggregate_nrl_profile.tsv
nucleosuite plot sample_aggregate_nrl_positive_regression.tsv
nucleosuite plot sample_aggregate_nrl_negative_regression.tsv
```

```text
--format {png,svg,pdf}
--width FLOAT
--height FLOAT
--dpi INT
--title TEXT
--no-title
--x-label TEXT
--y-label TEXT
--font-size FLOAT
--x-min FLOAT
--x-max FLOAT
--y-min FLOAT
--y-max FLOAT
--x-tick-rotation FLOAT
--y-tick-rotation FLOAT
--axes-facecolor COLOR
--transparent
--no-legend
```

For a source with one applicable figure, the default output is `<input_stem>_replot.png`. Multi-plot sources receive distinct names derived from the corresponding original figure or plot family. Use `--output` to choose the primary output name or `--plots` to select one figure explicitly.

## Direct Matplotlib customization

For settings that do not need dedicated NucleoSuite options, `plot` exposes two repeatable pass-through mechanisms.

Set any valid Matplotlib `rcParams` entry with:

```bash
--mpl-rc axes.linewidth=1.2 \
--mpl-rc font.family=Arial
```

Pass keyword arguments to figure artists with:

```bash
--mpl-kw raw.color=0.7 \
--mpl-kw smooth.linewidth=2.0 \
--mpl-kw points.s=25 \
--mpl-kw legend.loc='upper left'
```

Supported artist targets are `line`, `raw`, `smooth`, `points`, `bar`, `heatmap`, and `legend`. Values are parsed as Python literals when possible; ordinary strings are passed through unchanged.

For example, change a distance-regression replot from its default dotted fit and open circles to a solid fit and filled points:

```bash
nucleosuite plot distance_nrl_regression.tsv \
  --mpl-kw line.linestyle=- \
  --mpl-kw points.facecolors=black
```

This gives access to Matplotlib properties such as line styles, marker properties, colormaps, image interpolation, bar styling, and legend layout without requiring a dedicated NucleoSuite flag for every Matplotlib keyword.

## Generic tables

A tabular file that is not recognised can still be plotted:

```bash
nucleosuite plot custom.tsv \
  --plot-type generic-scatter \
  --x-column x_value \
  --y-column y_value
```

Available generic families are `generic-line`, `generic-scatter`, `generic-bar`, and `generic-heatmap`. Use `--group-column` to overlay several line groups when appropriate.

[Back to the command reference](../COMMAND_REFERENCE.md)
