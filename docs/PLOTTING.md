# Plot customization

NucleoSuite uses shared plot options across commands that generate figures. Each command supplies default settings; specify only the values you want to change.

## Command-line help

Analysis commands keep the shared plot-customization block out of ordinary help. Use:

```bash
nucleosuite COMMAND --help
```

for the analysis options, and:

```bash
nucleosuite COMMAND --help-plotting
```

to expand the plotting controls for that command. The dedicated `nucleosuite plot --help` command shows its plotting options directly.


## Replot an existing NucleoSuite result

Use [`nucleosuite plot`](commands/plot.md) when the analysis has already finished and you want deeper control over the figure without recalculating the underlying data:

```bash
nucleosuite plot sample_dac.tsv
```

The replot command auto-detects NucleoSuite output tables and reproduces each source command's default figure style and layout before applying user overrides. Distance distributions retain their solid, order-coloured lines and legend; distance and standalone NRL regressions use open circles, a dotted fit, and a square figure. It also provides independent major/minor tick and grid controls, heatmap saturation limits, optional DAC peak detection and NRL regression insets, generic line/scatter/bar/heatmap modes, and Matplotlib `rcParams`/artist keyword pass-through. See the [plot command page](commands/plot.md) for the complete interface.

For [`compare-positions`](commands/compare-positions.md), multi-callset outputs retain comparison labels in the plot tables. The combined distance distribution uses signed summit distance; correlation-by-distance bins and percentile-distance summaries use absolute distance. Replotting preserves the multiple comparison series/grouped boxes. Per-comparison matched-pair tables can be replotted as score-agreement figures, or explicitly as main-score-versus-distance figures with `--plot-type compare-positions-score-distance`.

## Choose PNG or SVG

```text
--plot-format {png,svg}
```

PNG is the default. SVG provides a vector figure for editing or publication.

```bash
nucleosuite distances ... --plot-format svg
```

## Set figure size and resolution

```text
--plot-width FLOAT
--plot-height FLOAT
--plot-dpi INT
--plot-transparent
```

Width and height are in inches. DPI mainly matters for PNG; SVG remains vector-based.

Example:

```bash
--plot-format svg --plot-width 10 --plot-height 5
```

## Change or remove titles and axis labels

```text
--plot-title TEXT
--no-plot-title
--plot-x-label TEXT
--plot-y-label TEXT
--plot-font-size FLOAT
```

Use `--plot-title` when the automatic title is not suitable for a figure panel. Use `--no-plot-title` when titles will be added later in a manuscript or presentation.

## Control tick labels

```text
--plot-x-tick-rotation FLOAT
--plot-y-tick-rotation FLOAT
```

These are useful when category labels or sample names overlap.

## Control grid lines

Grid lines are the lines extending from tick marks across the plotting area.

```text
--plot-grid {none,x,y,both}
--plot-grid-color COLOR
--plot-grid-alpha FLOAT
--plot-grid-width FLOAT
```

Use:

```bash
--plot-grid none
```

when you want no grid lines inside the plot. Use `x`, `y`, or `both` only when tick-aligned reference lines improve interpretation.

Color accepts Matplotlib color names or values such as hexadecimal colours. Alpha ranges from 0 to 1.

## Set displayed axis limits

```text
--plot-x-min FLOAT
--plot-x-max FLOAT
--plot-y-min FLOAT
--plot-y-max FLOAT
```

These options change the displayed range. Data-filtering options are documented separately on each command page.

## Change lines and fills

```text
--plot-line-width FLOAT
--plot-line-color COLOR
--plot-fill-color COLOR
```

For a single scientific series, the line/fill colour can be overridden directly. Multi-series figures keep their command-defined palettes by default so biological categories such as ChromHMM states or dinucleotide classes remain distinguishable.

## Show point markers on line plots

```text
--plot-points
--no-plot-points
--plot-point-size FLOAT
--plot-point-fill COLOR
--plot-point-edge COLOR
--plot-point-edge-width FLOAT
--plot-point-shape {circle,square,triangle,diamond}
```

Point markers and labels are controlled independently.

Example:

```bash
--plot-points \
--plot-point-shape circle \
--plot-point-size 4 \
--plot-point-fill white \
--plot-point-edge black
```

## Label selected points or called peaks

```text
--plot-label-points {none,peaks,all}
--plot-point-label-value {x,y,both}
--plot-point-label-offset FLOAT
```

Labels are placed immediately above the selected point, centred horizontally, and rotated vertically so sparse annotations use little horizontal space.

The default label content is the **x-axis value**.

- `none` — show no point labels.
- `peaks` — label only meaningful peak/local-maximum positions supplied by the plotting command.
- `all` — label every eligible point; intended for sparse plots.

### DAC and NRL defaults

DAC plots do not detect or label peaks by default. They show the raw DAC profile only. Enable DAC peak labels with:

```bash
--plot-label-points peaks
```

When DAC peak labels are enabled, DAC uses the same resolution-driven peak detection and smoothing as `nucleosuite nrl`, and the raw profile remains visible in grey behind the smoothed curves. NRL plots label their called peaks by default.

Change the text to y values or both coordinates with:

```bash
--plot-point-label-value y
--plot-point-label-value both
```

## Legends

```text
--plot-legend
--no-plot-legend
--plot-legend-position {best,upper-right,upper-left,lower-right,lower-left,outside-right}
```

Use `outside-right` when a legend would cover data inside the plotting area.

## A clean publication-style example

```bash
nucleosuite dac \
  --bigwig sample_dyad.bw \
  --chrom-sizes sample.bam \
  --scope combined_chromosomes \
  --out-prefix sample_dac \
  --plot-format svg \
  --plot-width 10 \
  --plot-height 5 \
  --plot-grid none \
  --plot-line-width 1.5 \
  --plot-points \
  --plot-point-shape circle \
  --plot-point-size 3 \
  --plot-point-fill white \
  --plot-point-edge black \
  --plot-point-label-value x
```

## Suite commands

The same shared options can be supplied to `cfdna-suite` and `mnase-suite`. The suite passes explicitly supplied plot settings to downstream plot-producing commands.

For example:

```bash
nucleosuite cfdna-suite ... \
  --plot-format svg \
  --plot-grid none \
  --no-plot-title
```

## Command-specific plot controls

Some commands provide additional controls for the scientific layout:

- `peak-states --plot-bar-gap` and `--plot-x-axis`;
- fragment-heatmap normalization/palette choices;
- distance-specific tick and state-overlay settings.
- `flank-spacing --top-categories`, `--ratio-x1`, `--ratio-x2`, and `--distribution`; highlighted lines are layered by rank so rank 1 is drawn on top.

Those options are explained on the relevant command page and work alongside the shared options above.

## Plot metadata sidecars

Every generated plot is accompanied by a `*_metadata.tsv` sidecar. The sidecar records the NucleoSuite version, plot path, complete command invocation, and the full parsed parameter set. Automatic filenames are intentionally limited to at most three central analysis-parameter tokens; use the metadata sidecar for the complete provenance of a figure.

