# Choosing a NucleoSuite command

Choose a command from your input and the result you need. Most analyses begin with paired-end fragments or a genomic signal or peak file.

## I have a paired-end BAM

Use [`fragments`](commands/fragments.md) if you first want an explicit fragment BED file. Use [`tracks`](commands/tracks.md) when you want several fragment-derived tracks from one pass through the data.

For nucleosome-oriented signal:

- [`pns`](commands/pns.md) — score positions with endpoint probability triangles (PNS), a balanced central boxcar (BNS), or a centred unit-mass triangle (TNS).
- [`wps`](commands/wps.md) — score whether fragments protect a fixed window or terminate inside it.
- [`dyads`](commands/dyads.md) — place signal at fragment centres.
- [`coverage`](commands/coverage.md) — count how many fragments cover each base.
- [`fragment-ends`](commands/fragment-ends.md) — count fragment starts and ends.

For a coordinated analysis with standard defaults, use [`cfdna-suite`](commands/cfdna-suite.md) or [`mnase-suite`](commands/mnase-suite.md).

## I already have nucleosome or other peak calls

Use [`distances`](commands/distances.md) to measure nearest-neighbour or higher-order spacing. Use [`compare-positions`](commands/compare-positions.md) to compare two callsets. Use [`peak-score-frequency`](commands/peak-score-frequency.md) to compare their score distributions.

Use [`peak-states`](commands/peak-states.md) to measure how peaks are distributed across chromatin states. Use [`filter-coverage`](commands/filter-coverage.md) when you want to remove peaks that do not meet a minimum BigWig coverage value at their summit or interval midpoint.

## I have a genomic signal and want periodicity

Use [`dac`](commands/dac.md) to measure whether **the same signal repeats at characteristic distances**. A dyad signal with regular nucleosome spacing is a common input.

Use [`dcc`](commands/dcc.md) to measure whether **one signal occurs at a characteristic offset from another**.

Use [`nrl`](commands/nrl.md) on a DAC or DCC distance profile when you want one recurring-period estimate from its repeated peaks.

## I want signal around genomic features

Use [`aggregate`](commands/aggregate.md) for a heatmap and average profile around reference sites. Use [`region-extract`](commands/region-extract.md) to export the signal vector for every region.

A bundled resource can be passed directly. For example, to aggregate signal around the bundled GM12878 CTCF sites:

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --out-prefix sample_ctcf
```

## I want chromatin-state analyses

The bundled GM12878 hg19 ChromHMM states can be used directly:

```bash
STATES="$(nucleosuite resources path gm12878-hg19-states)"
```

Then pass `$STATES` to commands that accept a state or region BED, for example:

```bash
nucleosuite distances peaks.bed --state-bed "$STATES" --output-prefix distances_by_state
```

or:

```bash
nucleosuite peak-states peaks.bed --state-bed "$STATES" --output-prefix peak_states
```

## I want fragment-length or sequence composition analyses

Use [`fragment-lengths`](commands/fragment-lengths.md) to count lengths and [`fragment-heatmap`](commands/fragment-heatmap.md) to compare those profiles across samples or region classes. Use [`dinuc-profile`](commands/dinuc-profile.md) for positional dinucleotide frequencies and [`ww-types`](commands/ww-types.md) for the centred WW/SS fragment classification.

## I want gene-centred analyses

Use [`gene-sets`](commands/gene-sets.md) to classify genes from chromatin-state overlaps, [`gene-expression`](commands/gene-expression.md) to relate expression to spacing or periodicity, and [`tss-expression-quintiles`](commands/tss-expression-quintiles.md) to compare signal around TSSs across expression groups.

The bundled genes and expression table can be addressed directly:

```bash
GENES="$(nucleosuite resources path hg19-genes)"
EXPR="$(nucleosuite resources path hpa-tissue-expression)"
```


## I already have NucleoSuite output and want a different figure

Use [`plot`](commands/plot.md) to recreate a figure from a NucleoSuite TSV or TSV.GZ without rerunning the underlying analysis. The command auto-detects major output families and provides independent major/minor tick and grid controls, heatmap saturation limits, DAC smoothing/peak labels/NRL insets, and Matplotlib pass-through options.

## I need a bundled annotation file

Use [`resources`](commands/resources.md). `resources list` shows installed resources. `resources path NAME` prints a path for use in another command.
