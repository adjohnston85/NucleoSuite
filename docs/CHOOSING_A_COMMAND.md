# Choosing a NucleoSuite command

Choose a command from your input and the result you need. Most analyses begin with paired-end fragments or a genomic signal or peak file.

## I have a paired-end BAM

Use [`fragments`](commands/fragments.md) if you first want an explicit fragment BED file. Use [`tracks`](commands/tracks.md) when you want several fragment-derived tracks from one pass through the data.

For nucleosome-oriented signal:

- [`pns`](commands/pns.md) — score positions with endpoint probability triangles (PNS), a balanced central boxcar (BNS), or a centred unit-mass triangle (TNS).
- [`wps`](commands/wps.md) — score whether fragments protect a mode-sized window or terminate inside it; the mode is estimated automatically unless fixed explicitly.
- [`dyads`](commands/dyads.md) — place signal at fragment centres.
- [`coverage`](commands/coverage.md) — count how many fragments cover each base.
- [`fragment-ends`](commands/fragment-ends.md) — count fragment starts and ends.

For a coordinated analysis with standard defaults, use [`cfdna-suite`](commands/cfdna-suite.md) or [`mnase-suite`](commands/mnase-suite.md).

For a ChIP-seq, CUT&RUN, or CUT&Tag treatment plus control, use [`chip-suite`](commands/chip-suite.md). It defaults to PNS nucleosome-peak discovery over the resolved mode ±30 bp and estimates fragment modes by bootstrap-stabilized random sampling. Each replicate’s method-specific score and matching positive-score track are generated together with broad 1–1,000 bp coverage in one `tracks` pass. The score is normalized by `posPNS`, `posBNS`, or `posTNS` before treatment replicates are averaged; coverage is independently scaled to a non-zero mean of 100 for peak measurements. Every treatment replicate > every control replicate (`all-controls`) is the default Stage 1 gate, with `mean` available as a less conservative alternative. Welch p-values and BH FDR remain optional exploratory filters. Clusters start at gate-passing p < 0.05 seeds; membership can include S+G peaks or significant-only S peaks, and `--cluster-max-non-member-gap` controls bridging. The selected PNS/BNS/TNS method is reused for cluster heatmaps, aggregate profiles, confidence bands and directional NRLs. Control and breakpoint peaks are not called. Replicate groups are independent and may differ in size unless `--bam-mode merged` is selected. Four supplied groups add cluster-only log-scale empirical-Bayes interaction tests plus descriptive Venn and base-occupancy overlap outputs. Supply an integer such as `--mode 167` when the analysis mode is already known.

Use [`chip-compare`](commands/chip-compare.md) when two conditions already have completed Stage 1 manifests. It compares overlap-connected cluster loci from saved scaled-coverage BigWigs, produces matched cluster-centred PNS aggregates, and does not read the BAMs again.

## I already have nucleosome or other peak calls

Use [`mean-scale`](commands/mean-scale.md) when BigWig signal or BED-family scores should be expressed relative to a reference mean, including non-zero BigWig signal, interval-score means, alternate region-score means, or a supplied reference value.

Use [`filter-peaks`](commands/filter-peaks.md) to create a reusable peak subset by score, score percentile, interval length, BigWig coverage, or combinations of these filters. Use [`peak-score-frequency`](commands/peak-score-frequency.md) to compare peak-score distributions and [`peak-states`](commands/peak-states.md) to measure how peaks are distributed across chromatin states.

Use [`pns-peak-fdr`](commands/pns-peak-fdr.md) when an observed PNS peak BED and one or more identically processed fragment-randomized peak BEDs are available. It preserves every observed BED field and appends empirical FDR. An optional `--fdr` cutoff adds a filtered BED without replacing the complete annotated output.

Use [`compare-positions`](commands/compare-positions.md) when one main nucleosome callset should be compared with one or more other callsets using one-to-one positional matching, main-score percentile groups, and optional within-percentile statistics.

Use [`distances`](commands/distances.md) to measure nearest-neighbour or higher-order spacing. Use [`flank-spacing`](commands/flank-spacing.md) when you have nucleosome calls plus a BED of reference sites divided into categories and want to compare the distance between the nearest upstream and downstream nucleosomes for each category. Column 4 is the default category column.

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
  --output-prefix sample_ctcf
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
