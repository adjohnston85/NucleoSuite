# Command reference

Each command page describes the calculation, a typical invocation, analysis options, and outputs. `nucleosuite COMMAND --help` shows the required inputs and core controls; `nucleosuite COMMAND --help-all` shows every command-specific option and current default. Commands that generate figures keep the shared plotting block under `nucleosuite COMMAND --help-plotting`.

## Prepare and inspect fragments

- [`fragments`](commands/fragments.md) — convert paired-end alignments to fragment intervals or combine fragment files.
- [`merge-bams`](commands/merge-bams.md) — merge BAMs while retaining alignment records and tags.
- [`randomize-fragments`](commands/randomize-fragments.md) — make coordinate-randomized control fragments.
- [`fragment-lengths`](commands/fragment-lengths.md) — count fragment lengths, optionally by genomic region class.
- [`fragment-heatmap`](commands/fragment-heatmap.md) — compare fragment-length profiles across samples or classes.

## Build genomic signals

- [`tracks`](commands/tracks.md) — make several fragment-derived tracks in one fragment pass.
- [`pns`](commands/pns.md) — calculate PNS, BNS or TNS and call positive/negative regions.
- [`wps`](commands/wps.md) — calculate WPS-family tracks and WPS peak calls.
- [`coverage`](commands/coverage.md) — count fragment coverage per base.
- [`mean-scale`](commands/mean-scale.md) — scale a BigWig relative to a supplied, region-score, or non-zero-signal reference mean.
- [`dyads`](commands/dyads.md) — place signal at fragment centres.
- [`fragment-ends`](commands/fragment-ends.md) — place signal at fragment starts and ends.
- [`dinuc-profile`](commands/dinuc-profile.md) — calculate dyad-aligned dinucleotide frequencies.
- [`ww-types`](commands/ww-types.md) — classify fragments from centred WW/SS patterns.

## Analyse peaks and spacing

- [`call-peaks`](commands/call-peaks.md) — call PNS or WPS features from an existing signal track.
- [`distances`](commands/distances.md) — measure adjacent and higher-order peak spacing.
- [`flank-spacing`](commands/flank-spacing.md) — compare category-wise spacing between nucleosomes flanking reference sites and rank the distributions.
- [`compare-positions`](commands/compare-positions.md) — compare one main callset with positional callsets and optional BigWig score comparators.
- [`peak-score-frequency`](commands/peak-score-frequency.md) — compare peak-score distributions.
- [`filter-coverage`](commands/filter-coverage.md) — retain peaks that meet a BigWig coverage threshold at their selected position.
- [`peak-states`](commands/peak-states.md) — count peaks and score-dependent enrichment across chromatin states.
- [`positive-runs`](commands/positive-runs.md) — measure continuous signal intervals above a threshold.

## Analyse periodicity and offsets

- [`dac`](commands/dac.md) — measure distances at which one signal repeats.
- [`dcc`](commands/dcc.md) — measure signed or absolute offsets between two signals.
- [`nrl`](commands/nrl.md) — estimate a recurring period from repeated DAC/DCC maxima.

## Analyse regions and genes

- [`aggregate`](commands/aggregate.md) — make a heatmap and average signal profile around reference sites.
- [`region-extract`](commands/region-extract.md) — retain per-region signal vectors and flanking peaks.
- [`gene-sets`](commands/gene-sets.md) — assign genes to chromatin-state-derived categories.
- [`gene-expression`](commands/gene-expression.md) — relate expression to spacing or signal periodicity.
- [`tss-expression-quintiles`](commands/tss-expression-quintiles.md) — compare TSS-centred signal across expression quintiles.

## Run coordinated workflows

- [`cfdna-suite`](commands/cfdna-suite.md) — coordinated cfDNA fragmentomics workflow.
- [`mnase-suite`](commands/mnase-suite.md) — coordinated MNase-seq workflow.

## Replot and customize figures

- [`plot`](commands/plot.md) — recreate and deeply customize all applicable figures from existing NucleoSuite output tables.

## Utilities

- [`resources`](commands/resources.md) — list bundled resources and print paths for direct reuse in other commands.
- [`chrom-sizes`](commands/chrom-sizes.md) — write contig names and lengths from an alignment header.
- [`validate-inputs`](commands/validate-inputs.md) — check input/reference compatibility before analysis.
- [`combine`](commands/combine.md) — rerun the combine stage of a multicontig analysis.
