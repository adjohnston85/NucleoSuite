# Choosing a NucleoSuite command

Choose a command from the input you have and the result you need. Most analyses begin with paired-end fragments, a genomic signal, or a peak file.

## I have a paired-end BAM or fragment file

Use [`fragments`](commands/fragments.md) to make an explicit fragment BED, or [`tracks`](commands/tracks.md) when several fragment-derived tracks and ranges should be generated in one pass.

For a single nucleosome-oriented signal, use [`pns`](commands/pns.md). It produces native PNS and `posPNS` BigWigs plus nucleosome and breakpoint calls, with automatic or explicit protected-DNA mode. Use [`wps`](commands/wps.md) for window protection scores, [`dyads`](commands/dyads.md) for fragment centres, [`coverage`](commands/coverage.md) for covered bases, and [`fragment-ends`](commands/fragment-ends.md) for fragment boundaries.

For complete cfDNA or MNase-seq analyses, use [`cfdna-suite`](commands/cfdna-suite.md) or [`mnase-suite`](commands/mnase-suite.md). For matched CUT&RUN or CUT&Tag treatment/control data, use [`cutn-suite`](commands/cutn-suite.md); it coordinates PNS discovery, broad coverage measurement, replicate-aware clustering, and optional condition comparison.

## I have a signal or peak file

Use [`call-peaks`](commands/call-peaks.md) for PNS/WPS peak calls, [`mean-scale`](commands/mean-scale.md) when a signal or interval score should be expressed relative to a reference mean, and [`filter-peaks`](commands/filter-peaks.md) to make a reusable peak subset. Use [`peak-score-frequency`](commands/peak-score-frequency.md) or [`peak-states`](commands/peak-states.md) to summarize peak scores and chromatin-state enrichment.

Use [`distances`](commands/distances.md) for adjacent or higher-order spacing and [`flank-spacing`](commands/flank-spacing.md) for spacing between nucleosomes flanking categorized reference sites. Use [`compare-positions`](commands/compare-positions.md) to compare one main callset with other positional callsets, and [`empirical-peak-fdr`](commands/empirical-peak-fdr.md) for observed versus fragment-randomized calls.

## I want periodicity, offsets, or regional profiles

Use [`dac`](commands/dac.md) when one signal repeats at characteristic distances, [`dcc`](commands/dcc.md) when one signal is offset from another, and [`nrl`](commands/nrl.md) to estimate a recurring period from a DAC/DCC profile. Use [`aggregate`](commands/aggregate.md) for signal around reference sites and [`region-extract`](commands/region-extract.md) for per-region signal vectors.

## I want sequence, chromatin-state, or gene analyses

Use [`dinuc-profile`](commands/dinuc-profile.md) for positional dinucleotide frequencies and [`ww-types`](commands/ww-types.md) for centred WW/SS classification. Use [`gene-sets`](commands/gene-sets.md) to classify genes from state overlaps, [`gene-expression`](commands/gene-expression.md) for expression relationships, and [`tss-expression-quintiles`](commands/tss-expression-quintiles.md) for TSS-centred signal across expression groups.

Bundled annotations can be addressed directly:

```bash
STATES="$(nucleosuite resources path gm12878-hg19-states)"
GENES="$(nucleosuite resources path hg19-genes)"
EXPR="$(nucleosuite resources path hpa-tissue-expression)"
```

## I want to reuse or inspect a run

Use [`combine`](commands/combine.md) for existing chromosome-wise outputs, [`plot`](commands/plot.md) to recreate or customize figures, [`cutn-compare`](commands/cutn-compare.md) to compare completed CUT&RUN/CUT&Tag Stage 1 manifests, and [`resources`](commands/resources.md) to list or locate bundled annotation files.
