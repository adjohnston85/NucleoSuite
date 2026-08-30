# `nucleosuite dac`

## What this command does

`dac` calculates distance autocorrelation (DAC) from one or more BigWig signals to identify distances at which the signal recurs.

A dyad track with regularly spaced nucleosomes is a common input. If dyad signal repeats about every 185 bp, the DAC profile should show peaks near 185 bp and its multiples.

## Why use it

Use DAC to compare periodicity across regions, fragment classes, or samples. The resulting distance profile can be passed to [`nrl`](nrl.md).

## How it works

For each distance `d`, NucleoSuite compares the signal with itself `d` bases away. Positions that have signal at both locations contribute their product. Those products are added within each analysis region.

The default DAC value then divides the raw product sum by the number of valid position pairs available at that distance. This opportunity correction matters because large distances naturally have fewer possible pairs.

A simple five-dyad example and the exact equations are shown in [Distance autocorrelation](../ALGORITHMS.md#distance-autocorrelation).

## Choose the regions to analyse

Supply one or more BigWigs and exactly one region source:

- `--chrom-sizes` for whole-genome or chromosome windows;
- `--regions-bed` for labelled genomic intervals; or
- `--genes-bed` together with gene selection options for gene-centred DAC.

Pairs never cross a region boundary.

## Basic usage

```bash
nucleosuite dac \
  --bigwig sample_dyad.bw \
  --chrom-sizes sample.bam \
  --scope combined_chromosomes \
  --dmax 2000 \
  --out-prefix sample_dac
```

This calculates one pooled DAC profile over the selected chromosome windows.

## DAC by chromatin state

The bundled GM12878 state annotation can be used directly:

```bash
nucleosuite dac \
  --bigwig sample_dyad.bw \
  --regions-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --state-column 4 \
  --dmax 2000 \
  --out-prefix sample_dac_by_state
```

Here each ChromHMM interval acts as its own region, and the state label in column 4 groups the resulting within-region pair products.

## Multiple BigWigs

When several BigWigs are supplied, NucleoSuite autocorrelates each track with itself and adds their raw DAC and opportunity counts. [`dcc`](dcc.md) calculates products between two signal collections.

## Normalization choices

Use `--no-normalize-dac` to report the uncorrected raw product sum.

The output also contains percentage and signal-depth-scaled columns derived from the raw and opportunity-normalized values. See [Algorithms](../ALGORITHMS.md#distance-autocorrelation) for the formulas.

## Sparse, FFT, and automatic calculation

- `sparse` works directly with non-zero signal positions;
- `fft` calculates the same autocorrelation from dense arrays; and
- `auto` chooses between them from signal density.

These are computational routes to the same pair-product definition.

## Blacklist and missing signal

Missing or non-finite BigWig values are treated as zero signal but are still included when counting the genomic positions available for comparison. Blacklisted bases are excluded entirely from both the signal calculation and the number of available comparisons.

## Outputs

Each analysed state, category, gene, or pooled profile receives:

- a DAC TSV containing distances, raw product sums, opportunities, the selected DAC value, and derived columns;
- a summary; and
- a DAC plot.

Distance 1 remains in the TSV but is omitted from the standard plot so longer-range structure is easier to see.

## What to do next

Estimate a recurring period from the DAC peaks with:

```bash
nucleosuite nrl sample_dac.tsv \
  --output-prefix sample_nrl
```

## Plot customization

DAC plots show the unsmoothed `DAC Value` column without peak detection or peak labels by default. This column is opportunity-normalized unless `--no-normalize-dac` was used. To annotate DAC peaks, use `--plot-label-points peaks`. Peak calling then uses the same resolution-driven smoothing and refinement method as `nucleosuite nrl`; control it with `--peak-resolution` (default 160 bp). The unsmoothed input profile remains visible in grey behind the local-maximum and detection-smoothed curves. See [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
