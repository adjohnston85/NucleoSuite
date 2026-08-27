# `nucleosuite combine`

## What this command does

`combine` reruns the combination stage of an existing NucleoSuite multicontig analysis without repeating the completed per-contig analyses.

## Why use it

Use this when:

- per-contig outputs completed but the final combine stage failed or was interrupted;
- you originally used `--skip-combine`;
- you want to rerun combination after restoring missing per-contig outputs; or
- a suite/run needs combined BigWig, bigBed, BED, or summary products regenerated.

## Basic usage

```bash
nucleosuite combine \
  --input-dir sample_multicontig \
  --cores 4
```

The input directory must contain the NucleoSuite multicontig manifest and the per-contig outputs required for the registered combination steps.

## How combination works

NucleoSuite combines the values used to calculate each result, then recalculates the final output. For example, it sums fragment-length counts before calculating percentages, sums DAC/DCC products and opportunities before normalization, and combines BED records before building a bigBed. The result represents the selected contigs as one analysis.

## BigWig combination

The default direct method streams per-contig BigWigs into the final combined BigWig in bounded chunks. Workflows can also use staged bedGraph combination when requested by the original run.

## Recovery and validation

Combination uses the run manifest/checkpoints to determine which outputs belong together and verifies output prerequisites before publication.

## Plot customization

When combination regenerates plots from sufficient statistics, the plotting settings recorded by the original workflow are used where applicable.

[Back to the command reference](../COMMAND_REFERENCE.md)
