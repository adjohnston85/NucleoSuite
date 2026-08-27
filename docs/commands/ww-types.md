# `nucleosuite ww-types`

## What this command does

`ww-types` classifies fragments into four WW/SS sequence-pattern classes using a centred 147 bp core.

## Why use it

Use it when you want to separate fragments by their rotationally patterned WW/SS sequence organization and then compare class abundance, fragment coordinates, or type-specific dyad positioning.

## How it works

NucleoSuite counts WW and SS dinucleotides at predefined minor-groove-associated and major-groove-associated positions in the centred 147 bp core. Because the two position sets contain different numbers of sites, major-groove counts are put on the same scale before the WW and SS enrichment comparisons are made.

The two enrichment results define four fragment classes. See [WW/SS fragment classes](../ALGORITHMS.md#dinucleotide-profiles-and-wwss-classes) for the exact position counts, scaling, and type definitions.

## Basic usage

```bash
nucleosuite ww-types \
  --bam sample.bam \
  --fasta genome.fa \
  --frag-lower 145 \
  --frag-upper 147 \
  --out-prefix sample_145_147
```

## Outputs

Selected options control whether the command writes a combined classified BED, type-specific BED files and dyad tracks, type counts and percentages, and summary figures. Multicontig percentages are calculated from the combined counts.

## Plot customization

WW/SS figures use the shared plotting interface described in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Wright GM, Cui F. (2019). The nucleosome position-encoding WW/SS sequence pattern is depleted in mammalian genes relative to other eukaryotes. *Nucleic Acids Research* 47, 7942–7954. https://doi.org/10.1093/nar/gkz544
