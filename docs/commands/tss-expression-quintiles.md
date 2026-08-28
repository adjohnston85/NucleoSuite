# `nucleosuite tss-expression-quintiles`

## What this command does

`tss-expression-quintiles` compares a genomic signal around transcription start sites across five gene-expression groups.

## Why use it

Use it to test whether a PNS, WPS, coverage, or other BigWig profile changes with expression level.

## Basic usage

The `hg19-gm12878` resource set supplies the bundled hg19 genes and HPA tissue-consensus expression table:

```bash
nucleosuite tss-expression-quintiles \
  --signal sample_PNS.bw \
  --sample sample \
  --signal-label SNS \
  --resource-set hg19-gm12878 \
  --tissue bone_marrow \
  --window 2000 \
  --output-prefix sample_PNS_bone_marrow
```

You can also resolve the underlying files explicitly:

```bash
GENES="$(nucleosuite resources path hg19-genes)"
EXPR="$(nucleosuite resources path hpa-tissue-expression)"
```

Explicit file options can replace either bundled resource.

## How the quintiles are formed

Genes are matched to the selected expression profile and sorted from lowest to highest expression. Ensembl gene ID is used as a deterministic secondary key. The ordered genes are divided into five groups whose sizes differ by at most one gene:

- `Q1_lowest`
- `Q2_20_40_percent`
- `Q3_middle`
- `Q4_60_80_percent`
- `Q5_highest`

Expression ties can be split between adjacent quintiles because the grouping is based on rank and equal group size.

## How the signal is aligned

Each gene contributes a window around its TSS. Minus-strand genes are reversed so negative relative positions are upstream and positive positions are downstream for every gene.

The mean at each relative position uses the genes that have a valid signal value at that position. Blacklisted positions remain missing and do not contribute to the mean.

## Choose another tissue

Underscores in the selector are interpreted as spaces. For example:

```bash
--tissue skeletal_muscle
```

selects `skeletal muscle` from the long-format tissue table.

## Outputs

The command writes mean signal at every relative TSS position for all five quintiles, a gene-count and expression-range summary, a combined profile figure, and run metadata.

## Suite integration

Both `mnase-suite` and `cfdna-suite` can run this analysis automatically. Suite options choose the tissue, window size, custom expression resource, or whether to skip the step.

## Blacklist handling

TSS anchors overlapping the blacklist are skipped. Blacklisted bases inside retained TSS windows remain missing even if ordinary absent BigWig values are otherwise treated as zero.

## Plot customization

Figures use the shared plotting interface described in [Plot customization](../PLOTTING.md).

## Automatic output naming

If `--output-prefix` is omitted, NucleoSuite derives the prefix from the signal-track basename and appends `_tss_expression_quintiles`.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Aken BL, Achuthan P, Akanni W, et al. (2017). Ensembl 2017. *Nucleic Acids Research* 45(D1), D635–D642. https://doi.org/10.1093/nar/gkw1104
- GTEx Consortium. (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* 369, 1318–1330. https://doi.org/10.1126/science.aaz1776
- Uhlén M, Fagerberg L, Hallström BM, et al. (2015). Tissue-based map of the human proteome. *Science* 347, 1260419. https://doi.org/10.1126/science.1260419
