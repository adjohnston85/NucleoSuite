# `nucleosuite gene-expression`

## What this command does

`gene-expression` relates gene expression to two different features of nucleosome organization:

- **peak spacing** — are nucleosomes farther apart or closer together in more highly expressed genes?
- **signal periodicity** — does the strength of a repeating PNS/WPS/other signal component change with expression?

The two analyses can be run separately or together.

## Why use it

Use this command when you want to compare chromatin organization with expression across genes, tissues, or cell-line expression profiles. Correlation describes association; it does not by itself establish a causal relationship.

## Use bundled genes and expression data

The bundled hg19 genes and HPA expression resources can be resolved directly:

```bash
GENES="$(nucleosuite resources path hg19-genes)"
TISSUE_EXPR="$(nucleosuite resources path hpa-tissue-expression)"
CELL_META="$(nucleosuite resources path hpa-cell-line-metadata)"
```

`--resource-set hg19-gm12878` supplies the standard bundled resources where supported. Explicit paths can replace individual components.

The default expression value is **nTPM**. `TPM` and `pTPM` remain selectable because they represent different expression measures.

## Peak-spacing analysis

Spacing analysis starts from a peak BED. For every gene, NucleoSuite measures adjacent peak distances and summarizes them with the median spacing. Gene body, upstream flank, and downstream flank are analysed separately in strand-aware orientation.

A gene needs at least two peaks in the analysed region to have a spacing estimate. The default flank is 10 kb.

Example:

```bash
nucleosuite gene-expression \
  --expression expression.tsv \
  --peaks sample=sample_PNS_nucleosome_regions.bed \
  --analysis spacing \
  --focus-profile NB-4 \
  --output-prefix sample_PNS_expression
```

NucleoSuite transforms expression with `log2(value + 1)` for this analysis and calculates the selected Pearson or Spearman correlation across eligible genes. A high-confidence summary is also reported for genes with at least 60 peaks in the selected region by default.

## FFT periodicity analysis

FFT analysis measures the association between expression and repeating signal components at particular periods.

For each gene, the default region is the first 10 kb in the direction of transcription. Minus-strand genes are reversed so every gene is analysed in the same transcriptional orientation.

The signal is centred/detrended, tapered, transformed to a periodogram, smoothed, and interpolated onto integer periods in the requested range. The default period range is 120–280 bp.

Example:

```bash
nucleosuite gene-expression \
  --expression expression.tsv \
  --signal sample=sample_PNS_smoothed.bw \
  --analysis fft \
  --expression-value-column nTPM \
  --output-prefix sample_PNS_expression
```

At each period, NucleoSuite correlates per-gene spectral intensity with expression. The default ranking score averages intensities at **193, 196, and 199 bp** and ranks expression profiles by their correlation with that combined periodic signal. Change those periods with `--fft-ranking-periods` when another period range is biologically relevant.

## Run both analyses together

```bash
nucleosuite gene-expression \
  --expression expression.tsv \
  --peaks sample=sample_PNS_nucleosome_regions.bed \
  --signal sample=sample_PNS_smoothed.bw \
  --analysis all \
  --focus-profile NB-4 \
  --output-prefix sample_PNS_expression
```

Multiple samples can be supplied by repeating `--peaks NAME=FILE` or `--signal NAME=FILE`. `--control-sample` allows rank changes to be reported relative to a nominated sample.

## Defaults

| Setting | Default |
|---|---:|
| Signal type | PNS |
| Expression value | nTPM |
| Gene flank | 10,000 bp |
| High-confidence spacing subset | at least 60 peaks |
| Spacing expression transform | `log2(value + 1)` |
| FFT region | first strand-aware 10,000 bp |
| FFT period range | 120–280 bp |
| FFT ranking periods | 193, 196, 199 bp |
| FFT expression transform | `log2(max(value, 0.04))` |
| Correlation | Pearson |
| Minimum matched genes | 30 |

## What it writes

Spacing analysis writes expression-correlation summaries and figures by default. A compressed plot-source table retains only the points used for the spacing scatter so it can be recreated with `nucleosuite plot`.

FFT analysis writes period-by-expression correlations, expression-profile rankings, optional rank changes, QC summaries, and figures by default.

Add `--write-detail-tables` to retain the large per-gene spacing table and, for FFT analysis, the per-gene spectral table. Metadata records the files, expression columns, filtering, and analysis settings used.

## Blacklist handling

A gene is excluded when its one-base TSS anchor overlaps the blacklist. Peak intervals overlapping the blacklist are removed. In retained FFT windows, blacklisted BigWig bases do not contribute to centring, detrending, filtering, or the Fourier transform.

## Plot customization

Figures use the shared plotting options in [Plot customization](../PLOTTING.md).

## Automatic output naming

If `--output-prefix` is omitted, NucleoSuite derives the prefix from the primary signal input, then the primary peak input, or finally the expression-table basename, and appends `_gene_expression`.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Aken BL, Achuthan P, Akanni W, et al. (2017). Ensembl 2017. *Nucleic Acids Research* 45(D1), D635–D642. https://doi.org/10.1093/nar/gkw1104
- Uhlén M, Fagerberg L, Hallström BM, et al. (2015). Tissue-based map of the human proteome. *Science* 347, 1260419. https://doi.org/10.1126/science.1260419
- GTEx Consortium. (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* 369, 1318–1330. https://doi.org/10.1126/science.aaz1776
- Jin H, Zhang C, Zwahlen M, et al. (2023). Systematic transcriptional analysis of human cell lines for gene expression landscape and tumor representation. *Nature Communications* 14, 5417. https://doi.org/10.1038/s41467-023-41132-w
