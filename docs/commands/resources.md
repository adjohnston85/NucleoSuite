# `nucleosuite resources`

## What this command does

`resources` lists bundled reference files, prints their installed paths, validates them, or copies them into a project.

## Why use it

Bundled resource names let workflows locate compatible annotations without hard-coding installation-specific paths.

## See what is installed

```bash
nucleosuite resources list
```

To show only the resources in one configured resource set:

```bash
nucleosuite resources list --resource-set hg19-gm12878
```

The `hg19-gm12878` set groups resources that are intended to work together on hg19/GRCh37 data.

## Get a resource path

`resources path` prints the installed filesystem path of one named resource:

```bash
nucleosuite resources path gm12878-hg19-states
```

Common resource names include:

```bash
nucleosuite resources path hg19-genes
nucleosuite resources path gm12878-hg19-ctcf
nucleosuite resources path hg19-blacklist-v2
nucleosuite resources path default-gene-sets
nucleosuite resources path hpa-tissue-expression
nucleosuite resources path hpa-cell-line-metadata
```

## Use a resource directly in another command

Because `resources path` prints only the path, shell command substitution can place it directly into another option.

For example, analyse nucleosome spacing by the bundled GM12878 chromatin states:

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --position-column 7 \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --output-prefix sample_distances_by_state
```

Aggregate signal around the bundled CTCF sites:

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --output-prefix sample_ctcf
```

Use the bundled genes in a gene-centred analysis:

```bash
nucleosuite gene-sets \
  --genes-bed "$(nucleosuite resources path hg19-genes)" \
  --states-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --config "$(nucleosuite resources path default-gene-sets)" \
  --output-dir gm12878_gene_sets \
  --output-prefix gm12878_gene_sets
```

If you will reuse a resource several times, store its path once:

```bash
STATES="$(nucleosuite resources path gm12878-hg19-states)"
CTCF="$(nucleosuite resources path gm12878-hg19-ctcf)"
GENES="$(nucleosuite resources path hg19-genes)"
```

Then use `$STATES`, `$CTCF`, or `$GENES` in later commands.

## Logical names inside a resource set

`resources show` resolves a logical name such as `states`, `genes`, or `ctcf` within a resource set:

```bash
nucleosuite resources show states --resource-set hg19-gm12878
nucleosuite resources show genes --resource-set hg19-gm12878
nucleosuite resources show ctcf --resource-set hg19-gm12878
```

Logical names allow a workflow to select resources from a named set.

## Validate the bundled files

```bash
nucleosuite resources validate --resource-set hg19-gm12878
```

Validation checks that manifest entries resolve to installed files and verifies checksums where the manifest provides them.

## Copy resources into a project

Copy all bundled resources:

```bash
nucleosuite resources copy --output-dir resources_copy
```

Copy selected files:

```bash
nucleosuite resources copy \
  --output-dir resources_copy \
  --name hg19-genes \
  --name gm12878-hg19-states
```

Copy every resource referenced by a resource set:

```bash
nucleosuite resources copy \
  --output-dir resources_copy \
  --resource-set hg19-gm12878
```

## What is in the hg19/GM12878 set

The bundled collection includes:

- **hg19 genes** derived from Ensembl release 87 (`Homo_sapiens.GRCh37.87.gtf`; Aken et al., 2017);
- **GM12878 15-state ChromHMM annotations** from Ernst et al. (2011);
- **GM12878 CTCF sites** generated from FIMO scanning of JASPAR MA0139.1 and intersection with ENCODE CTCF CUT&RUN/CUT&Tag accession ENCFF923ZBP (Grant et al., 2011; Barski et al., 2007; Davis et al., 2018);
- **hg19 blacklist v2** from the ENCODE/Boyle Lab blacklist resource (Amemiya et al., 2019);
- **Human Protein Atlas tissue-consensus expression** and cell-line metadata, with GTEx and Cellosaurus provenance where relevant; and
- the default active/weak/repressed gene-category rules.

Use hg19/GM12878 resources only when the genome assembly and biological annotation context are appropriate for your analysis.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Aken BL, Achuthan P, Akanni W, et al. (2017). Ensembl 2017. *Nucleic Acids Research* 45(D1), D635–D642. https://doi.org/10.1093/nar/gkw1104
- Amemiya HM, Kundaje A, Boyle AP. (2019). The ENCODE Blacklist: Identification of Problematic Regions of the Genome. *Scientific Reports* 9, 9354. https://doi.org/10.1038/s41598-019-45839-z
- Bairoch A. (2018). The Cellosaurus, a cell line knowledge resource. *Journal of Biomolecular Techniques* 29, 25–38. https://doi.org/10.7171/jbt.18-2902-002
- Barski A, Cuddapah S, Cui K, et al. (2007). High-resolution profiling of histone methylations in the human genome. *Cell* 129, 823–837. https://doi.org/10.1016/j.cell.2007.05.009
- Davis CA, Hitz BC, Sloan CA, et al. (2018). The Encyclopedia of DNA elements (ENCODE): data portal update. *Nucleic Acids Research* 46(D1), D794–D801. https://doi.org/10.1093/nar/gkx1081
- Ernst J, Kheradpour P, Mikkelsen TS, et al. (2011). Mapping and analysis of chromatin state dynamics in nine human cell types. *Nature* 473, 43–49. https://doi.org/10.1038/nature09906
- Grant CE, Bailey TL, Noble WS. (2011). FIMO: scanning for occurrences of a given motif. *Bioinformatics* 27, 1017–1018. https://doi.org/10.1093/bioinformatics/btr064
- GTEx Consortium. (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* 369, 1318–1330. https://doi.org/10.1126/science.aaz1776
- Jin H, Zhang C, Zwahlen M, et al. (2023). Systematic transcriptional analysis of human cell lines for gene expression landscape and tumor representation. *Nature Communications* 14, 5417. https://doi.org/10.1038/s41467-023-41132-w
- Uhlén M, Fagerberg L, Hallström BM, et al. (2015). Tissue-based map of the human proteome. *Science* 347, 1260419. https://doi.org/10.1126/science.1260419
