# `nucleosuite gene-sets`

## What this command does

`gene-sets` groups genes by overlapping chromatin states according to configured inclusion and exclusion rules.

## Why use it

Use this command when you want reproducible gene groups such as active, weakly active, or repressed genes for downstream DAC, aggregation, expression, or regional analyses.

## Typical use with bundled resources

NucleoSuite includes an hg19 gene BED, the GM12878 ChromHMM states, and default gene-set rules. They can be passed directly:

```bash
nucleosuite gene-sets \
  --genes-bed "$(nucleosuite resources path hg19-genes)" \
  --states-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --config "$(nucleosuite resources path default-gene-sets)" \
  --leftover-set-name leftover_genes \
  --output-dir gm12878_gene_sets \
  --venn-sets active_genes weak_genes repressed_genes
```

## How the rules work

`include_rule` uses:

```text
&    AND
|    OR
( )  grouping
```

For example:

```text
1_Active_Promoter & (9_Txn_Transition | 10_Txn_Elongation)
```

requires an active promoter plus at least one of the two transcription-associated states.

`exclude_if_candidate` makes final categories mutually exclusive by removing genes that also qualify for named competing candidate sets.

The rules produce two set types:

- **candidate set**: genes satisfying the inclusion rule;
- **final set**: candidate genes remaining after the configured exclusions.

A leftover category contains genes that did not enter **any** candidate set. Genes that entered multiple candidates and were subsequently excluded are not reclassified as leftover.

See [Gene-set assignment](../ALGORITHMS.md#gene-set-assignment) for the exact set definition.

## Bundled default categories

The bundled rules create active, weak, and repressed candidates. Their final categories are arranged so active/weak/repressed outputs are mutually exclusive, with an optional strict leftover group for genes in none of the candidate sets.

## What it writes

The selected options control which outputs are written:

- candidate and final gene BED6 files;
- final one-base TSS BED6 files;
- a complete gene assignment table showing candidate and final membership;
- overlap/shared-category files; and
- optional summary/Venn figures.

`<output-prefix>_final_states.bed` stores the final category in BED column 4 and is designed for pooled state-aware analyses.

For example, calculate DAC separately for the final gene categories:

```bash
nucleosuite dac \
  --bigwig sample_dyad.bw \
  --regions-bed gm12878_gene_sets_final_states.bed \
  --state-column 4 \
  --out-prefix sample_gene_category_dac
```

## Blacklist handling

If `--blacklist-bed` is supplied, genes whose one-base TSS anchor overlaps the blacklist are excluded before classification.

## Plot customization

Summary figures use the shared plotting interface described in [Plot customization](../PLOTTING.md).

## Automatic output naming

Outputs are written to the current directory by default. If `--output-prefix` is omitted, NucleoSuite combines the gene-annotation and state-annotation basenames and appends `_gene_sets`. `--output-dir` and `--output-prefix` can override these defaults.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Aken BL, Achuthan P, Akanni W, et al. (2017). Ensembl 2017. *Nucleic Acids Research* 45(D1), D635–D642. https://doi.org/10.1093/nar/gkw1104
- Ernst J, Kheradpour P, Mikkelsen TS, et al. (2011). Mapping and analysis of chromatin state dynamics in nine human cell types. *Nature* 473, 43–49. https://doi.org/10.1038/nature09906
