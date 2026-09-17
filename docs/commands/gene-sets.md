# `nucleosuite gene-sets`

## What this command does

`gene-sets` classifies genes using chromatin states at their transcript transcription start sites (TSSs) and across the gene body. Each final gene has one interval, starting at its selected promoter TSS when a promoter rule applies.

## Why use it

Use this command when you want reproducible gene groups such as active, weakly active, or repressed genes for downstream DAC, aggregation, expression, or regional analyses.

## Basic usage

NucleoSuite includes the hg19 gene BED, matching Ensembl GRCh37 release-87 transcript TSSs, GM12878 ChromHMM states and default gene-set rules. The command uses these bundled annotations directly:

```bash
nucleosuite gene-sets \
  --genes-bed "$(nucleosuite resources path hg19-genes)" \
  --states-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --config "$(nucleosuite resources path default-gene-sets)" \
  --leftover-set-name leftover_genes \
  --output-dir gm12878_gene_sets \
  --venn-sets active_genes weak_genes repressed_genes
```

### Transcript annotation

The bundled transcript TSS table is available with:

```bash
nucleosuite resources path hg19-ensembl87-transcript-tss
```

The complete matching Ensembl GTF is also bundled:

```bash
nucleosuite resources path hg19-ensembl87-gtf
```

For another compatible annotation, use `--transcript-gtf FILE` or `--transcript-tss-tsv FILE`.

## How the rules work

`include_rule` uses:

```text
&    AND
|    OR
( )  grouping
```

For example, an ordinary gene-overlap rule is:

```text
9_Txn_Transition | 10_Txn_Elongation
```

`required_tss_state` selects genes with **at least one transcript TSS** overlapping the specified state. NucleoSuite selects the most upstream qualifying TSS in transcriptional orientation and adjusts the gene start to that position. The gene-body inclusion and exclusion rules are evaluated over this adjusted interval.

`forbidden_tss_states` excludes a gene when **any of its transcript TSSs** overlaps a listed state. `forbidden_gene_states` excludes genes with a listed state in the interval being classified. `exclude_if_candidate` resolves competing candidate sets into mutually exclusive final categories.

The rules produce two set types:

- **candidate set**: genes satisfying the inclusion rule;
- **final set**: candidate genes remaining after the configured exclusions.

A leftover category contains genes that did not enter **any** candidate set. Genes that entered multiple candidates and were subsequently excluded are not reclassified as leftover.

For example, if a gene qualifies for both active and weak candidates and the weak rule excludes active candidates, it is removed from the weak final set. Its active final membership is decided by the active rule. It cannot become leftover because it qualified for at least one candidate category.

See [Gene-set assignment](../ALGORITHMS.md#gene-set-assignment) for the exact set definition.

## Bundled default categories

The bundled rules use the following candidate requirements:

- **Active:** At least one transcript TSS overlaps `1_Active_Promoter`; the adjusted gene overlaps `9_Txn_Transition` or `10_Txn_Elongation` and has no `12_Repressed` overlap.
- **Weak:** At least one transcript TSS overlaps `2_Weak_Promoter`; the adjusted gene overlaps `9_Txn_Transition`, `10_Txn_Elongation` or `11_Weak_Txn` and has no `12_Repressed` overlap.
- **Repressed:** The original gene overlaps `12_Repressed`, and none of its transcript TSSs overlaps `1_Active_Promoter` or `2_Weak_Promoter`.

Configured candidate exclusions then make the final Active, Weak, and Repressed outputs mutually exclusive. The optional strict leftover group contains genes that entered none of the candidate sets.

## Outputs

The selected options control which outputs are written:

- candidate and final gene BED6 files;
- final one-base selected-TSS BED6 files;
- an assignment table with original and adjusted gene coordinates, selected transcript IDs, transcript-TSS counts and final membership;
- a `<output-prefix>_selected_transcript_tss.tsv` table linking each selected transcript to its adjusted gene interval;
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

`--blacklist-bed` excludes genes whose original one-base gene TSS anchor overlaps the blacklist before classification.

## Plot customization

Summary figures use the shared plotting interface described in [Plot customization](../PLOTTING.md).

## Automatic output naming

Outputs are written to the current directory by default. If `--output-prefix` is omitted, NucleoSuite combines the gene-annotation and state-annotation basenames and appends `_gene_sets`. `--output-dir` and `--output-prefix` can override these defaults.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Aken BL, Achuthan P, Akanni W, et al. (2017). Ensembl 2017. *Nucleic Acids Research* 45(D1), D635–D642. https://doi.org/10.1093/nar/gkw1104
- Ernst J, Kheradpour P, Mikkelsen TS, et al. (2011). Mapping and analysis of chromatin state dynamics in nine human cell types. *Nature* 473, 43–49. https://doi.org/10.1038/nature09906
