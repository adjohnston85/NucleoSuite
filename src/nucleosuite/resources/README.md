# Bundled resources

The package includes an hg19/GM12878 resource collection for reproducible workflows:

- `hg19/hg19_ensembl_genes.bed`
- `hg19/wgEncodeBroadHmmGm12878HMM.bed`
- `hg19/CTCF_fimo_MA0139.1_INTERSECT_ENCFF923ZBP_GM12878_CTCF_hg19.bed`
- `hg19/hg19-blacklist.v2.bed.gz`
- `default_gene_sets.tsv`

Use `nucleosuite resources list`, `nucleosuite resources path NAME`, or `nucleosuite resources copy --output-dir DIR` to access these files.

These annotations are specific to hg19/GRCh37 and GM12878. Use annotations matched to the assembly and cell type being analysed.

## ChromHMM states

The GM12878 15-state ChromHMM segmentation is derived from Ernst et al. (2011).
The bundled BED uses distinct column-9 `itemRgb` colours for all 15 state labels,
including separate shades for paired enhancer, transcription, heterochromatin,
and repetitive/CNV states. `peak-states` displays these labels
in numeric order from `1_` through `15_`.

## Blacklist

The bundled ENCODE/Boyle Lab hg19 blacklist v2 follows Amemiya et al. (2019). The MNase and cfDNA suites enable it when selected primary-contig lengths match hg19/GRCh37. `--blacklist-bed FILE` selects another blacklist; `--no-blacklist` disables filtering.

## CTCF sites

The CTCF resource contains JASPAR MA0139.1 motif sites found with FIMO and intersected with GM12878 CTCF ENCODE peak set ENCFF923ZBP. It uses BED6 columns: chromosome, start, end, site identifier, score, and motif strand. CTCF-centred aggregation uses column 6 to orient all sites in the motif direction.

## Gene-set rules

The default gene-set configuration defines active, weak, and repressed candidate sets. Directed exclusions make the final categories mutually exclusive. `leftover_genes` contains genes that enter none of the candidate sets; genes removed because they enter competing candidates remain unassigned. Assignment outputs record both candidate and final membership.

## Gene annotation

The bundled gene file was generated from Ensembl release 87
`Homo_sapiens.GRCh37.87.gtf` (Aken et al., 2017) and contains six columns:
chromosome, start, end, Ensembl gene ID, gene name and strand. It is the default gene annotation for `gene-expression` and for expression analysis in both `mnase-suite` and `cfdna-suite`.

## Expression resources

- `hpa-tissue-expression`: compressed Human Protein Atlas tissue-consensus nTPM table used by TSS expression-quintile analyses. The HPA tissue resource is described by Uhlén et al. (2015); the consensus expression resource also incorporates GTEx tissue RNA-seq data (GTEx Consortium, 2020).
- `hpa-cell-line-metadata`: Human Protein Atlas cell-line disease metadata (Jin et al., 2023) with Cellosaurus identifiers/annotations (Bairoch, 2018), used in expression-profile ranking tables.

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
