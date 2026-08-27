# `nucleosuite mnase-suite`

## What this command does

`mnase-suite` runs the coordinated MNase-seq NucleoSuite workflow. With multiple contigs, tracks are produced per contig, combined, and normalization/downstream analyses are then performed once on the combined data.

## Why use it

Use the suite when the MNase track, sequence, spacing, periodicity and regional analyses should share the same fragment filtering, resources, provenance, combination step and post-combine normalization. Use standalone commands when only one analysis is needed or when different filters are required between analyses.

## Typical run

```bash
nucleosuite mnase-suite \
  --bam "merged_chr*.bam" \
  --fasta hg19.fa \
  --resource-set hg19-gm12878 \
  --contigs chr1-22,chrX \
  --cores 8 \
  --outdir sample_mnase_suite
```

## Workflow

```mermaid
flowchart TB
    A[MNase BAM or fragments] --> B[mnase-suite]
    B --> C[SNS, dyads, ends, sequence]
    C --> D[Combine chromosomes]
    D --> E[Scale SNS, posSNS, coverage]
    D --> F[DAC from 146-148 dyads]
    F --> G[NRL]
    E --> H[Aggregate and regional analyses]
    D --> I[SNS peak distances]
```

## Main defaults

| Setting | MNase default |
|---|---|
| SNS | 120–180 bp fragments; mode 147 bp |
| Ranged dyads/ends, DAC and WW/SS | 146–148 bp |
| Exact dyad and fragment ends | 147 bp |
| Exact dinucleotide profiles | 145 and 147 bp |
| SNS nucleosome distances | order 1 to 500 bp; orders 1–7 to 1500 bp |
| Long DAC-derived NRL | 1–1500 bp, resolution 160; first called peak excluded from regression |
| Short periodicity | 1–144 bp, resolution 1 |
| Intermediate periodicity | 150–220 bp, resolution 8 |


## Post-combine normalization

After chromosome combination, the suite retains the raw combined outputs and creates normalized analysis inputs:

- coverage mean-scaled to 100;
- posSNS mean-scaled to 100;
- SNS scaled to 100 relative to the mean column-5 score of the raw combined SNS nucleosome calls;
- combined SNS nucleosome-region BED scores mean-scaled to 100;
- combined SNS breakpoint-peak BED scores mean-scaled to 100.

The mean-scaled nucleosome and breakpoint BEDs are used by downstream peak-based suite analyses. SNS aggregate analyses use the scaled SNS track. Regional extraction uses the mean-scaled peak BEDs together with scaled SNS and scaled coverage.

## DAC, NRL and spacing

DAC is calculated only from the ranged 146–148 bp dyad track. Its DAC outputs feed:

1. NRL at 1–1500 bp, resolution 160, with the first called peak shown in the profile but omitted from regression;
2. short periodicity at 1–144 bp, resolution 1;
3. nucleosome-scale periodicity at 150–220 bp, resolution 8.

The combined SNS nucleosome calls are analysed separately for adjacent spacing (order 1, 1–500 bp) and for orders 1–7 (1–1500 bp), with the latter performing combined NRL regression.

## Other downstream analyses

The suite performs SNS peak calls, ChromHMM-stratified SNS spacing, CTCF/TSS aggregation, TSS expression quintiles, region extraction, fragment-length profiles and heatmaps, optional SNS gene-expression analysis, SNS positive runs, SNS peak-score-frequency analyses, dinucleotide profiles, WW/SS classification and WW/SS type-specific dyads.

`peak-score-frequency` uses the mean-scaled nucleosome and breakpoint BEDs directly with `--score-scale 1`; no additional display scaling is applied.

## Observed plus randomized execution and peak FDR

`--randomize` retains randomized-only execution. `--with-randomized-control` instead runs the complete observed workflow and complete randomized workflow with identical settings in one invocation. After both combined peak callsets are available, the suite writes observed combined nucleosome and breakpoint BEDs with `empirical_fdr` appended as the final column. Add `--fdr 0.05` to also write filtered combined BEDs.

See [Output layout](../OUTPUT_LAYOUT.md), [Workflows](../WORKFLOWS.md), and the command-line help for the complete accepted option set.

[Back to the command reference](../COMMAND_REFERENCE.md)
