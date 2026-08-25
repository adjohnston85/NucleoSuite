# Quick start

This page shows common NucleoSuite analyses as complete workflows. Replace filenames and contigs with those appropriate for your data.

## 1. Check the installation

```bash
nucleosuite --version
nucleosuite --help
```

Validate BAM/reference compatibility before a long run:

```bash
nucleosuite validate-inputs --bam sample.bam --fasta genome.fa
```

## 2. Call PNS nucleosome regions from a BAM

Use `pns` when you want a nucleosome-oriented signal and peak calls derived from fragment-end geometry:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1-22,chrX \
  --cores 8 \
  --out-prefix sample_pns
```

The most commonly reused outputs are the nucleosome-score BigWig and the nucleosome-region BED/bigBed. The protected-DNA mode is estimated automatically from the unsmoothed accepted-fragment histogram and printed during execution; use `--mode 167` to fix it explicitly. Select `--scoring-method bns` for Boxcar Nucleosome Score or `--scoring-method tns` for Triangular Nucleosome Score; both use the same peak caller as PNS. See [`pns`](commands/pns.md), [PNS in Algorithms](ALGORITHMS.md#probabilistic-nucleosome-scoring), [BNS in Algorithms](ALGORITHMS.md#boxcar-nucleosome-scoring), and [TNS in Algorithms](ALGORITHMS.md#triangular-nucleosome-scoring).


### Filter nucleosome peaks by coverage

When a coverage threshold should be applied during a PNS, BNS or TNS run:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --peak-coverage-threshold 2 \
  --out-prefix sample_pns
```

For an existing peak BED and coverage BigWig:

```bash
nucleosuite filter-peaks \
  sample_nucleosome_regions.bed \
  --coverage-bigwig sample_coverage.bw \
  --min-coverage 2 \
  --coverage-position-column 7
```

## 3. Assign empirical FDR from randomized peaks

Run the observed and randomized sample with identical PNS settings, then compare their peak BEDs:

```bash
nucleosuite pns-peak-fdr \
  sample_nucleosome_regions.bed \
  sample_randomized_control_nucleosome_regions.bed
```

This always writes every observed peak with `empirical_fdr` appended. Add `--fdr 0.05` to also write a filtered BED. [`cfdna-suite`](commands/cfdna-suite.md) and [`mnase-suite`](commands/mnase-suite.md) can generate both full workflows and perform this annotation with `--with-randomized-control`.

## 4. Run a matched ChIP-seq, CUT&RUN, or CUT&Tag analysis

The default `chip-suite` run uses TNS over 120–500 bp fragments and estimates target and control fragment modes independently before choosing an equal-weight pooled analysis mode:

```bash
nucleosuite chip-suite \
  --treatment1-bam target.bam \
  --control1-bam control.bam \
  --outdir target_chip_suite \
  --sample-name target \
  --cores 8
```

Use `--scoring-method bns` or `--scoring-method pns` to change the score. Use `--mode 167` to bypass automatic sampling and use exactly 167 bp for both target and control.

Each replicate TNS track is divided by its own mean `posTNS` before treatment tracks are averaged. This prevents sequencing-depth differences from weighting the consensus peak-discovery track toward one replicate. The score-track stage also writes raw fragment coverage and a separate coverage BigWig scaled to a non-zero mean of 100 for each sample. Maximum scaled coverage within each TNS-defined interval becomes the replicate peak score because it measures local abundance without making the score strongly dependent on peak width. A candidate passes the default Stage 1 selection when every treatment score exceeds every control score; p-values and BH FDR remain annotations unless optional cutoffs are supplied. Clusters use gate-passing peaks with p < 0.05. Treatment and control groups are independent, may contain different replicate counts, and are never paired by input order. Use `--bam-mode merged` to pool a group. Supplying `--treatment2-bam` and `--control2-bam` adds a log2 empirical-Bayes four-group interaction test. Two independently completed Stage 1 runs can instead be compared with [`chip-compare`](commands/chip-compare.md).

## 5. Measure nucleosome spacing

Use `--position-column 7` to measure from representative positions stored in BED column 7:

```bash
nucleosuite distances sample_methodpns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed \
  --position-column 7 \
  --min-distance 120 \
  --max-distance 250 \
  --scope combined_chromosomes \
  --output-prefix sample_spacing
```

The main histogram shows how often each peak-to-peak spacing occurs.

## 6. Compare one main nucleosome callset with multiple callsets

Use one main BED as the common reference and repeat `--compare-bed` for each additional callset:

```bash
nucleosuite compare-positions \
  --main-bed PNS=sample_pns_nucleosomes.bed \
  --compare-bed iNPS=sample_inps_nucleosomes.bed \
  --compare-bed DANPOS=sample_danpos_nucleosomes.bed \
  --stats \
  --output-prefix sample_position_compare
```

Each comparison is matched once with one-to-one unique pairs. The smaller callset is used as the query, but percentile grouping always uses the matched **main BED score**. Quartiles are the default, and the grouped percentile boxplot places the comparison methods side-by-side within each quartile. See [`compare-positions`](commands/compare-positions.md).

## 7. Measure spacing by chromatin state

Pass the installed path of the bundled GM12878 state BED directly:

```bash
nucleosuite distances sample_methodpns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed \
  --position-column 7 \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --min-distance 120 \
  --max-distance 250 \
  --scope combined_chromosomes \
  --output-prefix sample_spacing_by_state
```

## 8. Compare flanking spacing across categorized reference sites

If a BED contains reference-site categories in column 4, compare the spacing between the nearest upstream and downstream nucleosome calls for every category:

```bash
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosome_regions.bed \
  --region-bed categorized_sites.bed \
  --category-col 4 \
  --output-prefix sample_flank_spacing
```

Density curves are used by default. Categories are ranked by the default 190/260 bp density ratio, with the lowest ratio ranked first. The displayed x-axis extends to 500 bp by default. See [`flank-spacing`](commands/flank-spacing.md).

## 9. Aggregate signal around CTCF sites

```bash
nucleosuite aggregate \
  --bigwig sample_methodpns_mode167_lower137_upper197_smooth0x2_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --output-prefix sample_ctcf
```

This writes the average signal pattern and an individual-region heatmap.

## 10. Calculate DAC from a dyad signal

Use DAC when you want to detect distances at which the same signal repeats:

```bash
nucleosuite dac \
  --bigwig sample_dyad.bw \
  --chrom-sizes sample.bam \
  --scope combined_chromosomes \
  --dmax 2000 \
  --out-prefix sample_dac
```

Repeated nucleosome spacing appears as recurring DAC peaks. See [Distance autocorrelation](ALGORITHMS.md#distance-autocorrelation).


## 11. Replot an existing result

Once an analysis has finished, `plot` can recreate a figure directly from its TSV without rerunning the genomic calculation:

```bash
nucleosuite plot sample_dac.tsv \
  --x-major-tick 100 \
  --x-minor-tick 10 \
  --nrl-inset on
```

For a heatmap, the displayed colour range can be saturated explicitly:

```bash
nucleosuite plot sample_heatmap_matrix.tsv.gz \
  --vmin -10 \
  --vmax 10 \
  --x-minor-tick 10
```

See [`plot`](commands/plot.md) for automatic file detection, major/minor grid controls, and Matplotlib pass-through options.

## 12. Run the coordinated suites

For a standard cfDNA analysis:

```bash
nucleosuite cfdna-suite \
  --bam sample.bam \
  --fasta hg19.fa \
  --resource-set hg19-gm12878 \
  --contigs chr1-22,chrX \
  --cores 8 \
  --outdir sample_cfdna_suite
```

For MNase-seq:

```bash
nucleosuite mnase-suite \
  --bam sample.bam \
  --fasta hg19.fa \
  --resource-set hg19-gm12878 \
  --contigs chr1-22,chrX \
  --cores 8 \
  --outdir sample_mnase_suite
```

The suites coordinate fragment selection, tracks, calls, spacing, sequence profiles, and optional expression analyses in one output layout. Review the suite defaults before a large run.

To run a complete observed workflow and a complete fragment-randomized workflow, then annotate the observed combined peak BEDs with empirical FDR:

```bash
nucleosuite mnase-suite \
  --bam sample.bam \
  --fasta hg19.fa \
  --resource-set hg19-gm12878 \
  --outdir sample_mnase_suite \
  --with-randomized-control \
  --fdr 0.05
```

## 13. Find and reuse bundled resources

```bash
nucleosuite resources list
```

Common named resources include:

```bash
nucleosuite resources path hg19-genes
nucleosuite resources path gm12878-hg19-states
nucleosuite resources path gm12878-hg19-ctcf
nucleosuite resources path hg19-blacklist-v2
nucleosuite resources path hpa-tissue-expression
nucleosuite resources path default-gene-sets
```

Store a path in a shell variable when you will use it several times:

```bash
STATES="$(nucleosuite resources path gm12878-hg19-states)"
GENES="$(nucleosuite resources path hg19-genes)"
```

## Where to go next

- [Choosing a command](CHOOSING_A_COMMAND.md)
- [Workflows](WORKFLOWS.md)
- [Algorithms](ALGORITHMS.md)
- [Command reference](COMMAND_REFERENCE.md)
