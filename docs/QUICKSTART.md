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

The most commonly reused outputs are the nucleosome-score BigWig and the nucleosome-region BED/bigBed. Select `--scoring-method bns` for Boxcar Nucleosome Score or `--scoring-method tns` for Triangular Nucleosome Score; both use the same peak caller as PNS. See [`pns`](commands/pns.md), [PNS in Algorithms](ALGORITHMS.md#probabilistic-nucleosome-scoring), [BNS in Algorithms](ALGORITHMS.md#boxcar-nucleosome-scoring), and [TNS in Algorithms](ALGORITHMS.md#triangular-nucleosome-scoring).


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
nucleosuite filter-coverage \
  sample_nucleosome_regions.bed \
  --bigwig sample_coverage.bw \
  --coverage-threshold 2 \
  --position-column 7
```

## 3. Measure nucleosome spacing

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

## 4. Measure spacing by chromatin state

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

## 5. Aggregate signal around CTCF sites

```bash
nucleosuite aggregate \
  --bigwig sample_methodpns_mode167_lower137_upper197_smooth0x2_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --out-prefix sample_ctcf
```

This writes the average signal pattern and an individual-region heatmap.

## 6. Calculate DAC from a dyad signal

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


## 7. Replot an existing result

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

## 8. Run the coordinated suites

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

## 9. Find and reuse bundled resources

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
