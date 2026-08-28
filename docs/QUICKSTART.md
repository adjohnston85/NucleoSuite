# Quick start

This page shows common NucleoSuite entry points. Replace filenames and contigs with those appropriate for your data.

## 1. Check the installation

```bash
nucleosuite --version
nucleosuite --help
nucleosuite validate-inputs --bam sample.bam --fasta genome.fa
```

## 2. Generate a nucleosome-oriented signal and calls

Use `pns` for the length-adaptive probabilistic nucleosome score and its shared nucleosome/breakpoint peak caller:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1-22,chrX \
  --cores 8 \
  --out-prefix sample
```

The command estimates the protected-DNA mode from the accepted-fragment histogram unless an integer `--mode` is supplied. The default scoring bounds follow the resolved mode by `--frag-mode-padding` (30 bp). The PNS kernel has a positive mass of 100 and a negative mass of -100 per complete fragment; the positive distribution is therefore represented in percent. BigWig scores are written at their native values and are not rescaled. The `posPNS` track is a non-negative reference track for the same waveform.

For an existing PNS or WPS BigWig, use [`call-peaks`](commands/call-peaks.md). To require coverage at a PNS summit during generation:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --peak-coverage-threshold 2 \
  --out-prefix sample
```

Use [`tracks`](commands/tracks.md) when PNS, WPS, coverage, dyads, ends, sequence profiles, or multiple fragment ranges should be produced together.

## 3. Assign empirical FDR from randomized peaks

Run the observed and randomized samples with matching fragment and peak settings, then compare their peak BEDs:

```bash
nucleosuite empirical-peak-fdr \
  sample_nucleosome_regions.bed \
  sample_randomized_control_nucleosome_regions.bed
```

This writes every observed peak with empirical p-value and FDR columns. Add `--fdr 0.05` to also write a filtered BED. The cfDNA and MNase suites can generate randomized controls and perform this annotation with `--with-randomized-control`.

## 4. Run a matched CUT&RUN or CUT&Tag analysis

`cutn-suite` uses PNS discovery with a mode-centred fragment range and broad coverage for replicate measurement:

```bash
nucleosuite cutn-suite \
  --treatment1-bam target.bam \
  --control1-bam control.bam \
  --outdir target_cutn_suite \
  --sample-name target \
  --cores 8
```

Treatment and control modes are estimated independently and pooled for a compatible analysis geometry unless `--mode` is supplied. Replicate PNS BigWigs and the `posPNS` reference remain at native score values; coverage is normalized separately for Stage 1 measurements. Replicate-aware clustering, optional two-condition comparison, and downstream aggregates are recorded in the output manifests. See [`cutn-suite`](commands/cutn-suite.md).

## 5. Measure spacing and genomic context

Use column 7, the representative position in a peak BED, for distances:

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --position-column 7 \
  --min-distance 120 \
  --max-distance 250 \
  --scope combined_chromosomes \
  --output-prefix sample_spacing
```

Aggregate a signal around reference sites:

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --output-prefix sample_ctcf
```

For state-aware summaries, the bundled state BED can be inserted the same way: `--state-bed "$(nucleosuite resources path gm12878-hg19-states)"`. Use [`flank-spacing`](commands/flank-spacing.md) when the question is how nucleosome spacing changes on either side of a reference feature.

Use [`dac`](commands/dac.md), [`dcc`](commands/dcc.md), and [`nrl`](commands/nrl.md) for signal periodicity and offsets; [`region-extract`](commands/region-extract.md), [`peak-states`](commands/peak-states.md), and [`tss-expression-quintiles`](commands/tss-expression-quintiles.md) for genomic context; and [`gene-expression`](commands/gene-expression.md) for gene-centred relationships.

## 6. Reuse existing output

Use [`combine`](commands/combine.md) to rebuild combined chromosome outputs and [`plot`](commands/plot.md) to recreate or customize figures from saved NucleoSuite tables without rerunning the underlying analysis.
