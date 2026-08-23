# `nucleosuite chrom-sizes`

## What this command does

`chrom-sizes` writes the contig names and lengths stored in a BAM or CRAM header as a simple two-column table.

## Why use it

Many genomic formats and analyses need reference-sequence lengths. Use this command when you want a reusable chromosome-size file in the same contig order and naming convention as your alignment data.

## Typical use

```bash
nucleosuite chrom-sizes \
  --bam sample.bam \
  --output sample.chrom.sizes
```

The output looks like:

```text
chr1    249250621
chr2    243199373
```

Restrict the file to selected contigs with:

```bash
nucleosuite chrom-sizes \
  --bam sample.bam \
  --contigs chr1-22,chrX \
  --output sample.selected.chrom.sizes
```

## Use an alignment header directly

Commands that accept `--chrom-sizes` can read a BAM or CRAM path directly:

```bash
nucleosuite dac \
  --bigwig sample_dyad.bw \
  --chrom-sizes sample.bam \
  --scope combined_chromosomes \
  --out-prefix sample_dac
```

## Automatic output naming

If `--output` is omitted, NucleoSuite writes `<alignment-basename>.chrom.sizes` in the current directory. An explicit `--output` overrides this default.

[Back to the command reference](../COMMAND_REFERENCE.md)
