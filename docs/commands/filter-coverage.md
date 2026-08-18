# `nucleosuite filter-coverage`

## What this command does

`filter-coverage` removes BED peaks whose coverage at a selected peak position is below a user-defined threshold. It reads the coverage value from a BigWig track and writes the retained BED rows unchanged.

## Why use it

Use this command when you already have a nucleosome or other peak BED and want to require a minimum amount of supporting coverage at each peak. It is useful for removing calls that occur in very low-coverage regions without changing the original peak coordinates, names, scores, or other BED fields.

## How it works

For each BED record, NucleoSuite chooses one genomic position. By default it uses the midpoint between BED start and end. `--position-column` can instead specify a one-based BED column containing an explicit zero-based genomic position, such as column 7 in NucleoSuite PNS/BNS/TNS nucleosome BED8 output.

The BigWig coverage value is read at that single base. A peak is retained when

```math
coverage\ge threshold.
```

A threshold of 2 therefore retains coverage values of 2, 2.5, 3, and higher, while values below 2 are removed. Missing or non-finite BigWig values are treated as zero.

## Typical use

Use BED column 7 as the summit or representative position:

```bash
nucleosuite filter-coverage \
  sample_nucleosome_regions.bed \
  --bigwig sample_coverage.bw \
  --coverage-threshold 2 \
  --position-column 7
```

If `--position-column` is omitted, the interval midpoint is used:

```bash
nucleosuite filter-coverage \
  peaks.bed \
  --bigwig sample_coverage.bw \
  --coverage-threshold 2
```

## Output naming

The filtered BED filename is generated automatically from the input BED and threshold. For example:

```text
sample_nucleosome_regions.bed
```

with

```text
--coverage-threshold 2
```

writes

```text
sample_nucleosome_regions_coverage_ge2.bed
```

A decimal threshold is retained in the filename, so `--coverage-threshold 2.5` writes `sample_nucleosome_regions_coverage_ge2.5.bed`.

Use `--output` when you want to choose a different filename.

## What it writes

The filtered BED contains the original retained rows without rewriting their fields. BED comments, `track` lines, and `browser` lines are also preserved.

A summary TSV is written beside the filtered BED. It reports the input files, threshold, position source, total peaks, retained peaks, filtered peaks, retained percentage, and the number of missing BigWig values treated as zero.

## Relationship to PNS, BNS and TNS

For a new PNS, BNS or TNS run, the same filtering can be applied directly while peaks are being written:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --peak-coverage-threshold 2 \
  --out-prefix sample_pns
```

The direct option uses internally calculated fragment coverage at BED column 7 and filters only nucleosome peaks. Breakpoint peaks are unchanged.

[Back to the command reference](../COMMAND_REFERENCE.md)
