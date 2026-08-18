# `nucleosuite merge-bams`

## What this command does

`merge-bams` combines BAM files while retaining the original alignment records and tags.

## Why use it

Most NucleoSuite fragment-based commands accept multiple BAMs directly. Use `merge-bams` when a downstream tool or archive requires one alignment file. Use [`fragments`](fragments.md) to combine fragment coordinates into BED intervals.

## Typical use

```bash
nucleosuite merge-bams \
  --bam part1.bam part2.bam part3.bam \
  --output merged.bam
```

The input BAMs should be compatible with the same reference assembly and contig lengths.

## What it writes

The primary output is the merged BAM. Follow ordinary BAM indexing requirements before using it with commands that require indexed random access.

[Back to the command reference](../COMMAND_REFERENCE.md)
