# `nucleosuite fragments`

## What this command does

`fragments` converts paired-end alignments into complete fragment intervals, or filters/combines existing fragment interval files.

Each output BED interval represents one accepted fragment from its genomic start to its exclusive end.

## Why use it

Materialize fragments when you want to:

- reuse exactly the same filtered fragments in several analyses;
- inspect or archive fragment coordinates;
- combine fragment sets before downstream analysis; or
- use BED/BED.gz/bigBed instead of repeatedly reconstructing fragments from BAM.

## Choose one input family

From BAM:

```bash
--bam sample.bam [sample2.bam ...]
```

From existing fragment intervals:

```bash
--fragments sample.bed [sample2.bed.gz sample3.bb ...]
```

Fragment files need only:

```text
chrom    start    end
```

Coordinates are zero-based and half-open, so fragment length is `end - start`.

## Typical BAM extraction

```bash
nucleosuite fragments \
  --bam sample.bam \
  --frag-lower 120 \
  --frag-upper 180 \
  --max-duplicates 1 \
  --dedup-scope all_bams \
  --split-contigs \
  --output-format bed.gz \
  --output-prefix sample_120_180
```

## Combine existing fragment files

```bash
nucleosuite fragments \
  --fragments run1.bed.gz run2.bed.gz \
  --chrom-sizes sample.bam \
  --contigs chr1-22,chrX \
  --frag-lower 100 \
  --frag-upper 250 \
  --max-duplicates 1 \
  --output-format both \
  --output-prefix combined
```

## Duplicate handling

`--max-duplicates 1` retains at most one copy of each identical complete fragment coordinate. `0` disables the cap.

With multiple inputs:

- `--dedup-scope all_bams` applies the limit across the full input collection;
- `per_bam` applies it independently within each source input.

## Outputs

Depending on `--output-format`, the command writes BED, BED.gz, and/or bigBed fragment intervals plus summaries describing the retained fragment population.

Use the resulting fragment file directly with commands such as `nuc-score`, `wps`, `tracks`, `dyads`, `fragment-ends`, and `fragment-lengths`.

## Automatic output naming

If `--output-prefix` is omitted, NucleoSuite derives the output prefix from the primary BAM/CRAM or fragment-input basename in the current directory. An explicit `--output-prefix` overrides this default.

[Back to the command reference](../COMMAND_REFERENCE.md)
