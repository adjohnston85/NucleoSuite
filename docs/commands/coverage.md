# `nucleosuite coverage`

## What this command does

`coverage` counts how many accepted fragments span each genomic base.

## Why use it

Use coverage to inspect fragment depth in a genome browser, extract regional depth, assess data quality, or compare depth with PNS, WPS, or dyad signals.

## How it works

Every accepted fragment contributes 1 to each base it covers. Overlapping fragments add together, so the final value at a base is the number of retained fragments spanning that position. See [Coverage, dyads, and fragment ends](../ALGORITHMS.md#coverage-dyads-and-fragment-ends).

## Typical use

```bash
nucleosuite coverage \
  --bam sample.bam \
  --frag-lower 137 \
  --frag-upper 197 \
  --output-format bigwig \
  --out-prefix sample
```

The fragment-length range should match the biological population whose coverage you want to inspect.

## Inputs

Use indexed paired-end BAMs or materialized fragment BED/BED.gz/bigBed files. Fragment-input runs can use a chromosome-size table, BAM, or CRAM for complete contig lengths.

## What it writes

The primary output is a coverage BigWig or compressed WIG, accompanied by fragment summary and fragment-length count outputs.

The standalone command retains one copy of each identical complete fragment by default. `--max-duplicates 0` disables that coordinate deduplication.

[Back to the command reference](../COMMAND_REFERENCE.md)
