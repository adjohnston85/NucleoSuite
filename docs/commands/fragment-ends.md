# `nucleosuite fragment-ends`

## What this command does

`fragment-ends` counts fragment boundaries at genomic positions. It can write genomic-left ends, genomic-right ends, or both ends combined.

## Why use it

Fragment-end tracks are useful for studying cleavage/breakpoint enrichment, comparing end positions between fragment-length classes, or measuring offsets between ends and dyads with DCC.

## How it works

For a fragment stored as BED interval `[start,end)`, the left end is the first covered base (`start`) and the right end is the last covered base (`end - 1`). Each accepted fragment contributes one count to the selected endpoint tracks.

See [Coverage, dyads, and fragment ends](../ALGORITHMS.md#coverage-dyads-and-fragment-ends).

## Basic usage

```bash
nucleosuite fragment-ends \
  --bam sample.bam \
  --frag-lower 145 \
  --frag-upper 147 \
  --tracks left right \
  --output-format bigwig \
  --out-prefix sample_145_147
```

## Coordinate pile-up cap

`--max-per-coordinate` limits the maximum accumulated value that can be assigned to any single genomic base in the fragment-end track. It is applied after endpoint contributions are accumulated and is separate from `--max-duplicates`, which filters identical complete fragments.

## Outputs

Depending on `--tracks`, the command writes combined, left-end, and/or right-end signal plus fragment summaries.

[Back to the command reference](../COMMAND_REFERENCE.md)
