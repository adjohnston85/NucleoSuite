# `nucleosuite dyads`

## What this command does

`dyads` places one fragment-centre contribution into a sparse genomic signal track for every accepted fragment. The singular alias `nucleosuite dyad` is also accepted.

## Why use it

For nucleosome-protected fragments, the fragment centre can estimate the nucleosome dyad. Use dyad tracks to compare fragment-length classes, measure periodicity with DAC, or measure offsets with DCC.

## How even-length fragments are represented

Odd-length fragments have one central base. Even-length fragments have two central bases because their geometric centre lies between bases.

The default:

```text
--even-dyad split
```

places 0.5 on each central base. `left` or `right` places the full count 1 on the selected central base.

For example, `[100,267)` is 167 bp long and contributes 1 at position 183. `[100,268)` is 168 bp long and contributes 0.5 at 183 and 0.5 at 184 with the default split rule. See [Dyads](../ALGORITHMS.md#dyads) for the coordinate definitions.

## Basic usage

```bash
nucleosuite dyads \
  --bam sample.bam \
  --frag-lower 145 \
  --frag-upper 147 \
  --even-dyad split \
  --output-format bigwig \
  --out-prefix sample_145_147
```

## Duplicate fragments versus coordinate pile-ups

These options control different stages:

- `--max-duplicates` limits identical complete fragment intervals;
- `--max-per-coordinate` optionally caps the final accumulated dyad value at one genomic coordinate.

The default coordinate cap is 0, meaning unlimited.

## Outputs

The command writes a dyad BigWig or WIG plus fragment summaries and length counts.

## What to do next

Calculate repeating spacing:

```bash
nucleosuite dac \
  --bigwig sample_145_147_dyad.bw \
  --chrom-sizes sample.bam \
  --dmax 2000 \
  --out-prefix sample_145_147_dac
```

[Back to the command reference](../COMMAND_REFERENCE.md)
