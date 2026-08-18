# `nucleosuite compare-positions`

## What this command does

`compare-positions` compares the genomic positions in two BED-compatible callsets. It reports how close matched calls are, whether one caller is systematically shifted relative to the other, and how their scores relate.

## Why use it

Use it when matching calls by genomic coordinate is more informative than comparing overlap between their BED intervals.

## Typical use

For a BED file whose exact summit is stored in column 7:

```bash
nucleosuite compare-positions \
  --bed-a caller_A.bed \
  --bed-b caller_B.bed \
  --summit-column-a 7 \
  --summit-column-b 7 \
  --score-column-a 5 \
  --score-column-b 5 \
  --label-a A \
  --label-b B \
  --output-prefix A_vs_B
```

## Matching modes

### `unique` — default

One-to-one matching prevents the same target call from being reused. Use this when you want a concordance-style comparison where each call can participate in at most one accepted pair.

### many-to-one

Each query call is assigned its nearest target call on the same contig. A target can be reused for several query calls.

The suite PNS/WPS comparisons use many-to-one nearest matching in both directions so each callset is queried against the other independently.

## What the signed distance means

The output records the positional difference between the matched summits. Zero means the calls coincide. The sign shows which call lies upstream/downstream according to the command's documented A/B convention.

Absolute distance ignores direction and reports only the size of the separation.

## Compare scores as well as positions

When score columns are supplied, the command reports matched score relationships and can stratify directional comparisons by score rank.

Directional score-percentile groups target equal frequencies. Equal scores are ordered deterministically by genomic and input order and can be divided across groups.

## What it writes

Outputs include matched-pair tables, unmatched/summary counts, distance distributions, score-comparison summaries, and optional directional score-group plots/tables.

## Blacklist handling

`--blacklist-bed` removes complete overlapping intervals before matching.

## Plot customization

Comparison figures use the shared plotting interface described in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
