# `nucleosuite region-extract`

## What this command does

`region-extract` keeps the **individual region-level data** around BED features. It can extract per-base values from one or more BigWigs and report nearby upstream/downstream peaks for every input region.

## Why use it

Use this command to retain each region's values for downstream statistics, modelling, or inspection. [`aggregate`](aggregate.md) produces a mean profile across the regions, with an individual-region heatmap available through `--write-detail-tables`.

## Basic usage

```bash
nucleosuite region-extract \
  --bed regions.bed \
  --score-bw sample_pns.bw \
  --nucleosome-peaks sample_nucleosomes.bed \
  --peak-flank-bp 2000 \
  --out-prefix sample_regions
```

Multiple named signal tracks can be supplied with repeated `--signal-track NAME=BIGWIG` options.

## Region centres and peak positions

The input region's midpoint is the reference for upstream/downstream distances. Supplied peak files use an absolute centre from BED column 7 when it is present and numeric, falling back to the peak interval midpoint otherwise. `--peak-center-column` selects another one-based column; `--peak-center-column 0` always uses the peak midpoint.

`--peak-flank-bp` sets how far on either side of the input region centre to search for peaks. It does not change the region's per-base signal vector.

## Use bundled CTCF regions

```bash
nucleosuite region-extract \
  --bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --score-bw sample_pns.bw \
  --nucleosome-peaks sample_nucleosomes.bed \
  --out-prefix sample_ctcf_regions
```

This exports each CTCF site's PNS vector and nearest flanking nucleosome calls.

## Outputs

Depending on the supplied tracks:

- one signal table per named BigWig;
- a flanking-peak table; and
- skipped-line records for invalid or excluded BED entries.

## Multicontig processing

```bash
nucleosuite region-extract \
  --bed regions.bed \
  --score-bw sample_pns.bw \
  --cores 4 \
  --out-prefix sample_regions
```

The combined analysis is regenerated from the complete BED so region indices and final tables remain consistent across the full input.

## Blacklist handling

`--blacklist-bed` skips overlapping anchor records, excludes overlapping peak intervals, and preserves blacklisted signal positions as missing values inside retained vectors.

[Back to the command reference](../COMMAND_REFERENCE.md)
