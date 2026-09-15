# `nucleosuite flank-filter`

## What this command does

`flank-filter` filters one BED callset using the positions in a second BED callset. A record from `--regions` is retained only when `--flanks` contains at least one feature at a strictly smaller genomic position and at least one feature at a strictly larger genomic position on the same chromosome.

The main NucleoSuite use is to retain nucleosome calls that are flanked on both sides by breakpoint peaks.

For NucleoSuite BED8 peak files, the default representative position is BED column 7 (`thickStart`), which stores the retained peak centre. If column 7 is absent or non-numeric, the interval midpoint is used instead. Explicit one-based position columns can be selected independently for the region and flank BEDs.

## Why use it

Use `flank-filter` when a primary peak should only be accepted if a second feature type occurs on both sides. For PNS output, this provides a direct way to restrict nucleosome-region calls to nucleosomes bounded by upstream and downstream breakpoint peaks.

The command preserves every column of each retained `--regions` row unchanged, so the filtered BED can be passed directly into later NucleoSuite analyses.

## Basic usage

```bash
nucleosuite flank-filter \
  --regions sample_nucleosome_regions.bed \
  --flanks sample_breakpoint_peaks.bed \
  --out sample_nucleosome_regions_flanked.bed
```

A region is retained when the nearest available flank position on the left is strictly less than the region position and the nearest available flank position on the right is strictly greater. Flanks on another chromosome do not count.

## Limiting flank distance

By default there is no maximum distance between the retained region and its two flanking features. To require both nearest flanks to lie within a specified distance, use `--max-flank-distance`:

```bash
nucleosuite flank-filter \
  --regions sample_nucleosome_regions.bed \
  --flanks sample_breakpoint_peaks.bed \
  --max-flank-distance 150 \
  --out sample_nucleosome_regions_flanked_150bp.bed
```

The distance is measured between representative positions, independently on the upstream and downstream sides. Both distances must satisfy the limit.

## Choosing the representative position

The default is designed for NucleoSuite peak BED8 output:

1. use BED column 7 when it is present and numeric;
2. otherwise use the BED interval midpoint.

To force another absolute genomic-position column, use one-based column numbers:

```bash
nucleosuite flank-filter \
  --regions regions.bed \
  --flanks flanks.bed \
  --region-position-column 7 \
  --flank-position-column 7 \
  --out filtered.bed
```

## Inputs and outputs

`--regions` and `--flanks` accept BED3+, BED.gz, or bigBed inputs. The output is BED text; use an output name ending in `.gz` for gzip-compressed BED text.

Only rows from `--regions` are written. Their original columns and values are preserved exactly. Header/comment lines are not copied to the filtered output.

[Back to the command reference](../COMMAND_REFERENCE.md)
