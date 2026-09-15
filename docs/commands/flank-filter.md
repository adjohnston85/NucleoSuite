# `nucleosuite flank-filter`

## What this command does

`flank-filter` filters one BED callset using the positions in a second BED callset. A record from `--regions` is retained only when its surrounding events form a clean **flank → region → flank** arrangement on the same chromosome.

The main NucleoSuite use is to retain nucleosome calls that are directly bounded by breakpoint peaks. In that use, the accepted pattern is:

```text
breakpoint -- nucleosome -- breakpoint    KEEP
```

Two consecutive nucleosome-region calls between the same breakpoint peaks are not considered flanked:

```text
breakpoint -- nucleosome -- nucleosome -- breakpoint    DISCARD BOTH
```

Consecutive flank features are allowed. For example, the nucleosome below is retained because its immediately surrounding event types are still breakpoint and breakpoint:

```text
breakpoint -- breakpoint -- nucleosome -- breakpoint    KEEP
```

For NucleoSuite BED8 peak files, the default representative position is BED column 7 (`thickStart`), which stores the retained peak centre. If column 7 is absent or non-numeric, the interval midpoint is used instead. Explicit one-based position columns can be selected independently for the region and flank BEDs.

## Why use it

Use `flank-filter` when a primary peak is only valid when the nearest surrounding event type on both sides is the second feature type. For PNS output, this restricts nucleosome-region calls to isolated nucleosomes bounded directly by upstream and downstream breakpoint peaks.

The command preserves every column of each retained `--regions` row unchanged, so the filtered BED can be passed directly into later NucleoSuite analyses.

## Basic usage

```bash
nucleosuite flank-filter \
  --regions sample_nucleosome_regions.bed \
  --flanks sample_breakpoint_peaks.bed \
  --out sample_nucleosome_regions_flanked.bed
```

For each candidate region, NucleoSuite finds the nearest flank at a strictly smaller representative position and the nearest flank at a strictly larger representative position. The candidate is retained only when it is the **only `--regions` record between those two flanks**. This is equivalent to requiring the immediate event types around the candidate to be `FLANK -- REGION -- FLANK`.

Upstream and downstream flanks must have distinct representative positions on the same chromosome as the candidate.

## Limiting flank distance

By default there is no maximum distance between the retained region and its two directly bounding flank features. To require both flanks to lie within a specified distance, use `--max-flank-distance`:

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
