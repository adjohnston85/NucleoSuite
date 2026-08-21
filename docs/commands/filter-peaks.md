# `nucleosuite filter-peaks`

## What this command does

`filter-peaks` filters nucleosome or other peak interval files by score, score percentile, and BED interval length. It accepts BED, BED.gz, or bigBed input and writes the filtered intervals in the same format by default.

## Why use it

Use it when a downstream analysis should be restricted to a score-defined peak subset, a particular nucleosome-region length range, or both. It can also transform output scores before writing, including absolute-value conversion and multiplicative scaling.

## Basic use

Retain peaks with scores of at least 10:

```bash
nucleosuite filter-peaks peaks.bed --min-score 10
```

Retain peaks from the upper 10% of the score distribution:

```bash
nucleosuite filter-peaks peaks.bed --score-percentile 90
```

Retain regions between 120 and 180 bp long:

```bash
nucleosuite filter-peaks peaks.bed \
  --min-length 120 \
  --max-length 180
```

Filters can be combined across score and length. For example:

```bash
nucleosuite filter-peaks peaks.bed \
  --score-percentile 90 \
  --min-length 120 \
  --max-length 180
```

## Score filtering

Scores are read from BED column 5 by default. Use `--score-column` when the filtering score is stored in another column.

Absolute score bounds are inclusive:

```bash
nucleosuite filter-peaks peaks.bed \
  --min-score 5 \
  --max-score 25
```

`--score-percentile P` instead retains scores at or above percentile `P`. Absolute score bounds and `--score-percentile` are mutually exclusive.

When length bounds are also supplied, the length filter is applied first and the percentile threshold is calculated from the scores of the length-eligible peaks. This makes a percentile describe the peak population that is actually eligible for output.

## Negative scores

Negative scores are left unchanged and are compared as signed values by default.

Use:

```bash
--abs-score
```

to use `abs(score)` for filtering. When this option is enabled, the output score is also written as its absolute value.

## Score scaling

Use `--score-scale` to multiply retained output scores. For example:

```bash
nucleosuite filter-peaks peaks.bed \
  --min-score 0.1 \
  --score-scale 100
```

BED and BED.gz output retain floating-point scores after the requested transformation.

For bigBed output, BED column 5 must follow the UCSC score range. NucleoSuite therefore rounds the transformed column-5 score to an integer and clamps it to **0-1000** during bigBed conversion. Negative values become 0 unless `--abs-score` converts them before scaling.

## Region-length filtering

Region length is calculated as:

```text
end - start
```

`--min-length` and `--max-length` are inclusive. Either bound can be used independently.

## Output format

By default, the filtered output follows the input representation:

- BED input -> BED output;
- BED.gz input -> BED.gz output;
- bigBed input -> bigBed output.

Use `--output-format bed`, `--output-format bed.gz`, or `--output-format bigbed` to override this behaviour.

When producing bigBed from BED/BED.gz input, provide chromosome sizes with `--chrom-sizes`. For bigBed input, NucleoSuite attempts to inherit the chromosome sizes embedded in that file automatically.

Automatic output names include the central filtering parameters so changing a score or length threshold does not silently overwrite a previous filtered set.

## Summary output

A companion summary TSV records the requested filters, any percentile-derived score threshold, score transformation settings, and counts of valid, retained, length-filtered, score-filtered, and malformed records.

## Relationship to `distances`

`nucleosuite distances` provides its own score-percentile and target-peak selection controls and now also accepts `--min-length` and `--max-length`. Use `filter-peaks` when you want a reusable filtered interval file; use `distances` when the filtering is specific to a spacing analysis.

[Back to the command reference](../COMMAND_REFERENCE.md)
