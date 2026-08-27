# `nucleosuite peak-states`

## What this command does

`peak-states` assigns peaks to chromatin states and reports their state composition and coverage-adjusted enrichment. It can repeat those measurements across peak-score thresholds or non-overlapping score groups.

## Why use it

Use it to identify preferential peak occurrence in promoters, enhancers, insulators, transcribed regions, heterochromatin, or other annotated states, including changes among stronger peaks.

## Basic usage

```bash
nucleosuite peak-states sample_nucleosome_regions.bed \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --position-column 7 \
  --output-prefix sample_peak_states
```

BED column 5 is used as the peak score by default. For NucleoSuite PNS BED8 calls, column 7 contains the midpoint of the retained PNS region.

## How peaks are assigned

Each peak contributes one selected genomic position. By default that is the BED interval midpoint; `--position-column` can select an exact coordinate stored in the record.

The peak is assigned to the state interval containing that point. A point outside all supplied state intervals is reported as unassigned. If state intervals overlap, `--overlap-policy first` uses the first match in input order; `error` stops instead.

## What the main measurements mean

### Peak composition

The stacked-bar plot shows the percentage of **assigned peaks** in each state. Every populated bar therefore sums to 100%.

An additional column reports each state's percentage of **all retained peaks**, which includes unassigned peaks in the denominator.

### Peak density

Peak density reports the number of peaks per megabase of each state, allowing states with different genomic spans to be compared.

### Coverage-adjusted enrichment

Enrichment compares a state's fraction of assigned peaks with its fraction of the summed annotated state coverage.

- enrichment = 1: peak representation matches state coverage;
- enrichment > 1: more peaks than expected from coverage;
- enrichment < 1: fewer peaks than expected from coverage.

Coverage adjustment separates peak abundance from enrichment relative to annotated state coverage.

## Cumulative score thresholds

A score percentile threshold shows how state composition changes as progressively higher-scoring peaks are retained.

```bash
nucleosuite peak-states sample_nucleosome_regions.bed \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --pct-range \
  --pct-lower 0 \
  --pct-upper 99 \
  --pct-step 1 \
  --output-prefix sample_peak_states
```

A 90th-percentile threshold retains peaks at or above the 90th-percentile score, including complete score ties.

Exact nonuniform thresholds can be supplied with:

```bash
--pct-values 10,20,50,90,99
```

## Non-overlapping score groups

`--pct-bins` or `--pct-bin-size` assigns each peak to one non-overlapping score group.

```bash
nucleosuite peak-states sample_nucleosome_regions.bed \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --pct-values 0,10,30,50,90,100 \
  --pct-bins \
  --bin-tie-mode split \
  --pct-bin-seed 1 \
  --output-prefix sample_peak_state_bins
```

The tie mode controls what happens when many peaks have the same score at a group boundary:

- **`split`** (default) targets the requested group sizes and can divide a score tie between adjacent rank groups;
- **`keep`** keeps identical scores together, so group sizes can differ from the requested percentages.

`--pct-bin-size 1` creates 100 regular one-percent groups; `--pct-bin-size 5` creates 20 five-percent groups.

## Stacked-bar appearance

The default categorical x-axis gives each threshold or bin equal width. `--plot-x-axis continuous` uses numeric percentile spacing.

`--plot-bar-gap` controls the white space between adjacent stacked bars. The default is `0.18`; use:

```bash
--plot-bar-gap 0
```

to make adjacent bars touch.

## Outputs

- state coverage table;
- peak counts, percentages, density, and coverage-adjusted enrichment;
- threshold/bin summaries with retained, assigned, and unassigned counts;
- metadata describing assignment and coverage accounting; and
- a stacked state-composition figure.

## Blacklist handling

`--blacklist-bed` excludes complete overlapping peaks and removes blacklisted bases from state coverage before coverage-adjusted quantities are calculated.

## Plot customization

The stacked figure also accepts the shared plotting options in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Ernst J, Kheradpour P, Mikkelsen TS, et al. (2011). Mapping and analysis of chromatin state dynamics in nine human cell types. *Nature* 473, 43–49. https://doi.org/10.1038/nature09906
