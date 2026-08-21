# `nucleosuite distances`

## What this command does

`distances` measures adjacent and higher-order separations between called genomic positions.

## Why use it

Use it to compare spacing across score groups or chromatin states, or to estimate NRL from several neighbour orders.

## Which position is measured

The input can be BED, BED.gz, or bigBed. By default, NucleoSuite uses the interval midpoint. If your peak file stores an exact summit or dyad in another column, use `--position-column`.

PNS BED8 output stores the midpoint of the retained PNS region in column 7, so a typical PNS spacing command uses:

```bash
--position-column 7
```

## How it works

After positions are sorted within the relevant contig or state, order 1 measures each position to the next position. Order 2 measures to the position after one intervening call, and so on.

For example, if nucleosome calls occur at:

```text
1000, 1185, 1370, 1555
```

then the order-1 distances are:

```text
185, 185, 185
```

and order-2 distances are:

```text
370, 370
```

The exact definition is given in [Peak distances](../ALGORITHMS.md#peak-distances).

## Typical nearest-neighbour analysis

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --position-column 7 \
  --min-distance 120 \
  --max-distance 250 \
  --max-order 1 \
  --scope combined_chromosomes \
  --output-prefix sample_spacing
```

A histogram maximum near 185 bp means many adjacent calls in the supplied callset are approximately 185 bp apart. Peak density and caller behavior affect this distribution.

## Compare spacing by chromatin state

Use the bundled GM12878 state annotation directly:

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --position-column 7 \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --state-label-column 4 \
  --min-distance 120 \
  --max-distance 250 \
  --scope combined_chromosomes \
  --output-prefix sample_spacing_by_state
```

Standard state-stratified output assigns each peak to the state containing its selected position and groups pairs whose endpoint labels match.

`--state-overlay-plot` restarts adjacency at each state interval, requiring both peaks to lie within the same ChromHMM segment.

## Compare equal-sized score groups

`--pct-bins` and `--pct-bin-size` create non-overlapping score groups. `--bin-tie-mode` controls equal scores at a group boundary:

```text
--bin-tie-mode split   # default
--bin-tie-mode keep
```

### `split`

Use `split` when you want the requested group sizes to be as close as possible to the requested percentages. Peaks are randomly ordered within score ties using the reproducible `--pct-bin-seed`, then rank-sorted. A tie can therefore be divided between adjacent bins.

For example:

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --pct-values 0,10,30,50,90,100 \
  --pct-bins \
  --bin-tie-mode split \
  --pct-bin-seed 1 \
  --output-prefix sample_distance_bins
```

creates approximately 10%, 20%, 20%, 40%, and 10% groups.

### `keep`

Use `keep` when identical scores should never be divided. Bin boundaries are converted to score cutoffs and tied scores remain together, so the final group sizes can differ from the requested percentages.

### Regular-width bins

```bash
--pct-bin-size 1
```

creates 100 one-percent rank groups; `--pct-bin-size 5` creates 20 five-percent groups. The same tie mode applies.

## Cumulative score-percentile sweeps

Without binning, percentile analysis is cumulative. A 90th-percentile threshold retains peaks at or above the 90th-percentile score.

A regular sweep can be requested with:

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --pct-range \
  --pct-lower 0 \
  --pct-upper 90 \
  --pct-step 10 \
  --output-prefix sample_percentile_sweep
```

or exact thresholds with:

```bash
--pct-values 10,20,50,90,99
```

Sweep outputs include count and percentage overlays plus a table of the retained peak counts at each threshold.

## Zero-count distances

Distance distributions are written with explicit zero-count positions by default. For each order/group, the `*_distances.tsv` output contains every integer distance from `--min-distance` through `--max-distance`; distances that were not observed receive a count of zero. This makes the raw profile continuous and prevents plotting software from connecting separated non-zero observations directly across missing integer positions.

The complete positive-distance distribution is still retained internally for smoothing, summary statistics, and neighbour-order mode detection. Dense zero filling is restricted to the requested reporting range, so a rare multi-megabase neighbour distance cannot cause NucleoSuite to allocate or write millions of unnecessary zero rows when, for example, `--max-distance 1500` is requested.

Use the opt-out when a sparse observed-only table is specifically desired:

```bash
nucleosuite distances peaks.bed \
  --no-include-zero-distances
```

`nucleosuite plot` exposes the same `--include-zero-distances` / `--no-include-zero-distances` control for distance tables. This allows an older sparse distance TSV to be replotted with explicit zeros without modifying the TSV itself.

## Higher-order distances and NRL regression

Set `--max-order` above 1 when you want first-, second-, third-, or higher-neighbour distances in the same run. The default plotted/regression x-maximum is **1500 bp**; use `--max-distance` to change the visible range. The same bounds define the default dense zero-filled reporting support. `--max-distance` does **not** truncate distance calculation, smoothing, or the underlying per-order distribution used for mode calling.

Distances are always calculated **within contigs**. By default `--scope combined_chromosomes` pools the within-contig distance counts across all supplied contigs for the output distributions. The combined NRL regression is also the default:

```text
--regression-scope combined
```

Use `--regression-scope contig` for separate regressions per contig, or `--regression-scope both` to write both forms.

For each populated order, NucleoSuite first constructs the **complete positive-distance distribution**, regardless of the requested plot limits. The raw or smoothed mode is then determined from that complete distribution. Only after the mode is known are `--min-distance` and `--max-distance` applied to decide whether that order's mode is displayed and is eligible for the regression. If the true +8 mode is 1560 bp and `--max-distance 1500` is used, the +8 curve is simply clipped at 1500 bp; there is no marker at 1500 and +8 is omitted from the regression. Changing `--max-distance` therefore cannot change the underlying mode call.

In smoothed mode (`--nrl-mode smoothed`), the Savitzky-Golay window is controlled by `--count-smooth-window` and `--count-smooth-polyorder`, which default to 21 and 2. Smoothing is evaluated on the full distribution rather than on the plotted x-window, so the plot boundary cannot create a terminal smoothing peak. Use `--nrl-mode raw` to use the full raw-count mode instead.

The retained in-range mode distance for each order is fitted against neighbour order. The slope of that fit is reported as the NRL estimate. For a well-ordered array, order 1 might peak near 185 bp, order 2 near 370 bp, and order 3 near 555 bp.

When multiple neighbour orders are plotted together in smoothed mode, each smoothed order is a different colour and the corresponding raw distribution is retained in grey behind it. Peak markers represent the same full-distribution modes used for regression eligibility. Add `--label-peaks` to label displayed modes; `--peak-label-value x|y|both` controls whether the label shows the modal distance, count, or both.

See [Nucleosome repeat length](../ALGORITHMS.md#nucleosome-repeat-length) for the related peak-period fitting used by the standalone `nrl` command.

## What it writes

The requested options determine which outputs are written:

- raw and percentage distance tables, zero-filled by default across the requested `--min-distance` to `--max-distance` reporting range;
- per-state distance tables when `--state-bed` is used;
- score-threshold or score-bin outputs;
- percentile-sweep figures plus per-figure curve and retained-peak tables suitable for faithful replotting;
- NRL regression point tables, summaries, and plots for higher-order analysis; and
- optional state-relative percentage overlays.

## Filtered peak output

The command can also write filtered peak intervals corresponding to the score selection. Use bigBed output when you want an indexed peak set for downstream genomic tools.

## Blacklist handling

When `--blacklist-bed` is supplied, peaks overlapping the blacklist are excluded before distances are calculated.

## Plot customization

Distance figures use the shared plotting interface described in [Plot customization](../PLOTTING.md). Command-specific controls for distance ticks, state overlays, percentile overlays, and score bins remain available alongside those shared options.

[Back to the command reference](../COMMAND_REFERENCE.md)

## References

- Ernst J, Kheradpour P, Mikkelsen TS, et al. (2011). Mapping and analysis of chromatin state dynamics in nine human cell types. *Nature* 473, 43–49. https://doi.org/10.1038/nature09906
