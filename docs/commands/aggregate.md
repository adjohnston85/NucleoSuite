# `nucleosuite aggregate`

## What this command does

`aggregate` extracts equal-sized BigWig windows around reference sites, orients them by strand when requested, and writes an individual-region heatmap and mean profile.

## Why use it

Use `aggregate` to examine signal around BED-defined positions such as transcription start sites, transcription-factor sites, nucleosome calls, or chromatin-state features.

The heatmap shows whether the pattern is consistent across individual regions. The mean profile shows the average pattern across the complete retained set.

## How it works

Each accepted BED feature defines an aggregation centre. NucleoSuite extracts BigWig signal from `--window-half` bases on each side of that centre. Minus-strand regions are reversed so negative coordinates remain upstream and positive coordinates remain downstream.

The complete aggregate profile is the mean of the valid signal values at each relative position. The exact handling of missing and blacklisted positions is described in [Regional aggregation](../ALGORITHMS.md#regional-aggregation).

By default, `aggregate` also calls long-range peaks once across the complete negative-to-positive aggregate alignment. With the default 160 bp peak resolution, the complete profile is smoothed continuously across position 0 using 51 bp detection smoothing and 21 bp summit refinement. The resulting unified peak set is then divided by direction for two independent repeat-length regressions. The called peak nearest 0 within half the peak-calling resolution is assigned order 0 on both sides before regression filtering; remaining positive and negative peaks are ordered outward from it.

Every default x-axis is labelled `Distance from reference-site centre (bp)`. `--axis-label` remains available for an explicit alternative.

## Typical use around CTCF sites

Use the bundled GM12878 CTCF resource directly:

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --window-half 2500 \
  --output-dir sample_ctcf_aggregate \
  --output-prefix sample_ctcf
```

The CTCF BED stores motif strand in column 6, so the aggregated windows can be oriented consistently.

## Directional repeat length

The peak caller always operates across the complete aggregate alignment. It is not run separately on the two sides, and smoothing is not interrupted at position 0. After peak calling, positive and negative peaks are selected for separate regressions. Negative positions are converted to their absolute distance from 0 and ordered outward, so both fitted slopes are positive repeat lengths.

The called summit closest to 0 is treated as the shared order-0 peak when its absolute position is no greater than half `--nrl-peak-resolution` (±80 bp with the 160 bp default). This includes a central peak whose refined summit is slightly offset from exactly 0. If no called peak lies in that central interval, neither regression has an order-0 candidate.

Use `--nrl-regression-min` and `--nrl-regression-max` to control only the absolute-distance range entering both regressions. For example:

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed ctcf_sites.bed \
  --window-half 2500 \
  --nrl-regression-min 200 \
  --nrl-regression-max 1200 \
  --output-dir sample_ctcf_aggregate \
  --output-prefix sample_ctcf
```

This fits positive peaks from +200 to +1200 bp and negative peaks from −200 to −1200 bp. Peaks outside that range remain present in the unified peak table and plot.

By default, an inclusive regression-only exclusion interval spans half the peak resolution on either side of 0. It is therefore −80 to +80 bp at the default 160 bp resolution and changes automatically when `--nrl-peak-resolution` changes. This avoids making the reference-centred peak determine both outward fits. Use `--no-nrl-exclusion` to include an eligible central peak as order 0 in both regressions.

Explicit bounds replace the resolution-derived default:

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed ctcf_sites.bed \
  --nrl-regression-exclusion-start -250 \
  --nrl-regression-exclusion-end 250 \
  --output-dir sample_ctcf_aggregate \
  --output-prefix sample_ctcf
```

Both exclusion limits must be supplied. The shorter aliases `--nrl-exclusion-start` and `--nrl-exclusion-end` are also accepted. The explicit interval overrides the resolution-derived interval. Exclusion affects regression membership only: smoothing, unified peak calling and the complete peak plot remain unchanged. Excluded peaks retain their directional order numbers, so later peaks are not renumbered and the fitted repeat length is not compressed across the omitted interval. The unified profile plot shades the effective regression exclusion interval.

`--nrl-peak-resolution` controls the unified caller and defaults to 160 bp. Use `--no-nrl` to suppress aggregate peak calling and repeat-length outputs.

## Choosing the centre

By default, the interval midpoint is used. If a BED column contains the exact genomic coordinate you want to centre on, use `--point-col`.

`--nucleosome-bed` centres each reference region on a selected nearby nucleosome. `--nucleosome-offset` selects the strand-relative nucleosome; positive offsets are downstream and negative offsets are upstream.

## Heatmap rows versus the complete aggregate

The complete aggregate profile uses **all accepted regions**. `--max-heatmap-rows` limits only the plotted heatmap and its plotted-row mean.

`--subsample-mode` and `--seed` control how the plotted subset is selected.

## Missing signal and sparse tracks

By default, ordinary missing BigWig positions become zero, which supports sparse signals such as dyads. `--no-nan-to-zero` makes missing signal invalidate a window.

`--zero-thresh` can reject windows containing long zero runs. Set it to `0` when long zero stretches are expected and should not be treated as low-quality windows.

Recognized NucleoSuite BigWig suffixes set track-specific labels automatically. Unknown tracks use `Score` and `Mean score`. Dyad inputs ending in `_dyad.bw` also default to `--zero-thresh 0` and `--max-score inf`; explicit options always win. If every region is rejected, the error reports rejection counts and suggests these disabling values only when the corresponding filters rejected data.

## What it writes

`aggregate` writes:

- the heatmap matrix and row metadata;
- the complete aggregate profile from all accepted rows;
- the mean of the rows actually plotted in the heatmap;
- heatmap and mean-profile figures; and
- processing/parameter summaries.

With aggregate NRL enabled, it additionally writes:

- `_aggregate_nrl_profile.tsv` and `.png`, containing the complete unsmoothed profile, continuous 21 bp and 51 bp smoothed profiles, and every unified peak call;
- `_aggregate_nrl_peaks.tsv`, containing all called peaks, their signed positions, directional order numbers, shared-central status and regression inclusion or exclusion;
- `_aggregate_nrl_positive_regression.tsv` and `.png`, containing the positive-direction outward-distance fit;
- `_aggregate_nrl_negative_regression.tsv` and `.png`, containing the negative-direction outward-distance fit; and
- `_aggregate_nrl_summary.tsv`, containing both repeat lengths, fit statistics, caller settings and quality statuses.

The two regression plots are separate square figures with open circles and dotted fitted lines. All three figures can be recreated with [`nucleosuite plot`](plot.md).

The default half-window is 2500 bp, giving 5001 relative positions at base resolution.

Automatic output stems include the half-window, zero/maximum-score and missing-value filters, row sorting, NRL resolution, regression range, and effective exclusion interval. A supplied `--output-prefix` is a base prefix and receives the same tokens; plot-specific exact output options remain exact.

## Multicontig use

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed regions.bed \
  --contigs chr1 chr2 chr3 chr4 \
  --cores 4 \
  --output-dir aggregate \
  --output-prefix sample
```

The combined profile is calculated from the per-position sums and valid-row counts across all selected contigs.

## Blacklist handling

`--blacklist-bed` excludes reference anchors that overlap the blacklist. Blacklisted bases inside otherwise retained windows remain unavailable and do not contribute to the per-position average.

## Plot customization

Both the heatmap and mean-profile figures accept the shared plotting options described in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
