# `nucleosuite cutn-suite`

## What this command does

`cutn-suite` coordinates matched target/control CUT&RUN or CUT&Tag analyses. With condition 1 treatment and control inputs it runs Stage 1 discovery and measurement. Supplying treatment and control inputs for a second condition adds Stage 2 comparison, cluster overlap summaries, differential testing, and matched aggregates.

The suite uses PNS as its nucleosome-oriented discovery score. It estimates or accepts a protected-DNA mode, generates PNS and `posPNS` tracks together with broad coverage, calls treatment candidates, measures support in each replicate, and builds replicate-aware clusters. PNS BigWigs and peak scores remain native; coverage is normalized independently for measurement.

## Why use it

Use it when matched target/control replicates need shared discovery, control-aware gating, cluster formation, and optional two-condition comparison with reusable manifests.

## Basic invocation

```bash
nucleosuite cutn-suite \
  --treatment1-bam target.bam \
  --control1-bam control.bam \
  --condition1-name target_condition \
  --outdir target_cutn_suite \
  --cores 8
```

Multiple BAMs in each group are treated as independent biological replicates unless `--bam-mode merged` is selected. Treatment and control groups may contain different numbers of replicates.

## Stage 1: discovery and measurement

### 1. Resolve PNS geometry

With `--mode auto` (the default), treatment and control fragment histograms are estimated independently from seeded genomic-block samples. Their within-group probability distributions are pooled with equal group weight to select a compatible analysis mode. An integer `--mode` bypasses estimation. The discovery range defaults to the resolved mode ±30 bp; change this with `--frag-mode-padding` or override either bound with `--score-frag-lower` and `--score-frag-upper`.

Broad coverage uses 1–1,000 bp fragments by default and can be changed with `--coverage-frag-lower` and `--coverage-frag-upper`. The two ranges are generated in one `tracks` pass per replicate: the PNS range focuses positioning discovery, while broad coverage retains fragment abundance for treatment/control measurement.

### 2. Generate replicate tracks

Each replicate receives:

- a signed native `pns` BigWig;
- a native non-negative `posPNS` reference BigWig;
- optional smoothed PNS output and PNS nucleosome/breakpoint calls; and
- broad fragment coverage, normalized later for measurement.

Each complete fragment contributes positive PNS mass 100 and negative mass -100. The positive distribution is represented in percent. No score BigWig is divided by a reference mean or otherwise rescaled by the suite. The `posPNS` track is retained as a waveform reference, not as a score-normalization step.

### 3. Call a consensus treatment set

Treatment PNS tracks are averaged at native scale to form the discovery signal. Candidate nucleosome regions are called once from this condition-level signal so every treatment and control replicate is measured against the same intervals. The sign-inverted signal supplies breakpoint calls when requested.

### 4. Normalize coverage for Stage 1 measurement

Coverage is normalized independently for each replicate to a non-zero mean of 100:

```math
Cov_{100,i}(x)=100\frac{Cov_i(x)}{\mathrm{mean}(Cov_i(x)\mid Cov_i(x)>0)}.
```

For a candidate interval $R$, the default measurement is the mean normalized coverage over the complete interval:

```math
P_i(R)=\frac{1}{|R|}\sum_{x\in R}Cov_{100,i}(x).
```

Use `--stage1-coverage-statistic max` when the interval maximum is preferred. Treatment and control are never globally subtracted; each replicate is measured independently.

### 5. Gate peaks and form clusters

The default seed/member rules adapt to replicate count. When either group has fewer than three replicates, both seed and extension gates use every treatment replicate > every control replicate. With at least three replicates in both groups, seeds require a one-sided raw p-value below `--cluster-seed-p-value` (default 0.05) and mean treatment > mean control; extension members use the all-controls gate by default.

The main controls are:

- `--cluster-seed-mode pvalue|gated` selects the seed requirement;
- `--cluster-seed-gate-mode mean|all-controls` selects the seed gate;
- `--stage1-gate-mode mean|all-controls` selects the extension gate;
- `--cluster-member-mode seed-and-gated|significant-only` selects cluster membership;
- `--cluster-max-non-member-gap` controls bridging through non-members;
- `--max-cluster-gap` limits the distance between adjacent included members; and
- `--min-cluster-members` sets the minimum cluster size.

The complete peak statistics table retains candidate-level measurements, group summaries, gates, seed/member status, and p-values. Stage 1 manifests retain all paths and parameter choices for inspection, reruns, and Stage 2 comparison.

## Stage 2: compare two conditions

Supply both condition pairs:

```bash
nucleosuite cutn-suite \
  --treatment1-bam wt_target_R1.bam wt_target_R2.bam \
  --control1-bam wt_control_R1.bam wt_control_R2.bam \
  --treatment2-bam mutant_target_R1.bam mutant_target_R2.bam \
  --control2-bam mutant_control_R1.bam mutant_control_R2.bam \
  --outdir cutn_comparison
```

Stage 2 groups directly or transitively overlapping clusters into comparison loci. Shared loci use the actual genomic intersection for measurement by default; condition-specific loci use their complete interval. Saved broad coverage is measured in four independent groups: condition 1 treatment/control and condition 2 treatment/control.

The interaction model tests whether target-over-control enrichment changes between conditions. It reports raw and moderated p-values, BH-adjusted differential FDR, effect sizes, confidence intervals, direction, contributing clusters, measurement intervals, and measured bases. All loci remain in the complete table; `--differential-fdr` controls separate significant gain/loss BEDs.

Stage 2 also reports cluster-count and occupied-base overlap, and creates matched cluster-centred PNS aggregates using the native Stage 1 score tracks. It does not revisit the BAMs.

## Inspect and rerun

Inspect a completed run:

```bash
nucleosuite cutn-suite --inspect-run target_cutn_suite
```

The inspection report lists conditions, replicate inputs, PNS mode and fragment ranges, coverage ranges, gates, cluster controls, retained tracks, and available mode-estimation details.

Reuse retained per-replicate tracks for downstream changes:

```bash
nucleosuite cutn-suite \
  --rerun-from target_cutn_suite \
  --exclude-sample target_R2.bam \
  --cluster-member-mode significant-only
```

Reruns do not regenerate BAM-derived PNS or coverage tracks. They can change peak, measurement, clustering, Stage 2, and aggregate settings. Initial-track parameters—PNS geometry, fragment ranges, input grouping, contigs, blacklist, and duplicate handling—are inherited and locked so retained tracks cannot be mixed with incompatible settings. The source run is never overwritten; a numbered rerun directory is created.

## Outputs

The output directory includes setup and mode reports, per-replicate native PNS/`posPNS` and coverage tracks, combined tracks, treatment peak and replicate-statistics tables, cluster files, Stage 1 manifests, and optional Stage 2 comparison/aggregate outputs. The layout is described in [Output layout](../OUTPUT_LAYOUT.md).

Use [`cutn-compare`](cutn-compare.md) when two Stage 1 runs already exist and only the comparison needs to be repeated.

[Back to the command reference](../COMMAND_REFERENCE.md)
