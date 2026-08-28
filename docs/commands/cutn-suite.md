# `nucleosuite cutn-suite`

## What this command does

`cutn-suite` performs control-aware nucleosome-score analysis for CUT&RUN or CUT&Tag. A run with condition 1 treatment and control BAMs performs Stage 1 only. Supplying all four condition groups performs Stage 1 independently for both conditions and then compares their enrichments in Stage 2.

If condition labels are not supplied, they are named `condition1` and `condition2`. The output prefix remains independently controlled by `--sample-name` and is not used as a biological condition label.

SNS is the default discovery score. Its fragment range is resolved separately for treatment and control as the selected mode ±30 bp. Use `--scoring-method pns`, `--scoring-method bns`, or `--scoring-method tns` when another kernel is required. The mode-centred score/positive-score pair and broad 1–1,000 bp coverage are generated together by `tracks` in one fragment pass per replicate.

## Why use it

Use this workflow to distinguish antibody-target enrichment from a condition-matched control and, optionally, test whether that enrichment changes between conditions. SNS locates peaks by default. Every coverage BigWig is separately divided by its own finite non-zero mean and multiplied by 100 before peak strength is measured. The control is not globally subtracted from the treatment track.

## Stage 1: one condition

```bash
nucleosuite cutn-suite \
  --treatment1-bam H3K4me3.bam \
  --control1-bam H3.bam \
  --condition1-name wild_type \
  --outdir wild_type_stage1 \
  --cores 8
```


Stage 1 separates **peak discovery** from **peak measurement**. SNS, PNS, BNS, or TNS defines where candidate peaks occur. Scaled fragment coverage then provides the replicate values used for treatment-versus-control filtering and statistics. Keeping these roles separate avoids treating the height of a model-derived positioning score as if it were direct fragment abundance.

### 1. Generate one score and coverage set per replicate

For each treatment and control replicate, `cutn-suite` makes two deliberately different fragment selections:

- an SNS, PNS, BNS, or TNS discovery track and its matching non-negative `posSNS`, `posPNS`, `posBNS`, or `posTNS` track, using fragments from the resolved mode ±30 bp by default; and
- raw fragment coverage using all accepted fragments from 1–1,000 bp by default.

The narrow, mode-centred range focuses the discovery score on fragments most consistent with one protected nucleosome and prevents very short or long assay fragments from changing the positioning signal. Change the automatic ±30 bp distance with `--frag-mode-padding`. `--score-frag-lower` and `--score-frag-upper` can override the lower and upper bounds independently. The broader coverage range retains the fragment abundance associated with the enriched locus, including subnucleosomal and longer fragments that can be informative in CUT&RUN or CUT&Tag. Change it with `--coverage-frag-lower` and `--coverage-frag-upper`.

The centred score locates protected-DNA structure. The positive track measures the overall amount of method-specific score support and is used only as the normalization reference. Raw broad-range coverage is retained so the original sequencing-depth scale remains available.

The selected scoring method is used throughout the workflow. SNS uses `sns`/`posSNS`, PNS uses `pns`/`posPNS`, BNS uses `bns`/`posBNS`, and TNS uses `tns`/`posTNS`. The same method-specific normalized treatment tracks used for discovery are reused for cluster-centred profiles, heatmaps and directional NRLs; selecting SNS, BNS or TNS does not trigger an additional PNS pass. Score tracks from the mode-centred range and broad 1–1,000 bp coverage are generated together by `tracks` in one fragment pass. `cutn-suite` then calls treatment nucleosome candidates once from the consensus discovery track and does not produce per-replicate nucleosome or breakpoint callsets.

### 2. Normalize and average treatment score tracks for discovery

Each replicate score track is divided by the finite, non-zero mean of its own positive-score track:

```math
Z_{i,scaled}(x)=\frac{Z_i(x)}{\mathrm{mean}(Z_i^+(x)\mid Z_i^+(x)>0)}.
```

This normalization is performed **before averaging** because raw score magnitude depends on usable fragment depth. Without it, a deeper replicate would contribute more strongly to the average and could dominate which peaks are discovered. The normalization places the replicate score tracks on comparable scales while preserving their spatial peak structure.

When several treatment replicates are supplied, their normalized score tracks are averaged. Averaging reinforces signal shared across replicates and reduces the influence of replicate-specific fluctuations. Peaks are then called once on this consensus treatment track, giving every downstream replicate exactly the same candidate intervals. With one treatment replicate, its normalized score track is used directly. Control score tracks are retained, but control peaks are not called.

### 3. Normalize coverage for measurement

Every treatment and control coverage BigWig is independently divided by its finite, non-zero mean and multiplied by 100:

```math
Cov_{100,i}(x)=100\frac{Cov_i(x)}{\mathrm{mean}(Cov_i(x)\mid Cov_i(x)>0)}.
```

A value of 100 therefore represents the mean among covered bases in that replicate. This scaling compensates for overall coverage differences so local values can be compared between replicates. It is applied independently and does not subtract control signal from treatment signal.

### 4. Measure each treatment candidate in every replicate

For each treatment-defined candidate interval, `cutn-suite` measures the **mean scaled coverage across the complete interval** in every treatment and control replicate. Mean is the default because it represents support across the called nucleosome interval rather than allowing one unusually high base to determine the replicate value.

```math
P_i(R)=\frac{1}{|R|}\sum_{x\in R}Cov_{100,i}(x).
```

Use `--stage1-coverage-statistic max` when the interval maximum is preferred. Column 5 of the complete treatment-peak BED is the corresponding statistic from the condition-mean treatment coverage track.

### 5. Choose S seeds and G cluster members

`cutn-suite` treats the rule that starts a cluster (**S**) separately from the rule that allows neighbouring treatment peaks to extend it (**G**).

The automatic defaults depend on biological replicate count.

When either treatment or control has fewer than three replicates:

```text
S = all treatment replicates > all control replicates
G = all treatment replicates > all control replicates
```

Peak p-values are not used because replicate counts are insufficient for the default statistical seed rule.

When treatment and control each have at least three replicates:

```text
S = raw p < 0.05 AND mean treatment > mean control
G = all treatment replicates > all control replicates
```

The one-sided raw p-value uses Welch's test of treatment mean > control mean. The default threshold can be changed with `--cluster-seed-p-value`.

The automatic rule is printed when the run starts. Explicit controls are available when another design is required:

- `--cluster-seed-mode pvalue|gated` chooses whether S requires the raw p-value or only its gate;
- `--cluster-seed-gate-mode mean|all-controls` changes the S gate independently;
- `--stage1-gate-mode mean|all-controls` changes the G gate independently.

The complete `target_peaks_replicate_statistics.bed` retains every treatment candidate and appends its raw p-value when one is calculated. `target_peak_replicate_statistics.tsv` reports every replicate interval measurement, group means, mean treatment-minus-control difference, `min(T)-max(C)` excess, both gate results, S/G status, and the raw p-value.

### 6. Form Stage 1 clusters

With the default `--cluster-member-mode seed-and-gated`, both S and G peaks are cluster members. An S peak remains an included member even when it passes the default mean seed gate but not the stricter all-controls G gate. `--cluster-member-mode significant-only` restricts membership to S peaks.

One consecutive non-member can bridge included members by default (`--cluster-max-non-member-gap 1`). A longer run ends the current cluster. Adjacent included-member summits can be at most 1,000 bp apart by default (`--max-cluster-gap 1000`), and a cluster requires at least two included members (`--min-cluster-members 2`). Bridging candidates do not contribute to the cluster boundary or score.

Cluster coordinates extend from the first included member to the last. The aggregate anchor is the discovery summit of the included member with the strongest condition-mean Stage 1 coverage measurement.

### 7. Aggregate SNS around cluster anchors

Each replicate SNS track is normalized by the mean of its `posSNS` track before averaging. This places replicate score tracks on comparable scales so that differences in sequencing depth do not cause higher-depth libraries to contribute disproportionately to the condition-average discovery track.

Cluster-centred heatmaps and aggregate profiles use the selected discovery scorer. The default directional NRL peak resolution is **130 bp** and the default fitted orders are 0–3 on both sides.

## Replicates and merged input

Multiple BAMs in each group are independent biological replicates by default. Treatment and control groups may contain different numbers of BAMs. They are not paired by command-line order:

```bash
nucleosuite cutn-suite \
  --treatment1-bam wt_mark_r1.bam wt_mark_r2.bam wt_mark_r3.bam \
  --control1-bam wt_H3_r1.bam wt_H3_r2.bam wt_H3_r3.bam \
  --outdir wt_stage1
```

Each replicate SNS track is normalized by the mean of its `posSNS` track before the discovery tracks are averaged. This places replicate discovery tracks on comparable scales so that higher-depth libraries do not contribute disproportionately to the condition mean. Replicate-specific scaled coverage supplies the Stage 1 treatment/control measurements, and condition-mean treatment coverage supplies the reported BED score. All raw and normalized tracks are retained in `cutn_stage1_manifest.json` so downstream comparison can reuse them without returning to the BAM files.

Use `--bam-mode merged` to pass every treatment BAM as one logical treatment sample and every control BAM as one logical control sample. This matches the usual NucleoSuite multi-BAM pooling behaviour. Merged mode provides Stage 2 effect sizes and gain/loss direction from the pooled treatment/control groups.

## Stage 1 plus Stage 2

Supply treatment and control BAMs for both biological conditions to run Stage 1 independently and then compare their clusters:

```bash
nucleosuite cutn-suite \
  --treatment1-bam wt_target_R1.bam wt_target_R2.bam wt_target_R3.bam \
  --control1-bam wt_control_R1.bam wt_control_R2.bam wt_control_R3.bam \
  --treatment2-bam mutant_target_R1.bam mutant_target_R2.bam mutant_target_R3.bam \
  --control2-bam mutant_control_R1.bam mutant_control_R2.bam mutant_control_R3.bam \
  --outdir cutn_comparison
```

Stage 2 operates on Stage 1 clusters and retained **raw coverage** tracks. Overlap-connected clusters are grouped into comparison loci. When clusters from both conditions overlap, replicate measurements are taken from the **actual overlapping genomic section by default**. Condition-specific clusters are measured across their complete locus.

The Stage 2 statistic is mean raw coverage over the comparison interval. The four independent groups are condition 1 treatment, condition 1 control, condition 2 treatment, and condition 2 control. The values are transformed with `log2(mean raw coverage + 1)` and fitted with a condition-by-treatment interaction model.

The comparison table reports raw interaction p-values, empirical-Bayes moderated p-values, BH-adjusted differential FDR, effect sizes, confidence intervals, direction, contributing cluster IDs, the exact measurement interval(s), and the number of measured bases. Significant gain/loss BEDs use `--differential-fdr`; all gain/loss outputs remain available for ranking.

Stage 2 also reports cluster-count and occupied-base overlap summaries and creates matched cluster-centred aggregate plots using a shared anchor set. It uses the retained Stage 1 files and does not return to the BAMs.

## Inspect a completed run

Completed `cutn-suite` directories can be inspected without opening the manifests manually:

```bash
nucleosuite cutn-suite --inspect-run cutn_results_h3K4me3
```

The report lists each biological condition, its treatment and control BAMs, the retained replicate tracks, scoring method, resolved treatment/control modes, discovery and coverage fragment ranges, Stage 1 gate and cluster settings. For each retained replicate, `cutn-suite` also reports lightweight sample statistics from the saved track outputs, including the fragment-length mode within the nucleosome mode-search range, the number of fragments used in the broad coverage range, the positive-score normalization mean, and the pre-scaling non-zero coverage mean when those values are available. The original treatment/control/pooled bootstrap mode estimate and confidence interval are also shown when the mode report is present.


## Fast reruns from retained BigWigs

`--rerun-from` reuses the per-replicate normalized score and broad-coverage BigWigs from a completed run. It does **not** repeat mode estimation or the BAM-to-BigWig `tracks` stage. This is useful for leave-one-replicate-out checks and for changing downstream peak, statistics, clustering, Stage 2, or aggregate parameters.

Exclude one replicate by BAM path, filename, or an unambiguous filename stem:

```bash
nucleosuite cutn-suite \
  --rerun-from cutn_results_h3K4me3 \
  --exclude-sample wt_K4_R2_sort.bam
```

Repeat `--exclude-sample` to exclude several replicates. When a basename or stem matches more than one retained BAM, the command stops and asks for the full BAM path rather than excluding multiple samples silently. A BAM cannot be removed from a source run created with `--bam-mode merged`, because the retained BigWig already contains the merged group; sample exclusion therefore requires per-replicate source tracks.

The rerun is written as a subdirectory of the source run. One exclusion produces names such as:

```text
rerun_excluding_wt_K4_R2_sort_01/
rerun_excluding_wt_K4_R2_sort_02/
```

Several exclusions use a shorter form such as `rerun_excluding_2_samples_01/`. A rerun with no exclusions is named `rerun_01/`. Existing reruns are detected and the numeric suffix is incremented automatically. The source run is never overwritten.

The retained replicate BigWigs are filtered first, then treatment/control condition means are recalculated as required. Peak discovery, replicate statistics, seed peaks, clusters, Stage 2 and cluster aggregates are regenerated downstream of those new means. For example:

```bash
nucleosuite cutn-suite \
  --rerun-from cutn_results_h3K4me3 \
  --exclude-sample wt_K4_R2_sort.bam \
  --peak-min-region-length 60 \
  --peak-max-neg-run 2 \
  --cluster-seed-p-value 0.01 \
  --cluster-member-mode significant-only \
  --min-cluster-members 2
```

The peak caller can be changed downstream with `--peak-min-region-length`, `--peak-max-neg-run`, `--peak-smooth-window`, and `--peak-smooth-order`. Stage 1 coverage statistic, seed rule, S and G gates, seed threshold, cluster membership mode, maximum non-member gap, maximum adjacent-member distance, minimum cluster member count, Stage 2 `--differential-fdr`, and the cluster aggregate/NRL settings can likewise be changed because none alters the retained per-sample BigWigs. `--skip-cluster-aggregate` disables aggregate regeneration; `--run-cluster-aggregate` explicitly enables it when the source run had skipped it.

Parameters that define the initial BigWigs are inherited and locked during a reuse rerun. These include the scoring method, resolved scoring geometry/mode inputs, discovery and coverage fragment selections, BAM grouping, contigs, blacklist and duplicate handling. If one of these options is supplied with `--rerun-from`, the command stops with an explanatory error rather than silently mixing incompatible tracks. To change such a parameter, start a fresh `cutn-suite` run from the BAMs.

Every rerun writes its own `cutn_suite_run_manifest.json`. It records the source run, excluded BAMs, retained condition manifests, and downstream parameters explicitly changed for that rerun.

## Automatic fragment mode

`--mode auto` visits indexed genomic blocks in seeded random order and samples fragments that pass the shared BAM, duplicate, contig, blacklist, and default 1–1,000 bp coverage filters. It counts the nucleosome-sized 120–250 bp subset in one-base histogram bins and bootstraps that raw histogram until its mode stabilizes or the maximum sample size is reached. Random block order avoids making the estimate depend on chromosome input order; restricting the modal search prevents abundant assay-specific short or very long fragments from defining nucleosome geometry. The resolved mode then sets the default discovery range to mode ±30 bp.

Histogram smoothing is disabled by default. `--mode-histogram-smoothing binomial` explicitly enables the optional normalized `1,4,6,4,1` kernel. Smoothing remains opt-in because it can move or merge closely spaced modes.

The treatment, control, and pooled estimates are printed as soon as each calculation finishes. `00_setup/*_fragment_mode_estimation.tsv` records the resolved mode, bootstrap interval, fragment counts, search bounds, smoothing method, stability result, and checkpoint count.

For example, a two-condition run reports lines of the form:

```text
[cutn-suite] Condition 1 treatment fragment mode: 153 bp (...)
[cutn-suite] Condition 1 control fragment mode: 151 bp (...)
[cutn-suite] Condition 2 treatment fragment mode: 154 bp (...)
[cutn-suite] Condition 2 control fragment mode: 152 bp (...)
[cutn-suite] Resolved analysis modes: treatment=153 bp; control=153 bp; strategy=pooled
```

In a two-condition run, mode estimates are pooled across corresponding groups so both Stage 1 analyses use compatible scoring geometry.

The default `--mode-strategy pooled` gives the group histograms equal weight and uses one mode for every treatment and control. `separate` uses one pooled treatment mode and one pooled control mode; `target` and `control` apply the selected group mode to both.

Automatic estimation can be bypassed:

```bash
nucleosuite cutn-suite \
  --treatment1-bam target.bam \
  --control1-bam control.bam \
  --outdir cutn_results \
  --mode 167
```

For Stage 1 analyses that will later be compared with `cutn-compare`, using the same explicit mode is the simplest way to guarantee compatibility.

## Outputs

A one-condition run writes:

- `00_setup/`: mode and normalization reports;
- `01_score_tracks/`: method-specific mode-centred score/positive-score BigWigs plus unscaled 1–1,000 bp coverage BigWigs generated in the same `tracks` pass;
- `02_mean_scaled_tracks/`: method-specific positive-score-normalized score tracks and replicate/condition-mean coverage scaled to 100;
- `03_peak_calls/`: treatment-defined nucleosome candidate peaks;
- `04_peak_statistics/`: replicate statistics, annotated treatment peaks, seed/member classifications, and seeded clusters;
- `05_cluster_aggregate/`: strongest-member anchors, replicate and combined normalized-score profiles, heatmap, bootstrap confidence band, and directional NRL outputs;
- `cutn_stage1_manifest.json`: reusable Stage 1 metadata and scaled-track paths.
- `cutn_suite_run_manifest.json`: run-level condition, BAM, parameter, and reuse metadata used by `--inspect-run` and `--rerun-from`.

A two-condition run writes the two Stage 1 trees under `01_condition1_stage1/` and `02_condition2_stage1/`. `03_condition_comparison/` contains cluster-only differential tables and BEDs, overlap-component mapping, Venn and occupied-base summaries, and matched union-locus aggregate heatmaps. The root `cutn_suite_run_manifest.json` points to both condition manifests and the comparison manifest.

When contigs run in parallel, `cutn-suite` follows each native multicontig manifest to its combined BigWig or BED output before continuing.

[Back to the command reference](../COMMAND_REFERENCE.md)
