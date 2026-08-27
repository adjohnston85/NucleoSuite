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

### 4. Assign replicate peak scores

Every treatment-defined candidate interval is measured separately in every treatment and control scaled-coverage track. For interval $R$:

```math
T_i(R)=\max_{x\in R}Cov_{100,T_i}(x),
\qquad
C_j(R)=\max_{x\in R}Cov_{100,C_j}(x).
```

The maximum is used because it captures the strongest local fragment enrichment supporting the candidate without making the score strongly dependent on the width of the SNS/PNS/BNS/TNS-defined interval. These broad-range coverage maxima—not the SNS, PNS, BNS, or TNS peak heights—are the replicate peak scores used for filtering and inference.

Column 5 of the annotated BED is replaced with the maximum of the condition-mean treatment `Cov100` track over the same interval. This provides one convenient display score. The complete peak BED appends the raw Welch p-value and BH FDR as its final two columns, while the parameter-aware replicate-statistics TSV retains every individual replicate value used in the test.

### 5. Apply the treatment-control gate

The default gate is the conservative **all-controls** rule:

```math
\min_i T_i(R)>\max_j C_j(R),
```

which requires every treatment replicate to exceed every control replicate. Use `--stage1-gate-mode mean` to instead require

```math
\mathrm{mean}_i T_i(R)>\mathrm{mean}_j C_j(R).
```

The mean gate allows replicate-level variation while requiring average treatment enrichment. Treatment and control are independent groups rather than paired files, so input ordering does not affect either gate.

### 6. Select Stage 1 peaks and annotate statistical evidence

A one-sided Welch test compares the full treatment and control replicate vectors for every candidate. Welch's test is used because treatment and control replicates need not have equal variances or equal group sizes. The one-sided alternative matches the Stage 1 question: is treatment greater than control? Benjamini-Hochberg correction is calculated across all candidates and retained as an annotation.

By default, the selected treatment-control gate alone selects individual peaks for Stage 2. This default is intentional. A nucleosome-scale analysis may test hundreds of thousands of correlated candidates, while two or three biological replicates provide little power for a separate Welch test at every nucleosome. Requiring genome-wide FDR by default can therefore discard every peak even when treatment is consistently enriched over control. The p-value and FDR remain in the outputs so they can be inspected without being presented as stronger evidence than the replicate design supports.

`--stage1-p-value` optionally adds an exploratory raw-p cutoff, and `--peak-fdr` optionally adds an FDR cutoff. Neither is applied unless requested. The complete annotated BED is named `target_peaks_replicate_statistics_gate_<mode>.bed`; it preserves the discovery-track coordinates and other fields, uses the condition-mean treatment coverage maximum as column 5, and appends **raw p-value then BH FDR**. The selected-peak BED carries the same two statistical columns. `target_peak_replicate_statistics_gate_<mode>.tsv` reports every replicate maximum, group means, mean treatment-minus-control difference, the conservative `min(T)-max(C)` excess and fold/log2 enrichment with pseudocount 1, both gate results, the selected gate/excess, Stage 2 selection status, p-value, and FDR. At least two treatment and two control biological replicates are required to calculate p-values and FDR. Runs with fewer than three replicates in either group print a warning that these annotations are exploratory.

### 7. Form Stage 1 clusters

Clusters are seeded by the strongest statistical evidence and then extended through consistently treatment-enriched neighbouring peaks. A **seed** must pass the selected treatment-control gate and satisfy

```math
p(R)<0.05.
```

With the default all-controls gate this means every treatment replicate exceeds every control replicate; with `--stage1-gate-mode mean` it means mean treatment > mean control. Every seed is also written to `target_seed_peaks_gate_<mode>_seed_p<threshold>.bed`, with raw p-value and BH FDR appended as the final two columns. The nominal p-value supplies a reproducible point from which to begin the cluster. Change the default seed threshold with `--cluster-seed-p-value`.

After finding a seed, `cutn-suite` looks upstream and downstream through the ordered treatment candidates. Any neighbouring peak that passes the selected treatment-control gate can extend the cluster even when its own p-value is at least 0.05. This separates the evidence needed to **start** a domain from the evidence used to define its **extent**: local replicate variability should not cut a coherent run of consistently treatment-over-control nucleosomes into many small pieces.

The default permits one consecutive non-member candidate to bridge two included members (`--cluster-max-non-member-gap 1`). In `seed-and-gated` mode, non-members are gate-failing `x` peaks; in `significant-only` mode, both `G` and `x` are non-members. Bridging candidates never become endpoints or score contributors. A run of non-members longer than the configured limit ends the current cluster; later eligible peaks can start a new seeded cluster. Adjacent included-member summits must also be no more than 1,000 bp apart (`--max-cluster-gap 1000`).

The notation in these diagrams is:

```text
S = seed: selected treatment-control gate passes and p < 0.05
G = extension: selected treatment-control gate passes; no p-value requirement
x = selected treatment-control gate fails
. = not a cluster member
```

The smallest default cluster contains one seed and one additional gated member:

```text
state:    S G
cluster:  1 1
```

One non-gated candidate may bridge two gated members. It lies within the BED span but is not a member:

```text
state:    S x G
cluster:  1 . 1
span:     └───┘
```

Two consecutive non-gated candidates split the gated runs. Cluster boundaries remain on the outermost gated members, not on either `x`:

```text
state:    S G x x G S G
cluster:  1 1 . . 2 2 2
```

Gate-passing peaks without a seed do not form a cluster:

```text
state:    G S G x x G G
cluster:  1 1 1 . . . .
```

The right-hand `G G` is discarded because it contains no `S`. A lone `S` is also discarded because `--min-cluster-members` defaults to 2. Expansion from nearby seeds produces one cluster when the expansions connect; it does not emit overlapping duplicate clusters.

With `--cluster-member-mode significant-only`, `G` and `x` are both non-members for gap counting. For example, with the default `--cluster-max-non-member-gap 1`, `S G S` may form one cluster, but `S G x S` is split because two consecutive non-members separate the significant peaks.

Cluster coordinates run from the start of the first included member to the end of the last included member. With the default `--stage1-gate-mode all-controls`, each member contributes `minimum treatment - maximum control` to the cluster score. With `--stage1-gate-mode mean`, each member instead contributes `mean treatment - mean control`. The default `--cluster-member-mode seed-and-gated` treats both `S` and `G` peaks as members; `--cluster-member-mode significant-only` restricts membership and scoring to `S` peaks. `--cluster-max-non-member-gap` controls how many consecutive non-members may bridge included members; a longer run ends the current cluster and later eligible peaks are evaluated as a new cluster. Bridging candidates never contribute to the cluster boundary or score. The strongest peak is the included member with the largest maximum on the condition-mean treatment coverage track; selected treatment-over-control excess and genomic position break ties. `--cluster-fdr` can optionally filter the maximum seed FDR, but no cluster FDR cutoff is applied by default.

`cutn-suite` calls nucleosome peaks only. It does not call or retain breakpoint peaks because Stage 1 and Stage 2 use positive nucleosome-score candidates.

Because every cluster requires at least one p-value-defined seed, merged mode and groups lacking two biological replicates do not produce Stage 1 clusters. Gate-selected individual peaks are still written.

### 8. Aggregate the selected nucleosome score around the strongest peak in each cluster

For each treatment replicate, the selected score is divided by the finite, non-zero mean of its matching positive-score track:

```math
S_{scaled,i}(x)=\frac{S_i(x)}{\mathrm{mean}(posS_i(x)\mid posS_i(x)>0)}.
```

This normalization is done per replicate before averaging because raw score magnitude increases with usable fragment depth. Averaging independently normalized method-matched score tracks gives each replicate equal weight. This differs from coverage-to-100 scaling: scaled coverage measures peak abundance and supports treatment/control statistics, whereas the normalized SNS, PNS, BNS or TNS signal shows nucleosome positioning around the selected cluster anchor.

Each cluster is aligned at the SNS/PNS/BNS/TNS summit of its strongest coverage-scored member. Keeping the discovery summit rather than replacing it with the coordinate of the coverage maximum preserves the nucleosome-position estimate while using direct coverage only to decide which member is strongest.

The default aggregate window is ±1,000 bp. Outputs include a replicate-combined heatmap and mean profile, individual replicate mean profiles and an overlay, a cluster-bootstrap 95% confidence band, and positive- and negative-direction NRL fits. Directional NRL calling defaults to 130 bp resolution, includes the aligned central peak as order 0, uses peak orders 0 through 3, and disables the usual central regression exclusion. Missing peak orders are not renumbered. Change these settings with the `--cluster-aggregate-*` options or use `--skip-cluster-aggregate` when only Stage 1 peak and cluster tables are needed.

## Replicates and merged input

Multiple BAMs in each group are independent biological replicates by default. Treatment and control groups may contain different numbers of BAMs. They are not paired by command-line order:

```bash
nucleosuite cutn-suite \
  --treatment1-bam wt_mark_r1.bam wt_mark_r2.bam wt_mark_r3.bam \
  --control1-bam wt_H3_r1.bam wt_H3_r2.bam wt_H3_r3.bam \
  --outdir wt_stage1
```

Each replicate is scored and normalized separately before the discovery tracks are averaged. This gives the replicates equal footing during candidate discovery. Replicate-specific scaled coverage supplies the selected treatment-control gate and exploratory one-sided Welch annotations; condition-mean treatment coverage supplies the single reported BED score. All tracks are retained in `cutn_stage1_manifest.json` so Stage 2 can reuse the exact replicate measurements without returning to the BAM files.

Use `--bam-mode merged` to pass every treatment BAM as one logical treatment sample and every control BAM as one logical control sample. This matches the usual NucleoSuite multi-BAM pooling behaviour. Merged mode provides Stage 2 effect sizes and gain/loss direction, but not biological-replicate p-values or FDR.

## Stage 1 plus Stage 2

Supply condition 2 as a complete treatment/control pair:

```bash
nucleosuite cutn-suite \
  --treatment1-bam wt_mark_r1.bam wt_mark_r2.bam \
  --control1-bam wt_H3_r1.bam wt_H3_r2.bam \
  --condition1-name wild_type \
  --treatment2-bam mutant_mark_r1.bam mutant_mark_r2.bam \
  --control2-bam mutant_H3_r1.bam mutant_H3_r2.bam \
  --condition2-name mutant \
  --outdir mutant_vs_wild_type \
  --cores 8
```

Stage 2 compares **clusters only**. Individual Stage 1 peaks remain available for inspection and determine each cluster's aggregate anchor, but they are not separate differential tests.

All selected clusters from both conditions are combined into overlap-connected components called **cluster loci**. Any directly or transitively overlapping clusters form one locus. This preserves one-to-many and many-to-many relationships instead of forcing an arbitrary one-to-one match. Non-overlapping clusters remain separate even when they are close. A locus is labelled `overlap_union`, `condition1_only`, or `condition2_only`; the contributing cluster IDs and relationship class are retained in `cluster_overlap_components.tsv`.

Using the union prevents a cluster absent from one Stage 1 callset from being treated as zero signal. Its interval is still measured in every treatment and control replicate from both conditions. For each cluster locus $R$, every replicate contributes the **positive area** of its scaled-coverage track:

```math
A_i(R)=\sum_{x\in R}\max(Cov_{100,i}(x),0).
```

Positive area is used rather than the maximum used for individual peaks because a cluster represents an extended enriched domain: the statistic should retain both signal magnitude and occupied extent. The four independent vectors are condition 1 treatment ($T_1$), condition 1 control ($C_1$), condition 2 treatment ($T_2$), and condition 2 control ($C_2$). Each area is transformed as

```math
Y=\log_2(A(R)+1).
```

The logarithm makes multiplicative enrichment differences comparable across low- and high-coverage regions, while the pseudocount permits zero-valued regions. A four-group factorial model then tests the condition-by-treatment interaction

```math
\Delta_{log}(R)=\left(\overline{Y}_{T_2(R)}-\overline{Y}_{C_2(R)}\right)
-\left(\overline{Y}_{T_1(R)}-\overline{Y}_{C_1(R)}\right).
```

Subtracting control within each condition accounts for condition-specific background. Comparing the two treatment-minus-control enrichments then asks whether target-specific enrichment changed between conditions rather than merely whether treatment coverage differs.

For each cluster locus, the ordinary p-value comes from the t statistic for the interaction coefficient using that locus's residual variance and the factorial-model residual degrees of freedom. Because two or three replicates produce unstable locus-specific variances, an empirical-Bayes step estimates a shared variance prior across all cluster loci and combines it with each locus's residual variance. The moderated p-value uses this posterior variance and the combined prior-plus-residual degrees of freedom. It therefore stabilizes variance estimation without increasing the biological replicate count. Benjamini-Hochberg correction is applied across the moderated cluster-locus p-values. At least two replicates in each of the four groups are required; merged or undersampled groups receive descriptive effects without p-values or FDR.

Every row is retained whether or not it passes `--differential-fdr`. The table reports the contributing cluster IDs, raw and log-scale effects, ordinary and moderated p-values, standard errors, 95% confidence intervals, differential FDR, direction, and replicate consistency. Significant gain/loss BEDs contain only FDR-passing cluster loci; all-direction BEDs remain available for ranking.

For an additional assumption-free consistency annotation, Stage 2 calculates the lower and upper possible log enrichment within each condition from all treatment and control extrema. `robust_gain` means the entire condition 2 enrichment range exceeds the condition 1 range; `robust_loss` means the reverse. This stringent annotation describes complete replicate separation but does not replace the moderated FDR when making a formal differential claim.

Stage 2 also writes a descriptive Venn diagram of condition-only and shared cluster loci. The companion summary reports raw cluster counts, shared/condition-only locus counts, one-to-one and one-to-many topology counts, bases occupied in each condition, overlapping and union bases, condition-specific overlap percentages, and base-pair Jaccard percentage. These summaries describe the observed callsets; no genomic randomization overlap test is performed.

For a coordinate-matched visual comparison, Stage 2 aligns both conditions to the same union-locus anchor set and uses the same symmetric heatmap colour range. The anchor is the strongest coverage-scored Stage 1 member peak among clusters contributing to that locus. These matched heatmaps supplement the condition's own-cluster aggregates generated during Stage 1.

Stage 2 never returns to the BAM files. The saved scaled-coverage BigWigs supply differential measurements, and the saved method-matched normalized SNS, PNS, BNS or TNS tracks supply aggregate positioning signal.

## Inspect a completed run

Completed `cutn-suite` directories can be inspected without opening the manifests manually:

```bash
nucleosuite cutn-suite --inspect-run cutn_results_h3K4me3
```

The report lists each biological condition, its treatment and control BAMs, the retained replicate tracks, scoring method, resolved treatment/control modes, discovery and coverage fragment ranges, Stage 1 gate and cluster settings. For each retained replicate, `cutn-suite` also reports lightweight sample statistics from the saved track outputs, including the fragment-length mode within the nucleosome mode-search range, the number of fragments used in the broad coverage range, the positive-score normalization mean, and the pre-scaling non-zero coverage mean when those values are available. The original treatment/control/pooled bootstrap mode estimate and confidence interval are also shown when the mode report is present.

Runs created by 0.10.13 and later contain a root `cutn_suite_run_manifest.json` that records the condition manifests and run-level settings. `--inspect-run` also recognizes the 0.10.12 layout directly, so an existing run does not need to be repeated merely to inspect it.

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

The peak caller can be changed downstream with `--peak-min-region-length`, `--peak-max-neg-run`, `--peak-smooth-window`, and `--peak-smooth-order`. Stage 1 gate/p-value/FDR settings, seed threshold, cluster FDR, cluster membership mode, maximum non-member gap, maximum adjacent-member distance, minimum cluster member count, Stage 2 `--differential-fdr`, and the cluster aggregate/NRL settings can likewise be changed because none alters the retained per-sample BigWigs. `--skip-cluster-aggregate` disables aggregate regeneration; `--run-cluster-aggregate` explicitly enables it when the source run had skipped it.

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
- `04_peak_fdr/`: replicate statistics, annotated and gate-selected peaks, and seeded clusters;
- `05_cluster_aggregate/`: strongest-member anchors, replicate and combined normalized-score profiles, heatmap, bootstrap confidence band, and directional NRL outputs;
- `cutn_stage1_manifest.json`: reusable Stage 1 metadata and scaled-track paths.
- `cutn_suite_run_manifest.json`: run-level condition, BAM, parameter, and reuse metadata used by `--inspect-run` and `--rerun-from`.

A two-condition run writes the two Stage 1 trees under `01_condition1_stage1/` and `02_condition2_stage1/`. `03_condition_comparison/` contains cluster-only differential tables and BEDs, overlap-component mapping, Venn and occupied-base summaries, and matched union-locus aggregate heatmaps. The root `cutn_suite_run_manifest.json` points to both condition manifests and the comparison manifest.

When contigs run in parallel, `cutn-suite` follows each native multicontig manifest to its combined BigWig or BED output before continuing.

[Back to the command reference](../COMMAND_REFERENCE.md)
