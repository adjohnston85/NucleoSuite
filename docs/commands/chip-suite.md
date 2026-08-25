# `nucleosuite chip-suite`

## What this command does

`chip-suite` performs control-aware nucleosome-score analysis for ChIP-seq, CUT&RUN, or CUT&Tag. A run with condition 1 treatment and control BAMs performs Stage 1 only. Supplying all four condition groups performs Stage 1 independently for both conditions and then compares their enrichments in Stage 2.

TNS is the default score over 120–500 bp fragments. Use `--scoring-method bns` or `--scoring-method pns` when another kernel is required.

## Why use it

Use this workflow to distinguish antibody-target enrichment from a condition-matched control and, optionally, test whether that enrichment changes between conditions. TNS locates peaks by default. Every coverage BigWig is separately divided by its own finite non-zero mean and multiplied by 100 before peak strength is measured. The control is not globally subtracted from the treatment track.

## Stage 1: one condition

```bash
nucleosuite chip-suite \
  --treatment1-bam H3K4me3.bam \
  --control1-bam H3.bam \
  --condition1-name wild_type \
  --outdir wild_type_stage1 \
  --cores 8
```

The older names `--target-bam` and `--control-bam` remain aliases for `--treatment1-bam` and `--control1-bam`.

Stage 1 separates **peak discovery** from **peak measurement**. TNS, BNS, or PNS defines where candidate peaks occur. Scaled fragment coverage then provides the replicate values used for treatment-versus-control filtering and statistics. Keeping these roles separate avoids treating the height of a model-derived positioning score as if it were direct fragment abundance.

### 1. Generate one score and coverage set per replicate

For each treatment and control replicate, `chip-suite` generates:

- the centred TNS, BNS, or PNS track;
- its matching non-negative `posTNS`, `posBNS`, or `posPNS` track; and
- raw fragment coverage.

The centred score locates protected-DNA structure. The positive track measures the overall amount of method-specific score support and is used only as the normalization reference. Raw coverage is retained so the original sequencing-depth scale remains available.

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

The maximum is used because it captures the strongest local fragment enrichment supporting the candidate without making the score strongly dependent on the width of the TNS-defined interval. These coverage maxima—not the TNS, BNS, or PNS peak heights—are the replicate peak scores used for filtering and inference.

Column 5 of the annotated BED is replaced with the maximum of the condition-mean treatment `Cov100` track over the same interval. This provides one convenient display score, while `target_peak_replicate_statistics.tsv` retains every individual replicate value used in the test.

### 5. Require consistent treatment-over-control enrichment

A treatment peak is eligible to proceed only when

```math
\min_i T_i(R)>\max_j C_j(R).
```

This all-versus-all gate requires every treatment replicate to exceed every control replicate. It is intentionally conservative: one weak treatment or one strong control prevents the candidate from passing. Input ordering does not affect the decision because treatment and control are compared as independent groups rather than paired files.

### 6. Select Stage 1 peaks and annotate statistical evidence

A one-sided Welch test compares the full treatment and control replicate vectors for every candidate. Welch's test is used because treatment and control replicates need not have equal variances or equal group sizes. The one-sided alternative matches the Stage 1 question: is treatment greater than control? Benjamini-Hochberg correction is calculated across all candidates and retained as an annotation.

By default, the all-controls gate alone selects individual peaks for Stage 2. This default is intentional. A nucleosome-scale analysis may test hundreds of thousands of correlated candidates, while two or three biological replicates provide little power for a separate Welch test at every nucleosome. Requiring genome-wide FDR by default can therefore discard every peak even when its treatment values consistently exceed all controls. The p-value and FDR remain in the outputs so they can be inspected without being presented as stronger evidence than the replicate design supports.

`--stage1-p-value` optionally adds an exploratory raw-p cutoff, and `--peak-fdr` optionally adds an FDR cutoff. Neither is applied unless requested. The output BED preserves the discovery-track coordinates and other fields, uses the condition-mean treatment coverage maximum as column 5, and appends FDR. `target_peak_replicate_statistics.tsv` reports every replicate maximum, group means, the conservative `min(T)-max(C)` excess, conservative fold and log2 enrichment calculated with a pseudocount of 1, the gate and selection results, p-value, and FDR. At least two treatment and two control biological replicates are required to calculate p-values and FDR. Runs with fewer than three replicates in either group print a warning that these annotations are exploratory.

### 7. Form Stage 1 clusters

Clusters are built from stretches of treatment peaks that satisfy both requirements:

```math
\min_i T_i(R)>\max_j C_j(R)
\quad\text{and}\quad p(R)<0.05.
```

The first requirement ensures that every cluster member is consistently stronger in treatment than in every control. The second removes peaks whose treatment-control separation is too variable to support even nominal evidence. `--cluster-member-p-value` changes the default 0.05 member threshold. `--cluster-break`, `--max-cluster-gap`, and `--min-significant-peaks` control how qualifying members are grouped. The member p-value is deliberately not described as cluster-level FDR: it defines which nucleosomes may contribute to a cluster, while the resulting cluster is a broader feature for Stage 2 measurement.

`chip-suite` calls nucleosome peaks only. It does not call or retain breakpoint peaks because Stage 1 and Stage 2 use positive nucleosome-score candidates.

Because cluster membership requires a replicate p-value, merged mode and groups lacking two biological replicates do not produce Stage 1 clusters. Gate-selected individual peaks are still written.

## Replicates and merged input

Multiple BAMs in each group are independent biological replicates by default. Treatment and control groups may contain different numbers of BAMs. They are not paired by command-line order:

```bash
nucleosuite chip-suite \
  --treatment1-bam wt_mark_r1.bam wt_mark_r2.bam wt_mark_r3.bam \
  --control1-bam wt_H3_r1.bam wt_H3_r2.bam wt_H3_r3.bam \
  --outdir wt_stage1
```

Each replicate is scored and normalized separately before the discovery tracks are averaged. This gives the replicates equal footing during candidate discovery. Replicate-specific scaled coverage supplies the all-controls gate and exploratory one-sided Welch annotations; condition-mean treatment coverage supplies the single reported BED score. All tracks are retained in `chip_stage1_manifest.json` so Stage 2 can reuse the exact replicate measurements without returning to the BAM files.

Use `--bam-mode merged` to pass every treatment BAM as one logical treatment sample and every control BAM as one logical control sample. This matches the usual NucleoSuite multi-BAM pooling behaviour. Merged mode provides Stage 2 effect sizes and gain/loss direction, but not biological-replicate p-values or FDR.

## Stage 1 plus Stage 2

Supply condition 2 as a complete treatment/control pair:

```bash
nucleosuite chip-suite \
  --treatment1-bam wt_mark_r1.bam wt_mark_r2.bam \
  --control1-bam wt_H3_r1.bam wt_H3_r2.bam \
  --condition1-name wild_type \
  --treatment2-bam mutant_mark_r1.bam mutant_mark_r2.bam \
  --control2-bam mutant_H3_r1.bam mutant_H3_r2.bam \
  --condition2-name mutant \
  --outdir mutant_vs_wild_type \
  --cores 8
```

The Stage 2 candidate set is the union of gate-selected Stage 1 peaks, or the union of Stage 1 clusters whose members pass the gate and p < 0.05. Using the union prevents condition-specific discovery from being interpreted as zero signal in the other condition. Instead, every region is measured quantitatively in both conditions. Overlapping peaks are represented by the union of their coordinates so all four groups use one shared interval. A peak called in only one condition retains its Stage 1 interval and is still measured in every track from the other condition. `--peak-match-distance` can additionally combine nearby non-overlapping peaks.

Every differential row includes `region_origin`: `overlap_union` for matched overlapping peaks, `proximity_union` for non-overlapping peaks joined by `--peak-match-distance`, `condition1_only`, or `condition2_only`. The same value is appended to directional gain/loss BED records.

For every shared peak interval, Stage 2 retains four independent vectors of maximum scaled coverage: $T_1$, $C_1$, $T_2$ and $C_2$. Clusters use positive scaled-coverage area. Each replicate measurement is transformed as

```math
Y=\log_2(\mathrm{scaled\ coverage}+1).
```

The logarithm makes multiplicative enrichment differences comparable across low- and high-coverage regions, while the pseudocount permits zero-valued regions. A four-group factorial model then tests the condition-by-treatment interaction

```math
\Delta_{log}(R)=\left(\overline{Y}_{T_2(R)}-\overline{Y}_{C_2(R)}\right)
-\left(\overline{Y}_{T_1(R)}-\overline{Y}_{C_1(R)}\right).
```

Subtracting control within each condition accounts for condition-specific background. Comparing the two treatment-minus-control enrichments then asks whether target-specific enrichment changed between conditions rather than merely whether treatment coverage differs.

An empirical-Bayes moderated t statistic borrows information about residual variance across all Stage 2 regions. This stabilizes the very noisy region-specific variance estimates produced by only two or three replicates without treating regions as extra biological replicates. Benjamini-Hochberg correction is applied to the moderated p-values. At least two replicates in each of the four groups are required; merged or undersampled groups receive descriptive effects without p-values or FDR.

Every row is retained whether or not it passes `--differential-fdr`. The table reports raw and log-scale effects, ordinary and moderated p-values, differential FDR, direction, and replicate consistency. Significant gain/loss BEDs contain only FDR-passing regions; all-direction BEDs remain available for ranking.

For an additional assumption-free consistency annotation, Stage 2 calculates the lower and upper possible log enrichment within each condition from all treatment and control extrema. `robust_gain` means the entire condition 2 enrichment range exceeds the condition 1 range; `robust_loss` means the reverse. This stringent annotation describes complete replicate separation but does not replace the moderated FDR when making a formal differential claim.

Stage 2 never returns to the BAM files. TNS or the selected alternative score defines candidate locations; mean-scaled coverage defines their measured strength.

## Automatic fragment mode

`--mode auto` visits indexed genomic blocks in seeded random order, accumulates accepted fragment lengths, and bootstraps the raw 120–250 bp histogram until its mode stabilizes or the maximum sample size is reached. Random block order avoids making the estimate depend on chromosome input order. Only fragments passing the analysis filters enter the histogram, so the mode describes the fragments that actually contribute to scoring.

Histogram smoothing is disabled by default. `--mode-histogram-smoothing binomial` explicitly enables the optional normalized `1,4,6,4,1` kernel. Smoothing remains opt-in because it can move or merge closely spaced modes.

The treatment, control, and pooled estimates are printed as soon as each calculation finishes. `00_setup/*_fragment_mode_estimation.tsv` records the resolved mode, bootstrap interval, fragment counts, search bounds, smoothing method, stability result, and checkpoint count.

For example, a two-condition run reports lines of the form:

```text
[chip-suite] Condition 1 treatment fragment mode: 153 bp (...)
[chip-suite] Condition 1 control fragment mode: 151 bp (...)
[chip-suite] Condition 2 treatment fragment mode: 154 bp (...)
[chip-suite] Condition 2 control fragment mode: 152 bp (...)
[chip-suite] Resolved analysis modes: treatment=153 bp; control=153 bp; strategy=pooled
```

In a two-condition run, mode estimates are pooled across corresponding groups so both Stage 1 analyses use compatible scoring geometry.

The default `--mode-strategy pooled` gives the group histograms equal weight and uses one mode for every treatment and control. `separate` uses one pooled treatment mode and one pooled control mode; `target` and `control` apply the selected group mode to both.

Automatic estimation can be bypassed:

```bash
nucleosuite chip-suite \
  --target-bam target.bam \
  --control-bam control.bam \
  --outdir chip_results \
  --mode 167
```

For Stage 1 analyses that will later be compared with `chip-compare`, using the same explicit mode is the simplest way to guarantee compatibility.

## Output layout

A one-condition run writes:

- `00_setup/`: mode and normalization reports;
- `01_score_tracks/`: score, positive-score, and unscaled raw coverage BigWigs;
- `02_mean_scaled_tracks/`: positive-score-normalized discovery tracks plus replicate and condition-mean coverage scaled to 100;
- `03_peak_calls/`: treatment-defined nucleosome candidate peaks;
- `04_peak_fdr/`: replicate statistics, annotated and gate-selected peaks, and p-defined clusters;
- `chip_stage1_manifest.json`: reusable Stage 1 metadata and scaled-track paths.

A two-condition run writes the two Stage 1 trees under `01_condition1_stage1/` and `02_condition2_stage1/`, then writes differential peak and cluster tables under `03_condition_comparison/`.

When contigs run in parallel, `chip-suite` follows each native multicontig manifest to its combined BigWig or BED output before continuing.

[Back to the command reference](../COMMAND_REFERENCE.md)
