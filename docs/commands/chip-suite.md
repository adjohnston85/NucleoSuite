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

Stage 1 generates the score, corresponding positive-score track, and raw coverage for each logical sample. It divides the centred score by its positive-score mean for peak discovery. It also writes a second coverage BigWig with non-zero mean 100:

```math
Cov_{100}(x)=100\frac{Cov(x)}{\mathrm{mean}(Cov(x)\mid Cov(x)>0)}.
```

With the default `--stage1-control-mode all-controls`, peaks are called only on the condition-mean treatment score track. Control peaks are not called. Every treatment candidate interval is then measured separately in every treatment and control replicate. Define

```math
T_i(R)=\max_{x\in R}Cov_{100,T_i}(x),
\qquad
C_j(R)=\max_{x\in R}Cov_{100,C_j}(x).
```

A treatment peak is taken forward only when

```math
\min_i T_i(R)>\max_j C_j(R).
```

This is an all-versus-all rule: every treatment score must exceed every control score, and input ordering does not affect the Stage 1 decision. A one-sided Welch test compares the complete treatment and control replicate vectors for every treatment candidate. Benjamini-Hochberg correction is applied across all treatment candidates, including candidates that fail the all-controls gate. A peak enters the significant BED only when it passes both the gate and `--peak-fdr`.

The output BED preserves the TNS-defined coordinates and other fields, replaces column 5 with the maximum condition-mean treatment `Cov100` in the interval, and appends FDR. `target_peak_replicate_statistics.tsv` reports every replicate maximum, group means, the conservative `min(T)-max(C)` excess, gate result, p-value and FDR. At least two treatment and two control biological replicates are required for inferential p-values and FDR; otherwise the annotated BED contains `.` and no peak enters the FDR-filtered BED. `--stage1-control-mode condition-mean` remains as a legacy compatibility mode and retains the older control-peak empirical-decoy analysis.

## Replicates and merged input

Multiple BAMs in each group are independent biological replicates by default. Treatment and control groups may contain different numbers of BAMs. They are not paired by command-line order:

```bash
nucleosuite chip-suite \
  --treatment1-bam wt_mark_r1.bam wt_mark_r2.bam wt_mark_r3.bam \
  --control1-bam wt_H3_r1.bam wt_H3_r2.bam wt_H3_r3.bam \
  --outdir wt_stage1
```

Each replicate is scored and normalized separately. The condition-mean treatment score BigWig is used for Stage 1 peak discovery. Replicate-specific scaled coverage supplies the all-controls gate and one-sided Welch test; condition-mean treatment coverage supplies the single reported BED score. All tracks are retained in `chip_stage1_manifest.json` for later differential inference.

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

The Stage 2 candidate set is the union of significant Stage 1 features from either condition. Overlapping peaks are represented by the union of their coordinates. A peak called in only one condition retains its Stage 1 interval and is still measured in every track from the other condition. `--peak-match-distance` can additionally combine nearby non-overlapping peaks.

Every differential row includes `region_origin`: `overlap_union` for matched overlapping peaks, `proximity_union` for non-overlapping peaks joined by `--peak-match-distance`, `condition1_only`, or `condition2_only`. The same value is appended to directional gain/loss BED records.

For every shared interval, Stage 2 retains four independent vectors of maximum scaled coverage: $T_1$, $C_1$, $T_2$ and $C_2$. The reported interaction is

```math
\Delta(R)=\left(\overline{T_2(R)}-\overline{C_2(R)}\right)
-\left(\overline{T_1(R)}-\overline{C_1(R)}\right).
```

A Welch-style difference-of-differences test uses the variance and sample size of all four groups without pairing treatment and control BAMs. Benjamini-Hochberg correction is applied across all Stage 2 regions. At least two replicates in each of the four groups are required for inferential FDR. Peaks use maximum BigWig values; clusters use positive base-wise area.

Stage 2 never returns to the BAM files. TNS or the selected alternative score defines candidate locations; mean-scaled coverage defines their measured strength.

## Automatic fragment mode

`--mode auto` visits indexed genomic blocks in seeded random order, accumulates accepted fragment lengths, and bootstraps a lightly smoothed 120–250 bp histogram until its mode stabilizes or the maximum sample size is reached. In a two-condition run, mode estimates are pooled across corresponding groups so both Stage 1 analyses use compatible scoring geometry.

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
- `03_peak_calls/`: treatment candidates, plus control candidates only in legacy condition-mean mode;
- `04_peak_fdr/`: replicate statistics, annotated/significant peaks, and clusters;
- `chip_stage1_manifest.json`: reusable Stage 1 metadata and scaled-track paths.

A two-condition run writes the two Stage 1 trees under `01_condition1_stage1/` and `02_condition2_stage1/`, then writes differential peak and cluster tables under `03_condition_comparison/`.

When contigs run in parallel, `chip-suite` follows each native multicontig manifest to its combined BigWig or BED output before continuing.

[Back to the command reference](../COMMAND_REFERENCE.md)
