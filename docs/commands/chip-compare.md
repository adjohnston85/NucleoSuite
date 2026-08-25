# `nucleosuite chip-compare`

## What this command does

`chip-compare` performs Stage 2 separately from BAM processing. It reads two completed `chip-suite` Stage 1 manifests, validates their scoring compatibility, constructs shared peak and/or cluster intervals, and tests the treatment-by-condition interaction using the recorded coverage BigWigs scaled to a non-zero mean of 100.

## Why use it

Use this command when the two conditions were processed in separate Stage 1 runs, when Stage 2 needs to be repeated with a different feature set or cutoff, or when the original BAMs should not be revisited.

## Typical run

```bash
nucleosuite chip-compare \
  --condition1-results wild_type_stage1 \
  --condition2-results mutant_stage1 \
  --outdir mutant_vs_wild_type
```

Each argument may name a Stage 1 directory or its `chip_stage1_manifest.json` directly. The scoring method, treatment and control modes, fragment limits, contig selection, and BigWig chromosome definitions must match.

## Comparison statistics

### 1. Validate compatible Stage 1 analyses

The command first checks the scoring method, resolved modes, fragment limits, contigs, peak-measurement method, and BigWig chromosome definitions. These settings must agree because otherwise the two Stage 1 feature sets would have been discovered or measured using different geometries.

### 2. Construct one shared region set

The candidate set is the union of gate-selected Stage 1 peaks from either condition. At cluster level, it is the union of Stage 1 clusters made from peaks that passed both the all-controls gate and the configured member p-value. The union is used because failure to select a feature in one condition does not prove that its signal is zero. Retaining the region allows its actual scaled coverage to be measured in both conditions.

Overlapping peaks use the union of their Stage 1 coordinates so every replicate is summarized across one common interval. Condition-specific peaks retain their own coordinates. `--peak-match-distance` can additionally combine nearby non-overlapping peaks that plausibly represent the same shifted feature.

The `region_origin` field records how each comparison interval was formed:

| Value | Meaning |
| --- | --- |
| `overlap_union` | overlapping Stage 1 regions from both conditions were combined |
| `proximity_union` | non-overlapping peaks were combined using `--peak-match-distance` |
| `condition1_only` | the region came only from condition 1 Stage 1 |
| `condition2_only` | the region came only from condition 2 Stage 1 |

`--feature-level peaks`, `clusters`, or `both` controls which outputs are produced.

### 3. Measure all four replicate groups

Every region produces four independent replicate vectors: condition 1 treatment, condition 1 control, condition 2 treatment, and condition 2 control. Peaks use maximum scaled coverage because it captures local enrichment without strongly depending on interval width. Clusters use positive scaled-coverage area because a broader cluster represents both signal magnitude and enriched extent.

### 4. Transform replicate measurements

Each replicate measurement is transformed before inference:

```math
Y=\log_2(\mathrm{scaled\ coverage}+1).
```

Coverage enrichment is multiplicative and tends to become more variable as its magnitude increases. The log transformation makes fold-like changes comparable across the signal range. Adding 1 permits a region with zero measured coverage while having little effect near the mean-scaled value of 100.

### 5. Test the condition-by-treatment interaction

The factorial model contains condition, treatment/control status, and their interaction. Its tested coefficient is

```math
\Delta_{log}=(\bar Y_{T_2}-\bar Y_{C_2})-(\bar Y_{T_1}-\bar Y_{C_1}).
```

Control subtraction within each condition accounts for condition-specific background. Comparing those two enrichments tests whether target-specific enrichment changes between conditions rather than simply testing whether the treatment tracks differ.

Ordinary region-by-region variance estimates are unreliable with two or three replicates. `chip-compare` therefore estimates a shared variance prior across the full region set and combines it with each region's residual variance. The resulting empirical-Bayes moderated t test borrows precision across regions without treating them as biological replicates. Both the ordinary and moderated p-values are reported, but Benjamini-Hochberg correction uses the moderated p-values.

At least two replicates are required in every group. Merged or undersampled inputs receive descriptive log and raw effects without differential p-values or FDR. `--fdr 0.05` controls the separate significant gain and loss BEDs; the complete TSV and all-direction BEDs are always written.

### 6. Annotate replicate-consistent direction

For condition $k$, the possible log enrichment range is

```math
L_k=\min(Y_{T_k})-\max(Y_{C_k}),
\qquad
U_k=\max(Y_{T_k})-\min(Y_{C_k}).
```

`robust_gain` requires $L_2>U_1$, so every treatment-control contrast in condition 2 exceeds every contrast in condition 1. `robust_loss` requires $U_2<L_1$. Otherwise the region is `not_replicate_separated`. This stringent annotation is useful for 2×2 designs, but it is descriptive and does not replace moderated FDR for a formal differential call.

## Outputs

The output directory contains `differential_peaks.tsv` and/or `differential_clusters.tsv`, directional BED files, and `chip_comparison_manifest.json`. Tables include Stage 1 support, `region_origin`, every raw replicate score, raw and log group means, raw and log interaction effects, enrichment-range bounds, effect direction, replicate consistency, ordinary and moderated p-values, differential FDR, residual and posterior variances, and status.

For each feature level, the command writes all gains, all losses, robust gains, robust losses, FDR-significant gains, and FDR-significant losses as separate BEDs. Directional BED records append `region_origin` after differential FDR. The comparison manifest records the empirical-Bayes prior, model degrees of freedom, output paths, and counts by origin.

[Back to the command reference](../COMMAND_REFERENCE.md)
