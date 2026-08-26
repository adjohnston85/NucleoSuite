# `nucleosuite chip-compare`

## What this command does

`chip-compare` performs cluster-only Stage 2 analysis from two completed `chip-suite` Stage 1 manifests. It compares target-specific cluster enrichment between biological conditions, summarizes observed cluster overlap and occupied bases, and creates coordinate-matched cluster-centred method-specific score aggregates. It reads saved BigWigs and does not revisit BAM files.

## Why use it

Use it when the two conditions were processed separately or when Stage 2 needs to be repeated with another differential FDR cutoff.

```bash
nucleosuite chip-compare \
  --condition1-results wild_type_stage1 \
  --condition2-results mutant_stage1 \
  --outdir mutant_vs_wild_type \
  --fdr 0.05
```

Each results argument may name a Stage 1 directory or its `chip_stage1_manifest.json`.

## How and why the comparison works

### 1. Validate comparable Stage 1 analyses

The scoring method, resolved treatment/control modes, fragment limits, contig selection, Stage 1 measurement method, statistical method, and BigWig chromosome definitions must agree. This prevents differential results from mixing callsets discovered with different scoring geometries or measurements made on incompatible coordinate systems.

### 2. Construct overlap-connected cluster loci

Stage 2 compares clusters, not individual nucleosome peaks. All selected Stage 1 clusters from both conditions are sorted by coordinate. Directly or transitively overlapping clusters form one **cluster locus**. Non-overlapping clusters remain separate even if they are close.

This connected-component rule preserves genuine one-to-many and many-to-many relationships. For example:

```text
condition 1:  └──────── cluster A ────────┘
condition 2:    └─ cluster B ─┘  └─ C ─┘
Stage 2 locus: └────────── A+B+C ──────────┘   (1-to-many)
```

Forcing a one-to-one match would discard either B or C or make the result depend on an arbitrary matching order. `cluster_overlap_components.tsv` therefore retains every contributing cluster ID and labels shared loci as `1_to_1`, `1_to_many`, `many_to_1`, or `many_to_many`.

`region_origin` is:

- `overlap_union` when the component contains clusters from both conditions;
- `condition1_only` when it contains only condition 1 clusters; or
- `condition2_only` when it contains only condition 2 clusters.

A condition-specific cluster is still measured in every BigWig from the other condition. Failure to call a cluster is not treated as proof of zero signal.

### 3. Measure cluster enrichment in four replicate groups

Each Stage 1 coverage BigWig was independently scaled to a non-zero mean of 100. For every cluster locus $R$, `chip-compare` sums only positive scaled coverage in each replicate:

```math
A_i(R)=\sum_{x\in R}\max(Cov_{100,i}(x),0).
```

Positive area is used because a cluster is an extended domain: both signal magnitude and enriched span are relevant. The four independent replicate vectors are condition 1 treatment ($T_1$), condition 1 control ($C_1$), condition 2 treatment ($T_2$), and condition 2 control ($C_2$). Treatment and control files are not paired by command-line order and group sizes may differ.

Each area is transformed as

```math
Y=\log_2(A+1).
```

The logarithm makes multiplicative differences more comparable between weak and strong loci. The pseudocount permits loci with zero positive area.

### 4. Test the treatment-by-condition interaction

A four-group factorial model contains condition, treatment/control status, and their interaction. The tested coefficient is

```math
\Delta=(\bar Y_{T_2}-\bar Y_{C_2})-(\bar Y_{T_1}-\bar Y_{C_1}).
```

The first subtraction estimates target-specific enrichment within each condition. The second subtraction asks whether that target-over-control enrichment changes between conditions. This avoids mistaking a general background or coverage change for a target-specific biological change.

The ordinary p-value is a two-sided t test of the interaction coefficient using the locus-specific residual variance and factorial-model residual degrees of freedom. A positive coefficient is a gain in condition 2 and a negative coefficient is a loss.

With only two or three replicates, locus-specific variance estimates are noisy. NucleoSuite therefore estimates a scaled-inverse-chi-square variance prior across all cluster loci and calculates

```math
s_{post,R}^2=\frac{d_0s_0^2+d_Rs_R^2}{d_0+d_R}.
```

The moderated p-value uses the posterior variance $s_{post,R}^2$ and $d_0+d_R$ degrees of freedom. This borrows information about variance across loci without treating loci as biological replicates. Benjamini-Hochberg correction is applied across all moderated cluster-locus p-values. At least two biological replicates are required in every group; otherwise effect sizes and directions are reported without inferential p-values or FDR.

The complete table also reports ordinary and moderated standard errors and 95% confidence intervals. `--fdr 0.05` controls only the separate significant gain/loss BEDs; every locus remains in `differential_clusters.tsv` and the all-direction BEDs.

### 5. Report a conservative replicate-separation label

For condition $k$,

```math
L_k=\min(Y_{T_k})-\max(Y_{C_k}),
\qquad
U_k=\max(Y_{T_k})-\min(Y_{C_k}).
```

`robust_gain` requires $L_2>U_1$; `robust_loss` requires $U_2<L_1$. These labels mean every possible treatment-control contrast is ordered in the same direction. They are descriptive consistency annotations, not replacements for the moderated p-value or FDR.

### 6. Summarize observed cluster overlap

`cluster_locus_venn.png` shows condition 1-only, shared, and condition 2-only **overlap-connected loci**. The Venn uses loci rather than raw cluster counts because one cluster can overlap several clusters in the other condition. `cluster_overlap_summary.tsv` separately reports:

- raw cluster counts for each condition;
- clusters with any cross-condition overlap;
- shared and condition-only locus counts;
- one-to-one, one-to-many, many-to-one, and many-to-many shared loci;
- non-overlapping bases occupied by each condition's cluster union;
- overlapping, condition-unique, and union cluster bases;
- the percentage of each condition's occupied bases that overlaps; and
- base-pair Jaccard percentage.

These are descriptive summaries of the observed callsets. This build does not perform a genomic randomization test of whether overlap exceeds chance.

### 7. Create matched cluster-centred method-specific score aggregates

For each shared or condition-specific locus, the common anchor is the strongest coverage-scored member peak among all contributing Stage 1 clusters. Both conditions are aligned to the same ordered anchors, so heatmap rows refer to identical loci.

Each treatment replicate's selected score was independently divided by the finite, non-zero mean of its matching positive-score track before averaging. PNS uses `posPNS`, BNS uses `posBNS`, and TNS uses `posTNS`. Scaling before averaging prevents a deeper replicate from determining the condition mean. The selected Stage 1 scoring method is reused for this positioning view; scaled coverage remains the statistical measurement.

The matched outputs include condition-specific heatmaps with one common symmetric colour range, replicate and replicate-combined mean profiles, cluster-bootstrap 95% bands, and directional NRLs. Defaults are ±1,000 bp around the anchor, 140 bp peak resolution, central peak order 0 included, and regression through peak orders 0–3 with no central exclusion. Each Stage 1 directory also contains an own-cluster aggregate; the Stage 2 matched aggregate is intended for direct visual comparison between conditions.

## Outputs

The output directory contains:

- `differential_clusters.tsv`;
- all, robust, and FDR-significant gain/loss BEDs;
- `cluster_overlap_components.tsv`;
- `cluster_overlap_summary.tsv` and `cluster_locus_venn.png`;
- `cluster_aligned_aggregates/` with common anchors, heatmaps, profiles, confidence bands, and NRL outputs; and
- `chip_comparison_manifest.json` with all paths, model metadata, overlap counts, and aggregate settings.

[Back to the command reference](../COMMAND_REFERENCE.md)
