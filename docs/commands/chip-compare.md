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

The candidate set is the union of significant Stage 1 calls from either condition. Those calls were discovered with TNS by default. Overlapping peaks use the union of their Stage 1 coordinates. Condition-specific peaks retain their own coordinates and are measured in every replicate from both conditions; absence of a called peak is never treated as zero. `--peak-match-distance` can additionally combine nearby non-overlapping peaks. Peak scores use maximum scaled coverage, while clusters use scaled positive coverage area.

The `region_origin` field records how each comparison interval was formed:

| Value | Meaning |
| --- | --- |
| `overlap_union` | overlapping Stage 1 regions from both conditions were combined |
| `proximity_union` | non-overlapping peaks were combined using `--peak-match-distance` |
| `condition1_only` | the region came only from condition 1 Stage 1 |
| `condition2_only` | the region came only from condition 2 Stage 1 |

`--feature-level peaks`, `clusters`, or `both` controls which outputs are produced.

Every region produces four independent replicate vectors: condition 1 treatment, condition 1 control, condition 2 treatment and condition 2 control. The tested interaction is $(T_2-C_2)-(T_1-C_1)$. A Welch-style standard error and degrees of freedom use the variance and sample size of all four groups, so input-order pairing and equal group sizes are unnecessary. Benjamini-Hochberg correction is applied across all tested regions. `--fdr 0.05` sets the significant gain/loss cutoff. At least two replicates are required in every group; merged or undersampled inputs receive descriptive effects without differential p-values or FDR.

## Outputs

The output directory contains `differential_peaks.tsv` and/or `differential_clusters.tsv`, directional BED files, and `chip_comparison_manifest.json`. Tables include Stage 1 support, `region_origin`, every replicate score in all four groups, group means, within-condition mean enrichments, the interaction difference, p-value, differential FDR and status. Directional BED files append `region_origin` after differential FDR. The comparison manifest reports counts by origin for each feature level.

[Back to the command reference](../COMMAND_REFERENCE.md)
