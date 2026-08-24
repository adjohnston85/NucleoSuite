# `nucleosuite chip-suite`

## What this command does

`chip-suite` compares a ChIP-seq, CUT&RUN or CUT&Tag target with its matched control using nucleosome-position score tracks. It generates target and control scores independently, scales each centred score track by the finite non-zero mean of its corresponding positive score track, calls peaks independently, and uses control peaks as empirical decoys for target peak and cluster FDR.

TNS is the default score. BNS and PNS are available with `--scoring-method bns` or `--scoring-method pns`.

## Why use it

Use this workflow when target and control assays should be processed with identical fragment filters and score definitions, normalized without subtracting one track from the other, and evaluated with a control-derived empirical peak and cluster FDR.

## Typical run

```bash
nucleosuite chip-suite \
  --target-bam H3K4me3.bam \
  --control-bam H3_control.bam \
  --outdir H3K4me3_chip_suite \
  --sample-name H3K4me3 \
  --cores 8
```

The default accepted fragment range is 120–500 bp.

## Automatic fragment mode

`--mode auto` is the default. Target and control modes are estimated independently by visiting indexed genomic blocks from the selected analysis contigs in seeded random order, accumulating accepted fragment lengths, and bootstrapping a lightly smoothed 120–250-bp length histogram. Sampling stops when the mode is stable across three checkpoints and the bootstrap 95% interval is no wider than 4 bp, or after the maximum sample size is reached.

The default `--mode-strategy pooled` gives target and control histograms equal weight and uses one pooled analysis mode for both score tracks. Independent estimates remain in the QC report. Other strategies are `separate`, `target` and `control`.

Automatic estimation can be bypassed completely:

```bash
nucleosuite chip-suite \
  --target-bam target.bam \
  --control-bam control.bam \
  --outdir chip_results \
  --mode 167
```

An integer `--mode` is used exactly for both samples and the report records `mode_source=explicit`.

## Score-matched normalization

For default TNS scoring:

```math
T_{scaled}(x)=\frac{TNS_T(x)}{mean(posTNS_T)}
```

```math
C_{scaled}(x)=\frac{TNS_C(x)}{mean(posTNS_C)}
```

BNS uses mean posBNS and PNS uses mean posPNS. The mean is calculated across finite, non-zero bases in the corresponding positive-score BigWig. The control is not subtracted from the target. Peaks are called independently on the two normalized tracks, preserving their native peak shapes.

## Peak competition and FDR

Nearby target and control summits are matched one-to-one within half the analysis mode by default. The higher normalized summit score wins; ties conservatively belong to the control. Unmatched peaks remain target or control wins. Control winners act as decoys when cumulative target FDR is calculated. Target-winning peaks receive monotonic empirical q-values and are retained at `--peak-fdr 0.05` by default.

## Clusters

Significant target-winning peaks are clustered within each contig. A cluster ends after five consecutive nonsignificant candidate peaks or when significant summits are separated by more than 1000 bp. At least two significant peaks are required. Control-winning peaks are clustered with the same rules, and their cluster-score distribution supplies the empirical cluster FDR. The default cluster cutoff is `--cluster-fdr 0.05`.

## Output layout

- `00_setup/`: mode estimates and normalization QC;
- `01_score_tracks/`: raw target and control score/positive-score tracks;
- `02_mean_scaled_tracks/`: score divided by matching positive-score mean;
- `03_peak_calls/`: independent target and control peak calls;
- `04_peak_fdr/`: all target peaks with appended empirical FDR, significant peaks, cluster tables and significant clusters.

[Back to the command reference](../COMMAND_REFERENCE.md)
