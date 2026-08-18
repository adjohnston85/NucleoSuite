# `nucleosuite dcc`

## What this command does

`dcc` calculates distance cross-correlation (DCC) between two genomic signals to identify the offset of signal B relative to signal A.

## Why use it

Use DCC to measure upstream or downstream offsets between two signal collections, such as dyads from different fragment-length classes or dyads and fragment ends.

## How to read lag

DCC uses the convention:

```text
lag = position of B - position of A
```

Under this convention:

- a peak at **0 bp** means A and B tend to occur at the same position;
- a peak at **+10 bp** means B tends to occur 10 bp downstream of A; and
- a peak at **−10 bp** means B tends to occur 10 bp upstream of A.

For strand-oriented feature analyses, minus-strand regions are reversed so positive lag still means downstream in the feature's orientation.

The pair-product definition and a +10 bp example are shown in [Distance cross-correlation](../ALGORITHMS.md#distance-cross-correlation).

## Choose the input mode

### BigWig mode

Use `dcc bigwig` when A and B already exist as genomic signal tracks:

```bash
nucleosuite dcc bigwig \
  --bigwig-a short_fragment_dyad.bw \
  --bigwig-b long_fragment_dyad.bw \
  --chrom-sizes sample.bam \
  --signed-lags \
  --dmax 500 \
  --out-prefix short_vs_long
```

All A BigWigs are added base-by-base within each region to form one A signal; all B BigWigs are combined in the same way before DCC is calculated.

Missing or non-finite BigWig values are treated as zero signal but are still included when counting the genomic positions available for comparison. Blacklisted bases are excluded entirely from both the signal calculation and the number of available comparisons.

### Fragment/BAM mode

Use `dcc bam` when A and B should be constructed directly from fragments. Each side can select `dyad`, `left_end`, or `right_end` positions and its own fragment-length range.

```bash
nucleosuite dcc bam \
  --bam-a sample.bam --min-length-a 145 --max-length-a 147 --position-a dyad \
  --bam-b sample.bam --min-length-b 145 --max-length-b 147 --position-b right_end \
  --chrom-sizes sample.bam \
  --signed-lags \
  --out-prefix dyad_vs_right_end
```

## Analyse only selected regions

Both input modes can use a region BED. For example, restrict a BigWig comparison to the bundled GM12878 ChromHMM intervals:

```bash
nucleosuite dcc bigwig \
  --bigwig-a signal_A.bw \
  --bigwig-b signal_B.bw \
  --regions-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --state-column 4 \
  --signed-lags \
  --out-prefix A_vs_B_by_state
```

## Signed versus absolute distance

Use `--signed-lags` when upstream/downstream direction matters.

Without `--signed-lags`, NucleoSuite collapses `+d` and `−d` into one absolute distance. It combines the raw products and their opportunity counts **before** calculating the normalized value. This prevents unequal numbers of valid +d and −d comparisons from being weighted incorrectly.

## Normalization

The default DCC value divides the raw A×B product sum by the number of valid A/B position pairs at that lag. `--no-normalize-dcc` selects the raw product sum.

`--normalize-by-signal-totals` applies an additional signal-total normalization when required. Percentage and per-million signal-pair columns are derived after the principal profile is formed.

## What it writes

Each analysed state or category receives:

- a DCC TSV;
- a maximum-lag or maximum-distance summary;
- a run summary; and
- a DCC plot.

All calculated values remain in the TSV even when very short distances are omitted from the standard plot for readability.

## Plot customization

Use `--plot-label-points peaks` when you want retained DCC local maxima annotated. The shared plotting options are described in [Plot customization](../PLOTTING.md).

[Back to the command reference](../COMMAND_REFERENCE.md)
