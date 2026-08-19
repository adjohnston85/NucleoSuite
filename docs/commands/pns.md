# `nucleosuite pns`

## What this command does

`pns` converts paired-end fragments into a nucleosome-position signal and calls positive nucleosome regions and negative breakpoint peaks. The default scoring method is probabilistic nucleosome score (PNS). `--scoring-method bns` uses boxcar nucleosome score (BNS), while `--scoring-method tns` uses triangular nucleosome score (TNS).

## Why use it

Use PNS when you want endpoint-derived dyad support from nucleosome-sized fragments. Use BNS when you want a balanced central boxcar contribution. Use TNS when you want one centred triangular contribution that rises from the scoring boundaries to the fragment centre. All three methods use the same fragment filtering, genomic accumulation, optional smoothing, and peak caller, so their signals and peak calls can be compared directly.

## How it works

### PNS

Each fragment end contributes a triangular dyad-support distribution with mass 0.5. The two distributions are added, and their maxima coincide when the fragment length equals the selected protected-DNA mode. NucleoSuite subtracts the combined distribution's mean so that each complete fragment contributes values that sum to zero, then adds the fragment contributions across the genome.

### BNS

BNS uses the same scoring support as PNS for every accepted fragment length. Within that support, it builds a symmetric central boxcar with total mass 1 and zero outer flanks. The boxcar is mean-centred so that each fragment contributes a positive central region and negative flanks whose total contribution sums to zero. Discrete half-weight or zero transition positions are used where needed to keep odd and non-divisible support lengths symmetric.

BNS kernels are precomputed for every accepted fragment length in the same way that PNS kernels are precomputed, then reused for all fragments of that length.

### TNS

TNS uses the same fragment-length-dependent scoring support as PNS and BNS. It places one symmetric triangle across the entire support. The raw triangle begins at zero, rises toward the fragment centre, and returns to zero at the opposite boundary. Odd support lengths have one central maximum, while even support lengths have a two-base central plateau. The raw triangle is normalized to total mass 1 and then mean-centred so that each fragment contributes values that sum to zero.

For fragments shorter than the selected mode, the support expands in the same way as PNS and BNS. With a 167 bp mode, for example, a 137 bp fragment uses a 197 bp triangle because it is 30 bp shorter than the mode. A 197 bp fragment also uses a 197 bp triangle because it is already longer than the mode.

TNS kernels are precomputed for every accepted fragment length and reused for all fragments of that length.

The complete PNS, BNS and TNS constructions are described in [Nucleosome scoring](../ALGORITHMS.md#probabilistic-nucleosome-scoring).

## Typical use

Default PNS:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1-22,chrX \
  --cores 8 \
  --out-prefix sample_pns
```

BNS with the same fragment selection and peak caller:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1-22,chrX \
  --cores 8 \
  --scoring-method bns \
  --out-prefix sample
```

TNS with the same fragment selection and peak caller:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1-22,chrX \
  --cores 8 \
  --scoring-method tns \
  --out-prefix sample
```

## Defaults

Standalone `pns` uses fragment lengths **137–197 bp** with modal protected-DNA length **167 bp** by default. `mnase-suite` uses an MNase-specific PNS configuration of **120–180 bp** with mode **147 bp**.

`--scoring-method pns` is the default. Select `--scoring-method bns` for BNS or `--scoring-method tns` for TNS.

Raw PNS, BNS or TNS is the default peak-calling signal. Savitzky–Golay smoothing is available but disabled by default (`--smooth-window 0`).

## Peak calling

PNS, BNS and TNS use the same peak caller. The caller retains positive regions that meet `--min-region-length`. `--max-neg-run` sets the permitted number of consecutive zero-or-negative bases within a region; its default is 0.

Breakpoint calls apply the same region logic to the sign-inverted scoring signal, so negative regions become positive for the caller.

See [PNS peak calling](../ALGORITHMS.md#pns-peak-calling) for the exact definition.

## Optional coverage filtering of nucleosome peaks

`--peak-coverage-threshold` filters nucleosome peaks using the fragment coverage already calculated during the PNS, BNS or TNS run. Coverage is evaluated at the representative peak position written to BED column 7. A peak is retained when its coverage is greater than or equal to the selected threshold.

For example:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --peak-coverage-threshold 2 \
  --out-prefix sample_pns
```

Only nucleosome peaks with coverage of at least 2 at BED column 7 are written. Breakpoint peaks are not filtered by this option. The filter is off by default.

To apply the same type of filter to an existing BED and coverage BigWig, use [`filter-coverage`](filter-coverage.md).

## Peak scores and bigBed

Text BED peak scores remain floating-point values written to six decimal places after `--peak-score-scale`.

When peaks are converted to bigBed, `--bigbed-score-scale` controls the multiplier used before conversion to the required integer BED5 score. The default is **1000**, followed by rounding and clamping to the BED score range 0–1000.

This setting affects the bigBed score field, not the floating-point text BED score.

## What it writes

For PNS, the selected score tracks are `pns`, `posPNS`, and optionally `pns_smoothed`. For BNS they are `bns`, `posBNS`, and optionally `bns_smoothed`. For TNS they are `tns`, `posTNS`, and optionally `tns_smoothed`.

BNS and TNS output prefixes include the scoring method, mode, and accepted fragment range so they are distinct from PNS outputs. For example:

```text
sample_methodbns_mode167_lower137_upper197_smooth0x2_bns.bw
sample_methodbns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed
sample_methodtns_mode167_lower137_upper197_smooth0x2_tns.bw
sample_methodtns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed
```

BED column 7 stores the midpoint of the retained score region. Downstream `distances` can read that representative position with `--position-column 7`.

## What to do next

Measure spacing between called nucleosome regions:

```bash
nucleosuite distances sample_methodpns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed \
  --position-column 7 \
  --output-prefix sample_spacing
```

Or compare signal around bundled CTCF sites:

```bash
nucleosuite aggregate \
  --bigwig sample_methodpns_mode167_lower137_upper197_smooth0x2_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --output-prefix sample_ctcf
```

[Back to the command reference](../COMMAND_REFERENCE.md)
