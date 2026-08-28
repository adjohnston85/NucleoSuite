# `nucleosuite nuc-score`

## What this command does

`nuc-score` converts paired-end fragments into a nucleosome-position signal and calls positive nucleosome regions and negative breakpoint peaks. The command uses **Sinusoidal Nucleosome Scoring (SNS)** by default. The historical Probabilistic Nucleosome Score remains available with `--scoring-method pns`; `--scoring-method bns` selects Boxcar Nucleosome Scoring and `--scoring-method tns` selects Triangular Nucleosome Scoring.

## Why use it

Use SNS when you want each accepted fragment to contribute one smooth, symmetric nucleosome-centred wave whose spatial precision is greatest at the selected protected-DNA mode and broadens as fragment length departs from that mode. Use PNS for endpoint-derived dyad probability support, BNS for a balanced central boxcar, or TNS for a centred triangle. All four scoring methods use the same fragment filtering, genomic accumulation, optional smoothing, and PNS peak caller, so their signals and calls can be compared directly.

## How it works

Before constructing a score, `nuc-score` resolves the protected-DNA mode. The default `--mode auto` samples fragments from seeded randomly ordered genomic blocks and selects the most frequent integer length from the unsmoothed histogram. The automatic mode-search interval is 137–197 bp by default. After the mode is resolved, each omitted scoring bound is derived from the mode ± `--frag-mode-padding` (default 30 bp). Use an integer such as `--mode 167` when a fixed geometry is required.

### SNS — default

For a fragment of length `L` and selected mode `m`, SNS uses a support width

```math
W(L,m)=m+|L-m|.
```

Therefore a mode-length fragment has the narrowest wave. Longer fragments use their own length as the support width, while shorter fragments extend symmetrically beyond both observed fragment boundaries. With `m=167`, a 120 bp fragment uses a 214 bp wave, a 167 bp fragment uses a 167 bp wave, and a 180 bp fragment uses a 180 bp wave.

Across the discrete integer genomic positions of that support, SNS samples one inverted cosine cycle: it is lowest at both support boundaries and highest at the fragment centre. Positive and negative values are normalized separately so that every complete fragment contributes exactly **+50 positive mass** and **−50 negative mass**. Thus the signed contribution sums to zero and its total absolute mass is 100 regardless of fragment length.

Because `W(L,m)` always has the same odd/even parity as `L`, odd fragments have one central genomic bin and even fragments have two equal central bins around the half-base midpoint. This keeps the discrete wave exactly symmetric without shifting even-length fragments by half a base.

`posSNS` is the complete SNS waveform shifted upward so that its minimum is zero. It is not clipped or renormalized after the shift; it is an auxiliary non-negative reference track and is not the signed signal used for peak calling.

### PNS

Each fragment end contributes a triangular dyad-support distribution with mass 0.5. The two distributions are added, and their maxima coincide when the fragment length equals the selected protected-DNA mode. NucleoSuite subtracts the combined distribution's mean so that each complete fragment contributes values that sum to zero.

### BNS

BNS uses the same fragment-length-dependent support as SNS and PNS. It builds a symmetric central boxcar with total mass 1 and zero outer flanks, then mean-centres that distribution to produce a positive central contribution and negative flanks.

### TNS

TNS uses the same support geometry and places one symmetric triangle across the complete support. The raw triangle is normalized to total mass 1 and mean-centred so that the signed contribution sums to zero.

The full mathematical definitions are in [Algorithms](../ALGORITHMS.md), including [SNS](../ALGORITHMS.md#sinusoidal-nucleosome-scoring), [PNS](../ALGORITHMS.md#probabilistic-nucleosome-scoring), [BNS](../ALGORITHMS.md#boxcar-nucleosome-scoring), and [TNS](../ALGORITHMS.md#triangular-nucleosome-scoring).

## Basic usage

Default SNS scoring:

```bash
nucleosuite nuc-score \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1-22,chrX \
  --cores 8 \
  --out-prefix sample
```

Explicit PNS scoring with the same fragment selection and peak caller:

```bash
nucleosuite nuc-score \
  --bam sample.bam \
  --fasta genome.fa \
  --scoring-method pns \
  --out-prefix sample
```

BNS or TNS can be selected in the same way:

```bash
nucleosuite nuc-score --bam sample.bam --scoring-method bns --out-prefix sample
nucleosuite nuc-score --bam sample.bam --scoring-method tns --out-prefix sample
```

## Defaults

`nuc-score` estimates the modal protected-DNA length automatically and uses **SNS** as the scoring method. The default mode-search interval is **137–197 bp**. Once the mode is known, the scoring fragment interval defaults to **mode ±30 bp**. For example, an estimated mode of 165 bp gives a scoring range of **135–195 bp**. Change the distance with `--frag-mode-padding`; `--frag-lower` and `--frag-upper` override their corresponding automatic bounds independently. `--mode 167` bypasses mode estimation. `--scoring-method pns`, `bns`, or `tns` selects one of the alternative kernels.


### Controlling the mode-centred fragment range

The automatic scoring bounds are resolved only after the mode is known:

```math
L_{lower}=\max(1,m-p),
\qquad
L_{upper}=m+p,
```

where `m` is the resolved mode and `p` is `--frag-mode-padding` (default 30 bp). For example:

```bash
nucleosuite nuc-score \
  --bam sample.bam \
  --frag-mode-padding 25 \
  --out-prefix sample
```

If the estimated mode is 165 bp, this scores fragments from 140–190 bp. An explicit `--frag-lower` or `--frag-upper` replaces only that side of the automatic interval. `--mode-search-lower` and `--mode-search-upper` control which fragment lengths are considered when estimating an automatic mode; they do not fix the final scoring interval.

Raw SNS/PNS/BNS/TNS signal is used for peak calling by default. Savitzky–Golay smoothing is available but disabled (`--smooth-window 0`). Mode-histogram smoothing is also disabled by default; `--mode-histogram-smoothing binomial` is an independent option that affects only mode estimation.

The resolved mode is printed during execution and written to `*_fragment_mode_estimation.tsv`. Output prefixes include the scoring method, mode, fragment range, and smoothing parameters so analyses using different kernels or geometries do not silently overwrite one another.

## Peak calling

SNS, PNS, BNS, and TNS all use the same PNS peak caller. Positive score regions are retained when they meet `--min-region-length`; `--max-neg-run` controls how many consecutive zero-or-negative bases may be bridged inside a region. Breakpoint calls apply the same segmentation to the sign-inverted score.

Use `--no-peak-calling` when only score and auxiliary tracks are required. See [shared nucleosome-score peak calling](../ALGORITHMS.md#pns-peak-calling) for the exact shared peak-calling definition.

## Optional coverage filtering of nucleosome peaks

`--peak-coverage-threshold` filters nucleosome peaks using coverage at the representative position written to BED column 7. A peak is retained when coverage is greater than or equal to the selected threshold.

```bash
nucleosuite nuc-score \
  --bam sample.bam \
  --fasta genome.fa \
  --peak-coverage-threshold 2 \
  --out-prefix sample
```

Breakpoint peaks are not filtered by this option. Existing callsets can be filtered with [`filter-peaks`](filter-peaks.md).

## Peak scores and bigBed

Text BED peak scores remain floating-point values written to six decimal places after `--peak-score-scale`. `--bigbed-score-scale` controls conversion to the required integer BED5 score. SNS defaults to **1**, so native SNS peak scores are not rescaled before integer rounding/clamping. PNS, BNS and TNS default to **1000** because their native peak scores are fractional. An explicit value overrides the method-aware default.

## Outputs

For SNS, the score tracks are `sns`, `posSNS`, and optionally `sns_smoothed`. PNS uses `pns`, `posPNS`, and optionally `pns_smoothed`; BNS uses `bns`, `posBNS`, and optionally `bns_smoothed`; TNS uses `tns`, `posTNS`, and optionally `tns_smoothed`.

A default fixed-mode SNS run therefore produces names such as:

```text
sample_methodsns_mode167_lower137_upper197_smooth0x2_sns.bw
sample_methodsns_mode167_lower137_upper197_smooth0x2_posSNS.bw
sample_methodsns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed
```

BED column 7 stores the midpoint of the retained score region. Downstream `distances` can read that representative position with `--position-column 7`.

## What to do next

Measure spacing between called nucleosome regions:

```bash
nucleosuite distances sample_methodsns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed \
  --position-column 7 \
  --output-prefix sample_spacing
```

Or aggregate SNS around reference sites:

```bash
nucleosuite aggregate \
  --bigwig sample_methodsns_mode167_lower137_upper197_smooth0x2_sns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --output-prefix sample_ctcf
```

[Back to the command reference](../COMMAND_REFERENCE.md)
