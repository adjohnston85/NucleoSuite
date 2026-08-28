# `nucleosuite pns`

## What this command does

`pns` converts paired-end fragments into a probabilistic nucleosome score (PNS) and calls positive nucleosome regions and negative breakpoint peaks. Its length-adaptive sinusoidal kernels carry +100 positive and −100 negative mass, with positional probability represented in percent.

## Why use it

Use PNS to map recurrent protection and cleavage while allowing fragments near the selected protected-DNA mode to contribute more localized positional evidence. Each fragment has equal total mass; distance from the mode broadens its wave and lowers the per-base amplitude.

## How it works

Before constructing a score, `pns` resolves the protected-DNA mode. The default `--mode auto` samples fragments from seeded randomly ordered genomic blocks and selects the most frequent integer length from the unsmoothed histogram. The automatic mode-search interval is 137–197 bp by default. After the mode is resolved, each omitted scoring bound is derived from the mode ± `--frag-mode-padding` (default 30 bp). Use an integer such as `--mode 167` when a fixed geometry is required.

### PNS geometry and percent mass

For fragment length $L$ and protected-DNA mode $m$, the scoring width is

```math
W(L,m)=m+|L-m|.
```

For $L<m$, the wave extends $m-L$ bases beyond each fragment end. For $L\ge m$, it spans the fragment. An inverted cosine is sampled over that support and its positive and negative parts are separately normalized to +100 and −100. Every complete kernel sums to zero and has total absolute mass 200. Probability mass is expressed in percent; the peak height changes with width.

The `posPNS` reference is the whole signed kernel shifted upward by its minimum for each fragment. It retains the complete waveform and is not renormalized to mass 100. Genomic PNS and `posPNS` tracks are sums over fragments and retain their native values. Summed scores are not bounded by 100.

See [PNS in Algorithms](../ALGORITHMS.md#probabilistic-nucleosome-scoring) for the equations, discrete symmetry, worked geometry, explanatory plots, accumulation and interpretation.

## Basic usage

Estimate the mode and score accepted fragments within mode ±30 bp:

```bash
nucleosuite pns \
  --bam sample.bam \
  --contigs chr1-22,chrX \
  --cores 8 \
  --out-prefix sample
```

Use an explicit mode and scoring range:

```bash
nucleosuite pns \
  --bam sample.bam \
  --mode-length 152 \
  --frag-lower 122 \
  --frag-upper 182 \
  --out-prefix sample
```

For fragment BED/bigBed input, use `--fragments` and provide chromosome sizes through `--chrom-sizes` or an indexed FASTA. A FASTA is required for sequence-dependent operations, such as dinucleotide profiles or WW/SS classification, but not for PNS scoring alone.

## Defaults

`pns` estimates the modal protected-DNA length automatically and uses **PNS** as the scoring method. The default mode-search interval is **137–197 bp**. Once the mode is known, the scoring fragment interval defaults to **mode ±30 bp**. For example, an estimated mode of 165 bp gives a scoring range of **135–195 bp**. Change the distance with `--frag-mode-padding`; `--frag-lower` and `--frag-upper` override their corresponding automatic bounds independently. `--mode 167` bypasses mode estimation. 


### Controlling the mode-centred fragment range

The automatic scoring bounds are resolved only after the mode is known:

```math
L_{lower}=\max(1,m-p),
\qquad
L_{upper}=m+p,
```

where `m` is the resolved mode and `p` is `--frag-mode-padding` (default 30 bp). For example:

```bash
nucleosuite pns \
  --bam sample.bam \
  --frag-mode-padding 25 \
  --out-prefix sample
```

If the estimated mode is 165 bp, this scores fragments from 140–190 bp. An explicit `--frag-lower` or `--frag-upper` replaces only that side of the automatic interval. `--mode-search-lower` and `--mode-search-upper` control which fragment lengths are considered when estimating an automatic mode; they do not fix the final scoring interval.

Raw PNS signal is used for peak calling by default. Savitzky–Golay smoothing is available but disabled (`--smooth-window 0`). Mode-histogram smoothing is also disabled by default; `--mode-histogram-smoothing binomial` is an independent option that affects only mode estimation.

The resolved mode is printed during execution and written to `*_fragment_mode_estimation.tsv`. Output prefixes include the scoring method, mode, fragment range, and smoothing parameters so analyses using different geometries do not silently overwrite one another.

## Peak calling

PNS uses positive-region segmentation. Positive score regions are retained when they meet `--min-region-length`; `--max-neg-run` controls how many consecutive zero-or-negative bases may be bridged inside a region. Breakpoint calls apply the same segmentation to the sign-inverted score.

Use `--no-peak-calling` when only score and auxiliary tracks are required. See [shared nucleosome-score peak calling](../ALGORITHMS.md#pns-peak-calling) for the exact shared peak-calling definition.

## Optional coverage filtering of nucleosome peaks

`--peak-coverage-threshold` filters nucleosome peaks using coverage at the representative position written to BED column 7. A peak is retained when coverage is greater than or equal to the selected threshold.

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --peak-coverage-threshold 2 \
  --out-prefix sample
```

Breakpoint peaks are not filtered by this option. Existing callsets can be filtered with [`filter-peaks`](filter-peaks.md).

## Peak scores and bigBed

Text BED peak scores are floating-point values written to six decimal places after `--peak-score-scale` (default 1). `--bigbed-score-scale` also defaults to 1 and controls only conversion to the integer bigBed score field, which is rounded and clamped to 0–1000. Neither option scales BigWig values.

## Outputs

The score tracks are `pns`, `posPNS`, and optionally `pns_smoothed`. Coverage, dyads and combined/left/right fragment-end tracks are also written by default and can be selected with `--other-tracks`.

A default fixed-mode PNS run therefore produces names such as:

```text
sample_methodpns_mode167_lower137_upper197_smooth0x2_pns.bw
sample_methodpns_mode167_lower137_upper197_smooth0x2_posPNS.bw
sample_methodpns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed
```

BED column 7 stores the midpoint of the retained score region. Downstream `distances` can read that representative position with `--position-column 7`.

## What to do next

Measure spacing between called nucleosome regions:

```bash
nucleosuite distances sample_methodpns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed \
  --position-column 7 \
  --output-prefix sample_spacing
```

Or aggregate PNS around reference sites:

```bash
nucleosuite aggregate \
  --bigwig sample_methodpns_mode167_lower137_upper197_smooth0x2_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --output-prefix sample_ctcf
```

[Back to the command reference](../COMMAND_REFERENCE.md)
