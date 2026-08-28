# Glossary

## Breakpoint peak

A retained negative-score region from PNS, or a region called from the sign-inverted WPS signal. It marks a location compatible with fragment termination or exposed DNA.

## Chromosome-wise run

A run in which selected contigs are processed independently, often in parallel, before compatible outputs are combined into genome-level files.

## Coverage

The number of accepted fragments covering each genomic base. In cutn-suite, coverage—not the PNS score—is normalized to a non-zero mean of 100 for Stage 1 measurements.

## Dyad

The central position of a paired-end fragment, using the configured convention for even-length fragments.

## Fragment

A zero-based, half-open interval [start,end) inferred from a properly paired alignment or supplied directly as an interval. Its length is end-start.

## PNS

The probabilistic nucleosome score. PNS places a symmetric, length-adaptive inverted-cosine kernel across each accepted fragment’s scoring support. Each complete fragment contributes positive mass 100 and negative mass -100, so the positive distribution is represented in percent, the signed contribution sums to zero, and total absolute mass is 200. The native signed PNS BigWig is not score-scaled.

## posPNS

The non-negative PNS reference track made by shifting each signed PNS kernel upward until its minimum is zero. It retains the waveform and is not renormalized after shifting.

## PNS support

The genomic interval used by one fragment’s PNS kernel. For fragment length L and protected-DNA mode m, its width is W(L,m)=m+|L-m|. Short fragments extend around their dyad; fragments at or above the mode use their observed interval.

## Flank spacing

The [`flank-spacing`](commands/flank-spacing.md) analysis compares nucleosome spacing immediately upstream and downstream of reference sites, optionally by category.

## Protected-DNA mode

The dominant accepted fragment length used to define PNS geometry. pns and cutn-suite can estimate it from seeded samples and bootstrap stability, or an integer --mode can set it explicitly.

## Mean centring

Subtracting a vector mean from every value. PNS does not need a separate mean-centring step because its signed kernel is constructed with positive mass 100 and negative mass -100.

## Native score

A signal or interval score as produced by the scoring calculation, before a user-requested display multiplier or reference-mean normalization. PNS BigWigs and PNS peak scores remain native by default.

## Probability represented in percent

The PNS convention in which the positive part of each complete fragment distribution sums to 100 rather than 1. The value describes a percent-scale distribution; it does not imply that a genome-wide BigWig is itself a probability density.

## Window protection score (WPS)

A signal that rewards fragments enclosing a fixed protection window and penalizes fragments whose endpoints fall inside that window. NucleoSuite supports raw, smoothed, baseline-adjusted, and WPS peak-calling signals.

## Mean-scale

The mean-scale command divides a signal or interval score by a finite, non-zero reference mean and applies a requested multiplier. It is available for explicit downstream normalization; cutn-suite does not use it on PNS score BigWigs.

## Stage 1 and Stage 2

In cutn-suite, Stage 1 performs PNS discovery, coverage measurement, replicate gating, and cluster formation for one condition. Stage 2 compares two completed conditions using saved coverage, cluster overlap, interaction statistics, and matched aggregates.

## Seed and member gates

The S gate decides which candidate peaks seed clusters. The G gate decides which neighbouring candidates can extend a cluster. Their modes, p-value threshold, bridge gap, and cluster size are independently configurable.
