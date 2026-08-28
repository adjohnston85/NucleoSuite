# `nucleosuite pns`

`pns` converts accepted paired-end fragments into a length-adaptive probabilistic nucleosome score and calls positive nucleosome regions and negative breakpoint peaks. It accepts BAM, CRAM, BED, BED.gz, and bigBed fragment inputs through the shared fragment-processing options.

## What this command does

`pns` turns accepted paired-end fragments into a length-adaptive nucleosome-positioning signal and its positive-region/breakpoint calls.

## Why use it

Use `pns` when a single analysis needs the core PNS signal, native score tracks, and directly associated nucleosome or breakpoint features. Use [`tracks`](tracks.md) when the same fragment pass should also produce WPS, coverage, dyads, ends, or sequence-derived outputs.

## What it calculates

For a protected-DNA mode $m$ and fragment length $L$, PNS uses a symmetric support of

```math
W(L,m)=m+|L-m|.
```

An inverted cosine is sampled across that support. Positive and negative samples are normalized separately, giving each complete accepted fragment a positive mass of `+100` and a negative mass of `-100`. The positive distribution is therefore represented in percent. The signed contribution sums to zero and has total absolute mass 200. Shorter or longer fragments use wider support, while fragments at the mode use the narrowest support.

`posPNS` is the non-negative reference version of the same waveform: the signed kernel is shifted upward until its minimum is zero, without clipping or an additional normalization. The signed `pns` track is used for peak calling and downstream positioning analyses.

## Basic invocation

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1 chr2 chr3 chr4 \
  --cores 4 \
  --out-prefix sample
```

The command writes a mode-estimation report for automatic mode runs. To use a fixed protected-DNA mode and explicit scoring bounds:

```bash
nucleosuite pns \
  --fragments sample.fragments.bed.gz \
  --fasta genome.fa \
  --mode 167 \
  --frag-lower 137 \
  --frag-upper 197 \
  --out-prefix sample_pns
```

With `--mode auto` (the default), accepted fragments are sampled in seeded genomic-block order and the most frequent integer length in the mode-search interval is selected. The default search interval is 137–197 bp. If scoring bounds are omitted, they are derived from the resolved mode plus or minus `--frag-mode-padding` (default 30 bp). Supplying an integer mode bypasses estimation.

## Important options

| Option | Purpose |
|---|---|
| `--score-tracks pns posPNS pns_smoothed` | Select the signed score, non-negative reference, and optional smoothed score tracks. |
| `--smooth-window N` | Enable Savitzky–Golay score smoothing with an odd window of at least 3. |
| `--peak-calling/--no-peak-calling` | Write or suppress nucleosome and breakpoint peak calls. |
| `--min-region-length N` | Minimum positive-region length in bases. |
| `--max-neg-run N` | Allow up to N consecutive zero-or-negative bases within a positive region. |
| `--peak-coverage-threshold N` | Retain nucleosome peaks only when coverage at the representative position meets N. |
| `--score-format` | Write BigWig, WIG.GZ, both, or no PNS signal files. |
| `--other-tracks` | Select coverage, dyad, and fragment-end tracks produced in the same pass. |
| `--split-ww-types` | Write all-fragment and sequence-type-specific outputs; requires `--fasta`. |

## Peak scores and BigWig values

Text BED scores are floating-point PNS values after `--peak-score-scale`. BigBed scores are converted to the integer 0–1000 field using `--bigbed-score-scale`; the PNS default is `1`, so native PNS values are not rescaled. Set another multiplier only when a downstream integer-score convention requires it.

PNS BigWigs also retain their native values. NucleoSuite does not apply a post-generation score normalization. Coverage can be normalized independently by [`mean-scale`](mean-scale.md) where a coverage-based comparison requires it.

## Outputs

For output prefix `sample`, the principal files are:

- `sample_pns.bw` — signed PNS score;
- `sample_posPNS.bw` — non-negative PNS reference;
- `sample_pns_smoothed.bw` — optional smoothed score;
- `sample_nucleosome_regions.bed` or `.bb` — positive PNS regions;
- `sample_breakpoint_peaks.bed` or `.bb` — negative PNS regions after sign inversion; and
- `sample_fragment_mode_estimation.tsv` — automatic mode-estimation details, when applicable.

Coverage, dyad, fragment-end, dinucleotide, and WW/SS files are written when requested. On multiple contigs, per-contig outputs can be combined with [`combine`](../commands/combine.md). For several signal families or fragment ranges in one pass, use [`tracks`](tracks.md).

## Related analyses

Use [`call-peaks`](call-peaks.md) to call features from an existing PNS BigWig, [`distances`](distances.md) for peak spacing, [`aggregate`](aggregate.md) for signal around reference sites, and [`positive-runs`](positive-runs.md) for contiguous positive-signal intervals.

[Back to the command reference](../COMMAND_REFERENCE.md)
