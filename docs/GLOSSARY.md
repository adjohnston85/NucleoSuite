# Glossary

## Aggregate profile

The mean genomic signal obtained after centring multiple regions on a common reference position. For strand-aware analyses, negative-strand regions are reversed before the mean is calculated.

## BigWig

An indexed binary format for continuous genomic signal. BigWig files support efficient regional access and display in genome browsers.

## Breakpoint

A genomic feature associated with frequent fragment termination or exposed DNA. PNS, BNS and TNS breakpoint calls are retained negative-score regions, while WPS breakpoint calls are obtained by applying WPS segmentation to the sign-inverted calling signal.

## BNS

**Boxcar nucleosome score.** BNS uses the same fragment-length-dependent scoring support as PNS but places a symmetric unit-mass central boxcar within that support. Mean centring produces a positive central contribution and negative flanks with zero total contribution for every fragment. Half-weight or zero transition positions preserve symmetry when the discrete support cannot be divided into four equal integer blocks.

## Chromosome and contig

A contig is a named reference sequence in a genome assembly. Human chromosomes such as `chr1` and `chrX` are contigs, as are unplaced scaffolds. NucleoSuite can analyse contigs separately across multiple processor cores and combine compatible outputs afterward.

For BAM-based workflows, generated contig names follow the BAM headers. Equivalent `chr`/non-`chr` spellings are matched conservatively during analysis; when multiple BAMs mix both spellings, the canonical output uses the `chr` form.

## DAC

**Distance autocorrelation coefficient.** At distance $d$, DAC sums the products $S(x)S(x+d)$ for eligible pairs of positions within the same selected region. It is commonly applied to nucleosome dyad signal to detect regular nucleosome spacing. For binary dyad tracks, the raw DAC is the number of dyad pairs separated by $d$; for weighted tracks, each pair contributes the product of its values.

For example, if nucleosome dyads occur at positions 1,000, 1,185, 1,370 and 1,555 bp, the dyads are separated by 185 bp. The DAC profile will show peaks at 185, 370 and 555 bp, corresponding to one, two and three nucleosome spacings. The spacing between successive DAC peaks can be used to estimate nucleosome repeat length.

## DCC

**Distance cross-correlation coefficient.** At signed lag $\ell$, DCC sums the products $A(x)B(x+\ell)$ for eligible position pairs. NucleoSuite defines lag as `position_B - position_A`.

A DCC peak at 0 bp indicates that the signals tend to occur at the same positions. A peak at `+10 bp` indicates that the second signal is most strongly associated with positions 10 bp downstream of the first. DCC can compare dyads from different fragment-length classes or compare corresponding fragment-end signals.

## Dinucleotide profile

The frequency of each adjacent two-base sequence at each position across aligned fragments. Each retained fragment contributes its sequence once. Sequence-aware profiles use only fragments whose complete extracted sequence contains canonical bases (`A`, `C`, `G` or `T`); a fragment containing any non-ACGT base is skipped in full.

## Dyad

The central position of nucleosome-bound DNA. A dyad track uses the centre of each fragment as an estimate of this position. Odd-length fragments have one central base. Even-length fragments have two central bases; dyad tracks split the signal equally between them by default. `--even-dyad left` and `--even-dyad right` assign the full count to one central base.

Dinucleotide profiles do not split even-length fragments between alternative centres.

## Flanking spacing

The distance between the nearest nucleosome centre strictly upstream of a reference site and the nearest nucleosome centre strictly downstream. In `flank-spacing`, these distances are compared across categories defined in a reference BED.

## Fragment

The genomic interval spanned by a properly paired sequencing-read pair. Fragment coordinates are represented as zero-based, half-open intervals.

## Lag

The signed positional difference `position_B - position_A` used by DCC. A positive lag places signal B downstream of signal A; a negative lag places B upstream. For strand-aware regions, minus-strand inputs are reversed so this interpretation remains feature-oriented.

## L-WPS

**Long-window protection score.** NucleoSuite's WPS implementation was written to reproduce the L-WPS algorithm used by [Snyder et al.](https://doi.org/10.1016/j.cell.2015.11.050). It uses the same default parameters: 120–180 bp fragments, a 120 bp protection window, a 1,000 bp running-median adjustment, and 21 bp second-order Savitzky–Golay smoothing.

## Mean centring

Subtraction of a vector's arithmetic mean from every value. PNS, BNS and TNS mean-centre each fragment distribution before placing it on the genome, making the values contributed by each complete fragment sum to zero.

## NRL

**Nucleosome repeat length.** The average centre-to-centre spacing between successive nucleosomes, including nucleosomal DNA and linker DNA.

NucleoSuite estimates NRL by detecting regularly spaced DAC or DCC peaks and fitting a regression of peak position against nucleosome order. The fitted slope is the estimated repeat length.

## Nucleosome order

The sequential number assigned to recurring spacing peaks. The first peak represents approximately one nucleosome spacing, the second represents two spacings, and so on.

## Opportunity normalization

Adjustment for the number of position pairs that can contribute at each separation distance. Fewer comparisons may be possible at long distances or near region boundaries. Dividing the raw DAC or DCC value by the number of available comparisons makes distances with different numbers of opportunities comparable.

## PNS

**Probabilistic nucleosome score.** PNS constructs one dyad-support triangle from each fragment end. Each triangle has discrete probability mass 0.5, giving their combined distribution mass 1 for every fragment length. Odd mode lengths produce one maximum in each triangle; even mode lengths produce two equal maxima. Subtracting the combined distribution's mean makes the values contributed by each complete fragment sum to zero.

The fragment contributions are summed at each genomic position. Positive regions support nucleosome dyads, and negative regions are compatible with fragment boundaries or cleavage. A called nucleosome region receives its maximum PNS value as its score.

## `posPNS`

The non-negative PNS distribution before mean subtraction. `posPNS` adds the two endpoint-derived triangles from every accepted fragment. Each fragment contributes total mass 1. The name refers to this distribution, not to values clipped from the PNS track.

## `posBNS`

The non-negative BNS unit-mass boxcar before mean subtraction. Each accepted fragment contributes total mass 1.

## TNS

**Triangular nucleosome score.** TNS places one symmetric unit-mass triangle across the fragment scoring support. The triangle is zero at both support boundaries, has one central maximum for odd support lengths and a two-base central plateau for even support lengths, and is mean-centred before genome-wide accumulation.

## `posTNS`

The non-negative TNS unit-mass triangle before mean subtraction. Each accepted fragment contributes total mass 1.

## Positional offset

The separation between two genomic positions, measured in base pairs. A signed offset retains direction: negative values indicate upstream displacement and positive values indicate downstream displacement. An absolute offset reports only the separation.

## Summit

The representative coordinate of a called peak. It may be read from a specified interval column or derived from the interval midpoint.

## WPS

**Window protection score.** A window is centred at each genomic position. Fragments spanning the complete window increase the score, while fragments terminating within the window decrease it. High WPS values indicate regions frequently protected within intact fragments; low values indicate frequent fragmentation within the window.

NucleoSuite can also write Savitzky–Golay-smoothed WPS (`wps_smoothed`), raw WPS minus its running-median baseline (`mWPS`), and smoothed WPS minus the running median of raw WPS (`sm_mWPS`). The default peak caller evaluates `sm_mWPS`.

## WW and SS dinucleotides

WW dinucleotides contain two weak bases (`A` or `T`): `AA`, `AT`, `TA` and `TT`. SS dinucleotides contain two strong bases (`C` or `G`): `CC`, `CG`, `GC` and `GG`.

## WW/SS classification

Classification of fragments using a centred 147 bp reference core and predefined minor- and major-groove-associated positions. WW and SS enrichment are evaluated independently, and fragments are assigned to `type1`–`type4` according to the resulting enrichment pattern.
