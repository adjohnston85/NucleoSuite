# NucleoSuite workflows

This page shows how NucleoSuite commands connect. Each diagram is followed by a command-line example.

## cfDNA: from fragments to nucleosome organization

```mermaid
flowchart TB
    A[Paired-end cfDNA BAM or fragment BED] --> B[cfdna-suite]
    B --> C[PNS mode 167 and WPS]
    B --> D[Dyads, fragment ends, and fragment lengths]
    B --> E[Dinucleotide and WW/SS analyses]
    C --> F[Peak calls]
    D --> G[DAC and DCC]
    F --> H[Peak spacing and callset comparisons]
    G --> I[NRL and positional relationships]
    B --> J[CTCF, TSS, chromatin-state, and expression analyses]
```

[`cfdna-suite`](commands/cfdna-suite.md) applies cfDNA-oriented defaults and uses one filtered fragment set throughout the full analysis tree.

```bash
nucleosuite cfdna-suite \
  --bam sample.bam \
  --fasta hg19.fa \
  --resource-set hg19-gm12878 \
  --contigs chr1-22,chrX \
  --cores 8 \
  --outdir sample_cfdna_suite
```

Use standalone commands for individual analyses or different fragment-selection settings between outputs.

### PNS followed by spacing analysis

```mermaid
flowchart LR
    A[Paired-end BAM or fragment BED] --> B[PNS]
    B --> C[PNS BigWig]
    B --> D[Nucleosome-region calls]
    D --> E[distances]
    E --> F[Nearest-neighbour spacing profile]
```

The same workflow can use BNS with `--scoring-method bns` or TNS with `--scoring-method tns`; the nucleosome peak caller is unchanged.

For example:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta hg19.fa \
  --out-prefix sample_pns

nucleosuite distances sample_methodpns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed \
  --position-column 7 \
  --min-distance 120 \
  --max-distance 250 \
  --output-prefix sample_spacing
```

## MNase-seq: from protected fragments to nucleosome organization

```mermaid
flowchart TB
    A[MNase BAM or fragment BED] --> B[mnase-suite]
    B --> C[PNS mode 147 and WPS]
    B --> D[Coverage, dyads, and fragment ends]
    B --> E[Dinucleotide and WW/SS analyses]
    C --> F[Peak calls]
    D --> G[DAC, DCC, and NRL]
    F --> H[Peak spacing and callset comparisons]
    B --> I[CTCF, TSS, chromatin-state, and expression analyses]
```

[`mnase-suite`](commands/mnase-suite.md) applies MNase-oriented defaults and uses one retained fragment population across downstream outputs.

```bash
nucleosuite mnase-suite \
  --bam sample.bam \
  --fasta hg19.fa \
  --resource-set hg19-gm12878 \
  --contigs chr1-22,chrX \
  --cores 8 \
  --outdir sample_mnase_suite
```

## Generate several track classes in one pass

```mermaid
flowchart LR
    A[Paired-end BAM or fragment BED] --> B[tracks]
    B --> C[PNS and WPS]
    B --> D[Coverage and dyads]
    B --> E[Fragment ends]
    B --> F[Sequence profiles]
    C --> G[Peak calls]
```

Use [`tracks`](commands/tracks.md) when several outputs should share the same fragment filtering but you do not need the full suite. A fragment is read once per chunk and can contribute to every requested fragment-length range that contains it.

```bash
nucleosuite tracks \
  --bam sample.bam \
  --fasta hg19.fa \
  --output-dir sample_tracks \
  --output-prefix sample \
  --fragment-range "137-197=pns,posPNS,coverage,pns_peaks" \
  --fragment-range "120-180=wps,wps_smoothed,mWPS,sm_mWPS,wps_peaks" \
  --fragment-range "145-147=dyad,fragment_left_ends,fragment_right_ends"
```

## Generate a nucleosome signal and call peaks

```mermaid
flowchart LR
    A[Paired-end BAM or fragment BED] --> B[PNS or WPS]
    B --> C[BigWig signal]
    C --> D[call-peaks]
    D --> E[Nucleosome and breakpoint BED files]
```

`pns` and `wps` can call peaks during signal generation. [`call-peaks`](commands/call-peaks.md) applies the same callers to an existing compatible BigWig.

```bash
nucleosuite call-peaks \
  --input-bigwig sample_pns.bw \
  --method pns \
  --signal both \
  --out-prefix sample_peaks
```

## Compare nucleosome spacing between chromatin states

```mermaid
flowchart LR
    A[Nucleosome peak BED] --> C[distances]
    B[Chromatin-state BED] --> C
    C --> D[State-specific spacing tables]
    C --> E[State-specific spacing plots]
```

The bundled GM12878 hg19 ChromHMM annotation can be passed directly into the command. `resources path` prints the installed path, and `$(...)` substitutes that path into `--state-bed`.

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --position-column 7 \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --min-distance 120 \
  --max-distance 250 \
  --output-prefix spacing_by_state
```

[`peak-states`](commands/peak-states.md) counts peaks within each state and measures how state representation changes with peak score.

```bash
nucleosuite peak-states sample_nucleosome_regions.bed \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --output-prefix peak_state_counts
```

## Compare flanking nucleosome spacing across reference-site categories

```mermaid
flowchart LR
    A[Nucleosome call BED] --> C[flank-spacing]
    B[Categorized reference-site BED] --> C
    C --> D[Per-site flanking spacings]
    C --> E[Category distributions]
    C --> F[Ranked category summary]
```

Use [`flank-spacing`](commands/flank-spacing.md) when each reference site belongs to a category and the biological quantity of interest is the distance between the nearest nucleosome strictly upstream and the nearest nucleosome strictly downstream. Categories are read from BED column 4 by default.

```bash
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosome_regions.bed \
  --region-bed categorized_sites.bed \
  --category-col 4 \
  --output-prefix categorized_flank_spacing
```

The default figure uses density curves, evaluates each curve at 190 and 260 bp, ranks categories by `y(190) / y(260)` from lowest to highest, highlights the top seven categories, and shows 0-500 bp on the x-axis. Raw count curves are available with `--distribution count`.

## Detect repeating spacing with DAC and estimate NRL

```mermaid
flowchart LR
    A[Dyad or other positional BigWig] --> B[DAC]
    B --> C[Distance autocorrelation profile]
    C --> D[NRL]
    D --> E[Retained repeat peaks]
    D --> F[Recurring-period estimate]
```

[`dac`](commands/dac.md) measures recurrence of the same signal at each distance. A periodic signal produces DAC peaks at the repeat distance and its multiples. [`nrl`](commands/nrl.md) then fits those repeated maxima to estimate one recurring period.

```bash
nucleosuite dac \
  --bigwig sample_dyad.bw \
  --chrom-sizes sample.bam \
  --scope combined_chromosomes \
  --dmax 2000 \
  --out-prefix sample_dac

nucleosuite nrl sample_dac.tsv \
  --peak-resolution 160 \
  --output-prefix sample_nrl
```

See [Distance autocorrelation](ALGORITHMS.md#distance-autocorrelation) for a worked 185 bp example and figure.

## Compare two positional signals with DCC

```mermaid
flowchart LR
    A[Signal A] --> C[DCC]
    B[Signal B] --> C
    C --> D[Signed or absolute lag profile]
    D --> E[Preferred A-to-B separation]
```

[`dcc`](commands/dcc.md) measures where one signal occurs relative to another.

```bash
nucleosuite dcc bigwig \
  --bigwig-a short_fragment_dyad.bw \
  --bigwig-b long_fragment_dyad.bw \
  --chrom-sizes sample.bam \
  --signed-lags \
  --dmax 500 \
  --out-prefix short_vs_long
```

Positive signed lag places signal B downstream of signal A in the active coordinate orientation.

## Examine signal around CTCF sites or other annotations

```mermaid
flowchart LR
    A[BigWig signal] --> C[aggregate]
    B[Reference-site BED] --> C
    C --> D[Per-region matrix]
    C --> E[Heatmap]
    C --> F[Mean profile]
```

Use [`aggregate`](commands/aggregate.md) when you want a heatmap plus the average signal around many sites. The bundled GM12878 CTCF sites can be inserted directly:

```bash
nucleosuite aggregate \
  --bigwig sample_pns.bw \
  --region-bed "$(nucleosuite resources path gm12878-hg19-ctcf)" \
  --strand-col 6 \
  --window-half 2500 \
  --output-dir sample_ctcf_aggregate \
  --output-prefix sample_ctcf
```

[`region-extract`](commands/region-extract.md) exports every region's signal vector and nearby peak records.

## Compare one main peak callset with multiple callsets

```mermaid
flowchart LR
    A[Main nucleosome BED] --> C[compare-positions]
    B1[Comparison BED 1] --> C
    B2[Comparison BED 2] --> C
    B3[Comparison BED ...] --> C
    C --> D[One-to-one matched pairs per comparison]
    C --> E[Combined distance distributions]
    C --> F[Main-score percentile groups]
    F --> G[Grouped percentile boxplot]
    F --> H[Optional within-percentile pairwise tests]
    F --> I[1%-percentile median distance + IQR trend]
```

[`compare-positions`](commands/compare-positions.md) compares one main nucleosome BED with one or more comparison BEDs using one-to-one matching. Matched pairs are grouped by the main BED score, with quartiles used by default. `--stats` performs pairwise tests separately within each percentile group. Main and comparison labels can be supplied as `LABEL=path.bed`.

## Analyse fragment sequence periodicity

```mermaid
flowchart LR
    A[Fragment BED or BAM] --> C[dinuc-profile]
    B[Reference FASTA] --> C
    A --> D[ww-types]
    B --> D
    C --> E[All 16 dinucleotide profiles]
    D --> F[WW/SS classes and type-specific outputs]
```

Use [`dinuc-profile`](commands/dinuc-profile.md) to measure where each dinucleotide occurs relative to the fragment centre. Use [`ww-types`](commands/ww-types.md) when you want to classify fragments by the WW/SS pattern in the centred 147-bp reference core and generate type-specific outputs.

## Compare fragment-length profiles across chromatin states

```mermaid
flowchart LR
    A[Fragments] --> C[fragment-lengths]
    B[Chromatin-state BED] --> C
    C --> D[Length counts by state]
    D --> E[fragment-heatmap]
    E --> F[State-by-length heatmap]
```

Count fragment lengths within the bundled state annotation:

```bash
nucleosuite fragment-lengths \
  --bam sample.bam \
  --bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --bed-label-column 4 \
  --output sample_state_lengths.tsv
```

Then use [`fragment-heatmap`](commands/fragment-heatmap.md) to compare the resulting state-specific profiles. Row-percentage normalization shows the fragment-length distribution within each state; fragment-length z-scores highlight states that are relatively enriched or depleted for each length.

## Gene-centred analysis with bundled resources

```mermaid
flowchart LR
    A[Signal or distance profile] --> C[Gene-centred analysis]
    B[Bundled hg19 genes and expression] --> C
    C --> D[TSS-expression quintiles]
    C --> E[Gene-expression correlation or ranking]
```

The bundled hg19 gene BED and HPA tissue-expression table can be passed directly into commands:

```bash
GENES="$(nucleosuite resources path hg19-genes)"
EXPR="$(nucleosuite resources path hpa-tissue-expression)"
```

Use [`tss-expression-quintiles`](commands/tss-expression-quintiles.md) to group genes by expression and compare signal around their TSSs. Use [`gene-expression`](commands/gene-expression.md) to relate expression profiles to gene-level spacing or periodicity measures.

## Chromosome-wise execution and combination

```mermaid
flowchart TB
    A[Selected chromosomes or scaffolds] --> B[Per-contig worker jobs]
    B --> C[Per-contig outputs and sufficient statistics]
    C --> D[combine]
    D --> E[Combined tables and intervals]
    E --> F[Combined BigWigs and bigBeds]
```

Many commands can split indexed inputs by contig. Per-contig calculations are completed first, then raw counts, products, opportunities, and interval records are combined before derived percentages or normalized values are recalculated.

To defer the combine stage:

```bash
nucleosuite dyads ... --cores 4 --skip-combine
nucleosuite combine --input-dir sample_dyads_multicontig
```

If the per-contig work is already complete, [`combine`](commands/combine.md) can reconstruct the combined outputs without rerunning the analysis itself.

## Randomized controls

```mermaid
flowchart LR
    A[Observed fragments] --> B[randomize-fragments]
    B --> C[Materialized randomized fragment set]
    C --> D[Same downstream tracks or suite]
    D --> E[Randomized-control outputs]
```

[`randomize-fragments`](commands/randomize-fragments.md) changes fragment coordinates while preserving selected fragment properties and placement constraints. Process observed and randomized fragments with the same downstream settings for a positional null comparison.

The `cfdna-suite` and `mnase-suite` commands also provide randomized-only workflow execution with `--randomize`.
