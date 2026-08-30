# NucleoSuite workflows

This page shows how NucleoSuite commands connect. Each diagram is followed by a command-line example. Links alongside the diagrams open the corresponding command documentation.

## cfDNA: from fragments to nucleosome organization

```mermaid
flowchart TB
    A[Paired-end cfDNA BAM or fragment BED] --> B[cfdna-suite]
    B --> C[PNS mode 167]
    B --> D[Dyads and fragment ends]
    B --> E[Dinucleotide and WW/SS analyses]
    C --> F[Peak calls]
    D --> G[DAC]
    G --> I[NRL]
    F --> H[Peak spacing and NRL regression]
    C --> S[Native PNS after combination]
    S --> J[CTCF, TSS, chromatin-state, and expression analyses]
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

### Observed plus randomized suite execution

```mermaid
flowchart TB
    A[Observed fragments] --> B[Full observed suite]
    A --> C[Coordinate randomization]
    C --> D[Full randomized suite]
    B --> E[Observed combined peaks]
    D --> F[Randomized combined peaks]
    E --> G[empirical-peak-fdr]
    F --> G
    G --> H[All peaks plus empirical p and FDR]
```

Both coordinated suites accept `--with-randomized-control`. The randomized workflow uses the same filtering, track, scaling, and peak-calling settings. All observed combined peaks are retained with empirical p-value and FDR appended. Supplying `--fdr` also writes an additional filtered BED.

### Nucleosome scoring followed by spacing analysis

```mermaid
flowchart LR
    A[Paired-end BAM or fragment BED] --> B[pns command: PNS default]
    B --> C[PNS BigWig]
    B --> D[Nucleosome-region calls]
    D --> E[distances]
    E --> F[Nearest-neighbour spacing profile]
```

PNS assigns +100/−100 mass to each complete sinusoidal kernel, giving total absolute mass 200, and retains native genomic scores. The detailed geometry and explanatory plots are in [Algorithms](ALGORITHMS.md#probabilistic-nucleosome-scoring).

For example:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta hg19.fa \
  --mode 167 \
  --out-prefix sample_score

nucleosuite distances sample_score_methodpns_mode167_lower137_upper197_smooth0x2_nucleosome_regions.bed \
  --position-column 7 \
  --min-distance 120 \
  --max-distance 250 \
  --output-prefix sample_spacing
```

## Matched CUT&RUN or CUT&Tag

```mermaid
flowchart TB
    A[Condition 1 treatment/control BAMs] --> B[cutn-suite Stage 1]
    B --> C[PNS, coverage, gated peaks and clusters]
    D[Condition 2 treatment/control BAMs] --> E[cutn-suite Stage 1]
    E --> F[PNS, coverage, gated peaks and clusters]
    C --> G[cutn-compare Stage 2]
    F --> G
    G --> H[Differential cluster loci, overlap summaries and matched aggregates]
```

[`cutn-suite`](commands/cutn-suite.md) uses PNS discovery over the resolved mode ±30 bp (`--frag-mode-padding 30`) while broad coverage uses 1–1,000 bp fragments. Both ranges are generated together by `tracks` in one fragment pass. Native PNS replicate tracks are averaged directly; PNS and `posPNS` BigWigs are retained without score normalization. Sequencing depth can affect PNS amplitude. Coverage is independently scaled to a non-zero mean of 100 for Stage 1 measurement.

The default replicate statistic is mean coverage across each treatment-defined peak; `--stage1-coverage-statistic max` selects the maximum. With at least three replicates in each group, the automatic seed rule uses one-sided raw p-values plus mean treatment > mean control. With fewer replicates, a seed requires every treatment replicate to exceed every control replicate. Members extend seeded clusters subject to the selected gate, the non-member-gap allowance and the maximum interpeak distance. See the [CUT&RUN/CUT&Tag command guide](commands/cutn-suite.md) for all replicate-aware rules and the [S/G examples](ALGORITHMS.md#sg-clustering-examples) for how seeds, gated members, and intervening peaks form clusters.

Treatment and control groups can differ in size; `--bam-mode merged` pools each group. Stage 2 reads saved tracks and manifests without revisiting BAMs. Use an explicit mode to bypass estimation:

```bash
nucleosuite cutn-suite \
  --treatment1-bam target.bam \
  --control1-bam control.bam \
  --outdir target_cutn_suite \
  --mode 167
```

Provide all four treatment/control groups to run both stages together. Alternatively, run each Stage 1 independently and compare its manifest later:

```bash
nucleosuite cutn-compare \
  --condition1-results wild_type_stage1 \
  --condition2-results mutant_stage1 \
  --outdir mutant_vs_wild_type
```

## MNase-seq: from protected fragments to nucleosome organization

```mermaid
flowchart TB
    A[MNase BAM or fragment BED] --> B[mnase-suite]
    B --> C[PNS mode 147]
    B --> D[Coverage, dyads, and fragment ends]
    B --> E[Dinucleotide and WW/SS analyses]
    C --> F[Peak calls]
    D --> G[DAC]
    G --> N[NRL]
    F --> H[Peak spacing and NRL regression]
    C --> S[Native PNS after combination]
    S --> I[CTCF, TSS, chromatin-state, and expression analyses]
```

[`mnase-suite`](commands/mnase-suite.md) applies MNase-oriented defaults and uses one filtered fragment population across downstream outputs.

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
    A[Paired-end BAM or fragment BED] --> B[Nucleosome score or WPS]
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

At resolution 160 bp, NRL uses 61 bp smoothing to find broad peaks, 21 bp smoothing to refine their summits, and a minimum peak separation of 160 bp. The fitted repeat length still comes from the observed peak spacing. The [resolution table](commands/nrl.md#how-it-works) shows how another resolution changes these settings.

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

The `cfdna-suite` and `mnase-suite` commands provide randomized-only workflow execution with `--randomize` and paired full execution with `--with-randomized-control`.

When observed and randomized peak calls already exist, [`empirical-peak-fdr`](commands/empirical-peak-fdr.md) reports pooled empirical p-values and monotonic empirical FDR values without positional matching:

```bash
nucleosuite empirical-peak-fdr observed_peaks.bed randomized_peaks.bed --fdr 0.05
```
