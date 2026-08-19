# File formats

Use this page to decide which file type a command expects and what the columns mean.

## Practical summary

| File | Contains | Common NucleoSuite use |
|---|---|---|
| BAM + BAI | aligned paired-end reads | original input for fragment-based analyses |
| BED/BED.gz | genomic intervals in text form | fragments, regions, genes, states and peak calls |
| bigBed | indexed binary BED | large interval files and genome-browser viewing |
| BigWig | indexed numeric genomic signal | PNS, WPS, coverage, dyads, ends and downstream signal analysis |
| TSV/TSV.gz | ordinary tables | counts, summaries, correlations, matrices and statistics |
| FASTA + FAI | reference genome sequence | sequence-aware analyses and contig lengths |
| chromosome sizes, BAM or CRAM | contig names and lengths | BigWig/bigBed creation and analyses that require reference lengths |

BAM, FASTA, BED, and chromosome-size files must use the same genome assembly. NucleoSuite accepts conservative `chr`/non-`chr` aliases when equivalent contigs have identical lengths. Generated fragment and track outputs use the BAM-derived namespace.

## BAM

BAM inputs must be coordinate sorted and indexed. Commands that reconstruct paired-end fragments require mapped mates on the same contig.

With one BAM, output contig names retain the BAM-header spelling. With multiple BAMs, NucleoSuite forms the union of all header contigs. Equivalent names such as `20` and `chr20` are merged, and the `chr` spelling is preferred when both are present. Each BAM is still fetched using its own original header name. Conflicting equivalent lengths or both aliases within one BAM are rejected.

Identical-fragment deduplication is controlled by:

- `--max-duplicates 1`: retain one fragment for each identical `(contig, start, end)` coordinate.
- `--max-duplicates 0`: disable identical-fragment deduplication.
- `--dedup-scope all_bams`: apply the duplicate limit across all input BAMs.
- `--dedup-scope per_bam`: apply the limit independently to each BAM.

Commands that write sparse dyad or fragment-end tracks also accept `--max-per-coordinate`, which caps the final signal at each output coordinate. Its default is `0` (unlimited).

## Fragment-coordinate BED, BED.gz and bigBed

Fragment-aware commands accept materialised fragment intervals in BED, gzip-compressed BED, or bigBed format. Only the first three columns are required:

```text
chrom    start    end    [any additional columns ...]
```

Coordinates are zero-based and half-open, and fragment length is calculated as `end - start`. Additional columns are ignored unless a command documents a separate use for them. Invalid records (`start < 0`, `end <= start`, or non-integer coordinates) are rejected with the source filename and line number.

For plain BED/BED.gz input, supply complete chromosome lengths through `--chrom-sizes` as a two-column table, BAM, or CRAM, or through `--fasta`. If neither is supplied, NucleoSuite uses the largest observed fragment end on each contig. BigBed input is converted temporarily with UCSC `bigBedToBed`.

Multiple fragment files may be combined. `--dedup-scope all_bams` applies the coordinate limit across the complete input collection, while `per_bam` applies it independently to each input. The option values have the same meaning for BAM, BED, BED.gz, and bigBed inputs.

Large materialized fragment files are indexed temporarily in SQLite for chunk-based retrieval. The index stores both the source chromosome spelling and the canonical analysis spelling. It is deleted when the command finishes and is not a final output.

## BED and BED.gz

BED coordinates are zero-based and half-open. BED files should be tab-delimited.

```text
chrom    start    end    [name]    [score]    [strand]    ...
```

Blank lines, lines beginning with `#`, and UCSC `track` or `browser` lines are ignored by commands that support them.

### Region BED

Region-based commands require at least BED3:

```text
chrom    start    end
```

`aggregate` uses configurable one-based BED column numbers. Its defaults are chromosome column 1, start column 2, end column 3, and strand column 6. `--point-col 0` uses the interval midpoint; a positive column number reads an absolute genomic position from that column.

### Categorized reference-site BED

`flank-spacing` accepts an ordinary BED as its reference-site input. The first three columns are chromosome, start, and end. By default, column 4 supplies the category label used to group sites; select another one-based column with `--category-col`. The reference coordinate is the interval midpoint unless `--point-col` selects an exact coordinate column.

A typical BED6 input is:

```text
chr1	100000	100001	category_A	0	+
chr1	120000	120001	category_B	0	-
```

The nucleosome BED supplied to `--nucleosome-bed` is also interpreted from interval midpoints by default. `--nucleosome-center-col` can select an explicit centre column when the callset stores one.

### State BED

State or group annotations require at least BED4:

```text
chrom    start    end    state
```

DAC and DCC use column 4 as the state name by default. `distances` accepts a state BED through `--state-bed` and allows the label column to be changed with `--state-label-column`.

For ChromHMM-coloured relative-distance overlays, use BED9 and store the RGB colour in column 9:

```text
chrom  start  end  state  score  strand  thickStart  thickEnd  itemRgb
```

`itemRgb` must contain comma-separated red, green and blue values such as `10,190,254`.

### Gene BED

The bundled hg19 gene table was derived from Ensembl release 87 (`Homo_sapiens.GRCh37.87.gtf`; Aken et al., 2017), and `gene-expression` uses one unique record per gene with six columns:

```text
chrom    start    end    ensembl_gene_id    gene_name    strand
```

`gene-sets` uses column 4 as its default identifier. `gene-expression` uses columns 4, 5, and 6 as the Ensembl identifier, display name, and strand.

### Long-format expression TSV

`gene-expression` and the expression stage of both `mnase-suite` and `cfdna-suite` accept a tab-delimited long-format table. The suite stage runs automatically whenever `--expression` is supplied. Default columns are:

```text
Gene    Gene name    Cell line    TPM    pTPM    nTPM
```

Each row represents one gene/profile combination. The selected expression value is numeric and non-negative. `nTPM` is the default; `TPM` and `pTPM` remain selectable alternatives. Column names are configurable.

### Gene-set configuration TSV

Gene-set configurations require `set_name` and `include_rule`. The optional `exclude_if_candidate` column lists candidate sets that disqualify a gene from the current final category.

```text
set_name         include_rule                                                           exclude_if_candidate
active_genes     1_Active_Promoter & (9_Txn_Transition | 10_Txn_Elongation)            repressed_genes
weak_genes       2_Weak_Promoter & (9_Txn_Transition | 10_Txn_Elongation | 11_Weak_Txn) active_genes,repressed_genes
repressed_genes  12_Repressed                                                           active_genes,weak_genes
```

Rules support `&`, `|`, and parentheses. `--leftover-set-name leftover_genes` assigns genes that belong to none of the configured candidate sets. A candidate-overlap gene removed from named final categories remains recorded as unassigned and is excluded from the leftover set.

### Gene-set interval BED6

Candidate and final per-category files use standard BED6 with the Ensembl gene identifier in column 4:

```text
chrom  start  end  Ensembl_gene_ID  0  strand
```

`gene-sets` also writes a pooled-analysis BED6 with the final category in column 4:

```text
chrom  start  end  category  0  strand
```

Gene names, Ensembl identifiers, candidate memberships and final categories are retained together in `gene_sets_gene_assignments.tsv`.

The default category labels are `active_genes`, `weak_genes`, `repressed_genes`, and `leftover_genes`.

### Peak BED8

`pns`, `wps`, `call-peaks`, and filtered peak output from `distances` use BED8:

```text
chrom    chromStart    chromEnd    name    score    strand    thickStart    thickEnd
```

| Column | Field | NucleoSuite value |
|---:|---|---|
| 1 | `chrom` | Contig name |
| 2 | `chromStart` | Peak-region start |
| 3 | `chromEnd` | Peak-region end |
| 4 | `name` | Peak identifier |
| 5 | `score` | PNS text BED: six-decimal floating-point peak score; WPS and bigBed: integer BED score from 0 to 1000 |
| 6 | `strand` | `.` unless strand is defined by the input |
| 7 | `thickStart` | Representative call centre |
| 8 | `thickEnd` | Representative call centre plus one base |

For PNS, column 7 is the retained positive- or negative-region midpoint; for WPS it is the selected above-median subrun midpoint. PNS BED files preserve six decimal places in column 5. PNS peak bigBed output multiplies that score by `--bigbed-score-scale` (default 1000), rounds to the nearest integer, and clamps it to 0–1000. Commands that accept peak tracks, including `region-extract`, use column 7 by default.

### Fragment BED outputs

- Unclassified fragments are written as BED3: `chrom`, `start`, `end`.
- Combined WW-type fragments are written as BED4, with the WW type in the `name` field.
- Type-specific fragment files are written as BED3.

## Chromosome-wise and combined outputs

When a command is run on multiple chromosomes or scaffolds with `--cores`, NucleoSuite analyses each selected reference sequence separately and writes a manifest-backed directory tree:

```text
<parallel-root>/
├── per_contig/
│   ├── chr1/
│   ├── chr2/
│   └── ...
├── combined/
└── nucleosuite_multicontig_manifest.json
```

The manifest records the command, selected reference sequences, chromosome-specific output locations, and combination strategy. `nucleosuite combine --input-dir <parallel-root>` can rerun the combination stage without repeating completed chromosome-specific prerequisite generation.

Combined statistics are recalculated from their underlying counts, products, denominators, matrices, or records:

- DAC and DCC raw values and opportunity counts are summed before normalized values and percentages are recalculated.
- Dinucleotide counts and the numbers of fully canonical retained fragments spanning each profile position are summed before frequencies are recalculated.
- Fragment-length and distance-histogram counts are summed before percentages or plots are regenerated.
- Aggregate matrices retain per-region observations and recompute the combined mean profile.

Combined BigWigs are created after tabular and interval outputs. The default method reads per-contig BigWig intervals in 100,000 bp chunks and writes them to a new BigWig with pyBigWig. `--combine-bigwig-method bedgraph` writes validated per-contig bedGraphs under `combined/temporary_bedgraph_combine/per_contig/`, concatenates them in reference order, and converts them with `bedGraphToBigWig`. Suite progress is recorded in `combined/<sample>_combine_steps.log`. Temporary files are deleted after successful verification and retained after failure.

## BigWig

BigWig and bigBed are indexed binary genomic formats described by Kent et al. (2010). BigWig files use the `.bw` extension and are read or written with pyBigWig. Newly generated BigWigs use the BAM-derived canonical namespace. Analyses resolve exact names first and then conservative `chr`/non-`chr` and mitochondrial aliases when matching support files or existing tracks.

## Compressed WIG

Track-producing commands can write gzip-compressed WIG files with the `.wig.gz` extension when `wiggz` or `both` is selected. Dense tracks use fixed-step WIG records; sparse tracks use variable-step records.

## Chromosome sizes

Commands that accept `--chrom-sizes` can use either a two-column chromosome-size table or an alignment file:

```text
chr1    248956422
chr2    242193529
```

The first column is the reference-sequence name and the second is its length in bases. When a BAM or CRAM is supplied, NucleoSuite reads the names and lengths from the alignment header. Reference order is preserved.

In observed `mnase-suite` and `cfdna-suite` runs, `00_setup/analysis.chrom.sizes` records the complete canonical union from all BAM headers and `00_setup/selected.chrom.sizes` records only selected contigs. Randomized runs prefix both support filenames with the `_randomized_control` sample name.

A chromosome-size table can also be generated directly:

```bash
nucleosuite chrom-sizes \
  --bam sample.bam \
  --output sample.chrom.sizes
```

CRAM input may require the corresponding reference FASTA.

## Distance-order NRL regression outputs

When `nucleosuite distances` analyses more than one neighbour order, it writes regression outputs separately for the genome and each chromosome according to `--scope`.

Each regression points TSV contains:

```text
scope	chromosome	order	peak_distance_bp	peak_count	total_pairs	fitted_distance_bp	residual_bp
```

The accompanying summary TSV contains one row per completed regression and reports `nrl_bp` as the fitted slope, together with the intercept, R-squared, order range and paths to the corresponding points TSV and PNG.

## DAC and DCC profile TSVs

`nrl` accepts tab-delimited DAC or DCC output. It auto-detects:

```text
Distance    DAC Value
Lag         DCC Value
```

Use `--distance-column` or `--value-column` to select another numeric column, such as a percentage or raw signal column.

NRL outputs are headered TSV files containing the selected profile, called peaks and regression statistics. The profile table contains the unsmoothed signal, the finer local-maximum-smoothed signal, the broader peak-detection-smoothed signal, and flags for detected/refined peaks. The regression table reports `slope_bp_per_peak` as the estimated NRL or periodicity together with `peak_resolution_bp`, `detection_smoothing_window`, and `local_max_smoothing_window`.

## Position-comparison TSVs

`compare-positions` accepts one main BED, BED.gz, or bigBed plus one or more comparison interval files. The main summit/score columns and shared comparison summit/score columns are selected independently. When no summit column is supplied, the summit is the integer midpoint between BED start and end. Repeated comparison inputs can be written as `--compare-bed LABEL=path.bed` to set the plot and table label.

The `_pairs.tsv` output contains the original start, end, name, summit and score from methods A and B, together with signed and absolute summit distances and raw, z-score and percentile-rank score differences.

The `_summary.tsv` output uses two columns:

```text
metric    value
```

The `_distance_bins.tsv` output reports pair counts, summit-distance summaries, and score correlations for the selected score normalization. Distance-bin labels use discrete integer ranges such as `0-5`, `6-10`, and `11-20`.

The directional `_A_percentiles_vs_all_B_distances.tsv` and `_B_percentiles_vs_all_A_distances.tsv` outputs contain one row per matched query position. Each records the analysis direction, source percentile group, both scores and independently assigned percentiles, summits, signed distances, and absolute distance. Default group labels are `0-25`, `25-50`, `50-75`, and `75-100`; these represent lower-exclusive, upper-inclusive boundaries except at zero. The corresponding `_summary.tsv` files report source-group size, complete target-callset size, matched and unmatched counts, quartiles, mean, range, and standard deviation. Each directional distances TSV is the source data for its matching boxplot.

## Reference FASTA

Sequence-aware commands require an indexed FASTA:

```bash
samtools faidx genome.fa
```

A FASTA is required for:

- `dinuc-profile`
- `ww-types`
- `pns --dinuc-profile`
- `pns --split-ww-types`
- `--randomize-mode dinuc_anchor`
- `randomize-fragments --method dinucleotide`

## TSV

NucleoSuite tabular outputs use tab-separated values with a header row and the `.tsv` extension. Common examples include:

### Fragment-length counts

```text
fragment_length    count
```

When fragment lengths are grouped by BED labels:

```text
label    fragment_length    count
```

Fragment-size NRL output retains the data required to reproduce both plots. The profile table contains `fragment_length`, `count`, `unsmoothed_density`, `local_max_smoothed_density`, `detection_smoothed_density`, detection/refined flags and peak number. The peaks table contains `peak_number`, observed and fitted `fragment_length`, residual, and the three signal values. Regression and combined summary tables identify `nrl_method` as `fragment_size_distribution` and report the fit settings and diagnostics.

### Fragment-processing summary

```text
metric    value
```

### Fragment-heatmap metadata

```text
profile    group       condition
sampleA    Group_A     Control
sampleB    Group_B     Treatment
```

Use `--metadata-profile-column` and `--metadata-category-column` to select the identifier and grouping columns.

### Dinucleotide count tables

Dinucleotide-profile commands write `_dinuc_profile_counts.tsv` alongside the frequency profile. The count table retains per-position counts for all 16 dinucleotides, the number of fully canonical retained fragments spanning each position, and fragment-use totals. These values are used to reconstruct combined frequencies across chromosomes or scaffolds.

## Images and workbooks

- Heatmaps are written as PNG files.
- Aggregate-profile plots from `aggregate` are written as SVG files.
- Fragment-heatmap workbooks are written as XLSX files unless `--no-excel` is used.

## Common output suffixes

| Suffix | Contents |
|---|---|
| `_pns.bw` | Mean-centred PNS signal |
| `_pns_smoothed.bw` | Optional Savitzky–Golay-smoothed PNS signal; not written by default |
| `_posPNS.bw` | Endpoint-support distribution before mean subtraction |
| `_coverage.bw` | Fragment coverage |
| `_dyad.bw` | Fragment-centre signal |
| `_fragment_ends.bw` | Combined fragment-end signal |
| `_fragment_left_ends.bw` | Genomic left-end signal |
| `_fragment_right_ends.bw` | Genomic right-end signal |
| `_nucleosome_regions.bed` / `.bb` | Nucleosome-region peaks in BED8 or bigBed |
| `_breakpoint_peaks.bed` / `.bb` | Breakpoint peaks in BED8 or bigBed |
| `.fragments.bed` / `.fragments.bed.gz` / `.fragments.bb` | Materialised fragment coordinates |
| `.randomized.fragments.bed.gz` | Reusable randomized fragment control |
| `.randomization_qc.tsv` | Randomization matching, fallback, anchor and seed metrics |
| `_fragment_summary.tsv` | Fragment-processing metrics |
| `_fragment_length_counts.tsv` | Counts by fragment length |
| `_fragment_length_distribution.png` | Fragment-length count profile plotted from `_fragment_length_counts.tsv` |
| `_fragment_size_nrl_profile.tsv` / `.png` | Fragment-size density, both smoothing scales, called peaks and profile figure |
| `_fragment_size_nrl_peaks.tsv` | Called multinucleosome fragment-size peaks, fitted positions and residuals |
| `_fragment_size_nrl_regression.tsv` / `.png` | Fragment-size NRL fit settings, diagnostics and square regression figure |
| `_fragment_size_nrl_summary.tsv` | Fragment-size NRL results collected across region labels |
| `_dinuc_profile.tsv` | Dyad-aligned dinucleotide profile |
| `_dinuc_profile_counts.tsv` | Exact dinucleotide counts and contributing-fragment totals used for profile combination |
| `_dinuc_profile.png` | The 16 individual dinucleotide profiles |
| `_ww_ss_profile.png` | Aggregate WW and SS profiles |
| `_ww_type_summary.tsv` | WW-type counts and percentages |
| `_ww_type_summary.png` | WW-type count plot |
| `_ww_type_by_length.tsv` | WW/SS type counts and relative frequencies for each fragment length |
| `_ww_type_by_length_stacked.png` | Type1–type4 relative-frequency stacked bars by fragment length |
| `_score_frequency.tsv` | Exact score-frequency bins and values for peak callsets |
| `_score_frequency.png` | Peak score-frequency plot |
| `_nrl_*_profile.tsv` | Unsmoothed, local-max-smoothed, and detection-smoothed DAC/DCC profiles with detected/refined peak flags |
| `_nrl_*_peaks.tsv` | Refined periodic peaks with their finer-smoothed values and originating broad detection peaks |
| `_nrl_*_regression.tsv` | Estimated slope, intercept, R-squared and analysis settings |
| `_nrl_*_profile.png` | Unsmoothed profile plus the finer local-max and broader detection smoothing scales with called peaks |
| `_nrl_*_regression.png` | Peak-number versus distance regression plot |
| `_aggregate_nrl_profile.tsv` / `.png` | Complete signed aggregate profile, continuous smoothing layers, unified peak calls and optional shaded regression exclusion interval |
| `_aggregate_nrl_peaks.tsv` | Every aggregate peak with signed position, outward distance, directional order, shared-central status and regression inclusion/exclusion |
| `_aggregate_nrl_positive_regression.tsv` / `.png` | Positive-direction peak order versus distance from position 0, with an eligible central peak as order 0 |
| `_aggregate_nrl_negative_regression.tsv` / `.png` | Negative-direction peak order versus absolute distance from position 0, with the same eligible central peak as order 0 |
| `_aggregate_nrl_summary.tsv` | Positive and negative repeat lengths, central-peak and exclusion settings, fit diagnostics and unified caller settings |
| `_pairs.tsv` | Matched method A/B positions, scores, distances and normalized score differences |
| `_summary.tsv` | Position-comparison or other command summary metrics |
| `_distance_bins.tsv` | Score agreement and score differences by summit-distance bin |
| `_distance_histogram.tsv` | Exact histogram bin edges and matched-pair counts used for the distance histogram plot |
| `_score_correlation.png` | Method A versus method B score, coloured by summit distance |
| `_distance_histogram.png` | Distribution of matched summit distances |
| `_correlation_by_distance.png` | Score correlation within summit-distance bins |

Chromosome-specific filenames include the chromosome or scaffold identifier. Combined filenames omit that identifier because they represent the complete selected reference-sequence set.

## Parameterized output stems

Plot-producing analyses append their result-defining settings to automatic output
stems, including when `--output-prefix` or `--out-prefix` supplies the base name.
Examples include NRL resolution and fit bounds, DAC/DCC distance and normalization,
distance neighbour orders and score grouping, aggregate window/NRL/exclusion settings,
position-matching policy, fragment ranges, and signal-processing mode. Cosmetic
settings such as DPI, colours, fonts, and tick spacing are not included. An option
that names one exact output path remains exact.

## Indexed interval inputs and `--cores`

Fragment and region commands split work by contig when the primary input has a random-access index. Supported forms are bigBed and bgzip-compressed BED/TSV with a `.tbi` or `.csi` index. Tabix-indexed fragment files also require `--chrom-sizes` or `--fasta` because tabix does not store chromosome lengths. Plain BED, ordinary gzip, and unindexed TSV inputs run serially.


## References

- Aken BL, Achuthan P, Akanni W, et al. (2017). Ensembl 2017. *Nucleic Acids Research* 45(D1), D635–D642. https://doi.org/10.1093/nar/gkw1104
- Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. (2010). BigWig and BigBed: enabling browsing of large distributed datasets. *Bioinformatics* 26, 2204–2207. https://doi.org/10.1093/bioinformatics/btq351
