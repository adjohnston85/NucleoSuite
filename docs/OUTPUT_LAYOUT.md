# Full-suite output layout

`mnase-suite` and `cfdna-suite` use the same numbered analysis tree for their coordinated track, periodicity, spacing, regional and summary analyses.

```text
00_setup/
00_gene_sets/
01_combined_tracks/
02_dac/
04_nrl/
05_ctcf_aggregation/
06_tss_aggregation/
06_tss_expression_quintiles/
07_distances/
08_region_extract/
09_fragment_lengths/
10_fragment_heatmaps/
11_gene_expression/
12_positive_runs/
13_peak_analysis/
logs/
.done/
```

## Combined tracks and scaling

```text
01_combined_tracks/
├── pns/
├── scaled/
├── dyads/
│   ├── exact/<length>/
│   └── ranges/<lower-upper>/
├── fragment_ends/
│   ├── exact/<length>/
│   └── ranges/<lower-upper>/
├── sequence/
│   ├── dinucleotide_profiles/{exact,ranges}/
│   ├── ww_types/ranges/
│   ├── type_dyads/ranges/
│   └── summaries/
├── manifest.tsv
└── completion_report.tsv
```

Raw PNS, posPNS, coverage, nucleosome-region and breakpoint-peak outputs are written beneath `pns/`. After chromosome combination, `scaled/` receives mean-scaled coverage, mean-scaled posPNS, PNS scaled relative to the mean raw combined nucleosome-peak score, and mean-scaled nucleosome-region and breakpoint-peak BEDs. Downstream peak-based suite analyses use the mean-scaled peak BEDs, while PNS aggregate stages use the scaled PNS track.

MNase uses the 146–148 bp ranged class, exact 147 bp dyads/ends, and exact 145/147 bp dinucleotide profiles. cfDNA uses ranged classes 144–146, 160–162 and 166–168 bp plus exact 145, 161 and 167 bp dyads/ends.

## Downstream organisation

`02_dac/` contains DAC from ranged dyads. `04_nrl/from_dac/` mirrors the DAC range paths and stores the long, short, and intermediate periodicity fits.

```text
05_ctcf_aggregation/{pns,dyads,type_dyads}/
06_tss_aggregation/<signal>/<gene-set>/
07_distances/pns_peaks/
08_region_extract/ctcf/pns/
11_gene_expression/pns/
12_positive_runs/pns/
13_peak_analysis/score_frequencies/pns/
```

Fragment-length products remain under:

```text
09_fragment_lengths/combined_chromosomes/
09_fragment_lengths/chromhmm_states/
10_fragment_heatmaps/combined/
```

## Multicontig runs

With `--analysis-scope combined-only` (default), per-contig workers create combine prerequisites beneath `per_contig/<contig>/`; complete tracks are combined beneath `combined/`, then scaling and downstream analyses run once on the pooled selected chromosomes. `--resume` reuses matching completed work and `--force` reruns it.

Randomized runs use the same tree and mark their sample/output names with `_randomized_control`.

With `--with-randomized-control`, the observed and randomized trees are both completed before FDR annotation. Combined observed nucleosome and breakpoint BEDs with appended FDR are written beneath:

```text
combined/13_peak_analysis/pns/empirical_fdr/
```

If the suite is not using the multicontig wrapper layout, the same directory is created directly beneath the suite root.

## ChIP/CUT&RUN/CUT&Tag suite

`chip-suite` uses a separate target/control layout:

```text
00_setup/                 fragment-mode and score-scaling reports
01_score_tracks/          method-specific discovery score/positive score plus raw broad-range coverage
02_mean_scaled_tracks/    normalized method-specific score tracks and coverage scaled to mean 100
03_peak_calls/            treatment nucleosome candidates only
04_peak_fdr/              replicate statistics, gate-selected peaks and seeded clusters
05_cluster_aggregate/     strongest-member anchors, PNS heatmap/profiles, confidence band, NRLs
<sample>_chip_suite_summary.tsv
chip_stage1_manifest.json
```

With biological replicates, `01_score_tracks/` and `02_mean_scaled_tracks/` retain replicate-specific outputs. The resolved mode ±30 bp score/positive-score pair and broad 1–1,000 bp coverage are generated together in one `tracks` pass for each replicate. Each selected score is normalized by its matching positive-score mean before treatment tracks are averaged, preventing a deeper replicate from dominating consensus candidate discovery. Only nucleosome peaks from the mean treatment score track supply Stage 1 candidates. Replicate-specific broad-range coverage scaled to a non-zero mean of 100 supplies treatment/control measurements and exploratory Welch/BH annotations, while condition-mean treatment coverage supplies the reported peak score. Every treatment replicate > every control replicate (`all-controls`) is the default gate; mean mode is optional. Clusters start at gate-passing p < 0.05 seeds. The default member mode includes S and G peaks, one non-member bridge is allowed, at least two included members are required, and adjacent included-member summits cannot exceed 1,000 bp. The selected normalized PNS/BNS/TNS tracks are reused for cluster heatmaps, aggregate profiles, confidence bands, and directional NRLs. An explicit `--mode` is recorded with `mode_source=explicit`; automatic runs retain treatment, control, and pooled estimates together with the unsmoothed-by-default histogram setting.

A four-group run places the two layouts under `01_condition1_stage1/` and `02_condition2_stage1/`. `03_condition_comparison/` contains the complete cluster-only differential table, all-direction, robust-direction and FDR-significant BEDs, overlap-component mapping, Venn and occupied-base summaries, matched cluster-locus PNS aggregates, and `chip_comparison_manifest.json`. The manifest records the log-scale empirical-Bayes interaction model, variance prior, overlap topology, and aggregate paths. The standalone `chip-compare` command writes the same Stage 2 layout from two existing Stage 1 manifests.
