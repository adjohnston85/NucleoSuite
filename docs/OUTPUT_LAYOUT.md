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
01_score_tracks/          target/control score, positive-score, and raw coverage tracks
02_mean_scaled_tracks/    discovery-score normalization and coverage scaled to mean 100
03_peak_calls/            treatment candidates; legacy mode may also call control peaks
04_peak_fdr/              replicate statistics, gate-selected peaks and p-defined clusters
<sample>_chip_suite_summary.tsv
chip_stage1_manifest.json
```

With biological replicates, `01_score_tracks/` and `02_mean_scaled_tracks/` retain replicate-specific outputs. Each score track is normalized by its matching positive-score mean before treatment tracks are averaged, preventing a deeper replicate from dominating consensus candidate discovery. Only nucleosome peaks from the mean treatment score track supply Stage 1 candidates; `chip-suite` does not call breakpoint peaks. Replicate-specific coverage scaled to a non-zero mean of 100 supplies the all-treatment versus all-control gate and exploratory Welch/BH annotations, while condition-mean treatment coverage supplies the reported peak score. Gate-selected peaks proceed by default. Clusters contain gate-passing members with p < 0.05. An explicit `--mode` is recorded with `mode_source=explicit`; automatic runs retain treatment, control, and pooled estimates together with the unsmoothed-by-default histogram setting.

A four-group run places the two layouts under `01_condition1_stage1/` and `02_condition2_stage1/`. `03_condition_comparison/` contains complete differential peak and cluster tables, all-direction, robust-direction and FDR-significant BEDs, and `chip_comparison_manifest.json`. The manifest records the log-scale empirical-Bayes interaction model and variance prior. The standalone `chip-compare` command writes the same Stage 2 layout from two existing Stage 1 manifests.
