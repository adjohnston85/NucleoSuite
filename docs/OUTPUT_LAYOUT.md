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

Native PNS, `posPNS`, coverage, nucleosome-region and breakpoint-peak outputs are written beneath `pns/`. After chromosome combination, `scaled/` receives mean-scaled coverage. Downstream peak analyses use the native peak BEDs and aggregate analyses use native PNS. The score tracks and peak scores are not automatically rescaled.

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

## CUT&RUN/CUT&Tag suite

`cutn-suite` uses a separate target/control layout:

| Path | Contents and role |
|---|---|
| `00_setup/` | Fragment-mode and track-processing reports. |
| `01_score_tracks/` | Replicate PNS/`posPNS` for positioning and raw broad-range coverage for abundance. |
| `02_analysis_tracks/` | Native condition-mean PNS and replicate/condition-mean coverage scaled to a non-zero mean of 100. |
| `03_peak_calls/` | Treatment nucleosome candidates defined by the consensus PNS track. |
| `04_peak_statistics/` | Replicate interval measurements, gate results, and seeded clusters. |
| `05_cluster_aggregate/` | Strongest-member anchors, PNS heatmaps/profiles, confidence bands, and directional NRLs. |
| `<sample>_cutn_suite_summary.tsv` | Summary of the Stage 1 analysis. |
| `cutn_stage1_manifest.json` | Saved Stage 1 parameters and output paths needed by Stage 2 and reruns. |
| `cutn_suite_run_manifest.json` | Run-level condition membership, parameters, and Stage 1/Stage 2 manifests. |

With biological replicates, `01_score_tracks/` retains native PNS/`posPNS` and raw broad coverage. `02_analysis_tracks/` contains native condition-mean PNS plus normalized replicate/condition-mean coverage. PNS and `posPNS` are not divided by a reference mean. Broad coverage is scaled to a non-zero mean of 100 for Stage 1 interval measurement, using the mean across each peak by default. Automatic S/G clustering rules depend on replicate count and are recorded in the Stage 1 manifest. Native PNS tracks are reused for cluster heatmaps, aggregate profiles, confidence bands and directional NRLs. Explicit and estimated modes are recorded in the mode report.

A four-group run places the two layouts under `01_condition1_stage1/` and `02_condition2_stage1/`. `03_condition_comparison/` contains the complete cluster-only differential table, all-direction, robust-direction and FDR-significant BEDs, overlap-component mapping, Venn and occupied-base summaries, matched cluster-locus PNS aggregates, and `cutn_comparison_manifest.json`. The root `cutn_suite_run_manifest.json` records both biological conditions, their BAM membership, the Stage 1 manifests, run-level parameters and the Stage 2 manifest. The standalone `cutn-compare` command writes the same Stage 2 layout from two existing Stage 1 manifests.

Fast reruns created with `cutn-suite --rerun-from` are written inside the source run as `rerun_01/`, `rerun_excluding_<sample>_01/`, or `rerun_excluding_<N>_samples_01/`, with the numeric suffix incremented when a matching rerun already exists. These rerun trees deliberately omit `01_score_tracks/`: their Stage 1 manifests reference the retained per-replicate score and coverage BigWigs in the source run, while new condition-mean tracks, peak calls, replicate statistics, clusters, Stage 2 outputs and aggregates are written beneath the rerun directory. Each rerun has its own `cutn_suite_run_manifest.json` recording the source run, exclusions and downstream parameter changes.
