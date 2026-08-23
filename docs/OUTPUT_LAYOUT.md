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
