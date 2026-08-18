# Full-suite output layout

`mnase-suite` and `cfdna-suite` use the numbered analysis tree below. By default, they analyse observed fragments. `--randomize` writes a validated randomized fragment set in `00_setup/` and analyses that control set in the same tree.

```text
00_setup/
00_gene_sets/
01_combined_tracks/
02_dac/
03_dcc/
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

Randomized output names, logs, completion markers, and provenance files contain `_randomized_control`.

## Setup and provenance

`00_setup/` contains chromosome sizes, normalized resources, validation records, and parameter tables. Observed support files include `analysis.chrom.sizes`, `selected.chrom.sizes`, and `run_parameters.tsv`. Randomized mode prefixes support names with the marked sample and adds `*.randomized.fragments.bed.gz`, combined randomization QC, and relocation-distance outputs. Downstream work starts after the randomized BED and its QC pass validation.

Each run writes parameters, a suite report, and a JSON manifest containing the NucleoSuite version, run mode, input identities and provenance, parameter hash, blacklist selection, completion state, and outputs. `--resume` reuses matching completed steps, `--force` reruns them, and `--dry-run` validates inputs and prints the plan.

## Combined tracks

```text
01_combined_tracks/
├── pns/
├── wps/
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

In randomized runs, the manifest and completion-report filenames are prefixed with the `_randomized_control` sample name.

The MNase suite writes exact 145 and 147 bp tracks and a 145–147 bp range. Its PNS range is 120–180 bp with mode 147 bp.

The cfDNA suite writes exact 145, 161, and 167 bp tracks and the ranges 145–147, 160–162, and 166–168 bp. Its PNS range is 137–197 bp with mode 167 bp.

Both suites store primary coverage with the PNS outputs. WPS and its auxiliary coverage and dyad tracks use 120–180 bp fragments and a 120 bp protection window.

## Downstream organisation

DAC uses dyad and WW-type dyad signals. DCC is organised by signal family beneath `03_dcc/`; the cfDNA workflow calculates pairwise dyad and same-side fragment-end comparisons. NRL and periodicity outputs preserve the relative DAC or DCC path beneath `04_nrl/from_dac/` and `04_nrl/from_dcc/`.

Track-dependent stages are grouped by their input signal:

```text
05_ctcf_aggregation/{pns,wps,dyads,type_dyads}/
06_tss_aggregation/<signal>/<gene-set>/
07_distances/{pns_peaks,wps_peaks}/
08_region_extract/ctcf/{pns,wps}/
11_gene_expression/{pns,wps}/
12_positive_runs/{pns,wps}/
13_peak_analysis/{pns,wps,pns_vs_wps}/
```

Fragment-length products are grouped by the regions counted:

```text
09_fragment_lengths/combined_chromosomes/
09_fragment_lengths/chromhmm_states/
10_fragment_heatmaps/combined/
```

## Multicontig runs

With multiple selected contigs, the default `--analysis-scope combined-only` writes prerequisites beneath `per_contig/<contig>/`, combines complete contributions beneath `combined/`, and runs downstream analyses once on the pooled chromosomes. `combined/00_setup/combined_chromosomes.tsv` records included and skipped chromosomes; randomized runs prefix the filename with the marked sample name.

`--cores` controls per-contig processing, memory-light analyses, and streaming combines. `--analysis-cores` and `--streaming-combine-cores` override the corresponding stages. Indexed BigWig/bigBed creation uses `--indexed-combine-cores 1` by default. Whole-callset and other memory-heavy analyses use `--memory-intensive-analysis-cores 1` by default. Increase either value explicitly when sufficient memory is available.

Each top-level command writes a timestamped log containing the command, version, resolved parameters, working directory, console messages, exit status, and elapsed time. Output-adjacent logs are stored beneath `logs/commands/`; commands without an explicit destination use `nucleosuite_logs/` in the working directory. Suite steps also retain their dedicated logs.

Per-contig track directories contain completion checkpoints. Track writers use `.partial` files until an output is complete. Restarted commands report whether each checkpointed stage was completed, reused, or rerun.

[Back to the documentation index](README.md)
