# Output layout

NucleoSuite keeps source-derived outputs, combined outputs, reports, and figures in predictable directories. Exact filenames include the selected sample, mode, fragment bounds, smoothing, and requested output family.

## Single-command outputs

For a single-contig `pns` run, the requested prefix receives:

```text
sample_pns.bw
sample_posPNS.bw
sample_pns_smoothed.bw       # when smoothing is requested
sample_nucleosome_regions.bed|bb
sample_breakpoint_peaks.bed|bb
sample_coverage.bw
sample_dyad.bw
sample_fragment_ends.bw
sample_fragment_left_ends.bw
sample_fragment_right_ends.bw
sample_fragment_mode_estimation.tsv
```

Only requested tracks and reports are written. The signed PNS and `posPNS` values are native outputs; they are not divided by a reference mean or otherwise rescaled by `pns`.

## Multicontig commands

Commands that support chromosome-wise execution use a command-specific parallel directory and a combined output directory:

```text
<output>/
├── parallel/
│   ├── chr1/
│   ├── chr2/
│   └── nucleosuite_multicontig_manifest.json
└── combined/
    ├── <combined tracks and intervals>
    └── <completion reports>
```

The manifest records input signatures, resolved parameters, contigs, and per-contig output paths. `combine` validates the per-contig set before writing combined tracks and interval files. `--skip-combine` leaves the per-contig outputs available for a later combine stage.

## `tracks` outputs

`tracks` groups outputs by the requested fragment range and output prefix. A range may contain PNS, `posPNS`, WPS, coverage, dyads, fragment ends, dinucleotide profiles, WW/SS summaries, type-specific dyads, and PNS/WPS peak calls. Exact lengths and ranges can coexist in the same run.

## `cutn-suite` outputs

A Stage 1 run is organized as:

```text
cutn_results/
├── 00_setup/                  input, mode, and track specifications
├── 01_score_tracks/           per-replicate native PNS, posPNS, and coverage
├── 02_mean_scaled_tracks/     condition means and normalized coverage products
├── 03_peak_calls/             treatment candidates and breakpoint calls
├── 04_peak_statistics/        replicate measurements and gate results
├── 05_cluster_aggregate/      cluster-centred PNS aggregates and NRLs
├── cutn_stage1_manifest.json
└── cutn_suite_summary.tsv
```

The `02_mean_scaled_tracks` name identifies the coverage-normalization stage and is retained for layout compatibility. PNS score BigWigs and PNS peak scores are not mean-scaled. Coverage alone is normalized to a non-zero mean of 100 for Stage 1 treatment/control measurement.

With biological replicates, the manifest retains each replicate’s native PNS/`posPNS` paths, normalized coverage path, condition-level PNS score path, modes, fragment ranges, gate settings, cluster settings, and downstream output paths. This allows `--rerun-from` and `cutn-compare` to reuse compatible files without reopening the BAMs.

Stage 2 adds comparison tables, all/robust/significant gain and loss BEDs, cluster-overlap summaries, and matched condition aggregates beneath the comparison output directory.

## Coordinated cfDNA and MNase workflows

The cfDNA and MNase suites create workflow-specific directories for fragment preparation, PNS/coverage and coordinate tracks, sequence profiles, peak analysis, spacing, regional aggregation, gene analyses, and plots. Their PNS directories use `pns` and `posPNS` names. Raw PNS score BigWigs and PNS peak BEDs are passed to downstream analyses; coverage normalization remains a separate product.

## Replotting and reproducibility

Reports and TSV outputs are the inputs to `plot`, so figures can be recreated without rerunning fragment processing. Seeded mode estimation, randomized controls, and saved manifests preserve the settings needed to reproduce a run.
