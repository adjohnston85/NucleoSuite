# NucleoSuite documentation

Choose a command by the analysis to perform and the output you need.

## Start here

- [Quick start](QUICKSTART.md) — common analyses from BAM to interpretable output.
- [Choosing a command](CHOOSING_A_COMMAND.md) — choose a command for your analysis goal.
- [Workflows](WORKFLOWS.md) — see how commands fit together in cfDNA, MNase-seq, ChIP-seq/CUT&RUN/CUT&Tag, empirical peak FDR, peak spacing, categorized [`flank-spacing`](commands/flank-spacing.md), chromatin-state, and gene-centred analyses.
- [Command reference](COMMAND_REFERENCE.md) — links to every command page.
- [Version 0.10.1 release notes](RELEASE_NOTES_0.10.1.md) — fixes multicontig output discovery in `chip-suite`.
- [Version 0.10.0 release notes](RELEASE_NOTES_0.10.0.md) — summarizes the new workflows, interfaces, and compatibility changes.

## Understanding the calculations

- [Algorithms](ALGORITHMS.md) describes each calculation and gives its defining mathematics.
- [File formats](FILE_FORMATS.md) explains the coordinate and track formats NucleoSuite reads and writes.
- [Output layout](OUTPUT_LAYOUT.md) explains single-contig, multicontig, suite, and combined output directories.
- [Plot customization](PLOTTING.md) explains the shared plotting options.
- [Glossary](GLOSSARY.md) defines recurring NucleoSuite terms.

## Bundled resources

NucleoSuite includes reference files that can be used directly from the command line. List them with:

```bash
nucleosuite resources list
```

Print the installed path of a named resource with:

```bash
nucleosuite resources path gm12878-hg19-states
```

`resources path` prints a path that can be inserted into another command with shell command substitution:

```bash
nucleosuite distances sample_nucleosome_regions.bed \
  --position-column 7 \
  --state-bed "$(nucleosuite resources path gm12878-hg19-states)" \
  --output-prefix sample_distances
```

See [`nucleosuite resources`](commands/resources.md) for the complete bundled resource list and examples.

## CLI help

Command-line help is layered so routine options stay readable:

```bash
nucleosuite COMMAND --help       # core inputs and analysis controls
nucleosuite COMMAND --help-all   # all command-specific options and defaults
nucleosuite COMMAND --help-plotting  # shared plot controls, when available
```

- [`filter-peaks`](commands/filter-peaks.md) filters peak interval files by score, percentile, region length, and/or BigWig coverage.
