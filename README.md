# NucleoSuite

NucleoSuite is a command-line toolkit that converts paired-end alignments or fragment intervals into nucleosome-positioning signals, cfDNA and MNase-seq fragmentomic analyses, peak calls, chromatin-state profiles, sequence profiles, spacing measurements, and coordinated workflows.

```bash
nucleosuite COMMAND [options]
```

Verify an installation with:

```bash
nucleosuite --version
nucleosuite --help
```

## Commands

### Fragment preparation

| Command | Function |
|---|---|
| [`fragments`](docs/commands/fragments.md) | Convert paired-end alignments to BED fragment intervals or combine fragment files. |
| [`merge-bams`](docs/commands/merge-bams.md) | Merge BAM files while retaining alignment records and tags. |
| [`randomize-fragments`](docs/commands/randomize-fragments.md) | Generate coordinate-randomized control fragments. |
| [`fragment-lengths`](docs/commands/fragment-lengths.md) | Calculate fragment-length counts and percentages. |
| [`fragment-heatmap`](docs/commands/fragment-heatmap.md) | Compare fragment-length profiles across samples or region classes. |

### Signal and sequence profiles

| Command | Function |
|---|---|
| [`tracks`](docs/commands/tracks.md) | Generate multiple range-specific signal and coordinate tracks in one fragment pass. |
| [`pns`](docs/commands/pns.md) | Calculate PNS, BNS or TNS nucleosome score tracks and shared peak calls. |
| [`wps`](docs/commands/wps.md) | Calculate window protection score tracks and WPS peak calls. |
| [`coverage`](docs/commands/coverage.md) | Calculate per-base fragment coverage. |
| [`mean-scale`](docs/commands/mean-scale.md) | Scale a BigWig relative to a supplied, region-derived, or non-zero-signal reference mean. |
| [`dyads`](docs/commands/dyads.md) | Generate fragment-centre tracks. |
| [`fragment-ends`](docs/commands/fragment-ends.md) | Generate combined, left-end, and right-end tracks. |
| [`dinuc-profile`](docs/commands/dinuc-profile.md) | Calculate positional dinucleotide profiles. |
| [`ww-types`](docs/commands/ww-types.md) | Classify fragments by centred WW/SS sequence patterns. |

### Peaks, spacing, and comparisons

| Command | Function |
|---|---|
| [`call-peaks`](docs/commands/call-peaks.md) | Call nucleosome and breakpoint features from PNS or WPS BigWigs. |
| [`distances`](docs/commands/distances.md) | Calculate adjacent and higher-order distances between called positions. |
| [`flank-spacing`](docs/commands/flank-spacing.md) | Compare the spacing between nucleosomes flanking categorized reference sites and rank category-specific distributions. |
| [`dac`](docs/commands/dac.md) | Calculate distance autocorrelation within one signal. |
| [`dcc`](docs/commands/dcc.md) | Calculate distance cross-correlation between two signals. |
| [`nrl`](docs/commands/nrl.md) | Estimate nucleosome repeat length from recurring DAC or DCC peaks. |
| [`compare-positions`](docs/commands/compare-positions.md) | Compare one main callset with one or more positional callsets and optional BigWig score comparators. |
| [`positive-runs`](docs/commands/positive-runs.md) | Measure contiguous positive-signal intervals in a BigWig. |
| [`peak-score-frequency`](docs/commands/peak-score-frequency.md) | Compare peak-score distributions. |
| [`filter-coverage`](docs/commands/filter-coverage.md) | Filter BED peaks by BigWig coverage at the summit or interval midpoint. |
| [`filter-peaks`](docs/commands/filter-peaks.md) | Filter peak intervals by score, score percentile, or region length. |
| [`peak-states`](docs/commands/peak-states.md) | Measure peak abundance and score-dependent enrichment by chromatin state. |

### Regional and gene analyses

| Command | Function |
|---|---|
| [`aggregate`](docs/commands/aggregate.md) | Aggregate BigWig signal around genomic reference features. |
| [`region-extract`](docs/commands/region-extract.md) | Export region-level signal vectors and nearby peaks. |
| [`gene-sets`](docs/commands/gene-sets.md) | Define gene groups from chromatin-state overlaps. |
| [`gene-expression`](docs/commands/gene-expression.md) | Relate expression to peak spacing or signal periodicity. |
| [`tss-expression-quintiles`](docs/commands/tss-expression-quintiles.md) | Aggregate PNS/WPS around TSSs after splitting genes into tissue-expression quintiles. |

### Workflows and utilities

| Command | Function |
|---|---|
| [`mnase-suite`](docs/commands/mnase-suite.md) | Run the coordinated MNase-seq workflow. |
| [`cfdna-suite`](docs/commands/cfdna-suite.md) | Run the coordinated cfDNA fragmentomics workflow. |
| [`combine`](docs/commands/combine.md) | Combine outputs from an existing chromosome-wise run. |
| [`chrom-sizes`](docs/commands/chrom-sizes.md) | Write chromosome names and lengths from a BAM or CRAM header. |
| [`resources`](docs/commands/resources.md) | List, locate, validate, or copy bundled resources. |
| [`validate-inputs`](docs/commands/validate-inputs.md) | Validate input integrity and reference compatibility before a run. |
| [`plot`](docs/commands/plot.md) | Recreate and deeply customize all applicable figures from existing NucleoSuite output tables. |


## Typical workflows

For a single analysis, run the command that produces the required signal, profile, peak set, or spacing result. For coordinated analyses, use [`mnase-suite`](docs/commands/mnase-suite.md) or [`cfdna-suite`](docs/commands/cfdna-suite.md). Existing chromosome-wise runs can be recombined with [`combine`](docs/commands/combine.md), and generated figures can be recreated or customized with [`plot`](docs/commands/plot.md).

A minimal PNS example:

```bash
nucleosuite pns \
  --bam sample.bam \
  --fasta genome.fa \
  --contigs chr1 chr2 chr3 chr4 \
  --cores 4 \
  --out-prefix sample
```

Detailed command behaviour, advanced options, output layouts, resource handling, and workflow examples are documented in the pages linked above and in the guides below.

## Installation

### Recommended installation

NucleoSuite is intended to run in the supplied Conda environment. The recommended installation is to clone the repository, create that environment, build the package locally, and install the generated wheel.

Before starting, ensure that [Git](docs/INSTALLATION.md#installing-git) and [Conda or Mamba](docs/INSTALLATION.md#installing-conda-and-mamba) are available.

1. Choose a directory for the repository and clone NucleoSuite:

```bash
mkdir -p ~/software
cd ~/software
git clone https://github.com/adjohnston85/NucleoSuite.git
cd NucleoSuite
```

2. Create and activate the supplied environment. Mamba is recommended for environment creation:

```bash
mamba env create -f environment.yml
conda activate nucleosuite
```

If Mamba is not available, use Conda instead:

```bash
conda env create -f environment.yml
conda activate nucleosuite
```

3. Build the wheel and source distribution:

```bash
rm -rf dist/
python -m build
```

4. Install the wheel into the active environment:

```bash
python -m pip install --upgrade --force-reinstall --no-deps dist/*.whl
```

The supplied environment already provides NucleoSuite's dependencies, so `--no-deps` prevents pip from replacing Conda-managed packages.

5. Verify the installation:

```bash
nucleosuite --version
nucleosuite --help
```

The `dist/` directory is generated locally and is intentionally excluded from Git.

Other common installation approaches are also supported. See [Alternative installation methods](docs/INSTALLATION.md#alternative-installation-methods), including direct source installation, editable development installs, source distributions, and installation of an already-built wheel.

For complete setup, prerequisite, update, and troubleshooting instructions, see [Installation](docs/INSTALLATION.md).

## Documentation

- [Documentation index](docs/README.md)
- [Quick start](docs/QUICKSTART.md)
- [Choosing a command](docs/CHOOSING_A_COMMAND.md)
- [Workflows](docs/WORKFLOWS.md)
- [Command reference](docs/COMMAND_REFERENCE.md)
- [File formats](docs/FILE_FORMATS.md)
- [Glossary](docs/GLOSSARY.md)
- [Algorithms](docs/ALGORITHMS.md)
- [Plot customization](docs/PLOTTING.md)

The command-line help is the authoritative reference for accepted options and defaults:

```bash
nucleosuite COMMAND --help
nucleosuite COMMAND --help-all
nucleosuite COMMAND --help-plotting
```
