# Installation

## Install from a cloned repository

Create a Conda or Mamba environment for the Python packages and command-line bioinformatics tools, then install NucleoSuite in editable mode.

From the repository root, where `pyproject.toml` is located:

```bash
mamba env create -f environment.yml
conda activate nucleosuite
python -m pip install -e . --no-deps
```

### What `pip install -e .` means

- `.` means “install the project in the current directory”.
- `-e` means editable.
- The installed `nucleosuite` command points to the source code in this directory.
- Python code changes are therefore visible immediately without rebuilding and reinstalling a wheel.
- Packaging metadata changes, entry-point changes or dependency changes may still require reinstalling.

Check the installation:

```bash
which nucleosuite
nucleosuite --version
nucleosuite --help
```

The environment file includes every Python and command-line dependency used by
NucleoSuite, including `matplotlib-venn` for gene-set Venn diagrams, `openpyxl`
for spreadsheet output, `pysam`, `pyBigWig`, `samtools`, and the UCSC BigWig/
bigBed conversion tools.

To bring an existing environment into line with the current file:

```bash
mamba env update -n nucleosuite -f environment.yml --prune
conda activate nucleosuite
python -m pip install -e . --no-deps
```

## Install from a wheel

A wheel is a built Python package file:

```bash
python -m pip install nucleosuite-0.8.6-py3-none-any.whl
```

To reinstall the wheel explicitly:

```bash
python -m pip install --upgrade --force-reinstall --no-cache-dir \
  nucleosuite-0.8.6-py3-none-any.whl
hash -r
```

`hash -r` clears the shell’s remembered command location.

## Install from the source directory without editable mode

```bash
python -m pip install .
```

This installs a copy of the package. Reinstall after changing the source.

## Build a wheel and source archive

Install the build tool in the active environment:

```bash
python -m pip install build
```

Build:

```bash
python -m build
```

Expected outputs include:

```text
dist/nucleosuite-0.8.6-py3-none-any.whl
dist/nucleosuite-0.8.6.tar.gz
```

## Build a Conda package

The repository includes `recipe/meta.yaml`.

```bash
mamba install conda-build
conda build recipe
```

The recipe installs NucleoSuite and its run-time dependencies. Use `conda build purge` when a cached build needs to be cleared.

## Required input preparation

### BAM

Most BAM-based commands expect coordinate-sorted paired-end data and an index:

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

### FASTA

Sequence-aware commands need a FASTA index:

```bash
samtools faidx genome.fa
```

### Chromosome sizes

Create a two-column chromosome-size file from the FASTA index:

```bash
cut -f1,2 genome.fa.fai > genome.chrom.sizes
```

## WSL and Windows paths

A Windows `C:` path appears under WSL as `/mnt/c`:

```text
C:\NucleoSuite\NucleoSuite
```

becomes:

```text
/mnt/c/NucleoSuite/NucleoSuite
```

Run Linux bioinformatics tools from WSL-style paths. Large analyses can be faster and more reliable in the WSL Linux filesystem than directly under `/mnt/c`, especially when many temporary files are created.

[Back to the documentation index](README.md)
