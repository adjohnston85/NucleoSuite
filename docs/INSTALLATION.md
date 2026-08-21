# Installation

## Recommended installation

NucleoSuite is intended to run in the Conda environment supplied with the repository. The recommended installation is to clone the repository, create that environment, build the Python package locally, and install the generated wheel.

### Prerequisites

The recommended installation requires:

- [Git](#installing-git) to clone and update the repository.
- [Conda or Mamba](#installing-conda-and-mamba) to create the NucleoSuite software environment.

The supplied `environment.yml` includes the Python build package (`python-build`), so no separate installation of `build` is required after the environment has been created.

### 1. Clone NucleoSuite

Choose a directory in which to keep the repository. For example:

```bash
mkdir -p ~/software
cd ~/software
```

Clone the repository and enter it:

```bash
git clone https://github.com/adjohnston85/NucleoSuite.git
cd NucleoSuite
```

### 2. Create the NucleoSuite environment

Mamba is recommended for environment creation:

```bash
mamba env create -f environment.yml
conda activate nucleosuite
```

If Mamba is not available, Conda can be used instead:

```bash
conda env create -f environment.yml
conda activate nucleosuite
```

The environment provides the Python and command-line dependencies used by NucleoSuite, including `pysam`, `pyBigWig`, `matplotlib`, `matplotlib-venn`, `openpyxl`, `samtools`, and the UCSC BigWig/bigBed conversion tools.

### 3. Build NucleoSuite

From the repository root, remove any distributions left from an older build and build the current source tree:

```bash
rm -rf dist/
python -m build
```

This creates a wheel and source distribution under `dist/`, for example:

```text
dist/
├── nucleosuite-<version>-py3-none-any.whl
└── nucleosuite-<version>.tar.gz
```

The `dist/` directory is generated locally and is intentionally excluded from Git.

### 4. Install the wheel

Install the wheel into the active `nucleosuite` environment:

```bash
python -m pip install --upgrade --force-reinstall --no-deps dist/*.whl
```

`--no-deps` is used because the supplied Conda environment already provides the package dependencies. This avoids pip replacing Conda-managed scientific or bioinformatics packages.

### 5. Verify the installation

```bash
nucleosuite --version
nucleosuite --help
```

The Python package version can also be checked directly:

```bash
python -c "import nucleosuite; print(nucleosuite.__version__)"
```

If the shell has remembered an older executable location after replacing an existing installation, run:

```bash
hash -r
```

and check again.

## Alternative installation methods

The wheel-based installation above is recommended for normal use. The following methods are useful in other common Python workflows.

### Install directly from the source tree

After creating and activating the supplied environment, NucleoSuite can be installed without first building a wheel:

```bash
python -m pip install --no-deps .
```

This installs a copy of the package from the current source tree. Reinstall after changing the source.

### Editable installation for development

For development work, install the source tree in editable mode:

```bash
python -m pip install --no-deps -e .
```

The installed command then points to the repository source, so Python code changes are available without rebuilding and reinstalling the wheel after every edit. Packaging metadata, entry-point, or dependency changes may still require reinstalling.

### Install an already-built wheel

If a compatible NucleoSuite wheel has already been built or downloaded, install it into an activated NucleoSuite environment:

```bash
python -m pip install --upgrade --force-reinstall --no-deps \
  /path/to/nucleosuite-<version>-py3-none-any.whl
```

### Install the source distribution

After `python -m build`, the generated source distribution can also be installed:

```bash
python -m pip install --no-deps dist/nucleosuite-*.tar.gz
```

The wheel is preferred when both distributions are available because it installs the already-built Python package.

## Updating an existing installation

Move to the existing repository and pull the latest source:

```bash
cd /path/to/NucleoSuite
git pull
```

Activate the environment:

```bash
conda activate nucleosuite
```

Update the environment in case `environment.yml` has changed:

```bash
mamba env update -n nucleosuite -f environment.yml --prune
```

or with Conda:

```bash
conda env update -n nucleosuite -f environment.yml --prune
```

Rebuild and reinstall the wheel:

```bash
rm -rf dist/
python -m build
python -m pip install --upgrade --force-reinstall --no-deps dist/*.whl
hash -r
```

Then verify the installed version:

```bash
nucleosuite --version
```

## Installing Git

Git is required to clone and update the NucleoSuite repository.

### Ubuntu, Debian, and WSL2

```bash
sudo apt update
sudo apt install git
```

Verify the installation:

```bash
git --version
```

For other platforms, use the installation instructions from the official Git project:

- https://git-scm.com/downloads

Return to [Recommended installation](#recommended-installation).

## Installing Conda and Mamba

NucleoSuite uses Conda environments to provide both Python packages and external bioinformatics tools.

### Quick Conda installation on Linux or WSL2

If Conda is not already installed, Miniconda provides a small Conda installation. For a standard x86-64 Linux or WSL2 system:

```bash
curl -L -o Miniconda3-latest-Linux-x86_64.sh \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the installer prompts and allow Conda to initialise the shell. Close and reopen the terminal when installation is complete, then verify Conda:

```bash
conda --version
```

Official Conda installation instructions, including installers for other operating systems and architectures, are available at:

- https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

### Install Mamba into an existing Conda installation

Install Mamba into the Conda `base` environment from `conda-forge`:

```bash
conda install -n base -c conda-forge --override-channels mamba
```

Verify it:

```bash
mamba --version
```

Mamba can then be used for the NucleoSuite environment creation and update commands. The NucleoSuite `environment.yml` explicitly uses `conda-forge`, `bioconda`, and `nodefaults`.

Mamba documentation is available at:

- https://mamba.readthedocs.io/

If starting from scratch, Miniforge is another common option and provides a Conda installation configured for `conda-forge`; see the Conda installation documentation above.

Return to [Recommended installation](#recommended-installation).

## Build a Conda package

The repository includes `recipe/meta.yaml` for users who want to build a Conda package rather than use the recommended wheel installation.

Install `conda-build` if it is not already available:

```bash
mamba install conda-build
```

Then build the recipe from the repository root:

```bash
conda build recipe
```

Use `conda build purge` when a cached build needs to be cleared.

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
