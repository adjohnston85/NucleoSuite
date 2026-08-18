#!/usr/bin/env bash
set -euo pipefail

BAM=${1:?Usage: dyad_to_dac.sh sample.bam genome.chrom.sizes [prefix]}
CHROM_SIZES=${2:?Usage: dyad_to_dac.sh sample.bam genome.chrom.sizes [prefix]}
PREFIX=${3:-sample_145_147}

nucleosuite dyads \
  -b "$BAM" \
  -c autosomes \
  --frag-lower 145 \
  --frag-upper 147 \
  --even-dyad split \
  --output-format bigwig \
  -o "$PREFIX"

DYAD_BW="${PREFIX}_dyads_lower145_upper147_dyad.bw"

nucleosuite dac \
  --bigwig "$DYAD_BW" \
  --chrom-sizes "$CHROM_SIZES" \
  --scope combined_chromosomes \
  --dmax 1500 \
  --out-prefix "${PREFIX}_dyad"
