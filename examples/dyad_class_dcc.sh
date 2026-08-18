#!/usr/bin/env bash
set -euo pipefail

BAM=${1:?Usage: dyad_class_dcc.sh sample.bam genome.chrom.sizes [prefix]}
CHROM_SIZES=${2:?Usage: dyad_class_dcc.sh sample.bam genome.chrom.sizes [prefix]}
PREFIX=${3:-sample}

nucleosuite dyads -b "$BAM" -c autosomes \
  --frag-lower 145 --frag-upper 147 \
  --even-dyad split -o "${PREFIX}_145_147"

nucleosuite dyads -b "$BAM" -c autosomes \
  --frag-lower 166 --frag-upper 168 \
  --even-dyad split -o "${PREFIX}_166_168"

nucleosuite dcc bigwig \
  --bigwig-a "${PREFIX}_145_147_dyads_lower145_upper147_dyad.bw" \
  --bigwig-b "${PREFIX}_166_168_dyads_lower166_upper168_dyad.bw" \
  --chrom-sizes "$CHROM_SIZES" \
  --signed-lags \
  --dmax 50 \
  --label-a 145_147 \
  --label-b 166_168 \
  --out-prefix "${PREFIX}_145_147_vs_166_168"
