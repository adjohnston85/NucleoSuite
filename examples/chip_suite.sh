#!/usr/bin/env bash
set -euo pipefail

TARGET_BAM=${1:?usage: chip_suite.sh TARGET.bam CONTROL.bam OUTDIR [MODE]}
CONTROL_BAM=${2:?usage: chip_suite.sh TARGET.bam CONTROL.bam OUTDIR [MODE]}
OUTDIR=${3:?usage: chip_suite.sh TARGET.bam CONTROL.bam OUTDIR [MODE]}
MODE=${4:-auto}

nucleosuite chip-suite \
  --treatment1-bam "$TARGET_BAM" \
  --control1-bam "$CONTROL_BAM" \
  --outdir "$OUTDIR" \
  --scoring-method tns \
  --mode "$MODE" \
  --frag-lower 120 \
  --frag-upper 500 \
  --peak-fdr 0.05 \
  --cluster-fdr 0.05
