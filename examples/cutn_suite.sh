#!/usr/bin/env bash
set -euo pipefail

TARGET_BAM=${1:?usage: cutn_suite.sh TARGET.bam CONTROL.bam OUTDIR [MODE]}
CONTROL_BAM=${2:?usage: cutn_suite.sh TARGET.bam CONTROL.bam OUTDIR [MODE]}
OUTDIR=${3:?usage: cutn_suite.sh TARGET.bam CONTROL.bam OUTDIR [MODE]}
MODE=${4:-auto}

nucleosuite cutn-suite \
  --treatment1-bam "$TARGET_BAM" \
  --control1-bam "$CONTROL_BAM" \
  --outdir "$OUTDIR" \
  --mode "$MODE" \
  --frag-lower 120 \
  --frag-upper 500
