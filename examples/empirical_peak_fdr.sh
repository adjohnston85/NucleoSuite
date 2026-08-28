#!/usr/bin/env bash
set -euo pipefail

OBSERVED=${1:?usage: empirical_peak_fdr.sh OBSERVED.bed RANDOMIZED.bed OUTPUT_PREFIX [FDR]}
RANDOMIZED=${2:?usage: empirical_peak_fdr.sh OBSERVED.bed RANDOMIZED.bed OUTPUT_PREFIX [FDR]}
OUTPUT_PREFIX=${3:?usage: empirical_peak_fdr.sh OBSERVED.bed RANDOMIZED.bed OUTPUT_PREFIX [FDR]}
FDR=${4:-}

arguments=(
  "$OBSERVED"
  "$RANDOMIZED"
  --output-prefix "$OUTPUT_PREFIX"
)
if [[ -n "$FDR" ]]; then
  arguments+=(--fdr "$FDR")
fi

nucleosuite empirical-peak-fdr "${arguments[@]}"
