#!/usr/bin/env bash
set -euo pipefail

OBSERVED=${1:?usage: pns_peak_fdr.sh OBSERVED.bed RANDOMIZED.bed OUTPUT_PREFIX [FDR]}
RANDOMIZED=${2:?usage: pns_peak_fdr.sh OBSERVED.bed RANDOMIZED.bed OUTPUT_PREFIX [FDR]}
OUTPUT_PREFIX=${3:?usage: pns_peak_fdr.sh OBSERVED.bed RANDOMIZED.bed OUTPUT_PREFIX [FDR]}
FDR=${4:-}

arguments=(
  "$OBSERVED"
  "$RANDOMIZED"
  --output-prefix "$OUTPUT_PREFIX"
)
if [[ -n "$FDR" ]]; then
  arguments+=(--fdr "$FDR")
fi

nucleosuite pns-peak-fdr "${arguments[@]}"
