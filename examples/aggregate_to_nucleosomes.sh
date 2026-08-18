#!/usr/bin/env bash
set -euo pipefail

BIGWIG=${1:?Usage: aggregate_to_nucleosomes.sh signal.bw regions.bed nucleosomes.bed [output_dir]}
REGIONS=${2:?Usage: aggregate_to_nucleosomes.sh signal.bw regions.bed nucleosomes.bed [output_dir]}
NUCLEOSOMES=${3:?Usage: aggregate_to_nucleosomes.sh signal.bw regions.bed nucleosomes.bed [output_dir]}
OUTPUT_DIR=${4:-aggregated}

nucleosuite aggregate \
  --bigwig "$BIGWIG" \
  --region-bed "$REGIONS" \
  --nucleosome-bed "$NUCLEOSOMES" \
  --nucleosome-offset 1 \
  --window-half 3000 \
  --output-dir "$OUTPUT_DIR"
