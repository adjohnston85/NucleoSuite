#!/usr/bin/env bash
set -euo pipefail

OUT_PREFIX=${1:?Usage: fragment_length_heatmap.sh output_prefix table1.tsv [table2.tsv ...]}
shift

ARGS=()
for table in "$@"; do
  ARGS+=( -i "$table" )
done

nucleosuite fragment-heatmap \
  "${ARGS[@]}" \
  --normalization fragment-zscore \
  --min-frag 100 \
  --max-frag 250 \
  -o "$OUT_PREFIX"
