#!/usr/bin/env bash
set -euo pipefail

BAM=${1:?Usage: pns_full.sh sample.bam [prefix]}
PREFIX=${2:-sample}

nucleosuite pns \
  -b "$BAM" \
  -c autosomes \
  --mode-length 167 \
  --frag-lower 137 \
  --frag-upper 197 \
  --max-duplicates 1 \
  --dedup-scope all_bams \
  --pns-format bigwig \
  --other-format bigwig \
  -o "$PREFIX"
