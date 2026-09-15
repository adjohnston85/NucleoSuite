#!/usr/bin/env bash
set -euo pipefail

nucleosuite flank-filter \
  --regions sample_nucleosome_regions.bed \
  --flanks sample_breakpoint_peaks.bed \
  --out sample_nucleosome_regions_flanked.bed
