#!/usr/bin/env bash
set -euo pipefail

# Retain only region calls whose directly surrounding event types are
# FLANK -- REGION -- FLANK. Consecutive region calls between the same flanks
# are rejected.
nucleosuite flank-filter \
  --regions sample_nucleosome_regions.bed \
  --flanks sample_breakpoint_peaks.bed \
  --out sample_nucleosome_regions_flanked.bed
