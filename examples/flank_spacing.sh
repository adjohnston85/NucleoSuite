#!/usr/bin/env bash
set -euo pipefail

# Compare nucleosome spacing around categorized reference sites.
nucleosuite flank-spacing \
  --nucleosome-bed sample_nucleosome_regions.bed \
  --region-bed categorized_sites.bed \
  --category-col 4 \
  --output-prefix sample_flank_spacing
