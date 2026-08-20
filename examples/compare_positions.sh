#!/usr/bin/env bash
set -euo pipefail

nucleosuite compare-positions \
  --main-bed PNS=sample_pns_nucleosomes.bed \
  --compare-bed iNPS=sample_inps_nucleosomes.bed \
  --compare-bed DANPOS=sample_danpos_nucleosomes.bed \
  --main-score-column 5 \
  --compare-score-column 5 \
  --stats \
  --output-prefix sample_position_compare
