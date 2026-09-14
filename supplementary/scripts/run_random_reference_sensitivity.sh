#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

for dataset in "${OLD7_DATASETS[@]}"; do
  run_r redesign/run_older7_random_sensitivity.R "$dataset" all 30 "$GENESELECTR_BUDGET"
done
