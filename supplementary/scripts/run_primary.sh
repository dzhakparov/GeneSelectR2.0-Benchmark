#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

# Reported workflow: base seven-dataset benchmark plus saved-ranking QA.
bash redesign/run_corrected_benchmark.sh base
for dataset in "${PRIMARY_DATASETS[@]}"; do
  run_r redesign/reevaluate_saved_rankings.R "$dataset"
done
run_r redesign/validate_deterministic_results.R
