#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

run_r redesign/fetch_ot_seeds_dense.R
for dataset in "${OLD7_DATASETS[@]}"; do
  run_old7 redesign/run_older7_biology.R "$dataset"
done
run_r redesign/run_older7_biology_ot_dense.R
