#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

for dataset in "${OLD7_DATASETS[@]}"; do
  run_old7 redesign/run_older7_component_ablation.R "$dataset"
done
