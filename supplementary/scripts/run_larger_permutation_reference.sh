#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

for dataset in GSE107994 GSE13355; do
  run_r redesign/run_older7_null100_confirmation.R "$dataset" "$GENESELECTR_BUDGET"
done
