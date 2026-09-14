#!/usr/bin/env bash
set -euo pipefail

# Configuration block. These values are inherited by every supplementary entry point.
SUPPLEMENTARY_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd "$SUPPLEMENTARY_SCRIPT_DIR/../.." && pwd)"
GENESELECTR_WORKERS="${GENESELECTR_N_CORES:-2}"
GENESELECTR_BUDGET="${GENESELECTR_CALL_BUDGET:-604800}"
R_BIN="${R_BIN:-$(command -v Rscript)}"

export GENESELECTR_N_CORES="$GENESELECTR_WORKERS"
export GENESELECTR_CALL_BUDGET="$GENESELECTR_BUDGET"
export GENESELECTR_REDESIGN_RESULTS_ROOT="${GENESELECTR_REDESIGN_RESULTS_ROOT:-redesign/results_corrected}"

cd "$REPOSITORY_ROOT"

run_r() {
  "$R_BIN" "$@"
}

run_old7() {
  local script="$1"
  local dataset="$2"
  run_r "$script" "$dataset" all "$GENESELECTR_BUDGET"
}

OLD7_DATASETS=(GSE101794 GSE107994 GSE13355 GSE65682 GSE69683 imvigor210 sosall)
PRIMARY_DATASETS=(GSE101794 GSE107994 GSE13355 GSE65682 GSE69683 imvigor210 sosall)
