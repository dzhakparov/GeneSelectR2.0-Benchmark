#!/bin/zsh
set -euo pipefail

repo_root="${0:A:h:h}"
result_root="$repo_root/redesign/results_frozen_external_exact_2026-08-31"
log_root="$result_root/logs"
workers="${GENESELECTR_OUTER_WORKERS:-2}"
mkdir -p "$log_root"

export GENESELECTR_REDESIGN_RESULTS_ROOT="$result_root"
export GENESELECTR_N_CORES=1
export GENESELECTR_FIT_N_CORES=1
export GENESELECTR_VALIDATION_GS_ARMS=GS_full_ungrouped
export GENESELECTR_VALIDATION_COMP_ARMS=none
export R_USER_CACHE_DIR="$repo_root/data/r_user_cache"
export GENESELECTR_REPO_ROOT="$repo_root"
export GENESELECTR_LOG_ROOT="$log_root"

cd "$repo_root"

for dataset in GSE16879 GSE91061 GSE92415 GSE206285; do
  for repeat_idx in 1 2 3; do
    for fold_idx in 1 2 3 4 5; do
      print -- "$dataset r${repeat_idx}f${fold_idx}"
    done
  done
done | xargs -n 2 -P "$workers" sh -c '
  dataset="$1"
  split_key="$2"
  log_path="$GENESELECTR_LOG_ROOT/${dataset}_${split_key}_fit.log"
  cd "$GENESELECTR_REPO_ROOT"
  GENESELECTR_VALIDATION_SPLIT="$split_key" \
    Rscript redesign/run_validation_benchmark.R "$dataset" fit 604800 \
    >"$log_path" 2>&1
  status=$?
  printf "%s %s fit exit=%d\n" "$dataset" "$split_key" "$status"
  exit "$status"
' sh
