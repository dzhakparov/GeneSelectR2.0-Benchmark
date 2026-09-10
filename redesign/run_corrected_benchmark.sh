#!/usr/bin/env bash
# Run the corrected redesign benchmark from the project root.
#
# Every R runner retains its own provenance and checkpoint checks. This script
# supplies dependency order and verifies the expected number of completed rows
# before advancing. Two workers are used by default so other tasks retain CPU
# and memory capacity.

set -euo pipefail
cd "$(dirname "$0")/.."

scope="${1:-all}"
export GENESELECTR_N_CORES="${GENESELECTR_N_CORES:-2}"
export GENESELECTR_REDESIGN_RESULTS_ROOT="${GENESELECTR_REDESIGN_RESULTS_ROOT:-redesign/results_corrected}"

R_BIN="${R_BIN:-$(command -v Rscript)}"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3)}"
CALL_BUDGET="${GENESELECTR_CALL_BUDGET:-7000}"
RESULTS_ROOT="$GENESELECTR_REDESIGN_RESULTS_ROOT"
LOG_DIR="$RESULTS_ROOT/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/corrected_benchmark.log"

validation_datasets=(GSE107994 GSE13355 GSE101794 GSE69683 GSE65682)
all_datasets=(GSE107994 GSE13355 GSE101794 GSE69683 GSE65682 imvigor210 sosall)
cohorts=(imvigor210 sosall)

log_message() {
  printf '%s %s\n' "$(date '+%F %H:%M:%S')" "$*" | tee -a "$LOG_FILE"
}

run_r() {
  log_message "Rscript $*"
  "$R_BIN" "$@" 2>&1 | tee -a "$LOG_FILE"
}

dataset_dir() {
  case "$1" in
    imvigor210|sosall) printf '%s/full_recipe/%s\n' "$RESULTS_ROOT" "$1" ;;
    *) printf '%s/validation_benchmark/%s\n' "$RESULTS_ROOT" "$1" ;;
  esac
}

csv_rows() {
  if [[ -f "$1" ]]; then
    awk 'END {if (NR > 0) print NR - 1; else print 0}' "$1"
  else
    printf '0\n'
  fi
}

file_count() {
  local directory="$1"
  local pattern="$2"
  if [[ -d "$directory" ]]; then
    find "$directory" -maxdepth 1 -type f -name "$pattern" | wc -l | tr -d ' '
  else
    printf '0\n'
  fi
}

ensure_csv_rows() {
  local output="$1"
  local expected="$2"
  shift 2
  local attempts=0
  while [[ "$(csv_rows "$output")" -lt "$expected" ]]; do
    attempts=$((attempts + 1))
    if [[ "$attempts" -gt 300 ]]; then
      log_message "ERROR incomplete output after 300 passes: $output"
      return 1
    fi
    run_r "$@"
  done
  log_message "complete rows=$expected file=$output"
}

run_cohorts() {
  local cohort grouped_dir full_dir ranking_count
  for cohort in "${cohorts[@]}"; do
    grouped_dir="$RESULTS_ROOT/grouped_benchmark/$cohort"
    while [[ "$(file_count "$grouped_dir" 'split_r*_f*.rds')" -lt 15 ]]; do
      run_r redesign/run_grouped_benchmark.R "$cohort" prep "$CALL_BUDGET"
    done

    full_dir="$RESULTS_ROOT/full_recipe/$cohort"
    while true; do
      if [[ -d "$full_dir" ]]; then
        ranking_count=$(find "$full_dir" -maxdepth 1 -type f \
          -name 'ranking_r*_f*_full_*.csv' ! -name '*_meta.csv' | wc -l | tr -d ' ')
      else
        ranking_count=0
      fi
      [[ "$ranking_count" -ge 30 ]] && break
      run_r redesign/run_full_recipe.R "$cohort" fit "$CALL_BUDGET"
    done
    ensure_csv_rows "$full_dir/eval_results.csv" 270 \
      redesign/run_full_recipe.R "$cohort" eval "$CALL_BUDGET"
    run_r redesign/run_full_recipe.R "$cohort" report "$CALL_BUDGET"
  done
}

run_validation() {
  local dataset out_dir ranking_count
  for dataset in "${validation_datasets[@]}"; do
    out_dir="$RESULTS_ROOT/validation_benchmark/$dataset"
    while [[ "$(file_count "$out_dir" 'split_r*_f*.rds')" -lt 15 ]]; do
      run_r redesign/run_validation_benchmark.R "$dataset" prep "$CALL_BUDGET"
    done
    while true; do
      ranking_count=$(find "$out_dir" -maxdepth 1 -type f \
        -name 'ranking_r*_f*.csv' ! -name '*_meta.csv' | wc -l | tr -d ' ')
      [[ "$ranking_count" -ge 120 ]] && break
      run_r redesign/run_validation_benchmark.R "$dataset" fit "$CALL_BUDGET"
    done
    ensure_csv_rows "$out_dir/eval_results.csv" 810 \
      redesign/run_validation_benchmark.R "$dataset" eval "$CALL_BUDGET"
  done
}

run_predfirst() {
  local dataset expected out_dir
  for dataset in "${all_datasets[@]}"; do
    out_dir="$(dataset_dir "$dataset")"
    expected=180
    [[ "$dataset" == "imvigor210" || "$dataset" == "sosall" ]] && expected=720
    ensure_csv_rows "$out_dir/eval_predfirst.csv" "$expected" \
      redesign/run_predfirst_benchmark.R "$dataset" all "$CALL_BUDGET"
  done
}

run_extension() {
  local script="$1"
  local filename="$2"
  local expected="$3"
  local dataset out_dir
  for dataset in "${all_datasets[@]}"; do
    out_dir="$(dataset_dir "$dataset")"
    ensure_csv_rows "$out_dir/$filename" "$expected" \
      "$script" "$dataset" "$CALL_BUDGET"
  done
}

run_biology() {
  local dataset out_dir
  for dataset in "${all_datasets[@]}"; do
    out_dir="$(dataset_dir "$dataset")"
    ensure_csv_rows "$out_dir/biology_string.csv" 420 \
      redesign/run_biology_benchmark.R "$dataset" "$CALL_BUDGET"
    ensure_csv_rows "$out_dir/biology_disease.csv" 420 \
      redesign/run_biology_disease_benchmark.R "$dataset" "$CALL_BUDGET"
  done
}

run_report() {
  log_message "building corrected report"
  "$PYTHON_BIN" redesign/build_full_report_v10.py 2>&1 | tee -a "$LOG_FILE"
  "$PYTHON_BIN" redesign/build_comprehensive_report.py 2>&1 | tee -a "$LOG_FILE"
}

case "$scope" in
  cohorts) run_cohorts ;;
  validation) run_validation ;;
  base) run_cohorts; run_validation ;;
  predfirst) run_predfirst ;;
  ensemble) run_extension redesign/run_ensemble_eval.R eval_ensemble.csv 270 ;;
  softprior) run_extension redesign/run_softprior_benchmark.R eval_softprior.csv 180 ;;
  datadriven) run_extension redesign/run_datadriven_benchmark.R eval_datadriven.csv 360 ;;
  horseshoe) run_extension redesign/run_horseshoe_benchmark.R eval_horseshoe.csv 90 ;;
  ensemble2) run_extension redesign/run_ensemble2_eval.R eval_ensemble2.csv 180 ;;
  wsweep) run_extension redesign/run_wsweep.R eval_wsweep.csv 90 ;;
  hsstab) run_extension redesign/run_hsstab_benchmark.R eval_hsstab.csv 90 ;;
  prune) run_extension redesign/run_prune_benchmark.R eval_prune.csv 270 ;;
  adapt) run_extension redesign/run_adapt_benchmark.R eval_adapt.csv 90 ;;
  biology) run_biology ;;
  report) run_report ;;
  all)
    run_cohorts
    run_validation
    run_predfirst
    run_extension redesign/run_ensemble_eval.R eval_ensemble.csv 270
    run_extension redesign/run_softprior_benchmark.R eval_softprior.csv 180
    run_extension redesign/run_datadriven_benchmark.R eval_datadriven.csv 360
    run_extension redesign/run_horseshoe_benchmark.R eval_horseshoe.csv 90
    run_extension redesign/run_ensemble2_eval.R eval_ensemble2.csv 180
    run_extension redesign/run_wsweep.R eval_wsweep.csv 90
    run_extension redesign/run_hsstab_benchmark.R eval_hsstab.csv 90
    run_extension redesign/run_prune_benchmark.R eval_prune.csv 270
    run_extension redesign/run_adapt_benchmark.R eval_adapt.csv 90
    run_biology
    run_report
    ;;
  *)
    printf 'Unknown scope: %s\n' "$scope" >&2
    exit 2
    ;;
esac

log_message "scope complete: $scope"
