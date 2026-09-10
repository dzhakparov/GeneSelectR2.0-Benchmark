#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

scope="${1:-all}"
export GENESELECTR_N_CORES="${GENESELECTR_N_CORES:-2}"
export GENESELECTR_REDESIGN_RESULTS_ROOT="${GENESELECTR_REDESIGN_RESULTS_ROOT:-redesign/results_corrected}"
R_BIN="${R_BIN:-$(command -v Rscript)}"
CALL_BUDGET="${GENESELECTR_CALL_BUDGET:-604800}"

datasets=(GSE101794 GSE107994 GSE13355 GSE65682 GSE69683 imvigor210 sosall)
exploratory_datasets=(GSE16879 GSE91061 GSE92415 GSE206285)

run_tests() {
  bash analysis/check_repository_contents.sh
  "$R_BIN" analysis/check_inputs.R
  "$R_BIN" redesign/run_tests.R
}

run_primary() {
  local dataset
  bash redesign/run_corrected_benchmark.sh all
  for dataset in "${datasets[@]}"; do
    "$R_BIN" redesign/reevaluate_saved_rankings.R "$dataset"
  done
  "$R_BIN" redesign/validate_deterministic_results.R
}

run_ablations() {
  local dataset
  for dataset in "${datasets[@]}"; do
    "$R_BIN" redesign/run_older7_component_ablation.R "$dataset" all "$CALL_BUDGET"
    "$R_BIN" redesign/run_older7_calibration_diagnostics.R "$dataset" all "$CALL_BUDGET"
    "$R_BIN" redesign/run_older7_random_sensitivity.R "$dataset" all 30 "$CALL_BUDGET"
    "$R_BIN" redesign/run_older7_adaptive.R "$dataset" all "$CALL_BUDGET"
  done
  "$R_BIN" redesign/run_older7_calibration_onoff.R all "$CALL_BUDGET"
  "$R_BIN" redesign/run_older7_null100_confirmation.R GSE107994 "$CALL_BUDGET"
  "$R_BIN" redesign/run_older7_null100_confirmation.R GSE13355 "$CALL_BUDGET"
  "$R_BIN" redesign/run_older7_glmnet_convergence_check.R all
  "$R_BIN" redesign/summarise_older7_remaining_analyses.R
  "$R_BIN" redesign/summarise_older7_null_instability.R all
  "$R_BIN" redesign/summarise_older7_efficiency_biology.R
  "$R_BIN" redesign/summarise_older7_inference.R
}

run_biology() {
  local dataset
  "$R_BIN" redesign/fetch_ot_seeds_dense.R
  for dataset in "${datasets[@]}"; do
    "$R_BIN" redesign/run_older7_biology.R "$dataset" all "$CALL_BUDGET"
  done
  "$R_BIN" redesign/run_older7_biology_ot_dense.R
}

run_complementarity() {
  local dataset
  "$R_BIN" redesign/complementarity_analysis/01_gene_set_stability.R
  "$R_BIN" redesign/complementarity_analysis/02_dge_geneselectr_concordance.R
  "$R_BIN" redesign/complementarity_analysis/03_gene_group_characterization.R
  "$R_BIN" redesign/complementarity_analysis/04_gene_group_biology.R
  for dataset in "${datasets[@]}"; do
    "$R_BIN" redesign/complementarity_analysis/05_exploratory_dge_prefilter.R \
      run "$dataset"
  done
  "$R_BIN" redesign/complementarity_analysis/05_exploratory_dge_prefilter.R \
    assemble
  "$R_BIN" redesign/complementarity_analysis/06_make_figures.R
}

run_exploratory() (
  local dataset
  local result_root="redesign/results_frozen_external_exact_2026-08-31"
  export GENESELECTR_REDESIGN_RESULTS_ROOT="$result_root"
  export GENESELECTR_VALIDATION_GS_ARMS="GS_full_ungrouped"
  export GENESELECTR_VALIDATION_COMP_ARMS="none"
  for dataset in "${exploratory_datasets[@]}"; do
    "$R_BIN" redesign/run_validation_benchmark.R "$dataset" prep \
      "$CALL_BUDGET"
  done
  zsh redesign/run_frozen_external_parallel.sh
  for dataset in "${exploratory_datasets[@]}"; do
    "$R_BIN" redesign/run_validation_benchmark.R "$dataset" eval \
      "$CALL_BUDGET"
    "$R_BIN" redesign/run_frozen_external_competitors.R "$dataset" all \
      "$CALL_BUDGET"
    "$R_BIN" redesign/run_frozen_external_adaptive.R "$dataset" all \
      "$CALL_BUDGET"
    "$R_BIN" redesign/run_frozen_external_biology.R "$dataset" all \
      "$CALL_BUDGET"
  done
  "$R_BIN" redesign/summarise_frozen_external.R
  "$R_BIN" redesign/summarise_frozen_external_competitors.R
  "$R_BIN" redesign/summarise_frozen_external_efficiency_biology.R
)

run_review() {
  "$R_BIN" redesign/review_completed_analyses_2026_09_07.R
  python3 redesign/build_full_report_v10.py
  python3 redesign/build_comprehensive_report.py
}

bash analysis/bootstrap_package.sh

case "$scope" in
  tests) run_tests ;;
  primary) run_primary ;;
  ablations) run_ablations ;;
  biology) run_biology ;;
  complementarity) run_complementarity ;;
  exploratory) run_exploratory ;;
  review) run_review ;;
  all)
    run_tests
    run_primary
    run_ablations
    run_biology
    run_exploratory
    run_complementarity
    run_review
    ;;
  *)
    printf 'Unknown scope: %s\n' "$scope" >&2
    printf 'Use tests, primary, ablations, biology, complementarity, exploratory, review or all.\n' >&2
    exit 2
    ;;
esac
