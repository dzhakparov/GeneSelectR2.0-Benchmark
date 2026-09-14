#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

# Full paper-scope sequence. The exploratory four-dataset workflow is absent.
bash supplementary/scripts/check_inputs.sh
bash supplementary/scripts/run_tests.sh
bash supplementary/scripts/run_method_development.sh
bash supplementary/scripts/run_component_ablation.sh
bash supplementary/scripts/run_permutation_adjustment.sh
bash supplementary/scripts/run_permutation_diagnostics.sh
bash supplementary/scripts/run_larger_permutation_reference.sh
bash supplementary/scripts/run_random_reference_sensitivity.sh
bash supplementary/scripts/run_gene_set_size_sensitivity.sh
bash supplementary/scripts/run_biology.sh
bash supplementary/scripts/run_dge_complementarity.sh
bash supplementary/scripts/run_figures.sh
