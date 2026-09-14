#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

run_r redesign/complementarity_analysis/01_gene_set_stability.R
run_r redesign/complementarity_analysis/02_dge_geneselectr_concordance.R
run_r redesign/complementarity_analysis/03_gene_group_characterization.R
run_r redesign/complementarity_analysis/04_gene_group_biology.R
