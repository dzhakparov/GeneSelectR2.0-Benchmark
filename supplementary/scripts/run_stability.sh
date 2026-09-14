#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

run_r redesign/complementarity_analysis/01_gene_set_stability.R
