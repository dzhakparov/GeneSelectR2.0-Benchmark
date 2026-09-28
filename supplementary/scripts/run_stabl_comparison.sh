#!/usr/bin/env bash
# Stabl comparison: fixed-size predictive evaluation, stability summary and
# biological assessment. Requires GENESELECTR_STABL_PYTHON to point to a
# Python 3.11 environment with Stabl tag v1.0.1-lw installed (see
# redesign/stabl_comparison/README.md) and the saved outer splits from
# run_primary.sh.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

if [[ -z "${GENESELECTR_STABL_PYTHON:-}" ]]; then
  echo "Set GENESELECTR_STABL_PYTHON to the Python environment containing Stabl." >&2
  exit 1
fi

run_r redesign/stabl_comparison/run_comparison.R
run_r redesign/stabl_comparison/summarise.R
run_r redesign/stabl_comparison/run_biology.R
run_r redesign/stabl_comparison/summarise_biology.R
