#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

bash analysis/bootstrap_package.sh
run_r analysis/check_inputs.R
