#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

bash analysis/bootstrap_package.sh
bash analysis/check_repository_contents.sh
run_r redesign/run_tests.R
