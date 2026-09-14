#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

run_r redesign/run_older7_calibration_onoff.R all "$GENESELECTR_BUDGET"
