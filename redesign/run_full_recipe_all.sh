#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
exec bash redesign/run_corrected_benchmark.sh cohorts
