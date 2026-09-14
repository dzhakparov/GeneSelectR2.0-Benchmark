#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

# Runs the maintained seven-dataset extension sequence. This includes the
# recorded prediction-first, prior, ensemble, weight, horseshoe, pruning,
# adaptive, STRING, and disease-assessment extensions.
bash redesign/run_corrected_benchmark.sh all
