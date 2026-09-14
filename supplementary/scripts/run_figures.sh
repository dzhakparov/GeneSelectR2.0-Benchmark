#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/_common.sh"

run_r supplementary/scripts/generate_paper_figures.R
