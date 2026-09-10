#!/bin/zsh

# Run one outer split per process. GeneSelectR uses one core within each
# process, so the worker count is also the maximum parallel model count.
set -eu

workers="${1:-${GENESELECTR_OUTER_WORKERS:-2}}"
project_dir="${0:A:h:h}"
result_root="$project_dir/redesign/results_frozen_external_exact_2026-08-31"
log_dir="$result_root/logs/adaptive_biology"
job_file="$result_root/adaptive_jobs_v3.tsv"

mkdir -p "$log_dir"
: > "$job_file"
for dataset in GSE16879 GSE91061 GSE92415 GSE206285; do
  for repeat_idx in 1 2 3; do
    for fold_idx in 1 2 3 4 5; do
      print -r -- "${dataset}"$'\t'"r${repeat_idx}f${fold_idx}" >> "$job_file"
    done
  done
done

cd "$project_dir"
cat "$job_file" | xargs -P "$workers" -n 2 zsh -c '
  dataset="$1"
  split_name="$2"
  log_path="redesign/results_frozen_external_exact_2026-08-31/logs/adaptive_biology/${dataset}_${split_name}_adaptive_v3.log"
  R_USER_CACHE_DIR="$PWD/data/r_user_cache" \
    GENESELECTR_VALIDATION_SPLIT="$split_name" \
    Rscript redesign/run_frozen_external_adaptive.R \
      "$dataset" inner 604800 > "$log_path" 2>&1
' adaptive-worker

print -r -- "Adaptive inner-fold jobs complete."
