#!/usr/bin/env Rscript
# ==============================================================================
#  Horseshoe arm benchmark (added 2026-08-23).
#
#  One arm per split:
#    hs   horseshoe probit fit on the standardised training matrix
#         (n_iter=4000, burn=1000, seed=42, deterministic), ranked by
#         posterior mean |beta| x PIP(|beta| > 0.05)
#
#  Same cached splits, pools, evaluator as every other arm. Output:
#  eval_horseshoe.csv in the dataset's out dir, checkpointed per split.
#
#  Usage: Rscript redesign/run_horseshoe_benchmark.R <dataset> [budget_sec]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
ds_arg <- if (length(args) >= 1) args[1] else stop("dataset required")
budget <- if (length(args) >= 2) as.numeric(args[2]) else 270

validation_ds <- c("GSE65682", "GSE69683", "GSE13355", "GSE107994",
                   "GSE101794")
stopifnot(ds_arg %in% c(validation_ds, "imvigor210", "sosall"))

t_start <- proc.time()[["elapsed"]]
over_budget <- function() (proc.time()[["elapsed"]] - t_start) > budget

suppressPackageStartupMessages(library(glmnet))
for (f in c("bio_prior.R", "evaluator.R", "imvigor210_data.R",
            "horseshoe.R")) {
  source(file.path("redesign", "R", f))
}
source(file.path("redesign", "R", "run_provenance.R"))

is_validation <- ds_arg %in% validation_ds
results_root <- redesign_results_root()
split_dir <- if (is_validation)
  file.path(results_root, "validation_benchmark", ds_arg) else
  file.path(results_root, "grouped_benchmark", ds_arg)
out_dir <- if (is_validation) split_dir else
  file.path(results_root, "full_recipe", ds_arg)
require_redesign_run(split_dir)
if (!identical(out_dir, split_dir)) require_redesign_run(out_dir)
prepare_redesign_extension(
  out_dir, "horseshoe",
  redesign_extension_sources("redesign/run_horseshoe_benchmark.R"),
  config = list(dataset = ds_arg)
)

outcome <- if (is_validation)
  readRDS(file.path(split_dir, "base_data.rds"))$outcome else
  readRDS(file.path(out_dir, "base_outcome.rds"))

split_files <- sort(list.files(split_dir, pattern = "^split_r.*\\.rds$",
                               full.names = TRUE))
stopifnot(length(split_files) == 15)
panel_sizes <- c(10, 20, 50, 100, 200, 500)

eval_path <- file.path(out_dir, "eval_horseshoe.csv")
done <- if (file.exists(eval_path)) read.csv(eval_path) else data.frame()
done_key <- if (nrow(done) > 0)
  paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

for (sf in split_files) {
  m <- regmatches(basename(sf),
                  regexec("split_r(\\d+)_f(\\d+)\\.rds", basename(sf)))[[1]]
  rep_idx <- as.integer(m[2]); fold_idx <- as.integer(m[3])

  need <- any(!(paste(rep_idx, fold_idx, "hs", panel_sizes) %in% done_key))
  if (!need) next
  if (over_budget()) { cat("[budget] stop\n"); quit(save = "no") }

  sp <- readRDS(sf)
  y_train <- outcome[sp$train_idx]
  y_test  <- outcome[sp$test_idx]
  pool <- sp$pools$var2000
  std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                           sp$test_raw[, pool, drop = FALSE])

  fit <- hs_probit(std$train, as.integer(droplevels(y_train)) - 1L,
                   n_iter = 4000, burn = 1000, seed = 42)
  rk <- fit$ranking
  write.csv(
    data.frame(gene = rk),
    file.path(out_dir, sprintf("ranking_r%d_f%d_hs.csv", rep_idx, fold_idx)),
    row.names = FALSE
  )

  rows <- list()
  for (k in panel_sizes) {
    key <- paste(rep_idx, fold_idx, "hs", k)
    if (key %in% done_key) next
    panel <- head(rk, k)
    scores <- predict_with_ensemble(
      std$train[, panel, drop = FALSE], y_train,
      std$test[, panel, drop = FALSE])
    rows[[length(rows) + 1]] <- data.frame(
      repeat_idx = rep_idx, fold_idx = fold_idx, arm = "hs", k = k,
      AUC = bench_auc(y_test, scores), n_panel = length(panel))
    done_key <- c(done_key, key)
  }
  write.table(do.call(rbind, rows), eval_path,
              append = file.exists(eval_path), sep = ",",
              row.names = FALSE, col.names = !file.exists(eval_path))
  cat(sprintf("[%s r%d f%d] hs eval done\n", ds_arg, rep_idx, fold_idx))
}

cat(sprintf("\nDone dataset=%s (%.0f s)\n", ds_arg,
            proc.time()[["elapsed"]] - t_start))
