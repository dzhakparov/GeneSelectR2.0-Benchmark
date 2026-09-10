#!/usr/bin/env Rscript
# ==============================================================================
#  Post-result boost-weight extension (added 2026-08-23). This evaluates
#  cb_gs at w = 2.0 using cached GS_full_ungrouped scores and the
#  zscores/ddzscores caches. No model or null fit is repeated. Output:
#  eval_wsweep.csv in the dataset's output directory.
#
#  Usage: Rscript redesign/run_wsweep.R <dataset> [budget_sec]
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
for (f in c("bio_prior.R", "evaluator.R", "imvigor210_data.R")) {
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
require_redesign_extension(
  out_dir, "softprior",
  redesign_extension_sources("redesign/run_softprior_benchmark.R"),
  config = list(dataset = ds_arg)
)
require_redesign_extension(
  out_dir, "datadriven",
  redesign_extension_sources("redesign/run_datadriven_benchmark.R"),
  config = list(dataset = ds_arg)
)
prepare_redesign_extension(
  out_dir, "wsweep",
  redesign_extension_sources("redesign/run_wsweep.R"),
  config = list(dataset = ds_arg)
)

nm_ungrouped <- if (is_validation) "GS_full_ungrouped" else "full_ungrouped"

outcome <- if (is_validation)
  readRDS(file.path(split_dir, "base_data.rds"))$outcome else
  readRDS(file.path(out_dir, "base_outcome.rds"))

split_files <- sort(list.files(split_dir, pattern = "^split_r.*\\.rds$",
                               full.names = TRUE))
stopifnot(length(split_files) == 15)
panel_sizes <- c(10, 20, 50, 100, 200, 500)
ws <- c(2.0)
arms <- "cb_gs_w2"

eval_path <- file.path(out_dir, "eval_wsweep.csv")
done <- if (file.exists(eval_path)) read.csv(eval_path) else data.frame()
done_key <- if (nrow(done) > 0)
  paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

for (sf in split_files) {
  m <- regmatches(basename(sf),
                  regexec("split_r(\\d+)_f(\\d+)\\.rds", basename(sf)))[[1]]
  rep_idx <- as.integer(m[2]); fold_idx <- as.integer(m[3])

  need <- any(!(paste(rep_idx, fold_idx, rep(arms, each = 6),
                      panel_sizes) %in% done_key))
  if (!need) next
  if (over_budget()) { cat("[budget] stop\n"); quit(save = "no") }

  sp <- readRDS(sf)
  y_train <- outcome[sp$train_idx]
  y_test  <- outcome[sp$test_idx]
  pool <- sp$pools$var2000
  std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                           sp$test_raw[, pool, drop = FALSE])

  zs <- readRDS(file.path(out_dir, sprintf("zscores_r%d_f%d.rds", rep_idx,
                                           fold_idx)))
  dd <- readRDS(file.path(out_dir, sprintf("ddzscores_r%d_f%d.rds", rep_idx,
                                           fold_idx)))
  sets_cb <- c(zs$sets, dd$sets)
  z_cb <- c(zs$z, dd$z)

  gs_tab <- read.csv(file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv",
                               rep_idx, fold_idx, nm_ungrouped)))
  gs_scores <- setNames(gs_tab$final_score, gs_tab$gene)
  gs_scores <- gs_scores[intersect(names(gs_scores), pool)]

  rankings <- lapply(ws, function(w)
    soft_prior_rescore(gs_scores, sets_cb, z_cb, w = w))
  names(rankings) <- arms
  for (arm in names(rankings)) {
    write.csv(
      data.frame(gene = rankings[[arm]]),
      file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv", rep_idx,
                                 fold_idx, arm)),
      row.names = FALSE
    )
  }

  rows <- list()
  for (k in panel_sizes) {
    for (arm in arms) {
      key <- paste(rep_idx, fold_idx, arm, k)
      if (key %in% done_key) next
      panel <- head(rankings[[arm]], k)
      scores <- predict_with_ensemble(
        std$train[, panel, drop = FALSE], y_train,
        std$test[, panel, drop = FALSE])
      rows[[length(rows) + 1]] <- data.frame(
        repeat_idx = rep_idx, fold_idx = fold_idx, arm = arm, k = k,
        AUC = bench_auc(y_test, scores), n_panel = length(panel))
      done_key <- c(done_key, key)
    }
  }
  if (length(rows) > 0) {
    write.table(do.call(rbind, rows), eval_path,
                append = file.exists(eval_path), sep = ",",
                row.names = FALSE, col.names = !file.exists(eval_path))
  }
  cat(sprintf("[%s r%d f%d] wsweep done\n", ds_arg, rep_idx, fold_idx))
}

cat(sprintf("\nDone dataset=%s (%.0f s)\n", ds_arg,
            proc.time()[["elapsed"]] - t_start))
