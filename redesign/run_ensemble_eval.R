#!/usr/bin/env Rscript
# ==============================================================================
#  Ensemble evaluation: combines EXISTING rankings, no new model fits.
#
#  Arms per split:
#    ens_rank   rank-average of GS_full_grouped + predfirst_raw
#    ens_rank3  rank-average of those + RF_importance
#    kswitch    predfirst_raw for k <= 20, GS_full_grouped for k > 20
#
#  Rankings come from the same directories the earlier runs wrote, so all
#  splits/pools/evaluator stay identical. Output: eval_ensemble.csv next to
#  the existing eval files.
#
#  Usage: Rscript redesign/run_ensemble_eval.R <dataset> [budget_sec]
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
for (f in list.files("package/GeneSelectR/R", full.names = TRUE)) source(f)
for (f in c("bio_prior.R", "evaluator.R", "imvigor210_data.R",
            "ensemble.R")) {
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
  out_dir, "predfirst",
  redesign_extension_sources("redesign/run_predfirst_benchmark.R"),
  config = list(dataset = ds_arg)
)
prepare_redesign_extension(
  out_dir, "ensemble",
  redesign_extension_sources("redesign/run_ensemble_eval.R"),
  config = list(dataset = ds_arg)
)

#  Arm names differ between the two result trees.
nm_grouped   <- if (is_validation) "GS_full_grouped" else "full_grouped"
nm_ungrouped <- if (is_validation) "GS_full_ungrouped" else "full_ungrouped"

outcome <- if (is_validation)
  readRDS(file.path(split_dir, "base_data.rds"))$outcome else
  readRDS(file.path(out_dir, "base_outcome.rds"))

split_files <- sort(list.files(split_dir, pattern = "^split_r.*\\.rds$",
                               full.names = TRUE))
stopifnot(length(split_files) == 15)
panel_sizes <- c(10, 20, 50, 100, 200, 500)
k_cut <- 20

eval_path <- file.path(out_dir, "eval_ensemble.csv")
done <- if (file.exists(eval_path)) read.csv(eval_path) else data.frame()
done_key <- if (nrow(done) > 0)
  paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

get_ranking <- function(rep_idx, fold_idx, arm) {
  p <- file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv", rep_idx,
                                  fold_idx, arm))
  rk <- read.csv(p)$gene
  if (length(rk) == 0 && arm == nm_grouped) {
    #  Empty gate fallback: the ungrouped ranking, as in the main runs.
    rk <- read.csv(file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv",
                     rep_idx, fold_idx, nm_ungrouped)))$gene
  }
  rk
}

for (sf in split_files) {
  m <- regmatches(basename(sf),
                  regexec("split_r(\\d+)_f(\\d+)\\.rds", basename(sf)))[[1]]
  rep_idx <- as.integer(m[2]); fold_idx <- as.integer(m[3])

  need <- any(!(paste(rep_idx, fold_idx,
                      rep(c("ens_rank", "ens_rank3", "kswitch"), each = 6),
                      panel_sizes) %in% done_key))
  if (!need) next
  if (over_budget()) { cat("[budget] stop\n"); quit(save = "no") }

  sp <- readRDS(sf)
  y_train <- outcome[sp$train_idx]
  y_test  <- outcome[sp$test_idx]
  pool <- sp$pools$var2000
  std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                           sp$test_raw[, pool, drop = FALSE])

  rk_g <- get_ranking(rep_idx, fold_idx, nm_grouped)
  rk_p <- read.csv(file.path(out_dir, sprintf(
    "ranking_r%d_f%d_predfirst_raw.csv", rep_idx, fold_idx)))$gene
  rk_rf <- get_ranking(rep_idx, fold_idx, "RF_importance")

  ens2 <- rank_average(list(rk_g, rk_p), pool)
  ens3 <- rank_average(list(rk_g, rk_p, rk_rf), pool)

  # Persist the exact rankings used for evaluation. The biology assessment
  # reads these files so prediction and STRING connectivity use the same
  # panels. kswitch uses the grouped ranking at the biology panel size k=50.
  biology_rankings <- list(
    ens_rank = ens2,
    ens_rank3 = ens3,
    kswitch = kswitch_ranking(50, k_cut, rk_p, rk_g)
  )
  for (arm in names(biology_rankings)) {
    write.csv(
      data.frame(gene = biology_rankings[[arm]]),
      file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv", rep_idx,
                                 fold_idx, arm)),
      row.names = FALSE
    )
  }

  rows <- list()
  for (k in panel_sizes) {
    panels <- list(
      ens_rank  = head(ens2, k),
      ens_rank3 = head(ens3, k),
      kswitch   = head(kswitch_ranking(k, k_cut, rk_p, rk_g), k)
    )
    for (arm in names(panels)) {
      key <- paste(rep_idx, fold_idx, arm, k)
      if (key %in% done_key) next
      panel <- panels[[arm]]
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
  cat(sprintf("[%s r%d f%d] ensemble eval done\n", ds_arg, rep_idx,
              fold_idx))
}

cat(sprintf("\nDone dataset=%s (%.0f s)\n", ds_arg,
            proc.time()[["elapsed"]] - t_start))
