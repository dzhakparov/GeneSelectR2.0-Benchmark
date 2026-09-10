#!/usr/bin/env Rscript
# ==============================================================================
#  Soft module prior benchmark. Replaces the all-or-nothing Hallmark gate
#  with a continuous multiplicative boost:
#
#    module z     = (mean|t| observed - null mean) / null sd, B=1000 label
#                   permutations, same statistic/null as module_gate
#    multiplier   = 1 + 0.5 * max(0, best z among the gene's sets)
#    unannotated  = 1 (no set membership, no penalty, no veto)
#
#  Arms per split:
#    soft_gs  base = GS_full_ungrouped gene_scores (loaded from disk,
#             re-scored, no refit)
#    soft_pf  base = predfirst_raw (refit per split, deterministic via
#             fixed seed, same call as the original runs)
#
#  Module z is computed on the raw training matrix restricted to the split's
#  var2000 pool and cached as zscores_r*_f*.rds in the dataset's out dir.
#
#  Usage: Rscript redesign/run_softprior_benchmark.R <dataset> [budget_sec]
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
            "predfirst.R")) {
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
n_workers <- redesign_worker_count()
prepare_redesign_extension(
  out_dir, "softprior",
  redesign_extension_sources("redesign/run_softprior_benchmark.R"),
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

eval_path <- file.path(out_dir, "eval_softprior.csv")
done <- if (file.exists(eval_path)) read.csv(eval_path) else data.frame()
done_key <- if (nrow(done) > 0)
  paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

get_zscores <- function(rep_idx, fold_idx, X_train_pool, y_train) {
  zp <- file.path(out_dir, sprintf("zscores_r%d_f%d.rds", rep_idx, fold_idx))
  if (file.exists(zp)) return(readRDS(zp))
  ms <- module_scores(X_train_pool, as.integer(droplevels(y_train)) - 1L,
                      B = 1000, seed = 7)
  saveRDS(list(z = ms$z, sets = ms$sets, table = ms$table), zp)
  list(z = ms$z, sets = ms$sets, table = ms$table)
}

for (sf in split_files) {
  m <- regmatches(basename(sf),
                  regexec("split_r(\\d+)_f(\\d+)\\.rds", basename(sf)))[[1]]
  rep_idx <- as.integer(m[2]); fold_idx <- as.integer(m[3])

  need <- any(!(paste(rep_idx, fold_idx,
                      rep(c("soft_gs", "soft_pf"), each = 6),
                      panel_sizes) %in% done_key))
  if (!need) next
  if (over_budget()) { cat("[budget] stop\n"); quit(save = "no") }

  sp <- readRDS(sf)
  y_train <- outcome[sp$train_idx]
  y_test  <- outcome[sp$test_idx]
  pool <- sp$pools$var2000
  std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                           sp$test_raw[, pool, drop = FALSE])

  zs <- get_zscores(rep_idx, fold_idx, sp$train_raw[, pool, drop = FALSE],
                    y_train)

  #  soft_gs: re-score the saved GS_full_ungrouped gene_scores table.
  gs_file <- file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv", rep_idx,
                                        fold_idx, nm_ungrouped))
  gs_tab <- read.csv(gs_file)
  gs_scores <- setNames(gs_tab$final_score, gs_tab$gene)
  gs_scores <- gs_scores[intersect(names(gs_scores), pool)]
  rk_soft_gs <- soft_prior_rescore(gs_scores, zs$sets, zs$z, w = 0.5)

  #  soft_pf: refit predfirst_raw (deterministic), then re-score.
  pf <- fit_predfirst(std$train, y_train, B = 50, k_folds = 5,
                      alpha_grid = c(0.5, 1.0), gamma = 0.25,
                      random_seed = 42, n_cores = n_workers)
  rk_soft_pf <- soft_prior_rescore(pf$score, zs$sets, zs$z, w = 0.5)

  rankings <- list(soft_gs = rk_soft_gs, soft_pf = rk_soft_pf)
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
    panels <- lapply(rankings, head, n = k)
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
  cat(sprintf("[%s r%d f%d] soft-prior eval done\n", ds_arg, rep_idx,
              fold_idx))
}

cat(sprintf("\nDone dataset=%s (%.0f s)\n", ds_arg,
            proc.time()[["elapsed"]] - t_start))
