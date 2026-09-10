#!/usr/bin/env Rscript
# ==============================================================================
#  Prediction-first variant benchmark (changes 1+2+3), all 7 datasets.
#
#  Reuses the CACHED splits from the earlier runs, so every split is
#  identical to what GS_full_* and the competitors already saw:
#    validation datasets -> <corrected root>/validation_benchmark/<ds>/
#    imvigor210, sosall  -> <corrected root>/grouped_benchmark/<cohort>/
#                           (outcome from full_recipe/<cohort>/base_outcome.rds)
#
#  For imvigor210 and sosall this also produces FRESH competitor rankings
#  under the current protocol (the saved anchors were p=5000/n=192 on
#  IMvigor210 -- not like-for-like).
#
#  Arms written per split:
#    predfirst       score ranking + redundancy filter (the candidate)
#    predfirst_raw   score ranking only (ablation for change 2)
#    (+ DGE, LASSO, ElasticNet, mRMR, Boruta, RF_importance on imvigor/sosall)
#
#  Usage: Rscript redesign/run_predfirst_benchmark.R <dataset> <stage> [budget]
#    stage: fit | eval | all
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
ds_arg <- if (length(args) >= 1) args[1] else stop("dataset required")
stage  <- if (length(args) >= 2) args[2] else "all"
budget <- if (length(args) >= 3) as.numeric(args[3]) else 270

validation_ds <- c("GSE65682", "GSE69683", "GSE13355", "GSE107994",
                   "GSE101794")
cohort_ds     <- c("imvigor210", "sosall")
stopifnot(ds_arg %in% c(validation_ds, cohort_ds),
          stage %in% c("fit", "eval", "all"))

t_start <- proc.time()[["elapsed"]]
over_budget <- function() (proc.time()[["elapsed"]] - t_start) > budget

suppressPackageStartupMessages(library(glmnet))
for (f in list.files("package/GeneSelectR/R", full.names = TRUE)) source(f)
for (f in c("bio_prior.R", "evaluator.R", "imvigor210_data.R",
            "competitors.R", "predfirst.R")) {
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
  out_dir, "predfirst",
  redesign_extension_sources("redesign/run_predfirst_benchmark.R"),
  config = list(dataset = ds_arg)
)

#  outcome
if (is_validation) {
  outcome <- readRDS(file.path(split_dir, "base_data.rds"))$outcome
} else {
  outcome <- readRDS(file.path(out_dir, "base_outcome.rds"))
}

split_files <- sort(list.files(split_dir, pattern = "^split_r.*\\.rds$",
                               full.names = TRUE))
stopifnot(length(split_files) == 15)

random_seed <- 42
panel_sizes <- c(10, 20, 50, 100, 200, 500)
comp_arms <- c("DGE", "LASSO", "ElasticNet", "mRMR", "Boruta",
               "RF_importance")

split_ids <- function(path) {
  m <- regmatches(basename(path),
                  regexec("split_r(\\d+)_f(\\d+)\\.rds", basename(path)))[[1]]
  c(rep = as.integer(m[2]), fold = as.integer(m[3]))
}

# ---- fit --------------------------------------------------------------------
if (stage %in% c("fit", "all")) {
  for (sf in split_files) {
    ids <- split_ids(sf)
    rank_path <- file.path(out_dir, sprintf("ranking_r%d_f%d_predfirst.csv",
                                            ids["rep"], ids["fold"]))
    need_pf <- !file.exists(rank_path)
    need_comp <- !is_validation && !all(file.exists(file.path(
      out_dir, sprintf("ranking_r%d_f%d_%s.csv", ids["rep"], ids["fold"],
                       comp_arms))))
    if (!need_pf && !need_comp) next
    if (over_budget()) { cat("[budget] stop in fit\n"); quit(save = "no") }

    sp <- readRDS(sf)
    y_train <- outcome[sp$train_idx]
    std <- standardise_split(sp$train_raw[, sp$pools$var2000, drop = FALSE],
                             sp$test_raw[, sp$pools$var2000, drop = FALSE])

    if (need_pf) {
      t0 <- proc.time()[["elapsed"]]
      pf <- fit_predfirst(std$train, y_train, B = 50, k_folds = 5,
                          alpha_grid = c(0.5, 1.0), gamma = 0.25,
                          random_seed = random_seed, n_cores = n_workers)
      el <- proc.time()[["elapsed"]] - t0
      filt <- redundancy_filter(pf$ranking, std$train, tau = 0.7)
      write.csv(data.frame(gene = filt), rank_path, row.names = FALSE)
      write.csv(data.frame(gene = pf$ranking),
                sub("predfirst\\.csv", "predfirst_raw.csv", rank_path),
                row.names = FALSE)
      write.csv(data.frame(
        key = c("alpha", "cvm_a0.5", "cvm_a1.0", "fit_seconds"),
        value = c(pf$alpha, round(pf$cvm, 4), round(el, 1))),
        sub("\\.csv$", "_meta.csv", rank_path), row.names = FALSE)
      cat(sprintf("[%s r%d f%d] predfirst: %.0f s, alpha=%.1f, top: %s\n",
                  ds_arg, ids["rep"], ids["fold"], el, pf$alpha,
                  paste(head(pf$ranking, 3), collapse = ", ")))
    }

    if (need_comp) {
      for (arm in comp_arms) {
        cp <- file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv",
                                         ids["rep"], ids["fold"], arm))
        if (file.exists(cp)) next
        if (over_budget()) { cat("[budget] stop in comp\n")
          quit(save = "no") }
        rk <- switch(arm,
          DGE           = rank_by_differential_expression(std$train, y_train),
          LASSO         = rank_by_lasso(std$train, y_train),
          ElasticNet    = rank_by_elastic_net(std$train, y_train),
          mRMR          = rank_by_mrmr(std$train, y_train),
          Boruta        = rank_by_boruta(std$train, y_train),
          RF_importance = rank_by_random_forest(std$train, y_train,
                                                random_seed = random_seed))
        stopifnot(length(rk$ranked) == ncol(std$train))
        write.csv(data.frame(gene = rk$ranked), cp, row.names = FALSE)
        cat(sprintf("[%s r%d f%d] %s done\n", ds_arg, ids["rep"],
                    ids["fold"], arm))
      }
    }
  }
}

# ---- eval -------------------------------------------------------------------
if (stage %in% c("eval", "all")) {
  eval_path <- file.path(out_dir, "eval_predfirst.csv")
  done <- if (file.exists(eval_path)) read.csv(eval_path) else data.frame()
  done_key <- if (nrow(done) > 0)
    paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

  arms <- if (is_validation) c("predfirst", "predfirst_raw") else
    c("predfirst", "predfirst_raw", comp_arms)

  for (sf in split_files) {
    ids <- split_ids(sf)
    rank_files <- setNames(file.path(
      out_dir, sprintf("ranking_r%d_f%d_%s.csv", ids["rep"], ids["fold"],
                       arms)), arms)
    if (!all(file.exists(rank_files))) next

    sp <- readRDS(sf)
    y_train <- outcome[sp$train_idx]
    y_test  <- outcome[sp$test_idx]
    pool <- sp$pools$var2000
    std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                             sp$test_raw[, pool, drop = FALSE])

    rows <- list()
    flush_rows <- function() {
      if (length(rows) == 0) return(invisible(NULL))
      write.table(do.call(rbind, rows), eval_path,
                  append = file.exists(eval_path), sep = ",",
                  row.names = FALSE, col.names = !file.exists(eval_path))
      rows <<- list()
    }

    for (arm in arms) {
      ranking <- read.csv(rank_files[[arm]])$gene
      for (k in panel_sizes) {
        key <- paste(ids["rep"], ids["fold"], arm, k)
        if (key %in% done_key) next
        if (over_budget()) { flush_rows()
          cat("[budget] stop in eval\n"); quit(save = "no") }
        panel <- head(ranking, min(k, length(ranking)))
        scores <- predict_with_ensemble(
          std$train[, panel, drop = FALSE], y_train,
          std$test[, panel, drop = FALSE])
        rows[[length(rows) + 1]] <- data.frame(
          repeat_idx = ids["rep"], fold_idx = ids["fold"], arm = arm, k = k,
          AUC = bench_auc(y_test, scores), n_panel = length(panel))
        done_key <- c(done_key, key)
        flush_rows()
      }
    }

    #  Random baseline for the cohort datasets (validation ones already have
    #  it in eval_results.csv; the report merges).
    if (!is_validation) {
      for (k in panel_sizes) {
        key <- paste(ids["rep"], ids["fold"], "Random", k)
        if (key %in% done_key) next
        aucs <- random_panel_aucs(std$train, std$test, y_train, y_test,
                                  pool, k, n_draws = 3,
                                  seed = 99 + 1000 * ids["rep"] + ids["fold"])
        rows[[length(rows) + 1]] <- data.frame(
          repeat_idx = ids["rep"], fold_idx = ids["fold"], arm = "Random",
          k = k, AUC = mean(aucs), n_panel = min(k, length(pool)))
        done_key <- c(done_key, key)
        flush_rows()
      }
    }
  }
}

cat(sprintf("\nDone stage=%s dataset=%s (%.0f s)\n", stage, ds_arg,
            proc.time()[["elapsed"]] - t_start))
