#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
stage <- if (length(args) >= 2L) args[[2L]] else "all"
budget <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 604800
datasets <- c("GSE16879", "GSE91061", "GSE92415", "GSE206285")
stopifnot(dataset %in% datasets, stage %in% c("init", "fit", "eval", "all"))

started_at <- proc.time()[["elapsed"]]
over_budget <- function() proc.time()[["elapsed"]] - started_at > budget

suppressPackageStartupMessages({
  library(Boruta)
  library(glmnet)
  library(mRMRe)
  library(ranger)
  library(withr)
  library(xgboost)
})
options(warn = 1)

source(file.path("redesign", "R", "bio_prior.R"))
source(file.path("redesign", "R", "competitors.R"))
source(file.path("redesign", "R", "evaluator.R"))
source(file.path("redesign", "R", "run_provenance.R"))

result_root <- file.path(
  "redesign", "results_frozen_external_exact_2026-08-31"
)
dataset_dir <- file.path(result_root, "validation_benchmark", dataset)
methods <- c("DGE", "LASSO", "ElasticNet", "mRMR", "Boruta",
             "RF_importance")
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
extension <- "classical_comparators_v1"
evaluation_path <- file.path(dataset_dir, "competitor_eval_results.csv")

source_files <- c(
  "redesign/run_frozen_external_competitors.R",
  file.path("redesign", "R",
            c("bio_prior.R", "competitors.R", "evaluator.R",
              "run_provenance.R"))
)
prepare_redesign_extension(
  dataset_dir,
  extension = extension,
  source_files = source_files,
  config = list(
    dataset = dataset,
    methods = methods,
    panel_sizes = panel_sizes,
    ranking_seed_base = 520000L,
    evaluation_seed_base = 420000L,
    evaluator = "glmnet+xgboost+ranger"
  ),
  output_files = evaluation_path
)

if (stage == "init") {
  cat(sprintf("%s comparator extension initialized\n", dataset))
  quit(save = "no")
}

base <- readRDS(file.path(dataset_dir, "base_data.rds"))
split_grid <- expand.grid(repeat_idx = 1:3, fold_idx = 1:5)
split_override <- trimws(Sys.getenv("GENESELECTR_VALIDATION_SPLIT", ""))
if (nzchar(split_override)) {
  split_keys <- sprintf("r%df%d", split_grid$repeat_idx, split_grid$fold_idx)
  if (!(split_override %in% split_keys)) {
    stop("Unknown validation split: ", split_override, call. = FALSE)
  }
  split_grid <- split_grid[split_keys == split_override, , drop = FALSE]
}

ranking_path <- function(repeat_idx, fold_idx, method) {
  file.path(dataset_dir, sprintf(
    "competitor_r%d_f%d_%s.csv", repeat_idx, fold_idx, method
  ))
}

ranking_seed <- function(repeat_idx, fold_idx, method) {
  520000L + repeat_idx * 1000L + fold_idx * 100L + match(method, methods)
}

fit_competitor <- function(method, train_features, train_labels,
                           repeat_idx, fold_idx) {
  seed <- ranking_seed(repeat_idx, fold_idx, method)
  withr::with_seed(seed, switch(
    method,
    DGE = rank_by_differential_expression(train_features, train_labels),
    LASSO = rank_by_lasso(train_features, train_labels),
    ElasticNet = rank_by_elastic_net(train_features, train_labels),
    mRMR = rank_by_mrmr(train_features, train_labels),
    Boruta = rank_by_boruta(train_features, train_labels),
    RF_importance = rank_by_random_forest(
      train_features, train_labels, random_seed = seed
    ),
    stop("Unknown comparator: ", method, call. = FALSE)
  ))
}

if (stage %in% c("fit", "all")) {
  for (split_index in seq_len(nrow(split_grid))) {
    repeat_idx <- split_grid$repeat_idx[split_index]
    fold_idx <- split_grid$fold_idx[split_index]
    split <- readRDS(file.path(
      dataset_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
    ))
    pool <- split$pools$var2000
    standardised <- standardise_split(
      split$train_raw[, pool, drop = FALSE],
      split$test_raw[, pool, drop = FALSE]
    )
    train_labels <- base$outcome[split$train_idx]

    for (method in methods) {
      output_path <- ranking_path(repeat_idx, fold_idx, method)
      if (file.exists(output_path)) next
      if (over_budget()) {
        cat("[budget] stop in comparator fit\n")
        quit(save = "no")
      }

      method_seed <- ranking_seed(repeat_idx, fold_idx, method)
      method_started <- proc.time()[["elapsed"]]
      ranking <- fit_competitor(
        method, standardised$train, train_labels, repeat_idx, fold_idx
      )
      elapsed <- proc.time()[["elapsed"]] - method_started
      if (length(ranking$ranked) != length(pool) ||
          anyDuplicated(ranking$ranked) ||
          !setequal(ranking$ranked, pool)) {
        stop(sprintf("%s produced an invalid ranking in r%d f%d",
                     method, repeat_idx, fold_idx), call. = FALSE)
      }
      write.csv(data.frame(gene = ranking$ranked), output_path,
                row.names = FALSE)
      write.csv(data.frame(
        method = method,
        ranking_seed = method_seed,
        pool_size = length(pool),
        selected_size = length(ranking$selected),
        fit_seconds = elapsed
      ), sub("[.]csv$", "_meta.csv", output_path), row.names = FALSE)
      cat(sprintf("[%s r%d f%d] %s %.1f s\n", dataset, repeat_idx,
                  fold_idx, method, elapsed))
    }
  }
}

if (stage %in% c("eval", "all")) {
  completed <- if (file.exists(evaluation_path)) {
    read.csv(evaluation_path, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  completed_keys <- if (nrow(completed)) {
    paste(completed$repeat_idx, completed$fold_idx,
          completed$method, completed$k)
  } else {
    character(0)
  }

  for (split_index in seq_len(nrow(split_grid))) {
    repeat_idx <- split_grid$repeat_idx[split_index]
    fold_idx <- split_grid$fold_idx[split_index]
    split <- readRDS(file.path(
      dataset_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
    ))
    pool <- split$pools$var2000
    standardised <- standardise_split(
      split$train_raw[, pool, drop = FALSE],
      split$test_raw[, pool, drop = FALSE]
    )
    train_labels <- base$outcome[split$train_idx]
    test_labels <- base$outcome[split$test_idx]

    for (method in methods) {
      path <- ranking_path(repeat_idx, fold_idx, method)
      if (!file.exists(path)) {
        stop("Missing comparator ranking: ", path, call. = FALSE)
      }
      ranking <- read.csv(path, stringsAsFactors = FALSE)$gene
      for (panel_size in panel_sizes) {
        key <- paste(repeat_idx, fold_idx, method, panel_size)
        if (key %in% completed_keys) next
        if (over_budget()) {
          cat("[budget] stop in comparator eval\n")
          quit(save = "no")
        }
        panel <- head(ranking, panel_size)
        evaluation_seed <- 420000L + repeat_idx * 1000L + fold_idx * 100L +
          match(panel_size, panel_sizes)
        scores <- predict_with_ensemble(
          standardised$train[, panel, drop = FALSE], train_labels,
          standardised$test[, panel, drop = FALSE],
          random_seed = evaluation_seed
        )
        row <- data.frame(
          repeat_idx = repeat_idx,
          fold_idx = fold_idx,
          method = method,
          k = panel_size,
          AUC = bench_auc(test_labels, scores),
          n_panel = length(panel),
          ranking_seed = ranking_seed(repeat_idx, fold_idx, method),
          evaluation_seed = evaluation_seed,
          evaluator_components = paste(attr(scores, "components_used"),
                                       collapse = "+")
        )
        write.table(
          row, evaluation_path, append = file.exists(evaluation_path),
          sep = ",", row.names = FALSE,
          col.names = !file.exists(evaluation_path)
        )
        completed_keys <- c(completed_keys, key)
      }
    }
  }
}

cat(sprintf("%s comparator stage=%s complete\n", dataset, stage))
