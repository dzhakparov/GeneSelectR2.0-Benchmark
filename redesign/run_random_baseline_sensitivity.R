#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
n_draws <- if (length(args) >= 2L) as.integer(args[[2L]]) else 30L
datasets <- c("GSE16879", "GSE91061", "GSE92415", "GSE206285")
stopifnot(dataset %in% datasets, n_draws >= 3L)

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(withr)
  library(xgboost)
})
source(file.path("redesign", "R", "bio_prior.R"))
source(file.path("redesign", "R", "evaluator.R"))

result_root <- file.path(
  "redesign", "results_frozen_external_exact_2026-08-31"
)
dataset_dir <- file.path(result_root, "validation_benchmark", dataset)
output_path <- file.path(dataset_dir, sprintf(
  "random_baseline_sensitivity_%d_draws.csv", n_draws
))
base <- readRDS(file.path(dataset_dir, "base_data.rds"))
primary_panel_sizes <- c(10L, 20L, 50L)
all_panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)

rows <- list()
for (repeat_idx in 1:3) {
  for (fold_idx in 1:5) {
    split <- readRDS(file.path(
      dataset_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
    ))
    pool <- split$pools$var2000
    standardised <- standardise_split(
      split$train_raw[, pool, drop = FALSE],
      split$test_raw[, pool, drop = FALSE]
    )
    for (panel_size in primary_panel_sizes) {
      evaluation_seed <- 420000L + repeat_idx * 1000L + fold_idx * 100L +
        match(panel_size, all_panel_sizes)
      aucs <- random_panel_aucs(
        standardised$train, standardised$test,
        base$outcome[split$train_idx], base$outcome[split$test_idx],
        pool, panel_size, n_draws = n_draws,
        seed = 99L + 1000L * repeat_idx + fold_idx,
        model_seed = evaluation_seed
      )
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = dataset,
        repeat_idx = repeat_idx,
        fold_idx = fold_idx,
        k = panel_size,
        draw = seq_len(n_draws),
        AUC = aucs,
        evaluation_seed = evaluation_seed
      )
    }
  }
}

result <- do.call(rbind, rows)
keys <- paste(result$repeat_idx, result$fold_idx, result$k, result$draw)
stopifnot(
  nrow(result) == 15L * 3L * n_draws,
  !anyDuplicated(keys),
  all(is.finite(result$AUC)),
  all(result$AUC >= 0 & result$AUC <= 1)
)
write.csv(result, output_path, row.names = FALSE)
cat(sprintf("%s: wrote %d Random-panel evaluations to %s\n",
            dataset, nrow(result), output_path))
