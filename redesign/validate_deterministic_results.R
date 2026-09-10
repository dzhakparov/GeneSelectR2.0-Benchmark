#!/usr/bin/env Rscript
# Validate the deterministic benchmark tables and independently recompute one
# saved-ranking evaluation per dataset.

args <- commandArgs(trailingOnly = TRUE)
results_root <- if (length(args) >= 1L) args[[1L]] else
  file.path("redesign", "results_corrected")
output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(results_root, "full_benchmark_deterministic")

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(withr)
  library(xgboost)
})
source(file.path("redesign", "R", "bio_prior.R"))
source(file.path("redesign", "R", "evaluator.R"))

validation_datasets <- c("GSE65682", "GSE69683", "GSE13355", "GSE107994",
                         "GSE101794")
datasets <- c(validation_datasets, "imvigor210", "sosall")
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
expected_methods <- 27L
expected_rows <- (expected_methods + 1L) * 15L * length(panel_sizes)
qa_rows <- list()

for (dataset in datasets) {
  is_validation <- dataset %in% validation_datasets
  split_dir <- if (is_validation) {
    file.path(results_root, "validation_benchmark", dataset)
  } else {
    file.path(results_root, "grouped_benchmark", dataset)
  }
  ranking_dir <- if (is_validation) split_dir else
    file.path(results_root, "full_recipe", dataset)
  result <- read.csv(file.path(ranking_dir, "eval_deterministic.csv"),
                     stringsAsFactors = FALSE)
  keys <- paste(result$repeat_idx, result$fold_idx, result$arm, result$k)
  methods <- setdiff(unique(result$arm), "Random")
  control <- read.csv(file.path(ranking_dir, "reproducibility_control.csv"),
                      stringsAsFactors = FALSE)

  stopifnot(
    nrow(result) == expected_rows,
    !anyDuplicated(keys),
    length(methods) == expected_methods,
    !("dc_pf_ctrl" %in% methods),
    all(is.finite(result$AUC)),
    all(result$AUC >= 0 & result$AUC <= 1),
    nrow(control) == 15L,
    all(control$identical)
  )

  # Independent spot check: reconstruct dc_pf, k=50, repeat 1, fold 1 from
  # the stored outer split and compare its AUC with the written table.
  split <- readRDS(file.path(split_dir, "split_r1_f1.rds"))
  outcome <- if (is_validation) {
    readRDS(file.path(split_dir, "base_data.rds"))$outcome
  } else {
    readRDS(file.path(ranking_dir, "base_outcome.rds"))
  }
  pool <- split$pools$var2000
  standardised <- standardise_split(
    split$train_raw[, pool, drop = FALSE],
    split$test_raw[, pool, drop = FALSE]
  )
  ranking <- read.csv(file.path(ranking_dir, "ranking_r1_f1_dc_pf.csv"),
                      stringsAsFactors = FALSE)$gene
  panel <- head(ranking, 50L)
  evaluation_seed <- 420000L + 1000L + 100L + match(50L, panel_sizes)
  scores <- predict_with_ensemble(
    standardised$train[, panel, drop = FALSE], outcome[split$train_idx],
    standardised$test[, panel, drop = FALSE],
    random_seed = evaluation_seed
  )
  recomputed_auc <- bench_auc(outcome[split$test_idx], scores)
  written_auc <- result$AUC[
    result$repeat_idx == 1L & result$fold_idx == 1L &
      result$arm == "dc_pf" & result$k == 50L
  ]
  stopifnot(length(written_auc) == 1L,
            isTRUE(all.equal(recomputed_auc, written_auc, tolerance = 1e-12)))

  qa_rows[[length(qa_rows) + 1L]] <- data.frame(
    dataset = dataset,
    rows = nrow(result),
    methods = length(methods),
    duplicate_keys = anyDuplicated(keys),
    minimum_auc = min(result$AUC),
    maximum_auc = max(result$AUC),
    control_rankings_identical = all(control$identical),
    spot_check_written_auc = written_auc,
    spot_check_recomputed_auc = recomputed_auc,
    spot_check_absolute_difference = abs(written_auc - recomputed_auc)
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
qa <- do.call(rbind, qa_rows)
write.csv(qa, file.path(output_dir, "deterministic_validation.csv"),
          row.names = FALSE)
cat(sprintf("Validated %d datasets, %d rows, and %d unique methods per dataset.\n",
            nrow(qa), expected_rows, expected_methods))
print(qa, row.names = FALSE)
