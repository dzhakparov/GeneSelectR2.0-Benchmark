#!/usr/bin/env Rscript
# ==============================================================================
#  Deterministic re-evaluation of every saved benchmark ranking.
#
#  Feature selectors are not refitted. Each saved ranking is applied to the
#  original outer split and evaluated with explicit, matched inner folds and
#  model seeds. The dc_pf_ctrl duplicate is checked against dc_pf and excluded
#  from the method table.
#
#  Usage:
#    Rscript redesign/reevaluate_saved_rankings.R <dataset> [result_root]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
results_root <- if (length(args) >= 2L) args[[2L]] else
  file.path("redesign", "results_corrected")

validation_datasets <- c("GSE65682", "GSE69683", "GSE13355", "GSE107994",
                         "GSE101794")
stopifnot(dataset %in% c(validation_datasets, "imvigor210", "sosall"))

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(withr)
  library(xgboost)
})
source(file.path("redesign", "R", "bio_prior.R"))
source(file.path("redesign", "R", "evaluator.R"))

is_validation <- dataset %in% validation_datasets
split_dir <- if (is_validation) {
  file.path(results_root, "validation_benchmark", dataset)
} else {
  file.path(results_root, "grouped_benchmark", dataset)
}
ranking_dir <- if (is_validation) split_dir else
  file.path(results_root, "full_recipe", dataset)

outcome <- if (is_validation) {
  readRDS(file.path(split_dir, "base_data.rds"))$outcome
} else {
  readRDS(file.path(ranking_dir, "base_outcome.rds"))
}

panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
methods <- c(
  "Boruta", "DGE", "ElasticNet", "GS_full_grouped",
  "GS_full_ungrouped", "LASSO", "mRMR", "RF_importance",
  "predfirst", "predfirst_raw", "ens_rank", "ens_rank3", "kswitch",
  "soft_gs", "soft_pf", "dc_gs", "dc_pf", "cb_gs", "cb_pf", "hs",
  "ens2", "ens3", "cb_gs_w2", "hs_stab", "cbgs_prune", "dcpf_prune",
  "cbgs_adapt"
)

# Preserve method-selection metadata from the original tables. These fields
# describe the saved rankings and remain valid after deterministic prediction.
metadata <- list()
for (path in c(file.path(ranking_dir, "eval_results.csv"),
               file.path(ranking_dir, "eval_adapt.csv"))) {
  if (!file.exists(path)) next
  table <- read.csv(path, stringsAsFactors = FALSE)
  table$arm[table$arm == "full_grouped"] <- "GS_full_grouped"
  table$arm[table$arm == "full_ungrouped"] <- "GS_full_ungrouped"
  for (index in seq_len(nrow(table))) {
    key <- paste(table$repeat_idx[index], table$fold_idx[index],
                 table$arm[index], table$k[index])
    metadata[[key]] <- table[index, , drop = FALSE]
  }
}

metadata_value <- function(repeat_idx, fold_idx, arm, k, field) {
  row <- metadata[[paste(repeat_idx, fold_idx, arm, k)]]
  if (is.null(row) || !(field %in% names(row))) return(NA)
  row[[field]][1L]
}

output_path <- file.path(ranking_dir, "eval_deterministic.csv")
control_path <- file.path(ranking_dir, "reproducibility_control.csv")
manifest_path <- file.path(ranking_dir, "eval_deterministic_manifest.csv")
version_path <- file.path(ranking_dir, "eval_deterministic.version")
evaluation_version <- "2026-08-29-deterministic-v1"
manifest <- data.frame(
  evaluation_version = evaluation_version,
  evaluator_md5 = unname(tools::md5sum(file.path("redesign", "R",
                                                  "evaluator.R"))),
  script_md5 = unname(tools::md5sum(
    file.path("redesign", "reevaluate_saved_rankings.R"))),
  methods = length(methods),
  splits = 15L,
  panel_sizes = paste(panel_sizes, collapse = ";")
)
if (file.exists(output_path)) {
  if (!file.exists(manifest_path) || !file.exists(version_path)) {
    stop("deterministic checkpoint has no provenance manifest: ", output_path,
         call. = FALSE)
  }
  observed <- read.csv(manifest_path, stringsAsFactors = FALSE)
  if (!identical(observed, manifest) ||
      !identical(readLines(version_path, warn = FALSE), evaluation_version)) {
    stop("deterministic checkpoint provenance does not match current code: ",
         output_path, call. = FALSE)
  }
} else {
  write.csv(manifest, manifest_path, row.names = FALSE)
  writeLines(evaluation_version, version_path)
}
done <- if (file.exists(output_path)) read.csv(output_path) else data.frame()
done_keys <- if (nrow(done)) {
  paste(done$repeat_idx, done$fold_idx, done$arm, done$k)
} else character(0)

read_ranking <- function(repeat_idx, fold_idx, arm, k) {
  file_arm <- arm
  if (!is_validation) {
    file_arm <- switch(arm,
      GS_full_grouped = "full_grouped",
      GS_full_ungrouped = "full_ungrouped",
      arm
    )
  }
  # kswitch used prediction-first at k <= 20 and grouped GeneSelectR above
  # that threshold. Its saved k=50 ranking cannot represent both branches.
  if (arm == "kswitch") {
    file_arm <- if (k <= 20L) "predfirst_raw" else
      if (is_validation) "GS_full_grouped" else "full_grouped"
  }
  path <- file.path(ranking_dir, sprintf(
    "ranking_r%d_f%d_%s.csv", repeat_idx, fold_idx, file_arm))
  if (!file.exists(path)) stop("missing ranking: ", path, call. = FALSE)
  ranking <- read.csv(path, stringsAsFactors = FALSE)$gene
  if (!length(ranking) && arm %in% c("GS_full_grouped", "kswitch")) {
    fallback_arm <- if (is_validation) "GS_full_ungrouped" else
      "full_ungrouped"
    fallback_path <- file.path(ranking_dir, sprintf(
      "ranking_r%d_f%d_%s.csv", repeat_idx, fold_idx, fallback_arm))
    ranking <- read.csv(fallback_path, stringsAsFactors = FALSE)$gene
  }
  if (!length(ranking) || anyDuplicated(ranking)) {
    stop("invalid ranking: ", path, call. = FALSE)
  }
  ranking
}

append_rows <- function(rows, path) {
  if (!length(rows)) return(invisible(NULL))
  write.table(do.call(rbind, rows), path, append = file.exists(path),
              sep = ",", row.names = FALSE, col.names = !file.exists(path))
}

split_files <- sort(list.files(split_dir, pattern = "^split_r.*[.]rds$",
                               full.names = TRUE))
stopifnot(length(split_files) == 15L)

for (split_path in split_files) {
  parsed <- regmatches(basename(split_path),
    regexec("split_r([0-9]+)_f([0-9]+)[.]rds", basename(split_path)))[[1L]]
  repeat_idx <- as.integer(parsed[[2L]])
  fold_idx <- as.integer(parsed[[3L]])
  split <- readRDS(split_path)
  y_train <- outcome[split$train_idx]
  y_test <- outcome[split$test_idx]
  pool <- split$pools$var2000
  standardised <- standardise_split(
    split$train_raw[, pool, drop = FALSE],
    split$test_raw[, pool, drop = FALSE]
  )

  rows <- list()
  for (k in panel_sizes) {
    evaluation_seed <- 420000L + repeat_idx * 1000L + fold_idx * 100L +
      match(k, panel_sizes)
    for (arm in methods) {
      key <- paste(repeat_idx, fold_idx, arm, k)
      if (key %in% done_keys) next
      ranking <- read_ranking(repeat_idx, fold_idx, arm, k)
      panel <- head(ranking, min(k, length(ranking)))
      scores <- predict_with_ensemble(
        standardised$train[, panel, drop = FALSE], y_train,
        standardised$test[, panel, drop = FALSE],
        random_seed = evaluation_seed
      )
      rows[[length(rows) + 1L]] <- data.frame(
        repeat_idx = repeat_idx, fold_idx = fold_idx, arm = arm, k = k,
        AUC = bench_auc(y_test, scores), n_panel = length(panel),
        evaluation_seed = evaluation_seed,
        evaluator_components = paste(attr(scores, "components_used"),
                                     collapse = "+"),
        gate_fallback = metadata_value(repeat_idx, fold_idx, arm, k,
                                       "gate_fallback"),
        w_chosen = metadata_value(repeat_idx, fold_idx, arm, k, "w_chosen")
      )
      done_keys <- c(done_keys, key)
    }

    random_key <- paste(repeat_idx, fold_idx, "Random", k)
    if (!(random_key %in% done_keys)) {
      aucs <- random_panel_aucs(
        standardised$train, standardised$test, y_train, y_test, pool, k,
        n_draws = 3L, seed = 99L + 1000L * repeat_idx + fold_idx,
        model_seed = evaluation_seed
      )
      rows[[length(rows) + 1L]] <- data.frame(
        repeat_idx = repeat_idx, fold_idx = fold_idx, arm = "Random", k = k,
        AUC = mean(aucs), n_panel = min(k, length(pool)),
        evaluation_seed = evaluation_seed,
        evaluator_components = "glmnet+xgboost+ranger",
        gate_fallback = metadata_value(repeat_idx, fold_idx, "Random", k,
                                       "gate_fallback"),
        w_chosen = NA
      )
      done_keys <- c(done_keys, random_key)
    }
  }
  append_rows(rows, output_path)
  cat(sprintf("[%s r%d f%d] %d deterministic rows complete\n", dataset,
              repeat_idx, fold_idx, length(rows)))
}

# Verify the removed control by file content. dc_pf_ctrl must be an exact copy
# of dc_pf for every split. This check remains separate from method ranking.
control_rows <- list()
for (repeat_idx in 1:3) {
  for (fold_idx in 1:5) {
    main_path <- file.path(ranking_dir, sprintf(
      "ranking_r%d_f%d_dc_pf.csv", repeat_idx, fold_idx))
    control_copy <- file.path(ranking_dir, sprintf(
      "ranking_r%d_f%d_dc_pf_ctrl.csv", repeat_idx, fold_idx))
    control_rows[[length(control_rows) + 1L]] <- data.frame(
      repeat_idx = repeat_idx,
      fold_idx = fold_idx,
      dc_pf_md5 = unname(tools::md5sum(main_path)),
      dc_pf_ctrl_md5 = unname(tools::md5sum(control_copy)),
      identical = identical(readLines(main_path), readLines(control_copy))
    )
  }
}
write.csv(do.call(rbind, control_rows), control_path, row.names = FALSE)
if (!all(do.call(rbind, control_rows)$identical)) {
  stop("dc_pf_ctrl differs from dc_pf", call. = FALSE)
}

expected_rows <- (length(methods) + 1L) * length(panel_sizes) * 15L
result <- read.csv(output_path)
keys <- paste(result$repeat_idx, result$fold_idx, result$arm, result$k)
stopifnot(nrow(result) == expected_rows, !anyDuplicated(keys),
          all(is.finite(result$AUC)), all(result$AUC >= 0 & result$AUC <= 1))
cat(sprintf("[%s] complete: %d rows, %d methods plus Random\n", dataset,
            nrow(result), length(methods)))
