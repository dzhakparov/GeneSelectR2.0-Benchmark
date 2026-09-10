#!/usr/bin/env Rscript

# Repeated matched-Random baseline for the seven older benchmark datasets.
# Random panels are sampled from the same outer-split var2000 pool used by all
# methods. The first three draws reproduce the baseline in eval_deterministic;
# the additional draws quantify Monte Carlo uncertainty in AUC-minus-Random.

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
stage <- if (length(args) >= 2L) args[[2L]] else "all"
n_draws <- if (length(args) >= 3L) as.integer(args[[3L]]) else 30L
budget <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 604800
validation_datasets_old <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
)
datasets <- c(validation_datasets_old, "imvigor210", "sosall")
stopifnot(
  dataset %in% datasets,
  stage %in% c("init", "evaluate", "assemble", "all"),
  is.finite(n_draws), n_draws >= 3L,
  is.finite(budget), budget > 0
)

started_at <- proc.time()[["elapsed"]]
over_budget <- function() proc.time()[["elapsed"]] - started_at > budget

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(withr)
  library(xgboost)
})
options(warn = 1)
source(file.path("redesign", "R", "bio_prior.R"))
source(file.path("redesign", "R", "evaluator.R"))
source(file.path("redesign", "R", "run_provenance.R"))

results_root <- redesign_results_root()
is_validation <- dataset %in% validation_datasets_old
split_dir <- if (is_validation) {
  file.path(results_root, "validation_benchmark", dataset)
} else {
  file.path(results_root, "grouped_benchmark", dataset)
}
dataset_dir <- if (is_validation) split_dir else {
  file.path(results_root, "full_recipe", dataset)
}
require_redesign_run(split_dir)
if (!identical(split_dir, dataset_dir)) require_redesign_run(dataset_dir)

panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
split_grid <- expand.grid(
  repeat_idx = 1:3, fold_idx = 1:5,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
split_override <- Sys.getenv("GENESELECTR_VALIDATION_SPLIT", "")
if (nzchar(split_override)) {
  parsed <- regmatches(
    split_override, regexec("^r([1-3])f([1-5])$", split_override)
  )[[1L]]
  if (length(parsed) != 3L) {
    stop("GENESELECTR_VALIDATION_SPLIT must have the form r1f1.",
         call. = FALSE)
  }
  split_grid <- data.frame(
    repeat_idx = as.integer(parsed[[2L]]),
    fold_idx = as.integer(parsed[[3L]])
  )
}

split_path <- function(repeat_idx, fold_idx) file.path(
  split_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
)
checkpoint_path <- function(repeat_idx, fold_idx) file.path(
  dataset_dir, sprintf(
    "older7_random_%ddraw_r%d_f%d.rds", n_draws, repeat_idx, fold_idx
  )
)
all_split_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    split_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
all_checkpoint_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    checkpoint_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
evaluation_path <- file.path(dataset_dir, "eval_deterministic.csv")
outcome_path <- file.path(
  dataset_dir, if (is_validation) "base_data.rds" else "base_outcome.rds"
)
draw_output <- file.path(
  dataset_dir, sprintf("older7_random_baseline_%d_draws.csv", n_draws)
)
split_output <- file.path(
  dataset_dir, sprintf("older7_random_baseline_%d_by_split.csv", n_draws)
)
summary_output <- file.path(
  dataset_dir, sprintf("older7_random_baseline_%d_summary.csv", n_draws)
)
qa_output <- file.path(
  dataset_dir, sprintf("older7_random_baseline_%d_QA.csv", n_draws)
)
if (!all(file.exists(c(all_split_files, evaluation_path, outcome_path)))) {
  stop("Saved outer splits or deterministic evaluation are missing.",
       call. = FALSE)
}

extension <- sprintf("older7_random_%ddraw_v1", n_draws)
prepare_redesign_extension(
  dataset_dir,
  extension = extension,
  source_files = c(
    "redesign/run_older7_random_sensitivity.R",
    file.path("redesign", "R", c(
      "bio_prior.R", "evaluator.R", "run_provenance.R"
    )),
    file.path(split_dir, "run_manifest.rds"),
    all_split_files, evaluation_path, outcome_path
  ),
  config = list(
    dataset = dataset,
    n_draws = n_draws,
    panel_sizes = panel_sizes,
    pool = "outer_split_var2000",
    panel_seed = "99 + 1000 * repeat + fold",
    model_seed = "420000 + 1000 * repeat + 100 * fold + panel_index",
    evaluator = "glmnet+xgboost+ranger",
    first_three_draws_match_saved_baseline = TRUE
  ),
  output_files = c(
    all_checkpoint_files, draw_output, split_output, summary_output, qa_output
  )
)
if (stage == "init") {
  cat(sprintf(
    "%s older-seven %d-draw Random extension initialized\n",
    dataset, n_draws
  ))
  quit(save = "no")
}

outcome <- if (is_validation) {
  readRDS(outcome_path)$outcome
} else {
  readRDS(outcome_path)
}
outcome <- droplevels(as.factor(outcome))
if (nlevels(outcome) != 2L) {
  stop("Saved outcome must have exactly two levels.", call. = FALSE)
}

saved_evaluation <- read.csv(evaluation_path, stringsAsFactors = FALSE)
saved_random <- saved_evaluation[
  saved_evaluation$arm == "Random",
  c("repeat_idx", "fold_idx", "k", "AUC"), drop = FALSE
]
if (nrow(saved_random) != 15L * length(panel_sizes) ||
    anyDuplicated(saved_random[c("repeat_idx", "fold_idx", "k")])) {
  stop("eval_deterministic.csv has an invalid saved Random baseline.",
       call. = FALSE)
}

if (stage %in% c("evaluate", "all")) {
  for (job_idx in seq_len(nrow(split_grid))) {
    repeat_idx <- split_grid$repeat_idx[[job_idx]]
    fold_idx <- split_grid$fold_idx[[job_idx]]
    output_path <- checkpoint_path(repeat_idx, fold_idx)
    if (file.exists(output_path)) next
    if (over_budget()) {
      cat("[budget] stop in repeated Random evaluation\n")
      quit(save = "no")
    }

    split <- readRDS(split_path(repeat_idx, fold_idx))
    pool <- split$pools$var2000
    standardized <- standardise_split(
      split$train_raw[, pool, drop = FALSE],
      split$test_raw[, pool, drop = FALSE]
    )
    y_train <- droplevels(outcome[split$train_idx])
    y_test <- droplevels(outcome[split$test_idx])

    rows <- lapply(seq_along(panel_sizes), function(panel_idx) {
      panel_size <- panel_sizes[[panel_idx]]
      evaluation_seed <- 420000L + repeat_idx * 1000L +
        fold_idx * 100L + panel_idx
      panel_seed <- 99L + 1000L * repeat_idx + fold_idx
      aucs <- random_panel_aucs(
        standardized$train, standardized$test,
        y_train, y_test, pool, panel_size,
        n_draws = n_draws, seed = panel_seed,
        model_seed = evaluation_seed
      )
      data.frame(
        dataset = dataset,
        repeat_idx = repeat_idx,
        fold_idx = fold_idx,
        k = panel_size,
        draw = seq_len(n_draws),
        AUC = aucs,
        n_panel = min(panel_size, length(pool)),
        panel_seed = panel_seed,
        evaluation_seed = evaluation_seed,
        evaluator_components = "glmnet+xgboost+ranger",
        stringsAsFactors = FALSE
      )
    })
    result <- do.call(rbind, rows)
    keys <- paste(result$k, result$draw)
    if (nrow(result) != length(panel_sizes) * n_draws ||
        anyDuplicated(keys) || any(!is.finite(result$AUC)) ||
        any(result$AUC < 0 | result$AUC > 1)) {
      stop("Repeated Random checkpoint failed validation.", call. = FALSE)
    }

    saved_here <- saved_random[
      saved_random$repeat_idx == repeat_idx &
        saved_random$fold_idx == fold_idx, , drop = FALSE
    ]
    first_three <- aggregate(
      AUC ~ k, result[result$draw <= 3L, , drop = FALSE], mean
    )
    comparison <- merge(
      first_three, saved_here[c("k", "AUC")], by = "k",
      suffixes = c("_recomputed", "_saved"), sort = FALSE
    )
    if (nrow(comparison) != length(panel_sizes) ||
        max(abs(comparison$AUC_recomputed - comparison$AUC_saved)) > 1e-12) {
      stop(sprintf(
        "The first three Random draws do not reproduce r%d f%d.",
        repeat_idx, fold_idx
      ), call. = FALSE)
    }
    saveRDS(result, output_path, version = 3)
    cat(sprintf(
      "[%s r%d f%d] %d Random evaluations complete\n",
      dataset, repeat_idx, fold_idx, nrow(result)
    ))
  }
}

if (stage %in% c("assemble", "all")) {
  if (!all(file.exists(all_checkpoint_files))) {
    missing <- sum(!file.exists(all_checkpoint_files))
    stop(sprintf("%d repeated-Random checkpoints are incomplete.", missing),
         call. = FALSE)
  }
  result <- do.call(rbind, lapply(all_checkpoint_files, readRDS))
  keys <- paste(
    result$repeat_idx, result$fold_idx, result$k, result$draw
  )
  expected_rows <- 15L * length(panel_sizes) * n_draws
  if (nrow(result) != expected_rows || anyDuplicated(keys) ||
      any(!is.finite(result$AUC)) || any(result$AUC < 0 | result$AUC > 1)) {
    stop("Assembled repeated Random results failed validation.",
         call. = FALSE)
  }

  by_split <- aggregate(
    AUC ~ dataset + repeat_idx + fold_idx + k, result,
    function(values) c(
      mean = mean(values), sd = stats::sd(values),
      mc_se = stats::sd(values) / sqrt(length(values)),
      q025 = as.numeric(stats::quantile(values, 0.025)),
      q975 = as.numeric(stats::quantile(values, 0.975))
    )
  )
  value_matrix <- by_split$AUC
  by_split$AUC <- NULL
  by_split <- cbind(by_split, value_matrix)

  summary_rows <- lapply(panel_sizes, function(panel_size) {
    split_values <- by_split[by_split$k == panel_size, , drop = FALSE]
    data.frame(
      dataset = dataset,
      k = panel_size,
      mean_random_AUC = mean(split_values$mean),
      sd_across_splits = stats::sd(split_values$mean),
      se_across_splits = stats::sd(split_values$mean) /
        sqrt(nrow(split_values)),
      mean_within_split_mc_se = mean(split_values$mc_se),
      n_splits = nrow(split_values),
      n_draws_per_split = n_draws
    )
  })
  summary_table <- do.call(rbind, summary_rows)

  first_three <- aggregate(
    AUC ~ repeat_idx + fold_idx + k,
    result[result$draw <= 3L, , drop = FALSE], mean
  )
  qa <- merge(
    first_three, saved_random,
    by = c("repeat_idx", "fold_idx", "k"),
    suffixes = c("_recomputed", "_saved"), sort = FALSE
  )
  qa$absolute_difference <- abs(qa$AUC_recomputed - qa$AUC_saved)
  qa$passed <- qa$absolute_difference <= 1e-12
  if (nrow(qa) != 15L * length(panel_sizes) || !all(qa$passed)) {
    stop("Saved three-draw Random baseline was not reproduced.",
         call. = FALSE)
  }

  write.csv(result, draw_output, row.names = FALSE)
  write.csv(by_split, split_output, row.names = FALSE)
  write.csv(summary_table, summary_output, row.names = FALSE)
  write.csv(qa, qa_output, row.names = FALSE)
  cat(sprintf(
    "%s: wrote %d repeated Random evaluations; saved baseline reproduced\n",
    dataset, nrow(result)
  ))
}
