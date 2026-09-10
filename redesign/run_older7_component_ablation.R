#!/usr/bin/env Rscript

# Component ablation for the current GeneSelectR configuration on the seven
# older datasets. Rankings are reconstructed from each saved outer-training
# fit. No selector is refitted. Each ranking is evaluated on the untouched
# outer test fold with the benchmark's deterministic three-model ensemble.

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
stage <- if (length(args) >= 2L) args[[2L]] else "all"
budget <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 604800
validation_datasets_old <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
)
datasets <- c(validation_datasets_old, "imvigor210", "sosall")
stopifnot(
  dataset %in% datasets,
  stage %in% c("init", "evaluate", "assemble", "all"),
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
source(file.path("package", "GeneSelectR", "R", "utils.R"))
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
variants <- c(
  "GS_recurrence_only", "GS_SHAP_only", "GS_MI_only",
  "GS_SHAPxMI", "GS_current"
)
saved_method <- if (is_validation) "GS_full_ungrouped" else "full_ungrouped"
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
ranking_path <- function(repeat_idx, fold_idx) file.path(
  dataset_dir, sprintf(
    "ranking_r%d_f%d_%s.csv", repeat_idx, fold_idx, saved_method
  )
)
meta_path <- function(repeat_idx, fold_idx) file.path(
  dataset_dir, sprintf(
    "ranking_r%d_f%d_%s_meta.csv", repeat_idx, fold_idx, saved_method
  )
)
read_selected_alpha <- function(repeat_idx, fold_idx) {
  metadata <- read.csv(
    meta_path(repeat_idx, fold_idx), stringsAsFactors = FALSE
  )
  if (!identical(names(metadata), c("key", "value")) ||
      sum(metadata$key == "alpha") != 1L) {
    stop("Invalid GeneSelectR ranking metadata.", call. = FALSE)
  }
  alpha <- suppressWarnings(as.numeric(metadata$value[metadata$key == "alpha"]))
  if (length(alpha) != 1L || !alpha %in% c(0.5, 1.0)) {
    stop("Selected alpha must be 0.5 or 1.0.", call. = FALSE)
  }
  alpha
}
alpha_tag <- function(alpha) {
  if (identical(as.numeric(alpha), 0.5)) "0p5" else "1"
}
fit_path <- function(repeat_idx, fold_idx) {
  alpha <- read_selected_alpha(repeat_idx, fold_idx)
  file.path(dataset_dir, sprintf(
    "fit_r%d_f%d_%s_a%s.rds",
    repeat_idx, fold_idx, saved_method, alpha_tag(alpha)
  ))
}
checkpoint_path <- function(repeat_idx, fold_idx) file.path(
  dataset_dir, sprintf(
    "older7_component_ablation_r%d_f%d.rds", repeat_idx, fold_idx
  )
)
all_split_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    split_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
all_ranking_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    ranking_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
all_meta_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    meta_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
if (!all(file.exists(c(all_split_files, all_ranking_files, all_meta_files)))) {
  stop("Saved splits, rankings, or ranking metadata are missing.",
       call. = FALSE)
}
all_fit_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    fit_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
evaluation_path <- file.path(dataset_dir, "eval_deterministic.csv")
outcome_path <- file.path(
  dataset_dir, if (is_validation) "base_data.rds" else "base_outcome.rds"
)
all_checkpoint_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    checkpoint_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
evaluation_output <- file.path(
  dataset_dir, "older7_component_ablation_evaluation.csv"
)
panel_output <- file.path(
  dataset_dir, "older7_component_ablation_panels.csv"
)
summary_output <- file.path(
  dataset_dir, "older7_component_ablation_summary.csv"
)
qa_output <- file.path(dataset_dir, "older7_component_ablation_QA.csv")
if (!all(file.exists(c(all_fit_files, evaluation_path, outcome_path)))) {
  stop("Selected GeneSelectR fits or deterministic evaluation are missing.",
       call. = FALSE)
}

extension <- "older7_component_ablation_v1"
prepare_redesign_extension(
  dataset_dir,
  extension = extension,
  source_files = c(
    "redesign/run_older7_component_ablation.R",
    file.path("package", "GeneSelectR", "R", "utils.R"),
    file.path("redesign", "R", c(
      "bio_prior.R", "evaluator.R", "run_provenance.R"
    )),
    file.path(split_dir, "run_manifest.rds"),
    all_split_files, all_ranking_files, all_meta_files, all_fit_files,
    evaluation_path, outcome_path
  ),
  config = list(
    dataset = dataset,
    variants = variants,
    panel_sizes = panel_sizes,
    fit_source = "saved_outer_training_fit",
    recurrence_score = "calibrated_selection_frequency",
    SHAP_score = "held_out_instance_frequency",
    MI_score = "saved_mutual_information_percentile",
    SHAPxMI_score = "calibrated_combined_utility",
    current_score = "saved_stability_utility_geometric_mean",
    evaluator = "glmnet+xgboost+ranger",
    model_seed = "420000 + 1000 * repeat + 100 * fold + panel_index"
  ),
  output_files = c(
    all_checkpoint_files, evaluation_output, panel_output,
    summary_output, qa_output
  )
)
if (stage == "init") {
  cat(sprintf("%s older-seven component ablation initialized\n", dataset))
  quit(save = "no")
}

outcome <- if (is_validation) {
  readRDS(outcome_path)$outcome
} else {
  readRDS(outcome_path)
}
outcome <- droplevels(as.factor(outcome))
saved_evaluation <- read.csv(evaluation_path, stringsAsFactors = FALSE)

reconstruct_instance_frequency <- function(fit, genes) {
  shap_matrix <- fit$instance_importance
  if (is.null(shap_matrix) || is.null(colnames(shap_matrix)) ||
      anyDuplicated(colnames(shap_matrix)) ||
      !all(colnames(shap_matrix) %in% genes) ||
      any(!is.finite(shap_matrix)) || any(shap_matrix < 0)) {
    stop("Saved instance-level SHAP matrix is invalid.", call. = FALSE)
  }
  important <- matrix(
    FALSE, nrow = nrow(shap_matrix), ncol = ncol(shap_matrix),
    dimnames = dimnames(shap_matrix)
  )
  for (sample_idx in seq_len(nrow(shap_matrix))) {
    values <- shap_matrix[sample_idx, ]
    positive <- values[values > 0]
    if (length(positive) == 0L) next
    threshold <- as.numeric(stats::quantile(
      positive, 0.75, na.rm = TRUE
    ))
    important[sample_idx, ] <- values >= threshold
  }
  frequency <- stats::setNames(rep(0, length(genes)), genes)
  frequency[colnames(shap_matrix)] <- colMeans(important)
  frequency
}

make_rankings <- function(fit, saved_ranking) {
  scores <- fit$gene_scores
  required <- c(
    "gene", "raw_score", "pi_raw", "pi_scored", "u", "u_mi", "u_scored"
  )
  if (!all(required %in% names(scores)) || anyDuplicated(scores$gene) ||
      !setequal(scores$gene, saved_ranking$gene) ||
      !identical(as.character(scores$gene), as.character(saved_ranking$gene))) {
    stop("Saved fit and ranking are inconsistent.", call. = FALSE)
  }
  numeric_columns <- setdiff(required, "gene")
  if (any(!is.finite(as.matrix(scores[numeric_columns])))) {
    stop("GeneSelectR component scores contain non-finite values.",
         call. = FALSE)
  }

  shap_frequency <- reconstruct_instance_frequency(fit, scores$gene)
  shap_values <- unname(shap_frequency[scores$gene])
  shap_percentile <- percentile01(shap_values)
  reconstructed_u <- percentile01(sqrt(
    shap_percentile * scores$u_mi + 1e-10
  ))
  if (max(abs(reconstructed_u - scores$u)) > 1e-12) {
    stop("Held-out SHAP reconstruction does not reproduce saved utility.",
         call. = FALSE)
  }

  list(
    GS_recurrence_only = scores$gene[order(
      -scores$pi_scored, -scores$pi_raw, scores$gene, method = "radix"
    )],
    GS_SHAP_only = scores$gene[order(
      -shap_values, scores$gene, method = "radix"
    )],
    GS_MI_only = scores$gene[order(
      -scores$u_mi, scores$gene, method = "radix"
    )],
    GS_SHAPxMI = scores$gene[order(
      -scores$u_scored, -scores$u, scores$gene, method = "radix"
    )],
    GS_current = as.character(saved_ranking$gene)
  )
}

if (stage %in% c("evaluate", "all")) {
  for (job_idx in seq_len(nrow(split_grid))) {
    repeat_idx <- split_grid$repeat_idx[[job_idx]]
    fold_idx <- split_grid$fold_idx[[job_idx]]
    output_path <- checkpoint_path(repeat_idx, fold_idx)
    if (file.exists(output_path)) next
    if (over_budget()) {
      cat("[budget] stop in component ablation\n")
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
    wrapper <- readRDS(fit_path(repeat_idx, fold_idx))
    fit <- wrapper$fit
    selected_alpha <- read_selected_alpha(repeat_idx, fold_idx)
    if (is.null(fit) || !isTRUE(all.equal(wrapper$alpha, selected_alpha)) ||
        !identical(fit$parameters$components, c("stability", "utility")) ||
        !identical(fit$parameters$gate_method, "none") ||
        !identical(fit$parameters$calibration_mode, "evidence_ratio") ||
        !identical(fit$parameters$utility_method, "instance_shap")) {
      stop("Saved fit is not the current benchmark configuration.",
           call. = FALSE)
    }
    saved_ranking <- read.csv(
      ranking_path(repeat_idx, fold_idx), stringsAsFactors = FALSE
    )
    rankings <- make_rankings(fit, saved_ranking)
    if (!identical(names(rankings), variants) ||
        any(vapply(rankings, anyDuplicated, integer(1)) > 0L) ||
        any(!vapply(rankings, function(genes) {
          setequal(genes, pool)
        }, logical(1)))) {
      stop("A component ranking does not match the candidate pool.",
           call. = FALSE)
    }

    saved_random <- saved_evaluation[
      saved_evaluation$repeat_idx == repeat_idx &
        saved_evaluation$fold_idx == fold_idx &
        saved_evaluation$arm == "Random",
      c("k", "AUC"), drop = FALSE
    ]
    saved_current <- saved_evaluation[
      saved_evaluation$repeat_idx == repeat_idx &
        saved_evaluation$fold_idx == fold_idx &
        saved_evaluation$arm == "GS_full_ungrouped",
      c("k", "AUC"), drop = FALSE
    ]
    if (nrow(saved_random) != length(panel_sizes) ||
        nrow(saved_current) != length(panel_sizes)) {
      stop("Saved deterministic evaluation is incomplete.", call. = FALSE)
    }

    evaluation_rows <- list()
    panel_rows <- list()
    for (variant in variants) {
      ranking <- rankings[[variant]]
      panel_rows[[length(panel_rows) + 1L]] <- data.frame(
        dataset = dataset,
        repeat_idx = repeat_idx,
        fold_idx = fold_idx,
        variant = variant,
        rank = seq_len(max(panel_sizes)),
        gene = ranking[seq_len(max(panel_sizes))],
        stringsAsFactors = FALSE
      )
      for (panel_idx in seq_along(panel_sizes)) {
        panel_size <- panel_sizes[[panel_idx]]
        panel <- ranking[seq_len(min(panel_size, length(ranking)))]
        evaluation_seed <- 420000L + repeat_idx * 1000L +
          fold_idx * 100L + panel_idx
        predictions <- predict_with_ensemble(
          standardized$train[, panel, drop = FALSE], y_train,
          standardized$test[, panel, drop = FALSE],
          random_seed = evaluation_seed
        )
        random_auc <- saved_random$AUC[match(panel_size, saved_random$k)]
        auc <- bench_auc(y_test, predictions)
        evaluation_rows[[length(evaluation_rows) + 1L]] <- data.frame(
          dataset = dataset,
          repeat_idx = repeat_idx,
          fold_idx = fold_idx,
          variant = variant,
          k = panel_size,
          AUC = auc,
          random_AUC_3draw = random_auc,
          delta_random_3draw = auc - random_auc,
          n_panel = length(panel),
          selected_alpha = selected_alpha,
          evaluation_seed = evaluation_seed,
          evaluator_components = "glmnet+xgboost+ranger",
          stringsAsFactors = FALSE
        )
      }
    }
    evaluation <- do.call(rbind, evaluation_rows)
    panels <- do.call(rbind, panel_rows)
    keys <- paste(evaluation$variant, evaluation$k)
    current <- evaluation[evaluation$variant == "GS_current", ]
    current <- merge(
      current, saved_current, by = "k", suffixes = c("", "_saved"),
      sort = FALSE
    )
    current$absolute_difference <- abs(current$AUC - current$AUC_saved)
    qa <- data.frame(
      dataset = dataset,
      repeat_idx = repeat_idx,
      fold_idx = fold_idx,
      k = current$k,
      recomputed_AUC = current$AUC,
      saved_AUC = current$AUC_saved,
      absolute_difference = current$absolute_difference,
      passed = current$absolute_difference <= 1e-12
    )
    if (nrow(evaluation) != length(variants) * length(panel_sizes) ||
        anyDuplicated(keys) || any(!is.finite(evaluation$AUC)) ||
        any(evaluation$AUC < 0 | evaluation$AUC > 1) || !all(qa$passed)) {
      stop("Component-ablation checkpoint failed validation.",
           call. = FALSE)
    }
    saveRDS(
      list(evaluation = evaluation, panels = panels, qa = qa),
      output_path, version = 3
    )
    cat(sprintf("[%s r%d f%d] component ablation complete\n",
                dataset, repeat_idx, fold_idx))
  }
}

if (stage %in% c("assemble", "all")) {
  if (!all(file.exists(all_checkpoint_files))) {
    missing <- sum(!file.exists(all_checkpoint_files))
    stop(sprintf("%d component-ablation checkpoints are incomplete.", missing),
         call. = FALSE)
  }
  checkpoints <- lapply(all_checkpoint_files, readRDS)
  evaluation <- do.call(rbind, lapply(checkpoints, `[[`, "evaluation"))
  panels <- do.call(rbind, lapply(checkpoints, `[[`, "panels"))
  qa <- do.call(rbind, lapply(checkpoints, `[[`, "qa"))
  evaluation_keys <- paste(
    evaluation$repeat_idx, evaluation$fold_idx,
    evaluation$variant, evaluation$k
  )
  panel_keys <- paste(
    panels$repeat_idx, panels$fold_idx, panels$variant, panels$rank
  )
  if (nrow(evaluation) != 15L * length(variants) * length(panel_sizes) ||
      anyDuplicated(evaluation_keys) ||
      nrow(panels) != 15L * length(variants) * max(panel_sizes) ||
      anyDuplicated(panel_keys) || !all(qa$passed)) {
    stop("Assembled component ablation failed validation.", call. = FALSE)
  }

  summary_rows <- lapply(variants, function(variant) {
    do.call(rbind, lapply(panel_sizes, function(panel_size) {
      values <- evaluation[
        evaluation$variant == variant & evaluation$k == panel_size,
        , drop = FALSE
      ]
      data.frame(
        dataset = dataset,
        variant = variant,
        k = panel_size,
        mean_AUC = mean(values$AUC),
        sd_AUC = stats::sd(values$AUC),
        mean_random_AUC_3draw = mean(values$random_AUC_3draw),
        mean_delta_random_3draw = mean(values$delta_random_3draw),
        sd_delta_random_3draw = stats::sd(values$delta_random_3draw),
        n_splits = nrow(values)
      )
    }))
  })
  summary_table <- do.call(rbind, summary_rows)
  write.csv(evaluation, evaluation_output, row.names = FALSE)
  write.csv(panels, panel_output, row.names = FALSE)
  write.csv(summary_table, summary_output, row.names = FALSE)
  write.csv(qa, qa_output, row.names = FALSE)
  cat(sprintf(
    "%s: wrote %d component-ablation evaluations; current ranking reproduced\n",
    dataset, nrow(evaluation)
  ))
}
