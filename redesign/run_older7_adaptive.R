#!/usr/bin/env Rscript

# Leakage-safe panel-size selection for the seven datasets in the corrected
# benchmark. A selector is refitted in each inner-training fold. Panel size is
# selected by the one-standard-error rule and evaluated with the saved,
# deterministic outer-fold ranking and matched Random result.

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
stage <- if (length(args) >= 2L) args[[2L]] else "all"
budget <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 604800
validation_datasets_old <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
)
datasets <- c(validation_datasets_old, "imvigor210", "sosall")
stopifnot(dataset %in% datasets,
          stage %in% c("init", "inner", "assemble", "all"))

started_at <- proc.time()[["elapsed"]]
over_budget <- function() proc.time()[["elapsed"]] - started_at > budget

suppressPackageStartupMessages({
  library(Boruta)
  library(edgeR)
  library(glmnet)
  library(mRMRe)
  library(ranger)
  library(withr)
  library(xgboost)
})
options(warn = 1)
for (path in list.files("package/GeneSelectR/R", full.names = TRUE)) {
  source(path)
}
for (file_name in c(
    "bio_prior.R", "competitors.R", "evaluator.R", "imvigor210_data.R",
    "run_provenance.R", "sosall_data.R", "validation_data.R")) {
  source(file.path("redesign", "R", file_name))
}

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

methods <- c(
  "GS_full_ungrouped", "DGE", "LASSO", "ElasticNet", "mRMR", "Boruta",
  "RF_importance"
)
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
alpha_grid <- c(0.5, 1.0)
inner_folds_n <- 3L
extension <- "older7_nested_panel_size_v2"
outer_output <- file.path(dataset_dir, "older7_adaptive_nested_results.csv")
curve_output <- file.path(dataset_dir, "older7_adaptive_inner_curves.csv")
fold_output <- file.path(dataset_dir, "older7_adaptive_inner_fold_results.csv")
source_files <- c(
  "redesign/run_older7_adaptive.R",
  list.files("package/GeneSelectR/R", full.names = TRUE),
  file.path("redesign", "R", c(
    "bio_prior.R", "competitors.R", "evaluator.R", "imvigor210_data.R",
    "run_provenance.R", "sosall_data.R", "validation_data.R"
  )),
  file.path(split_dir, "run_manifest.rds")
)
prepare_redesign_extension(
  dataset_dir,
  extension = extension,
  source_files = source_files,
  config = list(
    dataset = dataset,
    methods = methods,
    panel_sizes = panel_sizes,
    inner_folds = inner_folds_n,
    panel_rule = "smallest_within_one_SE_of_best",
    pool_rule = "inner_training_variance_top_2000",
    GS_B = 50L,
    GS_subsample_scheme = "kfold_5",
    GS_calibration_permutations = 20L,
    GS_calibration_null_B = 20L,
    small_class_cv = "reduce_and_stratify_when_minority_count_below_5",
    competitor_cv = "class_stratified_fold_id",
    inner_grouping = "paired_patients_kept_together_when_spec_paired",
    alpha_grid = alpha_grid,
    evaluator = "glmnet+xgboost+ranger"
  ),
  output_files = c(outer_output, curve_output, fold_output)
)
if (stage == "init") {
  cat(sprintf("%s older-seven nested extension initialized\n", dataset))
  quit(save = "no")
}

base <- if (is_validation) {
  readRDS(file.path(dataset_dir, "base_data.rds"))
} else if (dataset == "imvigor210") {
  load_imvigor210()
} else {
  load_sosall()
}
outcome <- if (is_validation) base$outcome else base$outcome
saved_outcome <- if (!is_validation) {
  readRDS(file.path(dataset_dir, "base_outcome.rds"))
} else {
  outcome
}
if (!identical(as.character(outcome), as.character(saved_outcome))) {
  stop("Reconstructed outcome differs from the saved benchmark outcome.",
       call. = FALSE)
}

preprocess_inner_split <- function(training_indices, validation_indices) {
  if (is_validation) {
    if (base$spec$scale == "counts") {
      normalized <- normalise_count_split(
        t(base$expr), training_indices, validation_indices
      )
      train_raw <- normalized$train
      validation_raw <- normalized$test
    } else {
      train_raw <- base$expr[training_indices, , drop = FALSE]
      validation_raw <- base$expr[validation_indices, , drop = FALSE]
    }
    residualized <- residualise_split_generic(
      train_raw, validation_raw,
      base$confounders[training_indices, , drop = FALSE],
      base$confounders[validation_indices, , drop = FALSE]
    )
  } else if (dataset == "imvigor210") {
    residualized <- normalise_count_split(
      base$raw_counts, training_indices, validation_indices
    )
  } else {
    residualized <- sosall_residualise_split(
      base$expression_matrix, base$metadata,
      training_indices, validation_indices
    )
  }
  pool <- make_pools(
    residualized$train, character(0), top_var = 2000
  )$var2000
  standardized <- standardise_split(
    residualized$train[, pool, drop = FALSE],
    residualized$test[, pool, drop = FALSE]
  )
  list(train = standardized$train, validation = standardized$test,
       pool = pool)
}

fit_gs_alpha <- function(train_features, train_labels, alpha_value) {
  fit <- geneselectr2_fit(
    train_features, train_labels,
    gate_method = "none",
    B = 50,
    subsample_scheme = "kfold",
    subsample_k_folds = 5,
    utility_method = "instance_shap",
    components = c("stability", "utility"),
    score_formula = "geometric",
    calibration_mode = "evidence_ratio",
    calibration_n_permutations = 20,
    calibration_null_B = 20,
    alpha = alpha_value,
    n_cores = 1,
    random_seed = 42,
    use_cache = TRUE,
    verbose = FALSE
  )
  list(gene_scores = fit$gene_scores,
       internal_oob_auc = fit$cv_results$mean_auc,
       alpha = alpha_value)
}

fit_competitor <- function(method, train_features, train_labels, seed) {
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

method_seed <- function(repeat_idx, fold_idx, inner_idx, method) {
  910000L + repeat_idx * 10000L + fold_idx * 1000L +
    inner_idx * 100L + match(method, methods)
}
evaluation_seed <- function(repeat_idx, fold_idx, inner_idx, panel_size) {
  920000L + repeat_idx * 10000L + fold_idx * 1000L +
    inner_idx * 100L + match(panel_size, panel_sizes)
}
result_path <- function(repeat_idx, fold_idx, inner_idx, method) {
  file.path(dataset_dir, sprintf(
    "older7_adaptive_inner_r%d_f%d_i%d_%s.rds",
    repeat_idx, fold_idx, inner_idx, method
  ))
}
alpha_path <- function(repeat_idx, fold_idx, inner_idx, alpha_value) {
  file.path(dataset_dir, sprintf(
    "older7_adaptive_gsfit_r%d_f%d_i%d_a%s.rds",
    repeat_idx, fold_idx, inner_idx,
    gsub("[.]", "p", as.character(alpha_value))
  ))
}

if (stage %in% c("inner", "all")) {
  for (repeat_idx in 1:3) {
    for (fold_idx in 1:5) {
      outer_split <- readRDS(file.path(
        split_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
      ))
      outer_indices <- outer_split$train_idx
      outer_labels <- outcome[outer_indices]
      inner_seed <- 930000L + repeat_idx * 100L + fold_idx
      inner_folds <- if (is_validation && isTRUE(base$spec$paired)) {
        # GSE13355 has one lesional and one uninvolved sample per patient.
        # Both samples must remain in the same inner fold.
        make_paired_folds(
          base$groups[outer_indices], inner_folds_n, seed = inner_seed
        )
      } else {
        make_stratified_folds(
          outer_labels, inner_folds_n, seed = inner_seed
        )
      }
      for (inner_idx in seq_along(inner_folds)) {
        validation_positions <- inner_folds[[inner_idx]]
        training_positions <- setdiff(
          seq_along(outer_indices), validation_positions
        )
        training_indices <- outer_indices[training_positions]
        validation_indices <- outer_indices[validation_positions]
        processed <- preprocess_inner_split(
          training_indices, validation_indices
        )
        train_labels <- outcome[training_indices]
        validation_labels <- outcome[validation_indices]

        for (method in methods) {
          output_path <- result_path(
            repeat_idx, fold_idx, inner_idx, method
          )
          if (file.exists(output_path)) next
          if (over_budget()) {
            cat("[budget] stop in nested panel-size evaluation\n")
            quit(save = "no")
          }
          fit_started <- proc.time()[["elapsed"]]
          chosen_alpha <- NA_real_
          alpha_auc <- stats::setNames(
            rep(NA_real_, length(alpha_grid)), alpha_grid
          )
          if (method == "GS_full_ungrouped") {
            alpha_results <- vector("list", length(alpha_grid))
            for (alpha_index in seq_along(alpha_grid)) {
              alpha_value <- alpha_grid[alpha_index]
              saved_alpha <- alpha_path(
                repeat_idx, fold_idx, inner_idx, alpha_value
              )
              if (file.exists(saved_alpha)) {
                alpha_results[[alpha_index]] <- readRDS(saved_alpha)
              } else {
                alpha_results[[alpha_index]] <- fit_gs_alpha(
                  processed$train, train_labels, alpha_value
                )
                saveRDS(alpha_results[[alpha_index]], saved_alpha)
              }
              alpha_auc[alpha_index] <-
                alpha_results[[alpha_index]]$internal_oob_auc
            }
            best_alpha_index <- which.max(alpha_auc)
            chosen_alpha <- alpha_grid[best_alpha_index]
            ranking <- alpha_results[[best_alpha_index]]$gene_scores$gene
            clear_run_cache()
          } else {
            seed <- method_seed(repeat_idx, fold_idx, inner_idx, method)
            ranking <- fit_competitor(
              method, processed$train, train_labels, seed
            )$ranked
          }
          fit_seconds <- proc.time()[["elapsed"]] - fit_started
          if (length(ranking) != length(processed$pool) ||
              anyDuplicated(ranking) || !setequal(ranking, processed$pool)) {
            stop(sprintf(
              "%s produced an invalid ranking in r%d f%d i%d",
              method, repeat_idx, fold_idx, inner_idx
            ), call. = FALSE)
          }
          curve <- do.call(rbind, lapply(panel_sizes, function(panel_size) {
            panel <- head(ranking, panel_size)
            seed <- evaluation_seed(
              repeat_idx, fold_idx, inner_idx, panel_size
            )
            scores <- predict_with_ensemble(
              processed$train[, panel, drop = FALSE], train_labels,
              processed$validation[, panel, drop = FALSE],
              random_seed = seed
            )
            data.frame(
              k = panel_size,
              AUC = bench_auc(validation_labels, scores),
              evaluation_seed = seed,
              evaluator_components = paste(
                attr(scores, "components_used"), collapse = "+"
              )
            )
          }))
          saveRDS(list(
            dataset = dataset,
            repeat_idx = repeat_idx,
            fold_idx = fold_idx,
            inner_idx = inner_idx,
            method = method,
            training_indices = training_indices,
            validation_indices = validation_indices,
            pool = processed$pool,
            ranking = ranking,
            ranking_seed = method_seed(
              repeat_idx, fold_idx, inner_idx, method
            ),
            chosen_alpha = chosen_alpha,
            alpha_internal_oob_auc = alpha_auc,
            fit_seconds = fit_seconds,
            curve = curve
          ), output_path)
          cat(sprintf(
            "[%s r%d f%d i%d] %s %.1f s\n",
            dataset, repeat_idx, fold_idx, inner_idx, method, fit_seconds
          ))
        }
      }
    }
  }
}

if (stage %in% c("assemble", "all")) {
  outer_evaluation <- read.csv(
    file.path(dataset_dir, "eval_deterministic.csv"),
    stringsAsFactors = FALSE
  )
  fold_rows <- list()
  curve_rows <- list()
  outer_rows <- list()
  for (repeat_idx in 1:3) {
    for (fold_idx in 1:5) {
      for (method in methods) {
        method_results <- lapply(seq_len(inner_folds_n), function(inner_idx) {
          path <- result_path(repeat_idx, fold_idx, inner_idx, method)
          if (!file.exists(path)) stop("Missing nested result: ", path)
          readRDS(path)
        })
        for (result in method_results) {
          fold_rows[[length(fold_rows) + 1L]] <- data.frame(
            repeat_idx = repeat_idx,
            fold_idx = fold_idx,
            inner_idx = result$inner_idx,
            method = method,
            k = result$curve$k,
            AUC = result$curve$AUC,
            evaluation_seed = result$curve$evaluation_seed,
            chosen_alpha = result$chosen_alpha,
            fit_seconds = result$fit_seconds
          )
        }
        fold_table <- do.call(rbind, lapply(method_results, function(result) {
          data.frame(
            inner_idx = result$inner_idx,
            k = result$curve$k,
            AUC = result$curve$AUC
          )
        }))
        curve <- aggregate(AUC ~ k, fold_table, function(value) {
          c(mean = mean(value), sd = stats::sd(value),
            se = stats::sd(value) / sqrt(length(value)))
        })
        curve <- data.frame(
          k = curve$k,
          inner_mean_AUC = curve$AUC[, "mean"],
          inner_sd_AUC = curve$AUC[, "sd"],
          inner_se_AUC = curve$AUC[, "se"]
        )
        best_index <- which.max(curve$inner_mean_AUC)
        threshold <- curve$inner_mean_AUC[best_index] -
          curve$inner_se_AUC[best_index]
        chosen_k <- min(curve$k[curve$inner_mean_AUC >= threshold])
        curve$repeat_idx <- repeat_idx
        curve$fold_idx <- fold_idx
        curve$method <- method
        curve$best_inner_k <- curve$k[best_index]
        curve$one_se_threshold <- threshold
        curve$chosen_k <- chosen_k
        curve_rows[[length(curve_rows) + 1L]] <- curve

        observed <- outer_evaluation[
          outer_evaluation$repeat_idx == repeat_idx &
            outer_evaluation$fold_idx == fold_idx &
            outer_evaluation$arm == method &
            outer_evaluation$k == chosen_k, , drop = FALSE
        ]
        random <- outer_evaluation[
          outer_evaluation$repeat_idx == repeat_idx &
            outer_evaluation$fold_idx == fold_idx &
            outer_evaluation$arm == "Random" &
            outer_evaluation$k == chosen_k, , drop = FALSE
        ]
        stopifnot(nrow(observed) == 1L, nrow(random) == 1L,
                  observed$evaluation_seed == random$evaluation_seed)
        outer_rows[[length(outer_rows) + 1L]] <- data.frame(
          repeat_idx = repeat_idx,
          fold_idx = fold_idx,
          method = method,
          chosen_k = chosen_k,
          best_inner_k = curve$k[best_index],
          best_inner_mean_AUC = curve$inner_mean_AUC[best_index],
          best_inner_se_AUC = curve$inner_se_AUC[best_index],
          one_se_threshold = threshold,
          outer_AUC = observed$AUC,
          matched_random_AUC = random$AUC,
          outer_delta_random = observed$AUC - random$AUC,
          outer_evaluation_seed = observed$evaluation_seed
        )
      }
    }
  }
  fold_table <- do.call(rbind, fold_rows)
  curve_table <- do.call(rbind, curve_rows)
  outer_table <- do.call(rbind, outer_rows)
  stopifnot(
    nrow(fold_table) == 15L * inner_folds_n * length(methods) *
      length(panel_sizes),
    nrow(curve_table) == 15L * length(methods) * length(panel_sizes),
    nrow(outer_table) == 15L * length(methods),
    all(is.finite(fold_table$AUC)),
    all(is.finite(outer_table$outer_AUC)),
    all(outer_table$chosen_k %in% panel_sizes)
  )
  write.csv(fold_table, fold_output, row.names = FALSE)
  write.csv(curve_table, curve_output, row.names = FALSE)
  write.csv(outer_table, outer_output, row.names = FALSE)
  cat(sprintf("%s older-seven nested results assembled\n", dataset))
}

cat(sprintf("%s adaptive stage=%s complete\n", dataset, stage))
