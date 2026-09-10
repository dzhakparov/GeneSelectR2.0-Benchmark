#!/usr/bin/env Rscript

# Summarise prediction across all panel sizes, leakage-safe adaptive panel-size
# selection, and post-hoc biological evidence. Dataset-level estimates are
# averaged with equal weight so that the largest cohort does not dominate.

datasets <- c("GSE16879", "GSE91061", "GSE92415", "GSE206285")
methods <- c("GS_full_ungrouped", "DGE", "LASSO", "ElasticNet", "mRMR",
             "Boruta", "RF_importance")
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
result_root <- file.path(
  "redesign", "results_frozen_external_exact_2026-08-31"
)
benchmark_root <- file.path(result_root, "validation_benchmark")
output_dir <- file.path(result_root, "efficiency_biology")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bind_dataset_csv <- function(file_name) {
  do.call(rbind, lapply(datasets, function(dataset) {
    path <- file.path(benchmark_root, dataset, file_name)
    if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
    value <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    value$dataset <- dataset
    value
  }))
}

standard_error <- function(value) {
  stats::sd(value) / sqrt(length(value))
}

corrected_cv_interval <- function(value, test_train_ratio = 1 / 4) {
  corrected_se <- sqrt(
    (1 / length(value) + test_train_ratio) * stats::var(value)
  )
  critical_value <- stats::qt(0.975, df = length(value) - 1L)
  c(
    corrected_cv_se = corrected_se,
    corrected_cv_low = mean(value) - critical_value * corrected_se,
    corrected_cv_high = mean(value) + critical_value * corrected_se
  )
}

group_apply <- function(data, group_columns, function_name) {
  groups <- split(
    data,
    interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  )
  rows <- lapply(groups, function(value) {
    group <- value[1L, group_columns, drop = FALSE]
    cbind(group, function_name(value))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

trapezoid_area <- function(x, y) {
  order_index <- order(x)
  x <- x[order_index]
  y <- y[order_index]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2) /
    (max(x) - min(x))
}

original <- bind_dataset_csv("eval_results.csv")
names(original)[names(original) == "arm"] <- "method"
competitor <- bind_dataset_csv("competitor_eval_results.csv")
prediction <- rbind(
  original[, c("dataset", "repeat_idx", "fold_idx", "method", "k", "AUC",
               "evaluation_seed")],
  competitor[, c("dataset", "repeat_idx", "fold_idx", "method", "k", "AUC",
                 "evaluation_seed")]
)
stopifnot(
  nrow(prediction) == length(datasets) * 15L * 8L * length(panel_sizes),
  !anyDuplicated(prediction[c("dataset", "repeat_idx", "fold_idx", "method",
                              "k")]),
  all(is.finite(prediction$AUC)),
  setequal(unique(prediction$method), c(methods, "Random"))
)

random <- prediction[prediction$method == "Random", ]
names(random)[names(random) == "AUC"] <- "random_AUC"
observed <- prediction[prediction$method != "Random", ]
prediction_matched <- merge(
  observed,
  random[, c("dataset", "repeat_idx", "fold_idx", "k", "random_AUC",
             "evaluation_seed")],
  by = c("dataset", "repeat_idx", "fold_idx", "k"),
  suffixes = c("", "_random"), all.x = TRUE, sort = FALSE
)
stopifnot(
  nrow(prediction_matched) == length(datasets) * 15L * length(methods) *
    length(panel_sizes),
  prediction_matched$evaluation_seed ==
    prediction_matched$evaluation_seed_random
)
prediction_matched$delta_random <-
  prediction_matched$AUC - prediction_matched$random_AUC

fixed_panel_summary <- group_apply(
  prediction_matched, c("dataset", "method", "k"), function(value) {
    interval <- corrected_cv_interval(value$delta_random)
    data.frame(
      mean_AUC = mean(value$AUC),
      se_AUC = standard_error(value$AUC),
      mean_random_AUC = mean(value$random_AUC),
      mean_delta_random = mean(value$delta_random),
      se_delta_random = standard_error(value$delta_random),
      corrected_cv_se_delta = interval[["corrected_cv_se"]],
      corrected_cv_low_delta = interval[["corrected_cv_low"]],
      corrected_cv_high_delta = interval[["corrected_cv_high"]]
    )
  }
)
write.csv(fixed_panel_summary,
          file.path(output_dir, "fixed_panel_prediction_summary.csv"),
          row.names = FALSE)

efficiency_split <- group_apply(
  prediction_matched,
  c("dataset", "repeat_idx", "fold_idx", "method"),
  function(value) data.frame(
    integrated_delta_random = trapezoid_area(log(value$k), value$delta_random),
    integrated_AUC = trapezoid_area(log(value$k), value$AUC)
  )
)
efficiency_dataset <- group_apply(
  efficiency_split, c("dataset", "method"), function(value) data.frame(
    mean_integrated_delta_random = mean(value$integrated_delta_random),
    se_integrated_delta_random = standard_error(
      value$integrated_delta_random
    ),
    corrected_cv_se_delta = corrected_cv_interval(
      value$integrated_delta_random
    )[["corrected_cv_se"]],
    corrected_cv_low_delta = corrected_cv_interval(
      value$integrated_delta_random
    )[["corrected_cv_low"]],
    corrected_cv_high_delta = corrected_cv_interval(
      value$integrated_delta_random
    )[["corrected_cv_high"]],
    mean_integrated_AUC = mean(value$integrated_AUC),
    se_integrated_AUC = standard_error(value$integrated_AUC)
  )
)
efficiency_overall <- group_apply(
  efficiency_dataset, "method", function(value) data.frame(
    dataset_balanced_delta_random = mean(
      value$mean_integrated_delta_random
    ),
    dataset_sd_delta_random = stats::sd(
      value$mean_integrated_delta_random
    ),
    dataset_balanced_AUC = mean(value$mean_integrated_AUC)
  )
)
efficiency_overall <- efficiency_overall[
  order(-efficiency_overall$dataset_balanced_delta_random), ]
efficiency_overall$rank <- seq_len(nrow(efficiency_overall))
write.csv(efficiency_split,
          file.path(output_dir, "panel_efficiency_by_split.csv"),
          row.names = FALSE)
write.csv(efficiency_dataset,
          file.path(output_dir, "panel_efficiency_by_dataset.csv"),
          row.names = FALSE)
write.csv(efficiency_overall,
          file.path(output_dir, "panel_efficiency_overall.csv"),
          row.names = FALSE)

adaptive <- bind_dataset_csv("adaptive_nested_results_v3.csv")
adaptive_folds <- bind_dataset_csv("adaptive_inner_fold_results_v3.csv")
adaptive_curves <- bind_dataset_csv("adaptive_inner_curves_v3.csv")
stopifnot(
  nrow(adaptive) == length(datasets) * 15L * length(methods),
  nrow(adaptive_folds) == length(datasets) * 15L * 3L * length(methods) *
    length(panel_sizes),
  nrow(adaptive_curves) == length(datasets) * 15L * length(methods) *
    length(panel_sizes),
  !anyDuplicated(adaptive[c("dataset", "repeat_idx", "fold_idx", "method")]),
  all(is.finite(adaptive$outer_AUC)),
  all(adaptive$chosen_k %in% panel_sizes)
)

reconstructed_choices <- group_apply(
  adaptive_folds,
  c("dataset", "repeat_idx", "fold_idx", "method"),
  function(value) {
    mean_auc <- aggregate(AUC ~ k, value, mean)
    se_auc <- aggregate(AUC ~ k, value, standard_error)
    curve <- merge(mean_auc, se_auc, by = "k", suffixes = c("_mean", "_se"))
    best_candidates <- which(
      abs(curve$AUC_mean - max(curve$AUC_mean)) <= 1e-12
    )
    best <- best_candidates[which.min(curve$k[best_candidates])]
    threshold <- curve$AUC_mean[best] - curve$AUC_se[best]
    data.frame(
      reconstructed_k = min(curve$k[curve$AUC_mean >= threshold]),
      reconstructed_best_k = curve$k[best]
    )
  }
)
choice_check <- merge(
  adaptive, reconstructed_choices,
  by = c("dataset", "repeat_idx", "fold_idx", "method"),
  all.x = TRUE, sort = FALSE
)
stopifnot(
  choice_check$chosen_k == choice_check$reconstructed_k,
  choice_check$best_inner_k == choice_check$reconstructed_best_k
)

adaptive_dataset <- group_apply(
  adaptive, c("dataset", "method"), function(value) {
    interval <- corrected_cv_interval(value$outer_delta_random)
    data.frame(
      mean_outer_AUC = mean(value$outer_AUC),
      se_outer_AUC = standard_error(value$outer_AUC),
      mean_matched_random_AUC = mean(value$matched_random_AUC),
      mean_delta_random = mean(value$outer_delta_random),
      se_delta_random = standard_error(value$outer_delta_random),
      corrected_cv_se_delta = interval[["corrected_cv_se"]],
      corrected_cv_low_delta = interval[["corrected_cv_low"]],
      corrected_cv_high_delta = interval[["corrected_cv_high"]],
      median_chosen_k = stats::median(value$chosen_k),
      mean_chosen_k = mean(value$chosen_k),
      fraction_k_at_most_50 = mean(value$chosen_k <= 50L)
    )
  }
)
adaptive_overall <- group_apply(
  adaptive_dataset, "method", function(value) data.frame(
    dataset_balanced_outer_AUC = mean(value$mean_outer_AUC),
    dataset_sd_outer_AUC = stats::sd(value$mean_outer_AUC),
    dataset_balanced_delta_random = mean(value$mean_delta_random),
    dataset_sd_delta_random = stats::sd(value$mean_delta_random),
    median_of_dataset_median_k = stats::median(value$median_chosen_k),
    mean_of_dataset_mean_k = mean(value$mean_chosen_k)
  )
)
adaptive_overall <- adaptive_overall[
  order(-adaptive_overall$dataset_balanced_delta_random), ]
adaptive_overall$rank <- seq_len(nrow(adaptive_overall))
adaptive_k_counts <- as.data.frame(table(
  adaptive$dataset, adaptive$method, adaptive$chosen_k
), stringsAsFactors = FALSE)
names(adaptive_k_counts) <- c("dataset", "method", "chosen_k", "n_splits")
adaptive_k_counts <- adaptive_k_counts[adaptive_k_counts$n_splits > 0L, ]
write.csv(adaptive_dataset,
          file.path(output_dir, "adaptive_prediction_by_dataset.csv"),
          row.names = FALSE)
write.csv(adaptive_overall,
          file.path(output_dir, "adaptive_prediction_overall.csv"),
          row.names = FALSE)
write.csv(adaptive_k_counts,
          file.path(output_dir, "adaptive_panel_size_counts.csv"),
          row.names = FALSE)
write.csv(choice_check,
          file.path(output_dir, "adaptive_choice_reconstruction.csv"),
          row.names = FALSE)

fixed_10 <- prediction_matched[prediction_matched$k == 10L, c(
  "dataset", "repeat_idx", "fold_idx", "method", "AUC", "random_AUC",
  "delta_random"
)]
names(fixed_10)[names(fixed_10) == "AUC"] <- "fixed_10_AUC"
names(fixed_10)[names(fixed_10) == "random_AUC"] <- "fixed_10_random_AUC"
names(fixed_10)[names(fixed_10) == "delta_random"] <-
  "fixed_10_delta_random"
adaptive_fixed_10 <- merge(
  adaptive, fixed_10,
  by = c("dataset", "repeat_idx", "fold_idx", "method"),
  all.x = TRUE, sort = FALSE
)
stopifnot(nrow(adaptive_fixed_10) == nrow(adaptive))
adaptive_fixed_10$adaptive_minus_fixed_10_AUC <-
  adaptive_fixed_10$outer_AUC - adaptive_fixed_10$fixed_10_AUC
adaptive_fixed_10$adaptive_minus_fixed_10_delta <-
  adaptive_fixed_10$outer_delta_random -
  adaptive_fixed_10$fixed_10_delta_random
adaptive_fixed_10_dataset <- group_apply(
  adaptive_fixed_10, c("dataset", "method"), function(value) data.frame(
    mean_adaptive_minus_fixed_10_AUC = mean(
      value$adaptive_minus_fixed_10_AUC
    ),
    mean_adaptive_minus_fixed_10_delta = mean(
      value$adaptive_minus_fixed_10_delta
    ),
    fraction_adaptive_AUC_above_fixed_10 = mean(
      value$adaptive_minus_fixed_10_AUC > 0
    ),
    fraction_chosen_k_above_10 = mean(value$chosen_k > 10L)
  )
)
adaptive_fixed_10_overall <- group_apply(
  adaptive_fixed_10_dataset, "method", function(value) data.frame(
    dataset_balanced_adaptive_minus_fixed_10_AUC = mean(
      value$mean_adaptive_minus_fixed_10_AUC
    ),
    dataset_balanced_adaptive_minus_fixed_10_delta = mean(
      value$mean_adaptive_minus_fixed_10_delta
    ),
    mean_fraction_adaptive_AUC_above_fixed_10 = mean(
      value$fraction_adaptive_AUC_above_fixed_10
    )
  )
)
write.csv(adaptive_fixed_10_dataset,
          file.path(output_dir, "adaptive_vs_fixed_10_by_dataset.csv"),
          row.names = FALSE)
write.csv(adaptive_fixed_10_overall,
          file.path(output_dir, "adaptive_vs_fixed_10_overall.csv"),
          row.names = FALSE)

biology <- bind_dataset_csv("biology_multiaxis_v2.csv")
stopifnot(
  nrow(biology) == length(datasets) * 15L * length(methods) *
    length(panel_sizes),
  !anyDuplicated(biology[c("dataset", "repeat_idx", "fold_idx", "method",
                           "k")]),
  all(biology$n_null == 1000L),
  all(biology$pool_size == 2000L)
)

biology_axes <- data.frame(
  axis = c("Open Targets", "GO semantic", "Hallmark", "STRING"),
  enrichment = c("open_targets_enrichment", "GO_semantic_enrichment",
                 "hallmark_enrichment", "string_enrichment"),
  empirical_p = c("open_targets_empirical_p", "GO_semantic_empirical_p",
                  "hallmark_empirical_p", "string_empirical_p"),
  coverage = c("open_targets_overlap", "GO_annotated",
               "hallmark_annotated", "string_mapped"),
  stringsAsFactors = FALSE
)

biology_to_long <- function(value) {
  rows <- lapply(seq_len(nrow(biology_axes)), function(axis_index) {
    axis <- biology_axes[axis_index, ]
    data.frame(
      dataset = value$dataset,
      repeat_idx = value$repeat_idx,
      fold_idx = value$fold_idx,
      method = value$method,
      k = value$k,
      axis = axis$axis,
      enrichment = value[[axis$enrichment]],
      empirical_p = value[[axis$empirical_p]],
      coverage = value[[axis$coverage]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

biology_long <- biology_to_long(biology)
stopifnot(all(is.finite(biology_long$enrichment)),
          all(is.finite(biology_long$empirical_p)))
biology_panel_summary <- group_apply(
  biology_long, c("dataset", "method", "k", "axis"), function(value) {
    data.frame(
      median_enrichment = stats::median(value$enrichment),
      mean_enrichment = mean(value$enrichment),
      fraction_enrichment_above_1 = mean(value$enrichment > 1),
      fraction_empirical_p_at_most_0_05 = mean(value$empirical_p <= 0.05),
      median_coverage = stats::median(value$coverage)
    )
  }
)
write.csv(biology_panel_summary,
          file.path(output_dir, "biology_by_dataset_method_panel.csv"),
          row.names = FALSE)

adaptive_biology <- merge(
  adaptive[, c("dataset", "repeat_idx", "fold_idx", "method", "chosen_k",
               "outer_AUC", "outer_delta_random")],
  biology,
  by.x = c("dataset", "repeat_idx", "fold_idx", "method", "chosen_k"),
  by.y = c("dataset", "repeat_idx", "fold_idx", "method", "k"),
  all.x = TRUE, sort = FALSE
)
stopifnot(nrow(adaptive_biology) == nrow(adaptive),
          !anyNA(adaptive_biology$n_null))
names(adaptive_biology)[names(adaptive_biology) == "chosen_k"] <- "k"
adaptive_biology_long <- biology_to_long(adaptive_biology)
adaptive_biology_long$outer_AUC <- rep(
  adaptive_biology$outer_AUC, times = nrow(biology_axes)
)
adaptive_biology_long$outer_delta_random <- rep(
  adaptive_biology$outer_delta_random, times = nrow(biology_axes)
)
adaptive_biology_dataset <- group_apply(
  adaptive_biology_long, c("dataset", "method", "axis"), function(value) {
    data.frame(
      median_enrichment = stats::median(value$enrichment),
      fraction_enrichment_above_1 = mean(value$enrichment > 1),
      fraction_empirical_p_at_most_0_05 = mean(value$empirical_p <= 0.05),
      median_coverage = stats::median(value$coverage)
    )
  }
)
adaptive_biology_overall <- group_apply(
  adaptive_biology_dataset, c("method", "axis"), function(value) {
    data.frame(
      mean_of_dataset_median_enrichment = mean(value$median_enrichment),
      median_of_dataset_median_enrichment = stats::median(
        value$median_enrichment
      ),
      datasets_with_median_enrichment_above_1 = sum(
        value$median_enrichment > 1
      ),
      mean_fraction_significant = mean(
        value$fraction_empirical_p_at_most_0_05
      ),
      mean_median_coverage = mean(value$median_coverage)
    )
  }
)
write.csv(adaptive_biology_dataset,
          file.path(output_dir, "adaptive_biology_by_dataset.csv"),
          row.names = FALSE)
write.csv(adaptive_biology_overall,
          file.path(output_dir, "adaptive_biology_overall.csv"),
          row.names = FALSE)

prediction_biology <- merge(
  prediction_matched,
  biology_long,
  by = c("dataset", "repeat_idx", "fold_idx", "method", "k"),
  all.x = TRUE, sort = FALSE
)
stopifnot(nrow(prediction_biology) == nrow(prediction_matched) *
            nrow(biology_axes),
          !anyNA(prediction_biology$enrichment))
biology_prediction_correlation <- group_apply(
  prediction_biology, c("dataset", "axis"), function(value) {
    data.frame(
      spearman_enrichment_AUC = stats::cor(
        value$enrichment, value$AUC, method = "spearman"
      ),
      spearman_enrichment_delta_random = stats::cor(
        value$enrichment, value$delta_random, method = "spearman"
      ),
      n = nrow(value)
    )
  }
)
write.csv(biology_prediction_correlation,
          file.path(output_dir, "biology_prediction_correlations.csv"),
          row.names = FALSE)

qa <- data.frame(
  check = c(
    "prediction_rows", "adaptive_outer_rows", "adaptive_inner_rows",
    "adaptive_curve_rows", "adaptive_choices_reconstructed",
    "biology_rows", "biology_null_draws", "biology_candidate_pool"
  ),
  observed = c(
    nrow(prediction), nrow(adaptive), nrow(adaptive_folds),
    nrow(adaptive_curves), sum(
      choice_check$chosen_k == choice_check$reconstructed_k
    ), nrow(biology), unique(biology$n_null), unique(biology$pool_size)
  ),
  expected = c(
    2880L, 420L, 7560L, 2520L, 420L, 2520L, 1000L, 2000L
  ),
  stringsAsFactors = FALSE
)
qa$passed <- qa$observed == qa$expected
stopifnot(all(qa$passed))
write.csv(qa, file.path(output_dir, "validation_checks.csv"),
          row.names = FALSE)

cat("Efficiency, adaptive panel-size, and biology summaries complete.\n")
