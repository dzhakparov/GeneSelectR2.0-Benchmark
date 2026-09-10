#!/usr/bin/env Rscript

# Summaries for the seven older datasets using the same prediction-efficiency,
# nested panel-size, and biological-evidence definitions as the four-dataset
# frozen external analysis. Datasets receive equal weight in overall estimates.

datasets <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683",
  "imvigor210", "sosall"
)
validation_datasets_old <- datasets[1:5]
methods <- c(
  "GS_full_ungrouped", "DGE", "LASSO", "ElasticNet", "mRMR", "Boruta",
  "RF_importance"
)
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
results_root <- file.path("redesign", "results_corrected")
output_dir <- file.path(
  results_root, "older7_efficiency_biology_2026-09-01"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_dir <- function(dataset) {
  if (dataset %in% validation_datasets_old) {
    file.path(results_root, "validation_benchmark", dataset)
  } else {
    file.path(results_root, "full_recipe", dataset)
  }
}
bind_dataset_csv <- function(file_name) {
  do.call(rbind, lapply(datasets, function(dataset) {
    path <- file.path(dataset_dir(dataset), file_name)
    if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
    value <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    value$dataset <- dataset
    value
  }))
}
standard_error <- function(value) stats::sd(value) / sqrt(length(value))
corrected_cv_interval <- function(value, test_train_ratio = 1 / 4) {
  corrected_se <- sqrt(
    (1 / length(value) + test_train_ratio) * stats::var(value)
  )
  critical <- stats::qt(0.975, df = length(value) - 1L)
  c(se = corrected_se,
    low = mean(value) - critical * corrected_se,
    high = mean(value) + critical * corrected_se)
}
group_apply <- function(data, group_columns, function_name) {
  groups <- split(
    data, interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  )
  rows <- lapply(groups, function(value) {
    cbind(value[1L, group_columns, drop = FALSE], function_name(value))
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}
trapezoid_area <- function(x, y) {
  index <- order(x)
  x <- x[index]
  y <- y[index]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2) /
    (max(x) - min(x))
}

# Prediction across the complete panel-size range.
prediction <- bind_dataset_csv("eval_deterministic.csv")
names(prediction)[names(prediction) == "arm"] <- "method"
prediction <- prediction[
  prediction$method %in% c(methods, "Random") &
    prediction$k %in% panel_sizes,
  c("dataset", "repeat_idx", "fold_idx", "method", "k", "AUC",
    "evaluation_seed")
]
stopifnot(
  nrow(prediction) == length(datasets) * 15L * 8L * length(panel_sizes),
  !anyDuplicated(prediction[c(
    "dataset", "repeat_idx", "fold_idx", "method", "k"
  )]),
  all(is.finite(prediction$AUC))
)
random <- prediction[prediction$method == "Random", ]
names(random)[names(random) == "AUC"] <- "random_AUC"
observed <- prediction[prediction$method != "Random", ]
matched <- merge(
  observed,
  random[, c(
    "dataset", "repeat_idx", "fold_idx", "k", "random_AUC",
    "evaluation_seed"
  )],
  by = c("dataset", "repeat_idx", "fold_idx", "k"),
  suffixes = c("", "_random"), all.x = TRUE, sort = FALSE
)
stopifnot(matched$evaluation_seed == matched$evaluation_seed_random)
matched$delta_random <- matched$AUC - matched$random_AUC

fixed_summary <- group_apply(
  matched, c("dataset", "method", "k"), function(value) {
    interval <- corrected_cv_interval(value$delta_random)
    data.frame(
      mean_AUC = mean(value$AUC),
      mean_random_AUC = mean(value$random_AUC),
      mean_delta_random = mean(value$delta_random),
      se_delta_random = standard_error(value$delta_random),
      corrected_cv_se_delta = interval[["se"]],
      corrected_cv_low_delta = interval[["low"]],
      corrected_cv_high_delta = interval[["high"]],
      positive_split_fraction = mean(value$delta_random > 0)
    )
  }
)
write.csv(fixed_summary, file.path(
  output_dir, "fixed_panel_prediction_summary.csv"
), row.names = FALSE)

efficiency_split <- group_apply(
  matched, c("dataset", "repeat_idx", "fold_idx", "method"),
  function(value) data.frame(
    integrated_delta_random = trapezoid_area(
      log(value$k), value$delta_random
    ),
    integrated_AUC = trapezoid_area(log(value$k), value$AUC)
  )
)
efficiency_dataset <- group_apply(
  efficiency_split, c("dataset", "method"), function(value) {
    interval <- corrected_cv_interval(value$integrated_delta_random)
    data.frame(
      mean_integrated_delta_random = mean(value$integrated_delta_random),
      corrected_cv_low_delta = interval[["low"]],
      corrected_cv_high_delta = interval[["high"]],
      mean_integrated_AUC = mean(value$integrated_AUC)
    )
  }
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
  order(-efficiency_overall$dataset_balanced_delta_random),
]
efficiency_overall$rank <- seq_len(nrow(efficiency_overall))
write.csv(efficiency_split, file.path(
  output_dir, "panel_efficiency_by_split.csv"
), row.names = FALSE)
write.csv(efficiency_dataset, file.path(
  output_dir, "panel_efficiency_by_dataset.csv"
), row.names = FALSE)
write.csv(efficiency_overall, file.path(
  output_dir, "panel_efficiency_overall.csv"
), row.names = FALSE)

# Leakage-safe adaptive panel size.
adaptive <- bind_dataset_csv("older7_adaptive_nested_results.csv")
adaptive_folds <- bind_dataset_csv(
  "older7_adaptive_inner_fold_results.csv"
)
adaptive_curves <- bind_dataset_csv("older7_adaptive_inner_curves.csv")
stopifnot(
  nrow(adaptive) == length(datasets) * 15L * length(methods),
  nrow(adaptive_folds) == length(datasets) * 15L * 3L *
    length(methods) * length(panel_sizes),
  nrow(adaptive_curves) == length(datasets) * 15L * length(methods) *
    length(panel_sizes),
  all(adaptive$chosen_k %in% panel_sizes),
  all(is.finite(adaptive$outer_AUC))
)
adaptive_dataset <- group_apply(
  adaptive, c("dataset", "method"), function(value) {
    interval <- corrected_cv_interval(value$outer_delta_random)
    data.frame(
      mean_outer_AUC = mean(value$outer_AUC),
      mean_matched_random_AUC = mean(value$matched_random_AUC),
      mean_delta_random = mean(value$outer_delta_random),
      corrected_cv_low_delta = interval[["low"]],
      corrected_cv_high_delta = interval[["high"]],
      median_chosen_k = stats::median(value$chosen_k),
      mean_chosen_k = mean(value$chosen_k),
      fraction_k_at_most_50 = mean(value$chosen_k <= 50L)
    )
  }
)
adaptive_overall <- group_apply(
  adaptive_dataset, "method", function(value) data.frame(
    dataset_balanced_outer_AUC = mean(value$mean_outer_AUC),
    dataset_balanced_delta_random = mean(value$mean_delta_random),
    dataset_sd_delta_random = stats::sd(value$mean_delta_random),
    median_of_dataset_median_k = stats::median(value$median_chosen_k),
    mean_of_dataset_mean_k = mean(value$mean_chosen_k)
  )
)
adaptive_overall <- adaptive_overall[
  order(-adaptive_overall$dataset_balanced_delta_random),
]
adaptive_overall$rank <- seq_len(nrow(adaptive_overall))
adaptive_k_counts <- as.data.frame(table(
  adaptive$dataset, adaptive$method, adaptive$chosen_k
), stringsAsFactors = FALSE)
names(adaptive_k_counts) <- c("dataset", "method", "chosen_k", "n_splits")
adaptive_k_counts <- adaptive_k_counts[adaptive_k_counts$n_splits > 0L, ]
write.csv(adaptive_dataset, file.path(
  output_dir, "adaptive_prediction_by_dataset.csv"
), row.names = FALSE)
write.csv(adaptive_overall, file.path(
  output_dir, "adaptive_prediction_overall.csv"
), row.names = FALSE)
write.csv(adaptive_k_counts, file.path(
  output_dir, "adaptive_panel_size_counts.csv"
), row.names = FALSE)

fixed_10 <- matched[matched$k == 10L, c(
  "dataset", "repeat_idx", "fold_idx", "method", "AUC", "delta_random"
)]
names(fixed_10)[names(fixed_10) == "AUC"] <- "fixed_10_AUC"
names(fixed_10)[names(fixed_10) == "delta_random"] <-
  "fixed_10_delta_random"
adaptive_fixed_10 <- merge(
  adaptive, fixed_10,
  by = c("dataset", "repeat_idx", "fold_idx", "method"),
  all.x = TRUE, sort = FALSE
)
adaptive_fixed_10$adaptive_minus_fixed_10_AUC <-
  adaptive_fixed_10$outer_AUC - adaptive_fixed_10$fixed_10_AUC
adaptive_fixed_10$adaptive_minus_fixed_10_delta <-
  adaptive_fixed_10$outer_delta_random -
  adaptive_fixed_10$fixed_10_delta_random
adaptive_fixed_10_summary <- group_apply(
  adaptive_fixed_10, c("dataset", "method"), function(value) data.frame(
    mean_adaptive_minus_fixed_10_AUC = mean(
      value$adaptive_minus_fixed_10_AUC
    ),
    mean_adaptive_minus_fixed_10_delta = mean(
      value$adaptive_minus_fixed_10_delta
    ),
    fraction_adaptive_AUC_above_fixed_10 = mean(
      value$adaptive_minus_fixed_10_AUC > 0
    )
  )
)
write.csv(adaptive_fixed_10_summary, file.path(
  output_dir, "adaptive_vs_fixed_10_by_dataset.csv"
), row.names = FALSE)

# Four external biology axes.
biology <- bind_dataset_csv("older7_biology_multiaxis.csv")
stopifnot(
  nrow(biology) == length(datasets) * 15L * length(methods) *
    length(panel_sizes),
  !anyDuplicated(biology[c(
    "dataset", "repeat_idx", "fold_idx", "method", "k"
  )]),
  all(biology$n_null == 1000L),
  all(biology$pool_size == 2000L)
)
axis_columns <- data.frame(
  axis = c("Open Targets", "GO semantic", "Hallmark", "STRING"),
  enrichment = c(
    "open_targets_enrichment", "GO_semantic_enrichment",
    "hallmark_enrichment", "string_enrichment"
  ),
  empirical_p = c(
    "open_targets_empirical_p", "GO_semantic_empirical_p",
    "hallmark_empirical_p", "string_empirical_p"
  ),
  coverage = c(
    "open_targets_overlap", "GO_annotated", "hallmark_annotated",
    "string_mapped"
  ),
  stringsAsFactors = FALSE
)
biology_long <- do.call(rbind, lapply(seq_len(nrow(axis_columns)), function(i) {
  data.frame(
    dataset = biology$dataset,
    repeat_idx = biology$repeat_idx,
    fold_idx = biology$fold_idx,
    method = biology$method,
    k = biology$k,
    axis = axis_columns$axis[i],
    enrichment = biology[[axis_columns$enrichment[i]]],
    empirical_p = biology[[axis_columns$empirical_p[i]]],
    coverage = biology[[axis_columns$coverage[i]]]
  )
}))
biology_summary <- group_apply(
  biology_long, c("dataset", "method", "k", "axis"), function(value) {
    valid <- is.finite(value$enrichment)
    data.frame(
      median_enrichment = if (any(valid)) {
        stats::median(value$enrichment[valid])
      } else {
        NA_real_
      },
      fraction_enrichment_above_1 = if (any(valid)) {
        mean(value$enrichment[valid] > 1)
      } else {
        NA_real_
      },
      fraction_empirical_p_at_most_0_05 = mean(
        value$empirical_p <= 0.05, na.rm = TRUE
      ),
      median_coverage = stats::median(value$coverage),
      valid_splits = sum(valid)
    )
  }
)
write.csv(biology_summary, file.path(
  output_dir, "biology_by_dataset_method_panel.csv"
), row.names = FALSE)

prediction_biology <- merge(
  matched, biology_long,
  by = c("dataset", "repeat_idx", "fold_idx", "method", "k"),
  all.x = TRUE, sort = FALSE
)
biology_correlations <- group_apply(
  prediction_biology, c("dataset", "axis"), function(value) {
    valid <- is.finite(value$enrichment) & is.finite(value$delta_random)
    data.frame(
      spearman_enrichment_delta_random = if (sum(valid) >= 3L) {
        stats::cor(
          value$enrichment[valid], value$delta_random[valid],
          method = "spearman"
        )
      } else {
        NA_real_
      },
      n = sum(valid)
    )
  }
)
write.csv(biology_correlations, file.path(
  output_dir, "biology_prediction_correlations.csv"
), row.names = FALSE)

qa <- data.frame(
  check = c(
    "prediction_rows", "adaptive_outer_rows", "adaptive_inner_rows",
    "adaptive_curve_rows", "biology_rows", "biology_null_draws",
    "biology_candidate_pool"
  ),
  observed = c(
    nrow(prediction), nrow(adaptive), nrow(adaptive_folds),
    nrow(adaptive_curves), nrow(biology), unique(biology$n_null),
    unique(biology$pool_size)
  ),
  expected = c(5040L, 735L, 13230L, 4410L, 4410L, 1000L, 2000L)
)
qa$passed <- qa$observed == qa$expected
stopifnot(all(qa$passed))
write.csv(qa, file.path(output_dir, "validation_checks.csv"),
          row.names = FALSE)
cat("Older-seven efficiency, adaptive, and biology summaries complete.\n")
