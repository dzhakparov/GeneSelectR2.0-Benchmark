#!/usr/bin/env Rscript

# Summarise the classical-comparator extension to the exact frozen external
# benchmark. Dataset means receive equal weight so a large cohort cannot
# dominate the four-dataset result.

result_root <- file.path(
  "redesign", "results_frozen_external_exact_2026-08-31"
)
benchmark_root <- file.path(result_root, "validation_benchmark")
summary_dir <- file.path(result_root, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE16879", "GSE91061", "GSE92415", "GSE206285")
competitors <- c("DGE", "LASSO", "ElasticNet", "mRMR", "Boruta",
                 "RF_importance")
methods <- c("GS_full_ungrouped", competitors, "Random")
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
primary_sizes <- c(10L, 20L, 50L)

read_dataset_results <- function(dataset) {
  dataset_dir <- file.path(benchmark_root, dataset)
  original <- read.csv(file.path(dataset_dir, "eval_results.csv"),
                       stringsAsFactors = FALSE)
  names(original)[names(original) == "arm"] <- "method"
  comparator <- read.csv(
    file.path(dataset_dir, "competitor_eval_results.csv"),
    stringsAsFactors = FALSE
  )
  output <- rbind(
    original[, c("repeat_idx", "fold_idx", "method", "k", "AUC",
                 "n_panel", "evaluation_seed", "evaluator_components")],
    comparator[, c("repeat_idx", "fold_idx", "method", "k", "AUC",
                   "n_panel", "evaluation_seed", "evaluator_components")]
  )
  output$dataset <- dataset
  output
}

results <- do.call(rbind, lapply(datasets, read_dataset_results))
results$method <- factor(results$method, levels = methods)
results <- results[order(results$dataset, results$repeat_idx,
                         results$fold_idx, results$method, results$k), ]

expected_rows <- length(datasets) * 15L * length(methods) *
  length(panel_sizes)
result_keys <- with(results, paste(dataset, repeat_idx, fold_idx, method, k))
stopifnot(
  nrow(results) == expected_rows,
  !anyDuplicated(result_keys),
  setequal(as.character(unique(results$method)), methods),
  setequal(unique(results$k), panel_sizes),
  all(is.finite(results$AUC)),
  all(results$AUC >= 0 & results$AUC <= 1),
  all(results$n_panel == results$k),
  all(results$evaluation_seed ==
        420000L + results$repeat_idx * 1000L +
        results$fold_idx * 100L + match(results$k, panel_sizes)),
  all(results$evaluator_components == "glmnet+xgboost+ranger")
)

primary <- results[results$k %in% primary_sizes, ]
primary_by_dataset <- aggregate(
  AUC ~ dataset + method, primary, mean
)
random_by_dataset <- primary_by_dataset[
  primary_by_dataset$method == "Random", c("dataset", "AUC")
]
names(random_by_dataset)[2L] <- "random_AUC"
primary_by_dataset <- merge(primary_by_dataset, random_by_dataset,
                            by = "dataset")
primary_by_dataset$delta_random <- primary_by_dataset$AUC -
  primary_by_dataset$random_AUC
primary_by_dataset$rank <- ave(
  -primary_by_dataset$AUC, primary_by_dataset$dataset,
  FUN = function(value) rank(value, ties.method = "min")
)
primary_by_dataset <- primary_by_dataset[
  order(primary_by_dataset$dataset, primary_by_dataset$rank),
]

overall_primary <- aggregate(
  cbind(AUC, delta_random, rank) ~ method, primary_by_dataset, mean
)
overall_primary$datasets_ranked_first <- vapply(
  overall_primary$method,
  function(method) sum(primary_by_dataset$method == method &
                         primary_by_dataset$rank == 1L),
  integer(1)
)
overall_primary <- overall_primary[
  order(overall_primary$delta_random, decreasing = TRUE),
]

by_dataset_k <- aggregate(AUC ~ dataset + method + k, results, mean)
random_by_dataset_k <- by_dataset_k[
  by_dataset_k$method == "Random", c("dataset", "k", "AUC")
]
names(random_by_dataset_k)[3L] <- "random_AUC"
by_dataset_k <- merge(by_dataset_k, random_by_dataset_k,
                      by = c("dataset", "k"))
by_dataset_k$delta_random <- by_dataset_k$AUC - by_dataset_k$random_AUC
overall_by_k <- aggregate(
  cbind(AUC, delta_random) ~ method + k, by_dataset_k, mean
)
overall_by_k$rank <- ave(
  -overall_by_k$AUC, overall_by_k$k,
  FUN = function(value) rank(value, ties.method = "min")
)
overall_by_k <- overall_by_k[order(overall_by_k$k, overall_by_k$rank), ]

# Average the three primary panel sizes inside each outer split before paired
# uncertainty calculations. The corrected standard error uses the 1/5 test
# fraction of five-fold CV: n_test / n_train = 1/4.
split_primary <- aggregate(
  AUC ~ dataset + repeat_idx + fold_idx + method, primary, mean
)
gs_split <- split_primary[
  split_primary$method == "GS_full_ungrouped",
  c("dataset", "repeat_idx", "fold_idx", "AUC")
]
names(gs_split)[4L] <- "GS_AUC"
paired <- merge(split_primary, gs_split,
                by = c("dataset", "repeat_idx", "fold_idx"))
paired$GS_minus_method <- paired$GS_AUC - paired$AUC
paired <- paired[paired$method != "GS_full_ungrouped", ]

paired_summary <- do.call(rbind, lapply(
  split(paired, list(paired$dataset, paired$method), drop = TRUE),
  function(group) {
    difference <- group$GS_minus_method
    n <- length(difference)
    corrected_se <- sqrt((1 / n + 1 / 4) * stats::var(difference))
    critical_value <- stats::qt(0.975, df = n - 1L)
    data.frame(
      dataset = group$dataset[1L],
      comparator = as.character(group$method[1L]),
      GS_minus_comparator = mean(difference),
      corrected_95_low = mean(difference) - critical_value * corrected_se,
      corrected_95_high = mean(difference) + critical_value * corrected_se,
      GS_split_wins = sum(difference > 0),
      split_ties = sum(difference == 0),
      n_splits = n
    )
  }
))
paired_summary <- paired_summary[
  order(paired_summary$dataset, paired_summary$comparator),
]

meta_files <- list.files(
  benchmark_root,
  pattern = "competitor_r[0-9]+_f[0-9]+_.*_meta[.]csv$",
  recursive = TRUE, full.names = TRUE
)
meta <- do.call(rbind, lapply(meta_files, read.csv,
                             stringsAsFactors = FALSE))

valid_rankings <- 0L
for (dataset in datasets) {
  dataset_dir <- file.path(benchmark_root, dataset)
  for (repeat_idx in 1:3) {
    for (fold_idx in 1:5) {
      split_data <- readRDS(file.path(
        dataset_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
      ))
      pool <- split_data$pools$var2000
      for (method in competitors) {
        ranking <- read.csv(file.path(
          dataset_dir,
          sprintf("competitor_r%d_f%d_%s.csv", repeat_idx, fold_idx,
                  method)
        ), stringsAsFactors = FALSE)$gene
        if (length(ranking) == 2000L && !anyDuplicated(ranking) &&
            setequal(ranking, pool)) {
          valid_rankings <- valid_rankings + 1L
        }
      }
    }
  }
}

selected_size_summary <- do.call(rbind, lapply(
  split(meta$selected_size, meta$method),
  function(value) data.frame(
    min_selected = min(value),
    median_selected = stats::median(value),
    max_selected = max(value),
    splits_below_10 = sum(value < 10L),
    splits_below_20 = sum(value < 20L),
    splits_below_50 = sum(value < 50L),
    n_splits = length(value)
  )
))
selected_size_summary$method <- rownames(selected_size_summary)
rownames(selected_size_summary) <- NULL
selected_size_summary <- selected_size_summary[
  , c("method", setdiff(names(selected_size_summary), "method"))
]

qa_checks <- data.frame(
  check = c(
    "total_evaluation_rows",
    "unique_evaluation_keys",
    "finite_bounded_AUC",
    "expected_panel_sizes",
    "expected_evaluation_seeds",
    "expected_evaluator_components",
    "competitor_ranking_files",
    "valid_competitor_rankings",
    "competitor_meta_files"
  ),
  observed = c(
    nrow(results),
    length(unique(result_keys)),
    sum(is.finite(results$AUC) & results$AUC >= 0 & results$AUC <= 1),
    sum(results$n_panel == results$k),
    sum(results$evaluation_seed ==
          420000L + results$repeat_idx * 1000L +
          results$fold_idx * 100L + match(results$k, panel_sizes)),
    sum(results$evaluator_components == "glmnet+xgboost+ranger"),
    length(list.files(
      benchmark_root,
      pattern = "competitor_r[0-9]+_f[0-9]+_(DGE|LASSO|ElasticNet|mRMR|Boruta|RF_importance)[.]csv$",
      recursive = TRUE
    )),
    valid_rankings,
    length(meta_files)
  ),
  expected = c(expected_rows, expected_rows, expected_rows, expected_rows,
               expected_rows, expected_rows, 360L, 360L, 360L)
)
qa_checks$passed <- qa_checks$observed == qa_checks$expected
stopifnot(all(qa_checks$passed))

write.csv(primary_by_dataset,
          file.path(summary_dir, "competitor_primary_by_dataset.csv"),
          row.names = FALSE)
write.csv(overall_primary,
          file.path(summary_dir, "competitor_overall_primary.csv"),
          row.names = FALSE)
write.csv(overall_by_k,
          file.path(summary_dir, "competitor_overall_by_k.csv"),
          row.names = FALSE)
write.csv(paired_summary,
          file.path(summary_dir, "competitor_paired_vs_gs.csv"),
          row.names = FALSE)
write.csv(selected_size_summary,
          file.path(summary_dir, "competitor_selected_size_summary.csv"),
          row.names = FALSE)
write.csv(qa_checks,
          file.path(summary_dir, "competitor_qa_checks.csv"),
          row.names = FALSE)

cat(sprintf(
  "Wrote comparator summaries: %d evaluation rows, %d primary rows.\n",
  nrow(results), nrow(primary)
))
