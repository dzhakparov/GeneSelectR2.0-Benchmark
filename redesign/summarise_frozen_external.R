#!/usr/bin/env Rscript

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
validation_root <- file.path(result_root, "validation_benchmark")
summary_root <- file.path(result_root, "summary")
dir.create(summary_root, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE16879", "GSE91061", "GSE92415", "GSE206285")
primary_panel_sizes <- c(10L, 20L, 50L)
all_panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)

summary_by_k <- list()
primary_by_dataset <- list()
fold_deltas <- list()
alpha_rows <- list()
qa_rows <- list()
reconstruction_rows <- list()

append_row <- function(container, row) {
  container[[length(container) + 1L]] <- row
  container
}

hashes_match <- function(recorded) {
  paths <- names(recorded)
  all(file.exists(paths)) && identical(
    unname(as.character(tools::md5sum(paths))),
    unname(as.character(recorded))
  )
}

for (dataset in datasets) {
  dataset_dir <- file.path(validation_root, dataset)
  eval_path <- file.path(dataset_dir, "eval_results.csv")
  evaluation <- read.csv(eval_path, stringsAsFactors = FALSE)

  keys <- paste(evaluation$repeat_idx, evaluation$fold_idx,
                evaluation$arm, evaluation$k)
  expected_seed <- 420000L + evaluation$repeat_idx * 1000L +
    evaluation$fold_idx * 100L + match(evaluation$k, all_panel_sizes)

  stopifnot(
    nrow(evaluation) == 180L,
    !anyDuplicated(keys),
    setequal(unique(evaluation$arm), c("GS_full_ungrouped", "Random")),
    all(is.finite(evaluation$AUC)),
    all(evaluation$AUC >= 0 & evaluation$AUC <= 1),
    all(evaluation$evaluation_seed == expected_seed),
    all(evaluation$n_panel == evaluation$k),
    all(evaluation$evaluator_components == "glmnet+xgboost+ranger")
  )

  gs <- evaluation[evaluation$arm == "GS_full_ungrouped",
                   c("repeat_idx", "fold_idx", "k", "AUC",
                     "evaluation_seed")]
  random <- evaluation[evaluation$arm == "Random",
                       c("repeat_idx", "fold_idx", "k", "AUC",
                         "evaluation_seed")]
  names(gs)[4:5] <- c("GS_AUC", "GS_seed")
  names(random)[4:5] <- c("Random_AUC", "Random_seed")
  paired <- merge(gs, random, by = c("repeat_idx", "fold_idx", "k"),
                  sort = TRUE)
  stopifnot(nrow(paired) == 90L, all(paired$GS_seed == paired$Random_seed))
  paired$delta <- paired$GS_AUC - paired$Random_AUC
  paired$dataset <- dataset
  fold_deltas[[dataset]] <- paired

  for (panel_size in all_panel_sizes) {
    block <- paired[paired$k == panel_size, ]
    summary_by_k <- append_row(summary_by_k, data.frame(
      dataset = dataset,
      k = panel_size,
      n_outer_splits = nrow(block),
      mean_GS_AUC = mean(block$GS_AUC),
      mean_Random_AUC = mean(block$Random_AUC),
      mean_delta = mean(block$delta),
      median_delta = median(block$delta),
      positive_cells = sum(block$delta > 0),
      tied_cells = sum(block$delta == 0)
    ))
  }

  primary <- paired[paired$k %in% primary_panel_sizes, ]
  split_primary <- aggregate(
    delta ~ repeat_idx + fold_idx, data = primary, FUN = mean
  )
  repeat_primary <- aggregate(delta ~ repeat_idx, data = primary, FUN = mean)

  # Nadeau-Bengio corrected resampled standard error for repeated 5-fold CV.
  # This is a sensitivity interval. The primary estimate remains descriptive.
  corrected_se <- sqrt((1 / nrow(split_primary) + 1 / 4) *
                         stats::var(split_primary$delta))
  critical <- stats::qt(0.975, df = nrow(split_primary) - 1L)

  primary_by_dataset <- append_row(primary_by_dataset, data.frame(
    dataset = dataset,
    n_split_k_cells = nrow(primary),
    mean_GS_AUC = mean(primary$GS_AUC),
    mean_Random_AUC = mean(primary$Random_AUC),
    primary_delta = mean(primary$delta),
    median_cell_delta = median(primary$delta),
    positive_cells = sum(primary$delta > 0),
    tied_cells = sum(primary$delta == 0),
    positive_split_means = sum(split_primary$delta > 0),
    tied_split_means = sum(split_primary$delta == 0),
    repeat_mean_min = min(repeat_primary$delta),
    repeat_mean_max = max(repeat_primary$delta),
    corrected_CI_low = mean(split_primary$delta) - critical * corrected_se,
    corrected_CI_high = mean(split_primary$delta) + critical * corrected_se
  ))

  ranking_paths <- sort(list.files(
    dataset_dir,
    pattern = "^ranking_r[1-3]_f[1-5]_GS_full_ungrouped[.]csv$",
    full.names = TRUE
  ))
  meta_paths <- sub("[.]csv$", "_meta.csv", ranking_paths)
  stopifnot(length(ranking_paths) == 15L, all(file.exists(meta_paths)))

  ranking_valid <- TRUE
  alpha_valid <- TRUE
  for (index in seq_along(ranking_paths)) {
    ranking <- read.csv(ranking_paths[index], stringsAsFactors = FALSE)
    ranking_valid <- ranking_valid && nrow(ranking) == 2000L &&
      !anyDuplicated(ranking$gene) && all(nzchar(ranking$gene)) &&
      all(is.finite(ranking$final_score)) &&
      all(diff(ranking$final_score) <= 1e-12)

    metadata <- read.csv(meta_paths[index], stringsAsFactors = FALSE)
    values <- setNames(metadata$value, metadata$key)
    selected_alpha <- as.numeric(values[["alpha"]])
    alpha_05_auc <- as.numeric(values[["auc_alpha0.5"]])
    alpha_10_auc <- as.numeric(values[["auc_alpha1.0"]])
    expected_alpha <- if (alpha_05_auc >= alpha_10_auc) 0.5 else 1.0
    alpha_valid <- alpha_valid && identical(selected_alpha, expected_alpha)
    parsed <- regmatches(
      basename(ranking_paths[index]),
      regexec("ranking_r([1-3])_f([1-5])_", basename(ranking_paths[index]))
    )[[1L]]
    alpha_rows <- append_row(alpha_rows, data.frame(
      dataset = dataset,
      repeat_idx = as.integer(parsed[2L]),
      fold_idx = as.integer(parsed[3L]),
      selected_alpha = selected_alpha,
      OOB_AUC_alpha_0.5 = alpha_05_auc,
      OOB_AUC_alpha_1.0 = alpha_10_auc
    ))
  }
  stopifnot(ranking_valid, alpha_valid)

  manifest <- readRDS(file.path(dataset_dir, "run_manifest.rds"))
  config <- manifest$configuration
  config_valid <- identical(config$dataset, dataset) &&
    identical(config$panel_sizes, as.numeric(all_panel_sizes)) &&
    identical(config$alpha_grid, c(0.5, 1.0)) &&
    identical(config$B, 50L) &&
    identical(config$calibration_permutations, 20L) &&
    identical(config$calibration_null_B, 20L) &&
    identical(config$n_workers, 7L) &&
    identical(config$fit_n_cores, 1L) &&
    identical(config$gs_arms, "GS_full_ungrouped") &&
    length(config$comp_arms) == 0L
  source_hash_valid <- hashes_match(manifest$source_hashes)
  input_hash_valid <- hashes_match(manifest$input_hashes)
  stopifnot(config_valid, source_hash_valid, input_hash_valid)

  # Independently reconstruct one GS result and its matched Random baseline.
  repeat_idx <- 1L
  fold_idx <- 1L
  panel_size <- 20L
  split <- readRDS(file.path(dataset_dir, "split_r1_f1.rds"))
  base <- readRDS(file.path(dataset_dir, "base_data.rds"))
  pool <- split$pools$var2000
  standardised <- standardise_split(
    split$train_raw[, pool, drop = FALSE],
    split$test_raw[, pool, drop = FALSE]
  )
  ranking <- read.csv(
    file.path(dataset_dir, "ranking_r1_f1_GS_full_ungrouped.csv"),
    stringsAsFactors = FALSE
  )$gene
  panel <- head(ranking, panel_size)
  evaluation_seed <- 420000L + repeat_idx * 1000L + fold_idx * 100L +
    match(panel_size, all_panel_sizes)
  gs_scores <- predict_with_ensemble(
    standardised$train[, panel, drop = FALSE],
    base$outcome[split$train_idx],
    standardised$test[, panel, drop = FALSE],
    random_seed = evaluation_seed
  )
  reconstructed_gs <- bench_auc(base$outcome[split$test_idx], gs_scores)
  reconstructed_random <- mean(random_panel_aucs(
    standardised$train, standardised$test,
    base$outcome[split$train_idx], base$outcome[split$test_idx],
    pool, panel_size, n_draws = 3L,
    seed = 99L + 1000L * repeat_idx + fold_idx,
    model_seed = evaluation_seed
  ))
  written_gs <- paired$GS_AUC[
    paired$repeat_idx == repeat_idx & paired$fold_idx == fold_idx &
      paired$k == panel_size
  ]
  written_random <- paired$Random_AUC[
    paired$repeat_idx == repeat_idx & paired$fold_idx == fold_idx &
      paired$k == panel_size
  ]
  reconstruction_rows <- append_row(reconstruction_rows, data.frame(
    dataset = dataset,
    repeat_idx = repeat_idx,
    fold_idx = fold_idx,
    k = panel_size,
    written_GS_AUC = written_gs,
    reconstructed_GS_AUC = reconstructed_gs,
    absolute_GS_difference = abs(written_gs - reconstructed_gs),
    written_Random_AUC = written_random,
    reconstructed_Random_AUC = reconstructed_random,
    absolute_Random_difference = abs(written_random - reconstructed_random)
  ))

  log_paths <- list.files(file.path(result_root, "logs"),
                          pattern = paste0("^", dataset, ".*[.]log$"),
                          full.names = TRUE)
  log_error <- any(vapply(log_paths, function(path) {
    lines <- readLines(path, warn = FALSE)
    any(grepl("Execution halted|Error in|Error:", lines))
  }, logical(1)))

  qa_rows <- append_row(qa_rows, data.frame(
    dataset = dataset,
    evaluation_rows = nrow(evaluation),
    unique_evaluation_keys = length(unique(keys)),
    ranking_files = length(ranking_paths),
    ranking_valid = ranking_valid,
    alpha_choice_valid = alpha_valid,
    configuration_valid = config_valid,
    source_hashes_valid = source_hash_valid,
    input_hashes_valid = input_hash_valid,
    log_error_found = log_error
  ))
}

summary_by_k <- do.call(rbind, summary_by_k)
primary_by_dataset <- do.call(rbind, primary_by_dataset)
fold_deltas <- do.call(rbind, fold_deltas)
alpha_rows <- do.call(rbind, alpha_rows)
qa_rows <- do.call(rbind, qa_rows)
reconstruction_rows <- do.call(rbind, reconstruction_rows)

stopifnot(
  nrow(summary_by_k) == 24L,
  nrow(primary_by_dataset) == 4L,
  nrow(fold_deltas) == 360L,
  nrow(alpha_rows) == 60L,
  all(!qa_rows$log_error_found),
  max(reconstruction_rows$absolute_GS_difference) < 1e-12,
  max(reconstruction_rows$absolute_Random_difference) < 1e-12
)

overall <- data.frame(
  datasets = nrow(primary_by_dataset),
  mean_GS_AUC = mean(primary_by_dataset$mean_GS_AUC),
  mean_Random_AUC = mean(primary_by_dataset$mean_Random_AUC),
  mean_primary_delta = mean(primary_by_dataset$primary_delta),
  datasets_positive = sum(primary_by_dataset$primary_delta > 0),
  datasets_negative = sum(primary_by_dataset$primary_delta < 0)
)

development <- read.csv(file.path(
  "redesign", "results_corrected", "full_benchmark_deterministic",
  "primary_comparisons.csv"
), stringsAsFactors = FALSE)
development_primary <- development[
  development$arm == "GS_full_ungrouped",
  c("dataset", "primary_mean_auc", "primary_delta")
]
names(development_primary)[2L] <- "mean_GS_AUC"
development_primary$evaluation_set <- "configuration_selection_7"
external_primary <- primary_by_dataset[
  , c("dataset", "mean_GS_AUC", "primary_delta")
]
external_primary$evaluation_set <- "additional_4"
combined_primary <- rbind(development_primary, external_primary)
evaluation_set_summary <- do.call(rbind, lapply(
  split(combined_primary, combined_primary$evaluation_set),
  function(block) data.frame(
    evaluation_set = block$evaluation_set[1L],
    datasets = nrow(block),
    mean_primary_delta = mean(block$primary_delta)
  )
))
evaluation_set_summary <- rbind(
  evaluation_set_summary,
  data.frame(evaluation_set = "combined_11", datasets = 11L,
             mean_primary_delta = mean(combined_primary$primary_delta))
)

# These earlier files used the closest available configuration name and a
# different Random/evaluator protocol. The table is an audit of why those
# provisional results must not be substituted for the frozen evaluation.
previous_rows <- list()
for (dataset in datasets) {
  previous_path <- Sys.glob(file.path(
    "independent_benchmark_runs", "validation", dataset,
    "*_kfold_trainfilter", "data", "nested_results.csv"
  ))
  if (length(previous_path) != 1L) next
  previous <- read.csv(previous_path, stringsAsFactors = FALSE)
  previous <- previous[
    previous$Evaluator == "ensemble" &
      previous$Method %in% c("GS_stab_util", "Random") &
      previous$k %in% primary_panel_sizes,
  ]
  previous_means <- aggregate(AUC ~ Method, data = previous, FUN = mean)
  previous_delta <- previous_means$AUC[
    previous_means$Method == "GS_stab_util"
  ] - previous_means$AUC[previous_means$Method == "Random"]
  exact_delta <- primary_by_dataset$primary_delta[
    primary_by_dataset$dataset == dataset
  ]
  previous_rows[[length(previous_rows) + 1L]] <- data.frame(
    dataset = dataset,
    previous_configuration = "GS_stab_util approximation",
    previous_primary_delta = previous_delta,
    exact_configuration = "GS_full_ungrouped",
    exact_primary_delta = exact_delta,
    exact_minus_previous = exact_delta - previous_delta,
    protocols_comparable = FALSE
  )
}
previous_comparison <- do.call(rbind, previous_rows)

random_sensitivity_rows <- list()
random_sensitivity_by_k <- list()
for (dataset in datasets) {
  sensitivity <- read.csv(file.path(
    validation_root, dataset, "random_baseline_sensitivity_30_draws.csv"
  ), stringsAsFactors = FALSE)
  stopifnot(nrow(sensitivity) == 1350L,
            all(sensitivity$draw %in% 1:30))
  random_30 <- aggregate(
    AUC ~ repeat_idx + fold_idx + k, data = sensitivity, FUN = mean
  )
  names(random_30)[4L] <- "Random_30_AUC"
  random_3 <- aggregate(
    AUC ~ repeat_idx + fold_idx + k,
    data = sensitivity[sensitivity$draw <= 3L, ], FUN = mean
  )
  names(random_3)[4L] <- "Random_3_reconstructed_AUC"
  exact <- fold_deltas[
    fold_deltas$dataset == dataset &
      fold_deltas$k %in% primary_panel_sizes,
    c("repeat_idx", "fold_idx", "k", "GS_AUC", "Random_AUC")
  ]
  sensitivity_paired <- Reduce(
    function(left, right) merge(
      left, right, by = c("repeat_idx", "fold_idx", "k"), sort = TRUE
    ),
    list(exact, random_30, random_3)
  )
  stopifnot(
    nrow(sensitivity_paired) == 45L,
    max(abs(sensitivity_paired$Random_AUC -
              sensitivity_paired$Random_3_reconstructed_AUC)) < 1e-12
  )
  sensitivity_paired$delta_3 <- sensitivity_paired$GS_AUC -
    sensitivity_paired$Random_AUC
  sensitivity_paired$delta_30 <- sensitivity_paired$GS_AUC -
    sensitivity_paired$Random_30_AUC
  random_sensitivity_rows[[length(random_sensitivity_rows) + 1L]] <-
    data.frame(
      dataset = dataset,
      mean_GS_AUC = mean(sensitivity_paired$GS_AUC),
      mean_Random_3_AUC = mean(sensitivity_paired$Random_AUC),
      delta_3 = mean(sensitivity_paired$delta_3),
      mean_Random_30_AUC = mean(sensitivity_paired$Random_30_AUC),
      delta_30 = mean(sensitivity_paired$delta_30),
      delta_30_minus_delta_3 = mean(sensitivity_paired$delta_30) -
        mean(sensitivity_paired$delta_3)
    )
  for (panel_size in primary_panel_sizes) {
    block <- sensitivity_paired[sensitivity_paired$k == panel_size, ]
    random_sensitivity_by_k[[length(random_sensitivity_by_k) + 1L]] <-
      data.frame(
        dataset = dataset,
        k = panel_size,
        mean_GS_AUC = mean(block$GS_AUC),
        mean_Random_30_AUC = mean(block$Random_30_AUC),
        delta_30 = mean(block$delta_30)
      )
  }
}
random_sensitivity <- do.call(rbind, random_sensitivity_rows)
random_sensitivity_by_k <- do.call(rbind, random_sensitivity_by_k)
random_sensitivity <- rbind(
  random_sensitivity,
  data.frame(
    dataset = "mean_across_4",
    mean_GS_AUC = mean(random_sensitivity$mean_GS_AUC),
    mean_Random_3_AUC = mean(random_sensitivity$mean_Random_3_AUC),
    delta_3 = mean(random_sensitivity$delta_3),
    mean_Random_30_AUC = mean(random_sensitivity$mean_Random_30_AUC),
    delta_30 = mean(random_sensitivity$delta_30),
    delta_30_minus_delta_3 = mean(
      random_sensitivity$delta_30_minus_delta_3
    )
  )
)

write.csv(summary_by_k, file.path(summary_root, "summary_by_dataset_k.csv"),
          row.names = FALSE)
write.csv(primary_by_dataset,
          file.path(summary_root, "primary_by_dataset.csv"), row.names = FALSE)
write.csv(fold_deltas, file.path(summary_root, "paired_fold_results.csv"),
          row.names = FALSE)
write.csv(alpha_rows, file.path(summary_root, "alpha_selection.csv"),
          row.names = FALSE)
write.csv(qa_rows, file.path(summary_root, "qa_checks.csv"), row.names = FALSE)
write.csv(reconstruction_rows,
          file.path(summary_root, "independent_reconstruction.csv"),
          row.names = FALSE)
write.csv(overall, file.path(summary_root, "overall_primary.csv"),
          row.names = FALSE)
write.csv(combined_primary,
          file.path(summary_root, "development_and_additional_datasets.csv"),
          row.names = FALSE)
write.csv(evaluation_set_summary,
          file.path(summary_root, "evaluation_set_summary.csv"),
          row.names = FALSE)
write.csv(previous_comparison,
          file.path(summary_root, "previous_approximation_audit.csv"),
          row.names = FALSE)
write.csv(random_sensitivity,
          file.path(summary_root, "random_30_draw_sensitivity.csv"),
          row.names = FALSE)
write.csv(random_sensitivity_by_k,
          file.path(summary_root, "random_30_draw_sensitivity_by_k.csv"),
          row.names = FALSE)

cat("Frozen external summary complete\n")
print(primary_by_dataset, row.names = FALSE, digits = 4)
print(overall, row.names = FALSE, digits = 4)
