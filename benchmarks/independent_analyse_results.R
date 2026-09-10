# Independent analysis of newly generated benchmark rows.
#
# The input root must be the dedicated independent run directory. Historical
# result directories are not searched. The ensemble evaluator is primary;
# component evaluators are retained as a sensitivity analysis. Comparisons with
# Random are paired by repeat, fold, and panel size.

suppressPackageStartupMessages({
  library(dplyr)
})

run_root <- Sys.getenv("GS_OUTPUT_ROOT", "independent_benchmark_runs")
run_root <- normalizePath(run_root, mustWork = TRUE)
if (!grepl("independent_benchmark_runs$", run_root)) {
  stop("GS_OUTPUT_ROOT must identify the dedicated independent_benchmark_runs directory.")
}

analysis_root <- file.path(run_root, "independent_analysis")
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

result_files <- list.files(
  run_root, pattern = "^nested_results\\.csv$", recursive = TRUE,
  full.names = TRUE
)
result_files <- result_files[!grepl("/failed_runs/", result_files, fixed = TRUE)]
if (length(result_files) == 0L) {
  stop("No completed independent nested_results.csv files were found.")
}

analyse_one <- function(result_file) {
  data_dir <- dirname(result_file)
  run_dir <- dirname(data_dir)
  relative <- substring(run_dir, nchar(run_root) + 2L)
  dataset_id <- gsub("/", "__", relative, fixed = TRUE)
  output_dir <- file.path(analysis_root, dataset_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- read.csv(result_file, stringsAsFactors = FALSE)
  required <- c("Method", "Evaluator", "Repeat", "Fold", "k",
                "AUC", "BalAcc", "MCC")
  missing_columns <- setdiff(required, colnames(results))
  if (length(missing_columns) > 0L) {
    stop(dataset_id, ": missing columns: ", paste(missing_columns, collapse = ", "))
  }

  key <- results[c("Method", "Evaluator", "Repeat", "Fold", "k")]
  duplicate_rows <- duplicated(key)
  expected_split_keys <- expand.grid(
    Repeat = seq_len(3L), Fold = seq_len(5L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  observed_split_keys <- unique(results[c("Repeat", "Fold")])
  missing_split_keys <- dplyr::anti_join(
    expected_split_keys, observed_split_keys, by = c("Repeat", "Fold")
  )
  expected_grid <- expand.grid(
    Method = unique(results$Method),
    Evaluator = unique(results$Evaluator),
    Repeat = seq_len(3L), Fold = seq_len(5L),
    k = sort(unique(results$k)),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  missing_grid_rows <- nrow(dplyr::anti_join(
    expected_grid, key,
    by = c("Method", "Evaluator", "Repeat", "Fold", "k")
  ))
  valid_ranges <- all(results$AUC[is.finite(results$AUC)] >= 0 &
                        results$AUC[is.finite(results$AUC)] <= 1) &&
    all(results$BalAcc[is.finite(results$BalAcc)] >= 0 &
          results$BalAcc[is.finite(results$BalAcc)] <= 1) &&
    all(results$MCC[is.finite(results$MCC)] >= -1 &
          results$MCC[is.finite(results$MCC)] <= 1)

  integrity <- data.frame(
    dataset = dataset_id,
    rows = nrow(results),
    methods = n_distinct(results$Method),
    evaluators = n_distinct(results$Evaluator),
    repeats = n_distinct(results$Repeat),
    folds = n_distinct(paste(results$Repeat, results$Fold)),
    panel_sizes = paste(sort(unique(results$k)), collapse = ","),
    expected_rows = nrow(expected_grid),
    missing_grid_rows = missing_grid_rows,
    missing_split_keys = nrow(missing_split_keys),
    duplicate_keys = sum(duplicate_rows),
    valid_metric_ranges = valid_ranges,
    missing_auc = sum(!is.finite(results$AUC)),
    missing_balacc = sum(!is.finite(results$BalAcc)),
    missing_mcc = sum(!is.finite(results$MCC)),
    stringsAsFactors = FALSE
  )
  write.csv(integrity, file.path(output_dir, "integrity.csv"), row.names = FALSE)
  if (any(duplicate_rows) || nrow(missing_split_keys) > 0L ||
      missing_grid_rows > 0L || !valid_ranges) {
    stop(dataset_id, ": result-integrity checks failed.")
  }

  failures_file <- file.path(data_dir, "method_failures.csv")
  failures <- if (file.exists(failures_file)) {
    read.csv(failures_file, stringsAsFactors = FALSE)
  } else {
    data.frame(Method = character(0), error = character(0))
  }
  write.csv(failures, file.path(output_dir, "method_failures.csv"),
            row.names = FALSE)

  evaluator_sensitivity <- results %>%
    group_by(Method, Evaluator, k) %>%
    summarise(n = sum(is.finite(AUC)),
              AUC_mean = mean(AUC, na.rm = TRUE),
              AUC_sd = sd(AUC, na.rm = TRUE),
              BalAcc_mean = mean(BalAcc, na.rm = TRUE),
              MCC_mean = mean(MCC, na.rm = TRUE),
              .groups = "drop")
  write.csv(evaluator_sensitivity,
            file.path(output_dir, "evaluator_sensitivity.csv"), row.names = FALSE)

  primary <- results %>% filter(Evaluator == "ensemble")
  folds_per_repeat <- n_distinct(primary$Fold)
  if (folds_per_repeat < 2L) {
    stop(dataset_id, ": at least two folds per repeat are required.")
  }
  random <- primary %>%
    filter(Method == "Random") %>%
    select(Repeat, Fold, k,
           Random_AUC = AUC,
           Random_BalAcc = BalAcc,
           Random_MCC = MCC)
  if (nrow(random) == 0L) stop(dataset_id, ": Random ensemble rows are absent.")

  paired <- primary %>%
    inner_join(random, by = c("Repeat", "Fold", "k")) %>%
    mutate(AUC_minus_Random = AUC - Random_AUC,
           BalAcc_minus_Random = BalAcc - Random_BalAcc,
           MCC_minus_Random = MCC - Random_MCC)
  write.csv(paired, file.path(output_dir, "paired_ensemble_rows.csv"),
            row.names = FALSE)

  repeat_effects <- paired %>%
    group_by(Method, k, Repeat) %>%
    summarise(repeat_mean_delta_auc = mean(AUC_minus_Random, na.rm = TRUE),
              .groups = "drop")

  repeat_ci <- repeat_effects %>%
    group_by(Method, k) %>%
    summarise(n_repeats = sum(is.finite(repeat_mean_delta_auc)),
              repeat_delta_sd = sd(repeat_mean_delta_auc, na.rm = TRUE),
              repeat_delta_se = repeat_delta_sd / sqrt(n_repeats),
              repeat_ci_low = mean(repeat_mean_delta_auc, na.rm = TRUE) -
                qt(0.975, df = pmax(n_repeats - 1L, 1L)) * repeat_delta_se,
              repeat_ci_high = mean(repeat_mean_delta_auc, na.rm = TRUE) +
                qt(0.975, df = pmax(n_repeats - 1L, 1L)) * repeat_delta_se,
              .groups = "drop")

  primary_summary <- paired %>%
    group_by(Method, k) %>%
    summarise(n_splits = sum(is.finite(AUC)),
              AUC_mean = mean(AUC, na.rm = TRUE),
              AUC_sd = sd(AUC, na.rm = TRUE),
              BalAcc_mean = mean(BalAcc, na.rm = TRUE),
              MCC_mean = mean(MCC, na.rm = TRUE),
              Random_AUC_mean = mean(Random_AUC, na.rm = TRUE),
              delta_AUC_mean = mean(AUC_minus_Random, na.rm = TRUE),
              delta_AUC_sd = sd(AUC_minus_Random, na.rm = TRUE),
              delta_AUC_median = median(AUC_minus_Random, na.rm = TRUE),
              splits_above_random = sum(AUC_minus_Random > 0, na.rm = TRUE),
              splits_equal_random = sum(AUC_minus_Random == 0, na.rm = TRUE),
              splits_below_random = sum(AUC_minus_Random < 0, na.rm = TRUE),
              .groups = "drop") %>%
    # Repeated-CV folds share training observations. The corrected resampled
    # standard error inflates the ordinary split-level error by the approximate
    # test-to-training ratio, 1 / (K - 1), and is reported alongside the three-
    # repeat interval. Both summaries expose the limited effective sample size.
    mutate(corrected_se = sqrt(
             (1 / n_splits + 1 / (folds_per_repeat - 1)) * delta_AUC_sd^2
           ),
           corrected_ci_low = delta_AUC_mean -
             qt(0.975, df = pmax(n_splits - 1L, 1L)) * corrected_se,
           corrected_ci_high = delta_AUC_mean +
             qt(0.975, df = pmax(n_splits - 1L, 1L)) * corrected_se,
           corrected_p_value = ifelse(
             Method == "Random" | !is.finite(corrected_se) | corrected_se == 0,
             NA_real_,
             2 * pt(-abs(delta_AUC_mean / corrected_se),
                    df = pmax(n_splits - 1L, 1L))
           )) %>%
    group_by(k) %>%
    mutate(corrected_p_bh = p.adjust(corrected_p_value, method = "BH")) %>%
    ungroup() %>%
    left_join(repeat_ci, by = c("Method", "k")) %>%
    arrange(k, desc(delta_AUC_mean), desc(AUC_mean))
  write.csv(primary_summary, file.path(output_dir, "primary_summary.csv"),
            row.names = FALSE)

  # Signed-rank tests are retained as exploratory screens. Repeated CV folds
  # reuse observations, so effect sizes and repeat-level intervals are primary.
  exploratory_tests <- paired %>%
    group_by(Method, k) %>%
    summarise(
      n_pairs = sum(is.finite(AUC_minus_Random)),
      p_value = if (Method[1] == "Random" || n_pairs < 2L) NA_real_ else
        suppressWarnings(stats::wilcox.test(
          AUC_minus_Random, mu = 0, exact = FALSE
        )$p.value),
      .groups = "drop"
    ) %>%
    group_by(k) %>%
    mutate(p_adjusted_bh = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  write.csv(exploratory_tests,
            file.path(output_dir, "exploratory_vs_random.csv"), row.names = FALSE)

  headline <- primary %>%
    filter(Method == "GS_semantic") %>%
    select(Repeat, Fold, k,
           headline_AUC = AUC,
           headline_BalAcc = BalAcc,
           headline_MCC = MCC)
  headline_paired <- primary %>%
    filter(Method != "GS_semantic") %>%
    inner_join(headline, by = c("Repeat", "Fold", "k")) %>%
    mutate(delta_AUC = headline_AUC - AUC,
           delta_BalAcc = headline_BalAcc - BalAcc,
           delta_MCC = headline_MCC - MCC)

  headline_repeat <- headline_paired %>%
    group_by(Method, k, Repeat) %>%
    summarise(repeat_mean_delta_auc = mean(delta_AUC, na.rm = TRUE),
              .groups = "drop")
  headline_repeat_ci <- headline_repeat %>%
    group_by(Method, k) %>%
    summarise(n_repeats = sum(is.finite(repeat_mean_delta_auc)),
              repeat_delta_sd = sd(repeat_mean_delta_auc, na.rm = TRUE),
              repeat_delta_se = repeat_delta_sd / sqrt(n_repeats),
              repeat_ci_low = mean(repeat_mean_delta_auc, na.rm = TRUE) -
                qt(0.975, df = pmax(n_repeats - 1L, 1L)) * repeat_delta_se,
              repeat_ci_high = mean(repeat_mean_delta_auc, na.rm = TRUE) +
                qt(0.975, df = pmax(n_repeats - 1L, 1L)) * repeat_delta_se,
              .groups = "drop")

  headline_pairwise <- headline_paired %>%
    group_by(Method, k) %>%
    summarise(n_pairs = sum(is.finite(delta_AUC)),
              headline_AUC_mean = mean(headline_AUC, na.rm = TRUE),
              comparator_AUC_mean = mean(AUC, na.rm = TRUE),
              delta_AUC_mean = mean(delta_AUC, na.rm = TRUE),
              delta_AUC_sd = sd(delta_AUC, na.rm = TRUE),
              headline_wins = sum(delta_AUC > 0, na.rm = TRUE),
              ties = sum(delta_AUC == 0, na.rm = TRUE),
              headline_losses = sum(delta_AUC < 0, na.rm = TRUE),
              p_value = if (n_pairs < 2L) NA_real_ else suppressWarnings(
                stats::wilcox.test(delta_AUC, mu = 0, exact = FALSE)$p.value
              ),
              .groups = "drop") %>%
    mutate(corrected_se = sqrt(
             (1 / n_pairs + 1 / (folds_per_repeat - 1)) * delta_AUC_sd^2
           ),
           corrected_ci_low = delta_AUC_mean -
             qt(0.975, df = pmax(n_pairs - 1L, 1L)) * corrected_se,
           corrected_ci_high = delta_AUC_mean +
             qt(0.975, df = pmax(n_pairs - 1L, 1L)) * corrected_se,
           corrected_p_value = ifelse(
             !is.finite(corrected_se) | corrected_se == 0,
             NA_real_,
             2 * pt(-abs(delta_AUC_mean / corrected_se),
                    df = pmax(n_pairs - 1L, 1L))
           )) %>%
    left_join(headline_repeat_ci, by = c("Method", "k")) %>%
    group_by(k) %>%
    mutate(p_adjusted_bh = p.adjust(p_value, method = "BH"),
           corrected_p_bh = p.adjust(corrected_p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(k, desc(delta_AUC_mean))
  write.csv(headline_pairwise,
            file.path(output_dir, "headline_pairwise.csv"), row.names = FALSE)

  data.frame(dataset = dataset_id,
             completed = TRUE,
             methods = n_distinct(primary$Method),
             splits = n_distinct(paste(primary$Repeat, primary$Fold)),
             failures = n_distinct(failures$Method),
             stringsAsFactors = FALSE)
}

run_index <- bind_rows(lapply(result_files, analyse_one))
write.csv(run_index, file.path(analysis_root, "run_index.csv"), row.names = FALSE)
print(run_index, row.names = FALSE)
cat("\nIndependent analyses written to:", analysis_root, "\n")
