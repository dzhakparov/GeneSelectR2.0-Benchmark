# Compare frozen independent conclusions with selected historical runs.
#
# Historical runs are selected by recency and complete row grids. Configuration
# differences are retained in the output and prevent direct numerical
# equivalence claims.

suppressPackageStartupMessages(library(dplyr))

output_root <- file.path(
  "independent_benchmark_runs",
  "historical_comparison"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

run_specs <- data.frame(
  dataset = c("SOS-ALL", "IMvigor210"),
  historical_run = c(
    "results_sosall/2026-08-13_kfold",
    "results_imvigor210/2026-08-08_kfold"
  ),
  independent_run = c(
    "independent_benchmark_runs/sosall/2026-08-18_kfold_trainfilter",
    "independent_benchmark_runs/imvigor210/2026-08-18_kfold_trainfilter"
  ),
  historical_candidate_pool = c(5000L, 5000L),
  independent_candidate_pool = c(2000L, 2000L),
  historical_filter_scope = c("complete cohort", "complete cohort"),
  independent_filter_scope = c("outer training partition",
                               "outer training partition"),
  historical_count_normalization = c(
    "array expression; no count normalization",
    "complete-cohort TMM/logCPM before outer CV"
  ),
  independent_count_normalization = c(
    "array expression; no count normalization",
    "outer-training fixed-reference TMM/logCPM"
  ),
  stringsAsFactors = FALSE
)

analyse_nested <- function(path, dataset, source) {
  results <- read.csv(path, stringsAsFactors = FALSE)
  primary <- results %>% filter(Evaluator == "ensemble")
  random <- primary %>%
    filter(Method == "Random") %>%
    select(Repeat, Fold, k, Random_AUC = AUC)
  paired <- primary %>%
    inner_join(random, by = c("Repeat", "Fold", "k")) %>%
    mutate(delta_AUC = AUC - Random_AUC)
  folds_per_repeat <- n_distinct(primary$Fold)

  paired %>%
    group_by(Method, k) %>%
    summarise(
      AUC_mean = mean(AUC),
      Random_AUC_mean = mean(Random_AUC),
      delta_AUC_mean = mean(delta_AUC),
      delta_AUC_sd = sd(delta_AUC),
      n_splits = n(),
      .groups = "drop"
    ) %>%
    mutate(
      corrected_se = sqrt(
        (1 / n_splits + 1 / (folds_per_repeat - 1)) * delta_AUC_sd^2
      ),
      corrected_ci_low = delta_AUC_mean -
        qt(0.975, df = n_splits - 1L) * corrected_se,
      corrected_ci_high = delta_AUC_mean +
        qt(0.975, df = n_splits - 1L) * corrected_se,
      dataset = dataset,
      source = source
    )
}

all_summaries <- list()
for (idx in seq_len(nrow(run_specs))) {
  spec <- run_specs[idx, ]
  historical_file <- file.path(spec$historical_run, "data", "nested_results.csv")
  independent_file <- file.path(spec$independent_run, "data", "nested_results.csv")
  if (!file.exists(historical_file) || !file.exists(independent_file)) {
    stop("Missing comparison input for ", spec$dataset)
  }
  all_summaries[[paste0(spec$dataset, "_historical")]] <- analyse_nested(
    historical_file, spec$dataset, "historical"
  )
  all_summaries[[paste0(spec$dataset, "_independent")]] <- analyse_nested(
    independent_file, spec$dataset, "independent"
  )
}

summary_long <- bind_rows(all_summaries)
write.csv(summary_long, file.path(output_root, "selected_run_summaries.csv"),
          row.names = FALSE)

semantic <- summary_long %>%
  filter(Method == "GS_semantic") %>%
  select(dataset, source, k, AUC_mean, Random_AUC_mean, delta_AUC_mean,
         corrected_ci_low, corrected_ci_high) %>%
  tidyr::pivot_wider(
    names_from = source,
    values_from = c(AUC_mean, Random_AUC_mean, delta_AUC_mean,
                    corrected_ci_low, corrected_ci_high)
  ) %>%
  mutate(
    AUC_independent_minus_historical =
      AUC_mean_independent - AUC_mean_historical,
    delta_independent_minus_historical =
      delta_AUC_mean_independent - delta_AUC_mean_historical
  ) %>%
  left_join(run_specs, by = "dataset") %>%
  arrange(dataset, k)
write.csv(semantic, file.path(output_root, "semantic_comparison.csv"),
          row.names = FALSE)

best_methods <- summary_long %>%
  group_by(dataset, source, k) %>%
  slice_max(AUC_mean, n = 1L, with_ties = FALSE) %>%
  ungroup() %>%
  select(dataset, source, k, Method, AUC_mean, delta_AUC_mean) %>%
  arrange(dataset, k, source)
write.csv(best_methods, file.path(output_root, "best_method_comparison.csv"),
          row.names = FALSE)

print(semantic %>%
        select(dataset, k, AUC_mean_historical, AUC_mean_independent,
               AUC_independent_minus_historical,
               delta_AUC_mean_historical, delta_AUC_mean_independent,
               delta_independent_minus_historical),
      row.names = FALSE)
cat("\nBest methods by run:\n")
print(best_methods, row.names = FALSE)
