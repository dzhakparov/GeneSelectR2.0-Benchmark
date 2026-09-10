#!/usr/bin/env Rscript

# Confirmation run for the utility-null instability diagnosis
# (redesign/results_corrected/older7_null_instability_2026-09-03/diagnosis.md).
#
# Rebuilds the r1f1 calibration null for ONE dataset with 100 permutations
# instead of 20, keeping every other argument identical to the 2026-09-02
# rebuild in redesign/run_older7_calibration_diagnostics.R, including the
# withr::with_seed(42L, ...) wrapper. The permutation loop seeds each
# permutation with random_seed + permutation_idx and the ambient RNG (which
# drives the cv.glmnet fold assignment inside each fit) is consumed
# sequentially from the seed-42 state, so permutations 1..20 of this run MUST
# reproduce the saved 20-permutation checkpoint row for row. That equality is
# asserted before the new null is accepted -- it is the validity gate that
# makes the 20-perm and 100-perm LOPO summaries comparable.
#
# Outputs (new files only; the saved 20-permutation checkpoint is evidence
# and is never overwritten):
#   <dataset_dir>/older7_calibration_null100_r1f1.rds
#   <dataset_dir>/older7_calibration_null100_validation.csv
#
# Usage:
#   Rscript redesign/run_older7_null100_confirmation.R <dataset> [budget]
#   dataset: GSE107994 or GSE13355 (the two LOPO-unstable datasets).
#   budget: seconds, default 604800 (same convention as the driver).

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
budget <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 604800
stopifnot(
  dataset %in% c("GSE107994", "GSE13355"),
  is.finite(budget), budget > 0
)

started_at <- proc.time()[["elapsed"]]
over_budget <- function() proc.time()[["elapsed"]] - started_at > budget

suppressPackageStartupMessages({
  library(glmnet)
  library(withr)
})
options(warn = 1)
for (path in list.files(file.path("package", "GeneSelectR", "R"),
                        full.names = TRUE)) {
  source(path)
}
source(file.path("redesign", "R", "bio_prior.R"))

dataset_dir <- file.path(
  "redesign", "results_corrected", "validation_benchmark", dataset
)
null_split_name <- "r1f1"
n_permutations_new <- 100L

split_path <- file.path(dataset_dir, "split_r1_f1.rds")
meta_path <- file.path(
  dataset_dir, "ranking_r1_f1_GS_full_ungrouped_meta.csv"
)
outcome_path <- file.path(dataset_dir, "base_data.rds")
saved_null_path <- file.path(
  dataset_dir, "older7_calibration_null_r1f1.rds"
)
new_null_path <- file.path(
  dataset_dir, "older7_calibration_null100_r1f1.rds"
)
summary_path <- file.path(
  dataset_dir, "older7_calibration_null100_validation.csv"
)
for (path in c(split_path, meta_path, outcome_path, saved_null_path)) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
}
if (file.exists(new_null_path)) {
  # The checkpoint is only written after every QA gate passes, so an existing
  # file means this run already completed. Never overwrite it.
  cat(sprintf("%s: %s exists; nothing to do\n", dataset, new_null_path))
  quit(save = "no")
}

metadata <- read.csv(meta_path, stringsAsFactors = FALSE)
alpha <- as.numeric(metadata$value[metadata$key == "alpha"])
alpha_tag <- if (identical(alpha, 0.5)) "0p5" else "1"
fit <- readRDS(file.path(dataset_dir, sprintf(
  "fit_r1_f1_GS_full_ungrouped_a%s.rds", alpha_tag
)))$fit

saved_null <- readRDS(saved_null_path)
utility_epsilon <- saved_null$qa$utility_epsilon

split <- readRDS(split_path)
pool <- split$pools$var2000
standardized <- standardise_split(
  split$train_raw[, pool, drop = FALSE],
  split$test_raw[, pool, drop = FALSE]
)
X <- standardized$train
outcome <- droplevels(as.factor(readRDS(outcome_path)$outcome))
y <- droplevels(outcome[split$train_idx])
scores <- fit$gene_scores[match(colnames(X), fit$gene_scores$gene), ,
                          drop = FALSE]
stopifnot(identical(as.character(scores$gene), colnames(X)))

# ---------------------------------------------------------------------------
# Helpers copied verbatim from redesign/run_older7_calibration_diagnostics.R
# so this run stays on the identical code path.
# ---------------------------------------------------------------------------

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

reconstruct_raw_mi <- function(X, y, B = 50L, seed = 42L) {
  subsamples <- create_subsamples(
    y, B = B, random_seed = seed, scheme = "kfold", k_folds = 5L
  )
  mi_matrix <- vapply(subsamples, function(subsample) {
    compute_mi_vectorized(
      X[subsample$train, , drop = FALSE], y[subsample$train],
      method = "discrete", n_bins = 5L
    )
  }, numeric(ncol(X)))
  rowMeans(mi_matrix)
}

leave_one_permutation_out <- function(null_matrix, epsilon) {
  n_permutations <- nrow(null_matrix)
  t(vapply(seq_len(n_permutations), function(permutation_idx) {
    calibrate_by_null(
      null_matrix[permutation_idx, ],
      null_matrix[-permutation_idx, , drop = FALSE],
      epsilon = epsilon, winsorize_at = 4
    )
  }, numeric(ncol(null_matrix))))
}

summarise_ratio_matrix <- function(ratios, pillar, source) {
  data.frame(
    dataset = dataset,
    repeat_idx = 1L,
    fold_idx = 1L,
    pillar = pillar,
    source = source,
    n_values = length(ratios),
    median_ratio = stats::median(ratios),
    mean_log2_ratio = mean(log2(ratios)),
    sd_log2_ratio = stats::sd(log2(ratios)),
    fraction_between_half_and_two = mean(ratios >= 0.5 & ratios <= 2),
    fraction_above_two = mean(ratios > 2),
    fraction_below_half = mean(ratios < 0.5),
    stringsAsFactors = FALSE
  )
}

# Observed raw utility, re-derived on the driver's code path. The epsilon
# recomputed from it must equal the saved one; both are checked.
shap_frequency <- unname(
  reconstruct_instance_frequency(fit, colnames(X))[colnames(X)]
)
raw_mi <- reconstruct_raw_mi(X, y, B = 50L, seed = 42L)
observed_raw_utility <- sqrt(shap_frequency * raw_mi + 1e-10)
positive_utility <- observed_raw_utility[
  observed_raw_utility > 0 & is.finite(observed_raw_utility)
]
recomputed_epsilon <- if (length(positive_utility) > 0L) {
  as.numeric(stats::quantile(positive_utility, 0.95, na.rm = TRUE))
} else {
  1e-6
}
if (abs(recomputed_epsilon - utility_epsilon) > 1e-15) {
  stop(sprintf(
    "Recomputed utility epsilon %.6g != saved %.6g; inputs drifted.",
    recomputed_epsilon, utility_epsilon
  ), call. = FALSE)
}

if (over_budget()) {
  cat("[budget] stop before 100-permutation null build\n")
  quit(save = "no")
}

# WHY the seed wrapper: identical to the 2026-09-02 rebuild -- the null
# builder seeds the outcome permutations and the subsample design but NOT the
# cv.glmnet fold assignment, which draws from the ambient RNG. Only
# n_permutations differs from the saved rebuild (100 vs 20); every other
# argument is copied from run_older7_calibration_diagnostics.R lines 547-565.
clear_run_cache()
null_result <- withr::with_seed(42L, compute_null_selection_frequencies(
  X, y,
  B = 20L,
  n_permutations = n_permutations_new,
  regularization_method = "elastic_net",
  alpha = alpha,
  q_max = NULL,
  subsample_scheme = "kfold",
  k_folds = 5L,
  random_seed = 42L,
  compute_utility_null = TRUE,
  utility_method = "instance_shap",
  utility_mi = TRUE,
  mi_method = "discrete",
  mi_bins = 5L,
  max_failed_fraction = 0.1,
  use_cache = FALSE,
  verbose = TRUE
))
if (!identical(dim(null_result$null_frequencies), c(100L, ncol(X))) ||
    !identical(dim(null_result$null_utility), c(100L, ncol(X))) ||
    any(!is.finite(null_result$null_frequencies)) ||
    any(!is.finite(null_result$null_utility))) {
  stop("Recomputed 100-permutation calibration null is invalid.",
       call. = FALSE)
}

# Validity gate: permutations 1..20 must reproduce the saved 20-permutation
# checkpoint exactly. If this fails, the RNG design drifted and the 20-perm
# vs 100-perm comparison is meaningless -- stop before writing anything.
first20_utility_diff <- max(abs(
  null_result$null_utility[1:20, ] - saved_null$null_utility
))
first20_frequency_diff <- max(abs(
  null_result$null_frequencies[1:20, ] - saved_null$null_frequencies
))
if (first20_utility_diff > 1e-12 || first20_frequency_diff > 1e-12) {
  stop(sprintf(
    paste0("First 20 permutations do not reproduce the saved null ",
           "(utility diff %.3g, frequency diff %.3g)."),
    first20_utility_diff, first20_frequency_diff
  ), call. = FALSE)
}

null_utility_ratios <- leave_one_permutation_out(
  null_result$null_utility, epsilon = utility_epsilon
)
null_stability_ratios <- leave_one_permutation_out(
  null_result$null_frequencies, epsilon = 0.01
)
recomputed_utility <- calibrate_by_null(
  observed_raw_utility, null_result$null_utility,
  epsilon = utility_epsilon, winsorize_at = 4
)
recomputed_stability <- calibrate_by_null(
  scores$pi_raw, null_result$null_frequencies,
  epsilon = 0.01, winsorize_at = 4
)
null_summary <- rbind(
  summarise_ratio_matrix(
    null_stability_ratios, "stability", "leave_one_permutation_out_null100"
  ),
  summarise_ratio_matrix(
    null_utility_ratios, "utility", "leave_one_permutation_out_null100"
  ),
  summarise_ratio_matrix(
    recomputed_stability, "stability", "observed_outcome_null100"
  ),
  summarise_ratio_matrix(
    recomputed_utility, "utility", "observed_outcome_null100"
  )
)
qa <- data.frame(
  dataset = dataset,
  repeat_idx = 1L,
  fold_idx = 1L,
  selected_alpha = alpha,
  n_permutations = n_permutations_new,
  utility_epsilon = utility_epsilon,
  first20_utility_max_abs_difference = first20_utility_diff,
  first20_frequency_max_abs_difference = first20_frequency_diff,
  null_n_fits = null_result$diagnostics$n_fits,
  null_n_failed_fits = null_result$diagnostics$n_failed_fits,
  null_failure_fraction = null_result$diagnostics$failure_fraction,
  passed = first20_utility_diff <= 1e-12 &
    first20_frequency_diff <= 1e-12 &
    null_result$diagnostics$failure_fraction <= 0.1,
  stringsAsFactors = FALSE
)
if (!qa$passed) {
  stop("100-permutation calibration null failed validation.", call. = FALSE)
}
saveRDS(
  list(
    null_frequencies = null_result$null_frequencies,
    null_utility = null_result$null_utility,
    null_summary = null_summary,
    qa = qa,
    diagnostics = null_result$diagnostics
  ),
  new_null_path, version = 3
)
write.csv(null_summary, summary_path, row.names = FALSE)
cat(sprintf(
  paste0("[%s %s] 100-permutation null built; first-20 reproduction exact; ",
         "utility LOPO sd_log2 %.3f, median %.4g, %.1f%% within [0.5, 2]\n"),
  dataset, null_split_name,
  null_summary$sd_log2_ratio[null_summary$pillar == "utility" &
                             null_summary$source ==
                               "leave_one_permutation_out_null100"],
  null_summary$median_ratio[null_summary$pillar == "utility" &
                            null_summary$source ==
                              "leave_one_permutation_out_null100"],
  100 * null_summary$fraction_between_half_and_two[
    null_summary$pillar == "utility" &
      null_summary$source == "leave_one_permutation_out_null100"
  ]
))
