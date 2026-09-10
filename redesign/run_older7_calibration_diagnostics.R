#!/usr/bin/env Rscript

# Calibration and component-dependence diagnostics for the current
# GeneSelectR configuration on the seven older datasets. The extract stage
# checks every saved outer-fold fit. The null_check stage independently rebuilds
# the permutation null for the fixed r1f1 split and uses leave-one-permutation-
# out ratios to assess null centering without using the observed outcome fit as
# its own reference.

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
  stage %in% c("init", "extract", "null_check", "assemble", "all"),
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

saved_method <- if (is_validation) "GS_full_ungrouped" else "full_ungrouped"
extract_grid <- expand.grid(
  repeat_idx = 1:3, fold_idx = 1:5,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
extract_override <- Sys.getenv("GENESELECTR_VALIDATION_SPLIT", "")
if (nzchar(extract_override)) {
  parsed <- regmatches(
    extract_override, regexec("^r([1-3])f([1-5])$", extract_override)
  )[[1L]]
  if (length(parsed) != 3L) {
    stop("GENESELECTR_VALIDATION_SPLIT must have the form r1f1.",
         call. = FALSE)
  }
  extract_grid <- data.frame(
    repeat_idx = as.integer(parsed[[2L]]),
    fold_idx = as.integer(parsed[[3L]])
  )
}
null_split_name <- Sys.getenv("GENESELECTR_CALIBRATION_SPLIT", "r1f1")
null_parsed <- regmatches(
  null_split_name, regexec("^r([1-3])f([1-5])$", null_split_name)
)[[1L]]
if (length(null_parsed) != 3L) {
  stop("GENESELECTR_CALIBRATION_SPLIT must have the form r1f1.",
       call. = FALSE)
}
null_repeat <- as.integer(null_parsed[[2L]])
null_fold <- as.integer(null_parsed[[3L]])

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
extract_path <- function(repeat_idx, fold_idx) file.path(
  dataset_dir, sprintf(
    "older7_calibration_extract_r%d_f%d.rds", repeat_idx, fold_idx
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
all_extract_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) {
    extract_path(repeat_idx, fold_idx)
  }, character(1))
}), use.names = FALSE)
null_checkpoint <- file.path(dataset_dir, sprintf(
  "older7_calibration_null_%s.rds", null_split_name
))
outcome_path <- file.path(
  dataset_dir, if (is_validation) "base_data.rds" else "base_outcome.rds"
)
gene_output <- file.path(
  dataset_dir, "older7_calibration_component_gene_scores.csv"
)
dependence_output <- file.path(
  dataset_dir, "older7_calibration_dependence_by_split.csv"
)
summary_output <- file.path(
  dataset_dir, "older7_calibration_dependence_summary.csv"
)
null_output <- file.path(
  dataset_dir, "older7_calibration_null_validation.csv"
)
qa_output <- file.path(
  dataset_dir, "older7_calibration_diagnostics_QA.csv"
)
if (!all(file.exists(c(all_fit_files, outcome_path)))) {
  stop("Selected GeneSelectR fits are missing.", call. = FALSE)
}

extension <- sprintf("older7_calibration_diagnostics_%s_v1", null_split_name)
prepare_redesign_extension(
  dataset_dir,
  extension = extension,
  source_files = c(
    "redesign/run_older7_calibration_diagnostics.R",
    list.files(file.path("package", "GeneSelectR", "R"), full.names = TRUE),
    file.path("redesign", "R", c("bio_prior.R", "run_provenance.R")),
    file.path(split_dir, "run_manifest.rds"),
    all_split_files, all_ranking_files, all_meta_files, all_fit_files,
    outcome_path
  ),
  config = list(
    dataset = dataset,
    observed_splits = "3 repeats x 5 folds",
    null_validation_split = null_split_name,
    null_validation_rule = "same fixed split for every dataset",
    null_B = 20L,
    null_permutations = 20L,
    null_cross_calibration = "leave_one_permutation_out",
    stability_epsilon = 0.01,
    utility_epsilon = "95th percentile of positive observed raw utility",
    winsorized_absolute_log2_ratio = 4,
    dependence_measure = "Spearman"
  ),
  output_files = c(
    all_extract_files, null_checkpoint, gene_output, dependence_output,
    summary_output, null_output, qa_output
  )
)
if (stage == "init") {
  cat(sprintf(
    "%s older-seven calibration diagnostics initialized; null split=%s\n",
    dataset, null_split_name
  ))
  quit(save = "no")
}

outcome <- if (is_validation) {
  readRDS(outcome_path)$outcome
} else {
  readRDS(outcome_path)
}
outcome <- droplevels(as.factor(outcome))

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

safe_spearman <- function(x, y, keep = rep(TRUE, length(x))) {
  x <- x[keep]
  y <- y[keep]
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 3L || length(unique(x[valid])) < 2L ||
      length(unique(y[valid])) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[valid], y[valid], method = "spearman"))
}

top_jaccard <- function(x, y, k) {
  x_top <- head(x, min(k, length(x)))
  y_top <- head(y, min(k, length(y)))
  length(intersect(x_top, y_top)) / length(union(x_top, y_top))
}

calibration_summary_difference <- function(ratios, stored, pillar) {
  calculated <- summarise_calibration(ratios, pillar, verbose = FALSE)
  fields <- c(
    "median_ratio", "max_ratio", "frac_above_2", "frac_below_half",
    "frac_uninformative", "sd_log2_ratio"
  )
  max(abs(
    as.numeric(unlist(calculated[fields])) -
      as.numeric(unlist(stored[fields]))
  ))
}

read_current_fit <- function(repeat_idx, fold_idx) {
  wrapper <- readRDS(fit_path(repeat_idx, fold_idx))
  selected_alpha <- read_selected_alpha(repeat_idx, fold_idx)
  fit <- wrapper$fit
  if (is.null(fit) || !isTRUE(all.equal(wrapper$alpha, selected_alpha)) ||
      !identical(fit$parameters$components, c("stability", "utility")) ||
      !identical(fit$parameters$gate_method, "none") ||
      !identical(fit$parameters$calibration_mode, "evidence_ratio") ||
      !identical(fit$parameters$utility_method, "instance_shap") ||
      !identical(fit$parameters$subsample_scheme, "kfold") ||
      !identical(fit$parameters$subsample_k_folds, 5)) {
    stop("Saved fit is not the current benchmark configuration.",
         call. = FALSE)
  }
  list(wrapper = wrapper, fit = fit, alpha = selected_alpha)
}

if (stage %in% c("extract", "all")) {
  for (job_idx in seq_len(nrow(extract_grid))) {
    repeat_idx <- extract_grid$repeat_idx[[job_idx]]
    fold_idx <- extract_grid$fold_idx[[job_idx]]
    output_path <- extract_path(repeat_idx, fold_idx)
    if (file.exists(output_path)) next
    if (over_budget()) {
      cat("[budget] stop in calibration extraction\n")
      quit(save = "no")
    }

    split <- readRDS(split_path(repeat_idx, fold_idx))
    pool <- split$pools$var2000
    current <- read_current_fit(repeat_idx, fold_idx)
    fit <- current$fit
    scores <- fit$gene_scores
    saved_ranking <- read.csv(
      ranking_path(repeat_idx, fold_idx), stringsAsFactors = FALSE
    )
    required <- c(
      "gene", "raw_score", "pi_raw", "pi_scored", "u", "u_mi", "u_scored"
    )
    if (!all(required %in% names(scores)) || anyDuplicated(scores$gene) ||
        !identical(as.character(scores$gene), as.character(saved_ranking$gene)) ||
        !setequal(scores$gene, pool)) {
      stop("Saved fit, ranking, and candidate pool are inconsistent.",
           call. = FALSE)
    }
    if (any(!is.finite(as.matrix(scores[setdiff(required, "gene")]))) ||
        any(scores$pi_scored <= 0) || any(scores$u_scored <= 0)) {
      stop("Component scores are invalid.", call. = FALSE)
    }

    shap_frequency_named <- reconstruct_instance_frequency(fit, scores$gene)
    shap_frequency <- unname(shap_frequency_named[scores$gene])
    shap_percentile <- percentile01(shap_frequency)
    reconstructed_u <- percentile01(sqrt(
      shap_percentile * scores$u_mi + 1e-10
    ))
    shap_reconstruction_difference <- max(abs(reconstructed_u - scores$u))
    # Rows of the saved selection matrix follow the candidate-pool order used
    # to fit the model. The ranking table is sorted by the final score, so it
    # must be matched back to the pool before frequencies are compared.
    selection_frequency_difference <- max(abs(
      rowMeans(fit$stability$selection_matrix) -
        scores$pi_raw[match(pool, scores$gene)]
    ))
    stability_summary_difference <- calibration_summary_difference(
      scores$pi_scored, fit$calibration$stability, "stability"
    )
    utility_summary_difference <- calibration_summary_difference(
      scores$u_scored, fit$calibration$utility, "utility"
    )
    if (shap_reconstruction_difference > 1e-12 ||
        selection_frequency_difference > 1e-12 ||
        stability_summary_difference > 1e-12 ||
        utility_summary_difference > 1e-12) {
      stop("Saved component reconstruction failed.", call. = FALSE)
    }

    recurrence_order <- scores$gene[order(
      -scores$pi_scored, -scores$pi_raw, scores$gene, method = "radix"
    )]
    utility_order <- scores$gene[order(
      -scores$u_scored, -scores$u, scores$gene, method = "radix"
    )]
    selected <- scores$pi_raw > 0
    null_diagnostics <- fit$calibration$null_diagnostics
    dependence <- data.frame(
      dataset = dataset,
      repeat_idx = repeat_idx,
      fold_idx = fold_idx,
      selected_alpha = current$alpha,
      n_genes = nrow(scores),
      n_ever_selected = sum(selected),
      spearman_calibrated_recurrence_utility_all = safe_spearman(
        log2(scores$pi_scored), log2(scores$u_scored)
      ),
      spearman_calibrated_recurrence_utility_selected = safe_spearman(
        log2(scores$pi_scored), log2(scores$u_scored), selected
      ),
      spearman_raw_recurrence_SHAP_all = safe_spearman(
        scores$pi_raw, shap_frequency
      ),
      spearman_SHAP_MI_all = safe_spearman(
        shap_frequency, scores$u_mi
      ),
      spearman_calibrated_utility_MI_all = safe_spearman(
        scores$u_scored, scores$u_mi
      ),
      recurrence_utility_jaccard_k10 = top_jaccard(
        recurrence_order, utility_order, 10L
      ),
      recurrence_utility_jaccard_k50 = top_jaccard(
        recurrence_order, utility_order, 50L
      ),
      recurrence_utility_jaccard_k200 = top_jaccard(
        recurrence_order, utility_order, 200L
      ),
      stability_median_ratio = fit$calibration$stability$median_ratio,
      stability_sd_log2_ratio = fit$calibration$stability$sd_log2_ratio,
      stability_frac_above_2 = fit$calibration$stability$frac_above_2,
      stability_frac_below_half = fit$calibration$stability$frac_below_half,
      utility_median_ratio = fit$calibration$utility$median_ratio,
      utility_sd_log2_ratio = fit$calibration$utility$sd_log2_ratio,
      utility_frac_above_2 = fit$calibration$utility$frac_above_2,
      utility_frac_below_half = fit$calibration$utility$frac_below_half,
      null_n_fits = null_diagnostics$n_fits,
      null_n_failed_fits = null_diagnostics$n_failed_fits,
      null_failure_fraction = null_diagnostics$failure_fraction,
      stringsAsFactors = FALSE
    )
    gene_level <- data.frame(
      dataset = dataset,
      repeat_idx = repeat_idx,
      fold_idx = fold_idx,
      gene = scores$gene,
      pi_raw = scores$pi_raw,
      pi_scored = scores$pi_scored,
      SHAP_instance_frequency = shap_frequency,
      MI_percentile = scores$u_mi,
      SHAPxMI_percentile = scores$u,
      utility_scored = scores$u_scored,
      current_raw_score = scores$raw_score,
      stringsAsFactors = FALSE
    )
    qa <- data.frame(
      dataset = dataset,
      repeat_idx = repeat_idx,
      fold_idx = fold_idx,
      SHAP_utility_max_abs_difference = shap_reconstruction_difference,
      selection_frequency_max_abs_difference =
        selection_frequency_difference,
      stability_summary_max_abs_difference = stability_summary_difference,
      utility_summary_max_abs_difference = utility_summary_difference,
      ratios_finite_positive = all(
        is.finite(scores$pi_scored) & scores$pi_scored > 0 &
          is.finite(scores$u_scored) & scores$u_scored > 0
      ),
      null_failure_fraction = null_diagnostics$failure_fraction,
      passed = shap_reconstruction_difference <= 1e-12 &
        selection_frequency_difference <= 1e-12 &
        stability_summary_difference <= 1e-12 &
        utility_summary_difference <= 1e-12 &
        all(is.finite(scores$pi_scored) & scores$pi_scored > 0) &
        all(is.finite(scores$u_scored) & scores$u_scored > 0) &
        null_diagnostics$failure_fraction <= fit$parameters$max_failed_fraction,
      stringsAsFactors = FALSE
    )
    if (!qa$passed) {
      stop("Calibration extraction failed validation.", call. = FALSE)
    }
    saveRDS(
      list(gene_level = gene_level, dependence = dependence, qa = qa),
      output_path, version = 3
    )
    cat(sprintf("[%s r%d f%d] calibration diagnostics extracted\n",
                dataset, repeat_idx, fold_idx))
  }
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
    repeat_idx = null_repeat,
    fold_idx = null_fold,
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

if (stage %in% c("null_check", "all") && !file.exists(null_checkpoint)) {
  if (over_budget()) {
    cat("[budget] stop before calibration null check\n")
    quit(save = "no")
  }
  split <- readRDS(split_path(null_repeat, null_fold))
  pool <- split$pools$var2000
  standardized <- standardise_split(
    split$train_raw[, pool, drop = FALSE],
    split$test_raw[, pool, drop = FALSE]
  )
  X <- standardized$train
  y <- droplevels(outcome[split$train_idx])
  current <- read_current_fit(null_repeat, null_fold)
  fit <- current$fit
  scores <- fit$gene_scores
  if (!setequal(colnames(X), scores$gene)) {
    stop("Null-check matrix and saved fit use different candidate genes.",
         call. = FALSE)
  }
  scores <- scores[match(colnames(X), scores$gene), , drop = FALSE]
  shap_frequency <- reconstruct_instance_frequency(fit, scores$gene)
  shap_frequency <- unname(shap_frequency[scores$gene])
  raw_mi <- reconstruct_raw_mi(X, y, B = 50L, seed = 42L)
  saved_mi_reconstruction <- percentile01(raw_mi)
  mi_max_abs_difference <- max(abs(saved_mi_reconstruction - scores$u_mi))
  if (mi_max_abs_difference > 1e-12) {
    stop("Raw mutual-information reconstruction does not reproduce saved MI.",
         call. = FALSE)
  }
  observed_raw_utility <- sqrt(shap_frequency * raw_mi + 1e-10)
  positive_utility <- observed_raw_utility[
    observed_raw_utility > 0 & is.finite(observed_raw_utility)
  ]
  utility_epsilon <- if (length(positive_utility) > 0L) {
    as.numeric(stats::quantile(positive_utility, 0.95, na.rm = TRUE))
  } else {
    1e-6
  }

  clear_run_cache()
  # WHY the seed wrapper: the null builder seeds the outcome permutations and
  # the subsample design (random_seed + permutation_idx) but NOT the cv.glmnet
  # fold assignment inside each permuted fit, which draws from the ambient RNG.
  # The saved benchmark fits ran the whole pipeline under local_seed(42) with
  # the observed subsample fits on PSOCK workers, so the master RNG was still
  # at the seed-42 state when the null loop started. An independent rebuild
  # must recreate that state, otherwise the rebuilt null -- and every
  # calibrated ratio -- differs by Monte Carlo noise in the CV folds. Measured
  # on GSE101794 r1f1 without the wrapper: max abs difference 4.40 (stability)
  # and 1.29 (utility) against the saved ratios; with the wrapper: exactly 0.
  null_result <- withr::with_seed(42L, compute_null_selection_frequencies(
    X, y,
    B = 20L,
    n_permutations = 20L,
    regularization_method = "elastic_net",
    alpha = current$alpha,
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
  if (!identical(dim(null_result$null_frequencies), c(20L, ncol(X))) ||
      !identical(dim(null_result$null_utility), c(20L, ncol(X))) ||
      any(!is.finite(null_result$null_frequencies)) ||
      any(!is.finite(null_result$null_utility))) {
    stop("Recomputed calibration null is invalid.", call. = FALSE)
  }

  recomputed_stability <- calibrate_by_null(
    scores$pi_raw, null_result$null_frequencies,
    epsilon = 0.01, winsorize_at = 4
  )
  recomputed_utility <- calibrate_by_null(
    observed_raw_utility, null_result$null_utility,
    epsilon = utility_epsilon, winsorize_at = 4
  )
  stability_max_abs_difference <- max(abs(
    recomputed_stability - scores$pi_scored
  ))
  utility_max_abs_difference <- max(abs(
    recomputed_utility - scores$u_scored
  ))

  null_stability_ratios <- leave_one_permutation_out(
    null_result$null_frequencies, epsilon = 0.01
  )
  null_utility_ratios <- leave_one_permutation_out(
    null_result$null_utility, epsilon = utility_epsilon
  )
  null_summary <- rbind(
    summarise_ratio_matrix(
      null_stability_ratios, "stability", "leave_one_permutation_out_null"
    ),
    summarise_ratio_matrix(
      null_utility_ratios, "utility", "leave_one_permutation_out_null"
    ),
    summarise_ratio_matrix(
      scores$pi_scored, "stability", "observed_outcome"
    ),
    summarise_ratio_matrix(
      scores$u_scored, "utility", "observed_outcome"
    )
  )
  qa <- data.frame(
    dataset = dataset,
    repeat_idx = null_repeat,
    fold_idx = null_fold,
    selected_alpha = current$alpha,
    MI_percentile_max_abs_difference = mi_max_abs_difference,
    stability_ratio_max_abs_difference = stability_max_abs_difference,
    utility_ratio_max_abs_difference = utility_max_abs_difference,
    utility_epsilon = utility_epsilon,
    null_n_fits = null_result$diagnostics$n_fits,
    null_n_failed_fits = null_result$diagnostics$n_failed_fits,
    null_failure_fraction = null_result$diagnostics$failure_fraction,
    passed = mi_max_abs_difference <= 1e-12 &
      stability_max_abs_difference <= 1e-12 &
      utility_max_abs_difference <= 1e-12 &
      null_result$diagnostics$failure_fraction <= 0.1,
    stringsAsFactors = FALSE
  )
  if (!qa$passed) {
    stop("Independent calibration-null reconstruction failed.",
         call. = FALSE)
  }
  saveRDS(
    list(
      null_frequencies = null_result$null_frequencies,
      null_utility = null_result$null_utility,
      null_summary = null_summary,
      qa = qa,
      diagnostics = null_result$diagnostics
    ),
    null_checkpoint, version = 3
  )
  cat(sprintf(
    "[%s %s] calibration null reconstructed and saved mapping reproduced\n",
    dataset, null_split_name
  ))
}

if (stage %in% c("assemble", "all")) {
  if (!all(file.exists(all_extract_files))) {
    missing <- sum(!file.exists(all_extract_files))
    stop(sprintf("%d calibration-extraction checkpoints are incomplete.",
                 missing), call. = FALSE)
  }
  if (!file.exists(null_checkpoint)) {
    stop("The fixed-split calibration null check is incomplete.",
         call. = FALSE)
  }
  extracts <- lapply(all_extract_files, readRDS)
  gene_level <- do.call(rbind, lapply(extracts, `[[`, "gene_level"))
  dependence <- do.call(rbind, lapply(extracts, `[[`, "dependence"))
  extract_qa <- do.call(rbind, lapply(extracts, `[[`, "qa"))
  null_check <- readRDS(null_checkpoint)
  qa <- merge(
    extract_qa,
    data.frame(
      null_check_split = null_split_name,
      null_mapping_passed = null_check$qa$passed
    ),
    by = NULL
  )
  gene_keys <- paste(
    gene_level$repeat_idx, gene_level$fold_idx, gene_level$gene
  )
  split_keys <- paste(dependence$repeat_idx, dependence$fold_idx)
  if (nrow(gene_level) != 15L * 2000L || anyDuplicated(gene_keys) ||
      nrow(dependence) != 15L || anyDuplicated(split_keys) ||
      !all(qa$passed) || !all(qa$null_mapping_passed)) {
    stop("Assembled calibration diagnostics failed validation.",
         call. = FALSE)
  }

  metric_names <- c(
    "spearman_calibrated_recurrence_utility_all",
    "spearman_calibrated_recurrence_utility_selected",
    "spearman_raw_recurrence_SHAP_all",
    "spearman_SHAP_MI_all",
    "spearman_calibrated_utility_MI_all",
    "recurrence_utility_jaccard_k10",
    "recurrence_utility_jaccard_k50",
    "recurrence_utility_jaccard_k200",
    "stability_median_ratio", "stability_sd_log2_ratio",
    "utility_median_ratio", "utility_sd_log2_ratio",
    "null_failure_fraction"
  )
  summary_rows <- lapply(metric_names, function(metric) {
    values <- dependence[[metric]]
    finite <- values[is.finite(values)]
    data.frame(
      dataset = dataset,
      metric = metric,
      mean = if (length(finite) > 0L) mean(finite) else NA_real_,
      median = if (length(finite) > 0L) stats::median(finite) else NA_real_,
      minimum = if (length(finite) > 0L) min(finite) else NA_real_,
      maximum = if (length(finite) > 0L) max(finite) else NA_real_,
      n_splits = length(finite),
      stringsAsFactors = FALSE
    )
  })
  summary_table <- do.call(rbind, summary_rows)
  write.csv(gene_level, gene_output, row.names = FALSE)
  write.csv(dependence, dependence_output, row.names = FALSE)
  write.csv(summary_table, summary_output, row.names = FALSE)
  write.csv(null_check$null_summary, null_output, row.names = FALSE)
  write.csv(qa, qa_output, row.names = FALSE)
  cat(sprintf(
    "%s: wrote calibration diagnostics for 15 splits and one fixed null check\n",
    dataset
  ))
}
