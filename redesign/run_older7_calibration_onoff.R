#!/usr/bin/env Rscript

# Calibration on/off comparison at fixed components for the seven older
# datasets. The component ablation (redesign/run_older7_component_ablation.R)
# scored only the CALIBRATED pillars (pi_scored, u_scored, saved ranking); it
# never scored the calibration-OFF counterparts, so the contribution of
# evidence-ratio calibration itself was unmeasured. This script closes that
# gap. Rankings are reconstructed from each saved outer-training fit; no
# selector is refitted. Each ranking is evaluated on the untouched outer test
# fold with the benchmark's deterministic three-model ensemble.
#
# Variants per split (same selected alpha as the ablation, read from the
# ranking metadata):
#   GS_raw_recurrence  rank by pi_raw (raw selection frequency). Tiebreak:
#                      gene name. The ablation's recurrence_only used
#                      (-pi_scored, -pi_raw, gene); on the raw score the
#                      secondary key collapses to gene.
#   GS_raw_SHAPxMI     rank by u. APPROXIMATION: the exact raw product
#                      sqrt(u_instance_raw * u_mi_raw) is not stored in the
#                      fit; u is percentile01(sqrt(percentile01(SHAP instance
#                      frequency) * u_mi + 1e-10)), a percentile-of-percentiles
#                      construction. It is the closest saved calibration-off
#                      utility score. Ranking by u equals ranking by the
#                      saved u column, so results are exact for u and only
#                      approximate for the unsaved raw product.
#   GS_raw_combination rank by novelty_score, the saved calibration-off
#                      combined score percentile01(geomean(pi_final, u)) with
#                      pi_final = percentile01(pi_raw) (gate off, so
#                      pi_filtered == pi_raw).
#   GS_current         the saved ranking (calibration-on combined score).
#                      Included so every checkpoint re-verifies that this
#                      pipeline reproduces eval_deterministic.csv to 1e-12.
#
# Outputs go ONLY to redesign/results_corrected/older7_calibration_onoff_2026-09-03/.
# Nothing is written to the dataset directories.
#
# Run from the repository root:
#   Rscript redesign/run_older7_calibration_onoff.R [stage] [budget_seconds]

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1L) args[[1L]] else "all"
budget <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 604800
stopifnot(
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
validation_datasets <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
)
datasets <- c(validation_datasets, "imvigor210", "sosall")
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
headline_k <- c(10L, 20L, 50L)
variants <- c(
  "GS_raw_recurrence", "GS_raw_SHAPxMI", "GS_raw_combination", "GS_current"
)
# Calibrated counterparts from the component ablation, paired by construction:
# same split, same fit, same evaluator seed.
counterparts <- c(
  GS_raw_recurrence = "GS_recurrence_only",
  GS_raw_SHAPxMI = "GS_SHAPxMI",
  GS_raw_combination = "GS_current"
)

output_dir <- file.path(
  results_root, "older7_calibration_onoff_2026-09-03"
)
checkpoint_dir <- file.path(output_dir, "checkpoints")

split_dir_for <- function(dataset) {
  if (dataset %in% validation_datasets) {
    file.path(results_root, "validation_benchmark", dataset)
  } else {
    file.path(results_root, "grouped_benchmark", dataset)
  }
}
dataset_dir_for <- function(dataset) {
  if (dataset %in% validation_datasets) {
    split_dir_for(dataset)
  } else {
    file.path(results_root, "full_recipe", dataset)
  }
}
saved_method_for <- function(dataset) {
  if (dataset %in% validation_datasets) "GS_full_ungrouped" else "full_ungrouped"
}

split_path <- function(dataset, repeat_idx, fold_idx) file.path(
  split_dir_for(dataset), sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
)
ranking_path <- function(dataset, repeat_idx, fold_idx) file.path(
  dataset_dir_for(dataset), sprintf(
    "ranking_r%d_f%d_%s.csv", repeat_idx, fold_idx, saved_method_for(dataset)
  )
)
meta_path <- function(dataset, repeat_idx, fold_idx) file.path(
  dataset_dir_for(dataset), sprintf(
    "ranking_r%d_f%d_%s_meta.csv", repeat_idx, fold_idx,
    saved_method_for(dataset)
  )
)
# Alpha selection mirrors the ablation exactly: the ranking metadata records
# which of the two saved alpha fits (0.5, 1.0) the benchmark selected on
# inner performance, and the same selected fit is reused here.
read_selected_alpha <- function(dataset, repeat_idx, fold_idx) {
  metadata <- read.csv(
    meta_path(dataset, repeat_idx, fold_idx), stringsAsFactors = FALSE
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
fit_path <- function(dataset, repeat_idx, fold_idx) {
  alpha <- read_selected_alpha(dataset, repeat_idx, fold_idx)
  file.path(dataset_dir_for(dataset), sprintf(
    "fit_r%d_f%d_%s_a%s.rds",
    repeat_idx, fold_idx, saved_method_for(dataset), alpha_tag(alpha)
  ))
}
outcome_path_for <- function(dataset) file.path(
  dataset_dir_for(dataset),
  if (dataset %in% validation_datasets) "base_data.rds" else "base_outcome.rds"
)
evaluation_path_for <- function(dataset) file.path(
  dataset_dir_for(dataset), "eval_deterministic.csv"
)
ablation_path_for <- function(dataset) file.path(
  dataset_dir_for(dataset), "older7_component_ablation_evaluation.csv"
)
random30_path_for <- function(dataset) file.path(
  dataset_dir_for(dataset), "older7_random_baseline_30_draws.csv"
)
checkpoint_path <- function(dataset, repeat_idx, fold_idx) file.path(
  checkpoint_dir, sprintf(
    "older7_calibration_onoff_%s_r%d_f%d.rds", dataset, repeat_idx, fold_idx
  )
)

split_grid <- expand.grid(
  dataset = datasets, repeat_idx = 1:3, fold_idx = 1:5,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
split_override <- Sys.getenv("GENESELECTR_VALIDATION_SPLIT", "")
if (nzchar(split_override)) {
  parsed <- regmatches(
    split_override, regexec("^([A-Za-z0-9]+)_r([1-3])f([1-5])$", split_override)
  )[[1L]]
  if (length(parsed) != 4L || !parsed[[2L]] %in% datasets) {
    stop("GENESELECTR_VALIDATION_SPLIT must have the form <dataset>_r1f1.",
         call. = FALSE)
  }
  split_grid <- data.frame(
    dataset = parsed[[2L]],
    repeat_idx = as.integer(parsed[[3L]]),
    fold_idx = as.integer(parsed[[4L]]),
    stringsAsFactors = FALSE
  )
}

# --- Inputs and provenance ----------------------------------------------------
for (dataset in datasets) {
  require_redesign_run(split_dir_for(dataset))
  require_redesign_run(dataset_dir_for(dataset))
}
input_files <- unlist(lapply(datasets, function(dataset) {
  per_split <- unlist(lapply(1:3, function(repeat_idx) {
    unlist(lapply(1:5, function(fold_idx) {
      c(
        ranking_path(dataset, repeat_idx, fold_idx),
        meta_path(dataset, repeat_idx, fold_idx),
        fit_path(dataset, repeat_idx, fold_idx)
      )
    }))
  }))
  c(
    file.path(split_dir_for(dataset), "run_manifest.rds"),
    file.path(dataset_dir_for(dataset), "run_manifest.rds"),
    per_split,
    evaluation_path_for(dataset),
    ablation_path_for(dataset),
    random30_path_for(dataset),
    outcome_path_for(dataset)
  )
}), use.names = FALSE)
if (!all(file.exists(input_files))) {
  missing <- input_files[!file.exists(input_files)]
  stop(sprintf(
    "Required inputs are missing (%d), first: %s",
    length(missing), missing[[1L]]
  ), call. = FALSE)
}

by_split_output <- file.path(output_dir, "older7_calibration_onoff_by_split.csv")
summary_output <- file.path(output_dir, "older7_calibration_onoff_summary.csv")
qa_output <- file.path(output_dir, "older7_calibration_onoff_QA.csv")
panels_output <- file.path(output_dir, "older7_calibration_onoff_panels.csv")
report_output <- file.path(output_dir, "report.md")
all_checkpoint_files <- unlist(Map(
  checkpoint_path, split_grid$dataset, split_grid$repeat_idx,
  split_grid$fold_idx
), use.names = FALSE)

source_files <- c(
  "redesign/run_older7_calibration_onoff.R",
  file.path("package", "GeneSelectR", "R", "utils.R"),
  file.path("redesign", "R", c("bio_prior.R", "evaluator.R", "run_provenance.R"))
)
config <- list(
  datasets = datasets,
  variants = variants,
  calibrated_counterparts = counterparts,
  panel_sizes = panel_sizes,
  headline_k = headline_k,
  fit_source = "saved_outer_training_fit",
  raw_recurrence_score = "pi_raw",
  raw_SHAPxMI_score = paste0(
    "u (approximation: percentile-of-percentiles utility; the exact raw ",
    "product sqrt(u_instance_raw * u_mi_raw) is not saved)"
  ),
  raw_combination_score = "novelty_score (saved calibration-off geomean)",
  random_baseline = "mean of 30 matched Random draws per dataset/split/k",
  evaluator = "glmnet+xgboost+ranger",
  model_seed = "420000 + 1000 * repeat + 100 * fold + panel_index"
)

# A fresh output dir gets its own run manifest; the extension manifest then
# pins this script, its sources, and every input file it reads. Neither
# touches the dataset directories. Both calls are idempotent: on rerun they
# recompute the expected manifest and stop if code or inputs changed.
prepare_redesign_run(
  output_dir, config = config, source_files = source_files,
  input_files = input_files
)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
prepare_redesign_extension(
  output_dir,
  extension = "older7_calibration_onoff_v1",
  source_files = source_files,
  config = config,
  output_files = c(
    all_checkpoint_files, by_split_output, summary_output, qa_output,
    panels_output, report_output
  )
)
if (stage == "init") {
  cat("older-seven calibration on/off initialized\n")
  quit(save = "no")
}

make_rankings <- function(fit, saved_ranking) {
  scores <- fit$gene_scores
  required <- c(
    "gene", "pi_raw", "u", "u_mi", "pi_scored", "u_scored", "novelty_score"
  )
  if (!all(required %in% names(scores)) || anyDuplicated(scores$gene) ||
      !identical(as.character(scores$gene),
                 as.character(saved_ranking$gene))) {
    stop("Saved fit and ranking are inconsistent.", call. = FALSE)
  }
  numeric_columns <- setdiff(required, "gene")
  if (any(!is.finite(as.matrix(scores[numeric_columns])))) {
    stop("GeneSelectR component scores contain non-finite values.",
         call. = FALSE)
  }
  list(
    GS_raw_recurrence = scores$gene[order(
      -scores$pi_raw, scores$gene, method = "radix"
    )],
    GS_raw_SHAPxMI = scores$gene[order(
      -scores$u, scores$gene, method = "radix"
    )],
    GS_raw_combination = scores$gene[order(
      -scores$novelty_score, scores$gene, method = "radix"
    )],
    GS_current = as.character(saved_ranking$gene)
  )
}

if (stage %in% c("evaluate", "all")) {
  for (job_idx in seq_len(nrow(split_grid))) {
    dataset <- split_grid$dataset[[job_idx]]
    repeat_idx <- split_grid$repeat_idx[[job_idx]]
    fold_idx <- split_grid$fold_idx[[job_idx]]
    output_path <- checkpoint_path(dataset, repeat_idx, fold_idx)
    if (file.exists(output_path)) next
    if (over_budget()) {
      cat("[budget] stop in calibration on/off evaluation\n")
      quit(save = "no")
    }

    outcome <- if (dataset %in% validation_datasets) {
      readRDS(outcome_path_for(dataset))$outcome
    } else {
      readRDS(outcome_path_for(dataset))
    }
    outcome <- droplevels(as.factor(outcome))
    split <- readRDS(split_path(dataset, repeat_idx, fold_idx))
    pool <- split$pools$var2000
    standardized <- standardise_split(
      split$train_raw[, pool, drop = FALSE],
      split$test_raw[, pool, drop = FALSE]
    )
    y_train <- droplevels(outcome[split$train_idx])
    y_test <- droplevels(outcome[split$test_idx])
    wrapper <- readRDS(fit_path(dataset, repeat_idx, fold_idx))
    fit <- wrapper$fit
    selected_alpha <- read_selected_alpha(dataset, repeat_idx, fold_idx)
    if (is.null(fit) || !isTRUE(all.equal(wrapper$alpha, selected_alpha)) ||
        !identical(fit$parameters$components, c("stability", "utility")) ||
        !identical(fit$parameters$gate_method, "none") ||
        !identical(fit$parameters$calibration_mode, "evidence_ratio") ||
        !identical(fit$parameters$utility_method, "instance_shap")) {
      stop("Saved fit is not the current benchmark configuration.",
           call. = FALSE)
    }
    saved_ranking <- read.csv(
      ranking_path(dataset, repeat_idx, fold_idx), stringsAsFactors = FALSE
    )
    rankings <- make_rankings(fit, saved_ranking)
    if (!identical(names(rankings), variants) ||
        any(vapply(rankings, anyDuplicated, integer(1)) > 0L) ||
        any(!vapply(rankings, function(genes) {
          setequal(genes, pool)
        }, logical(1)))) {
      stop("A calibration on/off ranking does not match the candidate pool.",
           call. = FALSE)
    }

    saved_evaluation <- read.csv(
      evaluation_path_for(dataset), stringsAsFactors = FALSE
    )
    saved_current <- saved_evaluation[
      saved_evaluation$repeat_idx == repeat_idx &
        saved_evaluation$fold_idx == fold_idx &
        saved_evaluation$arm == "GS_full_ungrouped",
      c("k", "AUC"), drop = FALSE
    ]
    if (nrow(saved_current) != length(panel_sizes)) {
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
        evaluation_rows[[length(evaluation_rows) + 1L]] <- data.frame(
          dataset = dataset,
          repeat_idx = repeat_idx,
          fold_idx = fold_idx,
          variant = variant,
          k = panel_size,
          AUC = bench_auc(y_test, predictions),
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

    # Sanity check, per checkpoint: the recomputed GS_current AUCs must
    # reproduce eval_deterministic.csv to 1e-12. This validates the whole
    # pipeline (split handling, standardisation, panel extraction, evaluator,
    # seed) before any raw-variant number is trusted.
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
    keys <- paste(evaluation$variant, evaluation$k)
    if (nrow(evaluation) != length(variants) * length(panel_sizes) ||
        anyDuplicated(keys) || any(!is.finite(evaluation$AUC)) ||
        any(evaluation$AUC < 0 | evaluation$AUC > 1) || !all(qa$passed)) {
      stop("Calibration on/off checkpoint failed validation.", call. = FALSE)
    }
    saveRDS(
      list(evaluation = evaluation, panels = panels, qa = qa),
      output_path, version = 3
    )
    cat(sprintf(
      "[%s r%d f%d] calibration on/off complete (max |dAUC| vs saved %.3g)\n",
      dataset, repeat_idx, fold_idx, max(qa$absolute_difference)
    ))
  }
}

if (stage %in% c("assemble", "all")) {
  if (!all(file.exists(all_checkpoint_files))) {
    missing <- sum(!file.exists(all_checkpoint_files))
    stop(sprintf("%d calibration on/off checkpoints are incomplete.", missing),
         call. = FALSE)
  }
  checkpoints <- lapply(all_checkpoint_files, readRDS)
  evaluation <- do.call(rbind, lapply(checkpoints, `[[`, "evaluation"))
  panels <- do.call(rbind, lapply(checkpoints, `[[`, "panels"))
  qa <- do.call(rbind, lapply(checkpoints, `[[`, "qa"))
  evaluation_keys <- paste(
    evaluation$dataset, evaluation$repeat_idx, evaluation$fold_idx,
    evaluation$variant, evaluation$k
  )
  expected_rows <- nrow(split_grid) * length(variants) * length(panel_sizes)
  if (nrow(evaluation) != expected_rows || anyDuplicated(evaluation_keys) ||
      !all(qa$passed)) {
    stop("Assembled calibration on/off failed validation.", call. = FALSE)
  }

  # 30-draw Random baseline, matched per dataset/repeat/fold/k (same loading
  # as summarise_integrated_interpretation.R).
  random_draws <- do.call(rbind, lapply(datasets, function(dataset) {
    value <- read.csv(random30_path_for(dataset), stringsAsFactors = FALSE)
    value$dataset <- dataset
    value
  }))
  random_mean <- stats::aggregate(
    AUC ~ dataset + repeat_idx + fold_idx + k, data = random_draws, FUN = mean
  )
  names(random_mean)[names(random_mean) == "AUC"] <- "random_AUC_30draw"
  if (any(vapply(split(random_mean, random_mean$dataset), nrow,
                 integer(1)) != 15L * length(panel_sizes))) {
    stop("30-draw Random baseline is incomplete.", call. = FALSE)
  }
  split_keys <- c("dataset", "repeat_idx", "fold_idx", "k")
  by_split <- merge(
    evaluation, random_mean, by = split_keys, all.x = TRUE, sort = FALSE
  )
  by_split$delta_random_30draw <- by_split$AUC - by_split$random_AUC_30draw
  if (any(is.na(by_split$delta_random_30draw))) {
    stop("30-draw Random baseline does not cover every evaluation cell.",
         call. = FALSE)
  }

  # Calibrated counterparts from the component ablation, re-based on the same
  # 30-draw baseline so the pairing is exact.
  ablation <- do.call(rbind, lapply(datasets, function(dataset) {
    value <- read.csv(ablation_path_for(dataset), stringsAsFactors = FALSE)
    value[value$variant %in% counterparts,
          c("dataset", "repeat_idx", "fold_idx", "variant", "k", "AUC")]
  }))
  names(ablation)[names(ablation) == "variant"] <- "calibrated_variant"
  names(ablation)[names(ablation) == "AUC"] <- "calibrated_AUC"
  raw_long <- by_split[by_split$variant %in% names(counterparts), ]
  raw_long$calibrated_variant <- counterparts[raw_long$variant]
  paired <- merge(
    raw_long, ablation,
    by = c("dataset", "repeat_idx", "fold_idx", "k", "calibrated_variant"),
    all.x = TRUE, sort = FALSE
  )
  if (nrow(paired) != nrow(raw_long) || any(is.na(paired$calibrated_AUC))) {
    stop("Calibrated ablation counterparts do not cover every raw cell.",
         call. = FALSE)
  }
  # Same split, same evaluator seed, same Random baseline: the difference in
  # AUC equals the difference in delta_random_30draw.
  paired$raw_minus_calibrated <- paired$AUC - paired$calibrated_AUC

  # Summary: paired same-split differences, per dataset x pair x k, then
  # dataset-balanced; plus the primary endpoint (mean over k in {10,20,50}
  # per dataset, then dataset-balanced).
  pair_summary <- do.call(rbind, lapply(split(
    paired, interaction(paired$dataset, paired$variant, paired$k, drop = TRUE)
  ), function(value) {
    data.frame(
      dataset = value$dataset[1L],
      raw_variant = value$variant[1L],
      calibrated_variant = value$calibrated_variant[1L],
      k = value$k[1L],
      mean_raw_AUC = mean(value$AUC),
      mean_calibrated_AUC = mean(value$calibrated_AUC),
      mean_raw_delta_30draw = mean(value$delta_random_30draw),
      mean_raw_minus_calibrated = mean(value$raw_minus_calibrated),
      sd_raw_minus_calibrated = stats::sd(value$raw_minus_calibrated),
      n_splits = nrow(value),
      stringsAsFactors = FALSE
    )
  }))
  rownames(pair_summary) <- NULL
  balanced_summary <- do.call(rbind, lapply(split(
    pair_summary, interaction(pair_summary$raw_variant, pair_summary$k,
                              drop = TRUE)
  ), function(value) {
    data.frame(
      dataset = "dataset_balanced",
      raw_variant = value$raw_variant[1L],
      calibrated_variant = value$calibrated_variant[1L],
      k = value$k[1L],
      mean_raw_AUC = mean(value$mean_raw_AUC),
      mean_calibrated_AUC = mean(value$mean_calibrated_AUC),
      mean_raw_delta_30draw = mean(value$mean_raw_delta_30draw),
      mean_raw_minus_calibrated = mean(value$mean_raw_minus_calibrated),
      sd_raw_minus_calibrated = stats::sd(value$mean_raw_minus_calibrated),
      n_splits = nrow(value),
      stringsAsFactors = FALSE
    )
  }))
  rownames(balanced_summary) <- NULL
  # Primary endpoint: average the paired difference over k in {10,20,50}
  # WITHIN each split first, so the sd is a genuine split-level sd. Per
  # dataset: mean/sd over the 15 splits. Dataset-balanced: mean of the seven
  # per-dataset means (equal dataset weight), with the pooled split-level sd
  # over all 105 splits reported alongside.
  paired_headline <- paired[paired$k %in% headline_k, ]
  primary_split <- do.call(rbind, lapply(split(
    paired_headline,
    interaction(paired_headline$dataset, paired_headline$variant,
                paired_headline$repeat_idx, paired_headline$fold_idx,
                drop = TRUE)
  ), function(value) {
    if (nrow(value) != length(headline_k)) {
      stop("Primary endpoint does not cover every headline k.", call. = FALSE)
    }
    data.frame(
      dataset = value$dataset[1L],
      raw_variant = value$variant[1L],
      calibrated_variant = value$calibrated_variant[1L],
      repeat_idx = value$repeat_idx[1L],
      fold_idx = value$fold_idx[1L],
      primary_raw_AUC = mean(value$AUC),
      primary_calibrated_AUC = mean(value$calibrated_AUC),
      primary_raw_delta_30draw = mean(value$delta_random_30draw),
      primary_raw_minus_calibrated = mean(value$raw_minus_calibrated),
      stringsAsFactors = FALSE
    )
  }))
  rownames(primary_split) <- NULL
  primary <- do.call(rbind, lapply(split(
    primary_split,
    interaction(primary_split$dataset, primary_split$raw_variant, drop = TRUE)
  ), function(value) {
    data.frame(
      dataset = value$dataset[1L],
      raw_variant = value$raw_variant[1L],
      calibrated_variant = value$calibrated_variant[1L],
      k = "primary_k10_k20_k50",
      mean_raw_AUC = mean(value$primary_raw_AUC),
      mean_calibrated_AUC = mean(value$primary_calibrated_AUC),
      mean_raw_delta_30draw = mean(value$primary_raw_delta_30draw),
      mean_raw_minus_calibrated = mean(value$primary_raw_minus_calibrated),
      sd_raw_minus_calibrated = stats::sd(value$primary_raw_minus_calibrated),
      n_splits = nrow(value),
      stringsAsFactors = FALSE
    )
  }))
  rownames(primary) <- NULL
  primary_balanced <- do.call(rbind, lapply(split(
    primary_split, primary_split$raw_variant
  ), function(value) {
    per_dataset <- primary[primary$raw_variant == value$raw_variant[1L], ]
    data.frame(
      dataset = "dataset_balanced",
      raw_variant = value$raw_variant[1L],
      calibrated_variant = value$calibrated_variant[1L],
      k = "primary_k10_k20_k50",
      mean_raw_AUC = mean(per_dataset$mean_raw_AUC),
      mean_calibrated_AUC = mean(per_dataset$mean_calibrated_AUC),
      mean_raw_delta_30draw = mean(per_dataset$mean_raw_delta_30draw),
      mean_raw_minus_calibrated = mean(per_dataset$mean_raw_minus_calibrated),
      sd_raw_minus_calibrated = stats::sd(value$primary_raw_minus_calibrated),
      n_splits = nrow(value),
      stringsAsFactors = FALSE
    )
  }))
  rownames(primary_balanced) <- NULL
  summary_table <- rbind(pair_summary, balanced_summary, primary,
                         primary_balanced)
  summary_table <- summary_table[order(
    summary_table$raw_variant, summary_table$dataset,
    match(summary_table$k, c(as.character(panel_sizes), "primary_k10_k20_k50"))
  ), ]

  write.csv(by_split, by_split_output, row.names = FALSE)
  write.csv(panels, panels_output, row.names = FALSE)
  write.csv(summary_table, summary_output, row.names = FALSE)
  write.csv(qa, qa_output, row.names = FALSE)

  # Report: verdict from the primary-endpoint dataset-balanced differences.
  # Threshold for "neither": |mean difference| <= 0.005 AUC, smaller than the
  # split-level Monte Carlo noise of the evaluation itself.
  verdict_lines <- vapply(seq_len(nrow(primary_balanced)), function(row_idx) {
    value <- primary_balanced[row_idx, ]
    direction <- if (value$mean_raw_minus_calibrated > 0.005) {
      "calibration HURTS (raw ranks better)"
    } else if (value$mean_raw_minus_calibrated < -0.005) {
      "calibration HELPS (raw ranks worse)"
    } else {
      "NEITHER (difference within +/-0.005 AUC)"
    }
    sprintf(
      paste0(
        "- %s vs %s: raw-minus-calibrated = %+.4f AUC ",
        "(split-level sd %.4f, %d splits pooled over k=10/20/50) -> %s. ",
        "Raw mean delta vs 30-draw Random: %+.4f."
      ),
      value$raw_variant, value$calibrated_variant,
      value$mean_raw_minus_calibrated, value$sd_raw_minus_calibrated,
      value$n_splits, direction, value$mean_raw_delta_30draw
    )
  }, character(1))
  report <- c(
    "# Calibration on/off at fixed components — older seven datasets",
    "",
    sprintf("Assembled: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Question: does evidence-ratio calibration of the stability and utility",
    "pillars change predictive performance when the components are held fixed?",
    "Rankings are reconstructed from the saved outer-training fits; no selector",
    "is refitted. Deltas are vs the mean of 30 matched Random draws per",
    "dataset/split/k.",
    "",
    "Approximation note: GS_raw_SHAPxMI ranks by the saved `u`, which is",
    "percentile01(sqrt(percentile01(SHAP instance frequency) * u_mi + 1e-10)).",
    "The exact raw product sqrt(u_instance_raw * u_mi_raw) is not stored in",
    "the fits, so `u` is the closest saved calibration-off utility score.",
    "GS_raw_combination ranks by the saved `novelty_score` (calibration-off",
    "geomean of percentile-scored stability and `u`).",
    "",
    sprintf(
      paste0(
        "Sanity check: the recomputed current ranking reproduced the saved ",
        "eval_deterministic.csv AUCs to <=1e-12 on every split (%d of %d ",
        "dataset/split/k cells passed; max absolute difference %.3g)."
      ),
      sum(qa$passed), nrow(qa), max(qa$absolute_difference)
    ),
    "",
    "## Verdict (primary endpoint: mean over k = 10, 20, 50, dataset-balanced)",
    "",
    verdict_lines,
    "",
    "Per-dataset and per-k numbers are in",
    "`older7_calibration_onoff_summary.csv`; by-split AUCs and deltas in",
    "`older7_calibration_onoff_by_split.csv`."
  )
  writeLines(report, report_output)
  cat(sprintf(
    paste0(
      "wrote %d by-split rows (%d QA cells, max |dAUC| %.3g) and the paired ",
      "summary to %s\n"
    ),
    nrow(by_split), nrow(qa), max(qa$absolute_difference), output_dir
  ))
  cat(paste(verdict_lines, collapse = "\n"), "\n")
}
