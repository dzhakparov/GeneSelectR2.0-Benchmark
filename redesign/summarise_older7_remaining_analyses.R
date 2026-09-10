#!/usr/bin/env Rscript

# Combine the three remaining benchmark analyses (component ablation, 30-draw
# Random baseline sensitivity, calibration diagnostics) across the seven
# datasets. Run from the repository root:
#
#   Rscript redesign/summarise_older7_remaining_analyses.R
#
# All outputs are descriptive. The primary predictive metric is AUC minus the
# MEAN of the 30 Random draws matched on dataset/repeat/fold/k, because the
# earlier 3-draw Random reference is a noisy point estimate of what a random
# panel of the same size achieves (random signatures inherit dominant
# meta-genes, so raw AUC alone is misleading). Nothing here is a significance
# claim: repeated cross-validation splits share observations, so the Wilcoxon
# p-values in combined_ablation_vs_current.csv are indicative only and no
# multiple-testing correction is applied anywhere.

datasets <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683",
  "imvigor210", "sosall"
)
# The five GEO validation datasets live under validation_benchmark/; the two
# older benchmark datasets were run under full_recipe/. Same files, different
# parent directory, hence the helper.
validation_datasets <- datasets[1:5]
variants <- c(
  "GS_recurrence_only", "GS_SHAP_only", "GS_MI_only",
  "GS_SHAPxMI", "GS_current"
)
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
headline_k <- c(10L, 20L, 50L)

results_root <- file.path("redesign", "results_corrected")
output_dir <- file.path(
  results_root, "older7_remaining_analyses_2026-09-02"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_dir <- function(dataset) {
  if (dataset %in% validation_datasets) {
    file.path(results_root, "validation_benchmark", dataset)
  } else {
    file.path(results_root, "full_recipe", dataset)
  }
}

# Every file a complete dataset run must produce. Expected row counts are the
# contract of the three benchmark scripts; a mismatch means the run truncated
# silently and should be re-run rather than summarised.
expected_files <- c(
  ablation_evaluation = "older7_component_ablation_evaluation.csv",
  ablation_qa         = "older7_component_ablation_QA.csv",
  random_draws        = "older7_random_baseline_30_draws.csv",
  random_by_split     = "older7_random_baseline_30_by_split.csv",
  random_summary      = "older7_random_baseline_30_summary.csv",
  random_qa           = "older7_random_baseline_30_QA.csv",
  calibration_split   = "older7_calibration_dependence_by_split.csv",
  calibration_null    = "older7_calibration_null_validation.csv",
  calibration_qa      = "older7_calibration_diagnostics_QA.csv"
)
expected_rows <- c(
  ablation_evaluation = 15L * length(variants) * length(panel_sizes),  # 450
  random_draws        = 15L * length(panel_sizes) * 30L,               # 2700
  random_by_split     = 15L * length(panel_sizes),                     # 90
  random_summary      = length(panel_sizes),                           # 6
  calibration_split   = 15L,
  calibration_null    = 4L
)

# -----------------------------------------------------------------------------
# 1. Run status first: which files exist per dataset, and did QA pass.
# Written before anything else so an incomplete run is diagnosed from the
# output directory even though the script then stops.
# -----------------------------------------------------------------------------
qa_passed <- function(path, dataset, analysis) {
  if (!file.exists(path)) return(NA)
  value <- read.csv(path, stringsAsFactors = FALSE)
  if (!"passed" %in% names(value)) {
    # A QA file without a passed column is malformed; treat as failed rather
    # than silently assuming success.
    warning("QA file for ", dataset, " / ", analysis,
            " has no 'passed' column: ", path, call. = FALSE)
    return(NA)
  }
  isTRUE(all(value$passed))
}

run_status <- do.call(rbind, lapply(datasets, function(dataset) {
  # file.path() does not keep the names of expected_files, so re-attach them;
  # the qa_passed() calls below index paths by name.
  paths <- stats::setNames(
    file.path(dataset_dir(dataset), expected_files), names(expected_files)
  )
  present <- file.exists(paths)
  missing <- names(expected_files)[!present]
  data.frame(
    dataset = dataset,
    result_dir = dataset_dir(dataset),
    n_files_expected = length(expected_files),
    n_files_present = sum(present),
    missing_files = paste(missing, collapse = "; "),
    ablation_qa_passed = qa_passed(
      paths[["ablation_qa"]], dataset, "ablation"
    ),
    calibration_qa_passed = qa_passed(
      paths[["calibration_qa"]], dataset, "calibration"
    ),
    random30_qa_passed = qa_passed(
      paths[["random_qa"]], dataset, "random30"
    ),
    stringsAsFactors = FALSE
  )
}))
run_status$all_files_present <-
  run_status$n_files_present == run_status$n_files_expected
# Row-wise AND: a QA column is NA when its file was absent, and NA must not
# count as passed.
run_status$all_qa_passed <- mapply(
  function(present, a, c, r) {
    present && !is.na(a) && a && !is.na(c) && c && !is.na(r) && r
  },
  run_status$all_files_present,
  run_status$ablation_qa_passed,
  run_status$calibration_qa_passed,
  run_status$random30_qa_passed
)
write.csv(
  run_status, file.path(output_dir, "combined_run_status.csv"),
  row.names = FALSE
)
cat("Wrote", file.path(output_dir, "combined_run_status.csv"), "\n")

missing_report <- run_status[!run_status$all_files_present, ]
if (nrow(missing_report) > 0L) {
  details <- paste(
    sprintf("  %s (%s): missing %s",
            missing_report$dataset,
            missing_report$result_dir,
            missing_report$missing_files),
    collapse = "\n"
  )
  stop(
    "Incomplete benchmark runs detected; combined_run_status.csv was written ",
    "to ", output_dir, " before stopping. Missing files:\n", details,
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. Read everything. All datasets passed the existence check above, so a
# read failure here means a malformed file, which should be a loud error.
# -----------------------------------------------------------------------------
read_one <- function(dataset, key) {
  path <- file.path(dataset_dir(dataset), expected_files[[key]])
  value <- read.csv(path, stringsAsFactors = FALSE)
  # Overwrite rather than trust the stored dataset column: a file copied into
  # the wrong directory would otherwise contaminate another dataset silently.
  value$dataset <- dataset
  value
}
bind_all <- function(key) {
  do.call(rbind, lapply(datasets, read_one, key = key))
}

ablation <- bind_all("ablation_evaluation")
random_draws <- bind_all("random_draws")
random_by_split <- bind_all("random_by_split")
random_summary <- bind_all("random_summary")
calibration_split <- bind_all("calibration_split")
calibration_null <- bind_all("calibration_null")

# Row-count sanity checks: catches truncated runs that still wrote a header.
check_rows <- function(value, key) {
  per_dataset <- table(value$dataset)
  for (dataset in datasets) {
    if (is.na(per_dataset[[dataset]]) ||
        per_dataset[[dataset]] != expected_rows[[key]]) {
      stop(
        "Row-count check failed for ", dataset, " in ", expected_files[[key]],
        ": expected ", expected_rows[[key]], " rows, got ",
        ifelse(is.na(per_dataset[[dataset]]), 0L, per_dataset[[dataset]]),
        ". Re-run the benchmark for this dataset.",
        call. = FALSE
      )
    }
  }
}
check_rows(ablation, "ablation_evaluation")
check_rows(random_draws, "random_draws")
check_rows(calibration_split, "calibration_split")

# -----------------------------------------------------------------------------
# 3. Merge ablation AUCs with the MEAN of the 30 Random draws matched on
# dataset/repeat/fold/k. The mean is the reference because a single draw (or
# three) of random panels has large Monte Carlo error; the by-split files
# quantify that error explicitly.
# -----------------------------------------------------------------------------
random_mean <- stats::aggregate(
  AUC ~ dataset + repeat_idx + fold_idx + k,
  data = random_draws, FUN = mean
)
names(random_mean)[names(random_mean) == "AUC"] <- "mean_random_AUC_30draw"

split_keys <- c("dataset", "repeat_idx", "fold_idx", "k")
ablation_merged <- merge(
  ablation, random_mean, by = split_keys, all.x = TRUE, sort = FALSE
)
if (any(is.na(ablation_merged$mean_random_AUC_30draw))) {
  stop(
    "Merged ablation table has NA in mean_random_AUC_30draw: the 30-draw ",
    "Random baseline does not cover every dataset/repeat/fold/k cell of the ",
    "ablation evaluation. Check both runs for the affected dataset.",
    call. = FALSE
  )
}
ablation_merged$delta_random_30draw <-
  ablation_merged$AUC - ablation_merged$mean_random_AUC_30draw
if (any(is.na(ablation_merged$delta_random_30draw))) {
  stop(
    "NA in delta_random_30draw after merge; AUC or mean_random_AUC_30draw ",
    "contains missing values.",
    call. = FALSE
  )
}
write.csv(
  ablation_merged,
  file.path(output_dir, "combined_ablation_vs_random30_by_split.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Dataset x variant x k summary of AUC and delta vs the 30-draw Random
# mean, over the 15 splits.
# -----------------------------------------------------------------------------
summarise_delta <- function(value) {
  data.frame(
    mean_AUC = mean(value$AUC),
    sd_AUC = stats::sd(value$AUC),
    mean_delta_random_30draw = mean(value$delta_random_30draw),
    sd_delta_random_30draw = stats::sd(value$delta_random_30draw),
    n_splits = nrow(value),
    stringsAsFactors = FALSE
  )
}
group_summary <- function(data, group_columns, fn) {
  groups <- split(
    data, interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  )
  output <- do.call(rbind, lapply(groups, function(value) {
    cbind(value[1L, group_columns, drop = FALSE], fn(value))
  }))
  rownames(output) <- NULL
  output
}
ablation_summary <- group_summary(
  ablation_merged, c("dataset", "variant", "k"), summarise_delta
)
write.csv(
  ablation_summary,
  file.path(output_dir, "combined_ablation_vs_random30_summary.csv"),
  row.names = FALSE
)
headline_summary <- ablation_summary[ablation_summary$k %in% headline_k, ]
write.csv(
  headline_summary,
  file.path(output_dir, "combined_dataset_level_k10_k20_k50.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. Every ablation variant vs GS_current, per dataset x k. Paired on
# (repeat_idx, fold_idx) because both variants saw the same split. The
# Wilcoxon signed-rank p-value is descriptive: repeated CV splits share
# observations, so the splits are not independent and the p-value cannot be
# read as a significance test.
# -----------------------------------------------------------------------------
note_descriptive <- paste(
  "Descriptive only: repeated cross-validation splits share observations,",
  "so the paired Wilcoxon p-value is indicative, not a significance test.",
  "No multiple-testing correction applied."
)
paired_wilcox <- function(difference) {
  tryCatch(
    stats::wilcox.test(difference, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}
grid <- expand.grid(
  dataset = datasets, k = panel_sizes, variant = variants,
  stringsAsFactors = FALSE
)
grid <- grid[grid$variant != "GS_current", ]
vs_current <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  dataset <- grid$dataset[i]
  k <- grid$k[i]
  variant <- grid$variant[i]
  subset_key <- ablation_merged$dataset == dataset & ablation_merged$k == k
  cell <- ablation_merged[subset_key, ]
  current <- cell[cell$variant == "GS_current",
                  c("repeat_idx", "fold_idx", "delta_random_30draw")]
  other <- cell[cell$variant == variant,
                c("repeat_idx", "fold_idx", "delta_random_30draw")]
  paired <- merge(
    other, current, by = c("repeat_idx", "fold_idx"),
    suffixes = c("_variant", "_current"), sort = FALSE
  )
  difference <-
    paired$delta_random_30draw_variant - paired$delta_random_30draw_current
  data.frame(
    dataset = dataset,
    k = k,
    variant = variant,
    mean_delta_variant = mean(paired$delta_random_30draw_variant),
    mean_delta_current = mean(paired$delta_random_30draw_current),
    mean_difference_vs_current = mean(difference),
    sd_difference_vs_current = stats::sd(difference),
    n_paired_splits = nrow(paired),
    wilcoxon_p_descriptive = paired_wilcox(difference),
    note = note_descriptive,
    stringsAsFactors = FALSE
  )
}))
write.csv(
  vs_current, file.path(output_dir, "combined_ablation_vs_current.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Component contributions: mean delta_random_30draw of each variant side
# by side at k = 10/20/50, so the driver (recurrence, SHAP, MI, or the
# combination) is visible per dataset and panel size.
# -----------------------------------------------------------------------------
headline_merged <- ablation_merged[ablation_merged$k %in% headline_k, ]
contributions <- do.call(rbind, lapply(datasets, function(dataset) {
  do.call(rbind, lapply(headline_k, function(k) {
    cell <- headline_merged[
      headline_merged$dataset == dataset & headline_merged$k == k, ]
    row <- data.frame(
      dataset = dataset, k = k, stringsAsFactors = FALSE
    )
    for (variant in variants) {
      row[[paste0("mean_delta_", variant)]] <-
        mean(cell$delta_random_30draw[cell$variant == variant])
    }
    row
  }))
}))
write.csv(
  contributions,
  file.path(output_dir, "combined_component_contributions.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 7. Calibration behaviour. By-split table is a straight bind; the summary
# reports mean/median/min/max across the 15 splits for each diagnostic so a
# single pathological split is visible in min/max rather than hidden in a
# mean.
# -----------------------------------------------------------------------------
write.csv(
  calibration_split,
  file.path(output_dir, "combined_calibration_dependence_by_split.csv"),
  row.names = FALSE
)
calibration_metrics <- c(
  "spearman_calibrated_recurrence_utility_all",
  "spearman_calibrated_recurrence_utility_selected",
  "spearman_raw_recurrence_SHAP_all",
  "spearman_SHAP_MI_all",
  "spearman_calibrated_utility_MI_all",
  "recurrence_utility_jaccard_k10",
  "recurrence_utility_jaccard_k50",
  "recurrence_utility_jaccard_k200",
  "stability_median_ratio",
  "stability_sd_log2_ratio",
  "utility_median_ratio",
  "utility_sd_log2_ratio",
  "null_failure_fraction"
)
missing_metrics <- setdiff(calibration_metrics, names(calibration_split))
if (length(missing_metrics) > 0L) {
  stop(
    "calibration_dependence_by_split is missing expected columns: ",
    paste(missing_metrics, collapse = ", "),
    call. = FALSE
  )
}
calibration_summary <- do.call(rbind, lapply(datasets, function(dataset) {
  cell <- calibration_split[calibration_split$dataset == dataset, ]
  row <- data.frame(
    dataset = dataset, n_splits = nrow(cell), stringsAsFactors = FALSE
  )
  for (metric in calibration_metrics) {
    value <- cell[[metric]]
    row[[paste0(metric, "_mean")]] <- mean(value)
    row[[paste0(metric, "_median")]] <- stats::median(value)
    row[[paste0(metric, "_min")]] <- min(value)
    row[[paste0(metric, "_max")]] <- max(value)
  }
  row
}))
write.csv(
  calibration_summary,
  file.path(output_dir, "combined_calibration_dependence_summary.csv"),
  row.names = FALSE
)
write.csv(
  calibration_null,
  file.path(output_dir, "combined_calibration_null_validation.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8. Random baseline sensitivity: straight binds. The by-split table already
# carries the within-split Monte Carlo SE (mc_se) and the summary table the
# across-split sd; combining them is the whole job here.
# -----------------------------------------------------------------------------
write.csv(
  random_summary,
  file.path(output_dir, "combined_random30_summary.csv"),
  row.names = FALSE
)
write.csv(
  random_by_split,
  file.path(output_dir, "combined_random30_by_split.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 9. QA rollup: one row per dataset x analysis with pass counts.
# -----------------------------------------------------------------------------
qa_files <- c(
  ablation = "ablation_qa",
  calibration = "calibration_qa",
  random30 = "random_qa"
)
qa_rollup <- do.call(rbind, lapply(datasets, function(dataset) {
  do.call(rbind, lapply(names(qa_files), function(analysis) {
    path <- file.path(
      dataset_dir(dataset), expected_files[[qa_files[[analysis]]]]
    )
    value <- read.csv(path, stringsAsFactors = FALSE)
    data.frame(
      dataset = dataset,
      analysis = analysis,
      qa_file = basename(path),
      n_rows = nrow(value),
      n_passed = sum(value$passed),
      all_passed = isTRUE(all(value$passed)),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(
  qa_rollup, file.path(output_dir, "combined_QA_rollup.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Final report to stdout.
# -----------------------------------------------------------------------------
outputs <- c(
  "combined_run_status.csv",
  "combined_ablation_vs_random30_by_split.csv",
  "combined_ablation_vs_random30_summary.csv",
  "combined_dataset_level_k10_k20_k50.csv",
  "combined_ablation_vs_current.csv",
  "combined_component_contributions.csv",
  "combined_calibration_dependence_by_split.csv",
  "combined_calibration_dependence_summary.csv",
  "combined_calibration_null_validation.csv",
  "combined_random30_summary.csv",
  "combined_random30_by_split.csv",
  "combined_QA_rollup.csv"
)
cat("\nDone. Combined", length(datasets), "datasets into", output_dir, "\n")
cat("Files written:\n")
for (name in outputs) cat("  ", name, "\n")
cat("Rows: ablation-by-split", nrow(ablation_merged),
    "| dataset x variant x k summary", nrow(ablation_summary),
    "| variant-vs-current cells", nrow(vs_current),
    "| calibration by split", nrow(calibration_split), "\n")
cat("Primary metric: delta_random_30draw (AUC minus mean of 30 Random",
    "draws). All p-values descriptive only.\n")
