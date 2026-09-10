#!/usr/bin/env Rscript
# Reconstruct descriptive summaries from saved results. No models are fitted.
# Existing result files are read only. All new tables have a separate directory.
CONFIG <- list(
  root = normalizePath("."),
  output = "redesign/results_corrected/review_2026-09-07",
  benchmark = c("GSE101794", "GSE107994", "GSE13355", "GSE65682",
                "GSE69683", "imvigor210", "sosall"),
  additional = c("GSE16879", "GSE91061", "GSE92415", "GSE206285"),
  methods = c("GS_full_ungrouped", "DGE", "LASSO", "ElasticNet", "mRMR",
              "Boruta", "RF_importance"),
  k = c(10L, 20L, 50L, 100L, 200L, 500L), compact_k = c(10L, 20L, 50L),
  tolerance = 1e-12
)
out <- file.path(CONFIG$root, CONFIG$output)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
inputs <- character()
track <- function(path) {
  stopifnot(file.exists(path))
  inputs <<- unique(c(inputs, normalizePath(path)))
  path
}
track(file.path("analysis", "config.R"))
csv <- function(path) read.csv(track(path), stringsAsFactors = FALSE)
rds <- function(path) readRDS(track(path))
save_table <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")),
                                         row.names = FALSE)
root <- file.path(CONFIG$root, "redesign/results_corrected")
dataset_dir <- function(d) file.path(root, if (d %in% CONFIG$benchmark[1:5])
  "validation_benchmark" else "full_recipe", d)
split_dir <- function(d) file.path(root, if (d %in% CONFIG$benchmark[1:5])
  "validation_benchmark" else "grouped_benchmark", d)
key <- function(x, columns) do.call(paste, c(x[columns], sep = "|"))
keys <- c("dataset", "repeat_idx", "fold_idx", "k")
cell_keys <- c(keys, "method")
check_unique <- function(x, columns) stopifnot(!anyDuplicated(key(x, columns)))
means <- function(x, by, values) aggregate(x[values], x[by], mean)
load_dataset <- function(d, filename) {
  x <- csv(file.path(dataset_dir(d), filename)); x$dataset <- d; x
}

# Check all 30 draws in each matched cell before computing its reference mean.
random <- do.call(rbind, lapply(CONFIG$benchmark, load_dataset,
                              filename = "older7_random_baseline_30_draws.csv"))
check_unique(random, c(keys, "draw"))
stopifnot(nrow(random) == 7 * 15 * 6 * 30,
          all(table(key(random, keys)) == 30), all(is.finite(random$AUC)),
          all(random$draw %in% 1:30), all(random$AUC >= 0 & random$AUC <= 1))
random_mean <- means(random, keys, "AUC")
names(random_mean)[names(random_mean) == "AUC"] <- "random_AUC"
fixed <- do.call(rbind, lapply(CONFIG$benchmark, function(d) {
  x <- load_dataset(d, "eval_deterministic.csv")
  x$method <- x$arm
  x[x$method %in% CONFIG$methods, c(cell_keys, "AUC", "evaluation_seed")]
}))
check_unique(fixed, cell_keys)
stopifnot(nrow(fixed) == 7 * 15 * 6 * 7, all(is.finite(fixed$AUC)),
          all(fixed$AUC >= 0 & fixed$AUC <= 1))
fixed <- merge(fixed, random_mean, by = keys)
fixed$delta_random <- fixed$AUC - fixed$random_AUC
save_table(fixed, "benchmark7_fixed_by_split")
fixed_cells <- means(fixed, c("dataset", "method", "k"),
                     c("AUC", "random_AUC", "delta_random"))
save_table(fixed_cells, "benchmark7_fixed_by_dataset_k")
primary <- means(subset(fixed_cells, k %in% CONFIG$compact_k),
                 c("dataset", "method"), c("AUC", "random_AUC", "delta_random"))
primary$rank <- ave(-primary$AUC, primary$dataset, FUN = rank)
save_table(primary, "benchmark7_compact_by_dataset")
overall <- means(primary, "method", c("AUC", "random_AUC", "delta_random", "rank"))
overall$aggregate_AUC_rank <- rank(-overall$AUC)
save_table(overall, "benchmark7_compact_overall")

# Rebuild all nested panel-size choices from the three saved inner folds.
# Inspect patient separation for psoriasis and sample separation everywhere.
adaptive_rows <- list(); inner_rows <- list(); qa_rows <- list(); outer_pools <- list()
for (d in CONFIG$benchmark) {
  groups <- if (d == "GSE13355") rds(file.path(dataset_dir(d), "base_data.rds"))$groups else NULL
  local_rows <- list()
  for (r in 1:3) for (f in 1:5) {
    sp <- rds(file.path(split_dir(d), sprintf("split_r%d_f%d.rds", r, f)))
    outer_pools[[paste(d, r, f, sep = "|")]] <- sp$pools$var2000
    stopifnot(!length(intersect(sp$train_idx, sp$test_idx)))
    if (!is.null(groups)) stopifnot(!length(intersect(groups[sp$train_idx], groups[sp$test_idx])))
    reference_partitions <- NULL
    for (m in CONFIG$methods) {
      parts <- lapply(1:3, function(i) {
        x <- rds(file.path(dataset_dir(d), sprintf(
          "older7_adaptive_inner_r%d_f%d_i%d_%s.rds", r, f, i, m)))
        stopifnot(x$dataset == d, x$repeat_idx == r, x$fold_idx == f,
                  x$inner_idx == i, x$method == m,
                  !anyDuplicated(x$training_indices), !anyDuplicated(x$validation_indices),
                  !length(intersect(x$training_indices, x$validation_indices)),
                  setequal(c(x$training_indices, x$validation_indices), sp$train_idx),
                  length(x$pool) == 2000, !anyDuplicated(x$pool),
                  !anyDuplicated(x$ranking), all(x$ranking %in% x$pool),
                  length(x$ranking) >= max(CONFIG$k),
                  identical(as.integer(x$curve$k), CONFIG$k),
                  all(is.finite(x$curve$AUC)), all(x$curve$AUC >= 0 & x$curve$AUC <= 1))
        if (!is.null(groups)) stopifnot(!length(intersect(
          groups[x$training_indices], groups[x$validation_indices])))
        inner_rows[[length(inner_rows) + 1L]] <<- cbind(
          data.frame(dataset = d, repeat_idx = r, fold_idx = f, inner_idx = i, method = m), x$curve)
        x
      })
      partitions <- lapply(parts, function(x) sort(x$validation_indices))
      stopifnot(!anyDuplicated(unlist(partitions)),
                setequal(unlist(partitions), sp$train_idx))
      if (is.null(reference_partitions)) reference_partitions <- partitions
      stopifnot(identical(reference_partitions, partitions))
      mat <- do.call(cbind, lapply(parts, function(x) x$curve$AUC))
      # Match the original mean() reduction exactly: floating-point ties in
      # the best inner AUC can otherwise change the selected reference SE.
      mu <- apply(mat, 1, mean); se <- apply(mat, 1, sd) / sqrt(3)
      best <- which.max(mu); threshold <- mu[best] - se[best]
      chosen <- min(CONFIG$k[mu >= threshold])
      z <- fixed[fixed$dataset == d & fixed$repeat_idx == r &
                   fixed$fold_idx == f & fixed$method == m & fixed$k == chosen, ]
      stopifnot(nrow(z) == 1)
      local_rows[[length(local_rows) + 1L]] <- data.frame(
        dataset = d, repeat_idx = r, fold_idx = f, method = m, chosen_k = chosen,
        best_inner_k = CONFIG$k[best], best_inner_mean_AUC = mu[best],
        best_inner_se_AUC = se[best], one_se_threshold = threshold,
        outer_AUC = z$AUC, random_AUC = z$random_AUC,
        delta_random = z$delta_random, outer_evaluation_seed = z$evaluation_seed)
    }
  }
  local <- do.call(rbind, local_rows)
  oldpath <- file.path(dataset_dir(d), "older7_adaptive_nested_results.csv")
  err <- NA_real_
  if (file.exists(oldpath)) {
    old <- csv(oldpath)
    match_keys <- c("repeat_idx", "fold_idx", "method")
    check_unique(old, match_keys)
    old <- old[match(key(local, match_keys), key(old, match_keys)), ]
    cols <- c("chosen_k", "best_inner_k", "best_inner_mean_AUC", "best_inner_se_AUC",
              "one_se_threshold", "outer_AUC", "outer_evaluation_seed")
    err <- max(abs(as.matrix(local[cols]) - as.matrix(old[cols])))
    stopifnot(is.finite(err), err <= CONFIG$tolerance)
  }
  adaptive_rows[[d]] <- local
  qa_rows[[d]] <- data.frame(dataset = d, inner_checkpoints = 315,
    reconstructed_outer_rows = nrow(local), split_separation_passed = TRUE,
    psoriasis_patient_separation_checked = d == "GSE13355",
    previous_assembled_table = file.exists(oldpath), previous_max_absolute_error = err)
  cat("Validated adaptive results:", d, "\n"); gc(verbose = FALSE)
}
adaptive <- do.call(rbind, adaptive_rows)
save_table(adaptive, "benchmark7_adaptive_reconstructed")
save_table(do.call(rbind, inner_rows), "benchmark7_adaptive_inner_folds")
save_table(do.call(rbind, qa_rows), "benchmark7_adaptive_QA")
save_table(as.data.frame(xtabs(~ dataset + method + chosen_k, adaptive)), "benchmark7_adaptive_size_counts")
adataset <- means(adaptive, c("dataset", "method"), c("chosen_k", "outer_AUC", "random_AUC", "delta_random"))
save_table(adataset, "benchmark7_adaptive_by_dataset")
save_table(means(adataset, "method", c("chosen_k", "outer_AUC", "random_AUC", "delta_random")), "benchmark7_adaptive_overall")
acompare <- merge(adaptive, fixed, by = c("dataset", "repeat_idx", "fold_idx", "method"),
                  suffixes = c("_adaptive", "_fixed"))
acompare$AUC_change <- acompare$outer_AUC - acompare$AUC
acompare$random_change <- acompare$random_AUC_adaptive - acompare$random_AUC_fixed
acompare$adjusted_change <- acompare$delta_random_adaptive - acompare$delta_random_fixed
stopifnot(max(abs(acompare$adjusted_change - (acompare$AUC_change - acompare$random_change))) < 1e-12)
save_table(acompare, "benchmark7_adaptive_vs_fixed_by_split")
acompare_ds <- means(acompare, c("dataset", "method", "k"), c("AUC_change", "random_change", "adjusted_change"))
save_table(acompare_ds, "benchmark7_adaptive_vs_fixed_by_dataset")
save_table(means(acompare_ds, c("method", "k"), c("AUC_change", "random_change", "adjusted_change")), "benchmark7_adaptive_vs_fixed_overall")

# Component and calibration comparisons retain the saved selected alpha.
remaining <- file.path(root, "older7_remaining_analyses_2026-09-02")
ablation <- csv(file.path(remaining, "combined_ablation_vs_random30_by_split.csv"))
calibration <- csv(file.path(root, "older7_calibration_onoff_2026-09-03", "older7_calibration_onoff_by_split.csv"))
ab_current <- subset(ablation, variant == "GS_current")
cal_current <- subset(calibration, variant == "GS_current")
duplicate_check <- merge(ab_current, cal_current, by = keys)
stopifnot(nrow(duplicate_check) == 630,
          max(abs(duplicate_check$AUC.x - duplicate_check$AUC.y)) < 1e-12)
variants <- rbind(ablation[c(keys, "variant", "AUC")],
                  subset(calibration, variant != "GS_current")[c(keys, "variant", "AUC")])
check_unique(variants, c(keys, "variant"))
variants <- merge(variants, random_mean, by = keys)
variants$delta_random <- variants$AUC - variants$random_AUC
vds <- means(subset(variants, k %in% CONFIG$compact_k), c("dataset", "variant"), c("AUC", "delta_random"))
save_table(vds, "benchmark7_components_calibration_compact_by_dataset")
save_table(means(vds, "variant", c("AUC", "delta_random")), "benchmark7_components_calibration_compact_overall")
current <- subset(calibration, variant == "GS_current")
replay <- merge(current, subset(fixed, method == "GS_full_ungrouped"), by = keys)
stopifnot(nrow(replay) == 630, max(abs(replay$AUC.x - replay$AUC.y)) < 1e-12)
cal_pairs <- c(GS_raw_combination = "GS_current", GS_raw_recurrence = "GS_recurrence_only", GS_raw_SHAPxMI = "GS_SHAPxMI")
cp <- do.call(rbind, lapply(names(cal_pairs), function(raw) {
  x <- merge(subset(vds, variant == raw), subset(vds, variant == cal_pairs[[raw]]), by = "dataset")
  data.frame(dataset = x$dataset, raw_variant = raw, calibrated_variant = cal_pairs[[raw]],
             calibrated_minus_raw_AUC = x$AUC.y - x$AUC.x)
}))
save_table(cp, "benchmark7_calibration_paired_by_dataset")
dependence <- csv(file.path(remaining, "combined_calibration_dependence_by_split.csv"))
save_table(means(dependence, "dataset", c("spearman_calibrated_recurrence_utility_all",
  "spearman_calibrated_recurrence_utility_selected", "spearman_SHAP_MI_all")), "benchmark7_component_dependence")

# Report observed matched-null enrichment at two explicit database thresholds.
ot <- csv(file.path(root, "older7_biology_ot_dense_2026-09-03", "ot_dense_by_split.csv"))
check_unique(ot, c(cell_keys, "cutoff"))
stopifnot(nrow(ot) == 8820, all(is.finite(ot$open_targets_enrichment)))
stopifnot(max(abs(ot$open_targets_enrichment -
  ot$open_targets_sum / ot$open_targets_expected)) < 1e-10)
# Recompute every observed OT sum from the frozen rankings and dense annotation
# files. This checks the annotation calculation independently of its summary.
source(file.path("analysis", "config.R"))
bioconfig <- get_biology_config("older7")
ot_checks <- list(); examples <- list()
for (d in CONFIG$benchmark) {
  seedpath <- sub("_n100_s0.1", "_n3000_s0", bioconfig$association_files[bioconfig$dataset == d], fixed = TRUE)
  seeds <- rds(seedpath)
  stopifnot(!anyDuplicated(seeds$symbol))
  scores <- setNames(seeds$score, seeds$symbol)
  selected20 <- character(); max_error <- 0
  for (r in 1:3) for (f in 1:5) {
    pool <- outer_pools[[paste(d, r, f, sep = "|")]]
    ps <- scores[pool]; ps[is.na(ps)] <- 0; names(ps) <- pool
    for (m in CONFIG$methods) {
      saved_m <- if (d %in% c("imvigor210", "sosall") && m == "GS_full_ungrouped") "full_ungrouped" else m
      ranking <- csv(file.path(dataset_dir(d), sprintf("ranking_r%d_f%d_%s.csv", r, f, saved_m)))$gene
      stopifnot(!anyDuplicated(ranking), setequal(ranking, pool))
      if (m == "GS_full_ungrouped") selected20 <- c(selected20, head(ranking, 20))
      z <- ot[ot$dataset == d & ot$repeat_idx == r & ot$fold_idx == f & ot$method == m, ]
      for (j in seq_len(nrow(z))) {
        filtered <- ps; filtered[filtered < z$cutoff[j]] <- 0
        panel_scores <- filtered[head(ranking, z$k[j])]
        max_error <- max(max_error, abs(sum(panel_scores) - z$open_targets_sum[j]))
        stopifnot(sum(panel_scores > 0) == z$open_targets_overlap[j],
                  sum(filtered > 0) == z$open_targets_pool_targets[j])
      }
    }
  }
  stopifnot(max_error < 1e-12)
  ot_checks[[d]] <- data.frame(dataset = d, checked_rows = 1260,
    max_observed_sum_error = max_error, overlap_and_pool_counts_passed = TRUE)
  counts <- sort(table(selected20), decreasing = TRUE)
  ex <- data.frame(dataset = d, gene = names(counts), n_outer_top20 = as.integer(counts),
                   OT_association_score = unname(scores[names(counts)]))
  examples[[d]] <- ex
}
save_table(do.call(rbind, ot_checks), "benchmark7_dense_OT_reconstruction_QA")
save_table(do.call(rbind, examples), "benchmark7_GS_top20_gene_recurrence_OT")
otcells <- aggregate(open_targets_enrichment ~ dataset + method + k + cutoff,
                     subset(ot, k %in% CONFIG$compact_k), median)
otcells$within_cell_rank <- ave(-otcells$open_targets_enrichment,
  interaction(otcells$dataset, otcells$k, otcells$cutoff), FUN = rank)
otds <- means(otcells, c("dataset", "method", "cutoff"), c("open_targets_enrichment", "within_cell_rank"))
otoverall <- means(otds, c("method", "cutoff"), c("open_targets_enrichment", "within_cell_rank"))
otoverall$rank_of_mean_enrichment <- ave(-otoverall$open_targets_enrichment, otoverall$cutoff, FUN = rank)
save_table(otcells, "benchmark7_dense_OT_compact_cells")
save_table(otds, "benchmark7_dense_OT_by_dataset")
save_table(otoverall, "benchmark7_dense_OT_overall")
bio <- do.call(rbind, lapply(CONFIG$benchmark, load_dataset,
                            filename = "older7_biology_multiaxis.csv"))
check_unique(bio, cell_keys)
stopifnot(nrow(bio) == 4410)
bio_values <- c("GO_semantic_enrichment", "hallmark_enrichment", "string_enrichment")
bcells <- aggregate(bio[bio$k %in% CONFIG$compact_k, bio_values],
                     bio[bio$k %in% CONFIG$compact_k, c("dataset", "method", "k")], median)
bds <- means(bcells, c("dataset", "method"), bio_values)
save_table(bds, "benchmark7_other_biology_by_dataset")
boverall <- means(bds, "method", bio_values)
for (v in bio_values) boverall[[paste0(v, "_rank")]] <- rank(-boverall[[v]])
save_table(boverall, "benchmark7_other_biology_overall")

# Record the actual scope of the stored convergence audit.
convergence_dir <- file.path(root, "older7_glmnet_convergence_2026-09-03")
conv <- csv(file.path(convergence_dir, "summary_by_dataset.csv"))
conv <- subset(conv, dataset != "ALL")
stopifnot(sum(conv$n_fits) == 20400, sum(conv$n_replay_mismatched) == 0)
save_table(conv, "convergence_selector_summary")
save_table(csv(file.path(convergence_dir, "summary_evaluator.csv")), "convergence_evaluator_summary")
deep <- csv(file.path(convergence_dir, "partC_evaluator_deep_refit.csv"))
stopifnot(nrow(deep) == 3, max(abs(deep$auc_delta_deep_vs_saved)) < 1e-12)
save_table(deep, "convergence_deep_evaluator_refits")

# Independently check stored 100-permutation nulls and reproduce utility LOPO.
null_rows <- list(); null_qa <- list()
for (d in c("GSE107994", "GSE13355")) {
  x <- rds(file.path(dataset_dir(d), "older7_calibration_null100_r1f1.rds"))
  stopifnot(x$qa$passed, x$qa$n_permutations == 100)
  a <- x$null_utility; n <- nrow(a); eps <- x$qa$utility_epsilon
  refs <- (matrix(colSums(a), nrow = n, ncol = ncol(a), byrow = TRUE) - a) / (n - 1)
  ratios <- 2^pmax(-4, pmin(4, log2((a + eps) / (refs + eps))))
  z <- x$null_summary[x$null_summary$pillar == "utility" &
    x$null_summary$source == "leave_one_permutation_out_null100", ]
  stopifnot(abs(sd(as.vector(log2(ratios))) - z$sd_log2_ratio) < 1e-12,
    abs(mean(ratios >= .5 & ratios <= 2) - z$fraction_between_half_and_two) < 1e-12)
  null_rows[[d]] <- x$null_summary; null_qa[[d]] <- x$qa
}
save_table(do.call(rbind, null_rows), "null100_validation_rechecked")
save_table(do.call(rbind, null_qa), "null100_QA")
save_table(csv(file.path(remaining, "combined_calibration_null_validation.csv")),
           "null20_validation_reference")
save_table(csv(file.path(root, "older7_null_instability_2026-09-03", "combined_eps_floor_summary.csv")),
           "null_smoothing_diagnostics_reference")

# The four additional cohorts are reported separately with their 30-draw reference.
external_root <- file.path(CONFIG$root, "redesign/results_frozen_external_exact_2026-08-31/validation_benchmark")
external <- do.call(rbind, lapply(CONFIG$additional, function(d) {
  dir <- file.path(external_root, d)
  g <- csv(file.path(dir, "eval_results.csv")); g$method <- g$arm
  g <- subset(g, method == "GS_full_ungrouped")
  c <- csv(file.path(dir, "competitor_eval_results.csv"))
  cols <- c("repeat_idx", "fold_idx", "method", "k", "AUC")
  z <- rbind(g[cols], c[cols]); z$dataset <- d
  check_unique(z, cell_keys); stopifnot(nrow(z) == 630)
  r <- csv(file.path(dir, "random_baseline_sensitivity_30_draws.csv"))
  check_unique(r, c("repeat_idx", "fold_idx", "k", "draw"))
  # The additional-cohort sensitivity run covers only k=10,20,50.
  stopifnot(nrow(r) == 1350, setequal(r$k, CONFIG$compact_k),
            all(table(key(r, c("repeat_idx", "fold_idx", "k"))) == 30))
  r <- means(r, c("repeat_idx", "fold_idx", "k"), "AUC")
  names(r)[names(r) == "AUC"] <- "random_AUC"
  z <- merge(z, r, by = c("repeat_idx", "fold_idx", "k"))
  z$delta_random <- z$AUC - z$random_AUC; z
}))
save_table(external, "additional4_fixed_by_split")
exds <- means(subset(external, k %in% CONFIG$compact_k), c("dataset", "method"), c("AUC", "random_AUC", "delta_random"))
exds$rank <- ave(-exds$AUC, exds$dataset, FUN = rank)
save_table(exds, "additional4_compact_by_dataset")
exoverall <- means(exds, "method", c("AUC", "random_AUC", "delta_random", "rank"))
exoverall$aggregate_AUC_rank <- rank(-exoverall$AUC)
save_table(exoverall, "additional4_compact_overall")
save_table(means(external, c("dataset", "method", "k"), c("AUC", "random_AUC", "delta_random")), "additional4_by_dataset_k")

# Hash inputs used in reconstruction. No saved source file is changed.
source_texts <- c("redesign/run_older7_adaptive.R", "redesign/run_older7_calibration_onoff.R",
  "redesign/run_older7_null100_confirmation.R", "redesign/run_older7_biology.R",
  "redesign/run_older7_biology_ot_dense.R", "redesign/run_older7_glmnet_convergence_check.R",
  "package/GeneSelectR/R/calibration.R")
invisible(lapply(source_texts, track))
inputs <- unique(c(inputs, normalizePath("redesign/review_completed_analyses_2026_09_07.R")))
save_table(data.frame(path = inputs, bytes = file.info(inputs)$size,
                     md5 = unname(tools::md5sum(inputs))), "input_manifest")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("All reconstruction checks passed. Outputs:", out, "\n")
