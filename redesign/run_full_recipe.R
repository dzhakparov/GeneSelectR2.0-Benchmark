#!/usr/bin/env Rscript
# ==============================================================================
#  Package 1 of the full-benchmark plan: gate vs no gate at the FULL
#  incumbent recipe. Reuses the cached split matrices and module gates from
#  the corrected grouped-benchmark directory (same splits, same pools), but fits
#  with the two incumbent extras restored:
#    - calibration_mode = "evidence_ratio" with the incumbent's null size
#      (calibration_n_permutations = 20, calibration_null_B = 20)
#    - alpha grid c(0.5, 1.0), best alpha by GeneSelectR's internal OOB AUC
#      (fit$cv_results$mean_auc), exactly as the incumbent's
#      rank_with_geneselectr
#
#  Arms per split: full_ungrouped (var2000 pool) and full_grouped (module-
#  gate certified genes first; empty gate -> ungrouped ranking, flagged).
#
#  Checkpointing is PER (split, arm, alpha): one calibration fit takes most
#  of a 300 s tool call, so the two alphas are separate jobs. The price is
#  that each alpha recomputes the calibration null (the in-memory cache
#  cannot cross processes); that is accepted.
#
#  Usage: Rscript redesign/run_full_recipe.R <cohort> <stage> [budget_sec]
#    cohort: imvigor210 | sosall     stage: fit | eval | report | all
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
cohort <- if (length(args) >= 1) args[1] else stop("cohort required")
stage  <- if (length(args) >= 2) args[2] else "all"
budget <- if (length(args) >= 3) as.numeric(args[3]) else 270
stopifnot(cohort %in% c("imvigor210", "sosall"),
          stage %in% c("fit", "eval", "report", "all"))

t_start <- proc.time()[["elapsed"]]
over_budget <- function() (proc.time()[["elapsed"]] - t_start) > budget

suppressPackageStartupMessages(library(glmnet))
options(warn = 1)
for (f in list.files("package/GeneSelectR/R", full.names = TRUE)) source(f)
for (f in c("bio_prior.R", "evaluator.R", "sosall_data.R",
            "imvigor210_data.R", "run_provenance.R")) {
  source(file.path("redesign", "R", f))
}

src_dir  <- file.path(redesign_results_root(), "grouped_benchmark", cohort)
out_dir  <- file.path(redesign_results_root(), "full_recipe", cohort)
require_redesign_run(src_dir)

source_base_rds <- file.path(
  src_dir, if (cohort == "imvigor210") "imvigor210_base.rds" else
    "sosall_base.rds")
if (!file.exists(source_base_rds)) {
  stop("Missing grouped-benchmark cohort data: ", source_base_rds,
       call. = FALSE)
}
source_data <- readRDS(source_base_rds)
outcome <- source_data$outcome

split_files <- sort(list.files(src_dir, pattern = "^split_r.*\\.rds$",
                               full.names = TRUE))
if (length(split_files) != 15) {
  stop("Grouped benchmark must contain all 15 split files before the full ",
       "recipe is run.", call. = FALSE)
}

random_seed <- 42
panel_sizes <- c(10, 20, 50, 100, 200, 500)
alpha_grid <- c(0.5, 1.0)
n_workers <- redesign_worker_count()

prepare_redesign_run(
  out_dir,
  config = list(cohort = cohort, random_seed = random_seed,
                panel_sizes = panel_sizes, alpha_grid = alpha_grid,
                B = 50L, calibration_permutations = 20L,
                calibration_null_B = 20L, n_workers = n_workers,
                outcome_hash = redesign_object_hash(outcome)),
  source_files = c(list.files("package/GeneSelectR/R", full.names = TRUE),
                   file.path("redesign", "R",
                             c("bio_prior.R", "evaluator.R",
                               "sosall_data.R", "imvigor210_data.R",
                               "run_provenance.R")),
                   "redesign/run_full_recipe.R"),
  input_files = c(file.path(src_dir, "run_manifest.rds"), source_base_rds,
                  split_files))

fit_one_alpha <- function(X, y, alpha_value) {
  geneselectr2_fit(
    X, y, gate_method = "none", B = 50, subsample_scheme = "kfold",
    subsample_k_folds = 5, utility_method = "instance_shap",
    components = c("stability", "utility"), score_formula = "geometric",
    calibration_mode = "evidence_ratio", calibration_n_permutations = 20,
    calibration_null_B = 20,
    #  The default uses two workers so concurrent tasks retain capacity.
    alpha = alpha_value, n_cores = n_workers, random_seed = random_seed,
    use_cache = TRUE, verbose = FALSE)
}

#  The outcome is copied from the grouped benchmark and retained for the
#  downstream redesign variants.
base_rds <- file.path(out_dir, "base_outcome.rds")
if (file.exists(base_rds)) {
  if (!identical(readRDS(base_rds), outcome)) {
    stop("Cached full-recipe outcome does not match the grouped benchmark.",
         call. = FALSE)
  }
} else {
  saveRDS(outcome, base_rds)
}

# ---- fit ------------------------------------------------------------------
if (stage %in% c("fit", "all")) {
  for (sf in split_files) {
    sp <- readRDS(sf)
    m <- regmatches(basename(sf),
                    regexec("split_r(\\d+)_f(\\d+)\\.rds",
                            basename(sf)))[[1]]
    rep_idx <- as.integer(m[2]); fold_idx <- as.integer(m[3])
    y_train <- outcome[sp$train_idx]

    for (arm in c("full_ungrouped", "full_grouped")) {
      rank_path <- file.path(
        out_dir, sprintf("ranking_r%d_f%d_%s.csv", rep_idx, fold_idx, arm))
      if (file.exists(rank_path)) next

      pool_genes <- sp$pools$var2000
      if (arm == "full_grouped") {
        pool_genes <- intersect(pool_genes,
                                sp$gates$var2000$certified_genes$q20)
        if (length(pool_genes) < 10) {
          #  Defined fallback: empty gate -> ungrouped ranking, flagged.
          write.csv(data.frame(gene = character(0)), rank_path,
                    row.names = FALSE)
          write.csv(data.frame(key = c("gate_empty", "pool_size", "alpha"),
                               value = c(TRUE, length(pool_genes), NA)),
                    sub("\\.csv$", "_meta.csv", rank_path), row.names = FALSE)
          cat(sprintf("[%s r%d f%d] full_grouped: gate empty -> fallback\n",
                      cohort, rep_idx, fold_idx))
          next
        }
      }

      alpha_paths <- file.path(
        out_dir, sprintf("fit_r%d_f%d_%s_a%s.rds", rep_idx, fold_idx, arm,
                         gsub("\\.", "p", as.character(alpha_grid))))
      for (i in seq_along(alpha_grid)) {
        if (file.exists(alpha_paths[i])) next
        if (over_budget()) { cat("[budget] stop in fit\n"); quit(save = "no") }
        std <- standardise_split(sp$train_raw[, pool_genes, drop = FALSE],
                                 sp$test_raw[, pool_genes, drop = FALSE])
        t0 <- proc.time()[["elapsed"]]
        fit <- fit_one_alpha(std$train, y_train, alpha_grid[i])
        el <- proc.time()[["elapsed"]] - t0
        saveRDS(list(fit = fit, alpha = alpha_grid[i], seconds = el),
                alpha_paths[i])
        #  The calibration-null cache is keyed on the data matrix; the next
        #  job reuses it only if it is the same arm (same matrix).
        cat(sprintf("[%s r%d f%d] %s alpha=%.1f: %.0f s (OOB AUC %.3f)\n",
                    cohort, rep_idx, fold_idx, arm, alpha_grid[i], el,
                    fit$cv_results$mean_auc))
      }
      if (!all(file.exists(alpha_paths))) next   # budget stop mid-arm

      fits <- lapply(alpha_paths, readRDS)
      aucs <- vapply(fits, function(f) f$fit$cv_results$mean_auc,
                     numeric(1))
      best <- which.max(aucs)
      gs <- fits[[best]]$fit$gene_scores
      write.csv(gs, rank_path, row.names = FALSE)
      write.csv(data.frame(
        key = c("gate_empty", "pool_size", "alpha", "auc_alpha0.5",
                "auc_alpha1.0", "fit_seconds"),
        value = c(FALSE, length(pool_genes), alpha_grid[best],
                  round(aucs[1], 4), round(aucs[2], 4),
                  round(sum(vapply(fits, function(f) f$seconds,
                                   numeric(1))), 1))),
        sub("\\.csv$", "_meta.csv", rank_path), row.names = FALSE)
      cat(sprintf("[%s r%d f%d] %s done: %d genes, alpha=%.1f\n",
                  cohort, rep_idx, fold_idx, arm, length(pool_genes),
                  alpha_grid[best]))
      clear_run_cache()
    }
  }
}

# ---- eval -----------------------------------------------------------------
if (stage %in% c("eval", "all")) {
  eval_path <- file.path(out_dir, "eval_results.csv")
  done <- if (file.exists(eval_path)) read.csv(eval_path) else data.frame()
  done_key <- if (nrow(done) > 0)
    paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

  for (sf in split_files) {
    sp <- readRDS(sf)
    m <- regmatches(basename(sf),
                    regexec("split_r(\\d+)_f(\\d+)\\.rds",
                            basename(sf)))[[1]]
    rep_idx <- as.integer(m[2]); fold_idx <- as.integer(m[3])
    y_train <- outcome[sp$train_idx]
    y_test  <- outcome[sp$test_idx]

    rank_ungrouped <- file.path(
      out_dir, sprintf("ranking_r%d_f%d_full_ungrouped.csv",
                       rep_idx, fold_idx))
    rank_grouped <- file.path(
      out_dir, sprintf("ranking_r%d_f%d_full_grouped.csv", rep_idx, fold_idx))
    if (!file.exists(rank_ungrouped) || !file.exists(rank_grouped)) next

    std <- standardise_split(sp$train_raw[, sp$pools$var2000, drop = FALSE],
                             sp$test_raw[, sp$pools$var2000, drop = FALSE])

    for (arm in c("full_ungrouped", "full_grouped")) {
      rk <- read.csv(if (arm == "full_ungrouped") rank_ungrouped
                     else rank_grouped)
      fallback <- FALSE
      if (nrow(rk) == 0) {
        rk <- read.csv(rank_ungrouped)
        fallback <- TRUE
      }
      ranking <- rk$gene
      for (k in panel_sizes) {
        key <- paste(rep_idx, fold_idx, arm, k)
        if (key %in% done_key) next
        if (over_budget()) { cat("[budget] stop in eval\n")
          quit(save = "no") }
        panel <- head(ranking, min(k, length(ranking)))
        scores <- predict_with_ensemble(
          std$train[, panel, drop = FALSE], y_train,
          std$test[, panel, drop = FALSE])
        row <- data.frame(repeat_idx = rep_idx, fold_idx = fold_idx,
                          arm = arm, k = k,
                          AUC = bench_auc(y_test, scores),
                          n_panel = length(panel), gate_fallback = fallback,
                          evaluator_components = paste(
                            attr(scores, "components_used"), collapse = "+"))
        write.table(row, eval_path, append = file.exists(eval_path),
                    sep = ",", row.names = FALSE,
                    col.names = !file.exists(eval_path))
        done_key <- c(done_key, key)
        cat(sprintf("[%s r%d f%d] %s k=%d AUC=%.3f%s\n", cohort, rep_idx,
                    fold_idx, arm, k, row$AUC,
                    ifelse(fallback, " (fallback)", "")))
      }
    }


    #  Matched Random is generated in this corrected run. Keeping it in the
    #  same evaluation table makes the primary delta self-contained and keeps
    #  its provenance tied to the evaluator used for the GeneSelectR arms.
    for (k in panel_sizes) {
      key <- paste(rep_idx, fold_idx, "Random", k)
      if (key %in% done_key) next
      if (over_budget()) {
        cat("[budget] stop in eval\n")
        quit(save = "no")
      }
      pool_genes <- sp$pools$var2000
      aucs <- random_panel_aucs(
        std$train, std$test, y_train, y_test, pool_genes, k, n_draws = 3,
        seed = 99 + 1000 * rep_idx + fold_idx)
      row <- data.frame(
        repeat_idx = rep_idx, fold_idx = fold_idx, arm = "Random", k = k,
        AUC = mean(aucs), n_panel = min(k, length(pool_genes)),
        gate_fallback = FALSE,
        evaluator_components = "glmnet+xgboost+ranger")
      write.table(row, eval_path, append = file.exists(eval_path), sep = ",",
                  row.names = FALSE, col.names = !file.exists(eval_path))
      done_key <- c(done_key, key)
    }
  }
}

# ---- report ---------------------------------------------------------------
if (stage %in% c("report", "all")) {
  ev <- read.csv(file.path(out_dir, "eval_results.csv"))
  rnd <- ev[ev$arm == "Random", ]

  lines <- c(sprintf("# Full-recipe gate test — %s", cohort), "",
    "Recipe: evidence-ratio calibration (20 permutations x 20 fits) + alpha grid",
    "(0.5/1.0 by internal OOB AUC). Delta vs the matched-pool Random",
    "baseline from the grouped benchmark (same splits, same pool).", "",
    "| arm | k | mean AUC | Random | delta | splits | fallbacks |",
    "|---|---|---|---|---|---|---|")
  for (arm in c("full_ungrouped", "full_grouped")) {
    for (k in sort(unique(ev$k))) {
      a <- ev[ev$arm == arm & ev$k == k, ]
      r <- rnd[rnd$k == k, ]
      lines <- c(lines, sprintf(
        "| %s | %d | %.3f | %.3f | %+.3f | %d | %d |",
        arm, k, mean(a$AUC), mean(r$AUC), mean(a$AUC) - mean(r$AUC),
        nrow(a), sum(a$gate_fallback)))
    }
  }
  writeLines(lines, file.path(out_dir, "REPORT.md"))
  cat(sprintf("[report] wrote %s\n", file.path(out_dir, "REPORT.md")))
}

cat(sprintf("\nDone stage=%s cohort=%s (%.0f s)\n", stage, cohort,
            proc.time()[["elapsed"]] - t_start))
