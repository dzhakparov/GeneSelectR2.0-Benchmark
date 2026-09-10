#!/usr/bin/env Rscript
# ==============================================================================
#  Grouped-benchmark: biological priors (pool + module gate) on two cohorts.
#
#  Design (user requirements):
#    - Biology enters ONLY as a prior: (a) Hallmark candidate pool, (b) a
#      pre-selection module gate. NO post-hoc biological scoring anywhere.
#    - Ranking recipe R for every new arm: package geneselectr2_fit with
#      components = c("stability","utility"), geometric mean, instance_shap,
#      B = 50, kfold scheme, alpha = 0.5, calibration_mode = "percentile".
#      The percentile mode skips the permutation calibration null (97% of the
#      incumbent fit cost); incumbent numbers with full calibration are read
#      from saved results as external anchors, so like-for-like calibration
#      differences are stated in the report, not hidden.
#    - Same outer splits as the incumbents: 3 repeats x 5 folds, stratified,
#      seed = 42 + 1000 * repeat_idx.
#    - Same evaluator as the incumbents (redesign/R/evaluator.R).
#
#  Arms per split:
#    var2000_ungrouped   ranking R on the incumbent top-2000-variance pool
#    var2000_grouped     pool restricted to certified Hallmark modules first,
#                        then ranking R (empty gate -> ungrouped, flagged)
#    bio_ungrouped       ranking R on the Hallmark-union pool (no variance cap)
#    bio_grouped         bio pool + module gate, then ranking R
#    varlarge_ungrouped  ranking R on top-|bio| variance genes (size control)
#
#  Usage (from project root, terminal):
#    Rscript redesign/run_grouped_benchmark.R <cohort> <stage> [budget_sec]
#      cohort: imvigor210 | sosall | both
#      stage : prep | fit | eval | report | all
#      budget_sec: stop accepting new jobs after this many seconds (def 270)
#
#  Everything checkpoints per job and resumes; re-running never redoes work.
#  Outputs live under the corrected result root returned by
#  redesign_results_root().
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
cohort_arg <- if (length(args) >= 1) args[1] else "both"
stage_arg  <- if (length(args) >= 2) args[2] else "all"
budget     <- if (length(args) >= 3) as.numeric(args[3]) else 270
stopifnot(cohort_arg %in% c("imvigor210", "sosall", "both"),
          stage_arg %in% c("prep", "fit", "eval", "report", "all"))

t_start <- proc.time()[["elapsed"]]
over_budget <- function() (proc.time()[["elapsed"]] - t_start) > budget

suppressPackageStartupMessages({
  library(glmnet); library(edgeR)
})
options(warn = 1)
for (f in list.files("package/GeneSelectR/R", full.names = TRUE)) source(f)
for (f in c("bio_prior.R", "evaluator.R", "sosall_data.R",
            "imvigor210_data.R", "run_provenance.R")) {
  source(file.path("redesign", "R", f))
}

out_root <- file.path(redesign_results_root(), "grouped_benchmark")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

random_seed <- 42
k_folds <- 5
n_repeats <- 3
panel_sizes <- c(10, 20, 50, 100, 200, 500)
gate_B <- 1000
n_workers <- redesign_worker_count()
pool_jobs <- c("var2000", "bio", "varlarge", "var2000_gated", "bio_gated")

fit_recipe <- function(X, y, seed = random_seed) {
  geneselectr2_fit(
    X, y, gate_method = "none", B = 50, subsample_scheme = "kfold",
    subsample_k_folds = 5, utility_method = "instance_shap",
    components = c("stability", "utility"), score_formula = "geometric",
    calibration_mode = "percentile", alpha = 0.5,
    n_cores = n_workers, random_seed = seed, use_cache = FALSE,
    verbose = FALSE)
}

#  Split definitions identical to the incumbent benchmarks.
split_jobs <- function(outcome) {
  jobs <- list()
  for (rep in seq_len(n_repeats)) {
    folds <- make_stratified_folds(outcome, k_folds,
                                   seed = random_seed + 1000 * rep)
    for (f in seq_len(k_folds)) {
      jobs[[length(jobs) + 1]] <- list(repeat_idx = rep, fold_idx = f,
                                       test_indices = folds[[f]])
    }
  }
  jobs
}

cohorts <- if (cohort_arg == "both") c("imvigor210", "sosall") else cohort_arg

for (cohort in cohorts) {
  cohort_dir <- file.path(out_root, cohort)
  if (cohort == "imvigor210") {
    dat <- load_imvigor210()
    base_rds <- file.path(cohort_dir, "imvigor210_base.rds")
    input_files <- character(0)
  } else {
    dat <- load_sosall()
    base_rds <- file.path(cohort_dir, "sosall_base.rds")
    input_files <- c(file.path("data", "normalized_logcpm.csv"),
                     file.path("data", "metadata.csv"))
  }
  outcome <- dat$outcome
  data_hash <- redesign_object_hash(dat)
  prepare_redesign_run(
    cohort_dir,
    config = list(cohort = cohort, random_seed = random_seed,
                  k_folds = k_folds, n_repeats = n_repeats,
                  panel_sizes = panel_sizes, gate_B = gate_B,
                  pool_jobs = pool_jobs, calibration_mode = "percentile",
                  n_workers = n_workers, data_hash = data_hash),
    source_files = c(list.files("package/GeneSelectR/R", full.names = TRUE),
                     file.path("redesign", "R",
                               c("bio_prior.R", "evaluator.R",
                                 "sosall_data.R", "imvigor210_data.R",
                                 "run_provenance.R")),
                     "redesign/run_grouped_benchmark.R"),
    input_files = input_files)

  if (file.exists(base_rds)) {
    if (!identical(redesign_object_hash(readRDS(base_rds)), data_hash)) {
      stop("Cached cohort data do not match the current source: ", cohort,
           call. = FALSE)
    }
  } else {
    saveRDS(dat, base_rds)
  }
  jobs <- split_jobs(outcome)

  # ---- prep: per-split matrices + pools + module gates -------------------
  if (stage_arg %in% c("prep", "fit", "eval", "all")) {
    for (job in jobs) {
      split_path <- file.path(
        cohort_dir, sprintf("split_r%d_f%d.rds", job$repeat_idx,
                            job$fold_idx))
      if (file.exists(split_path)) next
      if (over_budget()) { cat("[budget] stop in prep\n"); quit(save = "no") }

      train_idx <- setdiff(seq_along(outcome), job$test_indices)
      if (cohort == "imvigor210") {
        norm <- normalise_count_split(dat$raw_counts, train_idx,
                                      job$test_indices)
        train_raw <- norm$train; test_raw <- norm$test
      } else {
        resid <- sosall_residualise_split(dat$expression_matrix,
                                          dat$metadata, train_idx,
                                          job$test_indices)
        train_raw <- resid$train; test_raw <- resid$test
      }

      pools <- make_pools(train_raw, hallmark_union(), top_var = 2000)

      ybin <- as.integer(outcome[train_idx]) - 1L
      gates <- list()
      for (pn in c("var2000", "bio")) {
        cat(sprintf("[%s r%d f%d] gate on %s (p=%d)\n", cohort,
                    job$repeat_idx, job$fold_idx, pn, length(pools[[pn]])))
        gates[[pn]] <- module_gate(train_raw[, pools[[pn]], drop = FALSE],
                                   ybin, B = gate_B, seed = 7)
      }

      saveRDS(list(train_raw = train_raw, test_raw = test_raw,
                   train_idx = train_idx, test_idx = job$test_indices,
                   pools = pools, gates = gates), split_path)
      cat(sprintf("[%s r%d f%d] prep done | pools: var2000=%d bio=%d ",
                  cohort, job$repeat_idx, job$fold_idx,
                  length(pools$var2000), length(pools$bio)))
      cat(sprintf("varlarge=%d | certified q20: var2000=%d bio=%d genes\n",
                  length(pools$varlarge),
                  length(gates$var2000$certified_genes$q20),
                  length(gates$bio$certified_genes$q20)))
    }
  }

  # ---- fit: ranking R per pool job ---------------------------------------
  if (stage_arg %in% c("fit", "all")) {
    for (job in jobs) {
      split_path <- file.path(
        cohort_dir, sprintf("split_r%d_f%d.rds", job$repeat_idx,
                            job$fold_idx))
      if (!file.exists(split_path)) next   # prep not done yet
      sp <- readRDS(split_path)
      y_train <- outcome[sp$train_idx]

      for (pj in pool_jobs) {
        rank_path <- file.path(
          cohort_dir, sprintf("ranking_r%d_f%d_%s.csv", job$repeat_idx,
                              job$fold_idx, pj))
        if (file.exists(rank_path)) next
        if (over_budget()) { cat("[budget] stop in fit\n"); quit(save = "no") }

        base_pool <- sub("_gated$", "", pj)
        gated <- grepl("_gated$", pj)
        pool_genes <- sp$pools[[base_pool]]
        gate_empty <- FALSE
        if (gated) {
          certified <- sp$gates[[base_pool]]$certified_genes$q20
          pool_genes <- intersect(pool_genes, certified)
          if (length(pool_genes) < 10) {
            #  Defined fallback: gate certified (almost) nothing -> the
            #  grouped arm uses the ungrouped ranking, flagged downstream.
            gate_empty <- TRUE
            meta <- data.frame(key = c("gate_empty", "pool_size"),
                               value = c(TRUE, length(pool_genes)))
            write.csv(meta, sub("\\.csv$", "_meta.csv", rank_path),
                      row.names = FALSE)
            write.csv(data.frame(gene = character(0)),
                      rank_path, row.names = FALSE)
            cat(sprintf("[%s r%d f%d] %s: gate empty (%d genes) -> fallback\n",
                        cohort, job$repeat_idx, job$fold_idx, pj,
                        length(pool_genes)))
            next
          }
        }

        std <- standardise_split(sp$train_raw[, pool_genes, drop = FALSE],
                                 sp$test_raw[, pool_genes, drop = FALSE])
        t0 <- proc.time()[["elapsed"]]
        fit <- fit_recipe(std$train, y_train)
        el <- proc.time()[["elapsed"]] - t0

        gs <- fit$gene_scores
        gs <- gs[order(gs$final_score, decreasing = TRUE), ]
        write.csv(gs, rank_path, row.names = FALSE)
        meta <- data.frame(
          key = c("gate_empty", "pool_size", "fit_seconds"),
          value = c(gate_empty, length(pool_genes), round(el, 1)))
        write.csv(meta, sub("\\.csv$", "_meta.csv", rank_path),
                  row.names = FALSE)
        cat(sprintf("[%s r%d f%d] %s: %d genes, %.1f s | top: %s\n",
                    cohort, job$repeat_idx, job$fold_idx, pj,
                    length(pool_genes), el,
                    paste(head(gs$gene, 3), collapse = ", ")))
      }
    }
  }

  # ---- eval: panels + ensemble + random baseline --------------------------
  if (stage_arg %in% c("eval", "all")) {
    eval_path <- file.path(cohort_dir, "eval_results.csv")
    done <- if (file.exists(eval_path)) read.csv(eval_path) else
      data.frame()
    done_key <- if (nrow(done) > 0)
      paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

    arms <- c("var2000_ungrouped", "var2000_grouped", "bio_ungrouped",
              "bio_grouped", "varlarge_ungrouped")

    for (job in jobs) {
      split_path <- file.path(
        cohort_dir, sprintf("split_r%d_f%d.rds", job$repeat_idx,
                            job$fold_idx))
      if (!file.exists(split_path)) next
      sp <- readRDS(split_path)
      y_train <- outcome[sp$train_idx]
      y_test  <- outcome[sp$test_idx]

      #  Rankings for all arms of this split; a missing ranking means the
      #  fit stage has not reached this split yet -> skip the split.
      rank_files <- setNames(file.path(
        cohort_dir, sprintf("ranking_r%d_f%d_%s.csv", job$repeat_idx,
                            job$fold_idx,
                            c("var2000", "var2000_gated", "bio",
                              "bio_gated", "varlarge"))), arms)
      if (!all(file.exists(rank_files))) next

      #  Standardised matrices per base pool (rebuilt from raw, cheap).
      std <- setNames(lapply(c("var2000", "bio", "varlarge"), function(pn) {
        standardise_split(sp$train_raw[, sp$pools[[pn]], drop = FALSE],
                          sp$test_raw[, sp$pools[[pn]], drop = FALSE])
      }), c("var2000", "bio", "varlarge"))

      rows <- list()
      add_row <- function(arm, k, auc, n_panel, fallback) {
        rows[[length(rows) + 1]] <<- data.frame(
          repeat_idx = job$repeat_idx, fold_idx = job$fold_idx,
          arm = arm, k = k, AUC = auc, n_panel = n_panel,
          gate_fallback = fallback,
          evaluator_components = "glmnet+xgboost+ranger")
      }

      flush_rows <- function() {
        #  Write accumulated rows immediately so a budget stop never loses
        #  completed work. done_key is updated alongside so a row is never
        #  written twice within one call.
        if (length(rows) == 0) return(invisible(NULL))
        write.table(do.call(rbind, rows), eval_path,
                    append = file.exists(eval_path),
                    sep = ",", row.names = FALSE,
                    col.names = !file.exists(eval_path))
        rows <<- list()
      }

      for (arm in arms) {
        base_pool <- sub("_ungrouped$|_grouped$", "", arm)
        grouped <- grepl("_grouped$", arm)
        rk <- read.csv(rank_files[[arm]])
        fallback <- FALSE
        if (grouped && nrow(rk) == 0) {
          #  Empty gate: defined fallback is the ungrouped ranking.
          rk <- read.csv(rank_files[[paste0(base_pool, "_ungrouped")]])
          fallback <- TRUE
        }
        ranking <- rk$gene
        for (k in panel_sizes) {
          key <- paste(job$repeat_idx, job$fold_idx, arm, k)
          if (key %in% done_key) next
          if (over_budget()) { flush_rows()
            cat("[budget] stop in eval\n"); quit(save = "no") }
          panel <- head(ranking, min(k, length(ranking)))
          scores <- predict_with_ensemble(
            std[[base_pool]]$train[, panel, drop = FALSE], y_train,
            std[[base_pool]]$test[, panel, drop = FALSE])
          add_row(arm, k, bench_auc(y_test, scores), length(panel), fallback)
          done_key <- c(done_key, key)
          flush_rows()
          cat(sprintf("[%s r%d f%d] %s k=%d AUC=%.3f%s\n", cohort,
                      job$repeat_idx, job$fold_idx, arm, k,
                      bench_auc(y_test, scores),
                      ifelse(fallback, " (fallback)", "")))
        }
      }

      #  Random baseline per base pool (3 draws; documented in the report).
      for (pn in c("var2000", "bio", "varlarge")) {
        for (k in panel_sizes) {
          arm <- paste0("Random_", pn)
          key <- paste(job$repeat_idx, job$fold_idx, arm, k)
          if (key %in% done_key) next
          if (over_budget()) { flush_rows()
            cat("[budget] stop in eval\n"); quit(save = "no") }
          aucs <- random_panel_aucs(
            std[[pn]]$train, std[[pn]]$test, y_train, y_test,
            sp$pools[[pn]], k, n_draws = 3,
            seed = 99 + 1000 * job$repeat_idx + job$fold_idx)
          add_row(arm, k, mean(aucs), min(k, length(sp$pools[[pn]])), FALSE)
          done_key <- c(done_key, key)
          flush_rows()
        }
      }
    }
  }
}

# ---- report ---------------------------------------------------------------
if (stage_arg %in% c("report", "all")) {
  source(file.path("redesign", "R", "grouped_report.R"))
  build_grouped_report(out_root)
}

cat(sprintf("\nDone stage=%s cohort=%s (%.0f s)\n", stage_arg, cohort_arg,
            proc.time()[["elapsed"]] - t_start))
