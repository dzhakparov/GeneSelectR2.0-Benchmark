#!/usr/bin/env Rscript
# ==============================================================================
#  Package 3 of the full-benchmark plan: the five validation datasets.
#
#  Per dataset x 15 splits (3 repeats x 5 folds, incumbent seeds
#  42 + 1000 * repeat_idx):
#    prep  platform branch (counts -> train-fitted TMM+logCPM; normalized ->
#          as-is), confounder residualisation fitted on train, variance pools,
#          Hallmark module gate on the var2000 pool
#    fit   GS_full_ungrouped + GS_full_grouped at the FULL recipe
#          (evidence-ratio calibration, 20 permutation nulls, alpha grid
#          c(0.5, 1.0) picked by internal OOB AUC) + competitor rankings
#          (DGE, LASSO, ElasticNet, mRMR, Boruta, RF) on the same pool
#    eval  incumbent 3-model ensemble at k = 10..500 + matched Random
#          baseline (3 draws)
#
#  Psoriasis (GSE13355) is PAIRED: outer folds move whole patients. The
#  package's inner stability subsamples do not know about patients; that
#  limitation is stated in the report. Psoriasis is report-only (positive
#  control) and is reported descriptively.
#
#  Checkpoints per (split, arm, alpha) / (split, competitor); resumes.
#
#  Usage: Rscript redesign/run_validation_benchmark.R <dataset> <stage> [budget]
#    dataset: any accession registered in redesign/R/validation_data.R
#    stage  : prep | fit | eval | all        budget: seconds (def 270)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
ds_arg <- if (length(args) >= 1) args[1] else stop("dataset required")
stage  <- if (length(args) >= 2) args[2] else "all"
budget <- if (length(args) >= 3) as.numeric(args[3]) else 270
stopifnot(stage %in% c("prep", "fit", "eval", "all"))

t_start <- proc.time()[["elapsed"]]
over_budget <- function() (proc.time()[["elapsed"]] - t_start) > budget

suppressPackageStartupMessages({
  library(glmnet); library(edgeR)
})
options(warn = 1)
for (f in list.files("package/GeneSelectR/R", full.names = TRUE)) source(f)
for (f in c("bio_prior.R", "evaluator.R", "imvigor210_data.R",
            "competitors.R", "validation_data.R", "run_provenance.R")) {
  source(file.path("redesign", "R", f))
}

if (!(ds_arg %in% names(validation_specs))) {
  stop("Unknown validation dataset: ", ds_arg, call. = FALSE)
}

random_seed <- 42
panel_sizes <- c(10, 20, 50, 100, 200, 500)
alpha_grid  <- c(0.5, 1.0)
gate_B      <- 1000
n_workers   <- redesign_worker_count()
fit_n_cores <- suppressWarnings(as.integer(Sys.getenv(
  "GENESELECTR_FIT_N_CORES", as.character(n_workers))))
if (is.na(fit_n_cores) || fit_n_cores < 1L) {
  stop("GENESELECTR_FIT_N_CORES must be a positive integer.", call. = FALSE)
}
parse_arm_override <- function(variable, default) {
  value <- trimws(Sys.getenv(variable, ""))
  if (!nzchar(value)) return(default)
  if (tolower(value) == "none") return(character(0))
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}
gs_arms <- parse_arm_override(
  "GENESELECTR_VALIDATION_GS_ARMS",
  c("GS_full_ungrouped", "GS_full_grouped")
)
comp_arms <- parse_arm_override(
  "GENESELECTR_VALIDATION_COMP_ARMS",
  c("DGE", "LASSO", "ElasticNet", "mRMR", "Boruta", "RF_importance")
)
stopifnot(length(gs_arms) > 0L,
          all(gs_arms %in% c("GS_full_ungrouped", "GS_full_grouped")),
          all(comp_arms %in% c("DGE", "LASSO", "ElasticNet", "mRMR",
                               "Boruta", "RF_importance")))

out_dir <- file.path(redesign_results_root(), "validation_benchmark", ds_arg)
source_files <- c(list.files("package/GeneSelectR/R", full.names = TRUE),
                  file.path("redesign", "R",
                            c("bio_prior.R", "evaluator.R",
                              "imvigor210_data.R", "competitors.R",
                              "validation_data.R", "run_provenance.R")),
                  "redesign/run_validation_benchmark.R")
prepare_redesign_run(
  out_dir,
  config = list(dataset = ds_arg, random_seed = random_seed,
                panel_sizes = panel_sizes, alpha_grid = alpha_grid,
                gate_B = gate_B, B = 50L, calibration_permutations = 20L,
                calibration_null_B = 20L, n_workers = n_workers,
                fit_n_cores = fit_n_cores, gs_arms = gs_arms,
                comp_arms = comp_arms),
  source_files = source_files,
  input_files = c(validation_specs[[ds_arg]]$expr_file,
                  validation_specs[[ds_arg]]$meta_file))

fit_one_alpha <- function(X, y, alpha_value) {
  geneselectr2_fit(
    X, y, gate_method = "none", B = 50, subsample_scheme = "kfold",
    subsample_k_folds = 5, utility_method = "instance_shap",
    components = c("stability", "utility"), score_formula = "geometric",
    calibration_mode = "evidence_ratio", calibration_n_permutations = 20,
    calibration_null_B = 20,
    #  The worker count defaults to two so concurrent benchmark work leaves
    #  memory and CPU capacity for other tasks. It can be changed explicitly.
    alpha = alpha_value, n_cores = fit_n_cores, random_seed = random_seed,
    use_cache = TRUE, verbose = FALSE)
}

# ---- data (cached: CSV parsing is slow at 20k x 480) ------------------------
base_rds <- file.path(out_dir, "base_data.rds")
if (file.exists(base_rds)) {
  dat <- readRDS(base_rds)
} else {
  dat <- load_validation_dataset(ds_arg)
  saveRDS(dat, base_rds)
}
outcome <- dat$outcome
jobs <- validation_split_jobs(dat)
stopifnot(length(jobs) == 15)
split_override <- trimws(Sys.getenv("GENESELECTR_VALIDATION_SPLIT", ""))
if (nzchar(split_override)) {
  split_keys <- vapply(jobs, function(job) sprintf("r%df%d", job$repeat_idx,
                                                    job$fold_idx),
                       character(1))
  if (!(split_override %in% split_keys)) {
    stop("Unknown validation split: ", split_override, call. = FALSE)
  }
  jobs <- jobs[split_keys == split_override]
}

# ---- prep -------------------------------------------------------------------
if (stage %in% c("prep", "all")) {
  for (job in jobs) {
    split_path <- file.path(out_dir, sprintf("split_r%d_f%d.rds",
                                             job$repeat_idx, job$fold_idx))
    if (file.exists(split_path)) next
    if (over_budget()) { cat("[budget] stop in prep\n"); quit(save = "no") }

    train_idx <- setdiff(seq_along(outcome), job$test_indices)

    if (dat$spec$scale == "counts") {
      #  normalise_count_split wants genes x samples.
      norm <- normalise_count_split(t(dat$expr), train_idx, job$test_indices)
      train_raw <- norm$train; test_raw <- norm$test
    } else {
      train_raw <- dat$expr[train_idx, , drop = FALSE]
      test_raw  <- dat$expr[job$test_indices, , drop = FALSE]
    }

    resid <- residualise_split_generic(
      train_raw, test_raw,
      dat$confounders[train_idx, , drop = FALSE],
      dat$confounders[job$test_indices, , drop = FALSE])
    train_raw <- resid$train; test_raw <- resid$test

    pools <- make_pools(train_raw, hallmark_union(), top_var = 2000)

    ybin <- as.integer(outcome[train_idx]) - 1L
    gate <- if ("GS_full_grouped" %in% gs_arms) {
      module_gate(train_raw[, pools$var2000, drop = FALSE],
                  ybin, B = gate_B, seed = 7)
    } else {
      list(certified_genes = list(q20 = character(0)))
    }

    saveRDS(list(train_raw = train_raw, test_raw = test_raw,
                 train_idx = train_idx, test_idx = job$test_indices,
                 pools = pools, gate = gate), split_path)
    cat(sprintf("[%s r%d f%d] prep: var2000=%d bio=%d | certified q20: %d genes\n",
                ds_arg, job$repeat_idx, job$fold_idx,
                length(pools$var2000), length(pools$bio),
                length(gate$certified_genes$q20)))
  }
}

# ---- fit --------------------------------------------------------------------
if (stage %in% c("fit", "all")) {
  for (job in jobs) {
    split_path <- file.path(out_dir, sprintf("split_r%d_f%d.rds",
                                             job$repeat_idx, job$fold_idx))
    if (!file.exists(split_path)) next
    sp <- readRDS(split_path)
    y_train <- outcome[sp$train_idx]

    #  -- GS full recipe, two arms -------------------------------------------
    for (arm in gs_arms) {
      rank_path <- file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv",
                                              job$repeat_idx, job$fold_idx,
                                              arm))
      if (file.exists(rank_path)) next

      pool_genes <- sp$pools$var2000
      if (arm == "GS_full_grouped") {
        pool_genes <- intersect(pool_genes, sp$gate$certified_genes$q20)
        if (length(pool_genes) < 10) {
          #  Defined fallback: empty gate -> ungrouped ranking, flagged.
          write.csv(data.frame(gene = character(0)), rank_path,
                    row.names = FALSE)
          write.csv(data.frame(key = c("gate_empty", "pool_size", "alpha"),
                               value = c(TRUE, length(pool_genes), NA)),
                    sub("\\.csv$", "_meta.csv", rank_path), row.names = FALSE)
          cat(sprintf("[%s r%d f%d] %s: gate empty -> fallback\n",
                      ds_arg, job$repeat_idx, job$fold_idx, arm))
          next
        }
      }

      alpha_paths <- file.path(
        out_dir, sprintf("fit_r%d_f%d_%s_a%s.rds", job$repeat_idx,
                         job$fold_idx, arm,
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
        cat(sprintf("[%s r%d f%d] %s alpha=%.1f: %.0f s (OOB AUC %.3f)\n",
                    ds_arg, job$repeat_idx, job$fold_idx, arm,
                    alpha_grid[i], el, fit$cv_results$mean_auc))
      }
      if (!all(file.exists(alpha_paths))) next   # budget stop mid-arm

      fits <- lapply(alpha_paths, readRDS)
      aucs <- vapply(fits, function(f) f$fit$cv_results$mean_auc, numeric(1))
      best <- which.max(aucs)
      write.csv(fits[[best]]$fit$gene_scores, rank_path, row.names = FALSE)
      write.csv(data.frame(
        key = c("gate_empty", "pool_size", "alpha", "auc_alpha0.5",
                "auc_alpha1.0", "fit_seconds"),
        value = c(FALSE, length(pool_genes), alpha_grid[best],
                  round(aucs[1], 4), round(aucs[2], 4),
                  round(sum(vapply(fits, function(f) f$seconds,
                                   numeric(1))), 1))),
        sub("\\.csv$", "_meta.csv", rank_path), row.names = FALSE)
      cat(sprintf("[%s r%d f%d] %s done: %d genes, alpha=%.1f\n",
                  ds_arg, job$repeat_idx, job$fold_idx, arm,
                  length(pool_genes), alpha_grid[best]))
      clear_run_cache()
    }

    #  -- competitors on the same var2000 pool --------------------------------
    for (arm in comp_arms) {
      rank_path <- file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv",
                                              job$repeat_idx, job$fold_idx,
                                              arm))
      if (file.exists(rank_path)) next
      if (over_budget()) { cat("[budget] stop in fit\n"); quit(save = "no") }

      pool_genes <- sp$pools$var2000
      std <- standardise_split(sp$train_raw[, pool_genes, drop = FALSE],
                               sp$test_raw[, pool_genes, drop = FALSE])
      t0 <- proc.time()[["elapsed"]]
      rk <- switch(arm,
        DGE           = rank_by_differential_expression(std$train, y_train),
        LASSO         = rank_by_lasso(std$train, y_train),
        ElasticNet    = rank_by_elastic_net(std$train, y_train),
        mRMR          = rank_by_mrmr(std$train, y_train),
        Boruta        = rank_by_boruta(std$train, y_train),
        RF_importance = rank_by_random_forest(std$train, y_train,
                                              random_seed = random_seed))
      el <- proc.time()[["elapsed"]] - t0
      stopifnot(length(rk$ranked) == length(pool_genes))
      write.csv(data.frame(gene = rk$ranked), rank_path, row.names = FALSE)
      write.csv(data.frame(key = c("pool_size", "fit_seconds"),
                           value = c(length(pool_genes), round(el, 1))),
                sub("\\.csv$", "_meta.csv", rank_path), row.names = FALSE)
      cat(sprintf("[%s r%d f%d] %s: %.1f s\n", ds_arg, job$repeat_idx,
                  job$fold_idx, arm, el))
    }
  }
}

# ---- eval -------------------------------------------------------------------
if (stage %in% c("eval", "all")) {
  eval_path <- file.path(out_dir, "eval_results.csv")
  done <- if (file.exists(eval_path)) read.csv(eval_path) else data.frame()
  done_key <- if (nrow(done) > 0)
    paste(done$repeat_idx, done$fold_idx, done$arm, done$k) else character(0)

  all_arms <- c(gs_arms, comp_arms)

  for (job in jobs) {
    split_path <- file.path(out_dir, sprintf("split_r%d_f%d.rds",
                                             job$repeat_idx, job$fold_idx))
    if (!file.exists(split_path)) next
    rank_files <- setNames(file.path(
      out_dir, sprintf("ranking_r%d_f%d_%s.csv", job$repeat_idx,
                       job$fold_idx, all_arms)), all_arms)
    if (!all(file.exists(rank_files))) next   # fit stage not there yet

    sp <- readRDS(split_path)
    y_train <- outcome[sp$train_idx]
    y_test  <- outcome[sp$test_idx]
    pool_genes <- sp$pools$var2000
    std <- standardise_split(sp$train_raw[, pool_genes, drop = FALSE],
                             sp$test_raw[, pool_genes, drop = FALSE])

    rows <- list()
    flush_rows <- function() {
      if (length(rows) == 0) return(invisible(NULL))
      write.table(do.call(rbind, rows), eval_path,
                  append = file.exists(eval_path), sep = ",",
                  row.names = FALSE, col.names = !file.exists(eval_path))
      rows <<- list()
    }

    for (arm in all_arms) {
      rk <- read.csv(rank_files[[arm]])
      fallback <- FALSE
      if (arm == "GS_full_grouped" && nrow(rk) == 0) {
        rk <- read.csv(rank_files[["GS_full_ungrouped"]])
        fallback <- TRUE
      }
      ranking <- rk$gene
      for (k in panel_sizes) {
        key <- paste(job$repeat_idx, job$fold_idx, arm, k)
        if (key %in% done_key) next
        if (over_budget()) { flush_rows()
          cat("[budget] stop in eval\n"); quit(save = "no") }
        panel <- head(ranking, min(k, length(ranking)))
        evaluation_seed <- 420000L + job$repeat_idx * 1000L +
          job$fold_idx * 100L + match(k, panel_sizes)
        scores <- predict_with_ensemble(
          std$train[, panel, drop = FALSE], y_train,
          std$test[, panel, drop = FALSE], random_seed = evaluation_seed)
        rows[[length(rows) + 1]] <- data.frame(
          repeat_idx = job$repeat_idx, fold_idx = job$fold_idx,
          arm = arm, k = k, AUC = bench_auc(y_test, scores),
          n_panel = length(panel), gate_fallback = fallback,
          evaluation_seed = evaluation_seed,
          evaluator_components = paste(attr(scores, "components_used"),
                                       collapse = "+"))
        done_key <- c(done_key, key)
        flush_rows()
      }
    }

    #  Matched Random baseline: 3 draws from the same var2000 pool.
    for (k in panel_sizes) {
      key <- paste(job$repeat_idx, job$fold_idx, "Random", k)
      if (key %in% done_key) next
      if (over_budget()) { flush_rows()
        cat("[budget] stop in eval\n"); quit(save = "no") }
      evaluation_seed <- 420000L + job$repeat_idx * 1000L +
        job$fold_idx * 100L + match(k, panel_sizes)
      aucs <- random_panel_aucs(
        std$train, std$test, y_train, y_test, pool_genes, k, n_draws = 3,
        seed = 99 + 1000 * job$repeat_idx + job$fold_idx,
        model_seed = evaluation_seed)
      rows[[length(rows) + 1]] <- data.frame(
        repeat_idx = job$repeat_idx, fold_idx = job$fold_idx,
        arm = "Random", k = k, AUC = mean(aucs),
        n_panel = min(k, length(pool_genes)), gate_fallback = FALSE,
        evaluation_seed = evaluation_seed,
        evaluator_components = "glmnet+xgboost+ranger")
      done_key <- c(done_key, key)
      flush_rows()
    }
  }
}

cat(sprintf("\nDone stage=%s dataset=%s (%.0f s)\n", stage, ds_arg,
            proc.time()[["elapsed"]] - t_start))
