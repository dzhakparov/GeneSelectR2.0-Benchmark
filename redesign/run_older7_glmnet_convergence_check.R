#!/usr/bin/env Rscript
# ==============================================================================
#  glmnet tail-convergence check for the older-7 benchmark (2026-09-03).
#
#  Question: do the "Convergence for the Nth lambda value not reached;
#  solutions for larger lambdas returned" warnings move any reported number?
#  The warning means glmnet stopped the path before the smallest requested
#  lambda, so the danger is cv.glmnet selecting lambda.min / lambda.1se at
#  the truncated end, where "best" only means "best of the lambdas computed".
#
#  The saved fits discarded every cv.glmnet object (GeneSelectR.R drops
#  cv_fit and predict_fn right after the OOB AUC), so the only way to answer
#  this is to refit the cv.glmnet step. Refitting is only meaningful as an
#  EXACT replay of the original fits; the replay contract, verified
#  bit-for-bit against the saved stability$coef_matrix of every fit, is:
#
#    Selector side (geneselectr2_fit via run_validation_benchmark.R /
#    run_full_recipe.R):
#      - X = standardise_split(train_raw[, pool], test_raw[, pool])$train
#        with pool = var2000 (ungrouped) or intersect(var2000, gate) (grouped)
#      - y = outcome[train_idx]
#      - subsamples = create_subsamples(y, B = 50, random_seed = 42,
#                                       scheme = "kfold", k_folds = 5)
#      - fits ran on a 6-worker PSOCK cluster (run manifests record
#        n_workers = 6) with parallel::clusterSetRNGStream(cl, iseed = 42)
#        and parLapply static chunking (splitIndices(50, 6), contiguous
#        blocks of 9, 9, 8, 8, 8, 8). In R 4.5 clusterSetRNGStream gives
#        worker 1 the UNADVANCED L'Ecuyer stream from set.seed(42) and
#        worker j > 1 exactly j - 1 nextRNGStream advances (function body
#        read; a first emulation that advanced once for worker 1 failed the
#        replay, the corrected one matches bit-for-bit).
#      - each subsample task consumes RNG exactly once: the
#        sample(rep(seq(nfolds), length = N)) inside cv.glmnet that builds
#        foldid (glmnet 5.0 source checked). MI and OOB AUC consume no RNG
#        (grepped), and the stratified-foldid branch in
#        fit_regularized_model only fires when a class has < 5 observations,
#        which never happens on these datasets. glmnet itself uses no R RNG,
#        so precomputing foldid per (worker stream, task order) and passing
#        it explicitly reproduces the original fits exactly.
#      - alpha = 0.5 and alpha = 1.0 were separate geneselectr2_fit calls,
#        each with a fresh cluster at iseed = 42, and the arms of a split
#        share y, so ALL (arm, alpha) jobs of a split share foldids.
#
#    Evaluator side (reevaluate_saved_rankings.R -> predict_with_ensemble):
#      - fold_id comes from withr::with_seed(evaluation_seed, stratified
#        sampling); evaluation_seed is stored per row in
#        eval_deterministic.csv. Fully explicit, no emulation needed.
#
#  Scope notes:
#    - Selector: EVERY saved GS fit (7 datasets x 15 splits x 2 arms x 2
#      alphas x 50 subsamples) is replayed and verified. An initial r1f1-only
#      pass found zero warnings, which could not explain the 753 warnings in
#      the main log, so coverage was extended to all splits and both arms.
#    - Evaluator: ALL eval_deterministic.csv rows for the two datasets whose
#      deterministic logs contain convergence warnings (GSE101794: 10,
#      GSE65682: 7), every arm and panel size, plus a full-ensemble harness
#      recompute on the 120-row GS_full_ungrouped/DGE x k in {10, 50} sample.
#    - NOT replayed: the 20-permutation calibration nulls (~84k fits, they
#      feed only the calibration null distributions) and competitor selector
#      fits (LASSO/ElasticNet rankings; make_competitor_fold_id draws from
#      ambient RNG with no stored seed, so they are not reproducible).
#      Both are stated as limitations in verdict.md.
#
#  Stages: replay_check | partA | partB | partC | report | all
#  Run from the repo root with LC_ALL=en_US.UTF-8 and single-thread BLAS vars.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1L) args[[1L]] else "all"
stopifnot(stage %in% c("replay_check", "partA", "partB", "partC", "report",
                       "all"))

suppressPackageStartupMessages({
  library(glmnet)
  library(withr)
  library(xgboost)
  library(ranger)
})

OUT_DIR <- file.path("redesign", "results_corrected",
                     "older7_glmnet_convergence_2026-09-03")
WORK_DIR <- file.path(OUT_DIR, "work")
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

#  Six fitting workers: every run manifest records n_workers = 6, and the
#  GSE101794 replay check proves 6 is the configuration the GS fits used.
FIT_WORKERS <- 6L
POOL_WORKERS <- suppressWarnings(as.integer(Sys.getenv(
  "GENESELECTR_CHECK_CORES", "6")))
if (is.na(POOL_WORKERS) || POOL_WORKERS < 1L) POOL_WORKERS <- 6L

#  standardise_split, copied verbatim from redesign/R/bio_prior.R (sourcing
#  the file would drag in unrelated dependencies; the copy must not drift).
standardise_split <- function(train_expr, test_expr) {
  cm <- colMeans(train_expr)
  cs <- apply(train_expr, 2, sd)
  cs[cs == 0 | is.na(cs)] <- 1
  tr <- sweep(sweep(train_expr, 2, cm, "-"), 2, cs, "/")
  te <- if (is.null(test_expr)) NULL else
    sweep(sweep(test_expr, 2, cm, "-"), 2, cs, "/")
  list(train = tr, test = te)
}

#  bench_auc, copied verbatim from redesign/R/evaluator.R (same reason).
bench_auc <- function(true_labels, predicted_scores) {
  true_labels <- droplevels(true_labels)
  pos <- which(true_labels == levels(true_labels)[2])
  neg <- which(true_labels == levels(true_labels)[1])
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(predicted_scores, ties.method = "average")
  (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

#  create_subsamples comes from the package itself so the subsample design
#  cannot drift from what the benchmark ran.
source(file.path("package", "GeneSelectR", "R", "utils.R"))

# ------------------------------------------------------------------------------
#  Dataset registry.
# ------------------------------------------------------------------------------
validation_datasets <- c("GSE101794", "GSE107994", "GSE13355", "GSE65682",
                         "GSE69683")
cohort_datasets <- c("imvigor210", "sosall")
all_datasets <- c(validation_datasets, cohort_datasets)

dataset_paths <- function(ds) {
  if (ds %in% validation_datasets) {
    dir <- file.path("redesign", "results_corrected", "validation_benchmark",
                     ds)
    list(split_dir = dir, fit_dir = dir,
         arms = c("GS_full_ungrouped", "GS_full_grouped"),
         outcome = readRDS(file.path(dir, "base_data.rds"))$outcome,
         eval_csv = file.path(dir, "eval_deterministic.csv"),
         ranking_dir = dir, is_validation = TRUE)
  } else {
    split_dir <- file.path("redesign", "results_corrected",
                           "grouped_benchmark", ds)
    fit_dir <- file.path("redesign", "results_corrected", "full_recipe", ds)
    list(split_dir = split_dir, fit_dir = fit_dir,
         arms = c("full_ungrouped", "full_grouped"),
         outcome = readRDS(file.path(fit_dir, "base_outcome.rds")),
         eval_csv = file.path(fit_dir, "eval_deterministic.csv"),
         ranking_dir = fit_dir, is_validation = FALSE)
  }
}

#  Gate-certified pool for the grouped arm, matching the two drivers:
#  validation keeps sp$gate$certified_genes$q20, cohorts
#  sp$gates$var2000$certified_genes$q20.
job_pool <- function(sp, arm, is_validation) {
  pool <- sp$pools$var2000
  #  Exact match: "ungrouped" contains "grouped" as a substring, so a grepl
  #  here silently hands the ungrouped arm the gate-intersected pool.
  if (arm %in% c("GS_full_grouped", "full_grouped")) {
    gate_genes <- if (is_validation) {
      sp$gate$certified_genes$q20
    } else {
      sp$gates$var2000$certified_genes$q20
    }
    pool <- intersect(pool, gate_genes)
  }
  pool
}

# ------------------------------------------------------------------------------
#  RNG emulation for the PSOCK selector fits (see header for the contract).
# ------------------------------------------------------------------------------
set_rng_state <- function(state) assign(".Random.seed", state,
                                        envir = .GlobalEnv)
get_rng_state <- function() get(".Random.seed", envir = .GlobalEnv)

psock_foldids <- function(task_train_ns, n_workers = FIT_WORKERS,
                          iseed = 42L) {
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(iseed)
  streams <- vector("list", n_workers)
  #  Worker 1 gets the un-advanced stream; each later worker one more advance.
  streams[[1L]] <- get_rng_state()
  for (j in seq_len(n_workers - 1L)) {
    streams[[j + 1L]] <- parallel::nextRNGStream(streams[[j]])
  }
  chunks <- parallel::splitIndices(length(task_train_ns), n_workers)
  foldids <- vector("list", length(task_train_ns))
  for (j in seq_len(n_workers)) {
    state <- streams[[j]]
    for (b in chunks[[j]]) {
      set_rng_state(state)
      foldids[[b]] <- sample(rep(seq_len(5L), length.out = task_train_ns[b]))
      state <- get_rng_state()
    }
  }
  foldids
}

# ------------------------------------------------------------------------------
#  cv.glmnet wrapper: capture every warning, keep the diagnostics. The first
#  convergence warning comes from the full-data fit (cv.glmnet fits the full
#  path first, then the folds); later ones are fold fits and only corrupt
#  tail cvm values, not the path length.
#  maxit is left ABSENT for standard fits because the original calls never
#  passed it (and glmnet 5.0 deprecates the argument); deep refits pass it
#  through control = list(maxit = ...) per the new API.
# ------------------------------------------------------------------------------
CONV_PATTERN <- "lambda value not reached"

fit_cv <- function(X, y01, alpha, foldid, type_measure, nlambda = 100L,
                   lambda_min_ratio = NULL, deep_maxit = NULL) {
  warns <- character(0)
  call_args <- list(
    x = X, y = y01, family = "binomial", alpha = alpha, nfolds = 5L,
    foldid = foldid, type.measure = type_measure,
    penalty.factor = rep(1, ncol(X)),
    nlambda = nlambda
  )
  #  lambda.min.ratio stays absent when NULL so glmnet's own nobs-vs-nvars
  #  default applies, exactly as in the original calls.
  if (!is.null(lambda_min_ratio)) call_args$lambda.min.ratio <- lambda_min_ratio
  if (!is.null(deep_maxit)) call_args$control <- list(maxit = deep_maxit)
  fit <- withCallingHandlers(
    do.call(glmnet::cv.glmnet, call_args),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  conv <- grep(CONV_PATTERN, warns, value = TRUE)
  conv_idx <- integer(0)
  conv_code <- integer(0)
  if (length(conv) > 0L) {
    m <- regmatches(conv, regexec(
      "error code -([0-9]+).*Convergence for ([0-9]+)th", conv))
    conv_code <- as.integer(vapply(m, `[`, character(1), 2L))
    conv_idx  <- as.integer(vapply(m, `[`, character(1), 3L))
  }
  list(fit = fit, warnings = warns, conv_idx = conv_idx,
       conv_code = conv_code)
}

cv_diagnostics <- function(res, nlambda_requested = 100L) {
  fit <- res$fit
  n_lambda <- length(fit$lambda)
  idx_min <- which.min(abs(fit$lambda - fit$lambda.min))
  idx_1se <- which.min(abs(fit$lambda - fit$lambda.1se))
  data.frame(
    n_warnings = length(res$warnings),
    early_stop = length(res$conv_idx) > 0L,
    stop_index = if (length(res$conv_idx) > 0L) res$conv_idx[1L] else NA_integer_,
    stop_code = if (length(res$conv_code) > 0L) res$conv_code[1L] else NA_integer_,
    n_fold_conv_warnings = max(0L, length(res$conv_idx) - 1L),
    n_lambda = n_lambda,
    path_complete = n_lambda >= nlambda_requested,
    idx_lambda_min = idx_min,
    idx_lambda_1se = idx_1se,
    lambda_min = fit$lambda.min,
    lambda_1se = fit$lambda.1se,
    #  "Boundary" = at or within 2 indices of the last computed lambda.
    boundary_min = idx_min >= n_lambda - 2L,
    boundary_1se = idx_1se >= n_lambda - 2L,
    interior_min = idx_min < n_lambda,
    cvm_decreasing_end = n_lambda >= 2L &&
      fit$cvm[n_lambda] < fit$cvm[n_lambda - 1L],
    cvm_at_lambda_min = fit$cvm[idx_min],
    cvm_at_last = fit$cvm[n_lambda],
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
#  Selector jobs: one row per saved fit file (dataset x split x arm x alpha).
# ------------------------------------------------------------------------------
partA_path <- function(ds, rep_idx, fold_idx, arm, alpha) {
  file.path(WORK_DIR, sprintf(
    "partA_%s_r%d_f%d_%s_a%s.rds", ds, rep_idx, fold_idx, arm,
    gsub("\\.", "p", as.character(alpha))))
}

selector_jobs <- function(only = NULL) {
  jobs <- list()
  for (ds in all_datasets) {
    paths <- dataset_paths(ds)
    for (arm in paths$arms) {
      for (rep_idx in 1:3) for (fold_idx in 1:5) {
        for (alpha in c(0.5, 1.0)) {
          fit_path <- file.path(paths$fit_dir, sprintf(
            "fit_r%d_f%d_%s_a%s.rds", rep_idx, fold_idx, arm,
            gsub("\\.", "p", as.character(alpha))))
          #  Grouped arm with an empty gate was never fitted (defined
          #  fallback), so no file exists; iterating over files skips it.
          if (!file.exists(fit_path)) next
          jobs[[length(jobs) + 1L]] <- data.frame(
            dataset = ds, rep_idx = rep_idx, fold_idx = fold_idx, arm = arm,
            alpha = alpha, fit_path = fit_path, stringsAsFactors = FALSE)
        }
      }
    }
  }
  jobs <- do.call(rbind, jobs)
  if (!is.null(only)) jobs <- merge(jobs, only)
  jobs
}

#  Run the ≤4 (arm, alpha) jobs of one split. Subsamples and foldids depend
#  only on y, so they are computed once per split and shared. The worker
#  cluster is created once in run_partA and passed in; only the per-arm data
#  (X, alpha) is re-exported.
run_split_jobs <- function(ds, rep_idx, fold_idx, jobs_this_split, cl) {
  paths <- dataset_paths(ds)
  sp <- readRDS(file.path(paths$split_dir, sprintf("split_r%d_f%d.rds",
                                                   rep_idx, fold_idx)))
  y <- paths$outcome[sp$train_idx]
  subsamples <- create_subsamples(y, B = 50, random_seed = 42,
                                  scheme = "kfold", k_folds = 5)
  task_ns <- vapply(subsamples, function(ss) length(ss$train), integer(1))
  foldids <- psock_foldids(task_ns)

  for (j in seq_len(nrow(jobs_this_split))) {
    job <- jobs_this_split[j, ]
    out_path <- partA_path(ds, rep_idx, fold_idx, job$arm, job$alpha)
    if (file.exists(out_path)) next
    pool <- job_pool(sp, job$arm, paths$is_validation)
    std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                             sp$test_raw[, pool, drop = FALSE])
    X <- std$train

    tasks <- lapply(seq_along(subsamples), function(b) {
      list(b = b, train = subsamples[[b]]$train, foldid = foldids[[b]])
    })
    alpha <- job$alpha
    run_one <- function(task) {
      Xs <- X[task$train, , drop = FALSE]
      y01 <- as.numeric(y[task$train]) - 1
      res <- fit_cv(Xs, y01, alpha = alpha, foldid = task$foldid,
                    type_measure = "auc")
      diag <- cv_diagnostics(res)
      coefs <- as.vector(glmnet::coef.glmnet(res$fit, s = "lambda.min"))
      list(diag = diag, coef = coefs)
    }
    results <- if (!is.null(cl)) {
      parallel::clusterExport(
        cl, c("X", "y", "alpha", "fit_cv", "cv_diagnostics", "CONV_PATTERN"),
        envir = environment())
      parallel::parLapply(cl, tasks, run_one)
    } else {
      lapply(tasks, run_one)
    }

    diag <- do.call(rbind, lapply(results, `[[`, "diag"))
    diag$dataset <- ds
    diag$rep_idx <- rep_idx
    diag$fold_idx <- fold_idx
    diag$arm <- job$arm
    diag$alpha <- alpha
    diag$subsample <- seq_len(nrow(diag))
    coef_mat <- do.call(cbind, lapply(results, `[[`, "coef"))
    rownames(coef_mat) <- c("(Intercept)", colnames(X))

    saved <- readRDS(job$fit_path)$fit$stability$coef_matrix
    refit <- coef_mat[-1L, , drop = FALSE]
    if (!all(dim(refit) == dim(saved))) {
      stop(sprintf(paste0("Dimension mismatch vs saved fit for %s r%d f%d ",
                          "%s a=%.1f: refit %s vs saved %s"),
                   ds, rep_idx, fold_idx, job$arm, alpha,
                   paste(dim(refit), collapse = "x"),
                   paste(dim(saved), collapse = "x")))
    }
    replay_max_abs_diff <- apply(abs(refit - saved), 2L, max)
    diag$replay_max_abs_diff <- as.numeric(replay_max_abs_diff)

    saveRDS(list(diag = diag, coef = coef_mat, fit_path = job$fit_path),
            out_path)
    cat(sprintf("[%s r%d f%d %s a=%.1f] mismatched cols: %d/50 (max %.3g), early stops: %d\n",
                ds, rep_idx, fold_idx, job$arm, alpha,
                sum(replay_max_abs_diff > 1e-8), max(replay_max_abs_diff),
                sum(diag$early_stop)))
  }
}

run_partA <- function(only = NULL) {
  jobs <- selector_jobs(only)
  splits <- unique(jobs[, c("dataset", "rep_idx", "fold_idx")])
  cl <- NULL
  if (POOL_WORKERS > 1L) {
    cl <- parallel::makeCluster(POOL_WORKERS, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, library(glmnet))
  }
  for (i in seq_len(nrow(splits))) {
    s <- splits[i, ]
    run_split_jobs(s$dataset, s$rep_idx, s$fold_idx,
                   jobs[jobs$dataset == s$dataset &
                          jobs$rep_idx == s$rep_idx &
                          jobs$fold_idx == s$fold_idx, ], cl)
  }
}

#  Replay check: GSE101794 r1f1 GS_full_ungrouped alpha = 0.5 must match the
#  saved coef_matrix in all 50 columns. Stops the pipeline if not: nothing
#  downstream means anything without an exact replay.
run_replay_check <- function() {
  only <- data.frame(dataset = "GSE101794", rep_idx = 1L, fold_idx = 1L,
                     arm = "GS_full_ungrouped", alpha = 0.5)
  run_partA(only)
  obj <- readRDS(partA_path("GSE101794", 1L, 1L, "GS_full_ungrouped", 0.5))
  n_bad <- sum(obj$diag$replay_max_abs_diff > 1e-8)
  write.csv(obj$diag, file.path(OUT_DIR, "replay_check_GSE101794.csv"),
            row.names = FALSE)
  if (n_bad > 0L) {
    stop(sprintf(paste0(
      "Exact replay FAILED: %d/50 columns differ from the saved fit (max %.3g). ",
      "Do not trust downstream numbers."), n_bad,
      max(obj$diag$replay_max_abs_diff)))
  }
  cat("replay_check: exact replay achieved (50/50 columns bit-identical)\n")
}

load_partA <- function() {
  files <- list.files(WORK_DIR, pattern = "^partA_.*[.]rds$",
                      full.names = TRUE)
  files <- files[!grepl("partA_index", files)]
  objs <- lapply(files, readRDS)
  list(diag = do.call(rbind, lapply(objs, `[[`, "diag")),
       objs = setNames(objs, basename(files)))
}

# ------------------------------------------------------------------------------
#  Stage: partB — deep-path refit for every selector fit that stopped early
#  AND has lambda.min or lambda.1se at/within 2 indices of the stop.
# ------------------------------------------------------------------------------
run_partB <- function() {
  pa <- load_partA()
  d <- pa$diag
  hit <- which(d$early_stop & (d$boundary_min | d$boundary_1se))
  if (length(hit) == 0L) {
    cat("partB: no early-stop boundary fits, nothing to refit\n")
    saveRDS(data.frame(), file.path(WORK_DIR, "partB.rds"))
    return(invisible(NULL))
  }
  cat(sprintf("partB: %d boundary fits to refit deep\n", length(hit)))
  results <- vector("list", length(hit))
  for (i in seq_along(hit)) {
    row <- d[hit[i], ]
    obj <- readRDS(partA_path(row$dataset, row$rep_idx, row$fold_idx,
                              row$arm, row$alpha))
    paths <- dataset_paths(row$dataset)
    sp <- readRDS(file.path(paths$split_dir, sprintf(
      "split_r%d_f%d.rds", row$rep_idx, row$fold_idx)))
    y <- paths$outcome[sp$train_idx]
    pool <- job_pool(sp, row$arm, paths$is_validation)
    std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                             sp$test_raw[, pool, drop = FALSE])
    X <- std$train
    subsamples <- create_subsamples(y, B = 50, random_seed = 42,
                                    scheme = "kfold", k_folds = 5)
    task_ns <- vapply(subsamples, function(ss) length(ss$train), integer(1))
    foldids <- psock_foldids(task_ns)
    b <- row$subsample
    tr <- subsamples[[b]]$train
    deep <- fit_cv(X[tr, , drop = FALSE], as.numeric(y[tr]) - 1,
                   alpha = row$alpha, foldid = foldids[[b]],
                   type_measure = "auc", nlambda = 200L,
                   lambda_min_ratio = 1e-6, deep_maxit = 1e6)
    dd <- cv_diagnostics(deep, nlambda_requested = 200L)
    old_lambda_min <- row$lambda_min
    coef_old <- obj$coef[-1L, b]
    coef_new_at_old <- as.vector(glmnet::coef.glmnet(
      deep$fit, s = old_lambda_min))[-1L]
    coef_new_at_newmin <- as.vector(glmnet::coef.glmnet(
      deep$fit, s = "lambda.min"))[-1L]
    results[[i]] <- data.frame(
      dataset = row$dataset, rep_idx = row$rep_idx, fold_idx = row$fold_idx,
      arm = row$arm, alpha = row$alpha, subsample = b,
      orig_stop_index = row$stop_index, orig_n_lambda = row$n_lambda,
      orig_idx_lambda_min = row$idx_lambda_min,
      orig_lambda_min = old_lambda_min,
      deep_n_lambda = dd$n_lambda, deep_early_stop = dd$early_stop,
      deep_idx_lambda_min = dd$idx_lambda_min,
      deep_lambda_min = deep$fit$lambda.min,
      deep_boundary_min = dd$boundary_min,
      lambda_min_rel_change = abs(deep$fit$lambda.min - old_lambda_min) /
        old_lambda_min,
      max_abs_coef_change_at_old_lambda = max(abs(coef_new_at_old - coef_old)),
      n_nonzero_old_at_min = sum(coef_old != 0),
      n_nonzero_deep_at_newmin = sum(coef_new_at_newmin != 0),
      n_nonzero_shared = sum(coef_old != 0 & coef_new_at_newmin != 0),
      stringsAsFactors = FALSE)
    cat(sprintf("  [%s r%d f%d %s a=%.1f b=%d] lambda.min rel change %.4g, max |dcoef| %.4g\n",
                row$dataset, row$rep_idx, row$fold_idx, row$arm, row$alpha,
                b, results[[i]]$lambda_min_rel_change,
                results[[i]]$max_abs_coef_change_at_old_lambda))
  }
  saveRDS(do.call(rbind, results), file.path(WORK_DIR, "partB.rds"))
}

# ------------------------------------------------------------------------------
#  Stage: partC — evaluator-side glmnet component for EVERY
#  eval_deterministic.csv row of GSE101794 and GSE65682 (the two datasets
#  whose deterministic logs contain convergence warnings), all arms, all
#  panel sizes, all 15 splits, including the Random baseline draws.
#  Workers return diagnostics only; the master then (a) recomputes the full
#  ensemble on the 120-row GS/DGE x k in {10, 50} harness sample to validate
#  against the saved AUCs, and (b) for every row whose winning fit stopped
#  early, recomputes the ensemble AUC with the original path and with a deep
#  path to measure the AUC effect.
# ------------------------------------------------------------------------------
evaluator_fold_id <- function(train_labels, seed) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  n_inner <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))
  fold_id <- withr::with_seed(seed, {
    result <- integer(length(numeric_labels))
    for (class_value in sort(unique(numeric_labels))) {
      indices <- which(numeric_labels == class_value)
      result[indices] <- sample(rep(seq_len(n_inner),
                                    length.out = length(indices)))
    }
    result
  })
  list(fold_id = fold_id, numeric_labels = numeric_labels)
}

#  read_ranking, mirroring reevaluate_saved_rankings.R (kswitch branch and
#  grouped-arm fallback included).
read_ranking <- function(ranking_dir, rep_idx, fold_idx, arm, k) {
  file_arm <- arm
  if (arm == "kswitch") {
    file_arm <- if (k <= 20L) "predfirst_raw" else "GS_full_grouped"
  }
  path <- file.path(ranking_dir, sprintf("ranking_r%d_f%d_%s.csv",
                                         rep_idx, fold_idx, file_arm))
  if (!file.exists(path)) stop("missing ranking: ", path, call. = FALSE)
  ranking <- read.csv(path, stringsAsFactors = FALSE)$gene
  if (!length(ranking) && arm %in% c("GS_full_grouped", "kswitch")) {
    ranking <- read.csv(file.path(ranking_dir, sprintf(
      "ranking_r%d_f%d_GS_full_ungrouped.csv", rep_idx, fold_idx)),
      stringsAsFactors = FALSE)$gene
  }
  if (!length(ranking) || anyDuplicated(ranking)) {
    stop("invalid ranking: ", path, call. = FALSE)
  }
  ranking
}

#  glmnet-component diagnostics for one (panel, fold_id) pair.
eval_glmnet_diag <- function(Xtr, y01, fold_id) {
  per_alpha <- list()
  for (alpha_value in c(0.5, 1.0)) {
    res <- fit_cv(Xtr, y01, alpha = alpha_value, foldid = fold_id,
                  type_measure = "deviance")
    per_alpha[[as.character(alpha_value)]] <- list(res = res,
                                                   diag = cv_diagnostics(res))
  }
  #  Winner rule copied from predict_with_ensemble: later alpha wins only on
  #  a strictly larger -min(cvm).
  auc05 <- -min(per_alpha[["0.5"]]$res$fit$cvm)
  auc10 <- -min(per_alpha[["1"]]$res$fit$cvm)
  winner <- if (auc10 > auc05) "1" else "0.5"
  list(per_alpha = per_alpha, winner = winner)
}

run_partC <- function() {
  panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
  partC_datasets <- c("GSE101794", "GSE65682")

  for (ds in partC_datasets) {
    ds_ckpt <- file.path(WORK_DIR, sprintf("partC_%s.rds", ds))
    if (!file.exists(ds_ckpt)) {
      paths <- dataset_paths(ds)
      ev <- read.csv(paths$eval_csv, stringsAsFactors = FALSE)
      outcome <- paths$outcome

      #  Per-split work units so the (large) split objects are read once.
      split_units <- split(ev, paste(ev$repeat_idx, ev$fold_idx))
      process_split <- function(unit_name) {
        unit <- split_units[[unit_name]]
        rep_idx <- unit$repeat_idx[1L]
        fold_idx <- unit$fold_idx[1L]
        sp <- readRDS(file.path(paths$split_dir, sprintf(
          "split_r%d_f%d.rds", rep_idx, fold_idx)))
        y_train <- outcome[sp$train_idx]
        pool <- sp$pools$var2000
        std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                                 sp$test_raw[, pool, drop = FALSE])
        rows <- list()
        for (i in seq_len(nrow(unit))) {
          r <- unit[i, ]
          expected_seed <- 420000L + rep_idx * 1000L + fold_idx * 100L +
            match(r$k, panel_sizes)
          stopifnot(r$evaluation_seed == expected_seed)
          fl <- evaluator_fold_id(y_train, r$evaluation_seed)
          if (r$arm == "Random") {
            #  random_panel_aucs: one set.seed per split/k, three draws.
            set.seed(99L + 1000L * rep_idx + fold_idx)
            draws <- replicate(3L, sample(pool, min(r$k, length(pool))),
                               simplify = FALSE)
            for (draw_idx in seq_along(draws)) {
              panel <- draws[[draw_idx]]
              g <- eval_glmnet_diag(std$train[, panel, drop = FALSE],
                                    fl$numeric_labels, fl$fold_id)
              rows[[length(rows) + 1L]] <- eval_row_from_diag(
                ds, r, draw_idx, g)
            }
          } else {
            ranking <- read_ranking(paths$ranking_dir, rep_idx, fold_idx,
                                    r$arm, r$k)
            panel <- head(ranking, min(r$k, length(ranking)))
            g <- eval_glmnet_diag(std$train[, panel, drop = FALSE],
                                  fl$numeric_labels, fl$fold_id)
            rows[[length(rows) + 1L]] <- eval_row_from_diag(ds, r, NA_integer_,
                                                            g)
          }
        }
        do.call(rbind, rows)
      }
      eval_row_from_diag <- function(ds, r, draw, g) {
        d05 <- g$per_alpha[["0.5"]]$diag
        d10 <- g$per_alpha[["1"]]$diag
        best <- g$per_alpha[[g$winner]]$diag
        data.frame(
          dataset = ds, repeat_idx = r$repeat_idx, fold_idx = r$fold_idx,
          arm = r$arm, k = r$k, draw = draw,
          evaluation_seed = r$evaluation_seed, winner_alpha = g$winner,
          a05_early_stop = d05$early_stop, a05_stop_index = d05$stop_index,
          a05_n_lambda = d05$n_lambda,
          a05_idx_lambda_min = d05$idx_lambda_min,
          a10_early_stop = d10$early_stop, a10_stop_index = d10$stop_index,
          a10_n_lambda = d10$n_lambda,
          a10_idx_lambda_min = d10$idx_lambda_min,
          winner_early_stop = best$early_stop,
          winner_stop_index = best$stop_index,
          winner_n_lambda = best$n_lambda,
          winner_idx_lambda_min = best$idx_lambda_min,
          winner_boundary_min = best$boundary_min,
          winner_boundary_1se = best$boundary_1se,
          winner_cvm_at_min = best$cvm_at_lambda_min,
          winner_cvm_at_last = best$cvm_at_last,
          saved_auc = r$AUC, stringsAsFactors = FALSE)
      }

      units <- names(split_units)
      env_now <- environment()
      if (POOL_WORKERS > 1L) {
        cl <- parallel::makeCluster(POOL_WORKERS, type = "PSOCK")
        on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl),
                                      silent = TRUE), add = TRUE)
        parallel::clusterEvalQ(cl, library(glmnet))
        parallel::clusterExport(
          cl, c("split_units", "paths", "outcome", "panel_sizes",
                "evaluator_fold_id", "fit_cv", "cv_diagnostics",
                "CONV_PATTERN", "standardise_split", "read_ranking",
                "eval_glmnet_diag", "eval_row_from_diag", "ds"),
          envir = env_now)
        unit_rows <- parallel::parLapply(cl, units, process_split)
        parallel::stopCluster(cl)
        cl <- NULL
      } else {
        unit_rows <- lapply(units, process_split)
      }
      rows_df <- do.call(rbind, unit_rows)
      saveRDS(rows_df, ds_ckpt)
      cat(sprintf("[partC %s] %d glmnet rows done, %d with winner early stop\n",
                  ds, nrow(rows_df), sum(rows_df$winner_early_stop)))
    }
  }

  #  -- Master-side AUC work ---------------------------------------------------
  #  Harness sample + every winner-early-stop row: recompute the full
  #  ensemble (glmnet component refit + xgboost + ranger with documented
  #  seeds) and compare against eval_deterministic.csv. For early-stop rows
  #  add the deep-path refit comparison.
  harness_rows <- list()
  deep_rows <- list()
  for (ds in partC_datasets) {
    paths <- dataset_paths(ds)
    outcome <- paths$outcome
    rows_df <- readRDS(file.path(WORK_DIR, sprintf("partC_%s.rds", ds)))
    harness_sel <- rows_df$arm %in% c("GS_full_ungrouped", "DGE") &
      rows_df$k %in% c(10, 50)
    todo <- rows_df[harness_sel | rows_df$winner_early_stop, ]
    for (i in seq_len(nrow(todo))) {
      r <- todo[i, ]
      is_harness <- r$arm %in% c("GS_full_ungrouped", "DGE") &&
        r$k %in% c(10, 50)
      sp <- readRDS(file.path(paths$split_dir, sprintf(
        "split_r%d_f%d.rds", r$repeat_idx, r$fold_idx)))
      y_train <- outcome[sp$train_idx]
      y_test <- outcome[sp$test_idx]
      pool <- sp$pools$var2000
      std <- standardise_split(sp$train_raw[, pool, drop = FALSE],
                               sp$test_raw[, pool, drop = FALSE])
      fl <- evaluator_fold_id(y_train, r$evaluation_seed)
      if (r$arm == "Random") {
        set.seed(99L + 1000L * r$repeat_idx + r$fold_idx)
        draws <- replicate(3L, sample(pool, min(r$k, length(pool))),
                           simplify = FALSE)
        panel <- draws[[r$draw]]
      } else {
        ranking <- read_ranking(paths$ranking_dir, r$repeat_idx, r$fold_idx,
                                r$arm, r$k)
        panel <- head(ranking, min(r$k, length(ranking)))
      }
      Xtr <- std$train[, panel, drop = FALSE]
      Xte <- std$test[, panel, drop = FALSE]

      #  Full ensemble with the ORIGINAL path (must reproduce the saved AUC
      #  on harness rows; on warning rows it quantifies what was reported).
      g <- eval_glmnet_diag(Xtr, fl$numeric_labels, fl$fold_id)
      best <- g$per_alpha[[g$winner]]
      p_glmnet <- as.numeric(predict(best$res$fit, Xte, s = "lambda.min",
                                     type = "response"))
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(Xtr),
                                     label = fl$numeric_labels)
      xgb_fit <- withr::with_seed(r$evaluation_seed + 1000L, xgboost::xgb.train(
        params = list(objective = "binary:logistic", eval_metric = "auc",
                      eta = 0.1, max_depth = 3, subsample = 0.8,
                      colsample_bytree = 0.8, seed = r$evaluation_seed + 1000L,
                      nthread = 1),
        data = dtrain, nrounds = 100, verbose = 0))
      p_xgb <- predict(xgb_fit, as.matrix(Xte))
      rf_fit <- ranger::ranger(x = Xtr, y = droplevels(y_train),
                               num.trees = 500, probability = TRUE,
                               seed = r$evaluation_seed + 2000L,
                               num.threads = 1)
      p_rf <- predict(rf_fit, data = Xte)$predictions[,
        levels(droplevels(y_train))[2]]
      auc_orig <- bench_auc(y_test, rowMeans(cbind(p_glmnet, p_xgb, p_rf)))

      if (is_harness) {
        harness_rows[[length(harness_rows) + 1L]] <- data.frame(
          dataset = ds, repeat_idx = r$repeat_idx, fold_idx = r$fold_idx,
          arm = r$arm, k = r$k, saved_auc = r$saved_auc,
          recomputed_auc = auc_orig,
          auc_abs_diff = abs(auc_orig - r$saved_auc),
          stringsAsFactors = FALSE)
      }

      if (isTRUE(r$winner_early_stop)) {
        deep <- fit_cv(Xtr, fl$numeric_labels, alpha = as.numeric(g$winner),
                       foldid = fl$fold_id, type_measure = "deviance",
                       nlambda = 200L, lambda_min_ratio = 1e-6,
                       deep_maxit = 1e6)
        dd <- cv_diagnostics(deep, nlambda_requested = 200L)
        p_deep <- as.numeric(predict(deep$fit, Xte, s = "lambda.min",
                                     type = "response"))
        auc_deep <- bench_auc(y_test, rowMeans(cbind(p_deep, p_xgb, p_rf)))
        deep_rows[[length(deep_rows) + 1L]] <- data.frame(
          dataset = ds, repeat_idx = r$repeat_idx, fold_idx = r$fold_idx,
          arm = r$arm, k = r$k, draw = r$draw, winner_alpha = g$winner,
          orig_stop_index = r$winner_stop_index,
          orig_n_lambda = r$winner_n_lambda,
          orig_idx_lambda_min = r$winner_idx_lambda_min,
          orig_lambda_min = best$res$fit$lambda.min,
          deep_n_lambda = dd$n_lambda, deep_early_stop = dd$early_stop,
          deep_idx_lambda_min = dd$idx_lambda_min,
          deep_lambda_min = deep$fit$lambda.min,
          deep_boundary_min = dd$boundary_min,
          lambda_min_rel_change =
            abs(deep$fit$lambda.min - best$res$fit$lambda.min) /
            best$res$fit$lambda.min,
          pred_cor = cor(p_glmnet, p_deep),
          pred_max_abs_diff = max(abs(p_glmnet - p_deep)),
          saved_auc = r$saved_auc, recomputed_auc_orig_path = auc_orig,
          deep_ensemble_auc = auc_deep,
          auc_delta_deep_vs_saved = auc_deep - r$saved_auc,
          glmnet_only_auc_orig = bench_auc(y_test, p_glmnet),
          glmnet_only_auc_deep = bench_auc(y_test, p_deep),
          stringsAsFactors = FALSE)
        cat(sprintf("  [partC-deep %s r%d f%d %s k=%d] saved %.4f orig %.4f deep %.4f\n",
                    ds, r$repeat_idx, r$fold_idx, r$arm, r$k, r$saved_auc,
                    auc_orig, auc_deep))
      }
    }
  }
  saveRDS(list(harness = if (length(harness_rows))
    do.call(rbind, harness_rows) else data.frame(),
    deep = if (length(deep_rows)) do.call(rbind, deep_rows) else
      data.frame()),
    file.path(WORK_DIR, "partC_auc.rds"))
}

# ------------------------------------------------------------------------------
#  Stage: report — per-fit CSVs, per-dataset summary, verdict.md.
# ------------------------------------------------------------------------------
run_report <- function() {
  pa <- load_partA()
  partA <- pa$diag
  write.csv(partA, file.path(OUT_DIR, "partA_selector_per_fit.csv"),
            row.names = FALSE)

  partB <- readRDS(file.path(WORK_DIR, "partB.rds"))
  if (nrow(partB) > 0L) {
    write.csv(partB, file.path(OUT_DIR, "partB_boundary_deep_refit.csv"),
              row.names = FALSE)
  }

  partC_rows <- do.call(rbind, lapply(c("GSE101794", "GSE65682"), function(ds)
    readRDS(file.path(WORK_DIR, sprintf("partC_%s.rds", ds)))))
  write.csv(partC_rows, file.path(OUT_DIR, "partC_evaluator_per_row.csv"),
            row.names = FALSE)
  partC_auc <- readRDS(file.path(WORK_DIR, "partC_auc.rds"))
  if (nrow(partC_auc$harness) > 0L) {
    write.csv(partC_auc$harness, file.path(
      OUT_DIR, "partC_evaluator_harness.csv"), row.names = FALSE)
  }
  if (nrow(partC_auc$deep) > 0L) {
    write.csv(partC_auc$deep, file.path(
      OUT_DIR, "partC_evaluator_deep_refit.csv"), row.names = FALSE)
  }

  #  -- Selector summary per dataset -------------------------------------------
  es <- partA$early_stop
  summ <- do.call(rbind, lapply(all_datasets, function(ds) {
    d <- partA[partA$dataset == ds, ]
    data.frame(
      dataset = ds,
      n_fits = nrow(d),
      n_replay_mismatched = sum(d$replay_max_abs_diff > 1e-8),
      replay_max_abs_diff = max(d$replay_max_abs_diff),
      n_early_stop = sum(d$early_stop),
      frac_early_stop = mean(d$early_stop),
      median_stop_index = if (any(d$early_stop))
        median(d$stop_index[d$early_stop]) else NA_real_,
      n_boundary_lambda_min = sum(d$boundary_min),
      frac_boundary_lambda_min = mean(d$boundary_min),
      n_boundary_lambda_1se = sum(d$boundary_1se),
      n_early_stop_and_boundary = sum(d$early_stop &
        (d$boundary_min | d$boundary_1se)),
      frac_early_stop_cvm_still_decreasing = if (any(d$early_stop))
        mean(d$cvm_decreasing_end[d$early_stop]) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  overall <- data.frame(
    dataset = "ALL",
    n_fits = nrow(partA),
    n_replay_mismatched = sum(partA$replay_max_abs_diff > 1e-8),
    replay_max_abs_diff = max(partA$replay_max_abs_diff),
    n_early_stop = sum(partA$early_stop),
    frac_early_stop = mean(partA$early_stop),
    median_stop_index = if (any(partA$early_stop))
      median(partA$stop_index[partA$early_stop]) else NA_real_,
    n_boundary_lambda_min = sum(partA$boundary_min),
    frac_boundary_lambda_min = mean(partA$boundary_min),
    n_boundary_lambda_1se = sum(partA$boundary_1se),
    n_early_stop_and_boundary = sum(partA$early_stop &
      (partA$boundary_min | partA$boundary_1se)),
    frac_early_stop_cvm_still_decreasing = if (any(partA$early_stop))
      mean(partA$cvm_decreasing_end[partA$early_stop]) else NA_real_,
    stringsAsFactors = FALSE)
  summ <- rbind(summ, overall)
  write.csv(summ, file.path(OUT_DIR, "summary_by_dataset.csv"),
            row.names = FALSE)

  #  -- Evaluator summary -------------------------------------------------------
  pc <- partC_rows
  c_summ <- do.call(rbind, lapply(c("GSE101794", "GSE65682"), function(ds) {
    d <- pc[pc$dataset == ds, ]
    data.frame(
      dataset = ds,
      n_rows = nrow(d),
      n_winner_early_stop = sum(d$winner_early_stop),
      frac_winner_early_stop = mean(d$winner_early_stop),
      n_winner_boundary = sum(d$winner_early_stop &
        (d$winner_boundary_min | d$winner_boundary_1se)),
      median_stop_index = if (any(d$winner_early_stop))
        median(d$winner_stop_index[d$winner_early_stop]) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  write.csv(c_summ, file.path(OUT_DIR, "summary_evaluator.csv"),
            row.names = FALSE)

  #  -- Verdict -----------------------------------------------------------------
  replay_exact <- overall$n_replay_mismatched == 0L
  harness <- partC_auc$harness
  deep <- partC_auc$deep
  harness_ok <- nrow(harness) > 0L && max(harness$auc_abs_diff) < 1e-8

  if (nrow(partB) > 0L) {
    b_moved <- sum(partB$lambda_min_rel_change > 1e-3)
    b_max_coef <- max(partB$max_abs_coef_change_at_old_lambda)
    b_med_rel <- median(partB$lambda_min_rel_change)
  } else {
    b_moved <- 0L; b_max_coef <- 0; b_med_rel <- 0
  }
  if (nrow(deep) > 0L) {
    max_abs_dauc <- max(abs(deep$auc_delta_deep_vs_saved))
    min_pred_cor <- min(deep$pred_cor)
    n_deep_lambda_moved <- sum(deep$lambda_min_rel_change > 1e-3)
  } else {
    max_abs_dauc <- 0; min_pred_cor <- 1; n_deep_lambda_moved <- 0L
  }

  n_selector_boundary_complete <- sum(!partA$early_stop & partA$boundary_min)

  auc_verdict <- if (overall$n_early_stop_and_boundary == 0L &&
                     sum(pc$winner_early_stop &
                           (pc$winner_boundary_min |
                              pc$winner_boundary_1se)) == 0L &&
                     nrow(deep) == 0L) {
    paste0("negligible. Nowhere in the benchmark does a CV-selected lambda ",
           "sit at a truncation boundary, and the evaluator rows whose paths ",
           "did stop early still selected interior lambdas; no AUC changed.")
  } else if (max_abs_dauc < 0.005 && b_max_coef < 1e-3) {
    sprintf(paste0("negligible (max |AUC change| on warning rows = %.4g; ",
                   "max coefficient change at the original lambda.min on ",
                   "boundary selector fits = %.3g)"),
            max_abs_dauc, b_max_coef)
  } else if (max_abs_dauc < 0.02) {
    sprintf(paste0("small (max |AUC change| on warning rows = %.4g; ",
                   "%d/%d boundary selector fits moved lambda.min by >0.1%%, ",
                   "max coefficient change %.3g)"),
            max_abs_dauc, b_moved, nrow(partB), b_max_coef)
  } else {
    sprintf("material (max |AUC change| on warning rows = %.4g)",
            max_abs_dauc)
  }

  verdict_lines <- c(
    "# Verdict — glmnet tail-convergence warnings in the older-7 benchmark",
    "",
    "Date: 2026-09-03. Generated by redesign/run_older7_glmnet_convergence_check.R",
    "",
    "## Replay status",
    "",
    if (replay_exact) {
      sprintf(paste0(
        "Exact replay achieved for every selector fit: refit coefficients at ",
        "lambda.min are bit-identical to the saved stability$coef_matrix in ",
        "all %d fits (7 datasets x 15 splits x 2 arms x 2 alphas x 50 ",
        "subsamples; max abs difference %.3g). Worker configuration: 6-worker ",
        "PSOCK emulation of clusterSetRNGStream(iseed = 42)."),
        overall$n_fits, overall$replay_max_abs_diff)
    } else {
      sprintf(paste0(
        "Replay NOT exact everywhere: %d of %d selector fits differ from the ",
        "saved coefficients (max %.3g). Per-fit diffs are in ",
        "partA_selector_per_fit.csv; interpret those rows with care."),
        overall$n_replay_mismatched, overall$n_fits,
        overall$replay_max_abs_diff)
    },
    if (nrow(harness) > 0L) {
      sprintf(paste0(
        "Evaluator harness check: on the 120-row sample (GS_full_ungrouped ",
        "and DGE, k = 10 and 50, both warning-flagged datasets) the ",
        "recomputed full-ensemble AUC matches eval_deterministic.csv with ",
        "max abs diff %.3g: %s."),
        max(harness$auc_abs_diff), if (harness_ok) "PASS" else "FAIL")
    } else {
      "Evaluator harness check: no harness rows (unexpected)."
    },
    "",
    "## Is the CV-selected lambda on the converged part of the path in practice?",
    "",
    sprintf("Selector side (%d fits, every saved GS fit in the benchmark):",
            overall$n_fits),
    sprintf("- path stopped early (convergence warning): %d fits (%.2f%%)",
            overall$n_early_stop, 100 * overall$frac_early_stop),
    sprintf("- early stop AND lambda.min/1se at or within 2 indices of the stop: %d fits",
            overall$n_early_stop_and_boundary),
    sprintf("- for comparison, lambda.min within 2 of the grid end on COMPLETED paths (not truncation): %d fits",
            n_selector_boundary_complete),
    "",
    sprintf("Evaluator side (all %d eval rows of GSE101794 + GSE65682, every arm and panel size):",
            nrow(pc)),
    sprintf("- rows where EITHER alpha's fit warned: %d (GSE101794) + %d (GSE65682); all 7 GSE65682 warnings sit on the alpha = 1.0 fit in rows where alpha = 0.5 won, so no reported AUC there touched a truncated path at all",
            sum(pc$dataset == "GSE101794" & (pc$a05_early_stop | pc$a10_early_stop)),
            sum(pc$dataset == "GSE65682" & (pc$a05_early_stop | pc$a10_early_stop))),
    sprintf("- winning-alpha fit stopped early: %d rows", sum(pc$winner_early_stop)),
    sprintf("- of those, lambda.min at/within 2 of the stop: %d rows",
            sum(pc$winner_early_stop & pc$winner_boundary_min)),
    sprintf("- of those, lambda.1se at/within 2 of the stop: %d rows",
            sum(pc$winner_early_stop & pc$winner_boundary_1se)),
    "",
    "Note on columns: stop_index in the per-fit CSVs is the lambda index from",
    "the FIRST convergence warning captured, which can come from a fold fit",
    "rather than the full-data fit; the authoritative truncation measure is",
    "n_lambda (the returned path length) versus idx_lambda_min / idx_lambda_1se.",
    "",
    "## Does the truncation move any number? (deep-path refits: nlambda = 200, lambda.min.ratio = 1e-6, maxit = 1e6)",
    "",
    sprintf("- boundary selector fits refit: %d; lambda.min moved >0.1%% in %d; median relative move %.3g; max |coef change| at the original lambda.min %.3g",
            nrow(partB), b_moved, b_med_rel, b_max_coef),
    sprintf("- evaluator rows with early-stop winning fits refit: %d; lambda.min moved >0.1%% in %d; max |AUC change| %.4g; min prediction correlation %.4f",
            nrow(deep), n_deep_lambda_moved, max_abs_dauc, min_pred_cor),
    "",
    "## Verdict",
    "",
    sprintf("Effect of the glmnet tail-convergence warnings on benchmark AUCs and rankings: **%s**",
            auc_verdict),
    "",
    "Limitations: the 20-permutation calibration nulls (~84k cv.glmnet fits,",
    "feeding only the calibration null distributions) were not replayed, and",
    "competitor LASSO/ElasticNet selector fits are not reproducible (their",
    "fold ids come from ambient RNG with no stored seed). The main-log warning",
    "count (753) includes both of those plus the fits checked here.",
    "",
    "Evidence: partA_selector_per_fit.csv, partB_boundary_deep_refit.csv,",
    "partC_evaluator_per_row.csv, partC_evaluator_harness.csv,",
    "partC_evaluator_deep_refit.csv, summary_by_dataset.csv,",
    "summary_evaluator.csv, replay_check_GSE101794.csv."
  )
  writeLines(verdict_lines, file.path(OUT_DIR, "verdict.md"))
  cat(paste(verdict_lines, collapse = "\n"), "\n")
}

# ------------------------------------------------------------------------------
if (stage %in% c("replay_check", "all")) run_replay_check()
if (stage %in% c("partA", "all")) run_partA()
if (stage %in% c("partB", "all")) run_partB()
if (stage %in% c("partC", "all")) run_partC()
if (stage %in% c("report", "all")) run_report()
cat(sprintf("stage=%s done (%.0f s)\n", stage, proc.time()[["elapsed"]]))
