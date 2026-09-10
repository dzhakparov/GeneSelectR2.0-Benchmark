# =============================================================================
# GeneSelectR Slim - self-contained driver
# =============================================================================
#
# The complete slim pipeline: repeated stratified K-fold elastic net, selection
# frequency pi, conditional OOB deviance utility u (utility_deviance.R), and
# the joint permutation-null combination (joint_null.R).
#
# This file deliberately does NOT depend on the GeneSelectR package. The slim
# redesign is being evaluated against the package's own variants; sharing code
# with the incumbent would make a bug in the incumbent invisible in the
# comparison. Everything here is small enough to audit at a glance.
#
# REPRODUCIBILITY CONTRACT
# ------------------------
# Every stochastic step is seeded EXPLICITLY and LOCALLY:
#   * subsample construction is seeded once per call;
#   * each cv.glmnet fit receives its own deterministic seed derived from
#     (random_seed, alpha index, subsample index);
#   * permutation m is seeded by random_seed + 7*m, and its fits live in a
#     disjoint seed range.
# Because no fit depends on ambient RNG state, serial and parallel execution
# produce IDENTICAL results at the same seed. There is a test for this.
#
# Seed ranges (all < 2^31 for the defaults):
#   real fits, alpha index a, subsample b : random_seed + a*1e4 + b
#   permutation m outcome shuffle         : random_seed + 7*m
#   null fits, permutation m, subsample b : random_seed + 1e6 + m*1e3 + b


# -----------------------------------------------------------------------------
#  Subsampling: stratified repeated K-fold
# -----------------------------------------------------------------------------
#
# Same construction as the incumbent ("kfold" scheme): B = repeats x k_folds
# training sets, each seeing (K-1)/K of the data. Kept identical so GS_slim
# differs from the incumbents ONLY in scoring, not in resampling.

slim_subsamples <- function(y, B = 50, k_folds = 5, random_seed = 123) {
  if (!is.factor(y) || nlevels(y) != 2L || anyNA(y)) {
    stop("y must be a non-missing factor with exactly two levels")
  }
  set.seed(random_seed)

  n <- length(y)
  class_indices <- split(seq_len(n), y)
  n_repeats <- ceiling(B / k_folds)
  subsamples <- vector("list", 0)

  for (repeat_idx in seq_len(n_repeats)) {
    fold_id <- integer(n)
    for (cls in names(class_indices)) {
      idx <- class_indices[[cls]]
      fold_id[idx] <- sample(rep(seq_len(k_folds), length.out = length(idx)))
    }
    for (fold in seq_len(k_folds)) {
      if (length(subsamples) >= B) break
      oob_idx <- which(fold_id == fold)
      subsamples[[length(subsamples) + 1]] <- list(
        train = setdiff(seq_len(n), oob_idx),
        oob = oob_idx
      )
    }
  }
  subsamples
}


# -----------------------------------------------------------------------------
#  AUC (rank-based, fixed direction)
# -----------------------------------------------------------------------------
#
# AUC with the second factor level as the positive class, computed from rank
# sums. No direction = "auto": a below-chance model must return AUC < 0.5,
# because the AUC-weighting in the utility relies on it (weight 0 for
# non-predictive subsample models).

slim_auc <- function(y_true, scores) {
  y_true <- droplevels(y_true)
  pos <- which(y_true == levels(y_true)[2])
  neg <- which(y_true == levels(y_true)[1])
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(scores, ties.method = "average")
  (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}


# -----------------------------------------------------------------------------
#  One subsample fit
# -----------------------------------------------------------------------------
#
# cv.glmnet at lambda.min, returning the FULL coefficient vector, the
# intercept, and the OOB AUC. The intercept is what the incumbent threw away
# with predict_fn; the deviance utility needs it to reconstruct eta exactly.
# The fit is seeded locally so it does not depend on ambient RNG state.

slim_fit_one <- function(X, y, train_idx, oob_idx, alpha, cv_folds = 5,
                         seed = 1) {
  y01 <- as.numeric(y[train_idx]) - 1
  #  glmnet changes AUC to deviance when a validation fold has fewer than ten
  #  observations. Select that objective explicitly so the recorded fit
  #  settings describe the computation that was performed.
  type_measure <- if (floor(length(train_idx) / cv_folds) < 10) {
    "deviance"
  } else {
    "auc"
  }

  set.seed(seed)
  fit <- tryCatch(
    glmnet::cv.glmnet(X[train_idx, , drop = FALSE], y01,
                      family = "binomial", alpha = alpha,
                      nfolds = cv_folds, type.measure = type_measure),
    error = function(e) {
      stop(sprintf("slim subsample glmnet fit failed: %s",
                   conditionMessage(e)), call. = FALSE)
    }
  )

  coefs <- glmnet::coef.glmnet(fit, s = "lambda.min")
  beta <- as.numeric(coefs)[-1]
  intercept <- as.numeric(coefs)[1]

  eta_oob <- as.numeric(X[oob_idx, , drop = FALSE] %*% beta) + intercept
  auc <- slim_auc(y[oob_idx], plogis(eta_oob))

  list(coef = beta, intercept = intercept, auc = auc)
}


# Fit all subsamples for one (X, y, alpha). fit_seeds holds one deterministic
# seed per subsample, so serial and parallel execution agree exactly.
slim_fit_all <- function(X, y, subsamples, alpha, fit_seeds,
                         n_cores = 1, cv_folds = 5) {
  if (length(fit_seeds) != length(subsamples)) {
    stop("fit_seeds must have one entry per subsample")
  }

  results <- if (n_cores > 1) {
    cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cluster))
    parallel::clusterEvalQ(cluster, library(glmnet))
    parallel::clusterExport(
      cluster, c("slim_fit_one", "slim_auc", "X", "y", "subsamples",
                 "alpha", "fit_seeds", "cv_folds"),
      envir = environment()
    )
    parallel::parLapply(cluster, seq_along(subsamples), function(i) {
      slim_fit_one(X, y, subsamples[[i]]$train, subsamples[[i]]$oob,
                   alpha, cv_folds = cv_folds, seed = fit_seeds[i])
    })
  } else {
    lapply(seq_along(subsamples), function(i) {
      slim_fit_one(X, y, subsamples[[i]]$train, subsamples[[i]]$oob,
                   alpha, cv_folds = cv_folds, seed = fit_seeds[i])
    })
  }

  p <- ncol(X)
  B <- length(subsamples)
  coef_matrix <- matrix(0, nrow = p, ncol = B)
  intercepts <- rep(NA_real_, B)
  auc_vec <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    coef_matrix[, b] <- results[[b]]$coef
    intercepts[b] <- results[[b]]$intercept
    auc_vec[b] <- results[[b]]$auc
  }
  colnames(coef_matrix) <- paste0("b", seq_len(B))

  list(coef_matrix = coef_matrix, intercepts = intercepts, auc_vec = auc_vec)
}


# -----------------------------------------------------------------------------
#  Permutation null (paired pi and u)
# -----------------------------------------------------------------------------
#
# M permutations of y; X is never touched, so the gene-gene correlation
# structure -- and with it whatever dependence pi and u actually have -- is
# preserved exactly. Each permutation yields one row of null_pi and one row of
# null_u from the SAME fits: the pairing is what makes the joint null valid.
#
# The utility for each permutation is computed with the same construction as
# the observed fit: candidates = genes selected at least once in that
# permutation, AUC-weighted with the permutation's own OOB AUCs.

slim_null_one_permutation <- function(X, y, m, B_null, k_folds, alpha,
                                      auc_weight, random_seed) {
  p <- ncol(X)

  set.seed(random_seed + 7 * m)
  y_perm <- sample(y)

  subsamples <- slim_subsamples(y_perm, B = B_null, k_folds = k_folds,
                                random_seed = random_seed + 11 * m)
  fit_seeds <- random_seed + 1e6 + m * 1e3 + seq_len(B_null)
  fits <- slim_fit_all(X, y_perm, subsamples, alpha, fit_seeds, n_cores = 1)

  null_pi <- rowMeans(fits$coef_matrix != 0)

  candidates <- which(null_pi > 0)
  null_u <- rep(0, p)
  if (length(candidates) > 0) {
    util <- compute_deviance_utility(
      X, y_perm, subsamples, fits$coef_matrix, fits$intercepts,
      candidate_idx = candidates,
      auc_vec = fits$auc_vec, auc_weight = auc_weight, verbose = FALSE
    )
    null_u <- util$u
  }
  list(pi = null_pi, u = null_u)
}


slim_permutation_null <- function(X, y, M = 100, B_null = 20, k_folds = 5,
                                  alpha = 0.5, auc_weight = TRUE,
                                  n_cores = 1, random_seed = 123,
                                  verbose = TRUE) {
  run_one <- function(m) {
    if (verbose && m %% 10 == 0) cat(sprintf("    null %d/%d\n", m, M))
    slim_null_one_permutation(X, y, m, B_null, k_folds, alpha,
                              auc_weight, random_seed)
  }

  null_list <- if (n_cores > 1) {
    cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cluster))
    parallel::clusterEvalQ(cluster, library(glmnet))
    parallel::clusterExport(
      cluster,
      c("slim_null_one_permutation", "slim_subsamples", "slim_fit_all",
        "slim_fit_one", "slim_auc", "compute_deviance_utility",
        ".binomial_deviance",
        "X", "y", "B_null", "k_folds", "alpha", "auc_weight", "random_seed"),
      envir = environment()
    )
    parallel::parLapply(cluster, seq_len(M), run_one)
  } else {
    lapply(seq_len(M), run_one)
  }

  null_pi <- do.call(rbind, lapply(null_list, `[[`, "pi"))
  null_u  <- do.call(rbind, lapply(null_list, `[[`, "u"))
  colnames(null_pi) <- colnames(X)
  colnames(null_u)  <- colnames(X)

  list(null_pi = null_pi, null_u = null_u)
}


# -----------------------------------------------------------------------------
#  The driver
# -----------------------------------------------------------------------------

#' GeneSelectR Slim
#'
#' @param X Numeric matrix, samples x genes, named columns
#' @param y Two-level factor; the second level is the modelled class
#' @param alpha_grid Elastic-net mixing parameters tried per fit; the winner is
#'   chosen by mean OOB AUC (the incumbent's internal tuning rule, unchanged)
#' @param B Subsamples per alpha
#' @param k_folds Folds for the subsampling scheme
#' @param M_null Permutations for the joint null
#' @param B_null Subsamples per null permutation
#' @param auc_weight Weight subsample utility contributions by OOB AUC - 0.5
#' @param n_cores PSOCK workers for the fitting loops
#' @param random_seed Seed controlling everything
#' @return List with ranked (gene names, best first), gene_table (one row per
#'   gene: pi_raw, u_dev, e_pi, e_u, p_combined, sign_consistency,
#'   n_selected_b), best_alpha, mean_oob_auc per alpha, and the null matrices
#'   (kept for diagnostics)
gs_slim_fit <- function(X, y,
                        alpha_grid = c(0.5, 1.0),
                        B = 50, k_folds = 5,
                        M_null = 100, B_null = 20,
                        auc_weight = TRUE,
                        n_cores = 1,
                        random_seed = 123,
                        verbose = TRUE) {

  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix")
  if (is.null(colnames(X))) stop("X must have gene-name columns")
  if (!is.factor(y) || nlevels(y) != 2L) stop("y must be a two-level factor")
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y)")
  if (any(!is.finite(X))) stop("X contains non-finite values")

  p <- ncol(X)
  subsamples <- slim_subsamples(y, B = B, k_folds = k_folds,
                                random_seed = random_seed)

  # --- Real fits per alpha; alpha chosen by mean OOB AUC --------------------
  fits_by_alpha <- list()
  mean_auc_by_alpha <- numeric(length(alpha_grid))
  for (a_idx in seq_along(alpha_grid)) {
    alpha <- alpha_grid[a_idx]
    if (verbose) cat(sprintf("  alpha = %.2f: %d subsample fits\n",
                             alpha, length(subsamples)))
    fit_seeds <- random_seed + a_idx * 1e4 + seq_len(length(subsamples))
    fits_by_alpha[[a_idx]] <- slim_fit_all(
      X, y, subsamples, alpha, fit_seeds, n_cores = n_cores
    )
    if (any(!is.finite(fits_by_alpha[[a_idx]]$auc_vec))) {
      stop(sprintf("Non-finite OOB AUC for alpha %.2f", alpha),
           call. = FALSE)
    }
    mean_auc_by_alpha[a_idx] <- mean(fits_by_alpha[[a_idx]]$auc_vec)
    if (verbose) {
      cat(sprintf("    mean OOB AUC: %.4f\n", mean_auc_by_alpha[a_idx]))
    }
  }
  best_idx <- which.max(mean_auc_by_alpha)
  best_alpha <- alpha_grid[best_idx]
  fits <- fits_by_alpha[[best_idx]]
  if (verbose) cat(sprintf("  best alpha: %.2f\n", best_alpha))

  # --- Stability -------------------------------------------------------------
  pi_raw <- rowMeans(fits$coef_matrix != 0)
  names(pi_raw) <- colnames(X)

  # --- Utility ----------------------------------------------------------------
  candidates <- which(pi_raw > 0)
  if (length(candidates) == 0) {
    stop("No gene was selected in any subsample; nothing to score.")
  }
  util <- compute_deviance_utility(
    X, y, subsamples, fits$coef_matrix, fits$intercepts,
    candidate_idx = candidates,
    auc_vec = fits$auc_vec, auc_weight = auc_weight, verbose = verbose
  )

  # --- Joint null --------------------------------------------------------------
  if (verbose) {
    cat(sprintf("  joint null: %d permutations x %d subsamples (alpha = %.2f)\n",
                M_null, B_null, best_alpha))
  }
  nulls <- slim_permutation_null(
    X, y, M = M_null, B_null = B_null, k_folds = k_folds,
    alpha = best_alpha, auc_weight = auc_weight,
    n_cores = n_cores, random_seed = random_seed, verbose = verbose
  )

  # --- Evidence combination ----------------------------------------------------
  evidence <- combine_evidence_joint(pi_raw, util$u,
                                     nulls$null_pi, nulls$null_u)

  gene_table <- data.frame(
    gene = colnames(X),
    pi_raw = pi_raw,
    n_selected_b = util$n_selected_b,
    sign_consistency = util$sign_consistency,
    u_dev = util$u,
    e_pi = evidence$e_pi,
    e_u = evidence$e_u,
    p_combined = evidence$p_combined,
    stringsAsFactors = FALSE
  )
  gene_table <- gene_table[order(gene_table$p_combined, -gene_table$pi_raw), ]
  rownames(gene_table) <- NULL

  list(
    ranked = gene_table$gene,
    gene_table = gene_table,
    best_alpha = best_alpha,
    mean_oob_auc = stats::setNames(mean_auc_by_alpha,
                                   paste0("alpha_", alpha_grid)),
    null_pi = nulls$null_pi,
    null_u = nulls$null_u,
    subsample_auc = fits$auc_vec
  )
}
