# =============================================================================
# GeneSelectR Slim - Conditional OOB Deviance Utility (reference implementation)
# =============================================================================
#
# WHAT THIS REPLACES
# ------------------
# The current utility pillar is geomean(instance_SHAP_frequency, MI). Two of
# its three ingredients are broken at the operating point of this package
# (n ~ 120-150, p = 2000):
#
#   * MI on ~120 samples with 3-5 quantile bins is mostly finite-sample bias,
#     and the per-subsample estimates are correlated under K-fold overlap, so
#     averaging them does not de-noise them.
#   * instance_SHAP_frequency divides by OOB APPEARANCES, including subsamples
#     where the gene was not selected (exact zero SHAP). The score is therefore
#     selection frequency x conditional magnitude -- the stability pillar
#     counted a second time, which is why the geometric mean double-counts and
#     why GS_no_stability / GS_harmonic win the benchmark.
#
# WHAT THIS IS
# ------------
# One held-out number per gene: the mean increase in OOB binomial deviance
# when the gene's coefficient is zeroed, averaged ONLY over subsamples that
# selected the gene. For a linear predictor this is exact and refit-free:
# eta^{-j} = eta - beta_j * x_j. pi owns "how often chosen", u owns "how much
# it helps when chosen" -- separable by construction, not just in expectation.
#
# Subsample contributions are weighted by max(OOB AUC - 0.5, 0): a gene's
# contribution inside a model that itself predicts nothing is not evidence.


#' Binomial Deviance of a Linear Predictor (internal)
#'
#' D(eta) = -2 * sum_i [ y_i * eta_i - softplus(eta_i) ]
#' softplus computed stably as log1p(exp(-|eta|)) + max(eta, 0).
#'
#' @keywords internal
.binomial_deviance <- function(eta, y01) {
  softplus <- log1p(exp(-abs(eta))) + pmax(eta, 0)
  -2 * sum(y01 * eta - softplus)
}


#' Conditional OOB Deviance Utility
#'
#' @param X Full expression matrix (n x p)
#' @param y Two-level factor outcome; the second level is the modelled class
#' @param subsamples List of subsample objects (train/oob indices)
#' @param coef_matrix p x B matrix of per-subsample coefficient vectors
#' @param intercepts Length-B numeric vector of per-subsample intercepts.
#'   Requires the Step-1 loop to retain fit$intercept (see spec section 6).
#' @param candidate_idx Column indices to score; genes never selected are not
#'   candidates because their conditional mean has no terms.
#' @param auc_vec Length-B vector of per-subsample OOB AUCs (may contain NA)
#' @param auc_weight Logical; weight each subsample's contribution by
#'   max(AUC_b - 0.5, 0). If every selecting subsample has weight 0, the gene
#'   gets u = 0: it was selected by models that do not predict.
#' @param verbose Print progress
#' @return List with u (length-p utility vector), n_selected_b (how many
#'   subsamples selected each gene), sign_consistency (fraction of selecting
#'   subsamples carrying the majority sign), and delta_matrix (n_candidates x B
#'   per-subsample contributions, NA where the gene was not selected -- kept
#'   for diagnostics and for the null construction, which must reuse the exact
#'   same statistic)
#' @keywords internal
compute_deviance_utility <- function(X, y, subsamples, coef_matrix, intercepts,
                                     candidate_idx,
                                     auc_vec = NULL,
                                     auc_weight = TRUE,
                                     verbose = TRUE) {

  n_genes <- ncol(X)
  B <- length(subsamples)
  y01 <- as.numeric(y) - 1   # second factor level = 1, matching the fit

  if (length(intercepts) != B) {
    stop("intercepts must have one entry per subsample")
  }
  if (is.null(auc_vec)) {
    auc_vec <- rep(1, B)     # uniform weights when AUCs are unavailable
  }
  if (auc_weight) {
    w_b <- pmax(ifelse(is.finite(auc_vec), auc_vec, 0.5) - 0.5, 0)
  } else {
    w_b <- rep(1, B)
  }

  n_cand <- length(candidate_idx)
  delta_matrix <- matrix(NA_real_, nrow = n_cand, ncol = B)
  colnames(delta_matrix) <- paste0("b", seq_len(B))
  rownames(delta_matrix) <- colnames(X)[candidate_idx]

  for (b in seq_len(B)) {
    beta <- coef_matrix[, b]
    a <- intercepts[b]
    active_cand <- which(beta[candidate_idx] != 0)
    if (length(active_cand) == 0 || !is.finite(a)) next

    oob_idx <- subsamples[[b]]$oob
    X_oob <- X[oob_idx, , drop = FALSE]
    y_oob <- y01[oob_idx]

    eta <- as.numeric(X_oob %*% beta) + a
    dev_full <- .binomial_deviance(eta, y_oob)

    # For all active candidates at once: eta^{-j} = eta - beta_j * x_j.
    # n_oob x n_active matrix, so the per-column deviance is vectorised.
    cols <- candidate_idx[active_cand]
    E <- matrix(eta, nrow = length(oob_idx), ncol = length(cols)) -
      sweep(X_oob[, cols, drop = FALSE], 2, beta[cols], "*")

    softplus <- log1p(exp(-abs(E))) + pmax(E, 0)
    dev_minus <- -2 * (colSums(y_oob * E) - colSums(softplus))

    delta_matrix[active_cand, b] <- dev_minus - dev_full

    if (verbose && b %% 25 == 0) {
      cat(sprintf("    deviance utility: %d/%d subsamples\n", b, B))
    }
  }

  # Conditional, weighted mean. A gene selected by zero-weight subsamples
  # only gets NA in the numerator AND denominator -> u = 0 below.
  w_mat <- matrix(w_b, nrow = n_cand, ncol = B, byrow = TRUE)
  w_mat[is.na(delta_matrix)] <- 0
  delta_zeroed <- delta_matrix
  delta_zeroed[is.na(delta_zeroed)] <- 0

  weight_sum <- rowSums(w_mat)
  u_cand <- ifelse(weight_sum > 0,
                   rowSums(delta_zeroed * w_mat) / weight_sum,
                   0)

  # Negative deltas (the gene HURT OOB prediction) are real evidence against;
  # clip at zero so "selected but harmful" ties "never useful" rather than
  # going below it. The joint null handles the asymmetry: under permuted y,
  # harmful genes get the same clip, so the comparison stays matched.
  u_cand <- pmax(u_cand, 0)

  u <- rep(0, n_genes)
  u[candidate_idx] <- u_cand
  names(u) <- colnames(X)

  n_selected_b <- rowSums(!is.na(delta_matrix))
  n_sel_full <- rep(0L, n_genes)
  n_sel_full[candidate_idx] <- n_selected_b
  names(n_sel_full) <- colnames(X)

  # Sign consistency over selecting subsamples. Not a score component -- a
  # diagnostic for correlated-cluster stand-ins, whose sign flips as the
  # elastic net re-picks cluster members.
  sign_consistency <- rep(NA_real_, n_genes)
  names(sign_consistency) <- colnames(X)
  for (k in seq_len(n_cand)) {
    j <- candidate_idx[k]
    signs <- sign(coef_matrix[j, ])
    signs <- signs[signs != 0]
    if (length(signs) > 0) {
      sign_consistency[j] <- max(mean(signs > 0), mean(signs < 0))
    }
  }

  list(
    u = u,
    n_selected_b = n_sel_full,
    sign_consistency = sign_consistency,
    delta_matrix = delta_matrix,
    weights = w_b
  )
}
