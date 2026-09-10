# =============================================================================
# GeneSelectR Slim - Joint-Null Evidence Combination (reference implementation)
# =============================================================================
#
# WHAT THIS REPLACES
# ------------------
# combine_scores() offered geometric / arithmetic / harmonic / minimum over
# three pillars with fixed weights. All four are fixed-formula aggregations of
# DEPENDENT inputs; the benchmark showed the winner is decided by which veto
# you remove, not by which formula is principled. The evidence-ratio
# calibration fixed the zero-point but kept two weaknesses: it divides by the
# MEAN of the null (one large permutation collapses a gene's ratio), and the
# likelihood-ratio reading of the geometric mean assumes pillar independence,
# which the pillars do not have.
#
# WHAT THIS IS
# ------------
# The permutation null produces PAIRED (null_pi, null_u) per permutation --
# same fits, same permuted y. Whatever dependence pi and u actually have is
# therefore present in the null. Combining against that joint null handles the
# dependence exactly, with no formula and no weights:
#
#   1. Per-gene exceedance of the observed statistics against their own nulls
#      (leave-one-out, add-one smoothing).
#   2. Fisher combination into T_j.
#   3. The SAME T computed for every permutation, treating each permutation as
#      the observation against the other M-1.
#   4. p_j = fraction of null T's at least as extreme as the observed T_j.
#
# There are no tuning constants in this file. The only choice is M, and it is
# a resolution choice (p-granularity 1/(M+1)), not a modelling choice.


#' Per-Gene Exceedance Against a Null Matrix
#'
#' e_j = (1 + #{m : null[m, j] >= obs_j}) / (M + 1)
#'
#' The add-one smoothing keeps e in (0, 1] so logs are finite, and is the
#' standard permutation-p-value correction (the observation itself is one
#' draw under the null hypothesis).
#'
#' @keywords internal
.exceedance <- function(obs, null_matrix) {
  M <- nrow(null_matrix)
  (1 + colSums(sweep(null_matrix, 2, obs, ">="), na.rm = TRUE)) / (M + 1)
}


#' Leave-One-Out Null Exceedances
#'
#' For each permutation m, the exceedance of its OWN statistics against the
#' other M-1 permutations. Needed to put the null replicates through the same
#' pipeline as the observation, so the null distribution of the combined
#' statistic is comparable to the observed one. Leave-one-out is required:
#' including m in its own null would bias every null exceedance downward and
#' make the combined null anti-conservative.
#'
#' @param null_matrix M x p matrix
#' @return M x p matrix of leave-one-out exceedances
#' @keywords internal
.loo_exceedance <- function(null_matrix) {
  M <- nrow(null_matrix)
  p <- ncol(null_matrix)
  out <- matrix(NA_real_, nrow = M, ncol = p)
  for (m in seq_len(M)) {
    others <- null_matrix[-m, , drop = FALSE]
    # (1 + #{others >= x}) / M keeps the same add-one scale as .exceedance,
    # where the denominator is M+1 for M null draws: here there are M-1
    # others, so the denominator is (M-1)+1 = M.
    out[m, ] <- (1 + colSums(sweep(others, 2, null_matrix[m, ], ">="),
                             na.rm = TRUE)) / M
  }
  out
}


#' Combine Stability and Utility Evidence Against the Joint Null
#'
#' @param pi_obs Length-p observed selection frequencies
#' @param u_obs Length-p observed deviance utilities
#' @param null_pi M x p null selection frequencies
#' @param null_u M x p null deviance utilities, PAIRED with null_pi by
#'   permutation row (same permuted fits)
#' @return List with per-gene exceedances (e_pi, e_u), the observed Fisher
#'   statistic T_obs, the joint-null p-value p_combined, and the null matrix
#'   of combined statistics T_null (M x p, for diagnostics)
#' @keywords internal
combine_evidence_joint <- function(pi_obs, u_obs, null_pi, null_u) {

  if (nrow(null_pi) != nrow(null_u) || ncol(null_pi) != ncol(null_u)) {
    stop("null_pi and null_u must have identical dimensions (paired rows)")
  }
  if (ncol(null_pi) != length(pi_obs) || ncol(null_u) != length(u_obs)) {
    stop("null matrices must have one column per gene")
  }
  M <- nrow(null_pi)
  if (M < 20) {
    warning(sprintf(
      paste0("M = %d permutations gives p-resolution ~%.3f. ",
             "Use M >= 100 for a ranking that separates the top of the list."),
      M, 1 / (M + 1)
    ))
  }

  # --- Observed ---
  e_pi <- .exceedance(pi_obs, null_pi)
  e_u  <- .exceedance(u_obs,  null_u)
  T_obs <- -2 * (log(e_pi) + log(e_u))

  # --- Null replicates through the identical pipeline ---
  e_pi_null <- .loo_exceedance(null_pi)
  e_u_null  <- .loo_exceedance(null_u)
  T_null <- -2 * (log(e_pi_null) + log(e_u_null))

  # --- Joint p-value ---
  # Per gene, the fraction of null combined statistics at least as extreme.
  # Add-one as in .exceedance: the observation is one draw under H0.
  p_combined <- (1 + colSums(sweep(T_null, 2, T_obs, ">="), na.rm = TRUE)) /
    (M + 1)

  names(p_combined) <- names(pi_obs)

  list(
    e_pi = e_pi,
    e_u = e_u,
    T_obs = T_obs,
    T_null = T_null,
    p_combined = p_combined
  )
}


#' Reportable Rank Uncertainty (optional, cheap)
#'
#' The panel-aggregation critique needs a RULE, not just a ranking. This is
#' the rule: refit the whole pipeline R times with different seeds and report
#' each gene's probability of landing in the top k. The panel is genes with
#' P(top-k) >= tau_membership. Deliberately a thin wrapper -- the point is
#' that the rule is pre-registrable in one sentence.
#'
#' @param fits List of R result objects, each with a p_combined vector
#'   (named by gene)
#' @param k Panel size
#' @return Named numeric vector: P(in top k) per gene
#' @keywords internal
panel_membership_probability <- function(fits, k) {
  genes <- names(fits[[1]]$p_combined)
  membership <- matrix(FALSE, nrow = length(fits), ncol = length(genes),
                       dimnames = list(NULL, genes))
  for (r in seq_along(fits)) {
    p <- fits[[r]]$p_combined[genes]
    top <- order(p, decreasing = FALSE)[seq_len(min(k, length(genes)))]
    membership[r, top] <- TRUE
  }
  colMeans(membership)
}
