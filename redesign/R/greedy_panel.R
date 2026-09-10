# =============================================================================
# Greedy forward panel selection
# =============================================================================
#
# Every method in this project scores genes ONE AT A TIME and hopes the top-k
# works as a team. The benchmark judges the team. This file picks the team
# directly: start with the best single gene, then repeatedly add whichever
# candidate most improves the panel's cross-validated AUC, and keep the full
# path so any top-k panel is a prefix.
#
# Methodological notes, stated because they are choices:
#
#   * The score is inner stratified K-fold CV AUC of a logistic regression on
#     the current panel + candidate. The inner folds are FIXED for the whole
#     greedy run, so every step compares candidates on identical splits --
#     without that, fold noise looks like improvement.
#   * Candidates are pre-filtered to the top n_candidates by univariate
#     |AUC - 0.5|. Without the filter, step 1 alone costs p inner-CV rounds.
#   * The scorer is plain logistic regression (stats::glm). It is fast at
#     k <= 50 and n ~ 120, and -- deliberately -- it is NOT one of the three
#     evaluator models, so greedy cannot win by sharing assumptions with the
#     benchmark's ensemble.
#   * Everything is deterministic: glm has no RNG, folds are fixed, ties go
#     to the earlier candidate in univariate order. n_cores only changes
#     wall-clock time.


#' Rank genes by univariate |AUC - 0.5|
#'
#' @param X Samples x genes matrix
#' @param y Two-level factor
#' @return Character vector of gene names, most discriminative first
univariate_auc_rank <- function(X, y) {
  scores <- apply(X, 2, function(x) {
    a <- slim_auc(y, x)
    if (is.na(a)) 0 else abs(a - 0.5)
  })
  names(sort(scores, decreasing = TRUE))
}


#' Inner-CV AUC of a logistic panel model
#'
#' @param X Training matrix (training samples only)
#' @param y Outcome factor
#' @param panel Character vector of genes in the panel
#' @param inner_folds List of test-index vectors (fixed before the greedy run)
#' @return Mean AUC across inner test folds
panel_cv_auc <- function(X, y, panel, inner_folds) {
  aucs <- vapply(inner_folds, function(test_idx) {
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    fit <- tryCatch(
      suppressWarnings(
        stats::glm(y ~ ., data = data.frame(y = y[train_idx],
                                            X[train_idx, panel, drop = FALSE]),
                   family = stats::binomial)
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    probs <- tryCatch(
      suppressWarnings(
        stats::predict(fit, newdata = data.frame(
          X[test_idx, panel, drop = FALSE]), type = "response")
      ),
      error = function(e) NULL
    )
    if (is.null(probs)) return(NA_real_)
    slim_auc(y[test_idx], probs)
  }, numeric(1))
  if (any(!is.finite(aucs))) {
    stop("Greedy panel evaluation failed in an inner fold.", call. = FALSE)
  }
  mean(aucs)
}


#' Greedy forward panel selection
#'
#' @param X Samples x genes matrix (TRAINING data only)
#' @param y Two-level factor outcome
#' @param n_candidates Size of the univariate pre-filter pool
#' @param k_max Path length; every prefix is a valid top-k panel
#' @param inner_k Folds for the inner CV scorer
#' @param n_cores PSOCK workers for candidate evaluation (results identical)
#' @param random_seed Seeds the inner fold assignment only
#' @return List with path (ordered gene vector, length k_max), pool (the
#'   candidate pool, univariate order), step_auc (inner-CV AUC after each
#'   addition), and inner_folds (for audit)
greedy_forward_panel <- function(X, y, n_candidates = 200, k_max = 50,
                                 inner_k = 5, n_cores = 1,
                                 random_seed = 123) {
  if (k_max > n_candidates) stop("k_max cannot exceed n_candidates")
  if (nrow(X) != length(y)) stop("nrow(X) must equal length(y)")

  pool <- head(univariate_auc_rank(X, y), n_candidates)

  # Fixed inner folds: every greedy step scores on identical splits.
  inner_folds <- (function() {
    set.seed(random_seed)
    class_indices <- split(seq_len(nrow(X)), y)
    fold_id <- integer(nrow(X))
    for (cls in names(class_indices)) {
      idx <- class_indices[[cls]]
      fold_id[idx] <- sample(rep(seq_len(inner_k), length.out = length(idx)))
    }
    lapply(seq_len(inner_k), function(f) which(fold_id == f))
  })()

  cluster <- NULL
  if (n_cores > 1) {
    cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cluster))
    parallel::clusterExport(
      cluster, c("panel_cv_auc", "slim_auc", "X", "y", "inner_folds"),
      envir = environment()
    )
  }

  panel <- character(0)
  step_auc <- numeric(0)

  for (step in seq_len(k_max)) {
    remaining <- setdiff(pool, panel)

    score_candidate <- function(gene) {
      panel_cv_auc(X, y, c(panel, gene), inner_folds)
    }
    scores <- if (!is.null(cluster)) {
      parallel::parLapply(cluster, remaining, score_candidate)
    } else {
      lapply(remaining, score_candidate)
    }
    scores <- unlist(scores)

    # Ties go to the earlier gene in univariate order: `remaining` inherits
    # pool's order, and which.max takes the first maximum.
    best <- remaining[which.max(scores)]
    panel <- c(panel, best)
    step_auc <- c(step_auc, max(scores))
  }

  list(path = panel, pool = pool, step_auc = step_auc,
       inner_folds = inner_folds)
}
