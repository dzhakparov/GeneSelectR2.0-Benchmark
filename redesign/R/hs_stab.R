# ==============================================================================
#  Subsampled horseshoe (added 2026-08-23). The single-fit horseshoe dies
#  on small-n sparse-signal data (IMvigor210): one posterior at n=123 puts
#  mass on nothing. This refits on B k-fold training subsamples and
#  averages the score across fits -- the stability trick applied to the
#  Bayesian score. Each fit uses fewer iterations (2000/500) since
#  averaging over B fits smooths Monte Carlo noise anyway.
#
#  Score_gene = mean over fits of (mean|beta| x PIP), restricted to the
#  genes kept in that subsample's training matrix (all of them here --
#  no filtering inside subsamples).
# ==============================================================================

hs_probit_stab <- function(X, y, B = 50, k_folds = 5, n_iter = 2000,
                           burn = 500, eps = 0.05, random_seed = 42,
                           n_cores = 2) {
  y_bin <- as.integer(droplevels(y)) - 1L
  subsamples <- create_subsamples(y, B = B, random_seed = random_seed,
                                  scheme = "kfold", k_folds = k_folds)

  fit_one <- function(idx) {
    tr <- subsamples[[idx]]$train
    fit <- hs_probit(X[tr, , drop = FALSE], y_bin[tr],
                     n_iter = n_iter, burn = burn, eps = eps,
                     seed = random_seed + idx * 31L)
    fit$score[colnames(X)]
  }

  scores <- parallel::mclapply(seq_along(subsamples), fit_one,
                               mc.cores = n_cores)
  smat <- do.call(rbind, scores)
  if (anyNA(smat)) stop("hs_probit_stab: a subsample fit failed")
  score <- colMeans(smat)
  names(score) <- colnames(X)

  list(ranking = names(sort(score, decreasing = TRUE)),
       score = sort(score, decreasing = TRUE))
}
