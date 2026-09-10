# ==============================================================================
#  "Prediction-first" variant (changes 1+2+3 from the 2026-08-22 review of
#  the full-benchmark loss):
#
#    1. Score = |mean coefficient over subsample fits| * freq^gamma,
#       gamma = 0.25. The coefficient magnitude (the predictive part) is the
#       backbone; selection frequency only modulates. The old geometric mean
#       gave stability veto power, which the benchmark showed costs AUC.
#    2. Redundancy filter at panel level: greedy walk, skip a gene if
#       |cor| > tau with anything already picked (correlations from the
#       TRAIN matrix only). This is the one thing mRMR -- which won two
#       datasets -- has that the pipeline lacked.
#    3. Alpha picked by CV deviance (min cvm) on the full training set, the
#       same criterion the ElasticNet competitor uses, instead of the
#       internal OOB AUC that mostly picked alpha = 0.5.
# ==============================================================================

#  Fit subsample elastic nets and return the prediction-first score ranking.
#  X: samples x genes (already standardised by the caller). y: 2-level factor.
fit_predfirst <- function(X, y, B = 50, k_folds = 5,
                          alpha_grid = c(0.5, 1.0), gamma = 0.25,
                          random_seed = 42, n_cores = 2) {
  y_bin <- as.integer(droplevels(y)) - 1L
  n_inner <- min(5, max(3, floor(min(table(y)) * 0.8)))

  #  Alpha by CV deviance on the full training set (change 3). cv.glmnet
  #  draws its folds from the RNG; a fixed seed per call makes the choice
  #  reproducible and gives every alpha the SAME folds (fair comparison).
  cvm <- vapply(alpha_grid, function(a) {
    withr::local_seed(random_seed + 7000L)
    fit <- glmnet::cv.glmnet(X, y_bin, family = "binomial", alpha = a,
                             nfolds = n_inner)
    min(fit$cvm)
  }, numeric(1))
  alpha <- alpha_grid[which.min(cvm)]

  subsamples <- create_subsamples(y, B = B, random_seed = random_seed,
                                  scheme = "kfold", k_folds = k_folds)

  fit_one <- function(idx) {
    tr <- subsamples[[idx]]$train
    #  Seed per subsample: cv.glmnet's internal folds must not depend on
    #  the (forked) global RNG stream.
    withr::local_seed(random_seed + idx)
    tryCatch({
      f <- glmnet::cv.glmnet(X[tr, , drop = FALSE], y_bin[tr],
                             family = "binomial", alpha = alpha,
                             nfolds = n_inner)
      as.numeric(coef(f, s = "lambda.min"))[-1]
    }, error = function(e) rep(NA_real_, ncol(X)))
  }
  coefs <- parallel::mclapply(seq_along(subsamples), fit_one,
                              mc.cores = n_cores)
  coef_mat <- do.call(rbind, coefs)
  if (anyNA(coef_mat)) stop("predfirst: a subsample fit failed")
  colnames(coef_mat) <- colnames(X)

  freq     <- colMeans(coef_mat != 0)
  mean_abs <- colMeans(abs(coef_mat))
  score    <- mean_abs * freq^gamma

  list(ranking = names(sort(score, decreasing = TRUE)),
       score = sort(score, decreasing = TRUE),
       alpha = alpha, cvm = setNames(cvm, alpha_grid))
}

#  Greedy redundancy filter over a ranking (change 2). Genes skipped for
#  redundancy are appended in original rank order after the clean prefix.
#  Correlations come from X_train only -- never the test half. The full
#  correlation matrix is computed once (vectorised) instead of per gene.
redundancy_filter <- function(ranked, X_train, tau = 0.7) {
  C <- suppressWarnings(abs(cor(X_train[, ranked, drop = FALSE])))
  C[is.na(C)] <- 0                       # zero-variance columns -> no link
  picked_idx <- integer(0)
  skipped_idx <- integer(0)
  for (i in seq_along(ranked)) {
    if (length(picked_idx) == 0 ||
        max(C[i, picked_idx]) <= tau) {
      picked_idx <- c(picked_idx, i)
    } else {
      skipped_idx <- c(skipped_idx, i)
    }
  }
  ranked[c(picked_idx, skipped_idx)]
}
