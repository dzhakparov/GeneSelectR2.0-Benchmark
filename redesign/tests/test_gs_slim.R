# =============================================================================
# GeneSelectR Slim - test suite
# =============================================================================
#
# Three layers, mirroring the code:
#   1. Deviance utility: hand-computed exactness, conditioning, weighting,
#      clipping, sign consistency, degenerate inputs.
#   2. Joint null: hand-computed exceedances, leave-one-out correctness,
#      Fisher statistic, uniformity under H0 (the "no invented results"
#      guard), extremes, determinism.
#   3. Driver: end-to-end signal recovery, reproducibility, serial/parallel
#      identity, no-signal calibration, class-convention handling, validation.
#
# The driver end-to-end fits are EXPENSIVE (~170 cv.glmnet fits each), so they
# are computed once at file level and shared across the test blocks that only
# assert on their outputs.
#
# Run with: Rscript redesign/run_tests.R

library(testthat)

# test_file() runs with the working directory set to THIS file's directory,
# so the sources are reached relative to it.
for (f in c("utility_deviance.R", "joint_null.R", "gs_slim.R")) {
  source(file.path("..", "R", f))
}


# =============================================================================
# 1. Deviance utility
# =============================================================================

test_that(".binomial_deviance matches the explicit Bernoulli formula", {
  eta <- c(-1, 0, 2, 0.5)
  y01 <- c(0, 1, 1, 0)
  manual <- -2 * sum(y01 * eta - log(1 + exp(eta)))
  expect_equal(.binomial_deviance(eta, y01), manual, tolerance = 1e-12)
})

test_that("deviance utility equals hand-computed contribution, single subsample", {
  # 6 OOB samples, 2 genes; only gene 1 has a nonzero coefficient.
  X <- matrix(c(0.2, -0.5, 1.1, 0.3, -0.8, 0.6,
                1.7,  0.4, -0.2, 0.9,  0.1, -1.3),
              nrow = 6, dimnames = list(NULL, c("g1", "g2")))
  y <- factor(c(0, 0, 1, 0, 1, 1))
  beta <- c(0.5, 0)
  intercept <- 0.2
  coef_matrix <- matrix(beta, ncol = 1)

  subsamples <- list(list(train = 1:6, oob = 1:6))
  util <- compute_deviance_utility(
    X, y, subsamples, coef_matrix, intercept,
    candidate_idx = 1:2, auc_weight = FALSE, verbose = FALSE
  )

  y01 <- as.numeric(y) - 1
  eta <- X %*% beta + intercept
  eta_minus <- eta - beta[1] * X[, 1]
  expected <- .binomial_deviance(eta_minus, y01) - .binomial_deviance(eta, y01)

  expect_equal(unname(util$u["g1"]), max(expected, 0), tolerance = 1e-10)
  expect_equal(unname(util$u["g2"]), 0)
  expect_equal(unname(util$n_selected_b["g1"]), 1)
  expect_equal(unname(util$n_selected_b["g2"]), 0)
})

test_that("utility averages ONLY over subsamples that selected the gene", {
  X <- matrix(rnorm(24), nrow = 6,
              dimnames = list(NULL, paste0("g", 1:4)))
  y <- factor(c(0, 0, 0, 1, 1, 1))
  # Gene 1 selected in subsamples 1 and 3, not 2.
  coef_matrix <- matrix(0, nrow = 4, ncol = 3)
  coef_matrix[1, ] <- c(0.4, 0, -0.3)
  intercepts <- c(0.1, 0.2, 0.3)
  subsamples <- lapply(1:3, function(b) list(train = 1:6, oob = 1:6))

  util <- compute_deviance_utility(
    X, y, subsamples, coef_matrix, intercepts,
    candidate_idx = 1:4, auc_weight = FALSE, verbose = FALSE
  )

  y01 <- as.numeric(y) - 1
  d <- numeric(3)
  for (b in c(1, 3)) {
    beta <- coef_matrix[, b]
    eta <- as.numeric(X %*% beta) + intercepts[b]
    eta_minus <- eta - beta[1] * X[, 1]
    d[b] <- .binomial_deviance(eta_minus, y01) - .binomial_deviance(eta, y01)
  }
  expect_equal(unname(util$u["g1"]), max(mean(c(d[1], d[3])), 0),
               tolerance = 1e-10)
  expect_equal(unname(util$n_selected_b["g1"]), 2)
})

test_that("AUC weighting excludes non-predictive subsample models", {
  X <- matrix(rnorm(18), nrow = 6, dimnames = list(NULL, paste0("g", 1:3)))
  y <- factor(c(0, 0, 0, 1, 1, 1))
  coef_matrix <- matrix(0, nrow = 3, ncol = 3)
  coef_matrix[1, ] <- c(0.4, 0.4, 0.4)   # selected everywhere
  intercepts <- c(0, 0, 0)
  subsamples <- lapply(1:3, function(b) list(train = 1:6, oob = 1:6))
  # Weights: max(auc - 0.5, 0) = c(0.4, 0, 0.4); subsample 2 is dead weight.
  auc_vec <- c(0.9, 0.4, 0.9)

  util <- compute_deviance_utility(
    X, y, subsamples, coef_matrix, intercepts,
    candidate_idx = 1:3, auc_vec = auc_vec, auc_weight = TRUE, verbose = FALSE
  )

  # With equal coefficients and intercepts, deltas are identical across b,
  # so the weighted mean equals the unweighted single-subsample delta here.
  y01 <- as.numeric(y) - 1
  beta <- coef_matrix[, 1]
  eta <- as.numeric(X %*% beta)
  eta_minus <- eta - beta[1] * X[, 1]
  d <- .binomial_deviance(eta_minus, y01) - .binomial_deviance(eta, y01)
  expect_equal(unname(util$u["g1"]), max(d, 0), tolerance = 1e-10)

  # If ALL selecting subsamples are non-predictive, the gene scores 0.
  util_dead <- compute_deviance_utility(
    X, y, subsamples, coef_matrix, intercepts,
    candidate_idx = 1:3, auc_vec = c(0.5, 0.4, 0.5), auc_weight = TRUE,
    verbose = FALSE
  )
  expect_equal(unname(util_dead$u["g1"]), 0)
})

test_that("harmful contributions are clipped to zero, not allowed negative", {
  # beta > 0 on a gene whose HIGH values mark the negative class: zeroing it
  # IMPROVES the fit, so the raw delta is negative.
  X <- matrix(c(2, 2, 2, -2, -2, -2), nrow = 6,
              dimnames = list(NULL, "g1"))
  y <- factor(c(0, 0, 0, 1, 1, 1))
  coef_matrix <- matrix(1.0, nrow = 1, ncol = 1)
  subsamples <- list(list(train = 1:6, oob = 1:6))

  y01 <- as.numeric(y) - 1
  eta <- X[, 1] * 1.0
  raw_delta <- .binomial_deviance(eta - X[, 1], y01) -
    .binomial_deviance(eta, y01)
  expect_true(raw_delta < 0)   # confirm the construction is really harmful

  util <- compute_deviance_utility(
    X, y, subsamples, coef_matrix, 0,
    candidate_idx = 1L, auc_weight = FALSE, verbose = FALSE
  )
  expect_equal(unname(util$u["g1"]), 0)
})

test_that("sign consistency reflects the majority sign over selecting fits", {
  X <- matrix(rnorm(12), nrow = 6, dimnames = list(NULL, c("g1", "g2")))
  y <- factor(c(0, 0, 0, 1, 1, 1))
  coef_matrix <- matrix(0, nrow = 2, ncol = 4)
  coef_matrix[1, ] <- c(0.3, 0.2, -0.5, 0.1)   # 3 positive, 1 negative
  subsamples <- lapply(1:4, function(b) list(train = 1:6, oob = 1:6))

  util <- compute_deviance_utility(
    X, y, subsamples, coef_matrix, rep(0, 4),
    candidate_idx = 1:2, auc_weight = FALSE, verbose = FALSE
  )
  expect_equal(unname(util$sign_consistency["g1"]), 0.75)
  expect_true(is.na(util$sign_consistency["g2"]))  # never selected
})

test_that("all-zero coefficient matrix yields all-zero utility, no errors", {
  X <- matrix(rnorm(12), nrow = 6, dimnames = list(NULL, c("g1", "g2")))
  y <- factor(c(0, 0, 0, 1, 1, 1))
  coef_matrix <- matrix(0, nrow = 2, ncol = 2)
  subsamples <- lapply(1:2, function(b) list(train = 1:6, oob = 1:6))

  util <- compute_deviance_utility(
    X, y, subsamples, coef_matrix, c(0, 0),
    candidate_idx = 1:2, auc_weight = FALSE, verbose = FALSE
  )
  expect_true(all(util$u == 0))
  expect_true(all(util$n_selected_b == 0))
  expect_true(all(is.na(util$delta_matrix)))
})

test_that("utility validates its inputs", {
  X <- matrix(rnorm(12), nrow = 6, dimnames = list(NULL, c("g1", "g2")))
  y <- factor(c(0, 0, 0, 1, 1, 1))
  expect_error(
    compute_deviance_utility(X, y, list(list(train = 1:6, oob = 1:6)),
                             matrix(0, 2, 1), c(0, 0),  # 2 intercepts, 1 ss
                             candidate_idx = 1:2, verbose = FALSE),
    "intercepts"
  )
})


# =============================================================================
# 2. Joint null
# =============================================================================

test_that(".exceedance is the add-one permutation p-value", {
  null_mat <- matrix(c(1, 2, 3, 4,   5, 5, 5, 5), nrow = 4, byrow = FALSE)
  e <- get(".exceedance")(c(2.5, 1), null_mat)
  expect_equal(e[1], (1 + 2) / 5)   # 2 of 4 nulls >= 2.5
  expect_equal(e[2], (1 + 4) / 5)   # all 4 nulls >= 1
})

test_that(".loo_exceedance leaves the query row out", {
  null_mat <- matrix(c(1, 2, 3), nrow = 3)
  loo <- get(".loo_exceedance")(null_mat)
  expect_equal(loo[1, 1], (1 + 2) / 3)  # others 2,3 both >= 1
  expect_equal(loo[2, 1], (1 + 1) / 3)  # only 3 >= 2
  expect_equal(loo[3, 1], (1 + 0) / 3)  # none >= 3
})

test_that("observed Fisher statistic equals -2 * sum of log exceedances", {
  set.seed(11)
  null_pi <- matrix(runif(40 * 30), 40)
  null_u  <- matrix(rexp(40 * 30), 40)
  res <- combine_evidence_joint(runif(30), rexp(30), null_pi, null_u)
  expect_equal(res$T_obs, -2 * (log(res$e_pi) + log(res$e_u)),
               tolerance = 1e-12)
})

test_that("p-values are finite, in (0, 1], and inputs are validated", {
  set.seed(12)
  null_pi <- matrix(runif(40 * 25), 40)
  null_u  <- matrix(rexp(40 * 25), 40)
  res <- combine_evidence_joint(runif(25), rexp(25), null_pi, null_u)
  expect_true(all(is.finite(res$p_combined)))
  expect_true(all(res$p_combined > 0 & res$p_combined <= 1))

  expect_error(combine_evidence_joint(runif(25), rexp(24), null_pi, null_u),
               "one column per gene")
  expect_error(combine_evidence_joint(runif(25), rexp(25),
                                      null_pi, null_u[1:30, ]),
               "identical dimensions")
  expect_warning(combine_evidence_joint(runif(25), rexp(25),
                                        null_pi[1:10, ], null_u[1:10, ]),
                 "p-resolution")
})

test_that("under the global null, p_combined is not enriched for small values", {
  # Observed statistics drawn from the SAME distributions as the null rows:
  # any deviation from uniformity here is a bug in the construction, not a
  # property of the data. Fixed seed -> deterministic bounds.
  set.seed(13)
  M <- 99; p <- 400
  null_pi <- matrix(runif(M * p), M)
  null_u  <- matrix(rexp(M * p), M)
  res <- combine_evidence_joint(runif(p), rexp(p), null_pi, null_u)

  expect_gt(mean(res$p_combined), 0.35)
  expect_lt(mean(res$p_combined), 0.65)
  expect_lt(mean(res$p_combined < 0.1), 0.20)
})

test_that("overwhelming signal attains the resolution floor", {
  set.seed(14)
  M <- 99; p <- 50
  null_pi <- matrix(runif(M * p, max = 0.5), M)
  null_u  <- matrix(rexp(M * p), M)
  res <- combine_evidence_joint(rep(1, p), rep(1e6, p), null_pi, null_u)
  expect_true(all(res$p_combined == 1 / (M + 1)))
})

test_that("combination is deterministic", {
  set.seed(15)
  null_pi <- matrix(runif(30 * 20), 30, 20)
  null_u  <- matrix(rexp(30 * 20), 30, 20)
  pi_obs <- runif(20); u_obs <- rexp(20)
  r1 <- combine_evidence_joint(pi_obs, u_obs, null_pi, null_u)
  r2 <- combine_evidence_joint(pi_obs, u_obs, null_pi, null_u)
  expect_identical(r1$p_combined, r2$p_combined)
})


# =============================================================================
# 3. Driver: building blocks
# =============================================================================

test_that("slim_subsamples is stratified, exhaustive per repeat, seeded", {
  y <- factor(rep(c("a", "b"), each = 20))
  ss <- slim_subsamples(y, B = 50, k_folds = 5, random_seed = 1)
  expect_equal(length(ss), 50)
  # First repeat (5 folds) partitions all 40 samples.
  first_repeat_oob <- unlist(lapply(ss[1:5], `[[`, "oob"))
  expect_equal(sort(first_repeat_oob), 1:40)
  # Every OOB fold carries both classes.
  for (s in ss) {
    expect_true(all(c("a", "b") %in% y[s$oob]))
  }
  # Same seed, same splits.
  ss2 <- slim_subsamples(y, B = 50, k_folds = 5, random_seed = 1)
  expect_identical(ss, ss2)
})

test_that("slim_auc has a fixed direction", {
  y <- factor(rep(c("neg", "pos"), each = 5))
  expect_equal(slim_auc(y, c(1:5, 6:10)), 1)        # perfect
  expect_equal(slim_auc(y, c(6:10, 1:5)), 0)        # perfectly wrong
  expect_true(is.na(slim_auc(factor(rep("pos", 10)), 1:10)))
})


# =============================================================================
# 4. Driver: end-to-end (shared fits, computed once)
# =============================================================================

make_synthetic <- function(n = 64, p = 40, n_signal = 5, seed = 101) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, paste0("g", 1:p)))
  # Correlated signal block, as expression data actually looks.
  latent <- rnorm(n)
  for (j in seq_len(n_signal)) X[, j] <- 0.7 * latent + 0.7 * X[, j]
  eta <- 1.5 * latent + rnorm(n, sd = 0.5)
  y <- factor(rbinom(n, 1, plogis(eta)), levels = c(0, 1),
              labels = c("ctrl", "case"))
  list(X = X, y = y, signal = paste0("g", 1:n_signal))
}

# --- Shared fit WITH signal (one driver run, many assertions) ----------------
dat_signal <- make_synthetic()
fit_signal <- gs_slim_fit(dat_signal$X, dat_signal$y, alpha_grid = 0.5,
                          B = 12, k_folds = 4, M_null = 20, B_null = 8,
                          n_cores = 1, random_seed = 7, verbose = FALSE)

# --- Shared fit WITHOUT signal (the "no invented results" guard) -------------
set.seed(104)
n0 <- 64; p0 <- 40
X_nosig <- matrix(rnorm(n0 * p0), nrow = n0,
                  dimnames = list(NULL, paste0("g", 1:p0)))
y_nosig <- factor(rbinom(n0, 1, 0.5), levels = c(0, 1),
                  labels = c("ctrl", "case"))
fit_nosig <- gs_slim_fit(X_nosig, y_nosig, alpha_grid = 0.5,
                         B = 12, k_folds = 4, M_null = 20, B_null = 8,
                         n_cores = 1, random_seed = 13, verbose = FALSE)


test_that("driver output is complete, sorted, finite, and bounded", {
  tab <- fit_signal$gene_table
  expect_equal(sort(tab$gene), sort(colnames(dat_signal$X)))
  expect_true(all(diff(tab$p_combined) >= 0))
  expect_true(all(is.finite(tab$p_combined)))
  expect_true(all(tab$p_combined > 0 & tab$p_combined <= 1))
  expect_true(all(tab$pi_raw >= 0 & tab$pi_raw <= 1))
  expect_true(all(tab$u_dev >= 0))
  expect_equal(length(fit_signal$ranked), ncol(dat_signal$X))
  expect_equal(fit_signal$ranked, tab$gene)
})

test_that("driver ranks true signal above noise end-to-end", {
  tab <- fit_signal$gene_table
  is_signal <- tab$gene %in% dat_signal$signal
  expect_lt(mean(tab$p_combined[is_signal]),
            mean(tab$p_combined[!is_signal]))
})

test_that("on signal-free data, small p-values are not overproduced", {
  # A method that reports discoveries on data with no signal is manufacturing
  # them. Bounds are lenient (one random dataset, not a uniformity proof) but
  # any systematic bug -- inverted exceedances, self-referential nulls,
  # leakage -- blows straight past them.
  frac_small <- mean(fit_nosig$gene_table$p_combined < 0.1)
  expect_lt(frac_small, 0.35)
  expect_gt(median(fit_nosig$gene_table$p_combined), 0.25)
})

test_that("driver is exactly reproducible at the same seed", {
  dat <- make_synthetic(seed = 102)
  run <- function() gs_slim_fit(dat$X, dat$y, alpha_grid = c(0.5, 1.0),
                                B = 8, k_folds = 4, M_null = 20, B_null = 6,
                                n_cores = 1, random_seed = 9, verbose = FALSE)
  r1 <- run()
  r2 <- run()
  expect_identical(r1$gene_table, r2$gene_table)
  expect_identical(r1$best_alpha, r2$best_alpha)
  expect_identical(r1$null_pi, r2$null_pi)
})

test_that("serial and parallel execution give identical results", {
  skip_if(parallel::detectCores(logical = FALSE) < 2,
          "needs at least 2 cores")
  dat <- make_synthetic(seed = 103)
  run <- function(cores) gs_slim_fit(dat$X, dat$y, alpha_grid = 0.5,
                                     B = 8, k_folds = 4, M_null = 20,
                                     B_null = 6,
                                     n_cores = cores, random_seed = 11,
                                     verbose = FALSE)
  serial <- run(1)
  par <- run(2)
  expect_identical(serial$gene_table, par$gene_table)
  expect_identical(serial$null_pi, par$null_pi)
  expect_identical(serial$null_u, par$null_u)
})

test_that("class convention: the SECOND factor level is modelled", {
  # Signal gene is HIGH in "Responder", the FIRST level. The pipeline must
  # still find it -- it models P(NonResponder), so gene 1 enters with a
  # negative coefficient, and both the AUC direction and the deviance must
  # follow that convention consistently.
  set.seed(105)
  n <- 64; p <- 30
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, paste0("g", 1:p)))
  y <- factor(rep(c("NonResponder", "Responder"), each = n / 2),
              levels = c("Responder", "NonResponder"))
  X[y == "Responder", 1] <- X[y == "Responder", 1] + 2

  fit <- gs_slim_fit(X, y, alpha_grid = 0.5,
                     B = 8, k_folds = 4, M_null = 20, B_null = 6,
                     n_cores = 1, random_seed = 15, verbose = FALSE)
  expect_true("g1" %in% fit$gene_table$gene[1:5])
  expect_gt(fit$gene_table$pi_raw[fit$gene_table$gene == "g1"], 0.5)
})

test_that("driver validates its inputs", {
  dat <- make_synthetic(n = 40, p = 20, seed = 106)
  # Fewer rows in X than outcomes: the nrow/length mismatch must be caught.
  expect_error(gs_slim_fit(dat$X[1:30, ], dat$y, verbose = FALSE),
               "nrow")
  expect_error(gs_slim_fit(unname(dat$X), dat$y, verbose = FALSE),
               "gene-name")
  expect_error(gs_slim_fit(dat$X, as.character(dat$y), verbose = FALSE),
               "two-level factor")
  X_bad <- dat$X; X_bad[1, 1] <- NA
  expect_error(gs_slim_fit(X_bad, dat$y, verbose = FALSE), "non-finite")
})
