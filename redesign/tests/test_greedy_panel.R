# =============================================================================
# Greedy forward panel selection - tests
# =============================================================================
#
# Run with: Rscript redesign/run_tests_greedy.R

library(testthat)

# test_file() runs with the working directory set to THIS file's directory.
source(file.path("..", "R", "gs_slim.R"))       # for slim_auc, slim_subsamples
source(file.path("..", "R", "greedy_panel.R"))


make_panel_data <- function(n = 90, p = 40, seed = 201) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, paste0("g", 1:p)))
  # Three independent signal genes with different signs, plus a correlated
  # decoy that shadows g1 (the kind of gene a univariate ranker also likes).
  lp <- 1.6 * X[, 1] - 1.4 * X[, 2] + 1.2 * X[, 3]
  X[, 4] <- 0.9 * X[, 1] + rnorm(n, sd = 0.4)
  y <- factor(rbinom(n, 1, plogis(lp)), levels = c(0, 1),
              labels = c("ctrl", "case"))
  list(X = X, y = y, signal = c("g1", "g2", "g3"))
}


test_that("univariate ranking finds both directions of signal", {
  dat <- make_panel_data()
  ranked <- univariate_auc_rank(dat$X, dat$y)
  # g1, g2, g3 and the decoy g4 should all beat pure noise genes.
  expect_true(all(c("g1", "g2", "g3", "g4") %in% ranked[1:8]))
  expect_equal(length(ranked), ncol(dat$X))
})

test_that("a signal panel scores higher than a noise panel", {
  dat <- make_panel_data()
  folds <- local({
    set.seed(1)
    idx <- sample(rep(1:5, length.out = nrow(dat$X)))
    lapply(1:5, function(f) which(idx == f))
  })
  good <- panel_cv_auc(dat$X, dat$y, c("g1", "g2", "g3"), folds)
  bad  <- panel_cv_auc(dat$X, dat$y, c("g10", "g20", "g30"), folds)
  expect_gt(good, bad + 0.1)
  expect_gt(good, 0.7)
})

test_that("greedy picks the signal genes first and stays deterministic", {
  dat <- make_panel_data()
  run <- function() greedy_forward_panel(dat$X, dat$y, n_candidates = 15,
                                         k_max = 5, inner_k = 5,
                                         n_cores = 1, random_seed = 3)
  r1 <- run()
  r2 <- run()

  expect_equal(length(r1$path), 5)
  expect_false(any(duplicated(r1$path)))
  expect_true(all(r1$path %in% r1$pool))

  # The three true signal genes should appear in the first four picks.
  # (g4 the decoy may take one slot -- it genuinely helps prediction.)
  expect_gte(sum(dat$signal %in% r1$path[1:4]), 3)

  # Deterministic: same seed, same path, same scores.
  expect_identical(r1$path, r2$path)
  expect_identical(r1$step_auc, r2$step_auc)
})

test_that("greedy is identical serial and parallel", {
  skip_if(parallel::detectCores(logical = FALSE) < 2,
          "needs at least 2 cores")
  dat <- make_panel_data(seed = 202)
  serial <- greedy_forward_panel(dat$X, dat$y, n_candidates = 12, k_max = 4,
                                 inner_k = 5, n_cores = 1, random_seed = 5)
  par <- greedy_forward_panel(dat$X, dat$y, n_candidates = 12, k_max = 4,
                              inner_k = 5, n_cores = 2, random_seed = 5)
  expect_identical(serial$path, par$path)
  expect_identical(serial$step_auc, par$step_auc)
})

test_that("inputs are validated", {
  dat <- make_panel_data()
  expect_error(greedy_forward_panel(dat$X, dat$y, n_candidates = 5,
                                    k_max = 10),
               "k_max")
  expect_error(greedy_forward_panel(dat$X, dat$y[1:30]),
               "nrow")
})
