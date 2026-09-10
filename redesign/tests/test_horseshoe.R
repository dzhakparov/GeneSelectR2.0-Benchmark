#  Tests for the horseshoe probit sampler (redesign/R/horseshoe.R).

suppressPackageStartupMessages(library(testthat))

proj_root <- normalizePath(file.path(getwd(), "..", ".."))
setwd(proj_root)
source(file.path(proj_root, "redesign", "R", "horseshoe.R"))

make_toy <- function(n = 100, p = 60, n_signal = 5, effect = 1.0,
                     seed = 9) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, paste0("G", seq_len(p))))
  sig <- paste0("G", seq_len(n_signal))
  beta_true <- setNames(rep(0, p), colnames(X))
  beta_true[sig] <- effect
  eta <- as.vector(X %*% beta_true)
  ybin <- rbinom(n, 1, pnorm(eta))
  list(X = X, ybin = ybin, sig = sig)
}

test_that("return structure and sizes", {
  toy <- make_toy()
  fit <- hs_probit(toy$X, toy$ybin, n_iter = 300, burn = 100, seed = 1)
  expect_named(fit, c("ranking", "score", "mean_abs", "pip", "draws"))
  expect_equal(nrow(fit$draws), 200)
  expect_equal(ncol(fit$draws), 60)
  expect_equal(length(fit$ranking), 60)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_true(all(fit$mean_abs >= 0))
})

test_that("planted signal genes rank at the top", {
  toy <- make_toy(n = 150, p = 100, n_signal = 5, effect = 1.2)
  fit <- hs_probit(toy$X, toy$ybin, n_iter = 1500, burn = 500, seed = 1)
  #  All 5 true genes should be in the top 15 of 100.
  expect_true(all(toy$sig %in% fit$ranking[1:15]))
  #  Their PIPs should clearly exceed the noise genes' median.
  expect_gt(median(fit$pip[toy$sig]), median(fit$pip) + 0.2)
})

test_that("null data produces no dominant genes", {
  set.seed(4)
  X <- matrix(rnorm(100 * 60), nrow = 100,
              dimnames = list(NULL, paste0("G", 1:60)))
  ybin <- rbinom(100, 1, 0.5)
  fit <- hs_probit(X, ybin, n_iter = 1000, burn = 500, seed = 1)
  #  Max PIP should stay low; nothing is really included.
  expect_lt(max(fit$pip), 0.6)
  #  Score spread should be small relative to signal case.
  expect_lt(max(fit$mean_abs), 0.5)
})

test_that("deterministic for fixed seed", {
  toy <- make_toy()
  a <- hs_probit(toy$X, toy$ybin, n_iter = 300, burn = 100, seed = 7)
  b <- hs_probit(toy$X, toy$ybin, n_iter = 300, burn = 100, seed = 7)
  expect_identical(a$score, b$score)
  c <- hs_probit(toy$X, toy$ybin, n_iter = 300, burn = 100, seed = 8)
  expect_false(identical(a$score, c$score))
})

test_that("no NA/NaN in scores on near-separable data", {
  toy <- make_toy(n = 80, p = 40, n_signal = 3, effect = 3)
  fit <- hs_probit(toy$X, toy$ybin, n_iter = 500, burn = 200, seed = 1)
  expect_true(all(is.finite(fit$score)))
})
