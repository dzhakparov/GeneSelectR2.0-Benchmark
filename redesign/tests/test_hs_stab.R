#  Tests for the subsampled horseshoe (redesign/R/hs_stab.R).

suppressPackageStartupMessages(library(testthat))

proj_root <- normalizePath(file.path(getwd(), "..", ".."))
setwd(proj_root)
for (f in list.files(file.path(proj_root, "package", "GeneSelectR", "R"),
                     full.names = TRUE)) source(f)
source(file.path(proj_root, "redesign", "R", "horseshoe.R"))
source(file.path(proj_root, "redesign", "R", "hs_stab.R"))

make_toy <- function(n = 120, p = 60, n_signal = 5, effect = 1.2,
                     seed = 9) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, paste0("G", seq_len(p))))
  sig <- paste0("G", seq_len(n_signal))
  eta <- rowSums(X[, sig, drop = FALSE]) * effect / sqrt(n_signal)
  ybin <- rbinom(n, 1, pnorm(eta))
  list(X = X, ybin = factor(ybin), sig = sig)
}

test_that("structure and determinism", {
  toy <- make_toy()
  a <- hs_probit_stab(toy$X, toy$ybin, B = 10, k_folds = 5, n_iter = 500,
                      burn = 200, n_cores = 2)
  b <- hs_probit_stab(toy$X, toy$ybin, B = 10, k_folds = 5, n_iter = 500,
                      burn = 200, n_cores = 2)
  expect_named(a, c("ranking", "score"))
  expect_equal(length(a$ranking), 60)
  expect_true(all(is.finite(a$score)))
  expect_identical(a$score, b$score)
})

test_that("planted signal ranks at top", {
  toy <- make_toy(n = 150, p = 80, n_signal = 5)
  fit <- hs_probit_stab(toy$X, toy$ybin, B = 10, k_folds = 5, n_iter = 800,
                        burn = 300, n_cores = 2)
  expect_true(all(toy$sig %in% fit$ranking[1:20]))
})

test_that("null data: no dominant genes", {
  set.seed(4)
  X <- matrix(rnorm(120 * 60), nrow = 120,
              dimnames = list(NULL, paste0("G", 1:60)))
  y <- factor(rbinom(120, 1, 0.5))
  fit <- hs_probit_stab(X, y, B = 10, k_folds = 5, n_iter = 500,
                        burn = 200, n_cores = 2)
  #  Absolute scale: with no signal, no gene's mean|beta| x PIP should be
  #  large (signal toy genes score ~0.3+; null max is ~0.07).
  expect_lt(max(fit$score), 0.15)
})
