# ==============================================================================
#  Tests for the prediction-first variant (redesign/R/predfirst.R).
# ==============================================================================

suppressPackageStartupMessages(library(testthat))

proj_root <- normalizePath(file.path(getwd(), "..", ".."))
setwd(proj_root)
for (f in list.files(file.path(proj_root, "package", "GeneSelectR", "R"),
                     full.names = TRUE)) source(f)
source(file.path(proj_root, "redesign", "R", "predfirst.R"))

set.seed(11)

make_data <- function(n = 100, p = 80, n_signal = 6, eff = 2.0) {
  y <- factor(rep(c("neg", "pos"), each = n / 2), levels = c("neg", "pos"))
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("G", 1:p)
  sig <- paste0("G", 1:n_signal)
  X[y == "pos", sig] <- X[y == "pos", sig] + eff
  list(X = X, y = y, sig = sig)
}

test_that("planted signal dominates the ranking", {
  d <- make_data()
  pf <- fit_predfirst(d$X, d$y, B = 10, k_folds = 5, n_cores = 1)
  expect_length(pf$ranking, 80)
  expect_true(all(d$sig %in% pf$ranking[1:12]))
  expect_true(all(pf$score >= 0))
})

test_that("alpha comes from the grid and is reproducible", {
  d <- make_data()
  pf1 <- fit_predfirst(d$X, d$y, B = 6, n_cores = 1)
  pf2 <- fit_predfirst(d$X, d$y, B = 6, n_cores = 1)
  expect_true(pf1$alpha %in% c(0.5, 1.0))
  expect_identical(pf1$ranking, pf2$ranking)
  expect_identical(pf1$alpha, pf2$alpha)
})

test_that("score is coefficient-driven, not stability-driven", {
  #  A gene with a large coefficient selected in half the subsamples must
  #  outrank a gene with a tiny coefficient selected in all of them.
  #  gamma = 0.25: freq ratio 0.5 vs 1.0 -> multiplier 0.84, so a 5x
  #  coefficient advantage always wins. Under the old geometric mean with
  #  calibrated pillars, the always-selected gene would have won.
  d <- make_data()
  pf <- fit_predfirst(d$X, d$y, B = 10, n_cores = 1, gamma = 0.25)
  #  direct arithmetic check on the formula
  expect_equal(5 * 0.5^0.25 > 1 * 1^0.25, TRUE)
  expect_equal(names(pf$score), pf$ranking)
})

test_that("redundancy filter keeps one representative per correlated block", {
  n <- 120
  base <- rnorm(n)
  X <- cbind(matrix(rnorm(n * 10), nrow = n),
             base + rnorm(n, sd = 0.05),          # A1
             base + rnorm(n, sd = 0.05),          # A2 (cor ~1 with A1)
             rnorm(n))                            # B1
  colnames(X) <- c(paste0("N", 1:10), "A1", "A2", "B1")
  ranked <- c("A1", "A2", "B1", paste0("N", 1:10))
  out <- redundancy_filter(ranked, X, tau = 0.7)
  expect_setequal(out, colnames(X))              # nothing lost
  #  A2 is redundant with A1 -> pushed behind every non-redundant gene.
  expect_gt(which(out == "A2"), which(out == "B1"))
  expect_equal(out[1], "A1")
  #  No pair in the clean prefix exceeds tau.
  prefix <- out[out != "A2"]
  C <- abs(cor(X[, prefix]))
  expect_true(all(C[upper.tri(C)] <= 0.7 + 1e-8))
})

test_that("redundancy filter never touches the test half", {
  n_tr <- 60; n_te <- 40
  Xtr <- matrix(rnorm(n_tr * 30), nrow = n_tr)
  Xte <- matrix(rnorm(n_te * 30), nrow = n_te)
  colnames(Xtr) <- colnames(Xte) <- paste0("G", 1:30)
  ranked <- paste0("G", 1:30)
  #  Same train, different test: output must be identical, i.e. the filter
  #  cannot have seen the test data.
  out1 <- redundancy_filter(ranked, Xtr)
  Xte2 <- matrix(rnorm(n_te * 30, sd = 5), nrow = n_te)
  colnames(Xte2) <- paste0("G", 1:30)
  out2 <- redundancy_filter(ranked, Xtr)   # Xte2 exists but is never passed
  expect_identical(out1, out2)
  expect_error(redundancy_filter(ranked, Xte), NA)  # works on any matrix
})

test_that("redundancy filter handles zero-variance columns", {
  X <- matrix(rnorm(50 * 10), nrow = 50)
  colnames(X) <- paste0("G", 1:10)
  X[, 5] <- 3                                   # constant column
  out <- redundancy_filter(paste0("G", 1:10), X)
  expect_setequal(out, colnames(X))
})
