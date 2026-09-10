#  Tests for the soft module prior (module_scores + soft_prior_rescore).
#  Design fixed 2026-08-22: z = (observed - null mean) / null sd under label
#  permutations; gene multiplier = 1 + w * max(0, best z); unannotated = 1.

suppressPackageStartupMessages(library(testthat))

proj_root <- normalizePath(file.path(getwd(), "..", ".."))
setwd(proj_root)
source(file.path(proj_root, "redesign", "R", "bio_prior.R"))

make_toy <- function(n = 120, p = 200, n_signal = 15, effect = 1.2,
                     seed = 11) {
  set.seed(seed)
  ybin <- rep(0:1, length.out = n)
  X <- matrix(rnorm(n * p), nrow = n,
              dimnames = list(NULL, paste0("G", seq_len(p))))
  sig <- paste0("G", seq_len(n_signal))
  X[ybin == 1, sig] <- X[ybin == 1, sig] + effect
  list(X = X, ybin = ybin, sig = sig)
}

toy_sets <- function(sig, p) {
  list(SIG_SET = sig,
       NOISE_SET = paste0("G", (p - 19):p),
       TINY_SET = paste0("G", 1:5))  # < min_members -> dropped
}

test_that("module_scores: planted module scores high z, null data does not", {
  toy <- make_toy()
  sets <- toy_sets(toy$sig, ncol(toy$X))
  ms <- module_scores(toy$X, toy$ybin, B = 200, seed = 7, sets = sets)
  expect_named(ms, c("z", "sets", "table"))
  expect_true(ms$z["SIG_SET"] > 3)
  expect_true(abs(ms$z["NOISE_SET"]) < 3)
  expect_false("TINY_SET" %in% names(ms$z))  # too small, dropped

  # Pure noise: no set should reach a large z
  set.seed(3)
  Xn <- matrix(rnorm(120 * 200), nrow = 120,
               dimnames = list(NULL, paste0("G", 1:200)))
  yn <- rep(0:1, length.out = 120)
  mn <- module_scores(Xn, yn, B = 200, seed = 7,
                      sets = list(ANY1 = paste0("G", 1:30),
                                  ANY2 = paste0("G", 31:60)))
  expect_true(all(abs(mn$z) < 4))
})

test_that("module_scores is deterministic for fixed seed", {
  toy <- make_toy()
  sets <- toy_sets(toy$sig, ncol(toy$X))
  a <- module_scores(toy$X, toy$ybin, B = 100, seed = 42, sets = sets)
  b <- module_scores(toy$X, toy$ybin, B = 100, seed = 42, sets = sets)
  expect_identical(a$z, b$z)
})

test_that("soft_prior_rescore: boost math and ranking are correct", {
  scores <- c(a = 1.0, b = 0.9, c = 0.8, orphan = 1.0)
  sets <- list(S1 = c("b", "c"))
  z <- c(S1 = 2)
  ranked <- soft_prior_rescore(scores, sets, z, w = 0.5)
  # b: 0.9*2 = 1.8 -> first; c: 0.8*2 = 1.6 -> second; a/orphan tie at 1.0
  expect_equal(ranked[1], "b")
  expect_equal(ranked[2], "c")
  expect_equal(ranked[3:4], c("a", "orphan"))  # tie order: original order
})

test_that("soft_prior_rescore: exact multiplier arithmetic", {
  scores <- c(x = 0.8, y = 0.5, z0 = 2.0)
  sets <- list(S1 = c("x", "y"))
  boosted <- soft_prior_rescore(scores, sets, z = c(S1 = 2), w = 0.5)
  # x -> 0.8*2 = 1.6, y -> 0.5*2 = 1.0, z0 untouched = 2.0
  expect_equal(boosted, c("z0", "x", "y"))
})

test_that("soft_prior_rescore: negative z gives no boost; best-of-sets wins", {
  scores <- c(g = 1.0, h = 1.0)
  sets <- list(NEG = "g", POS = "g", ALSO = "h")
  ranked <- soft_prior_rescore(scores, sets, z = c(NEG = -5, POS = 3, ALSO = 1),
                               w = 0.5)
  # g gets max(1, 1+0.5*3) = 2.5; h gets 1.5 -> g first
  expect_equal(ranked, c("g", "h"))
  # all-negative z: order unchanged
  ranked2 <- soft_prior_rescore(scores, sets, z = c(NEG = -5, POS = -1, ALSO = -2))
  expect_equal(ranked2, c("g", "h"))
})

test_that("soft_prior_rescore: genes in no set keep multiplier 1", {
  scores <- c(anno = 0.1, plain = 0.2)
  ranked <- soft_prior_rescore(scores, sets = list(S1 = "anno"),
                               z = c(S1 = 1), w = 0.5)
  # anno -> 0.15, plain stays 0.2 -> plain still first
  expect_equal(ranked[1], "plain")
})
