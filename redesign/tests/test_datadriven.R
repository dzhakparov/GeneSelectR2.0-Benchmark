#  Tests for data-driven modules + combined (Hallmark + data-driven) soft prior.

suppressPackageStartupMessages(library(testthat))

proj_root <- normalizePath(file.path(getwd(), "..", ".."))
setwd(proj_root)
source(file.path(proj_root, "redesign", "R", "bio_prior.R"))
source(file.path(proj_root, "redesign", "R", "datadriven_modules.R"))

make_blocks <- function(n = 120, seed = 5) {
  #  Three correlated blocks of 20 genes; block 1 also carries class signal.
  set.seed(seed)
  ybin <- rep(0:1, length.out = n)
  X <- matrix(rnorm(n * 200), nrow = n,
              dimnames = list(NULL, paste0("G", seq_len(200))))
  for (b in 0:2) {
    base <- rnorm(n)
    idx <- (b * 20 + 1):(b * 20 + 20)
    for (j in idx) X[, j] <- base + rnorm(n, sd = 0.3)
  }
  X[ybin == 1, 1:20] <- X[ybin == 1, 1:20] + 1.0
  list(X = X, ybin = ybin)
}

test_that("datadriven_modules recovers planted correlation blocks", {
  toy <- make_blocks()
  sets <- datadriven_modules(toy$X, min_size = 10, max_size = 200)
  expect_true(length(sets) >= 2)
  #  The three planted blocks should appear (roughly) as modules.
  blocks <- list(1:20, 21:40, 41:60)
  found <- vapply(blocks, function(b) {
    any(vapply(sets, function(s)
      length(intersect(s, paste0("G", b))) >= 15, logical(1)))
  }, logical(1))
  expect_true(all(found))
})

test_that("datadriven_modules is deterministic", {
  toy <- make_blocks()
  a <- datadriven_modules(toy$X)
  b <- datadriven_modules(toy$X)
  expect_identical(a, b)
})

test_that("size filters drop tiny and huge modules", {
  toy <- make_blocks()
  sets <- datadriven_modules(toy$X, min_size = 10, max_size = 25)
  sizes <- vapply(sets, length, integer(1))
  expect_true(all(sizes >= 10 & sizes <= 25))
})

test_that("signal block gets high z under module_scores", {
  toy <- make_blocks()
  sets <- datadriven_modules(toy$X, min_size = 10, max_size = 200)
  ms <- module_scores(toy$X, toy$ybin, B = 200, seed = 7, sets = sets)
  sig_set <- names(which.max(vapply(sets, function(s)
    length(intersect(s, paste0("G", 1:20))), integer(1))))
  expect_true(ms$z[sig_set] > 3)
  #  Non-signal modules should sit near zero.
  other <- setdiff(names(sets), sig_set)
  expect_true(all(abs(ms$z[other]) < 4))
})

test_that("combined prior takes the max boost across both sources", {
  scores <- c(g_hall = 1.0, g_data = 1.0, g_both = 1.0, g_none = 1.0)
  hall_sets <- list(H1 = c("g_hall", "g_both"))
  dd_sets <- list(DC1 = c("g_data", "g_both"))
  z_hall <- c(H1 = 2)      # boost 1+0.5*2 = 2
  z_dd <- c(DC1 = 4)       # boost 1+0.5*4 = 3

  r_hall <- soft_prior_rescore(scores, hall_sets, z_hall, w = 0.5)
  r_dd <- soft_prior_rescore(scores, dd_sets, z_dd, w = 0.5)
  r_cb <- soft_prior_rescore(scores, c(hall_sets, dd_sets),
                             c(z_hall, z_dd), w = 0.5)
  #  g_data and g_both both reach 3.0 (tie -> original order); g_hall 2.0;
  #  g_none 1.0.
  expect_equal(r_cb[1:2], c("g_data", "g_both"))
  expect_equal(r_dd[1:2], c("g_data", "g_both"))
  expect_equal(r_hall[1:2], c("g_hall", "g_both"))  # tie 2.0, original order
  expect_equal(r_cb[4], "g_none")
})
