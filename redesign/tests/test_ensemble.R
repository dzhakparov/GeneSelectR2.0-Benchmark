#  Tests for redesign/R/ensemble.R

suppressPackageStartupMessages(library(testthat))
proj_root <- normalizePath(file.path(getwd(), "..", ".."))
setwd(proj_root)
source(file.path(proj_root, "redesign", "R", "ensemble.R"))

test_that("rank_average averages ranks and keeps the full pool", {
  pool <- paste0("G", 1:6)
  r1 <- c("G1", "G2", "G3", "G4", "G5", "G6")
  r2 <- c("G1", "G3", "G2", "G4", "G5", "G6")
  out <- rank_average(list(r1, r2), pool)
  expect_setequal(out, pool)
  expect_equal(out[1], "G1")              # mean rank 1.0
  #  G2 and G3 tie at mean 2.5; ties follow pool order.
  expect_setequal(out[2:3], c("G2", "G3"))
  expect_equal(out[4:6], c("G4", "G5", "G6"))
})

test_that("rank_average penalises genes missing from a short ranking", {
  pool <- paste0("G", 1:5)
  full <- pool
  gated <- c("G3", "G1")                  # covers only 2 genes
  out <- rank_average(list(full, gated), pool)
  expect_setequal(out, pool)
  #  G3 (rank 3, 1) and G1 (rank 1, 2) must beat genes missing from gated.
  expect_true(all(c("G3", "G1") %in% out[1:2]))
})

test_that("kswitch picks the right ranking at the cut", {
  small <- c("a", "b"); large <- c("c", "d")
  expect_equal(kswitch_ranking(10, 20, small, large), small)
  expect_equal(kswitch_ranking(20, 20, small, large), small)
  expect_equal(kswitch_ranking(50, 20, small, large), large)
})
