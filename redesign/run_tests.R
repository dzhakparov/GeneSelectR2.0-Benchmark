#!/usr/bin/env Rscript
# Test runner for the complete GeneSelectR redesign.
# Usage: Rscript redesign/run_tests.R   (from the project root)

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_dir(
  "redesign/tests",
  reporter = "progress"
)

df <- as.data.frame(results)
if (any(df$failed > 0) || any(df$error)) {
  quit(status = 1L, save = "no")
}
cat("\nAll redesign tests passed.\n")
