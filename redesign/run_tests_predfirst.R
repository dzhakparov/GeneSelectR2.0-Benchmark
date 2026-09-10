#!/usr/bin/env Rscript
# Test runner for the prediction-first variant.
# Usage: Rscript redesign/run_tests_predfirst.R   (from the project root)

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_file(
  "redesign/tests/test_predfirst.R",
  reporter = "progress"
)

df <- as.data.frame(results)
if (any(df$failed > 0) || any(df$error)) {
  quit(status = 1L, save = "no")
}
cat("\nAll predfirst tests passed.\n")
