#!/usr/bin/env Rscript
# Test runner for data-driven modules + combined soft prior.
# Usage: Rscript redesign/run_tests_datadriven.R   (from the project root)

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_file(
  "redesign/tests/test_datadriven.R",
  reporter = "progress"
)

df <- as.data.frame(results)
if (any(df$failed > 0) || any(df$error)) {
  quit(status = 1L, save = "no")
}
cat("\nAll data-driven module tests passed.\n")
