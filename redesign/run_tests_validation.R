#!/usr/bin/env Rscript
# Test runner for the validation benchmark.
# Usage: Rscript redesign/run_tests_validation.R   (from the project root)

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_file(
  "redesign/tests/test_validation.R",
  reporter = "progress"
)

df <- as.data.frame(results)
if (any(df$failed > 0) || any(df$error)) {
  quit(status = 1L, save = "no")
}
cat("\nAll validation-benchmark tests passed.\n")
