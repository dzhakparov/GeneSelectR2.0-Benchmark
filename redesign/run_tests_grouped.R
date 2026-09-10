#!/usr/bin/env Rscript
# Test runner for the grouped-benchmark redesign.
# Usage: Rscript redesign/run_tests_grouped.R   (from the project root)

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_file(
  "redesign/tests/test_grouped_benchmark.R",
  reporter = "progress"
)

df <- as.data.frame(results)
if (any(df$failed > 0) || any(df$error)) {
  quit(status = 1L, save = "no")
}
cat("\nAll grouped-benchmark tests passed.\n")
