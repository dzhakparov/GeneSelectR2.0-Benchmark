#!/usr/bin/env Rscript
# Test runner for the greedy panel selector.
# Usage: Rscript redesign/run_tests_greedy.R   (from the project root)

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_file(
  "redesign/tests/test_greedy_panel.R",
  reporter = "progress"
)

df <- as.data.frame(results)
if (any(df$failed > 0) || any(df$error)) {
  quit(status = 1L, save = "no")
}
cat("\nAll greedy-panel tests passed.\n")
