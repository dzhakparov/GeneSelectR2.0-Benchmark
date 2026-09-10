#!/usr/bin/env Rscript
# Test runner for the soft module prior.
# Usage: Rscript redesign/run_tests_softprior.R   (from the project root)

suppressPackageStartupMessages(library(testthat))

results <- testthat::test_file(
  "redesign/tests/test_softprior.R",
  reporter = "progress"
)

df <- as.data.frame(results)
if (any(df$failed > 0) || any(df$error)) {
  quit(status = 1L, save = "no")
}
cat("\nAll soft-prior tests passed.\n")
