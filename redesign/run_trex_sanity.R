#!/usr/bin/env Rscript
# ==============================================================================
#  Phase 0 sanity check: does T-Rex run at all on IMvigor210 at p=2000, n~122?
#
#  Builds repeat-1/fold-1 exactly as the pilot runners do (same seed, same
#  preprocessing), then runs trex() at tFDR 0.1 and 0.2 and reports:
#    - wall time
#    - how many genes selected (pass criterion: not 0, not ~1000)
#    - whether 3 repeated runs at the same seed give the same panel
#
#  Usage (from the project root, terminal only):
#    Rscript redesign/run_trex_sanity.R
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(TRexSelector)
})
source(file.path("redesign", "R", "imvigor210_data.R"))

out_dir <- file.path("redesign", "results", "trex_sanity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- load_imvigor210()
folds <- make_stratified_folds(dat$outcome, k_folds = 5, seed = imv_random_seed)
test_idx <- folds[[1]]
train_idx <- setdiff(seq_along(dat$outcome), test_idx)

pp <- preprocess_split(dat$raw_counts, train_idx, test_idx, top_genes = 2000)
X <- pp$train
y <- as.integer(dat$outcome[train_idx]) - 1L   # 0/1 for LARS
cat(sprintf("Train matrix: n=%d p=%d | responders=%d\n\n",
            nrow(X), ncol(X), sum(y)))

run_one <- function(tFDR, seed) {
  t0 <- proc.time()[["elapsed"]]
  fit <- trex(X, y, tFDR = tFDR, K = 20, seed = seed,
              parallel_process = TRUE, parallel_max_cores = 7,
              verbose = FALSE)
  elapsed <- proc.time()[["elapsed"]] - t0
  sel <- sort(which(fit$selected_var > 0))
  list(tFDR = tFDR, seed = seed, elapsed = elapsed, selected = sel)
}

results <- list()
for (tfdr in c(0.1, 0.2)) {
  for (seed in c(1, 1, 1, 2)) {   # 3x same seed (determinism) + 1x different
    r <- run_one(tfdr, seed)
    key <- sprintf("tFDR%.2f_seed%d", tfdr, seed)
    cat(sprintf("[%s] %6.1f s | selected %d genes\n",
                key, r$elapsed, length(r$selected)))
    results[[length(results) + 1]] <- r
  }
}

#  Determinism: the three seed=1 runs at each tFDR must be identical.
for (tfdr in c(0.1, 0.2)) {
  same_seed <- Filter(function(r) r$tFDR == tfdr && r$seed == 1, results)
  panels <- lapply(same_seed, function(r) r$selected)
  identical_all <- all(vapply(panels[-1], function(p)
    identical(p, panels[[1]]), logical(1)))
  cat(sprintf("\ntFDR=%.1f determinism across 3 same-seed runs: %s\n",
              tfdr, ifelse(identical_all, "IDENTICAL", "DIFFERS")))
}

#  Cross-seed Jaccard: how much does the panel move with the dummy RNG?
for (tfdr in c(0.1, 0.2)) {
  p1 <- Filter(function(r) r$tFDR == tfdr && r$seed == 1, results)[[1]]$selected
  p2 <- Filter(function(r) r$tFDR == tfdr && r$seed == 2, results)[[1]]$selected
  j <- length(intersect(p1, p2)) / max(1L, length(union(p1, p2)))
  cat(sprintf("tFDR=%.1f cross-seed Jaccard (seed1 vs seed2): %.3f\n", tfdr, j))
}

#  Save panels for later comparison with GS_harmonic top-k.
for (r in results[unique(vapply(results, function(x)
  paste(x$tFDR, x$seed), character(1)))]) {
  genes <- colnames(X)[r$selected]
  write.csv(data.frame(index = r$selected, gene = genes),
            file.path(out_dir, sprintf("panel_tFDR%.2f_seed%d.csv",
                                       r$tFDR, r$seed)),
            row.names = FALSE)
}
cat(sprintf("\nPanels written to %s\n", out_dir))
