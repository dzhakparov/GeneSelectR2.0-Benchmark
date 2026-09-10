#!/usr/bin/env Rscript
# ==============================================================================
#  Phase 0 follow-up: WHY did trex() return empty panels on IMvigor210?
#
#  Diagnostics on repeat-1/fold-1 (same data construction as the pilot):
#    1. default call (max_T_stop=TRUE  -> T_stop = ceil(n/2)) at tFDR 0.2
#    2. paper setting (max_T_stop=FALSE -> T_stop = 1)        at tFDR 0.2
#    3. default call at tFDR 0.5 (how far must FDR go before anything passes?)
#  For each: T_stop, calibrated v_thresh, top-15 occurrence probabilities
#  Phi_prime, number selected. This separates "signal too weak" (Phi low for
#  all genes) from "calibration too strict" (Phi high but below v_thresh).
#
#  Usage: Rscript redesign/run_trex_diagnose.R
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(TRexSelector)
})
source(file.path("redesign", "R", "imvigor210_data.R"))

dat <- load_imvigor210()
folds <- make_stratified_folds(dat$outcome, k_folds = 5, seed = imv_random_seed)
test_idx <- folds[[1]]
train_idx <- setdiff(seq_along(dat$outcome), test_idx)
pp <- preprocess_split(dat$raw_counts, train_idx, test_idx, top_genes = 2000)
X <- pp$train
y <- as.integer(dat$outcome[train_idx]) - 1L
cat(sprintf("Train matrix: n=%d p=%d | responders=%d\n\n",
            nrow(X), ncol(X), sum(y)))

diagnose <- function(label, ...) {
  fit <- trex(X, y, K = 20, seed = 1, parallel_process = TRUE,
              parallel_max_cores = 7, verbose = FALSE, ...)
  sel <- which(fit$selected_var > 0)
  top_phi <- round(sort(fit$Phi_prime, decreasing = TRUE)[1:15], 3)
  cat(sprintf("--- %s\n", label))
  cat(sprintf("    T_stop=%s | v_thresh=%.3f | selected=%d\n",
              fit$T_stop, fit$v_thresh, length(sel)))
  cat(sprintf("    top Phi_prime: %s\n", paste(top_phi, collapse = " ")))
  if (length(sel) > 0) {
    cat(sprintf("    genes: %s\n",
                paste(colnames(X)[sel][1:min(10, length(sel))],
                      collapse = ", ")))
  }
  invisible(fit)
}

diagnose("default (max_T_stop=TRUE),  tFDR=0.2",
         tFDR = 0.2, max_T_stop = TRUE)
diagnose("paper    (max_T_stop=FALSE), tFDR=0.2",
         tFDR = 0.2, max_T_stop = FALSE)
diagnose("default (max_T_stop=TRUE),  tFDR=0.5",
         tFDR = 0.5, max_T_stop = TRUE)
