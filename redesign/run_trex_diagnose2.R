#!/usr/bin/env Rscript
# ==============================================================================
#  Phase 0 follow-up 2: is the empty selection a K-granularity artifact?
#
#  With K=20 experiments the occurrence frequencies Phi have resolution 1/20
#  and the calibrated threshold v_thresh was forced to 1.000 at tFDR=0.2.
#  The paper's real-data examples use larger K. Test:
#    A. tFDR=0.2, K=100, max_T_stop=TRUE
#    B. tFDR=0.2, K=100, max_T_stop=FALSE (paper T_stop=1 setting)
#    C. tFDR=0.3, K=100, max_T_stop=TRUE
#    D. tFDR=0.2, K=20, method='trex+GVS' (group version, module-level power)
#
#  Usage: Rscript redesign/run_trex_diagnose2.R
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
  t0 <- proc.time()[["elapsed"]]
  fit <- trex(X, y, seed = 1, parallel_process = TRUE,
              parallel_max_cores = 7, verbose = FALSE, ...)
  elapsed <- proc.time()[["elapsed"]] - t0
  sel <- which(fit$selected_var > 0)
  top_phi <- round(sort(fit$Phi_prime, decreasing = TRUE)[1:10], 3)
  cat(sprintf("--- %s  (%.0f s)\n", label, elapsed))
  cat(sprintf("    T_stop=%s | v_thresh=%.3f | selected=%d\n",
              fit$T_stop, fit$v_thresh, length(sel)))
  cat(sprintf("    top Phi_prime: %s\n", paste(top_phi, collapse = " ")))
  if (length(sel) > 0) {
    cat(sprintf("    genes: %s\n",
                paste(colnames(X)[sel][1:min(15, length(sel))],
                      collapse = ", ")))
  }
  invisible(fit)
}

diagnose("A: tFDR=0.2 K=100 max_T_stop=TRUE",  tFDR = 0.2, K = 100,
         max_T_stop = TRUE)
diagnose("B: tFDR=0.2 K=100 max_T_stop=FALSE", tFDR = 0.2, K = 100,
         max_T_stop = FALSE)
diagnose("C: tFDR=0.3 K=100 max_T_stop=TRUE",  tFDR = 0.3, K = 100,
         max_T_stop = TRUE)
diagnose("D: tFDR=0.2 K=20  trex+GVS",         tFDR = 0.2, K = 20,
         method = "trex+GVS")
