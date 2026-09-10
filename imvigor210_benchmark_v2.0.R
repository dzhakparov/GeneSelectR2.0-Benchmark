# ==============================================================================
#  GeneSelectR 2.0 — IMvigor210 benchmark (bladder cancer immunotherapy)
# ==============================================================================
#
#  Dataset: Mariathasan et al. 2018, anti-PD-L1 (atezolizumab) in metastatic
#  urothelial carcinoma. Binary outcome: Responder (CR/PR) vs Non-responder
#  (SD/PD).
#
#  This script mirrors sosall_benchmark.R exactly — same methods, same nested
#  CV, same three-model ensemble evaluator, same native-set stability, same
#  biology ablation. Only the data layer and the disease context differ, so the
#  two benchmarks are directly comparable.
#
#  It produces:
#    * Nested CV performance (AUC, Balanced Accuracy, MCC) for each method,
#      evaluated with a three-model soft-voting ensemble (elastic net +
#      XGBoost + random forest) to avoid linear-selector / linear-classifier bias
#    * Native-set stability: Nogueira index (headline, with CIs) and mean
#      pairwise Jaccard (robustness), both on each method's actually-selected set
#    * Nogueira top-k curves (supplementary)
#    * STRING protein-protein interaction (PPI) enrichment
#    * Biology-method ablation (multilayer vs network vs semantic)
#    * Headline figures with a single GeneSelectR point, plus detailed sweeps
#    * A full text log of the run
#
#  Data access (auto-detected, in order):
#    A) easierData (Bioconductor): BiocManager::install("easierData")
#    B) GitHub RDS: data/IMvigor210.all.rds (downloaded if absent)
#
# ==============================================================================


# Null-coalescing operator: returns y if x is NULL or empty, else x.
# Defined early because it's used in the method configuration below.
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x


# ------------------------------------------------------------------------------
#  Error surfacing
# ------------------------------------------------------------------------------
#
#  Ranking methods are wrapped in tryCatch so one failure doesn't kill the whole
#  benchmark. But a silently-swallowed error looks identical to a method that
#  genuinely produces a poor ranking — mRMR failing into the random fallback is
#  indistinguishable from mRMR performing like Random. This helper reports the
#  failure loudly (and records it) so the distinction is visible.

# Records every caught failure: method name -> vector of error messages.
method_failures <- new.env(parent = emptyenv())

report_failure <- function(method_name, stage, error_object) {
  message_text <- conditionMessage(error_object)
  cat(sprintf("  !! %s FAILED at %s: %s\n", method_name, stage, message_text))

  existing <- if (exists(method_name, envir = method_failures)) {
    get(method_name, envir = method_failures)
  } else {
    character(0)
  }
  assign(method_name,
         c(existing, sprintf("[%s] %s", stage, message_text)),
         envir = method_failures)
  NULL
}


# ------------------------------------------------------------------------------
#  Package loading
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GeneSelectR)
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# Install and load auxiliary packages if missing. Some dependencies here are
# Bioconductor packages (edgeR, STRINGdb), which install.packages() cannot
# fetch, so we route those through BiocManager.
bioconductor_packages <- c("edgeR", "STRINGdb", "easierData",
                           "SummarizedExperiment", "ExperimentHub")

ensure_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (package_name %in% bioconductor_packages) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org",
                         quiet = TRUE)
      }
      BiocManager::install(package_name, ask = FALSE, update = FALSE)
    } else {
      install.packages(package_name,
                       repos = "https://cloud.r-project.org",
                       quiet = TRUE)
    }
  }
  suppressPackageStartupMessages(library(package_name, character.only = TRUE))
}

ensure_package("mRMRe")        # mRMR feature selection
ensure_package("Boruta")       # Boruta wrapper around random forest
ensure_package("ranger")       # fast random forest for importance ranking
ensure_package("stabm")        # Nogueira stability index
ensure_package("STRINGdb")     # STRING PPI database client (Bioconductor)
ensure_package("xgboost")      # gradient boosting for the ensemble classifier
ensure_package("edgeR")        # TMM + log-CPM normalisation (Bioconductor)
ensure_package("knockoff")     # model-X knockoff filter (FDR-controlled gate)
ensure_package("ExperimentHub")
ensure_package("SummarizedExperiment")
ensure_package("easierData")   # primary IMvigor210 source (EH6677)


# ------------------------------------------------------------------------------
#  Configuration
# ------------------------------------------------------------------------------
#
#  All tunable parameters are gathered here so the script can be re-run with
#  different settings without searching through the code.

random_seed <- 42
set.seed(random_seed)

# --- Data source ------------------------------------------------------------
# Resolved automatically: easierData if installed, else the GitHub RDS.
imvigor_rds_path <- "data/IMvigor210.all.rds"
imvigor_rds_url  <- "https://github.com/snijeshvp/IMvigor210/raw/main/IMvigor210.all.rds"

# Response coding in the source metadata.
responder_code     <- "R"
nonresponder_code  <- "NR"

# Count filter applied before edgeR normalisation: keep genes with >= 10 counts
# in >= 10 samples.
min_count_per_gene    <- 10
min_samples_per_gene  <- 10

# --- Metadata schema --------------------------------------------------------
# IMvigor210 has no tissue-location confounder analogous to SOS-ALL, so
# residualisation is off by default. If you want to adjust for a covariate
# (e.g. "Immune phenotype" or tissue site), add its column name here and set
# do_residualisation <- TRUE; the residualisation machinery is dataset-agnostic.
confounders_categorical <- character(0)
do_residualisation      <- FALSE

# --- Outer cross-validation -------------------------------------------------
# Stratified K-fold repeated R times. K_OUTER folds per repeat.
k_outer_folds  <- 5
n_outer_repeats <- 3

# Top-k panel sizes evaluated in the parsimony curves.
panel_sizes_evaluated <- c(10, 20, 50, 100, 200, 500)

# --- GeneSelectR settings ---------------------------------------------------
# Number of complementary half-samples for stability selection.
gs_n_subsamples       <- 50
# --- Selection gate ---------------------------------------------------------
# The gate used by EVERY GS variant unless its config overrides gate_method.
#
# PFER is NOT the default and is not used. At n << p the Meinshausen-Buhlmann
# threshold pi_thr = 0.5 + q_bar^2/(2 p PFER) clamps to 1.0, meaning a gene
# would have to be selected in EVERY subsample to pass. None are, so the
# guarantee-unmet fallback fired on every fit and "PFER-controlled" described a
# case this data never reaches.
#
# The model-X knockoff filter replaces it: active at n << p, needs no q_max,
# targets FDR (the right error rate for discovery), and decouples error control
# from the ranking statistic.
#
# Caveats, stated because they belong in the methods section: second-order
# knockoffs have no proof of FDR control (expression is not Gaussian and Sigma
# is estimated), and knockoff power is weakest with weak signal and small n.
# Selecting nothing is a possible honest outcome, not a tuning failure.
#
# Set to "pfer" only to reproduce the legacy behaviour.
gs_gate_method        <- "none"   # no gate: see the note above; panel is top-k
# fdr = 0 is unsatisfiable BY CONSTRUCTION (the FDP numerator is floored by the
# offset, so every threshold returns Inf and nothing passes — seen in the run
# log as "target FDR = 0.00 ... threshold = Inf"). fdr >= 1 accepts every
# discovery as false and provides no control (the package warns as much).
# 0.1–0.2 is the usable range for a discovery panel.
gs_knockoff_fdr       <- 0.2
# 2 draws is too few for the e-value aggregation to stabilise the selected
# set; the package default of 5 exists for that reason.
gs_knockoff_draws     <- 5

# Retained for the legacy gate only; ignored unless gs_gate_method == "pfer".
gs_pfer_bound         <- 2

# q_max caps how many genes each subsample fit may select. It exists ONLY to
# make the PFER threshold attainable -- a hyperparameter chosen to satisfy a
# formula rather than the biology. Under the knockoff gate it serves no purpose
# and actively discards information from the subsample fits that feed the
# utility pillar, so it is disabled.
gs_q_max              <- if (gs_gate_method == "pfer") 50 else NULL
# Permutation GATE settings. No variant currently uses gate_method =
# "permutation" -- the gate tests the global null (Y independent of all of X)
# rather than the conditional null that variable selection needs, so a gene
# merely co-expressed with a true signal can pass it. Knockoffs test the
# correct null and are the gate in use. These knobs are retained only so a
# permutation-gate variant can be re-enabled without re-plumbing; the package
# still supports it.
gs_permutation_n      <- 10
gs_permutation_fdr    <- 0.1
# --- Pillar calibration -----------------------------------------------------
# Every GS variant calibrates its pillars to evidence ratios unless it
# overrides calibration_mode.
#
# Under the default percentile scale each pillar is a RANK within the candidate
# set, so its zero point is "worst gene here" -- a property of the set, not the
# gene. Two consequences, both visible in the SOS-ALL/IMvigor output:
#
#   * Unannotated genes get b = 0 exactly (percentile01 excludes zeros from the
#     ranking), which the geometric mean reads as evidence AGAINST rather than
#     as no evidence. Novel biology cannot surface.
#   * Deeply-annotated genes score high on any similarity measure simply by
#     being close to everything, so the pillar rewards fame. This is how MYH6
#     (cardiac myosin) scored b = 0.85 in a bladder-cancer run.
#
# Evidence ratios fix both: 1 = no evidence (multiplicative identity, no veto),
# and biology is compared against a depth-matched null so annotation depth
# stops being rewarded.
#
# Cost: one permutation loop per (backbone, alpha, fold). Memoised, so variants
# sharing a backbone reuse it rather than each rebuilding an identical null.
gs_calibration_mode        <- "evidence_ratio"
gs_calibration_permutations <- 10
# Minimum panel size for the LEGACY PFER gate only: if fewer than this pass the
# threshold, GeneSelectR falls back to the top gs_min_selected by frequency and
# records pfer_passed = FALSE. The knockoff gate has no such fallback -- it
# selects what clears the FDR threshold, and selecting nothing is a valid
# answer rather than something to pad.
gs_min_selected       <- 20
# --- Subsampling scheme (selects the whole run's regime) ---------------------
# Passed on the command line so the two arms can run CONCURRENTLY into separate
# output folders, letting you inspect one while the other is still going:
#
#     Rscript imvigor210_benchmark.R half     ->  results_imvigor210/<date>_half/
#     Rscript imvigor210_benchmark.R kfold    ->  results_imvigor210/<date>_kfold/
#
# Defaults to "kfold" when sourced interactively with no argument.
#
#   "half"  = stratified 50% draws. Each fit sees half the data and training
#             sets overlap only by chance. This is the regime the
#             Meinshausen-Buhlmann PFER bound is derived under, so it is the
#             only scheme where the error-control guarantee holds.
#
#   "kfold" = stratified K-fold. Each fit sees (K-1)/K of the data (80% at
#             K=5) and is therefore much better conditioned at small n. BUT
#             training sets share >=60% of their samples, so selections are
#             correlated by construction: selection frequency rises partly
#             because fits improve and partly for a trivial reason. The PFER
#             bound is NOT valid here.
#
# Comparing the two runs, the diagnostic is stability AND AUC together:
#   stability up + AUC up   -> real gain from the extra training data
#   stability up + AUC flat -> overlap artifact; do not adopt
command_args <- commandArgs(trailingOnly = TRUE)
gs_subsample_scheme <- if (length(command_args) >= 1 &&
                           nzchar(command_args[1])) {
  command_args[1]
} else {
  "kfold"
}
if (!gs_subsample_scheme %in% c("half", "kfold")) {
  stop("Subsampling scheme must be 'half' or 'kfold', got: ",
       gs_subsample_scheme)
}
gs_subsample_k_folds  <- 5

# Print GeneSelectR's internal step-by-step output (including Open Targets seed
# retrieval and STRING network diagnostics for bio_mode = "network"). Set TRUE
# to debug the biology methods; FALSE keeps the benchmark console readable.
gs_internal_verbose   <- TRUE
# Elastic net mixing parameter grid (tuned within each GS call).
gs_alpha_grid         <- c(0.5, 1.0)
# Parallel cores. The step-1 subsample fits, the knockoff draws, and the
# calibration permutation null all parallelise across cores. Capped at 8 so
# cluster startup does not dominate the small jobs. NOTE: if you run the
# "half" and "kfold" arms CONCURRENTLY on one machine, halve this value so
# the two processes do not oversubscribe the CPU. Set back to 1 when
# debugging (clusters make tracebacks harder to read).
# HARD CAP ON WORKERS -- 7, set by the user, not negotiable by any formula.
# This script has no split-level parallelism; n_parallel_cores is handed to
# geneselectr2_fit, which spawns that many PSOCK workers INSIDE each fit. Those
# are real R sessions and they are what fills RAM, so the same cap applies.
# It was previously 8. Raise only on explicit instruction, after measuring
# per-worker RSS with `ps`.
MAX_PARALLEL_WORKERS <- 7L

.detected_cores <- parallel::detectCores(logical = FALSE)
n_parallel_cores <- if (is.na(.detected_cores)) 1L else {
  max(1L, min(MAX_PARALLEL_WORKERS, .detected_cores - 1L))
}
rm(.detected_cores)

# The cap above is per-process, so it does not constrain two runs started in
# separate terminals -- three concurrent runs took 7 each and crashed the
# machine. gs_claim_workers() takes only the share of the 7 that other live
# runs are not already holding, and refuses to start when none is left.
source("benchmarks/worker_budget.R")
n_parallel_cores <- gs_claim_workers(n_parallel_cores, label = "IMvigor210")

# --- Variance pre-filter ----------------------------------------------------
# Keep the top-N most variable genes before any analysis to bound dimensionality.
.p_override <- Sys.getenv("GS_TOP_VARIABLE_GENES", "")
top_variable_genes <- if (nzchar(.p_override)) as.integer(.p_override) else 2000L
if (length(top_variable_genes) != 1L || !is.finite(top_variable_genes) ||
    top_variable_genes < 1L) {
  stop("GS_TOP_VARIABLE_GENES must be one positive integer.")
}
rm(.p_override)

# Determine the candidate genes within each training fold. A global variance
# filter uses test-fold expression values and is retained only for a labelled
# sensitivity analysis of historical runs.
variance_filter_scope <- Sys.getenv("GS_VARIANCE_FILTER_SCOPE", "train")
if (!variance_filter_scope %in% c("train", "global")) {
  stop("GS_VARIANCE_FILTER_SCOPE must be 'train' or 'global'.")
}
if (!identical(variance_filter_scope, "train")) {
  stop("IMvigor210 requires GS_VARIANCE_FILTER_SCOPE=train because count ",
       "filtering and TMM normalization are fitted within each training fold.")
}

# --- Biology settings -------------------------------------------------------
# GO terms used by supervised biology modes (immune response and related).
# Terms centre on T-cell mediated immunity and immune checkpoint biology, which
# is the mechanism anti-PD-L1 therapy acts through.
target_go_terms <- c(
  "GO:0006955",  # immune response
  "GO:0002376",  # immune system process
  "GO:0002250",  # adaptive immune response
  "GO:0042110",  # T cell activation
  "GO:0050863",  # regulation of T cell activation
  "GO:0002682",  # regulation of immune system process
  "GO:0050776",  # regulation of immune response
  "GO:0031341"   # regulation of cell killing
)
disease_term <- "bladder carcinoma"

# Glmnet alpha values searched by the downstream classifier (independent of GS).
glmnet_alpha_grid <- c(0.5, 1.0)

# --- Stability / quality analyses -------------------------------------------
# Panel sizes at which the Nogueira index is reported.
nogueira_panel_sizes  <- c(10, 20, 50, 100, 200)
# Subsamples used for the shared Nogueira / Jaccard cache.
shared_n_subsamples   <- 50
# Maximum top-k stored per subsample (must cover the largest Nogueira k).
shared_max_k          <- max(nogueira_panel_sizes)

# --- STRING database --------------------------------------------------------
# Permutations for the STRING enrichment null. 1000 gives a p-value floor of
# ~0.001; 50 gave ~0.02, which every real method saturated at.
string_n_permutations <- 1000
string_version    <- "12.0"
string_species_id <- 9606  # Homo sapiens

# --- STRING download timeout and cache --------------------------------------
# R's default download timeout is 60s and STRING's human interaction file is
# 79 MB. A cut-off download both fails the analysis and leaves a truncated file
# behind that poisons the next attempt. Raise it, and keep the ~100 MB of
# reference files in one place instead of refetching into a fresh tempdir.
options(timeout = max(3600, getOption("timeout")))
string_cache_dir <- file.path("data", "string_db_cache")
dir.create(string_cache_dir, recursive = TRUE, showWarnings = FALSE)

# Keep the primary Bioconductor input in the project so its downloaded object
# and metadata database are retained with the benchmark provenance.
experimenthub_cache_dir <- file.path("data", "experimenthub_cache")
dir.create(experimenthub_cache_dir, recursive = TRUE, showWarnings = FALSE)

# Point GeneSelectR's own network-biology layer at the SAME cache. Without
# this, bio_mode = "network" re-downloads STRING on its own and, before the
# package was fixed, silently returned a constant pillar when that download
# was cut off.
options(GeneSelectR.string_cache = normalizePath(string_cache_dir))


# --- Output paths -----------------------------------------------------------
run_date   <- format(Sys.Date(), "%Y-%m-%d")
run_stamp  <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
.output_root <- Sys.getenv("GS_OUTPUT_ROOT", "")
output_dir <- if (nzchar(.output_root)) {
  file.path(.output_root, "imvigor210",
            sprintf("%s_%s_%sfilter", run_date, gs_subsample_scheme,
                    variance_filter_scope))
} else {
  file.path("results_imvigor210",
            sprintf("%s_%s_%sfilter", run_date, gs_subsample_scheme,
                    variance_filter_scope))
}
rm(.output_root)
dir.create(file.path(output_dir, "data"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "logs"),    recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
#  Run logging
# ------------------------------------------------------------------------------
#
#  Everything printed to the console (cat, print, message via the warning
#  handler below) is also teed to a timestamped log file, so each run leaves a
#  complete record: configuration, per-step diagnostics, timings, results
#  tables, and any method failures. sink(split = TRUE) mirrors output rather
#  than redirecting it, so the console stays live while the file fills.

log_file <- file.path(output_dir, "logs",
                      sprintf("run_%s.log", run_stamp))
log_connection <- file(log_file, open = "wt")

# Tee normal output (stdout) to the log while keeping it on screen.
sink(log_connection, split = TRUE)
# Also capture warnings/messages (stderr) into the same file. This is NOT
# split, so warnings go to the file; they still surface on the console at the
# end of the run via R's normal warning collection.
sink(log_connection, type = "message")

# Guarantee the sinks are released and the file closed even if the script
# errors out partway through — otherwise the log is left locked and truncated.
close_log <- function() {
  # CAREFUL: the two sink.number() calls return DIFFERENT KINDS of value.
  #
  #   sink.number()                 -> a COUNT of diverted output connections.
  #                                    Reaches 0 when none are active.
  #   sink.number(type = "message") -> the CONNECTION NUMBER handling messages.
  #                                    It is 2 (stderr) when NO message sink is
  #                                    active, and never returns 0.
  #
  # Looping "while (sink.number(type = 'message') > 0)" therefore spins
  # forever: popping the sink leaves it at 2, which is still > 0. Check
  # against 2, not 0.
  if (sink.number(type = "message") != 2) {
    sink(type = "message")
  }
  while (sink.number() > 0) sink()

  # tryCatch: the connection may already be closed if this runs twice (normal
  # exit plus the error handler), and close() on a closed connection errors.
  tryCatch({
    if (isOpen(log_connection)) close(log_connection)
  }, error = function(e) invisible(NULL))
}

# on.exit() only works inside a function; at the top level of a sourced script
# it has no effect. Instead register a global error option so that if the run
# aborts, the log is flushed and closed before R prints the error. The normal
# completion path calls close_log() explicitly at the very end.
options(GeneSelectR_prev_error = getOption("error"))
options(error = function() {
  message("!! Run aborted by an error — flushing and closing log.")
  try(close_log(), silent = TRUE)
})

# Wall-clock start, used for total-runtime reporting at the end.
benchmark_start_time <- Sys.time()

cat(paste(rep("=", 78), collapse = ""), "\n")
cat("GeneSelectR IMvigor210 benchmark run\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
cat(sprintf("Timestamp:        %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Log file:         %s\n", log_file))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("R version:        %s\n", R.version.string))
cat(sprintf("Platform:         %s\n", R.version$platform))
cat(sprintf("Hostname:         %s\n",
            tryCatch(Sys.info()[["nodename"]], error = function(e) "unknown")))
cat("\n--- Configuration ---\n")
cat(sprintf("  random_seed        = %d\n", random_seed))
cat(sprintf("  top_variable_genes = %d\n", top_variable_genes))
cat(sprintf("  gs_n_subsamples    = %d\n", gs_n_subsamples))
cat(sprintf("  gs_pfer_bound      = %g\n", gs_pfer_bound))
cat(sprintf("  gs_q_max           = %d\n", gs_q_max))
cat(sprintf("  gs_min_selected    = %d\n", gs_min_selected))
cat(sprintf("  SUBSAMPLE SCHEME   = %s%s\n", gs_subsample_scheme,
            if (gs_subsample_scheme == "kfold")
              sprintf(" (K=%d; PFER bound VOID under this scheme)",
                      gs_subsample_k_folds)
            else " (50%% draws; PFER bound valid)"))
cat(sprintf("  gs_internal_verbose= %s\n", gs_internal_verbose))
cat(sprintf("  n_parallel_cores   = %d\n", n_parallel_cores))
cat(sprintf("  panel sizes        = %s\n",
            paste(panel_sizes_evaluated, collapse = ", ")))
cat(sprintf("  disease term       = %s\n", disease_term))
cat("  (GS configurations listed once defined, below)\n")
cat("\n")


# ------------------------------------------------------------------------------
#  Performance metric helpers
# ------------------------------------------------------------------------------

# Fast AUC computation using rank-based formula. Avoids the overhead of pROC.
compute_auc <- function(true_labels, predicted_scores) {
  true_labels   <- droplevels(true_labels)
  positive_idx  <- which(true_labels == levels(true_labels)[2])
  negative_idx  <- which(true_labels == levels(true_labels)[1])

  if (length(positive_idx) == 0 || length(negative_idx) == 0) {
    return(NA_real_)
  }

  # AUC = (sum of ranks for positives - tie correction) / (n_pos * n_neg)
  rank_scores <- rank(predicted_scores, ties.method = "average")
  numerator   <- sum(rank_scores[positive_idx]) -
    length(positive_idx) * (length(positive_idx) + 1) / 2

  numerator / (length(positive_idx) * length(negative_idx))
}

# Threshold-based metrics that complement AUC.
# Balanced accuracy handles class imbalance better than raw accuracy.
# MCC (Matthews correlation coefficient) is informative even with imbalance.
compute_classification_metrics <- function(true_labels,
                                           predicted_probabilities,
                                           threshold = 0.5) {
  true_labels    <- droplevels(true_labels)
  positive_class <- levels(true_labels)[2]
  negative_class <- levels(true_labels)[1]

  predicted_labels <- ifelse(predicted_probabilities >= threshold,
                             positive_class, negative_class)
  predicted_labels <- factor(predicted_labels, levels = levels(true_labels))

  # Confusion matrix counts
  true_positives  <- sum(predicted_labels == positive_class & true_labels == positive_class)
  true_negatives  <- sum(predicted_labels == negative_class & true_labels == negative_class)
  false_positives <- sum(predicted_labels == positive_class & true_labels == negative_class)
  false_negatives <- sum(predicted_labels == negative_class & true_labels == positive_class)

  sensitivity <- if ((true_positives + false_negatives) > 0)
    true_positives / (true_positives + false_negatives) else 0
  specificity <- if ((true_negatives + false_positives) > 0)
    true_negatives / (true_negatives + false_positives) else 0

  balanced_accuracy <- (sensitivity + specificity) / 2

  # MCC denominator can underflow; guard against it.
  mcc_denominator <- sqrt(as.numeric(true_positives + false_positives) *
                            (true_positives + false_negatives) *
                            (true_negatives + false_positives) *
                            (true_negatives + false_negatives))
  mcc <- if (mcc_denominator > 0)
    (true_positives * true_negatives - false_positives * false_negatives) /
    mcc_denominator
  else 0

  list(balanced_accuracy = balanced_accuracy, mcc = mcc)
}


# ------------------------------------------------------------------------------
#  Confounder residualisation
# ------------------------------------------------------------------------------
#
#  To remove the effect of tissue location (skin vs blood etc.), we regress
#  the expression matrix on a design matrix of confounders, fit on the training
#  data only, then apply the same fit to the held-out test samples. This
#  prevents leakage that would occur if confounders were residualised globally.

# Build a one-hot design matrix from categorical columns of a metadata data frame.
# Returns both the matrix and the reference levels so the same encoding can be
# applied to held-out samples.
build_design_matrix <- function(metadata, reference_levels = NULL) {
  encoded_matrices  <- list()
  recorded_levels   <- list()

  for (confounder in confounders_categorical) {
    raw_values <- metadata[[confounder]]
    raw_values[is.na(raw_values) | raw_values == ""] <- "NA_level"

    if (is.null(reference_levels)) {
      # Training: derive levels from the data.
      levels_here <- levels(factor(raw_values))
      encoded     <- model.matrix(~ factor(raw_values, levels = levels_here))
      colnames(encoded) <- gsub("^factor\\(raw_values, levels = levels_here\\)",
                                confounder, colnames(encoded))
      if (ncol(encoded) > 1) encoded <- encoded[, -1, drop = FALSE]  # drop intercept

      encoded_matrices[[confounder]] <- encoded
      recorded_levels[[confounder]]  <- levels_here

    } else {
      # Test: enforce the training levels so columns align.
      encoded <- model.matrix(~ factor(raw_values,
                                       levels = reference_levels[[confounder]]))
      colnames(encoded) <- gsub("^factor\\(raw_values, levels = reference_levels\\[\\[confounder\\]\\]\\)",
                                confounder, colnames(encoded))
      if (ncol(encoded) > 1) encoded <- encoded[, -1, drop = FALSE]
      encoded_matrices[[confounder]] <- encoded
    }
  }

  combined <- if (length(encoded_matrices) > 0) {
    do.call(cbind, encoded_matrices)
  } else {
    matrix(0, nrow = nrow(metadata), ncol = 0)
  }

  # Always include an intercept column.
  combined <- cbind("(Intercept)" = 1, combined)

  if (is.null(reference_levels)) {
    list(design_matrix = combined, reference_levels = recorded_levels)
  } else {
    combined
  }
}

# Fit and apply: subtract the projection of expression onto the confounder space.
apply_residualisation <- function(design_matrix, expression_matrix,
                                  coefficients) {
  expression_matrix - design_matrix %*% coefficients
}


# ------------------------------------------------------------------------------
#  Cross-validation fold construction
# ------------------------------------------------------------------------------

# Stratified K-fold splits that preserve class balance in every fold.
make_stratified_folds <- function(outcome_factor, k_folds, seed) {
  set.seed(seed)
  outcome_factor <- droplevels(outcome_factor)

  class1_indices <- which(outcome_factor == levels(outcome_factor)[1])
  class2_indices <- which(outcome_factor == levels(outcome_factor)[2])

  class1_folds <- split(sample(class1_indices),
                        rep(1:k_folds, length.out = length(class1_indices)))
  class2_folds <- split(sample(class2_indices),
                        rep(1:k_folds, length.out = length(class2_indices)))

  # Test indices per fold are the union of both class fold assignments.
  lapply(1:k_folds, function(fold) {
    sort(c(class1_folds[[fold]], class2_folds[[fold]]))
  })
}


# ------------------------------------------------------------------------------
#  Downstream classifier: three-model soft-voting ensemble
# ------------------------------------------------------------------------------
#
#  A reviewer noted that evaluating linearly-selected features with a linear
#  classifier biases the comparison: genes chosen for linear separability
#  look good on a linear model by construction. To neutralise this, we
#  evaluate every gene panel with an ensemble of three classifiers drawn
#  from three different model families, each with a distinct inductive bias:
#
#    * Elastic net   (linear, L1/L2-regularised)
#    * XGBoost       (gradient-boosted trees)
#    * Random forest (bagged trees)
#
#  Their predicted class probabilities are averaged ("soft voting"). No single
#  family's assumptions dominate, so no feature-selection method gains an
#  advantage from sharing assumptions with the evaluator. The same ensemble
#  is applied identically to every method's panel.
#
#  Each component returns a probability vector for the positive class, or a
#  vector of NA if that component fails to train. The ensemble averages over
#  whichever components succeeded.

# --- Component 1: elastic net with alpha tuning -----------------------------
predict_component_glmnet <- function(train_features, train_labels,
                                     test_features) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  n_inner_folds  <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))

  best_fit <- NULL
  best_auc <- -Inf

  for (alpha_value in glmnet_alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(train_features, numeric_labels,
                        family = "binomial",
                        alpha = alpha_value,
                        nfolds = n_inner_folds),
      error = function(e) NULL
    )
    # cv.glmnet stores -AUC in $cvm; -min(cvm) is the best AUC across alphas.
    if (!is.null(fit) && -min(fit$cvm) > best_auc) {
      best_auc <- -min(fit$cvm)
      best_fit <- fit
    }
  }

  if (is.null(best_fit)) return(rep(NA_real_, nrow(test_features)))

  as.numeric(predict(best_fit, test_features,
                     s = "lambda.min", type = "response"))
}

# --- Component 2: XGBoost (gradient-boosted trees) --------------------------
predict_component_xgboost <- function(train_features, train_labels,
                                      test_features) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])

  train_matrix <- xgboost::xgb.DMatrix(data = as.matrix(train_features),
                                       label = numeric_labels)

  # Fixed, sensible defaults. We avoid heavy tuning so XGBoost is neither
  # advantaged nor handicapped relative to the other components; shallow
  # trees and a modest round count suit the small sample sizes here.
  params <- list(
    objective        = "binary:logistic",
    eval_metric      = "auc",
    eta              = 0.1,
    max_depth        = 3,
    subsample        = 0.8,
    colsample_bytree = 0.8
  )

  fit <- tryCatch(
    xgboost::xgb.train(params = params, data = train_matrix,
                       nrounds = 100, verbose = 0),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(test_features)))

  predictions <- tryCatch(
    predict(fit, as.matrix(test_features)),
    error = function(e) rep(NA_real_, nrow(test_features))
  )
  predictions
}

# --- Component 3: random forest (bagged trees) ------------------------------
predict_component_random_forest <- function(train_features, train_labels,
                                            test_features) {
  train_labels   <- droplevels(train_labels)
  positive_class <- levels(train_labels)[2]

  fit <- tryCatch(
    ranger::ranger(x = train_features, y = train_labels,
                   num.trees = 500, probability = TRUE,
                   seed = random_seed, num.threads = 1),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(test_features)))

  predictions <- tryCatch({
    probability_matrix <- predict(fit, data = test_features)$predictions
    probability_matrix[, positive_class]
  }, error = function(e) rep(NA_real_, nrow(test_features)))

  predictions
}

# --- Ensemble: average the probabilities of all successful components -------
#
#  Returns EVERY component's probabilities alongside the soft-voting average,
#  not just the average. The three components are trained regardless, so their
#  individual predictions are already paid for -- discarding them and reporting
#  only the ensemble would throw away the ablation for free.
#
#  This makes the evaluator itself testable. Two questions it answers:
#
#    1. Does soft voting actually beat its own components? An ensemble is
#       assumed to help; on n ~ 100 with k = 10 features, XGBoost with 100
#       rounds is likely overfitting, and averaging a bad model into a good one
#       can hurt. Worth checking rather than assuming.
#
#    2. Do the method rankings depend on the evaluator? If GeneSelectR looks
#       better under a linear evaluator and worse under a tree evaluator, that
#       is a selector/evaluator interaction, not a property of the selector --
#       exactly the bias the ensemble exists to avoid. Reporting per-evaluator
#       metrics makes that visible instead of hiding it inside an average.
#
#  Cost: none. Same three fits, more outputs.
#
#  @return Named list of probability vectors: glmnet, xgboost, rf, ensemble.
predict_with_ensemble <- function(train_features, train_labels, test_features) {

  glmnet_probabilities <- predict_component_glmnet(
    train_features, train_labels, test_features)
  xgboost_probabilities <- predict_component_xgboost(
    train_features, train_labels, test_features)
  rf_probabilities <- predict_component_random_forest(
    train_features, train_labels, test_features)

  # Stack the three probability vectors as columns; average row-wise over
  # whichever components produced non-NA predictions.
  probability_components <- cbind(glmnet_probabilities,
                                  xgboost_probabilities,
                                  rf_probabilities)

  ensemble_probabilities <- rowMeans(probability_components, na.rm = TRUE)

  # If every component failed for a sample, rowMeans returns NaN; convert
  # those back to NA so the metric functions treat them as missing.
  ensemble_probabilities[is.nan(ensemble_probabilities)] <- NA_real_

  list(
    glmnet   = glmnet_probabilities,
    xgboost  = xgboost_probabilities,
    rf       = rf_probabilities,
    ensemble = ensemble_probabilities
  )
}


# ==============================================================================
#  GeneSelectR ranking wrapper
# ==============================================================================
#
#  All four GS configurations differ only in their bio_mode and active
#  components. We define a single parameterised ranker and instantiate the
#  variants with different argument lists.

rank_with_geneselectr <- function(train_features, train_labels, gs_config,
                                  gs_method_label = "GeneSelectR") {

  best_result <- NULL
  best_auc    <- -Inf

  # Try each alpha in the grid and keep the configuration with the best
  # internal CV AUC. This is GeneSelectR's own out-of-bag estimate, not the
  # outer CV evaluation.
  for (alpha_value in gs_config$alpha_grid) {

    fit_attempt <- tryCatch(
      geneselectr2_fit(
        train_features, train_labels,
        # --- Selection backbone -------------------------------------------
        selection_method      = "stability_selection",
        B                     = gs_n_subsamples,
        gate_method           = gs_config$gate_method %||% gs_gate_method,
        knockoff_fdr          = gs_knockoff_fdr,
        knockoff_draws        = gs_knockoff_draws,
        permutation_n         = gs_permutation_n,
        permutation_fdr       = gs_permutation_fdr,
        pfer                  = gs_pfer_bound,
        q_max                 = gs_q_max,
        min_selected          = gs_min_selected,
        subsample_scheme      = gs_config$subsample_scheme %||% gs_subsample_scheme,
        subsample_k_folds     = gs_subsample_k_folds,
        regularization_method = gs_config$regularization_method %||% "elastic_net",
        alpha                 = alpha_value,
        # --- Score components --------------------------------------------
        utility_method        = gs_config$utility_method,
        components            = gs_config$components,
        score_formula         = gs_config$score_formula %||% "geometric",
        calibration_mode      = gs_config$calibration_mode %||% gs_calibration_mode,
        calibration_n_permutations = gs_calibration_permutations,
        # --- Biology -----------------------------------------------------
        bio_mode              = gs_config$bio_mode,
        bio_layers            = gs_config$bio_layers %||% c("go"),
        target_terms          = target_go_terms,
        bio_disease_term      = disease_term,
        bio_string_threshold  = gs_bio_defaults$string_threshold,
        bio_rwr_restart       = gs_bio_defaults$rwr_restart,
        bio_max_seeds         = gs_bio_defaults$max_seeds,
        bio_seed_min_score    = gs_bio_defaults$seed_min_score,
        # --- Infrastructure ----------------------------------------------
        n_cores               = n_parallel_cores,
        random_seed           = random_seed,
        verbose               = gs_internal_verbose
      ),
      error = function(e) report_failure(gs_method_label, "geneselectr2_fit", e)
    )

    if (!is.null(fit_attempt) &&
        fit_attempt$cv_results$mean_auc > best_auc) {
      best_auc    <- fit_attempt$cv_results$mean_auc
      best_result <- fit_attempt
    }
  }

  if (is.null(best_result)) {
    # Fallback: random ranking so downstream code doesn't crash.
    return(list(ranked   = sample(colnames(train_features)),
                selected = character(0),
                gs_object = NULL))
  }

  # GeneSelectR's natural "selected" set is genes that pass RENT.
  # These are the genes the framework considers reliable given the
  # stability, coefficient significance, and sign consistency criteria.
  # GeneSelectR's selected set is now the PFER-thresholded genes, exposed as
  # the logical `selected` column (RENT was removed from the method).
  selected_flag <- best_result$gene_scores$selected
  if (is.null(selected_flag)) selected_flag <- rep(FALSE, nrow(best_result$gene_scores))

  selected_genes <- best_result$gene_scores$gene[selected_flag]

  list(ranked    = best_result$gene_scores$gene,
       selected  = selected_genes,
       gs_object = best_result)
}


# ==============================================================================
#  Competitor ranking functions
# ==============================================================================
#
#  Each ranker takes train features and labels and returns a list with two
#  elements: `ranked` (character vector of genes ordered from best to worst)
#  and `gs_object` (NULL for non-GS methods — kept for interface uniformity).

# Differential expression (paper review fix: use -log10(p) * |t-statistic|
# instead of just |t-statistic| to properly weight by significance).
# Selected set = genes with BH-adjusted p < 0.05 (the natural sparsity
# criterion for differential expression analysis).
rank_by_differential_expression <- function(train_features, train_labels) {
  train_labels <- droplevels(train_labels)

  # Compute scores and p-values in two explicit numeric vectors. We do NOT
  # use apply() here: t.test()$statistic carries a "t" name, and when that
  # named value flows into c(score = ...), R builds the name "score.t",
  # which corrupts the row names of the assembled matrix and breaks the
  # downstream per_gene_results["score", ] lookup. Explicit vectors avoid it.
  n_genes      <- ncol(train_features)
  scores       <- numeric(n_genes)
  raw_p_values <- numeric(n_genes)
  names(scores)       <- colnames(train_features)
  names(raw_p_values) <- colnames(train_features)

  for (gene_idx in seq_len(n_genes)) {
    t_test <- tryCatch(t.test(train_features[, gene_idx] ~ train_labels),
                       error = function(e) NULL)
    if (is.null(t_test)) {
      scores[gene_idx]       <- 0
      raw_p_values[gene_idx] <- 1
    } else {
      # Strip the "t" name off the statistic before arithmetic.
      t_statistic            <- unname(t_test$statistic)
      scores[gene_idx]       <- -log10(max(t_test$p.value, 1e-300)) *
        abs(t_statistic)
      raw_p_values[gene_idx] <- unname(t_test$p.value)
    }
  }

  # Selected set = genes with BH-adjusted p < 0.05 (the natural sparsity
  # criterion for differential expression analysis).
  adjusted_p     <- p.adjust(raw_p_values, method = "BH")
  selected_genes <- names(raw_p_values)[adjusted_p < 0.05]

  list(ranked   = names(sort(scores, decreasing = TRUE)),
       selected = selected_genes,
       gs_object = NULL)
}

# LASSO with internal CV for lambda selection.
# Selected set = genes with non-zero coefficient (the natural sparsity
# criterion for L1-regularised regression).
rank_by_lasso <- function(train_features, train_labels) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  n_inner_folds  <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))

  fit <- tryCatch(
    glmnet::cv.glmnet(train_features, numeric_labels,
                      family = "binomial", alpha = 1, nfolds = n_inner_folds),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(ranked   = sample(colnames(train_features)),
                selected = character(0), gs_object = NULL))
  }

  # First coefficient is the intercept; drop it.
  coefficients <- as.numeric(coef(fit, s = "lambda.min"))[-1]
  names(coefficients) <- colnames(train_features)

  selected_genes <- names(coefficients)[coefficients != 0]

  list(ranked   = names(sort(abs(coefficients), decreasing = TRUE)),
       selected = selected_genes,
       gs_object = NULL)
}

# Elastic net with alpha tuning across the standard grid.
# Selected set = genes with non-zero coefficient at the best alpha/lambda.
rank_by_elastic_net <- function(train_features, train_labels) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  n_inner_folds  <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))

  best_fit       <- NULL
  best_cv_error  <- Inf

  for (alpha_value in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    fit <- tryCatch(
      glmnet::cv.glmnet(train_features, numeric_labels,
                        family = "binomial", alpha = alpha_value,
                        nfolds = n_inner_folds),
      error = function(e) NULL
    )
    if (!is.null(fit) && min(fit$cvm) < best_cv_error) {
      best_cv_error <- min(fit$cvm)
      best_fit      <- fit
    }
  }

  if (is.null(best_fit)) {
    return(list(ranked   = sample(colnames(train_features)),
                selected = character(0), gs_object = NULL))
  }

  coefficients <- as.numeric(coef(best_fit, s = "lambda.min"))[-1]
  names(coefficients) <- colnames(train_features)

  selected_genes <- names(coefficients)[coefficients != 0]

  list(ranked   = names(sort(abs(coefficients), decreasing = TRUE)),
       selected = selected_genes,
       gs_object = NULL)
}

# Minimum-redundancy maximum-relevance feature selection.
# Capped at 5000 features to keep mRMRe tractable.
rank_by_mrmr <- function(train_features, train_labels) {
  numeric_labels <- as.numeric(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  max_features <- min(ncol(train_features), 5000)

  feature_subset <- train_features[, 1:max_features, drop = FALSE]
  original_names <- colnames(feature_subset)

  # mRMRe requires every column to be numeric (not integer) and will silently
  # mangle non-syntactic column names via make.names(). We therefore build the
  # data frame with safe placeholder names and map back afterwards, and coerce
  # everything to double.
  safe_names <- paste0("V", seq_len(max_features))
  input_df <- data.frame(
    outcome = as.numeric(numeric_labels),
    matrix(as.numeric(feature_subset), nrow = nrow(feature_subset)),
    stringsAsFactors = FALSE
  )
  colnames(input_df) <- c("outcome", safe_names)

  random_fallback <- list(ranked = sample(colnames(train_features)),
                          selected = NULL, gs_object = NULL)

  mrmr_data <- tryCatch(mRMRe::mRMR.data(data = input_df),
                        error = function(e) report_failure("mRMR", "mRMR.data", e))
  if (is.null(mrmr_data)) return(random_fallback)

  mrmr_result <- tryCatch(
    mRMRe::mRMR.classic(data = mrmr_data,
                        target_indices = 1,
                        feature_count = min(200, max_features)),
    error = function(e) report_failure("mRMR", "mRMR.classic", e)
  )
  if (is.null(mrmr_result)) return(random_fallback)

  selected_indices <- tryCatch(
    as.integer(mRMRe::solutions(mrmr_result)[[1]]),
    error = function(e) report_failure("mRMR", "solutions", e)
  )
  if (is.null(selected_indices) || length(selected_indices) == 0) {
    return(random_fallback)
  }

  # solutions() indexes into input_df, whose column 1 is the outcome.
  # Drop the outcome and anything out of range, then map back to gene names.
  selected_indices <- selected_indices[!is.na(selected_indices) &
                                         selected_indices >= 2 &
                                         selected_indices <= (max_features + 1)]
  if (length(selected_indices) == 0) return(random_fallback)

  selected_genes <- original_names[selected_indices - 1L]

  remaining_genes <- setdiff(colnames(train_features), selected_genes)
  list(ranked   = c(selected_genes, remaining_genes),
       selected = selected_genes,
       gs_object = NULL)
}

# Boruta wrapper around random forest (full strength settings).
# Selected set = features confirmed by Boruta (the natural selection criterion).
rank_by_boruta <- function(train_features, train_labels) {
  train_labels <- droplevels(train_labels)
  max_features <- min(ncol(train_features), 5000)

  if (ncol(train_features) > max_features) {
    variances    <- apply(train_features, 2, var)
    top_indices  <- order(variances, decreasing = TRUE)[1:max_features]
    train_features <- train_features[, top_indices, drop = FALSE]
  }

  boruta_result <- tryCatch(
    # Boruta forwards ... to ranger; passing num.trees here collides with
    # Boruta's own argument handling and errors. Boruta's ranger default is
    # already 500 trees, so we simply omit it. maxRuns controls the number of
    # Boruta iterations, which is the knob that actually matters here.
    Boruta::Boruta(x = train_features, y = train_labels,
                   doTrace = 0, maxRuns = 50),
    error = function(e) report_failure("Boruta", "Boruta", e)
  )
  if (is.null(boruta_result)) {
    return(list(ranked   = sample(colnames(train_features)),
                selected = NULL, gs_object = NULL))
  }

  boruta_result <- tryCatch(Boruta::TentativeRoughFix(boruta_result),
                            error = function(e) boruta_result)
  importance_stats <- Boruta::attStats(boruta_result)
  importance_stats <- importance_stats[order(-importance_stats$meanImp), ]

  # Boruta classifies each feature as Confirmed, Tentative, or Rejected.
  # "Confirmed" is the natural selection set. On small subsamples with many
  # features Boruta often confirms NOTHING — that is a genuine property of the
  # method here, not an error, and it means Boruta has no native selected set
  # to compute stability on. We return an empty character vector (not NULL) so
  # downstream code can distinguish "ran, selected nothing" from "never ran".
  selected_genes <- rownames(importance_stats)[
    importance_stats$decision == "Confirmed"
  ]

  remaining_genes <- setdiff(colnames(train_features),
                             rownames(importance_stats))
  list(ranked   = c(rownames(importance_stats), remaining_genes),
       selected = selected_genes,
       gs_object = NULL)
}

# Random forest variable importance via ranger (faster than randomForest).
# No natural sparsity threshold — selected = NULL means the downstream
# Nogueira analysis will fall back to top-k slicing of the ranking.
rank_by_random_forest <- function(train_features, train_labels) {
  train_labels <- droplevels(train_labels)
  max_features <- min(ncol(train_features), 5000)

  if (ncol(train_features) > max_features) {
    variances     <- apply(train_features, 2, var)
    top_indices   <- order(variances, decreasing = TRUE)[1:max_features]
    train_features <- train_features[, top_indices, drop = FALSE]
  }

  rf_fit <- tryCatch(
    ranger::ranger(x = train_features, y = train_labels,
                   num.trees = 1000, importance = "impurity",
                   seed = random_seed),
    error = function(e) NULL
  )
  if (is.null(rf_fit)) {
    return(list(ranked   = sample(colnames(train_features)),
                selected = NULL, gs_object = NULL))
  }

  sorted_importance <- sort(rf_fit$variable.importance, decreasing = TRUE)
  remaining_genes   <- setdiff(colnames(train_features),
                               names(sorted_importance))
  list(ranked   = c(names(sorted_importance), remaining_genes),
       selected = NULL,
       gs_object = NULL)
}

# Random baseline. No natural sparsity — selected = NULL.
rank_at_random <- function(train_features, train_labels) {
  list(ranked = sample(colnames(train_features)),
       selected = NULL, gs_object = NULL)
}


# ==============================================================================
#  Method registry
# ==============================================================================
#
#  GS configurations are kept here as a list so adding a new variant requires
#  only one new entry. The variants span the biology methods we want to
#  compare (multilayer, network propagation, semantic coherence x2) plus the
#  no-biology baselines for ablation. All share the identical stability +
#  utility backbone, so they differ only in the biology pillar and therefore
#  in the final ranking — not in the selected set (which the backbone fixes).

# Default biology knobs reused across configs.
gs_bio_defaults <- list(
  string_threshold = 400,
  rwr_restart      = 0.5,
  max_seeds        = 100,
  seed_min_score   = 0.1
)

gs_configurations <- list(

  # --- Biology methods (the comparison of interest) ---------------------

  GS_multilayer = list(
    bio_mode       = "multilayer",
    bio_layers     = c("go", "msigdb"),
    components     = c("stability", "utility", "bio"),
    utility_method = "instance_shap",
    # No source-matched evidence-ratio null exists for this composite source.
    calibration_mode = "percentile",
    alpha_grid     = gs_alpha_grid
  ),

  GS_network = list(
    bio_mode       = "network",
    components     = c("stability", "utility", "bio"),
    utility_method = "instance_shap",
    # STRING network scores currently support within-run percentile scaling.
    calibration_mode = "percentile",
    alpha_grid     = gs_alpha_grid
  ),

  # Semantic similarity to the disease's target GO terms. (The former
  # "to_set" variant was removed: scoring genes against the selected set is
  # circular. Panel coherence is now a post-hoc evaluation metric instead.)
  GS_semantic = list(
    bio_mode       = "semantic",
    components     = c("stability", "utility", "bio"),
    utility_method = "instance_shap",
    alpha_grid     = gs_alpha_grid
  ),

  # --- No-biology baselines (ablation) ----------------------------------

  GS_stab_util = list(
    bio_mode       = "none",
    components     = c("stability", "utility"),
    utility_method = "instance_shap",
    alpha_grid     = gs_alpha_grid
  ),

  GS_stab_only = list(
    bio_mode       = "none",
    components     = "stability",
    utility_method = "global_coef",
    alpha_grid     = gs_alpha_grid
  ),

  # --- Knockoff gate (fixes the inactive-PFER problem) ------------------
  # GS_knockoff scores utility x biology over an FDR-controlled gate. Stability
  # is deliberately NOT a scoring component here: it was anti-predictive in the
  # ablation, and selection frequency is correlated with SHAP by construction
  # (both derive from the same fits), so multiplying them double-counted one
  # evidence source. Under this gate the two scores come from genuinely
  # distinct sources -- held-out prediction and external databases.
  GS_no_stability = list(
    bio_mode       = "semantic",
    components     = c("utility", "bio"),
    utility_method = "instance_shap",
    alpha_grid     = gs_alpha_grid
  ),

  # Same gate, no biology: isolates what the biology pillar contributes once
  # the gate is sound.
  GS_utility_only = list(
    bio_mode       = "none",
    components     = c("utility"),
    utility_method = "instance_shap",
    alpha_grid     = gs_alpha_grid
  ),

  # --- Uncalibrated comparator ------------------------------------------
  # Identical to GS_semantic except it keeps the legacy percentile scale, so
  # the pair isolates what calibration changes. Expect two differences:
  # unannotated genes stop being vetoed, and deeply-annotated genes lose the
  # advantage they had from being close to everything.
  GS_uncalibrated = list(
    bio_mode         = "semantic",
    components       = c("stability", "utility", "bio"),
    utility_method   = "instance_shap",
    calibration_mode = "percentile",
    alpha_grid       = gs_alpha_grid
  ),

  # --- Pillar combination rule (veto vs compensatory) -------------------
  # Same three pillars, different combination.
  #
  #   geometric (default): a logical AND. A gene at stability=1, utility=1,
  #     biology=0 scores ~0.0005 -- below a gene that is mediocre at all three.
  #     Any pillar can veto. This makes the method VALIDATION-oriented by
  #     construction: a genuinely novel gene with no annotation cannot surface,
  #     because percentile normalisation puts the least-annotated gene at 0.
  #
  #   arithmetic (soft voting): compensatory. The same gene scores 0.667 and
  #     survives. Novel biology is discoverable. The cost is that biology can
  #     no longer veto -- but also that biology alone can CARRY a gene: zero
  #     stability, zero utility, perfect annotation still scores 0.333. A
  #     famous gene that predicts nothing can outrank a real signal.
  #
  # Which is right depends on the intended use, and that should be a stated
  # design choice rather than an accident of the formula. Run both.
  GS_softvote = list(
    bio_mode       = "semantic",
    components     = c("stability", "utility", "bio"),
    utility_method = "instance_shap",
    score_formula  = "arithmetic",
    alpha_grid     = gs_alpha_grid
  ),

  # Harmonic sits between the two: it punishes low scores harder than the
  # arithmetic mean but does not annihilate on a single zero the way the
  # geometric mean effectively does.
  GS_harmonic = list(
    bio_mode       = "semantic",
    components     = c("stability", "utility", "bio"),
    utility_method = "instance_shap",
    score_formula  = "harmonic",
    alpha_grid     = gs_alpha_grid
  )

  # --- Backbone comparison (experimental) -- DISABLED at p = 5000 ---------
  # GS_randomized (randomized_lasso) and GS_sparsegroup (sparse_group_lasso)
  # are removed from this roster, not deleted, because they are intractable at
  # this problem size rather than wrong.
  #
  # Measured: the 2026-08-05 run reached GS_sparsegroup and then ran for over
  # 30 hours on that one method without finishing it. Ten other GS variants and
  # the entire competitor set had completed in about 75 minutes. The cost is
  # the backbone itself -- sparse group lasso is far slower than elastic net at
  # p = 5000 -- multiplied by the evidence-ratio calibration null, which is
  # 10 permutations x 20 subsamples = 200 additional fits per configuration on
  # top of the B = 50 real fits.
  #
  # They are labelled experimental in this config and contribute only the
  # backbone comparison, not the headline result, so the run proceeds without
  # them. To restore: re-add the two entries below and drop top_variable_genes
  # to 2000, where they were previously tractable.
  #
  # GS_randomized = list(
  #   bio_mode              = "none",
  #   components            = c("stability", "utility"),
  #   utility_method        = "instance_shap",
  #   regularization_method = "randomized_lasso",
  #   alpha_grid            = gs_alpha_grid
  # ),
  #
  # GS_sparsegroup = list(
  #   bio_mode              = "none",
  #   components            = c("stability", "utility"),
  #   utility_method        = "instance_shap",
  #   regularization_method = "sparse_group_lasso",
  #   alpha_grid            = c(1.0)
  # )
)

# All methods to benchmark.
all_methods <- c(names(gs_configurations),
                 "DGE", "LASSO", "ElasticNet", "mRMR",
                 "Boruta", "RF_importance", "Random")

# Log the full roster and each GS variant's backbone/biology, now that the
# configurations are defined.
cat("--- Methods to benchmark ---\n")
cat(sprintf("  Competitors: %s\n",
            paste(setdiff(all_methods, names(gs_configurations)),
                  collapse = ", ")))
cat("  GeneSelectR variants:\n")
for (gs_name in names(gs_configurations)) {
  cfg <- gs_configurations[[gs_name]]
  cat(sprintf("    %-20s backbone=%-18s bio=%-10s components=[%s]\n",
              gs_name,
              cfg$regularization_method %||% "elastic_net",
              cfg$bio_mode,
              paste(cfg$components, collapse = "+")))
}
cat(sprintf("  Headline GS config: %s\n", "GS_semantic"))
cat("\n")

# Dispatch a method name to its ranking function.
get_ranker_function <- function(method_name) {
  if (method_name %in% names(gs_configurations)) {
    gs_config <- gs_configurations[[method_name]]
    return(function(train_features, train_labels) {
      rank_with_geneselectr(train_features, train_labels, gs_config,
                            gs_method_label = method_name)
    })
  }

  switch(method_name,
         DGE           = rank_by_differential_expression,
         LASSO         = rank_by_lasso,
         ElasticNet    = rank_by_elastic_net,
         mRMR          = rank_by_mrmr,
         Boruta        = rank_by_boruta,
         RF_importance = rank_by_random_forest,
         Random        = rank_at_random,
         stop("Unknown method: ", method_name)
  )
}


# ==============================================================================
#  Nested cross-validation
# ==============================================================================
#
#  Splits are the outer loop. All methods therefore use one preprocessed split
#  and one GeneSelectR cache. The split-specific cache is released before the
#  next split so memory use does not grow with the number of folds.

normalise_count_split <- function(train_indices, test_indices = NULL) {
  train_counts <- raw_count_matrix[, train_indices, drop = FALSE]

  # Abundance eligibility is learned from the training fold. A gene cannot
  # enter the model because it happened to be expressed in held-out samples.
  keep <- rowSums(train_counts >= min_count_per_gene) >=
    min(min_samples_per_gene, length(train_indices))
  if (!any(keep)) {
    stop("No genes passed the training-fold count filter.")
  }
  train_counts <- train_counts[keep, , drop = FALSE]
  test_counts <- if (is.null(test_indices)) NULL else
    raw_count_matrix[keep, test_indices, drop = FALSE]

  # edgeR's default TMM reference is based on upper-quartile CPM. Select that
  # reference from training samples only and retain it for every held-out
  # sample. This prevents the composition of the test fold from changing the
  # normalized training matrix.
  library_sizes <- colSums(train_counts)
  if (any(library_sizes <= 0)) {
    stop("A training sample has a zero count-library size.")
  }
  upper_quartile_cpm <- apply(train_counts, 2, function(counts) {
    stats::quantile(counts / sum(counts) * 1e6, 0.75,
                    names = FALSE, type = 7)
  })
  reference_index <- which.min(
    abs(upper_quartile_cpm - mean(upper_quartile_cpm))
  )

  train_dge <- edgeR::DGEList(counts = train_counts)
  train_dge <- edgeR::calcNormFactors(train_dge,
                                      refColumn = reference_index)
  train_factors <- train_dge$samples$norm.factors
  train_logcpm <- edgeR::cpm(train_dge, log = TRUE, prior.count = 1)

  if (is.null(test_counts)) {
    return(list(train = t(train_logcpm), test = NULL))
  }

  # Normalize each test sample independently against the fixed training
  # reference. The ratio of the two pairwise TMM factors is invariant to
  # edgeR's product-one rescaling, so it can be placed on the training-factor
  # scale without using any other held-out sample.
  test_logcpm <- vapply(seq_len(ncol(test_counts)), function(j) {
    pair_dge <- edgeR::DGEList(counts = cbind(
      reference = train_counts[, reference_index],
      held_out = test_counts[, j]
    ))
    pair_dge <- edgeR::calcNormFactors(pair_dge, refColumn = 1L)
    pair_ratio <- pair_dge$samples$norm.factors[2] /
      pair_dge$samples$norm.factors[1]
    held_out_factor <- train_factors[reference_index] * pair_ratio
    pair_dge$samples$norm.factors <- c(
      train_factors[reference_index], held_out_factor
    )
    edgeR::cpm(pair_dge, log = TRUE, prior.count = 1)[, 2]
  }, numeric(nrow(train_counts)))
  rownames(test_logcpm) <- rownames(train_counts)
  colnames(test_logcpm) <- colnames(test_counts)

  list(
    train = t(train_logcpm),
    test = t(test_logcpm)
  )
}

preprocess_split <- function(train_indices, test_indices = NULL) {
  normalized <- normalise_count_split(train_indices, test_indices)
  train_expression <- normalized$train
  test_expression <- normalized$test

  if (identical(variance_filter_scope, "train")) {
    training_variances <- apply(train_expression, 2, var)
    eligible <- which(is.finite(training_variances) & training_variances > 0)
    if (length(eligible) == 0L) {
      stop("No positive-variance genes remain in this training split.")
    }
    if (length(eligible) > top_variable_genes) {
      eligible <- eligible[
        order(training_variances[eligible], decreasing = TRUE)[
          seq_len(top_variable_genes)
        ]
      ]
    }
    train_expression <- train_expression[, eligible, drop = FALSE]
    if (!is.null(test_expression)) {
      test_expression <- test_expression[, eligible, drop = FALSE]
    }
  }

  column_means <- colMeans(train_expression)
  column_sds   <- apply(train_expression, 2, sd)
  column_sds[column_sds == 0 | is.na(column_sds)] <- 1

  train_expression <- sweep(sweep(train_expression, 2, column_means, "-"),
                            2, column_sds, "/")
  if (!is.null(test_expression)) {
    test_expression <- sweep(sweep(test_expression, 2, column_means, "-"),
                             2, column_sds, "/")
  }

  list(train = train_expression, test = test_expression)
}

build_cv_jobs <- function(outcome_factor) {
  jobs <- list()
  for (repeat_idx in seq_len(n_outer_repeats)) {
    folds <- make_stratified_folds(
      outcome_factor, k_outer_folds,
      seed = random_seed + 1000 * repeat_idx
    )
    for (fold_idx in seq_len(k_outer_folds)) {
      jobs[[length(jobs) + 1L]] <- list(
        repeat_idx = repeat_idx,
        fold_idx = fold_idx,
        test_indices = folds[[fold_idx]]
      )
    }
  }
  jobs
}

run_split_all_methods <- function(job) {
  test_indices  <- job$test_indices
  train_indices <- setdiff(seq_along(outcome_factor), test_indices)
  train_labels  <- outcome_factor[train_indices]
  test_labels   <- outcome_factor[test_indices]
  split_data    <- preprocess_split(train_indices, test_indices)
  evaluation_rows <- list()

  for (method_name in all_methods) {
    cat(sprintf("  [r%d f%d] %s\n", job$repeat_idx, job$fold_idx,
                method_name))
    ranking_function <- get_ranker_function(method_name)

    start_time <- proc.time()
    ranking_result <- tryCatch(
      ranking_function(split_data$train, train_labels),
      error = function(e) {
        report_failure(method_name, "nested_cv", e)
        list(ranked = sample(colnames(split_data$train)), gs_object = NULL)
      }
    )
    elapsed_seconds <- (proc.time() - start_time)["elapsed"]

    for (panel_size in panel_sizes_evaluated) {
      genes_in_panel <- head(
        ranking_result$ranked[ranking_result$ranked %in%
                                colnames(split_data$train)],
        panel_size
      )
      if (length(genes_in_panel) < 5L) next

      probability_set <- tryCatch(
        predict_with_ensemble(
          split_data$train[, genes_in_panel, drop = FALSE], train_labels,
          split_data$test[, genes_in_panel, drop = FALSE]
        ),
        error = function(e) {
          report_failure(method_name,
                         sprintf("classifier_k%d", panel_size), e)
          na_vector <- rep(NA_real_, length(test_labels))
          list(glmnet = na_vector, xgboost = na_vector,
               rf = na_vector, ensemble = na_vector)
        }
      )

      for (evaluator_name in names(probability_set)) {
        probabilities <- probability_set[[evaluator_name]]
        auc_score <- if (all(is.na(probabilities))) NA_real_ else
          compute_auc(test_labels, probabilities)
        extra_metrics <- if (all(is.na(probabilities))) {
          list(balanced_accuracy = NA_real_, mcc = NA_real_)
        } else {
          compute_classification_metrics(test_labels, probabilities)
        }
        evaluation_rows[[length(evaluation_rows) + 1L]] <- data.frame(
          Method = method_name,
          Evaluator = evaluator_name,
          Repeat = job$repeat_idx,
          Fold = job$fold_idx,
          k = panel_size,
          AUC = auc_score,
          BalAcc = extra_metrics$balanced_accuracy,
          MCC = extra_metrics$mcc,
          Time = elapsed_seconds,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  clear_run_cache()
  do.call(rbind, evaluation_rows)
}


# ==============================================================================
#  Data loading and preprocessing
# ==============================================================================

# IMvigor210 ships as raw counts plus a clinical table. We resolve the source
# automatically, normalise with edgeR (TMM + log-CPM), and end up with the same
# objects the rest of the script expects: `expression_matrix` (samples x genes,
# log-CPM, gene-symbol columns), `metadata`, and `outcome_factor`.

raw_counts   <- NULL
clinical     <- NULL
data_source  <- "none"
response_col <- NA_character_

# --- Option A: easierData (Bioconductor) ------------------------------------
if (requireNamespace("easierData", quietly = TRUE)) {
  cat("Loading IMvigor210 via easierData...\n")
  ok <- tryCatch({
    suppressPackageStartupMessages({
      library(ExperimentHub); library(easierData); library(SummarizedExperiment)
    })
    dat        <- ExperimentHub(cache = experimenthub_cache_dir)[["EH6677"]]
    raw_counts <- as.matrix(SummarizedExperiment::assay(dat, "counts"))
    clinical   <- as.data.frame(SummarizedExperiment::colData(dat))
    data_source <- "easierData"
    TRUE
  }, error = function(e) {
    cat(sprintf("  easierData failed (%s); checking for a local RDS.\n",
                conditionMessage(e)))
    raw_counts <<- NULL
    FALSE
  })
}

# --- Option B: GitHub RDS ---------------------------------------------------
if (is.null(raw_counts)) {
  if (!file.exists(imvigor_rds_path)) {
    stop(paste0(
      "The primary easierData resource failed and no local IMvigor210 RDS ",
      "is available at ", imvigor_rds_path, ". The former GitHub fallback ",
      "currently returns HTTP 404, so it is not used automatically."
    ))
  }
  imv <- readRDS(imvigor_rds_path)

  count_slot <- intersect(c("rawcounts", "counts"), names(imv))
  if (length(count_slot) == 0) {
    stop("Unexpected IMvigor210 RDS structure: no 'rawcounts' or 'counts' slot. ",
         "Found: ", paste(names(imv), collapse = ", "))
  }
  raw_counts  <- as.matrix(imv[[count_slot[1]]])
  clinical    <- as.data.frame(imv$clinical)
  data_source <- "github_rds"
}

# --- Locate the response column ---------------------------------------------
candidate_response_cols <- c("binaryResponse", "BOR",
                             "Best.Confirmed.Overall.Response")
found_cols <- intersect(candidate_response_cols, names(clinical))
if (length(found_cols) == 0) {
  stop("Cannot find a response column. Looked for: ",
       paste(candidate_response_cols, collapse = ", "),
       "\n  Available: ", paste(names(clinical), collapse = ", "))
}
response_col <- found_cols[1]

cat(sprintf("Source: %s | %d genes x %d samples | response column: '%s'\n",
            data_source, nrow(raw_counts), ncol(raw_counts), response_col))

# --- Restrict to evaluable samples and build the outcome --------------------
response_values <- as.character(clinical[[response_col]])
evaluable <- response_values %in% c(responder_code, nonresponder_code)

if (sum(evaluable) == 0) {
  stop("No samples matched the responder/non-responder codes ('",
       responder_code, "'/'", nonresponder_code, "').\n",
       "  Observed values: ",
       paste(utils::head(unique(response_values), 10), collapse = ", "))
}

raw_counts <- raw_counts[, evaluable, drop = FALSE]
clinical   <- clinical[evaluable, , drop = FALSE]

clinical$outcome <- ifelse(
  as.character(clinical[[response_col]]) == responder_code,
  "Responder", "NonResponder"
)
# Responder is the second level, so it is the positive class for AUC.
outcome_factor <- factor(clinical$outcome,
                         levels = c("NonResponder", "Responder"))

cat(sprintf("Evaluable samples: %d\n", length(outcome_factor)))
print(table(outcome_factor))

# Preserve raw integer counts through the outer split. Count filtering, TMM,
# and log-CPM conversion are fitted inside preprocess_split().
if (any(!is.finite(raw_counts)) || any(raw_counts < 0) ||
    any(abs(raw_counts - round(raw_counts)) > .Machine$double.eps^0.5)) {
  stop("IMvigor210 count assay must contain finite non-negative integers.")
}

# Duplicate symbols represent the same count feature and are summed before
# normalization so every split uses one fixed gene universe.
gene_ids <- rownames(raw_counts)
if (any(duplicated(gene_ids))) {
  cat(sprintf("Collapsing %d duplicate gene symbols by summed counts\n",
              sum(duplicated(gene_ids))))
  raw_counts <- rowsum(raw_counts, group = gene_ids, reorder = FALSE)
}

# expression_matrix supplies the fixed sample and gene names used by the
# downstream bookkeeping. Modeling values always come from preprocess_split().
raw_count_matrix <- raw_counts
expression_matrix <- t(raw_count_matrix)
storage.mode(expression_matrix) <- "numeric"

# `metadata` is used by the residualisation machinery and the CV loop. With
# do_residualisation = FALSE it is carried through untouched, but it must exist
# and be row-aligned with the expression matrix.
metadata <- clinical
rownames(expression_matrix) <- rownames(metadata)

cat(sprintf("Loaded %d samples x %d genes\n",
            nrow(expression_matrix), ncol(expression_matrix)))
cat(sprintf("Count filter fitted per training fold: >= %d counts in >= %d samples\n",
            min_count_per_gene, min_samples_per_gene))

if (identical(variance_filter_scope, "global")) {
  gene_variances <- apply(expression_matrix, 2, var)
  expression_matrix <- expression_matrix[
    , is.finite(gene_variances) & gene_variances > 0, drop = FALSE
  ]
  if (ncol(expression_matrix) > top_variable_genes) {
    top_indices <- order(apply(expression_matrix, 2, var),
                         decreasing = TRUE)[seq_len(top_variable_genes)]
    expression_matrix <- expression_matrix[, top_indices, drop = FALSE]
  }
}
cat(sprintf("Candidate filter: %s-fold, p = %d (input genes = %d)\n",
            variance_filter_scope, top_variable_genes,
            ncol(expression_matrix)))


# ==============================================================================
#  Main run: nested CV across all methods
# ==============================================================================

cv_jobs <- build_cv_jobs(outcome_factor)
cat(sprintf("Nested CV: %d splits x %d methods\n",
            length(cv_jobs), length(all_methods)))
nested_start <- proc.time()
nested_results_df <- bind_rows(lapply(cv_jobs, run_split_all_methods))
cat(sprintf("Nested CV done in %.1f seconds\n",
            (proc.time() - nested_start)["elapsed"]))
write.csv(nested_results_df,
          file.path(output_dir, "data", "nested_results.csv"),
          row.names = FALSE)


# ==============================================================================
#  Shared per-subsample ranking cache
# ==============================================================================
#
#  Both the Nogueira and the RENT analyses need each method's gene rankings
#  for every subsample. Computing them once and slicing avoids re-running
#  every method five times for Nogueira and once more for RENT.
#
#  IMPORTANT: this cache mirrors the preprocessing pipeline used by the
#  nested CV — confounder residualisation and standardisation are fit on
#  the subsample's training half and applied to it. Without this, the
#  cached rankings would be computed on different data than the AUC
#  evaluations, and the stability/AUC comparison would be apples-to-oranges.
#
#  The cache uses the same configured subsampling scheme as GeneSelectR. The
#  current benchmark default is repeated K-fold; Nogueira is descriptive under
#  the overlapping K-fold training sets.
#
#  Each entry stores TWO pieces of information per subsample:
#    * `ranked`:   the full ranked gene list (used for top-k slicing)
#    * `selected`: the genes the method natively considers selected
#                  (non-zero for LASSO/EN, p<0.05 BH for DGE, confirmed
#                  for Boruta, passes_rent for GS variants, NULL for
#                  RF and Random which have no natural sparsity).

cat("\n=== Caching per-subsample rankings (with residualisation) ===\n")

subsample_assignments <- create_subsamples(outcome_factor,
                                           B = shared_n_subsamples,
                                           random_seed = random_seed,
                                           scheme = gs_subsample_scheme,
                                           k_folds = gs_subsample_k_folds)

cached_rankings <- stats::setNames(
  lapply(all_methods, function(x) vector("list", shared_n_subsamples)),
  all_methods
)

# Subsamples form the outer loop so all GeneSelectR variants for one matrix
# reuse the same fitted-model and calibration caches. Clearing those entries
# before the next subsample bounds memory independently of B.
for (subsample_idx in seq_len(shared_n_subsamples)) {
  cat(sprintf("  Subsample %d/%d\n", subsample_idx, shared_n_subsamples))

  train_indices <- subsample_assignments[[subsample_idx]]$train

  # Use the identical training-fold filter and scaling path as nested CV.
  train_expression <- preprocess_split(train_indices)$train
  train_labels <- outcome_factor[train_indices]

  for (method_name in all_methods) {
    ranking_function <- get_ranker_function(method_name)
    # --- Run the ranking method on preprocessed data ----------------------
    ranking_attempt <- tryCatch(
      ranking_function(train_expression, train_labels),
      error = function(e) {
        report_failure(method_name, "stability_cache", e)
        NULL
      }
    )

    if (is.null(ranking_attempt)) {
      cached_rankings[[method_name]][[subsample_idx]] <- list(
        ranked   = character(0),
        selected = character(0)
      )
    } else {
      # Filter ranking to genes that exist in the original matrix and
      # cap at shared_max_k to control memory.
      ranked_genes <- head(
        ranking_attempt$ranked[ranking_attempt$ranked %in%
                                 colnames(expression_matrix)],
        shared_max_k
      )

      # selected may be NULL (for methods without natural sparsity) or
      # a character vector. Keep both possibilities clean.
      selected_genes <- if (is.null(ranking_attempt$selected)) {
        NULL
      } else {
        intersect(ranking_attempt$selected, colnames(expression_matrix))
      }

      cached_rankings[[method_name]][[subsample_idx]] <- list(
        ranked   = ranked_genes,
        selected = selected_genes
      )
    }
  }
  clear_run_cache()
}

# Persist the cache so we can re-run analyses without redoing this step.
saveRDS(cached_rankings,
        file.path(output_dir, "data", "cached_rankings.rds"))


# ==============================================================================
#  Analysis 1: Nogueira stability
# ==============================================================================
#
#  We report stability in two complementary ways:
#
#  (A) Native selected sets — the gene set each method actually claims as
#      selected (non-zero coefficient for LASSO/EN, BH-significant for DGE,
#      confirmed for Boruta, passes_rent for GS variants, mRMR's chosen set).
#      This is the apples-to-apples comparison: each method's stability is
#      measured on the genes it actually selects, regardless of count.
#      Methods without natural sparsity (RF, Random) are excluded from
#      this comparison since they have no selected/not-selected distinction.
#
#  (B) Top-k by ranking — for completeness, we also report Nogueira at
#      fixed top-k panel sizes. This is the conventional approach but
#      systematically inflates stability for sparse methods because their
#      lower-ranked positions are dominated by zero-coefficient tie-broken
#      orderings. (A) is the headline number; (B) is supplementary context.

cat("\n=== Nogueira stability: native selected sets ===\n")

nogueira_native_rows <- list()

for (method_name in all_methods) {

  # Check whether this method has a native selected set defined.
  # Methods that returned selected = NULL are skipped for this analysis.
  any_selected_set_defined <- any(
    sapply(cached_rankings[[method_name]],
           function(entry) !is.null(entry$selected))
  )

  if (!any_selected_set_defined) {
    cat(sprintf("  %s: no natural sparsity criterion, skipped\n", method_name))
    next
  }

  # Build a selection matrix where cell (i, b) is TRUE iff gene i is in
  # the method's natively selected set for subsample b.
  selection_matrix <- matrix(FALSE,
                             nrow = ncol(expression_matrix),
                             ncol = shared_n_subsamples)

  selected_set_sizes <- integer(shared_n_subsamples)

  for (subsample_idx in 1:shared_n_subsamples) {
    selected_genes <- cached_rankings[[method_name]][[subsample_idx]]$selected
    if (!is.null(selected_genes) && length(selected_genes) > 0) {
      selection_matrix[match(selected_genes,
                             colnames(expression_matrix)),
                       subsample_idx] <- TRUE
    }
    selected_set_sizes[subsample_idx] <- length(selected_genes %||% character(0))
  }

  stability_result <- compute_nogueira_stability(
    selection_matrix, colnames(expression_matrix)
  )

  nogueira_native_rows[[method_name]] <- data.frame(
    Method            = method_name,
    Nogueira          = stability_result$nogueira_index,
    CI_lower          = stability_result$nogueira_ci[1],
    CI_upper          = stability_result$nogueira_ci[2],
    mean_set_size     = mean(selected_set_sizes),
    median_set_size   = median(selected_set_sizes),
    min_set_size      = min(selected_set_sizes),
    max_set_size      = max(selected_set_sizes),
    stringsAsFactors  = FALSE
  )

  cat(sprintf("  %s: Nogueira = %.3f [%.3f, %.3f] | set size = %.0f (median)\n",
              method_name,
              stability_result$nogueira_index,
              stability_result$nogueira_ci[1],
              stability_result$nogueira_ci[2],
              median(selected_set_sizes)))
}

nogueira_native_df <- bind_rows(nogueira_native_rows)
write.csv(nogueira_native_df,
          file.path(output_dir, "data", "nogueira_native_selected.csv"),
          row.names = FALSE)


cat("\n=== Nogueira stability curves (top-k by ranking, supplementary) ===\n")

nogueira_topk_rows <- list()

for (method_name in all_methods) {
  for (panel_size in nogueira_panel_sizes) {

    # Build top-k selection matrix by slicing each subsample's ranking.
    selection_matrix <- matrix(FALSE,
                               nrow = ncol(expression_matrix),
                               ncol = shared_n_subsamples)

    for (subsample_idx in 1:shared_n_subsamples) {
      top_k_genes <- head(cached_rankings[[method_name]][[subsample_idx]]$ranked,
                          panel_size)
      if (length(top_k_genes) > 0) {
        selection_matrix[match(top_k_genes,
                               colnames(expression_matrix)),
                         subsample_idx] <- TRUE
      }
    }

    stability_result <- compute_nogueira_stability(
      selection_matrix, colnames(expression_matrix)
    )

    nogueira_topk_rows[[length(nogueira_topk_rows) + 1]] <- data.frame(
      Method    = method_name,
      k         = panel_size,
      Nogueira  = stability_result$nogueira_index,
      CI_lower  = stability_result$nogueira_ci[1],
      CI_upper  = stability_result$nogueira_ci[2],
      stringsAsFactors = FALSE
    )
  }
}

# Keep the original variable name `nogueira_df` for downstream figures
# and tradeoff plot — they expect this name.
nogueira_df <- bind_rows(nogueira_topk_rows)
write.csv(nogueira_df,
          file.path(output_dir, "data", "nogueira_curves.csv"),
          row.names = FALSE)


# ==============================================================================
#  Analysis 2: Native-set stability — Jaccard (supplements Nogueira)
# ==============================================================================
#
#  Nogueira on native selected sets is the headline stability metric (it has
#  an analytical variance, giving CIs). As a robustness check we also report
#  mean pairwise Jaccard similarity of the native selected sets across
#  subsamples. Both are computed on each method's actually-selected genes, so
#  neither is inflated by top-k padding of sparse methods.

cat("\n=== Native-set stability: Jaccard (robustness check) ===\n")

mean_pairwise_jaccard <- function(selection_matrix) {
  # selection_matrix: genes x subsamples logical. Mean Jaccard over all
  # subsample pairs. Returns NA if fewer than two non-empty subsamples.
  n_sub <- ncol(selection_matrix)
  sets  <- lapply(seq_len(n_sub), function(b) which(selection_matrix[, b]))
  non_empty <- which(lengths(sets) > 0)
  if (length(non_empty) < 2) return(NA_real_)

  pair_vals <- c()
  for (i in seq_along(non_empty)[-length(non_empty)]) {
    for (j in (i + 1):length(non_empty)) {
      a <- sets[[non_empty[i]]]; b <- sets[[non_empty[j]]]
      union_len <- length(union(a, b))
      jac <- if (union_len == 0) NA_real_ else length(intersect(a, b)) / union_len
      pair_vals <- c(pair_vals, jac)
    }
  }
  mean(pair_vals, na.rm = TRUE)
}

jaccard_rows <- list()
for (method_name in all_methods) {

  any_selected_set_defined <- any(
    vapply(cached_rankings[[method_name]],
           function(entry) !is.null(entry$selected), logical(1))
  )
  if (!any_selected_set_defined) next

  selection_matrix <- matrix(FALSE,
                             nrow = ncol(expression_matrix),
                             ncol = shared_n_subsamples)
  for (subsample_idx in seq_len(shared_n_subsamples)) {
    selected_genes <- cached_rankings[[method_name]][[subsample_idx]]$selected
    if (!is.null(selected_genes) && length(selected_genes) > 0) {
      selection_matrix[match(selected_genes, colnames(expression_matrix)),
                       subsample_idx] <- TRUE
    }
  }

  jaccard_rows[[method_name]] <- data.frame(
    Method = method_name,
    mean_jaccard = mean_pairwise_jaccard(selection_matrix),
    stringsAsFactors = FALSE
  )
}
jaccard_df <- bind_rows(jaccard_rows)
write.csv(jaccard_df,
          file.path(output_dir, "data", "jaccard_native_selected.csv"),
          row.names = FALSE)
if (nrow(jaccard_df) > 0) print(jaccard_df)


# ==============================================================================
#  Analysis 2b: Biology-method ablation (ranked-output comparison)
# ==============================================================================
#
#  The GS variants share the selection backbone, so they differ only in the
#  biology pillar and therefore in their RANKED output — not the selected set.
#  This ablation compares the biology methods on the metrics that actually
#  reflect the ranking: AUC at low k, and STRING coherence of the top-50
#  (computed later). Here we assemble the AUC side; STRING coherence is joined
#  in after the STRING analysis below.
#
#  The comparison is: GS_multilayer vs GS_network vs GS_semantic vs
#  GS_stab_util (no biology) vs GS_stab_only (stability alone).

cat("\n=== Biology-method ablation (AUC at low k) ===\n")

bio_ablation_methods <- intersect(
  c("GS_multilayer", "GS_network", "GS_semantic",
    "GS_stab_util", "GS_stab_only"),
  unique(nested_results_df$Method)
)

bio_ablation_auc <- nested_results_df %>%
  filter(Evaluator == "ensemble", Method %in% bio_ablation_methods,
         k %in% c(10, 20, 50)) %>%
  group_by(Method, k) %>%
  summarise(AUC_mean = mean(AUC, na.rm = TRUE),
            AUC_sd   = sd(AUC,  na.rm = TRUE),
            BalAcc_mean = mean(BalAcc, na.rm = TRUE),
            MCC_mean    = mean(MCC,    na.rm = TRUE),
            .groups = "drop")
write.csv(bio_ablation_auc,
          file.path(output_dir, "data", "biology_ablation_auc.csv"),
          row.names = FALSE)
if (nrow(bio_ablation_auc) > 0) print(bio_ablation_auc)



# ==============================================================================
#  Analysis 3: STRING protein-protein interaction enrichment
# ==============================================================================
#
#  Run each method once on the full data to get its consensus top-50 genes,
#  then query STRING to see how connected those genes are in the PPI network.
#  Biologically coherent gene sets form denser subnetworks than random sets.

cat("\n=== STRING PPI enrichment ===\n")

compute_string_coherence <- function(gene_symbols, string_db_handle,
                                     background_graph, symbol_to_id,
                                     n_permutations = 1000) {

  empty_result <- list(n_mapped = 0L, n_edges = 0L, expected_edges = 0,
                       enrichment_ratio = 0, ppi_p_value = 1,
                       edge_density = 0)

  # --- Map the panel's symbols onto the prebuilt background graph -----------
  panel_ids <- unique(stats::na.omit(symbol_to_id[gene_symbols]))
  panel_ids <- panel_ids[panel_ids %in% igraph::V(background_graph)$name]
  n_nodes   <- length(panel_ids)
  if (n_nodes < 5) {
    empty_result$n_mapped <- n_nodes
    return(empty_result)
  }

  # --- Observed edges among the panel --------------------------------------
  observed_edges <- igraph::ecount(
    igraph::induced_subgraph(background_graph, vids = panel_ids)
  )

  # --- Null: random panels drawn from the CANDIDATE POOL, not the proteome --
  # Sampling the null from all of STRING is too permissive: the candidate pool
  # is the top-N most variable genes in disease-relevant tissue, which are
  # enriched for well-studied, highly-connected genes. Any real panel beats a
  # proteome-wide draw, so every method saturates at the same p-value floor.
  # Drawing the null from the same pool the method selected from asks the right
  # question: is this panel more connected than a random panel of comparably
  # studied genes?
  background_ids <- igraph::V(background_graph)$name
  permuted_edge_counts <- vapply(seq_len(n_permutations), function(i) {
    random_ids <- sample(background_ids, n_nodes)
    igraph::ecount(igraph::induced_subgraph(background_graph,
                                            vids = random_ids))
  }, numeric(1))

  expected_edges <- mean(permuted_edge_counts)

  # Effect size: how many times more connected than a comparable random panel.
  # This is the discriminating quantity -- unlike the p-value, it does not
  # saturate when a method is far above chance.
  enrichment_ratio <- if (expected_edges > 0) {
    observed_edges / expected_edges
  } else {
    0
  }

  # Empirical p-value, add-one smoothed. With n_permutations = 1000 the floor
  # is ~0.001 rather than the ~0.02 that 50 permutations imposed.
  ppi_p_value <- (sum(permuted_edge_counts >= observed_edges) + 1) /
    (n_permutations + 1)

  max_possible_edges <- n_nodes * (n_nodes - 1) / 2
  edge_density <- if (max_possible_edges > 0) {
    observed_edges / max_possible_edges
  } else {
    0
  }

  list(
    n_mapped         = n_nodes,
    n_edges          = observed_edges,
    expected_edges   = expected_edges,
    enrichment_ratio = enrichment_ratio,
    ppi_p_value      = ppi_p_value,
    edge_density     = edge_density
  )
}

string_db <- tryCatch(
  STRINGdb$new(version = string_version,
               species = string_species_id,
               input_directory = string_cache_dir,
               score_threshold = 400),
  error = function(e) {
    cat(sprintf("  STRING initialisation failed: %s\n", e$message))
    NULL
  }
)

# --- Build the background graph ONCE ----------------------------------------
# Every method's null is drawn from the same candidate pool, and the pool never
# changes, so we map it and pull its interactions a single time. This replaces
# n_methods x n_permutations get_interactions() calls with one call plus fast
# in-memory igraph subsetting.
background_graph <- NULL
symbol_to_id     <- NULL

if (!is.null(string_db)) {
  full_data_expression <- preprocess_split(
    seq_len(nrow(expression_matrix))
  )$train
  candidate_genes <- colnames(full_data_expression)
  if (length(candidate_genes) != min(top_variable_genes,
                                     ncol(expression_matrix))) {
    stop("Full-data candidate pool has ", length(candidate_genes),
         " genes; expected ", min(top_variable_genes,
                                   ncol(expression_matrix)), ".")
  }
  cat("  Mapping candidate pool to STRING (once)...\n")
  background_mapping <- tryCatch(
    string_db$map(data.frame(gene = candidate_genes,
                             stringsAsFactors = FALSE),
                  "gene", removeUnmappedRows = TRUE),
    error = function(e) NULL
  )

  if (!is.null(background_mapping) && nrow(background_mapping) > 50) {
    symbol_to_id   <- stats::setNames(background_mapping$STRING_id,
                                      background_mapping$gene)
    background_ids <- unique(background_mapping$STRING_id)

    background_interactions <- tryCatch(
      string_db$get_interactions(background_ids),
      error = function(e) NULL
    )

    if (!is.null(background_interactions) && nrow(background_interactions) > 0) {
      # Column names vary across STRINGdb versions; resolve defensively.
      int_cols <- colnames(background_interactions)
      from_col <- intersect(c("from", "STRING_id.a", "protein1"), int_cols)[1]
      to_col   <- intersect(c("to", "STRING_id.b", "protein2"), int_cols)[1]

      if (!is.na(from_col) && !is.na(to_col)) {
        edge_df <- data.frame(
          from = background_interactions[[from_col]],
          to   = background_interactions[[to_col]],
          stringsAsFactors = FALSE
        )
        edge_df <- edge_df[stats::complete.cases(edge_df), , drop = FALSE]

        background_graph <- igraph::simplify(
          igraph::graph_from_data_frame(edge_df, directed = FALSE)
        )
        # Isolated candidate genes must still be samplable for the null.
        missing_nodes <- setdiff(background_ids,
                                 igraph::V(background_graph)$name)
        if (length(missing_nodes) > 0) {
          background_graph <- igraph::add_vertices(
            background_graph, length(missing_nodes), name = missing_nodes
          )
        }
        cat(sprintf("  Background graph: %d nodes, %d edges\n",
                    igraph::vcount(background_graph),
                    igraph::ecount(background_graph)))
      } else {
        cat(sprintf("  Unexpected get_interactions columns: [%s]\n",
                    paste(int_cols, collapse = ", ")))
      }
    }
  }
  if (is.null(background_graph)) {
    cat("  Background graph unavailable; STRING analysis will be skipped.\n")
  }
}

string_rows <- list()

if (!is.null(string_db)) {
  for (method_name in all_methods) {
    cat(sprintf("  STRING: %s\n", method_name))

    ranking_function <- get_ranker_function(method_name)
    full_data_result <- tryCatch(
      ranking_function(full_data_expression, outcome_factor),
      error = function(e) NULL
    )

    if (is.null(full_data_result)) {
      string_rows[[method_name]] <- data.frame(
        Method = method_name, n_mapped = 0, n_edges = 0,
        expected_edges = 0, enrichment_ratio = 0, ppi_p_value = 1,
        edge_density = 0, neg_log_p = 0, stringsAsFactors = FALSE
      )
      next
    }

    top_50_genes <- head(
      full_data_result$ranked[full_data_result$ranked %in%
                                candidate_genes],
      50
    )

    coherence_result <- if (is.null(background_graph)) {
      list(n_mapped = 0L, n_edges = 0L, expected_edges = 0,
           enrichment_ratio = 0, ppi_p_value = 1, edge_density = 0)
    } else {
      compute_string_coherence(top_50_genes, string_db,
                               background_graph, symbol_to_id,
                               n_permutations = string_n_permutations)
    }

    # Coerce every field to a guaranteed length-1 scalar before building the
    # data frame. This is belt-and-suspenders against any zero-length field
    # that might slip through from the STRING client.
    as_scalar <- function(x, default) {
      if (is.null(x) || length(x) == 0 || is.na(x[1])) default else x[1]
    }

    n_mapped_value       <- as_scalar(coherence_result$n_mapped, 0)
    n_edges_value        <- as_scalar(coherence_result$n_edges, 0)
    expected_edges_value <- as_scalar(coherence_result$expected_edges, 0)
    ppi_p_value_value    <- as_scalar(coherence_result$ppi_p_value, 1)
    edge_density_value   <- as_scalar(coherence_result$edge_density, 0)

    string_rows[[method_name]] <- data.frame(
      Method           = method_name,
      n_mapped         = n_mapped_value,
      n_edges          = n_edges_value,
      expected_edges   = expected_edges_value,
      enrichment_ratio = as_scalar(coherence_result$enrichment_ratio, 0),
      ppi_p_value      = ppi_p_value_value,
      edge_density     = edge_density_value,
      neg_log_p        = -log10(max(ppi_p_value_value, 1e-300)),
      stringsAsFactors = FALSE
    )
  }
}

string_df <- bind_rows(string_rows)
write.csv(string_df,
          file.path(output_dir, "data", "string_coherence.csv"),
          row.names = FALSE)


# ==============================================================================
#  Figures
# ==============================================================================

# Consistent palette: GS variants in reds/oranges, competitors in other hues.
method_colors <- c(
  GS_multilayer       = "#E41A1C",
  GS_network          = "#FF7F00",
  GS_semantic         = "#A6761D",
  GS_stab_util        = "#FB9A99",
  GS_stab_only        = "#A65628",
  GS_no_stability     = "#33A02C",
  GS_utility_only     = "#B2DF8A",
  GS_softvote         = "#CAB2D6",
  GS_harmonic         = "#FDBF6F",
  GS_uncalibrated     = "#E31A1C",
  GS_randomized       = "#6A3D9A",
  GS_sparsegroup      = "#B15928",
  DGE                 = "#4DAF4A",
  LASSO               = "#984EA3",
  ElasticNet          = "#F781BF",
  mRMR                = "#999999",
  Boruta              = "#66C2A5",
  RF_importance       = "#8DA0CB",
  Random              = "grey50"
)

# ------------------------------------------------------------------------------
#  Headline vs supplementary figures
# ------------------------------------------------------------------------------
#  The configuration sweep (GS_multilayer / GS_network / GS_semantic / the
#  ablations) is a configuration study, not a set of competing methods. Headline
#  figures therefore show ONE GeneSelectR point/line -- the configuration named
#  here -- against the competitors. Choose it on the ablation evidence and state
#  the choice in the paper. Detailed all-configuration figures are still written.
headline_gs_config <- "GS_semantic"
competitor_methods <- setdiff(all_methods, names(gs_configurations))

# Preserve evaluator-specific estimates as a sensitivity analysis. The
# pre-specified ensemble supplies the primary summary.
evaluator_sensitivity_summary <- nested_results_df %>%
  group_by(Method, Evaluator, k) %>%
  summarise(AUC_mean = mean(AUC, na.rm = TRUE),
            AUC_sd = sd(AUC, na.rm = TRUE),
            BalAcc_mean = mean(BalAcc, na.rm = TRUE),
            MCC_mean = mean(MCC, na.rm = TRUE),
            .groups = "drop")
write.csv(evaluator_sensitivity_summary,
          file.path(output_dir, "data", "summary_by_evaluator.csv"),
          row.names = FALSE)

performance_summary <- nested_results_df %>%
  filter(Evaluator == "ensemble") %>%
  group_by(Method, k) %>%
  summarise(
    AUC_mean    = mean(AUC, na.rm = TRUE),
    AUC_sd      = sd(AUC, na.rm = TRUE),
    BalAcc_mean = mean(BalAcc, na.rm = TRUE),
    BalAcc_sd   = sd(BalAcc, na.rm = TRUE),
    MCC_mean    = mean(MCC, na.rm = TRUE),
    MCC_sd      = sd(MCC, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(performance_summary,
          file.path(output_dir, "data", "summary.csv"),
          row.names = FALSE)

# --- Figure: three-panel parsimony curves ----------------------------------
plot_auc <- ggplot(performance_summary,
                   aes(x = k, y = AUC_mean, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd),
                width = 0.08, alpha = 0.3) +
  scale_x_log10(breaks = panel_sizes_evaluated) +
  scale_color_manual(values = method_colors) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  theme_bw(base_size = 11) +
  labs(title = "AUC", x = "Top-k genes", y = "AUC")

plot_balanced_accuracy <- ggplot(performance_summary,
                                 aes(x = k, y = BalAcc_mean, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = BalAcc_mean - BalAcc_sd,
                    ymax = BalAcc_mean + BalAcc_sd),
                width = 0.08, alpha = 0.3) +
  scale_x_log10(breaks = panel_sizes_evaluated) +
  scale_color_manual(values = method_colors) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  theme_bw(base_size = 11) +
  labs(title = "Balanced Accuracy", x = "Top-k genes", y = "Balanced Accuracy")

plot_mcc <- ggplot(performance_summary,
                   aes(x = k, y = MCC_mean, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = MCC_mean - MCC_sd, ymax = MCC_mean + MCC_sd),
                width = 0.08, alpha = 0.3) +
  scale_x_log10(breaks = panel_sizes_evaluated) +
  scale_color_manual(values = method_colors) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw(base_size = 11) +
  labs(title = "MCC", x = "Top-k genes", y = "MCC")

ggsave(file.path(output_dir, "figures", "parsimony_3panel.pdf"),
       plot_auc / plot_balanced_accuracy / plot_mcc +
         plot_layout(guides = "collect") &
         theme(legend.position = "bottom"),
       width = 12, height = 14)

# ==============================================================================
#  Performance heatmaps (method x panel size)
# ==============================================================================
#
#  Thirteen methods x six panel sizes x three metrics is unreadable as
#  overlapping curves. A heatmap shows the whole method x k grid at a glance.
#
#  Two versions are produced:
#
#    1. RAW metric values. Useful for absolute reference, but misleading on its
#       own: in bulk expression the Random baseline is NOT 0.5, and it climbs
#       steadily with k as more of the dominant transcriptome axes (immune
#       infiltration, proliferation) get captured by any gene set. Much of the
#       left-to-right gradient in this heatmap is that background, not method
#       quality (cf. Venet, Dutoit & Delorenzi, PLoS Comput Biol 2011: most
#       random gene signatures are significantly associated with outcome).
#
#    2. DELTA vs Random at matched k. This subtracts the background and is the
#       honest comparison: it asks whether a method beats a random gene set of
#       the SAME size. Zero (white) means "no better than random"; negative
#       (blue) means "worse than a random gene set of equal size".

build_metric_heatmap <- function(plot_df, title_text, subtitle_text,
                                 fill_label, diverging = FALSE,
                                 digits = 3) {

  # Order methods by mean value across k, best at the top of the plot.
  method_order <- plot_df %>%
    group_by(Method) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_value) %>%
    pull(Method)

  base_plot <- ggplot(plot_df,
                      aes(x = factor(k),
                          y = factor(Method, levels = method_order),
                          fill = value)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf(paste0("%+.", digits, "f"), value)),
              size = 2.9) +
    theme_minimal(base_size = 11) +
    labs(title = title_text, subtitle = subtitle_text,
         x = "Top-k genes", y = NULL) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 9))

  if (diverging) {
    # Symmetric limits so the white midpoint sits exactly at zero.
    max_abs <- max(abs(plot_df$value), na.rm = TRUE)
    base_plot +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-max_abs, max_abs),
                           name = fill_label)
  } else {
    base_plot +
      scale_fill_viridis_c(name = fill_label, option = "viridis")
  }
}

# --- Heatmap 1: raw metric values -------------------------------------------
raw_heatmaps <- lapply(
  list(list(col = "AUC_mean",    label = "AUC"),
       list(col = "BalAcc_mean", label = "Balanced Accuracy"),
       list(col = "MCC_mean",    label = "MCC")),
  function(spec) {
    plot_df <- performance_summary %>%
      select(Method, k, value = all_of(spec$col))
    build_metric_heatmap(
      plot_df,
      title_text    = spec$label,
      subtitle_text = "Raw value — note the Random row: the baseline is not 0.5 and rises with k",
      fill_label    = spec$label,
      diverging     = FALSE
    )
  }
)

ggsave(file.path(output_dir, "figures", "heatmap_raw_metrics.pdf"),
       raw_heatmaps[[1]] / raw_heatmaps[[2]] / raw_heatmaps[[3]],
       width = 11, height = 16)


# --- Heatmap 2: delta vs Random at matched k (the honest view) ---------------
if ("Random" %in% performance_summary$Method) {

  random_baseline <- performance_summary %>%
    filter(Method == "Random") %>%
    select(k,
           AUC_random    = AUC_mean,
           BalAcc_random = BalAcc_mean,
           MCC_random    = MCC_mean)

  delta_summary <- performance_summary %>%
    filter(Method != "Random") %>%
    left_join(random_baseline, by = "k") %>%
    mutate(AUC_delta    = AUC_mean    - AUC_random,
           BalAcc_delta = BalAcc_mean - BalAcc_random,
           MCC_delta    = MCC_mean    - MCC_random)

  write.csv(delta_summary %>%
              select(Method, k, AUC_delta, BalAcc_delta, MCC_delta),
            file.path(output_dir, "data", "delta_vs_random.csv"),
            row.names = FALSE)

  delta_heatmaps <- lapply(
    list(list(col = "AUC_delta",    label = "AUC"),
         list(col = "BalAcc_delta", label = "Balanced Accuracy"),
         list(col = "MCC_delta",    label = "MCC")),
    function(spec) {
      plot_df <- delta_summary %>%
        select(Method, k, value = all_of(spec$col))
      build_metric_heatmap(
        plot_df,
        title_text    = sprintf("%s minus Random (matched k)", spec$label),
        subtitle_text = "Blue = worse than a random gene set of the same size; white = no better than random",
        fill_label    = sprintf("delta %s", spec$label),
        diverging     = TRUE
      )
    }
  )

  ggsave(file.path(output_dir, "figures", "heatmap_delta_vs_random.pdf"),
         delta_heatmaps[[1]] / delta_heatmaps[[2]] / delta_heatmaps[[3]],
         width = 11, height = 16)

  cat("\n--- AUC minus Random at matched k ---\n")
  print(delta_summary %>%
          select(Method, k, AUC_delta) %>%
          tidyr::pivot_wider(names_from = k, values_from = AUC_delta,
                             names_prefix = "k=") %>%
          arrange(desc(`k=50`)))
}



# --- Figure: parsimony curves, headline (one GeneSelectR line) --------------
# Same data, but showing only the designated GeneSelectR configuration against
# the competitors. The full configuration sweep stays in parsimony_3panel.pdf.
headline_parsimony <- performance_summary %>%
  filter(Method %in% c(headline_gs_config, competitor_methods)) %>%
  mutate(Method = ifelse(Method == headline_gs_config, "GeneSelectR", Method))

if (nrow(headline_parsimony) > 0) {

  headline_curve_colors <- c(method_colors[competitor_methods],
                             GeneSelectR = unname(method_colors[headline_gs_config]))

  # Bold the GeneSelectR line so it reads as the method, not one of many.
  line_widths <- ifelse(levels(factor(headline_parsimony$Method)) == "GeneSelectR",
                        1.6, 0.8)

  build_headline_panel <- function(y_var, y_sd, y_label, hline) {
    ggplot(headline_parsimony,
           aes(x = k, y = .data[[y_var]], color = Method)) +
      geom_line(aes(linewidth = Method == "GeneSelectR")) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = .data[[y_var]] - .data[[y_sd]],
                        ymax = .data[[y_var]] + .data[[y_sd]]),
                    width = 0.08, alpha = 0.25) +
      scale_x_log10(breaks = panel_sizes_evaluated) +
      scale_color_manual(values = headline_curve_colors) +
      scale_linewidth_manual(values = c(`FALSE` = 0.7, `TRUE` = 1.7),
                             guide = "none") +
      geom_hline(yintercept = hline, linetype = "dotted") +
      theme_bw(base_size = 11) +
      labs(title = y_label, x = "Top-k genes", y = y_label)
  }

  ggsave(
    file.path(output_dir, "figures", "headline_parsimony.pdf"),
    build_headline_panel("AUC_mean", "AUC_sd", "AUC", 0.5) /
      build_headline_panel("BalAcc_mean", "BalAcc_sd", "Balanced Accuracy", 0.5) /
      build_headline_panel("MCC_mean", "MCC_sd", "MCC", 0) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom"),
    width = 12, height = 14
  )
}

# --- Figure: Nogueira on native selected sets (HEADLINE FIGURE) ------------
#
# This is the apples-to-apples comparison. Each bar shows the Nogueira
# stability of a method's actual selected set, measured on the genes the
# method itself claims as selected. Methods with no natural sparsity
# criterion (RF, Random) are excluded.

if (nrow(nogueira_native_df) > 0) {

  # Sort methods by stability descending for the bar chart.
  ordered_methods <- nogueira_native_df %>%
    arrange(desc(Nogueira)) %>%
    pull(Method)

  ggsave(
    file.path(output_dir, "figures", "nogueira_native.pdf"),
    ggplot(nogueira_native_df,
           aes(x = factor(Method, levels = rev(ordered_methods)),
               y = Nogueira, fill = Method)) +
      geom_col(alpha = 0.85) +
      geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                    width = 0.3, color = "black") +
      geom_text(aes(label = sprintf("set size: %.0f",
                                    median_set_size)),
                hjust = -0.1, size = 3.5, color = "grey30") +
      scale_fill_manual(values = method_colors) +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(title = "Selection Stability on Natively-Selected Gene Sets",
           subtitle = "Nogueira index with 95% CI | each method's actual selected set",
           x = NULL,
           y = "Nogueira Stability Index") +
      theme(legend.position = "none"),
    width = 10, height = 6
  )
}


# --- Figure: Nogueira stability curves (supplementary) ---------------------
ggsave(
  file.path(output_dir, "figures", "nogueira_curves.pdf"),
  ggplot(nogueira_df, aes(x = k, y = Nogueira, color = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = Method),
                alpha = 0.1, color = NA) +
    scale_x_log10(breaks = nogueira_panel_sizes) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    theme_bw(base_size = 12) +
    labs(title = "Selection Stability at Fixed Top-k (supplementary)",
         subtitle = "Top-k by ranking — may inflate stability for sparse methods due to tie-broken zero tails",
         x = "Top-k genes",
         y = "Nogueira Stability Index (95% CI)"),
  width = 12, height = 7
)

# --- Figure: Biology-method ablation (AUC at low k) ------------------------
# The GS variants share the backbone and differ only in the biology pillar,
# so this compares the biology methods on the ranked output (AUC at the
# clinically-relevant low-k panel sizes).
if (exists("bio_ablation_auc") && nrow(bio_ablation_auc) > 0) {
  ggsave(
    file.path(output_dir, "figures", "biology_ablation.pdf"),
    ggplot(bio_ablation_auc,
           aes(x = factor(k), y = AUC_mean, fill = Method)) +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd),
                    position = position_dodge(width = 0.9), width = 0.2,
                    alpha = 0.4) +
      scale_fill_manual(values = method_colors) +
      geom_hline(yintercept = 0.5, linetype = "dotted") +
      theme_bw(base_size = 12) +
      labs(title = "Biology-Method Ablation",
           subtitle = "Same backbone, different biology pillar | AUC at low k",
           x = "Top-k genes", y = "AUC") +
      theme(legend.position = "bottom"),
    width = 10, height = 6
  )
}

# --- Figure: STRING coherence ----------------------------------------------
# Plotted as fold-enrichment over a random panel of equal size drawn from the
# same candidate pool. The permutation p-value is NOT used as the axis: every
# method that beats all permutations lands on the same floor (1/(nperm+1)), so
# the p-value cannot rank methods that are all far above chance. The ratio can.
# A ratio of 1 (dashed line) means "no more connected than a random panel of
# comparably studied genes".
if (nrow(string_df) > 0 && any(!is.na(string_df$enrichment_ratio))) {
  ggsave(
    file.path(output_dir, "figures", "string_coherence.pdf"),
    ggplot(string_df %>% filter(!is.na(enrichment_ratio)),
           aes(x = reorder(Method, enrichment_ratio),
               y = enrichment_ratio, fill = Method)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("%d obs / %.1f exp  (p=%.3f)",
                                    n_edges, expected_edges, ppi_p_value)),
                hjust = -0.05, size = 3, color = "grey30") +
      geom_hline(yintercept = 1, linetype = "dashed",
                 color = "red", alpha = 0.6) +
      scale_fill_manual(values = method_colors) +
      coord_flip(ylim = c(0, max(string_df$enrichment_ratio,
                                 na.rm = TRUE) * 1.45)) +
      theme_bw(base_size = 12) +
      labs(title = "Biological Coherence: STRING PPI Enrichment",
           subtitle = paste("Top-50 genes | fold-enrichment vs random panels",
                            "from the same candidate pool | dashed = random"),
           x = NULL,
           y = "Observed / expected PPI edges") +
      theme(legend.position = "none"),
    width = 11, height = 6
  )
}

# ==============================================================================
#  HEADLINE FIGURE: one GeneSelectR point vs competitors
# ==============================================================================
#
#  The configuration sweep (GS_multilayer / GS_network / GS_semantic / the
#  ablations) belongs in the supplementary ablation figures. The headline
#  tradeoff plot should show THE METHOD -- a single GeneSelectR point -- against
#  the competing selectors, or readers cannot tell which point is "the method".
#
#  `headline_gs_config` names the configuration that IS GeneSelectR. Choose it
#  on the ablation evidence, and state the choice in the paper. Everything else
#  is reported as a configuration study, not as competing methods.

headline_data <- performance_summary %>%
  filter(k == 50,
         Method %in% c(headline_gs_config, competitor_methods)) %>%
  select(Method, AUC_mean) %>%
  inner_join(nogueira_native_df %>%
               select(Method, Nogueira, median_set_size),
             by = "Method") %>%
  left_join(string_df %>% select(Method, enrichment_ratio), by = "Method") %>%
  mutate(
    enrichment_ratio = ifelse(is.na(enrichment_ratio), 0, enrichment_ratio),
    # Collapse the chosen config's name to the method name.
    Method    = ifelse(Method == headline_gs_config, "GeneSelectR", Method),
    is_ours   = Method == "GeneSelectR"
  )

if (nrow(headline_data) > 0) {

  headline_colors <- c(method_colors[competitor_methods],
                       GeneSelectR = unname(method_colors[headline_gs_config]))

  ggsave(
    file.path(output_dir, "figures", "headline_tradeoff.pdf"),
    ggplot(headline_data,
           aes(x = Nogueira, y = AUC_mean,
               size = enrichment_ratio, color = Method)) +
      geom_point(aes(shape = is_ours), alpha = 0.85) +
      geom_text(aes(label = sprintf("%s (n=%.0f)", Method, median_set_size),
                    fontface = ifelse(is_ours, "bold", "plain")),
                vjust = -1.3, size = 3.2, show.legend = FALSE) +
      scale_color_manual(values = headline_colors) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
      scale_size_continuous(range = c(3, 12), name = "STRING edges\nvs random panel\n(fold enrichment)") +
      scale_x_continuous(expand = expansion(mult = 0.2)) +
      geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.4) +
      theme_bw(base_size = 12) +
      labs(title = "GeneSelectR vs established feature-selection methods",
           subtitle = paste0("k=50 | stability on native selected sets | ",
                             "point size = STRING PPI fold-enrichment vs a random panel of equal size"),
           x = "Nogueira Stability (native selected set)",
           y = "AUC") +
      theme(legend.position = "right"),
    width = 10, height = 7
  )

  cat(sprintf("\nHeadline figure uses '%s' as GeneSelectR.\n", headline_gs_config))
}


# --- Figure: tradeoff scatter, ALL configurations (supplementary) -----------
# Uses the native-selected-set Nogueira (the honest comparison), not the
# top-k version which inflates stability for sparse methods. The inner_join
# drops methods with no native selected set (RF, Random) — plotting them at
# an inflated top-k value would reintroduce the artifact we removed.
tradeoff_data <- performance_summary %>%
  filter(k == 50) %>%
  select(Method, AUC_mean) %>%
  inner_join(nogueira_native_df %>%
               select(Method, Nogueira, median_set_size),
             by = "Method") %>%
  left_join(string_df %>% select(Method, enrichment_ratio),
            by = "Method")

# Methods without STRING coherence are placed at zero on the size axis.
tradeoff_data$enrichment_ratio[is.na(tradeoff_data$enrichment_ratio)] <- 0

ggsave(
  file.path(output_dir, "figures", "tradeoff_scatter.pdf"),
  ggplot(tradeoff_data,
         aes(x = Nogueira, y = AUC_mean,
             size = enrichment_ratio, color = Method)) +
    geom_point(alpha = 0.8) +
    # Label shows method name and its median selected-set size, since a
    # high Nogueira on a 5-gene set is a weaker claim than on a 100-gene set.
    geom_text(aes(label = sprintf("%s (n=%.0f)", Method, median_set_size)),
              vjust = -1.2, size = 3, show.legend = FALSE) +
    scale_color_manual(values = method_colors) +
    scale_size_continuous(range = c(3, 12),
                          name = "STRING edges\nvs random panel\n(fold enrichment)") +
    scale_x_continuous(expand = expansion(mult = 0.18)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.4) +
    theme_bw(base_size = 12) +
    labs(title = "Stability x Prediction x Biology (k=50)",
         subtitle = "Stability = Nogueira on each method's native selected set | point size = STRING PPI fold-enrichment vs a random panel of equal size",
         x = "Nogueira Stability (native selected set)",
         y = "AUC"),
  width = 11, height = 8
)

# --- Figure: AUC boxplots --------------------------------------------------
ggsave(
  file.path(output_dir, "figures", "boxplots_auc.pdf"),
  ggplot(nested_results_df %>% filter(Evaluator == "ensemble"),
         aes(x = factor(k), y = AUC, fill = Method)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    scale_fill_manual(values = method_colors) +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme_bw() +
    labs(x = "Top-k", y = "AUC"),
  width = 14, height = 7
)


# ==============================================================================
#  Pairwise statistical tests
# ==============================================================================
#
#  Paired Wilcoxon signed-rank tests on the ensemble evaluator. Pairing uses
#  identical repeat/fold splits. Repeated CV splits remain dependent, so the
#  p-values are exploratory and paired effect sizes are retained.

for (metric_name in c("AUC", "BalAcc", "MCC")) {
  for (panel_size in c(50, 200)) {

    metric_subset <- nested_results_df %>%
      filter(Evaluator == "ensemble", k == panel_size,
             !is.na(.data[[metric_name]])) %>%
      select(Method, Repeat, Fold, value = all_of(metric_name))

    method_pairs <- combn(sort(unique(metric_subset$Method)), 2,
                          simplify = FALSE)
    pairwise_long <- bind_rows(lapply(method_pairs, function(method_pair) {
      paired <- inner_join(
        filter(metric_subset, Method == method_pair[1]) %>%
          select(Repeat, Fold, value_1 = value),
        filter(metric_subset, Method == method_pair[2]) %>%
          select(Repeat, Fold, value_2 = value),
        by = c("Repeat", "Fold")
      )
      p_value <- if (nrow(paired) < 2L) NA_real_ else suppressWarnings(
        stats::wilcox.test(paired$value_1, paired$value_2,
                           paired = TRUE, exact = FALSE)$p.value
      )
      data.frame(Method1 = method_pair[1], Method2 = method_pair[2],
                 mean_paired_difference = mean(paired$value_1 - paired$value_2),
                 n_pairs = nrow(paired), p_value = p_value)
    })) %>%
      mutate(p_adjusted_bh = p.adjust(p_value, method = "BH"),
             k = panel_size, metric = metric_name)

    write.csv(pairwise_long,
              file.path(output_dir, "data",
                        sprintf("wilcoxon_%s_k%d.csv",
                                metric_name, panel_size)),
              row.names = FALSE)
  }
}


# ==============================================================================
#  Runtime and configuration record
# ==============================================================================

runtime_summary <- nested_results_df %>%
  group_by(Method) %>%
  summarise(mean_time_seconds  = mean(Time, na.rm = TRUE),
            total_time_seconds = sum(Time, na.rm = TRUE),
            .groups = "drop")

write.csv(runtime_summary,
          file.path(output_dir, "data", "runtime.csv"),
          row.names = FALSE)

# Record every parameter so the run is fully reproducible.
saveRDS(list(
  run_date           = run_date,
  run_stamp          = run_stamp,
  geneselectr_version = as.character(utils::packageVersion("GeneSelectR")),
  geneselectr_library = find.package("GeneSelectR"),
  data_source        = data_source,
  experimenthub_resource = if (identical(data_source, "easierData"))
    "EH6677" else NA_character_,
  random_seed        = random_seed,
  n_samples          = nrow(expression_matrix),
  n_genes            = ncol(expression_matrix),
  methods            = all_methods,
  target_go_terms    = target_go_terms,
  disease_term       = disease_term,
  gs_n_subsamples    = gs_n_subsamples,
  gs_gate_method     = gs_gate_method,
  gs_calibration_mode = gs_calibration_mode,
  gs_pfer_bound      = gs_pfer_bound,
  top_variable_genes = top_variable_genes,
  variance_filter_scope = variance_filter_scope,
  count_normalization_scope = "training_fold_fixed_reference_tmm",
  min_count_per_gene = min_count_per_gene,
  min_samples_per_gene = min_samples_per_gene,
  n_parallel_cores  = n_parallel_cores,
  shared_n_subsamples = shared_n_subsamples,
  k_outer_folds      = k_outer_folds,
  n_outer_repeats    = n_outer_repeats
), file.path(output_dir, "data", "config.rds"))


# ==============================================================================
#  Final console summary
# ==============================================================================

cat("\n\n========== RESULTS ==========\n")

# --- Method failures: surfaced loudly ---------------------------------------
# A method that errored into its random fallback looks identical to a method
# that genuinely performs like Random. This block makes that distinction
# unmissable. Any method listed here has results that are NOT its own.
failed_methods <- ls(envir = method_failures)
if (length(failed_methods) > 0) {
  cat("\n!!! METHOD FAILURES DETECTED !!!\n")
  cat("The following methods errored and fell back to a random ranking.\n")
  cat("Their reported metrics are MEANINGLESS — do not interpret them.\n\n")
  for (failed_method in failed_methods) {
    messages <- unique(get(failed_method, envir = method_failures))
    cat(sprintf("  %s (%d failures):\n", failed_method,
                length(get(failed_method, envir = method_failures))))
    for (msg in head(messages, 3)) cat(sprintf("      %s\n", msg))
  }
  cat("\n")

  # Persist for inspection.
  failure_df <- bind_rows(lapply(failed_methods, function(m) {
    data.frame(Method = m,
               error  = unique(get(m, envir = method_failures)),
               stringsAsFactors = FALSE)
  }))
  write.csv(failure_df, file.path(output_dir, "data", "method_failures.csv"),
            row.names = FALSE)
} else {
  cat("\nNo method failures detected.\n")
}

cat("\n--- AUC / Balanced Accuracy / MCC at k=50 ---\n")
print(performance_summary %>%
        filter(k == 50) %>%
        arrange(desc(AUC_mean)) %>%
        select(Method, AUC_mean, AUC_sd, BalAcc_mean, MCC_mean))

cat("\n--- Nogueira on native selected sets (HEADLINE) ---\n")
if (nrow(nogueira_native_df) > 0) {
  print(nogueira_native_df %>% arrange(desc(Nogueira)))
} else {
  cat("(no methods reported a native selected set)\n")
}

cat("\n--- Nogueira at k=50 (top-k by ranking, supplementary) ---\n")
print(nogueira_df %>%
        filter(k == 50) %>%
        arrange(desc(Nogueira)))

cat("\n--- Native-set Jaccard (robustness check) ---\n")
if (exists("jaccard_df") && nrow(jaccard_df) > 0) {
  print(jaccard_df %>% arrange(desc(mean_jaccard)))
}

cat("\n--- Biology-method ablation (AUC at k=50) ---\n")
if (exists("bio_ablation_auc") && nrow(bio_ablation_auc) > 0) {
  print(bio_ablation_auc %>% filter(k == 50) %>% arrange(desc(AUC_mean)))
}

cat("\n--- STRING ---\n")
if (nrow(string_df) > 0) {
  print(string_df %>% arrange(desc(enrichment_ratio)))
}

# --- Backbone comparison: did any regulariser attain the PFER guarantee? -----
# This is the key experimental question for the weak-signal data: does
# randomized lasso or sparse group lasso stabilise selection where elastic net
# could not. We surface the PFER-attainment and native set size per GS variant.
cat("\n--- Backbone comparison (PFER attainment & selection stability) ---\n")
backbone_variants <- intersect(
  c("GS_stab_util", "GS_randomized", "GS_sparsegroup"),
  nogueira_native_df$Method
)
if (length(backbone_variants) > 0) {
  backbone_summary <- nogueira_native_df %>%
    filter(Method %in% backbone_variants) %>%
    select(Method, Nogueira, median_set_size, mean_set_size)
  print(backbone_summary)
  cat("\n  (Higher Nogueira on the native selected set = more stable selection.\n")
  cat("   If GS_randomized or GS_sparsegroup beats GS_stab_util here, the\n")
  cat("   instability was the elastic-net backbone, not the data.)\n")
} else {
  cat("  (backbone variants not present in native Nogueira results)\n")
}

# --- Total runtime -----------------------------------------------------------
benchmark_end_time <- Sys.time()
total_runtime_secs <- as.numeric(difftime(benchmark_end_time,
                                          benchmark_start_time, units = "secs"))
cat(sprintf("\n%s\n", paste(rep("=", 78), collapse = "")))
cat(sprintf("Run complete: %s\n", format(benchmark_end_time, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Total runtime: %.1f seconds (%.2f minutes / %.2f hours)\n",
            total_runtime_secs, total_runtime_secs / 60,
            total_runtime_secs / 3600))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Log file: %s\n", log_file))
cat(paste(rep("=", 78), collapse = ""), "\n")

cat("\nDone. Output directory:", output_dir, "\n")

# Release the log sinks and close the file. The global error handler (set at
# setup) also calls this defensively if the run aborts; here we handle the
# normal completion path and restore the previous error option.
close_log()
options(error = getOption("GeneSelectR_prev_error"))
