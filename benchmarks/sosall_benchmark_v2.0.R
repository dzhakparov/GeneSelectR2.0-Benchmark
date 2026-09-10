# ==============================================================================
#  GeneSelectR 2.0 — SOS-ALL Benchmark (Atopic Dermatitis)
# ==============================================================================
#
#  This script benchmarks GeneSelectR against established feature-selection
#  methods on the SOS-ALL atopic dermatitis cohort. It produces:
#
#    * Nested CV performance (AUC, Balanced Accuracy, MCC) for each method,
#      evaluated with a three-model soft-voting ensemble (elastic net +
#      XGBoost + random forest) to avoid linear-selector / linear-classifier bias
#    * Native-set stability: Nogueira index (headline, with CIs) and mean
#      pairwise Jaccard (robustness), both on each method's actually-selected set
#    * Nogueira top-k curves (supplementary)
#    * STRING protein-protein interaction (PPI) enrichment for biological
#      coherence of the selected gene sets
#    * Biology-method ablation (multilayer vs network vs semantic) on the
#      GS ranked output
#    * A tradeoff scatter plot summarising stability x prediction x biology
#
#  The Nogueira and Jaccard analyses both depend on per-subsample rankings, so
#  the script computes those rankings once (the shared cache) and reuses them.
#
#  Required input files:
#    data/normalized_logcpm.csv  (genes x samples, log2-CPM expression)
#    data/metadata.csv           (sample metadata with treatment column)
#
# ==============================================================================


# ------------------------------------------------------------------------------
#  macOS fork safety -- resolved by NOT FORKING
# ------------------------------------------------------------------------------
#
#  This benchmark used to fork its workers and died doing it:
#
#    objc[NNN]: +[NSCharacterSet initialize] may have been in progress in
#    another thread when fork() was called. ... Crashing instead.
#
#  The trigger is curl, which the biology pillar uses to reach Open Targets,
#  STRING and MSigDB inside every worker.
#
#  OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES DOES NOT FIX IT. That was tested
#  directly -- via ~/.Renviron and exported in the shell before R launches --
#  and forked children doing a curl fetch died 0/4 in both cases. The variable
#  is therefore NOT set here, because setting it would imply a protection that
#  does not exist.
#
#  The workers are PSOCK processes instead. See the backend note in the
#  configuration section for the measurements behind that choice.


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

# Install and load auxiliary packages if missing.
ensure_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name,
                     repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  }
  suppressPackageStartupMessages(library(package_name, character.only = TRUE))
}

ensure_package("mRMRe")        # mRMR feature selection
ensure_package("Boruta")       # Boruta wrapper around random forest
ensure_package("ranger")       # fast random forest for importance ranking
ensure_package("stabm")        # Nogueira stability index
ensure_package("STRINGdb")     # STRING PPI database client
ensure_package("xgboost")      # gradient boosting for the ensemble classifier
ensure_package("knockoff")     # model-X knockoff filter (FDR-controlled gate)


# ------------------------------------------------------------------------------
#  Configuration
# ------------------------------------------------------------------------------
#
#  All tunable parameters are gathered here so the script can be re-run with
#  different settings without searching through the code.

random_seed <- 42
set.seed(random_seed)

# --- File paths -------------------------------------------------------------
expression_file <- "data/normalized_logcpm.csv"
metadata_file   <- "data/metadata.csv"

# --- Metadata schema --------------------------------------------------------
# The treatment column encodes both tissue location and diagnosis in a single
# string like "skin_AD" or "blood_Healthy". We split on the first underscore.
treatment_column <- "treatment"
treatment_split  <- "_"
location_column  <- "location"
diagnosis_column <- "diagnosis"
sample_id_column <- "X"
outcome_column   <- "diagnosis"

# Confounders that will be residualised out before any selection runs.
confounders_categorical <- c(location_column)
do_residualisation      <- TRUE

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
gs_knockoff_fdr       <- 0.1
gs_knockoff_draws     <- 5
# --- knockoff+ vs plain knockoffs: offset = 1 (FDR), TESTED ------------------
#
# offset = 1 is knockoff+ and controls FDR. offset = 0 is plain knockoffs and
# controls a MODIFIED FDR -- a weaker guarantee. This stays at 1.
#
# It was briefly set to 0. The reasoning was that knockoff+ floors the estimated
# FDP at 1/k, so at fdr = 0.1 it cannot certify fewer than 10 genes, while the
# SOS-ALL log showed best thresholds picking a median of 4 -- and that 26.6% of
# individual knockoff draws would clear fdr = 0.10 under offset = 0.
#
# THAT REASONING WAS WRONG, and the mistake is worth keeping written down: the
# 26.6% was a PER-DRAW threshold calculation, but the filter does not certify
# per draw. It derandomises e-values across all knockoff_draws draws first, so a
# quarter of draws individually clearing the threshold does not make the
# aggregate clear it.
#
# Measured directly, full production settings (p = 2000, n = 149, B = 50,
# 5 draws), one fit per cell, ~22 min each:
#
#     offset = 1, fdr = 0.10  ->  0 genes certified
#     offset = 0, fdr = 0.10  ->  0 genes certified
#     offset = 1, fdr = 0.20  ->  0 genes certified
#     (offset = 0, fdr = 0.20 not run; the first three settled it)
#
# The offset makes no difference on this data and neither does doubling the
# target. So there is no case for accepting the weaker mFDR guarantee, and the
# stricter setting stands.
#
# What that means for SOS-ALL: the gate certifies nothing, every GS variant has
# an empty selected set, Nogueira is NA, and the GS variants drop out of the
# native-stability and tradeoff figures. That is consistent with SOS-ALL being a
# dead dataset (everything within 0.08 of Random) and is a result to report, not
# a knob to turn. It has NOT been shown to hold on cohorts with real signal --
# check the "Knockoff filter selected N genes" lines per dataset.
gs_knockoff_offset    <- 1

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
gs_permutation_n      <- 20
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
gs_calibration_permutations <- 5
# Subsamples per permutation when building the calibration null. The null needs
# each gene's EXPECTED selection frequency; B drives the VARIANCE of that
# estimate, not its expectation, and it is averaged over permutations on top.
# Running the null at a smaller B than the real fits is therefore a
# variance-for-speed trade that leaves the evidence ratio unbiased.
gs_calibration_null_B       <- 20
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
#     Rscript sosall_benchmark.R half     ->  results_sosall/<date>_half/
#     Rscript sosall_benchmark.R kfold    ->  results_sosall/<date>_kfold/
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
# One alpha, not two. Every value in this grid doubles the whole benchmark:
# each GS variant is fitted once per alpha and the best is kept by internal CV
# AUC. On weak-signal n << p data the two settings were selecting near-identical
# panels, so the second fit was buying very little. Restore c(0.5, 1.0) if you
# want the sweep back.
gs_alpha_grid         <- c(1.0)
# --- Variance pre-filter ----------------------------------------------------
# Keep the top-N most variable genes before any analysis to bound dimensionality.
#
# Reduced from 5000 to 2000 for runtime. The binding cost is the knockoff gate,
# which builds a p x p covariance matrix, so nearly everything here scales with
# p^2: dropping p by 2.5x cuts that work by ~6x AND shrinks per-worker memory
# from ~1.0 GB to ~0.5 GB, which roughly doubles the number of workers that fit
# in the same RAM. The two effects compound.
#
# Two things make this a fair change rather than a thumb on the scale:
#   * The filter is UNSUPERVISED (variance only, y is never consulted), so it
#     cannot leak outcome information or favour any selector.
#   * It is applied identically to every method, GeneSelectR and competitors
#     alike, before anything else runs.
#
# What it does change is the CANDIDATE POOL, and therefore the null that the
# STRING enrichment is measured against, so the biological-coherence numbers are
# not directly comparable to a p = 5000 run. Re-run SOS-ALL and IMvigor210 at
# this setting before comparing against the earlier results, and state the value
# in the paper -- it is a real analysis choice, not an implementation detail.
#
# 2000 still comfortably covers the signal: published expression signatures are
# tens of genes, and the panels evaluated here top out at k = 500.
.p_override <- Sys.getenv("GS_TOP_VARIABLE_GENES", "")
top_variable_genes <- if (nzchar(.p_override)) as.integer(.p_override) else 2000L
if (length(top_variable_genes) != 1L || !is.finite(top_variable_genes) ||
    top_variable_genes < 1L) {
  stop("GS_TOP_VARIABLE_GENES must be one positive integer.")
}
rm(.p_override)

# Feature eligibility is determined from the training observations in each
# outer split. Although variance filtering is outcome-free, using test-fold
# expression values makes the evaluation transductive and can change which
# genes enter the candidate pool. "global" is retained only for a labelled
# sensitivity analysis against runs produced before this correction.
variance_filter_scope <- Sys.getenv("GS_VARIANCE_FILTER_SCOPE", "train")
if (!variance_filter_scope %in% c("train", "global")) {
  stop("GS_VARIANCE_FILTER_SCOPE must be 'train' or 'global'.")
}

# (Defined above Parallelism: worker memory scales with p^2.)
# --- Parallelism ------------------------------------------------------------
# Applied at the OUTER level: one worker per data split, each running EVERY
# method for that split. Not inside a single GeneSelectR fit -- GeneSelectR
# caches its expensive intermediates keyed on the data matrix, so all methods
# looking at one split share them. Splitting a fold's methods across workers
# would give each its own cache and discard that sharing.
#
# WORKER COUNT IS BOUND BY MEMORY, NOT CORES.
#
# The knockoff gate constructs second-order knockoffs, which requires a p x p
# covariance matrix. At p = 5000 that single object is 5000^2 x 8 bytes = 200 MB,
# and shrinkage estimation holds several copies, so each worker peaks near 1 GB
# during Step 3. Sizing workers by detectCores() therefore exhausts RAM on any
# laptop and the machine swaps to death long before the CPUs saturate.
#
# Workers are capped by an explicit memory budget instead. Raise
# memory_budget_gb only if you know the machine has the headroom to spare.
# How much of the machine's RAM this benchmark may occupy. Detected rather than
# assumed, so the same script sizes itself correctly on a 16 GB laptop and a
# 128 GB server.
#
# Raised from 0.55 to 0.8 deliberately, to get 7 workers instead of 4 on the
# 24 GB machine this is developed on:
#
#     24 GB x 0.8 = 19.2 GB budget / 2.70 GB per worker = 7 workers
#     measured:     7 x 2.34 GB   = 16.4 GB, leaving ~7.6 GB for macOS and the
#                                   parent R session
#
# That is a deliberate move toward the edge, not a safe default. The measured
# per-worker figure is 2.34 GB and the budget assumes 2.70 GB, so the margin is
# the difference between those plus whatever macOS is not already using. If you
# see the machine start swapping -- beachballs, memory pressure in Activity
# Monitor, wall-clock per split climbing -- put this back to 0.55 (4 workers)
# rather than trying to ride it out. A swapping run is slower than a smaller one.
#
# On a machine with different RAM this fraction does NOT give 7 workers; the
# count is printed at startup, so read it rather than assuming.
.memory_fraction_override <- Sys.getenv("GS_MEMORY_FRACTION", "0.55")
memory_fraction <- suppressWarnings(as.numeric(.memory_fraction_override))
if (length(memory_fraction) != 1L || !is.finite(memory_fraction) ||
    memory_fraction <= 0 || memory_fraction > 1) {
  stop("GS_MEMORY_FRACTION must be a number in (0, 1].")
}
rm(.memory_fraction_override)

detect_total_ram_gb <- function() {
  bytes <- tryCatch({
    if (Sys.info()[["sysname"]] == "Darwin") {
      as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
    } else if (file.exists("/proc/meminfo")) {
      kb <- as.numeric(sub("\\D+", "",
                           grep("^MemTotal", readLines("/proc/meminfo"),
                                value = TRUE)[1]))
      kb * 1024
    } else NA_real_
  }, error = function(e) NA_real_, warning = function(w) NA_real_)

  # Fall back to a deliberately pessimistic 8 GB: under-guessing costs a few
  # workers, over-guessing costs the whole run to an out-of-memory crash.
  if (is.na(bytes) || bytes <= 0) 8 else bytes / 1024^3
}

total_ram_gb     <- detect_total_ram_gb()
memory_budget_gb <- total_ram_gb * memory_fraction

# Peak memory per PSOCK worker. MEASURED, not modelled from first principles --
# the previous two formulas here were both wrong, in opposite directions, and the
# second one filled a 24 GB machine and had to be killed.
#
# Measured with `ps -o rss=` on a real worker: fresh R session, the ten packages
# this script loads, the full exported globals, then one knockoff-gated
# geneselectr2_fit on the real SOS-ALL data:
#
#     bare PSOCK worker ........................ 0.07 GB
#     + packages ............................... 0.52 GB
#     + exported globals (55 MB payload) ....... 0.57 GB
#     + ONE knockoff GS fit at p = 2000 ........ 2.40 GB   <- peak
#                            at p = 1000 ........ 1.74 GB
#
# TWO THINGS THAT MODEL GOT WRONG, both worth stating:
#
#   1. It is NOT p^2. Doubling p from 1000 to 2000 raised the fit's own cost
#      from 1.20 GB to 1.83 GB -- a factor of 1.53, where p^2 predicts 4. The
#      p x p knockoff covariance is only 32 MB at p = 2000; it was never the
#      dominant term, so scaling everything by p^2 was wrong in SHAPE, not just
#      in magnitude. The two points fit ~0.57 + 0.00063 * p.
#
#   2. The fixed cost dominates at these p. A PSOCK worker carries its own
#      packages (0.5 GB) where a forked child shared them copy-on-write, so the
#      old fork-era base of 0.45 GB was never right for this backend either.
#
# The constants below are the measured line plus deliberate headroom, because
# the measurement was ONE fit and a real worker runs nineteen methods in
# sequence without R reliably returning memory to the OS between them. Erring
# high costs wall-clock; erring low costs the whole run.
#
# Only two data points, both at n = 149. Re-measure before trusting this at a
# markedly larger n or p -- the script prints its prediction at startup so you
# can compare it against `ps` while a run is warm.
gb_per_worker <- 1.5   # flat: the p^2 knockoff covariance is gone with the gate

.detected_cores  <- parallel::detectCores(logical = FALSE)
.usable_cores    <- if (is.na(.detected_cores)) 1L else max(1L, .detected_cores - 1L)
.memory_limit    <- max(1L, as.integer(floor(memory_budget_gb / gb_per_worker)))

# ------------------------------------------------------------------------------
#  HARD CAP ON WORKERS -- 7, set by the user, not negotiable by any formula
# ------------------------------------------------------------------------------
#  This exists because the memory-derived worker count has silently exceeded the
#  agreed limit twice. The count is computed from RAM / gb_per_worker, so ANY
#  change to gb_per_worker or top_variable_genes changes it -- flattening
#  gb_per_worker to 1.5 and raising p to 5000 produced 11 workers, filled RAM,
#  and crashed the machine.
#
#  Every worker count below is min()'d against this. Raise it only on explicit
#  instruction, and re-measure per-worker RSS with `ps` before you do.
MAX_PARALLEL_WORKERS <- 7L

n_parallel_jobs  <- min(.usable_cores, .memory_limit, MAX_PARALLEL_WORKERS)

# The cap above is per-process, so it does not constrain two runs started in
# separate terminals -- three concurrent runs took 7 each and crashed the
# machine. gs_claim_workers() takes only the share of the 7 that other live
# runs are not already holding, and refuses to start when none is left.
source("benchmarks/worker_budget.R")
n_parallel_jobs  <- gs_claim_workers(n_parallel_jobs, label = "SOS-ALL")

# Inner cores stay at 1: nesting a cluster inside each worker oversubscribes the
# CPU and cluster startup dominates these small jobs.
n_parallel_cores <- 1L

cat(sprintf(paste0("Parallelism: %d worker(s)\n",
                   "  RAM %.0f GB total, %.0f%% budget = %.1f GB",
                   " | ~%.2f GB/worker at p=%d\n",
                   "  limits: %d by cores, %d by memory -> %d used\n"),
            n_parallel_jobs, total_ram_gb, 100 * memory_fraction,
            memory_budget_gb, gb_per_worker, top_variable_genes,
            .usable_cores, .memory_limit, n_parallel_jobs))
rm(.detected_cores, .usable_cores, .memory_limit)

# Stop the BLAS inside each forked worker from spawning its own threads. Without
# this, N workers x M BLAS threads oversubscribe the CPU and thrash.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

# --- Worker backend: PSOCK, NOT fork ------------------------------------------
#
#  THE WORKERS ARE NOT FORKED. This is not a style preference; forking is
#  actively broken for this benchmark on macOS and cost the project a nine-hour
#  run that reported success while computing nothing.
#
#  What happens under fork():
#
#    objc[NNN]: +[NSCharacterSet initialize] may have been in progress in
#    another thread when fork() was called. ... Crashing instead.
#
#  Every worker dies. mclapply raises no condition for a killed child -- it just
#  yields NULL -- so the run continues and writes empty results.
#
#  THE TRIGGER IS curl, not model fitting. Measured, by forking children that do
#  one thing each:
#
#      arithmetic only ................ 4/4 survived
#      DNS resolution ................. 4/4 survived
#      iconv / locale ................. 4/4 survived
#      curl HTTPS fetch ............... 0/4 survived
#
#  The biology pillar queries Open Targets, STRING and MSigDB over HTTPS inside
#  every worker, so every worker touches curl, so every worker dies. This is
#  also why the previous "fork probe" was useless: it fitted glmnet, ranger and
#  xgboost, all of which fork perfectly happily, and reported OK.
#
#  OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES DOES NOT FIX THIS. Verified both via
#  ~/.Renviron and exported in the shell before R launches: 0/4 either way. Any
#  note elsewhere in this project recommending that variable is wrong for this
#  failure and should not be trusted.
#
#  PSOCK spawns fresh R processes instead of forking, so no Obj-C state is
#  inherited and there is nothing for the guard to object to. Measured 4/4
#  survived with curl in the worker, ~1.3s cluster startup. It is also what
#  GeneSelectR itself uses internally (makeCluster(type = "PSOCK")), so the
#  whole stack is already known to work this way.
#
#  What PSOCK costs, stated plainly:
#    * Each worker loads its own packages and receives its own copy of the data
#      (~2.4 MB at p = 2000 -- negligible).
#    * Annotation caches are NOT inherited copy-on-write. Each worker builds its
#      own GO/STRING/Open Targets cache, so first-fetch cost is paid once per
#      worker rather than once per run. If you see API rate limiting, that is
#      the cause; lower n_parallel_jobs.
#
#  What PSOCK does NOT change: the split-outer / method-inner ordering, or
#  GeneSelectR's within-process memoisation. Each worker still runs EVERY method
#  for ONE split in ONE process, so every cache hit that made this fast is
#  intact.
use_psock_cluster <- n_parallel_jobs > 1 &&
  requireNamespace("parallel", quietly = TRUE)

if (use_psock_cluster) {
  cat(sprintf("  Worker backend: PSOCK (%d processes; fork is unusable here)\n",
              n_parallel_jobs))
} else {
  cat("  Worker backend: serial\n")
}

# --- Biology settings -------------------------------------------------------
# GO terms used by supervised biology modes (immune response and related).
target_go_terms <- c(
  "GO:0006955",  # immune response
  "GO:0002376",  # immune system process
  "GO:0006952",  # defense response
  "GO:0002250",  # adaptive immune response
  "GO:0045087",  # innate immune response
  "GO:0002682",  # regulation of immune system process
  "GO:0050776"   # regulation of immune response
)
disease_term <- "atopic dermatitis"

# Glmnet alpha values searched by the downstream classifier (independent of GS).
glmnet_alpha_grid <- c(0.5, 1.0)

# --- Stability / quality analyses -------------------------------------------
# Panel sizes at which the Nogueira index is reported.
nogueira_panel_sizes  <- c(10, 20, 50, 100, 200)
# Subsamples used for the shared Nogueira / Jaccard cache.
# Subsamples for the stability cache. Halved from 50: this feeds a
# supplementary stability estimate, and the Nogueira index has an analytical
# variance, so the cost of 25 is a modestly wider CI rather than a bias.
shared_n_subsamples   <- 25
# Maximum top-k stored per subsample (must cover the largest Nogueira k).
shared_max_k          <- max(nogueira_panel_sizes)

# --- STRING database --------------------------------------------------------
# Permutations for the STRING enrichment null. 1000 gives a p-value floor of
# ~0.001; 50 gave ~0.02, which every real method saturated at.
string_n_permutations <- 1000
string_version    <- "12.0"
string_species_id <- 9606  # Homo sapiens

# --- STRING download timeout and cache --------------------------------------
#
# R's default download timeout is 60 SECONDS, and STRING's human interaction
# file is 79 MB. On the 2026-08-03 run it reached 76 MB and then died:
#
#   downloaded length 79740657 != reported length 83164437
#   URL '.../9606.protein.links.v12.0.txt.gz': Timeout of 60 seconds was reached
#
# get_interactions() therefore returned nothing, background_graph stayed NULL,
# and every method was written to string_coherence.csv with enrichment_ratio = 0
# -- indistinguishable, in the CSV, from a real measurement of zero coherence.
# An 11.5 hour run produced a biology column that was pure artifact.
#
# 3600s is not a considered number, it is simply far more than any plausible
# download needs; the failure mode being guarded against is a hard cutoff
# mid-transfer, not slowness per se.
options(timeout = max(3600, getOption("timeout")))

# Persist the ~100 MB STRING reference files instead of re-fetching them into a
# fresh tempdir() every run. Without input_directory STRINGdb downloads to a
# temporary directory that dies with the session, so every run re-downloads all
# of it and re-rolls the same timeout dice. Downloaded once, this directory
# makes the STRING section start instantly.
string_cache_dir <- file.path("data", "string_db_cache")
dir.create(string_cache_dir, recursive = TRUE, showWarnings = FALSE)

# Point GeneSelectR's own network-biology layer at the SAME cache, so
# bio_mode = "network" reuses these files instead of refetching STRING itself.
options(GeneSelectR.string_cache = normalizePath(string_cache_dir))


# --- Output paths -----------------------------------------------------------
run_date   <- format(Sys.Date(), "%Y-%m-%d")
run_stamp  <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
.output_root <- Sys.getenv("GS_OUTPUT_ROOT", "")
output_dir <- if (nzchar(.output_root)) {
  file.path(.output_root, "sosall",
            sprintf("%s_%s_%sfilter", run_date, gs_subsample_scheme,
                    variance_filter_scope))
} else {
  file.path("results_sosall",
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
cat("GeneSelectR SOS-ALL benchmark run\n")
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
cat(sprintf("  n_parallel_jobs    = %d (outer: splits)\n", n_parallel_jobs))
cat(sprintf("  n_parallel_cores   = %d (inner: within one fit)\n", n_parallel_cores))
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
      #
      # A level present in the test split but NOT in training becomes NA under
      # factor(levels = training_levels), and model.matrix SILENTLY DROPS that
      # row -- which then breaks cbind (two or more confounders) or silently
      # desynchronises the row count from the expression matrix (one
      # confounder). On GSE107994 this killed 13 of 15 splits.
      #
      # Unseen levels are folded into the training reference level: the sample
      # keeps its row and contributes only the intercept, which is the honest
      # handling since no coefficient was estimated for that level.
      unseen <- !(raw_values %in% reference_levels[[confounder]])
      if (any(unseen)) {
        warning(sprintf(paste0("%s: %d test sample(s) have level(s) not seen in ",
                               "training (%s); folded into the reference level '%s'."),
                        confounder, sum(unseen),
                        paste(unique(raw_values[unseen]), collapse = ", "),
                        reference_levels[[confounder]][1]), call. = FALSE)
        raw_values[unseen] <- reference_levels[[confounder]][1]
      }

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


# ------------------------------------------------------------------------------
#  Parallel execution helper
# ------------------------------------------------------------------------------
#
#  IMPORTANT — why we parallelise over DATA SPLITS and not over methods.
#
#  GeneSelectR memoises its expensive intermediates (subsample fits, the
#  knockoff filter, the calibration permutation null) keyed on the data matrix
#  and the backbone. Every one of those keys contains X, so a cache entry can
#  only ever be hit by another method looking at the SAME split. Running all
#  methods for one split in one process therefore turns N scoring variants into
#  roughly one variant's worth of model fitting.
#
#  Parallelising over methods would put each method in its own worker with its
#  own cache and destroy that sharing entirely -- it would look faster per job
#  while doing far more total work. Parallelising over splits keeps every cache
#  hit intact and is embarrassingly parallel on top.
#
#  fork() is used where available so workers inherit the loaded annotation
#  caches copy-on-write rather than rebuilding them. Windows has no fork, so it
#  falls back to serial.

parallel_available <- function() isTRUE(use_psock_cluster)

# Names to hand to a PSOCK worker. Workers get a fresh, empty global
# environment, so every function and configuration value the job touches has to
# be shipped explicitly -- unlike fork, which inherits everything for free.
#
# Exporting the whole global environment is the right call here: the job
# functions reach into dozens of config values (all_methods, gs_configurations,
# top_variable_genes, target_go_terms, expression_matrix, ...) and an
# omission would surface as a confusing "object not found" inside a worker
# hours into a run. The payload is small -- the expression matrix is ~2.4 MB.
#
# Connections are the one thing that must NOT go: the log file connection is
# meaningless in another process, and serialising it is at best useless.
exportable_globals <- function() {
  global_names <- ls(globalenv())
  keep <- vapply(global_names, function(nm) {
    !inherits(get(nm, envir = globalenv()), "connection")
  }, logical(1))
  global_names[keep]
}


# ------------------------------------------------------------------------------
#  Progress reporting
# ------------------------------------------------------------------------------
#
#  A forked worker cannot write to the parent's memory, so an in-process counter
#  is invisible from outside. A file append is the only channel that survives
#  fork(), so every worker appends its own timestamped lines to one progress
#  file and the parent tells you how to watch it.
#
#  This matters because worker stdout is buffered and replayed only after the
#  whole batch finishes: without the file you would see nothing for hours.

progress_file <- file.path(output_dir, "logs", "progress.log")
cat(sprintf("# GeneSelectR progress log -- started %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    file = progress_file)

log_progress <- function(fmt, ...) {
  # Open / write / flush / close on EVERY call, and tag the line with the PID.
  #
  # cat(file = <path>, append = TRUE) looks atomic but is not: each call opens
  # its own handle, and when several forked children write concurrently the
  # buffered writes interleave and lines are lost. In the previous run that
  # silently swallowed methods 2-19 of every split -- the loop was running fine,
  # but only the first method of each worker ever reached the file, which made a
  # perfectly healthy multi-hour run look frozen.
  #
  # Opening in "a" mode and flushing immediately keeps each line whole, and the
  # PID makes it obvious which worker produced it.
  handle <- file(progress_file, open = "a")
  on.exit(close(handle), add = TRUE)
  writeLines(sprintf("[%s][%d] %s", format(Sys.time(), "%H:%M:%S"), Sys.getpid(),
                     sprintf(fmt, ...)), con = handle)
  flush(handle)
  invisible(NULL)
}

# Human-readable duration.
format_duration <- function(seconds) {
  if (is.na(seconds)) return("?")
  if (seconds < 90) return(sprintf("%.0fs", seconds))
  if (seconds < 5400) return(sprintf("%.1fm", seconds / 60))
  sprintf("%.1fh", seconds / 3600)
}


# Runs `fun` over `jobs`, reporting progress.
#
# In PARALLEL mode worker output is captured and replayed in job order once the
# batch completes (so the log does not interleave), and live progress goes to
# progress_file. In SERIAL mode output streams live and a progress bar with an
# ETA is drawn instead.
run_jobs <- function(jobs, fun, label = "job") {

  total      <- length(jobs)
  is_par     <- parallel_available()
  batch_start <- Sys.time()

  log_progress("=== %s: %d jobs, %d worker(s) ===", label, total,
               if (is_par) n_parallel_jobs else 1L)
  cat(sprintf("  %d %ss, %s\n", total, label,
              if (is_par) sprintf("%d parallel workers", n_parallel_jobs)
              else "serial"))
  cat(sprintf("  Live progress:  tail -f %s\n", progress_file))

  run_one <- function(job_index) {
    job_start <- Sys.time()
    log_progress("%s %d/%d START", label, job_index, total)

    execute <- function() {
      tryCatch(fun(jobs[[job_index]]), error = function(e) {
        log_progress("%s %d/%d ERROR: %s", label, job_index, total,
                     conditionMessage(e))
        cat(sprintf("  !! %s %d FAILED: %s\n", label, job_index,
                    conditionMessage(e)))
        NULL
      })
    }

    # Capture only when parallel; serial runs stream live so you can watch the
    # GeneSelectR step output directly.
    if (is_par) {
      job_log <- utils::capture.output({ job_value <- execute() })
    } else {
      job_log <- character(0)
      job_value <- execute()
    }

    elapsed <- as.numeric(difftime(Sys.time(), job_start, units = "secs"))
    log_progress("%s %d/%d DONE in %s", label, job_index, total,
                 format_duration(elapsed))

    list(value = job_value, log = job_log, seconds = elapsed)
  }

  outputs <- if (is_par) {

    # PSOCK, not mclapply -- see the backend note in the configuration section.
    cluster <- parallel::makeCluster(n_parallel_jobs, type = "PSOCK")
    on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)

    # L'Ecuyer streams, so the RNG-dependent methods (Random, subsampling,
    # Boruta) are reproducible across runs. Forked workers were seeded by
    # mclapply's own scheme, which was never reproducible run-to-run.
    parallel::clusterSetRNGStream(cluster, random_seed)

    parallel::clusterExport(cluster, exportable_globals(), envir = globalenv())
    parallel::clusterEvalQ(cluster, suppressPackageStartupMessages({
      library(GeneSelectR); library(glmnet); library(dplyr)
      library(mRMRe); library(Boruta); library(ranger); library(stabm)
      library(STRINGdb); library(xgboost); library(knockoff)
    }))

    # parLapply propagates a worker-level error by aborting the whole call,
    # where mclapply returned a partial result. run_one already catches errors
    # from `fun`, so reaching here means the worker process itself failed --
    # report that as itself rather than as a mysterious batch failure.
    tryCatch(
      parallel::parLapply(cluster, seq_along(jobs), run_one),
      error = function(e) stop(sprintf(paste0(
        "PSOCK cluster failed during the %s batch: %s\n",
        "  This is a worker PROCESS failure, not a method failure.\n",
        "  Common causes: a package missing from the clusterEvalQ list above, ",
        "or memory\n  exhaustion -- each PSOCK worker is a full R session, so ",
        "lower memory_fraction\n  if the machine started swapping."),
        label, conditionMessage(e)), call. = FALSE)
    )

  } else {
    progress_bar <- utils::txtProgressBar(min = 0, max = total, style = 3)
    collected <- vector("list", total)
    for (i in seq_along(jobs)) {
      collected[[i]] <- run_one(i)
      utils::setTxtProgressBar(progress_bar, i)

      # ETA from the mean of completed jobs.
      done_secs <- vapply(collected[seq_len(i)],
                          function(o) o$seconds %||% NA_real_, numeric(1))
      eta <- mean(done_secs, na.rm = TRUE) * (total - i)
      cat(sprintf("   %s %d/%d done | elapsed %s | ETA %s\n", label, i, total,
                  format_duration(as.numeric(difftime(Sys.time(), batch_start,
                                                      units = "secs"))),
                  format_duration(eta)))
    }
    close(progress_bar)
    collected
  }

  # Replay captured worker logs in job order.
  for (out in outputs) {
    if (!is.null(out$log) && length(out$log) > 0) {
      cat(paste(out$log, collapse = "\n"), "\n", sep = "")
    }
  }

  # --- Did the workers actually come back? ---------------------------------
  #
  # THIS IS THE CHECK THAT MATTERS. A forked child killed by the Obj-C fork
  # guard does not raise an R condition in the parent -- mclapply simply yields
  # NULL for that job. Without this block the batch "completes", every
  # downstream table is built from nothing, and the run reports success while
  # having computed zero results.
  #
  # That is not hypothetical. The 2026-07-31 run spent 3.4h on the CV batch and
  # 6.1h on the subsample batch with every worker dead on arrival, wrote a
  # nested_results.csv containing only its header, printed "Nogueira = NA" for
  # all nineteen methods, and exited normally. Nine and a half hours to discover
  # something the first dead worker knew.
  #
  # The tell in that log was "per-job median ?, max -Infs" -- the summary line
  # was already reporting that no job had a recorded duration. Now it stops.
  dead_jobs <- vapply(outputs,
                      function(o) is.null(o) || is.null(o$value), logical(1))
  n_dead <- sum(dead_jobs)

  if (n_dead == total) {
    log_progress("=== %s batch: ALL %d jobs died ===", label, total)
    stop(sprintf(paste0(
      "Every %s (%d/%d) returned nothing -- the workers died rather than ",
      "failed.\n",
      "  Nothing was computed, so the run is being stopped here instead of ",
      "writing empty results.\n\n",
      "  On macOS this is almost always the Objective-C fork guard, and the ",
      "usual cause is\n  running in the RStudio console, where ~/.Renviron is ",
      "read too late to help.\n  Look for lines like:\n",
      "      objc[NNN]: +[NSCharacterSet initialize] ... Crashing instead.\n\n",
      "  Fix: run from a terminal, not the RStudio console --\n",
      "      Rscript benchmarks/sosall_benchmark_v2.0.R %s\n\n",
      "  If it must run in this session, force serial by setting ",
      "memory_fraction low\n  enough that n_parallel_jobs is 1."),
      label, n_dead, total, gs_subsample_scheme))
  }

  if (n_dead > 0) {
    # A partial loss is still a corrupted result -- the surviving splits are a
    # non-random subset (whichever forked at a quiet moment), so the metrics
    # are computed on a biased sample of the CV. Loud, and recorded.
    warning(sprintf(paste0(
      "%d of %d %ss returned nothing (workers died, indices: %s). Results ",
      "below are computed from the survivors only and are NOT the full ",
      "cross-validation."),
      n_dead, total, label, paste(which(dead_jobs), collapse = ", ")),
      call. = FALSE)
    cat(sprintf("\n  !! WARNING: %d of %d %ss died. Results are INCOMPLETE.\n\n",
                n_dead, total, label))
    log_progress("=== %s batch: %d/%d jobs died ===", label, n_dead, total)
  }

  total_secs <- as.numeric(difftime(Sys.time(), batch_start, units = "secs"))
  job_secs <- vapply(outputs, function(o) o$seconds %||% NA_real_, numeric(1))
  cat(sprintf("  %s batch complete: %s wall | %d/%d ok | per-job median %s, max %s\n",
              label, format_duration(total_secs), total - n_dead, total,
              format_duration(stats::median(job_secs, na.rm = TRUE)),
              format_duration(suppressWarnings(max(job_secs, na.rm = TRUE)))))
  log_progress("=== %s batch complete in %s (%d/%d ok) ===", label,
               format_duration(total_secs), total - n_dead, total)

  lapply(outputs, function(out) out$value)
}


# ------------------------------------------------------------------------------
#  Fold preprocessing (shared by every method)
# ------------------------------------------------------------------------------
#
#  Residualisation and standardisation depend only on the split, not on the
#  method, so they are done ONCE per split and reused. Previously this ran once
#  per method per split.

preprocess_split <- function(train_indices, test_indices = NULL) {

  train_expression_raw <- expression_matrix[train_indices, , drop = FALSE]
  train_metadata       <- metadata[train_indices, , drop = FALSE]

  test_expression_raw <- if (is.null(test_indices)) NULL else
    expression_matrix[test_indices, , drop = FALSE]
  test_metadata <- if (is.null(test_indices)) NULL else
    metadata[test_indices, , drop = FALSE]

  # --- Residualise confounders (fit on train, applied to test) -------------
  if (do_residualisation) {
    train_design <- build_design_matrix(train_metadata)
    residualisation_coefficients <- qr.coef(qr(train_design$design_matrix),
                                            train_expression_raw)
    train_expression <- apply_residualisation(
      train_design$design_matrix, train_expression_raw,
      residualisation_coefficients
    )
    if (!is.null(test_expression_raw)) {
      test_design <- build_design_matrix(test_metadata,
                                         train_design$reference_levels)
      test_expression <- apply_residualisation(
        test_design, test_expression_raw, residualisation_coefficients
      )
    } else {
      test_expression <- NULL
    }
  } else {
    train_expression <- train_expression_raw
    test_expression  <- test_expression_raw
  }

  # Fit the unsupervised candidate-pool filter on the training observations.
  # The selected column names are then applied unchanged to the test fold.
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

  # --- Standardise using TRAINING statistics only --------------------------
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
        knockoff_offset       = gs_knockoff_offset,
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
        calibration_null_B    = gs_calibration_null_B,
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
  max_features <- min(ncol(train_features), top_variable_genes)

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
  max_features <- min(ncol(train_features), top_variable_genes)

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
  max_features <- min(ncol(train_features), top_variable_genes)

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
  # Measured twice. On IMvigor210 the 2026-08-05 run reached GS_sparsegroup and
  # spent over 30 hours on that single method without finishing it, while every
  # other GS variant and the whole competitor set completed in about 75
  # minutes. On SOS-ALL the 2026-08-09 run then did the same thing: all 11
  # workers sat on "method 12/19: GS_sparsegroup" for 33 hours and completed
  # ZERO splits.
  #
  # Two costs compound. sparse_group_lasso is far slower than elastic net at
  # p = 5000, and the evidence-ratio calibration adds a permutation null of
  # 10 x 20 = 200 fits per configuration on top of the B = 50 real fits.
  #
  # They contribute only the backbone comparison, not the headline result. To
  # restore them, re-add the entries below AND drop top_variable_genes to 2000,
  # which is where they were last tractable.
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
#  Stratified K-fold CV repeated R times. The loop is organised SPLIT-OUTER,
#  METHOD-INNER: for each (repeat, fold) we preprocess once and then run every
#  method against that same preprocessed split.
#
#  This ordering is what makes GeneSelectR's memoisation effective. Its
#  expensive intermediates -- the B subsample fits, the knockoff filter, the
#  calibration permutation null -- depend only on (X, y, backbone, alpha), not
#  on bio_mode, components or score_formula. All the scoring variants sharing a
#  backbone therefore need the SAME fits and differ only in how they rank them.
#  With method-outer ordering each variant recomputed them from scratch; with
#  split-outer ordering the first variant pays and the rest are re-ranking.
#
#  The numbers are unchanged -- this is caching, not approximation.

build_cv_jobs <- function(outcome_factor) {
  jobs <- list()
  for (repeat_idx in 1:n_outer_repeats) {
    fold_assignments <- make_stratified_folds(
      outcome_factor, k_outer_folds,
      seed = random_seed + 1000 * repeat_idx
    )
    for (fold_idx in 1:k_outer_folds) {
      jobs[[length(jobs) + 1]] <- list(
        repeat_idx   = repeat_idx,
        fold_idx     = fold_idx,
        test_indices = fold_assignments[[fold_idx]]
      )
    }
  }
  jobs
}

# Runs EVERY method against one split. Returns the evaluation rows plus any
# failures, which the parent merges (a forked worker cannot write to the
# parent's failure environment).
run_split_all_methods <- function(job) {

  test_indices  <- job$test_indices
  train_indices <- setdiff(seq_along(outcome_factor), test_indices)

  train_labels <- outcome_factor[train_indices]
  test_labels  <- outcome_factor[test_indices]

  split_data       <- preprocess_split(train_indices, test_indices)
  train_expression <- split_data$train
  test_expression  <- split_data$test

  evaluation_rows <- list()
  local_failures  <- list()

  for (method_name in all_methods) {

    cat(sprintf("  [r%d f%d] %s\n", job$repeat_idx, job$fold_idx, method_name))
    log_progress("  r%d f%d | method %d/%d: %s", job$repeat_idx, job$fold_idx,
                 which(all_methods == method_name), length(all_methods),
                 method_name)

    ranking_function <- get_ranker_function(method_name)

    start_time <- proc.time()
    ranking_result <- tryCatch(
      ranking_function(train_expression, train_labels),
      error = function(e) {
        local_failures[[method_name]] <<- c(
          local_failures[[method_name]],
          sprintf("[nested_cv] %s", conditionMessage(e))
        )
        list(ranked = sample(colnames(train_expression)), gs_object = NULL)
      }
    )
    elapsed_seconds <- (proc.time() - start_time)["elapsed"]

    # --- Evaluate panels at each k ----------------------------------------
    for (panel_size in panel_sizes_evaluated) {

      genes_in_panel <- head(
        ranking_result$ranked[ranking_result$ranked %in%
                                colnames(train_expression)],
        panel_size
      )
      if (length(genes_in_panel) < 5) next

      probability_set <- tryCatch(
        predict_with_ensemble(
          train_expression[, genes_in_panel, drop = FALSE], train_labels,
          test_expression[,  genes_in_panel, drop = FALSE]
        ),
        error = function(e) {
          na_vector <- rep(NA_real_, length(test_labels))
          list(glmnet = na_vector, xgboost = na_vector,
               rf = na_vector, ensemble = na_vector)
        }
      )

      for (evaluator_name in names(probability_set)) {

        predicted_probabilities <- probability_set[[evaluator_name]]

        auc_score <- if (all(is.na(predicted_probabilities))) NA_real_
        else compute_auc(test_labels, predicted_probabilities)

        extra_metrics <- if (all(is.na(predicted_probabilities))) {
          list(balanced_accuracy = NA, mcc = NA)
        } else {
          compute_classification_metrics(test_labels, predicted_probabilities)
        }

        evaluation_rows[[length(evaluation_rows) + 1]] <- data.frame(
          Method    = method_name,
          Evaluator = evaluator_name,
          Repeat    = job$repeat_idx,
          Fold      = job$fold_idx,
          k         = panel_size,
          AUC       = auc_score,
          BalAcc    = extra_metrics$balanced_accuracy,
          MCC       = extra_metrics$mcc,
          Time      = elapsed_seconds,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Release this split's memoised fits. They are keyed on X, so no later split
  # can hit them; holding them would accumulate every split's fits in memory.
  # Annotation caches (GO, PubTator) are deliberately preserved.
  clear_run_cache()

  list(results  = do.call(rbind, evaluation_rows),
       failures = local_failures)
}


# ==============================================================================
#  Data loading and preprocessing
# ==============================================================================

# Read log-CPM expression. Rows are genes (formatted "ENSG__SYMBOL"), columns
# are samples. We collapse to gene symbol by averaging identical symbols.
expression_long <- read.csv(expression_file,
                            row.names = 1, check.names = FALSE,
                            stringsAsFactors = FALSE)

gene_symbols <- sub(".*__", "", rownames(expression_long))
expression_long$gene_symbol <- gene_symbols

expression_collapsed <- aggregate(. ~ gene_symbol,
                                  data = expression_long, FUN = mean)
rownames(expression_collapsed) <- expression_collapsed$gene_symbol
expression_collapsed$gene_symbol <- NULL

# Transpose so rows are samples, columns are genes.
expression_matrix <- t(as.matrix(expression_collapsed))
storage.mode(expression_matrix) <- "numeric"
expression_matrix[is.na(expression_matrix)] <- 0

# Load metadata and parse the treatment column into location + diagnosis.
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
treatment_parts <- strsplit(metadata[[treatment_column]],
                            treatment_split, fixed = TRUE)
metadata[[location_column]]  <- vapply(treatment_parts, `[`,
                                       character(1), 1)
metadata[[diagnosis_column]] <- vapply(treatment_parts, function(parts) {
  if (length(parts) >= 2) {
    paste(parts[-1], collapse = treatment_split)
  } else {
    NA_character_
  }
}, character(1))

# Align expression and metadata by sample ID.
common_samples    <- intersect(rownames(expression_matrix),
                               metadata[[sample_id_column]])
expression_matrix <- expression_matrix[common_samples, , drop = FALSE]
metadata          <- metadata[match(common_samples,
                                    metadata[[sample_id_column]]), ,
                              drop = FALSE]

# Fix the event direction explicitly. Alphabetical factor ordering makes
# Healthy the second level, while every metric treats level 2 as the event.
outcome_factor <- factor(metadata[[outcome_column]],
                         levels = c("Healthy", "AD"))
if (any(is.na(outcome_factor))) {
  stop("SOS-ALL diagnosis contains values outside Healthy and AD.")
}

cat(sprintf("Loaded %d samples x %d genes\n",
            nrow(expression_matrix), ncol(expression_matrix)))
print(table(outcome_factor))

# The historical global filter remains available as an explicitly labelled
# sensitivity setting. Primary runs fit this filter inside preprocess_split().
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

resume_after_stability <- identical(
  Sys.getenv("GS_RESUME_AFTER_STABILITY", "0"), "1"
)
nested_results_path <- file.path(output_dir, "data", "nested_results.csv")
cached_rankings_path <- file.path(output_dir, "data", "cached_rankings.rds")

if (resume_after_stability) {
  if (!file.exists(nested_results_path) || !file.exists(cached_rankings_path)) {
    stop("GS_RESUME_AFTER_STABILITY=1 requires nested_results.csv and ",
         "cached_rankings.rds in the selected output directory.")
  }
  nested_results_df <- read.csv(nested_results_path,
                                stringsAsFactors = FALSE)
  expected_rows <- length(all_methods) * 4L *
    n_outer_repeats * k_outer_folds * length(panel_sizes_evaluated)
  if (nrow(nested_results_df) != expected_rows) {
    stop("Saved nested_results.csv has ", nrow(nested_results_df),
         " rows; expected ", expected_rows, ".")
  }
  cat(sprintf("Resuming from %d saved nested-CV rows.\n",
              nrow(nested_results_df)))
} else {
cv_jobs <- build_cv_jobs(outcome_factor)
cat(sprintf("Nested CV: %d splits x %d methods%s\n",
            length(cv_jobs), length(all_methods),
            if (parallel_available())
              sprintf(" | %d parallel workers", n_parallel_jobs)
            else " | serial"))

nested_start <- proc.time()
split_outputs <- run_jobs(cv_jobs, run_split_all_methods, label = "split")
cat(sprintf("Nested CV done in %.1f seconds\n",
            (proc.time() - nested_start)["elapsed"]))

# Merge failures reported by workers back into the parent's registry.
for (out in split_outputs) {
  if (is.null(out$failures)) next
  for (method_name in names(out$failures)) {
    existing <- if (exists(method_name, envir = method_failures)) {
      get(method_name, envir = method_failures)
    } else character(0)
    assign(method_name, c(existing, out$failures[[method_name]]),
           envir = method_failures)
  }
}

nested_results_df <- bind_rows(lapply(split_outputs, function(o) o$results))
write.csv(nested_results_df,
          nested_results_path,
          row.names = FALSE)
}


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
#  The cache uses repeated K-fold training sets under the current default.
#  Their overlap inflates Nogueira, so stability is descriptive.
#
#  Each entry stores TWO pieces of information per subsample:
#    * `ranked`:   the full ranked gene list (used for top-k slicing)
#    * `selected`: the genes the method natively considers selected
#                  (non-zero for LASSO/EN, p<0.05 BH for DGE, confirmed
#                  for Boruta, passes_rent for GS variants, NULL for
#                  RF and Random which have no natural sparsity).

if (resume_after_stability) {
  cached_rankings <- readRDS(cached_rankings_path)
  if (!setequal(names(cached_rankings), all_methods) ||
      any(lengths(cached_rankings) != shared_n_subsamples)) {
    stop("Saved cached_rankings.rds does not match the current methods or B.")
  }
  cat(sprintf("Resuming from saved rankings for %d methods x %d subsamples.\n",
              length(cached_rankings), shared_n_subsamples))
} else {
cat("\n=== Caching per-subsample rankings ===\n")

subsample_assignments <- create_subsamples(outcome_factor,
                                           B = shared_n_subsamples,
                                           random_seed = random_seed,
                                           scheme = gs_subsample_scheme,
                                           k_folds = gs_subsample_k_folds)

# Same split-outer / method-inner ordering as the nested CV, and for the same
# reason: every method on a given subsample shares one set of cached fits.
# Previously this block ran method-outer, which recomputed the whole GeneSelectR
# pipeline once per variant per subsample and made it the single most expensive
# section of the benchmark -- several times the cost of the nested CV it
# supports.
run_subsample_all_methods <- function(job) {

  subsample_idx <- job$subsample_idx
  train_indices <- job$train_indices

  train_labels     <- outcome_factor[train_indices]
  train_expression <- preprocess_split(train_indices)$train

  per_method <- list()

  for (method_name in all_methods) {

    log_progress("  subsample %d | method %d/%d: %s", subsample_idx,
                 which(all_methods == method_name), length(all_methods),
                 method_name)

    ranking_attempt <- tryCatch(
      get_ranker_function(method_name)(train_expression, train_labels),
      error = function(e) NULL
    )

    if (is.null(ranking_attempt)) {
      per_method[[method_name]] <- list(ranked   = character(0),
                                        selected = character(0))
      next
    }

    ranked_genes <- head(
      ranking_attempt$ranked[ranking_attempt$ranked %in%
                               colnames(expression_matrix)],
      shared_max_k
    )

    # selected may be NULL for methods with no natural sparsity criterion;
    # that is distinct from "ran and selected nothing" and must stay NULL.
    selected_genes <- if (is.null(ranking_attempt$selected)) {
      NULL
    } else {
      intersect(ranking_attempt$selected, colnames(expression_matrix))
    }

    per_method[[method_name]] <- list(ranked   = ranked_genes,
                                      selected = selected_genes)
  }

  clear_run_cache()

  list(subsample_idx = subsample_idx, per_method = per_method)
}

cache_jobs <- lapply(seq_len(shared_n_subsamples), function(i) {
  list(subsample_idx = i,
       train_indices = subsample_assignments[[i]]$train)
})

cat(sprintf("  %d subsamples x %d methods%s\n",
            shared_n_subsamples, length(all_methods),
            if (parallel_available())
              sprintf(" | %d parallel workers", n_parallel_jobs)
            else " | serial"))

cache_start   <- proc.time()
cache_outputs <- run_jobs(cache_jobs, run_subsample_all_methods,
                          label = "subsample")
cat(sprintf("  Caching done in %.1f seconds\n",
            (proc.time() - cache_start)["elapsed"]))

# Re-shape from [subsample][method] to [method][subsample], which is what the
# stability analyses below expect.
cached_rankings <- list()
for (method_name in all_methods) {
  cached_rankings[[method_name]] <- lapply(
    seq_len(shared_n_subsamples),
    function(i) {
      out <- cache_outputs[[i]]
      if (is.null(out) || is.null(out$per_method[[method_name]])) {
        list(ranked = character(0), selected = character(0))
      } else {
        out$per_method[[method_name]]
      }
    }
  )
}


# Persist the cache so we can re-run analyses without redoing this step.
saveRDS(cached_rankings,
        cached_rankings_path)
}


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

cat("\n=== Nogueira stability: ever-selected sets (diagnostic only) ===\n")

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
#  Stability reported at MATCHED TOP-K, not on "native selected sets"
# ==============================================================================
#
#  The native-selected-set comparison has been retired. It never measured what
#  its label claimed. Under the PFER gate the guarantee was unattainable at
#  n << p, so the fallback padded GeneSelectR's set to exactly min_selected
#  genes on every fit -- a constant 20, in both the SOS-ALL and IMvigor210
#  results -- while LASSO, DGE and Boruta contributed genuinely variable native
#  sets. Nogueira is sensitive to set size, so a fixed-size set was being
#  plotted against variable-size ones as though they were the same quantity.
#  Under the knockoff gate that replaced it, the set was empty instead.
#
#  With no gate there is no native set to speak of: pi_raw > 0 means "ever
#  selected", which is hundreds of genes and not comparable across methods with
#  different sparsity.
#
#  So stability is reported where every method is on equal footing: the same
#  top-k. k = 50 matches the panel size used for the tradeoff figures and for
#  the headline AUC comparison.
#
#  nogueira_native_df keeps its name because the figures and joins below
#  consume it, but it now holds top-k values. The CSV written above
#  (nogueira_curves.csv) is the full curve across k.
stability_k <- 50
nogueira_native_df <- nogueira_df %>%
  filter(k == stability_k) %>%
  mutate(mean_set_size   = stability_k,
         median_set_size = stability_k,
         min_set_size    = stability_k,
         max_set_size    = stability_k) %>%
  select(Method, Nogueira, CI_lower, CI_upper,
         mean_set_size, median_set_size, min_set_size, max_set_size)

write.csv(nogueira_native_df,
          file.path(output_dir, "data",
                    sprintf("nogueira_at_k%d.csv", stability_k)),
          row.names = FALSE)
cat(sprintf("\nStability reported at matched top-%d for %d methods\n",
            stability_k, nrow(nogueira_native_df)))


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
               score_threshold = 400,
               # Persistent, so the 100 MB of reference files is fetched once
               # across all runs rather than once per run into a tempdir.
               input_directory = string_cache_dir),
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
  # For the full-data descriptive biology analysis, the entire cohort is the
  # training set. This applies the same residualisation, variance filter, and
  # scaling path used inside each outer-CV training fold.
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
    # Loud, and a real warning() so it resurfaces at the end of the run. The
    # previous behaviour was a single quiet cat() followed by a CSV full of
    # zeros, which is how a failed download became a table of results.
    cat(paste(rep("!", 78), collapse = ""), "\n")
    cat("STRING BACKGROUND GRAPH UNAVAILABLE -- biology column is NOT measured.\n")
    cat("  Every method will be written as NA, not 0. Do not report this column.\n")
    cat("  Usual cause: the 79 MB protein.links download was cut off. Check the\n")
    cat("  lines above for 'Timeout ... was reached' or a length mismatch, then\n")
    cat(sprintf("  delete %s and re-run to re-fetch.\n", string_cache_dir))
    cat(paste(rep("!", 78), collapse = ""), "\n")
    warning("STRING background graph unavailable -- the biology column in ",
            "string_coherence.csv is NA, not measured.", call. = FALSE)
  }
}

string_rows <- list()

if (!is.null(string_db)) {
  # NOTE: this section re-runs EVERY method on the FULL data, serially, to get a
  # consensus top-50 for the PPI test. That is 19 more complete GeneSelectR fits
  # with no parallelism -- it is slow, and previously it logged nothing at all,
  # so a run that had finished both parallel batches looked like it had hung
  # here. It now reports position and elapsed time per method.
  log_progress("=== STRING: %d methods on full data (serial) ===",
               length(all_methods))
  string_section_start <- Sys.time()

  for (method_name in all_methods) {
    method_index <- which(all_methods == method_name)
    cat(sprintf("  STRING: %s (%d/%d)\n", method_name, method_index,
                length(all_methods)))
    log_progress("STRING %d/%d: %s", method_index, length(all_methods),
                 method_name)
    method_start <- Sys.time()

    ranking_function <- get_ranker_function(method_name)
    full_data_result <- tryCatch(
      ranking_function(full_data_expression, outcome_factor),
      error = function(e) NULL
    )

    log_progress("STRING %d/%d: %s done in %s", method_index,
                 length(all_methods), method_name,
                 format_duration(as.numeric(difftime(Sys.time(), method_start,
                                                     units = "secs"))))

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
      # NA, not 0. A zero here reads as "measured, and this panel is no more
      # connected than chance", which is a scientific claim. NA says the
      # measurement did not happen, and the figures already skip NA.
      list(n_mapped = NA_integer_, n_edges = NA_integer_,
           expected_edges = NA_real_, enrichment_ratio = NA_real_,
           ppi_p_value = NA_real_, edge_density = NA_real_)
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
      labs(title = sprintf("Selection Stability at Matched Top-%d", stability_k),
           subtitle = sprintf("Nogueira index with 95%% CI | all methods at k = %d", stability_k),
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
           subtitle = paste0("k=50 | stability = Nogueira at matched top-k | ",
                             "point size = STRING PPI fold-enrichment vs a random panel of equal size"),
           x = "Nogueira Stability (matched top-k)",
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
         subtitle = "Stability = Nogueira at matched top-k | point size = STRING PPI fold-enrichment vs a random panel of equal size",
         x = "Nogueira Stability (matched top-k)",
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
  memory_fraction   = memory_fraction,
  n_parallel_jobs   = n_parallel_jobs,
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

cat("\n--- Nogueira at matched top-50 (HEADLINE) ---\n")
if (nrow(nogueira_native_df) > 0) {
  print(nogueira_native_df %>% arrange(desc(Nogueira)))
} else {
  cat("(no stability values available)\n")
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
  cat("\n  (Higher Nogueira at matched top-k = more stable selection.\n")
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
