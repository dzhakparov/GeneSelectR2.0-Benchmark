# ==============================================================================
#  GeneSelectR 2.0 — STRING rescore (repair a run whose STRING section failed)
# ==============================================================================
#
#  Recomputes ONLY the STRING biological-coherence section of an existing run
#  and rewrites the three outputs that depend on it, leaving everything else in
#  the results directory untouched.
#
#  USAGE
#      Rscript benchmarks/string_rescore.R results_sosall/2026-08-03_kfold
#      Rscript benchmarks/string_rescore.R results_validation/GSE69683/2026-08-05_kfold
#      Rscript benchmarks/string_rescore.R <dir> --dry-run     # print the plan only
#
#  WHY THIS EXISTS
#  ---------------
#  R's default download timeout is 60 seconds and STRING's human interaction
#  file is 79 MB. The 2026-08-03 SOS-ALL run reached 76 MB, timed out, and wrote
#  enrichment_ratio = 0 for all nineteen methods -- a failed download that read
#  as a measurement. The nested CV (5.3h, 15/15 splits) and the subsample cache
#  (4.6h, 25/25) in that run were completely fine.
#
#  Re-running the whole benchmark to repair a 1.5h section would cost ~10h. This
#  repairs the section.
#
#  WHAT IT RECOMPUTES, AND WHAT IT CANNOT REUSE
#  --------------------------------------------
#  The STRING section needs each method's consensus top-50 from a fit on the
#  FULL data. That is not in cached_rankings.rds -- those are per-subsample
#  rankings, a different quantity. So the fits are redone. That is the honest
#  1.5h, and substituting the subsample rankings to save time would silently
#  change what the biology column means.
#
#  HOW IT AVOIDS DRIFTING FROM THE BENCHMARK
#  -----------------------------------------
#  It does not reimplement anything. It evaluates the parent benchmark script
#  (the same one that produced the run) up to just before the nested CV, which
#  yields the identical data, preprocessing, method registry and rankers; then
#  it evaluates that script's OWN STRING section and OWN figure blocks, by
#  matching them in the parsed source.
#
#  The consequence worth stating: if the parent script changes, this follows it
#  automatically -- but a rescore run therefore uses the CURRENT script, not the
#  one that produced the original results. If you have edited the method
#  registry or the preprocessing since, the repaired STRING column is not
#  strictly comparable to the rest of that directory. It is recorded in the
#  rescore manifest so this is visible rather than assumed.
# ==============================================================================


`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

command_args <- commandArgs(trailingOnly = TRUE)
dry_run      <- any(command_args == "--dry-run")
smoke_test   <- any(command_args == "--smoke-test")
positional   <- command_args[!startsWith(command_args, "--")]

if (length(positional) < 1) {
  cat("Usage: Rscript benchmarks/string_rescore.R <results_dir> [--dry-run|--smoke-test]\n\n")
  cat("  <results_dir> is a completed run directory, e.g.\n")
  cat("      results_sosall/2026-08-03_kfold\n")
  cat("      results_validation/GSE69683/2026-08-05_kfold\n\n")
  cat("  --dry-run     print the plan and stop\n")
  cat("  --smoke-test  run the FULL pipeline with three cheap methods into a\n")
  cat("                throwaway directory, to prove it works before you commit\n")
  cat("                ~1.5h to the real thing. Touches nothing.\n")
  quit(status = 1)
}

results_dir <- sub("/+$", "", positional[1])

# The real rescore refits nineteen methods on the full data and takes about an
# hour and a half. Enough long runs in this project have died at the very end
# that it is worth being able to exercise the whole path -- prefix evaluation,
# STRING section, CSV write, all three figures -- in a few minutes first.
#
# Smoke mode keeps every step and shrinks only the method roster, and works on a
# copy so the run being repaired cannot be touched.
smoke_methods <- c("DGE", "LASSO", "Random")


# ------------------------------------------------------------------------------
#  Work out which benchmark produced this directory
# ------------------------------------------------------------------------------
#
#  The layout carries the identity, and it has to, because the two benchmarks
#  record different things in config.rds -- the validation script stores the GEO
#  accession, the SOS-ALL script has no dataset field at all. Parsing the path
#  is what makes this work for BOTH without guessing.
#
#      results_sosall/<date>_<scheme>                -> sosall_benchmark_v2.0.R <scheme>
#      results_validation/<ACCESSION>/<date>_<scheme> -> validation_benchmark.R <ACC> <scheme>
#
#  Anything else stops, rather than being handled by a default that would run
#  the wrong dataset against the right-looking directory.
identify_run <- function(dir) {
  parts <- strsplit(normalizePath(dir, mustWork = FALSE), .Platform$file.sep)[[1]]
  n     <- length(parts)
  leaf  <- parts[n]

  scheme <- sub("^\\d{4}-\\d{2}-\\d{2}_", "", leaf)
  if (!scheme %in% c("half", "kfold")) {
    stop(sprintf(paste0(
      "Cannot read a subsampling scheme from '%s'.\n",
      "  Expected a run directory named <date>_half or <date>_kfold."), leaf))
  }

  if (n >= 2 && parts[n - 1] == "results_sosall") {
    return(list(kind = "sosall", scheme = scheme, accession = NA_character_,
                script = "benchmarks/sosall_benchmark_v2.0.R",
                args = scheme, label = "SOS-ALL (atopic dermatitis)"))
  }
  if (n >= 3 && parts[n - 2] == "results_validation") {
    acc <- parts[n - 1]
    return(list(kind = "validation", scheme = scheme, accession = acc,
                script = "benchmarks/validation_benchmark.R",
                args = c(acc, scheme), label = acc))
  }

  stop(sprintf(paste0(
    "Unrecognised results layout: %s\n",
    "  Expected results_sosall/<date>_<scheme> or ",
    "results_validation/<ACCESSION>/<date>_<scheme>."), dir))
}

run <- identify_run(results_dir)

# Everything the rescore needs from the original run, checked up front so a
# missing input fails now rather than after the fits.
required_inputs <- c(
  config  = file.path(results_dir, "data", "config.rds"),
  summary = file.path(results_dir, "data", "summary.csv"),
  nog     = file.path(results_dir, "data", "nogueira_native_selected.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

cat(paste(rep("=", 78), collapse = ""), "\n")
cat("STRING rescore\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
cat(sprintf("  results dir : %s\n", results_dir))
cat(sprintf("  dataset     : %s\n", run$label))
cat(sprintf("  benchmark   : %s %s\n", run$script, paste(run$args, collapse = " ")))
cat(sprintf("  will rewrite: data/string_coherence.csv\n"))
cat(sprintf("                figures/string_coherence.pdf\n"))
cat(sprintf("                figures/tradeoff_scatter.pdf\n"))
cat(sprintf("                figures/headline_tradeoff.pdf\n"))
cat(sprintf("  untouched   : nested_results.csv, summary.csv, nogueira_*, jaccard_*,\n"))
cat(sprintf("                cached_rankings.rds, and every other figure\n"))

if (length(missing_inputs) > 0) {
  cat("\n  !! MISSING INPUTS:\n")
  for (nm in names(missing_inputs)) cat(sprintf("      %s\n", missing_inputs[[nm]]))
  stop("This does not look like a completed run directory.")
}
if (!file.exists(run$script)) {
  stop(sprintf("Benchmark script not found: %s (run from the project root)",
               run$script))
}

if (dry_run) {
  cat("\n  --dry-run: stopping before any work.\n")
  quit(status = 0)
}

# In smoke mode, redirect everything at a throwaway copy. The inputs the figures
# read (summary.csv, nogueira_native_selected.csv) are copied across so the
# figure code exercises real data rather than placeholders.
if (smoke_test) {
  smoke_dir <- file.path(tempdir(), paste0("string_smoke_", Sys.getpid()))
  dir.create(file.path(smoke_dir, "data"),    recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(smoke_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  file.copy(required_inputs, file.path(smoke_dir, "data"), overwrite = TRUE)
  real_results_dir <- results_dir
  results_dir      <- smoke_dir
  required_inputs  <- file.path(smoke_dir, "data", basename(required_inputs))
  names(required_inputs) <- c("config", "summary", "nog")

  cat("\n")
  cat(paste(rep("~", 78), collapse = ""), "\n")
  cat("SMOKE TEST -- proving the pipeline, not producing results.\n")
  cat(sprintf("  methods reduced to: %s\n", paste(smoke_methods, collapse = ", ")))
  cat(sprintf("  writing to:         %s\n", smoke_dir))
  cat(sprintf("  NOT touching:       %s\n", real_results_dir))
  cat(paste(rep("~", 78), collapse = ""), "\n")
}


# ------------------------------------------------------------------------------
#  Source-level helpers
# ------------------------------------------------------------------------------
#
#  These evaluate selected top-level expressions OUT OF the parent script, in
#  file order, into the global environment. That is what keeps this repair
#  byte-identical to what a real run would have computed.

deparse1 <- function(e) paste(deparse(e), collapse = " ")

# Evaluate from the top of the script until `stop_pattern` matches -- i.e. all
# configuration, helpers, the method registry and the data, but not the CV.
eval_prefix <- function(path, stop_pattern) {
  for (e in parse(path)) {
    if (grepl(stop_pattern, deparse1(e))) return(invisible(TRUE))
    eval(e, envir = globalenv())
  }
  stop("Never reached the nested-CV marker; the benchmark script has changed shape.")
}

# Evaluate every top-level expression whose deparsed source matches any pattern,
# in file order (so assignments land before the blocks that consume them).
eval_matching <- function(path, patterns, label) {
  hits <- 0L
  for (e in parse(path)) {
    t <- deparse1(e)
    if (any(vapply(patterns, function(p) grepl(p, t), logical(1)))) {
      eval(e, envir = globalenv())
      hits <- hits + 1L
    }
  }
  cat(sprintf("  %s: evaluated %d block(s)\n", label, hits))
  hits
}


# ------------------------------------------------------------------------------
#  Rebuild the run's world
# ------------------------------------------------------------------------------

# The parent script reads its own command line; shadow it. R resolves this from
# globalenv before base, so the script sees the arguments we want.
commandArgs <- function(trailingOnly = FALSE) run$args

# Sourcing the parent creates a NEW output directory for today. Note whether it
# already existed, so we only remove one we created ourselves -- deleting a real
# run that happens to share today's date would be unforgivable.
todays_dir <- if (run$kind == "sosall") {
  file.path("results_sosall", sprintf("%s_%s", format(Sys.Date(), "%Y-%m-%d"),
                                      run$scheme))
} else {
  file.path("results_validation", run$accession,
            sprintf("%s_%s", format(Sys.Date(), "%Y-%m-%d"), run$scheme))
}
todays_dir_preexisted <- dir.exists(todays_dir)

cat("\n--- Rebuilding the run's configuration and data ---\n")
eval_prefix(run$script, "^cv_jobs <- build_cv_jobs")

# The parent installed log sinks and an error handler; release both so this
# script's output goes to the console and errors surface normally.
while (sink.number() > 0) sink()
if (sink.number(type = "message") != 2) sink(type = "message")
options(error = NULL)
try(close(log_connection), silent = TRUE)

# Point every output path at the directory being repaired.
output_dir    <<- results_dir
progress_file <<- file.path(results_dir, "logs", "string_rescore_progress.log")
dir.create(file.path(results_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
cat(sprintf("# STRING rescore -- started %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = progress_file)

# Remove the empty directory the parent script created on the way past.
if (!todays_dir_preexisted && dir.exists(todays_dir) &&
    !identical(normalizePath(todays_dir), normalizePath(results_dir))) {
  unlink(todays_dir, recursive = TRUE)
}

if (smoke_test) {
  all_methods <<- intersect(smoke_methods, all_methods)
}

cat(sprintf("  data: %d samples x %d genes | %d methods\n",
            nrow(expression_matrix), ncol(expression_matrix), length(all_methods)))
cat(sprintf("  STRING cache: %s | download timeout: %ds\n",
            string_cache_dir, getOption("timeout")))


# ------------------------------------------------------------------------------
#  Preserve the failed table before overwriting it
# ------------------------------------------------------------------------------

string_csv <- file.path(results_dir, "data", "string_coherence.csv")
if (file.exists(string_csv)) {
  backup <- file.path(results_dir, "data",
                      sprintf("string_coherence.superseded_%s.csv",
                              format(Sys.time(), "%Y%m%d_%H%M%S")))
  file.copy(string_csv, backup, overwrite = FALSE)
  cat(sprintf("\n  previous table kept as %s\n", basename(backup)))
}


# ------------------------------------------------------------------------------
#  Run the parent script's own STRING section
# ------------------------------------------------------------------------------
#
#  Everything from compute_string_coherence() through the write.csv of
#  string_coherence.csv, evaluated verbatim out of the benchmark script.

cat("\n--- STRING section (re-fitting every method on the full data) ---\n")
string_start <- Sys.time()

string_section_hits <- eval_matching(
  run$script,
  c("^compute_string_coherence <- function",
    "^string_db <- tryCatch",
    "^background_graph <- NULL",
    "^symbol_to_id <- NULL",
    "Mapping candidate pool to STRING",
    "^string_rows <- list\\(\\)",
    "STRING: %s \\(%d/%d\\)",
    "string_df <- bind_rows\\(string_rows\\)",
    "\"string_coherence\\.csv\""),
  "STRING section")

if (string_section_hits < 7) {
  stop("Recognised too few STRING blocks in the benchmark script; it has been ",
       "restructured and this rescore's patterns need updating.")
}

cat(sprintf("  STRING section took %.1f minutes\n",
            as.numeric(difftime(Sys.time(), string_start, units = "mins"))))

# Did it actually work this time? The whole point of this script is that a
# failed STRING section must not look like a result.
if (!exists("string_df") || nrow(string_df) == 0) {
  stop("STRING produced no rows at all.")
}
measured <- sum(is.finite(string_df$enrichment_ratio) &
                  string_df$enrichment_ratio > 0)
if (measured == 0) {
  stop(sprintf(paste0(
    "STRING FAILED AGAIN -- 0 of %d methods have a usable enrichment ratio.\n",
    "  The repaired table has NOT been accepted. Check above for a download\n",
    "  timeout or length mismatch, delete %s, and re-run."),
    nrow(string_df), string_cache_dir))
}
cat(sprintf("  %d of %d methods scored\n", measured, nrow(string_df)))


# ------------------------------------------------------------------------------
#  Rebuild the three STRING-dependent figures
# ------------------------------------------------------------------------------
#
#  These read performance_summary and nogueira_native_df, which the original run
#  already computed and wrote -- they are loaded from the results directory
#  rather than recomputed, so the repaired figures use the SAME prediction and
#  stability numbers as the rest of the run.

cat("\n--- Figures ---\n")
performance_summary <- utils::read.csv(required_inputs[["summary"]],
                                       stringsAsFactors = FALSE)
nogueira_native_df  <- utils::read.csv(required_inputs[["nog"]],
                                       stringsAsFactors = FALSE)
cat(sprintf("  loaded summary.csv (%d rows) and nogueira_native_selected.csv (%d rows)\n",
            nrow(performance_summary), nrow(nogueira_native_df)))

figure_hits <- eval_matching(
  run$script,
  c("^method_colors <- c\\(",
    "^headline_gs_config <- ",
    "^competitor_methods <- ",
    "^figure_subtitle <- ",        # validation only; absent in SOS-ALL
    "\"string_coherence\\.pdf\"",
    "^headline_data <- performance_summary",
    "\"headline_tradeoff\\.pdf\"",
    "^tradeoff_data <- performance_summary",
    "tradeoff_data\\$enrichment_ratio\\[is\\.na",
    "\"tradeoff_scatter\\.pdf\""),
  "figure blocks")

if (figure_hits < 8) {
  warning("Fewer figure blocks matched than expected -- check that all three ",
          "PDFs below were actually rewritten.", call. = FALSE)
}


# ------------------------------------------------------------------------------
#  Record what was done
# ------------------------------------------------------------------------------
#
#  A repaired directory is a mixture: most of it from the original run, the
#  STRING column from today, possibly under a since-edited script. That mixture
#  has to be written down or it becomes invisible.

manifest <- list(
  rescored_on        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  results_dir        = results_dir,
  dataset            = run$label,
  accession          = run$accession,
  benchmark_script   = run$script,
  scheme             = run$scheme,
  methods_scored     = measured,
  methods_total      = nrow(string_df),
  top_variable_genes = top_variable_genes,
  string_version     = string_version,
  string_cache_dir   = string_cache_dir,
  n_permutations     = string_n_permutations,
  note = paste("STRING section recomputed after a download timeout voided the",
               "original. Prediction and stability results in this directory",
               "are from the original run and were not touched.")
)
saveRDS(manifest, file.path(results_dir, "data", "string_rescore_manifest.rds"))

cat("\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
cat("Rescore complete\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
print(string_df[order(-string_df$enrichment_ratio),
                c("Method", "n_mapped", "n_edges", "expected_edges",
                  "enrichment_ratio", "ppi_p_value")])
cat(sprintf("\n  wrote %s\n", string_csv))
made <- file.path(results_dir, "figures",
                  c("string_coherence.pdf", "tradeoff_scatter.pdf",
                    "headline_tradeoff.pdf"))
for (p in made) {
  cat(sprintf("  %s %s\n", if (file.exists(p)) "wrote" else "MISSING ->", p))
}
if (!all(file.exists(made))) {
  stop("At least one STRING figure was not written -- see MISSING above.")
}
cat(sprintf("  manifest: %s/data/string_rescore_manifest.rds\n", results_dir))

if (smoke_test) {
  cat("\n")
  cat(paste(rep("~", 78), collapse = ""), "\n")
  cat("SMOKE TEST PASSED -- the full path works end to end.\n")
  cat(sprintf("  These numbers are from %d methods and are NOT results.\n",
              length(all_methods)))
  cat(sprintf("  %s was not touched. Run for real with:\n", real_results_dir))
  cat(sprintf("      Rscript benchmarks/string_rescore.R %s\n", real_results_dir))
  cat(paste(rep("~", 78), collapse = ""), "\n")
  unlink(results_dir, recursive = TRUE)
}
cat("\n  Reminder: STRING connectivity is NOT a quality axis on its own. On\n")
cat("  IMvigor210 DGE reached 42x enrichment while predicting worst.\n")
