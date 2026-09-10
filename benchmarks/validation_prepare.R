# ==============================================================================
#  GeneSelectR 2.0 — validation cohort preparation
# ==============================================================================
#
#  Turns what GEO shipped into what the benchmark reads, for one of the five
#  pre-registered validation cohorts:
#
#      data/<ACC>/expression.tsv       ->  data/<ACC>/expression_prepared.csv
#                                        or counts_prepared.csv for RNA-seq
#      data/<ACC>/metadata.tsv         ->  data/<ACC>/metadata_prepared.csv
#                                          data/<ACC>/prep_manifest.json
#
#  USAGE
#      Rscript benchmarks/validation_prepare.R GSE69683
#      Rscript benchmarks/validation_prepare.R all
#
#  The analysis plan -- filters, contrast, grouping, confounders -- comes from
#  benchmarks/validation_datasets.R and is not repeated here.
#
#  WHAT THIS SCRIPT IS FOR
#  -----------------------
#  Two cohorts are RNA-seq counts and three are normalised microarray. They
#  cannot share a preprocessing path: TMM-normalising an already-logged array
#  matrix produces numbers that are not wrong in any way the run would notice.
#  So the branch is explicit, and the declared scale is VERIFIED against the
#  matrix rather than trusted:
#
#      counts            -> verified integer counts retained for fold-local
#                           filtering, TMM and log2-CPM in the benchmark
#      log_normalized    -> used as-is (variance filtering happens downstream)
#      linear_normalized -> log2(x + 1)
#      linear_background_corrected -> floor negative background estimates at
#                           zero, then log2(x + 1)
#
#  The verification is the important half. A matrix distributed as "counts"
#  that actually contains TPM or RPKM is the classic silent disaster in this
#  corner of the field: TMM assumes integer library-size-scaled counts, applying
#  it to length-normalised values is wrong rather than merely suboptimal, and
#  every downstream step accepts the result without complaint. GSE57945 in
#  particular distributes both a counts archive and an RPKM table, so this is a
#  live hazard and not a hypothetical one. Mismatches STOP the run.
#
#  ORDER OF OPERATIONS, AND WHY
#  ----------------------------
#  Samples are filtered before files are written. Count normalization remains
#  deferred because its abundance threshold and TMM reference must be estimated
#  independently in every training fold. Pointwise array transformations are
#  applied during preparation because they do not depend on other samples.
# ==============================================================================


suppressPackageStartupMessages({
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

source("benchmarks/validation_datasets.R")


# ------------------------------------------------------------------------------
#  Scale detection
# ------------------------------------------------------------------------------
#
#  Classifies a matrix from its own values, with no reference to what the file
#  was called. The three signatures are far enough apart that this is reliable:
#
#    counts            all values whole numbers, none negative, large maximum
#    log_normalized    small dynamic range (log2 of expression tops out in the
#                      tens)
#    linear_background_corrected non-integer array signal with negative values
#                      and a large maximum
#    linear_normalized non-negative, non-integer, large maximum -- this is what
#                      TPM/RPKM/FPKM look like, and it is the case that must
#                      never be mistaken for counts
#
#  `all_integers` tolerates floating-point representation of whole numbers
#  (2.0000000001 from a round-trip through a text file is still a count) but
#  not genuine fractions.
detect_expression_scale <- function(expression_matrix) {

  finite_values <- expression_matrix[is.finite(expression_matrix)]
  if (length(finite_values) == 0) stop("Expression matrix has no finite values.")

  # Sample rather than test all 20k x 500 values; the signature is a property of
  # the distribution, not of any particular cell.
  probe <- if (length(finite_values) > 2e6) {
    sample(finite_values, 2e6)
  } else {
    finite_values
  }

  all_integers   <- all(abs(probe - round(probe)) < 1e-8)
  has_negatives  <- any(probe < 0)
  maximum_value  <- max(probe)

  detected <- if (all_integers && !has_negatives && maximum_value > 100) {
    "counts"
  } else if (has_negatives && maximum_value >= 100) {
    "linear_background_corrected"
  } else if (maximum_value < 30) {
    "log_normalized"
  } else {
    "linear_normalized"
  }

  list(detected = detected, all_integers = all_integers,
       has_negatives = has_negatives, maximum = maximum_value,
       minimum = min(probe))
}


# Compares the registry's declaration against the matrix and stops on a
# mismatch. The message names the specific hazard rather than saying "mismatch",
# because the correct response differs: a counts/linear swap means you
# downloaded the wrong supplementary file, while a log/linear swap means the
# submitter logged (or did not log) against the series documentation.
verify_expression_scale <- function(declared, evidence, accession) {

  cat(sprintf("  scale declared: %-18s detected: %-18s\n",
              declared, evidence$detected))
  cat(sprintf("    range [%.3f, %.3f] | all integers: %s | negatives: %s\n",
              evidence$minimum, evidence$maximum,
              evidence$all_integers, evidence$has_negatives))

  if (identical(declared, evidence$detected)) return(invisible(TRUE))

  if (declared == "counts" && evidence$detected == "linear_normalized") {
    stop(sprintf(paste0(
      "%s: declared 'counts' but the matrix is NOT integer (max %.2f, ",
      "fractional values present).\n",
      "  This is almost certainly a TPM/RPKM/FPKM table, not counts.\n",
      "  edgeR TMM on length-normalised values is WRONG, not merely ",
      "suboptimal, and nothing downstream would report it.\n",
      "  Download the raw counts matrix instead, or change `scale` in the ",
      "registry to 'linear_normalized' and accept that this cohort has no ",
      "counts-based normalisation."), accession, evidence$maximum))
  }

  if (declared == "counts" && evidence$detected == "log_normalized") {
    stop(sprintf(paste0(
      "%s: declared 'counts' but the matrix looks already log-transformed ",
      "(max %.2f%s).\n",
      "  Re-normalising it would compound two transformations. Point the ",
      "registry at the raw counts, or set `scale` to 'log_normalized'."),
      accession, evidence$maximum,
      if (evidence$has_negatives) ", contains negatives" else ""))
  }

  stop(sprintf(paste0(
    "%s: scale mismatch -- registry declares '%s', the matrix looks like '%s' ",
    "(range [%.2f, %.2f], all integers: %s).\n",
    "  Resolve this before running: the two paths are not interchangeable."),
    accession, declared, evidence$detected,
    evidence$minimum, evidence$maximum, evidence$all_integers))
}


# ------------------------------------------------------------------------------
#  Platform branch
# ------------------------------------------------------------------------------
#
#  Returns the prepared matrix plus a record of the transformation. Counts stay
#  on their integer scale for fold-local processing in validation_benchmark.R.

normalize_expression <- function(expression_matrix, scale, accession) {

  if (scale == "counts") {
    if (any(!is.finite(expression_matrix)) || any(expression_matrix < 0) ||
        any(abs(expression_matrix - round(expression_matrix)) > 1e-8)) {
      stop(sprintf(
        "%s: count preparation requires finite non-negative integers.",
        accession))
    }
    cat("  retaining integer counts for fold-local normalization\n")
    return(list(
      matrix = expression_matrix,
      method = paste(
        "integer counts retained; training-fold abundance filter,",
        "fixed-reference TMM and log2-CPM deferred to benchmark"
      ),
      n_genes_filtered = 0L
    ))
  }

  if (scale == "linear_normalized") {
    cat("  applying log2(x + 1)\n")
    return(list(matrix = log2(expression_matrix + 1),
                method = "log2(x + 1)",
                n_genes_filtered = 0L))
  }

  if (scale == "linear_background_corrected") {
    n_nonpositive <- sum(expression_matrix <= 0)
    cat(sprintf(
      "  flooring %d non-positive background estimates at zero; applying log2(x + 1)\n",
      n_nonpositive
    ))
    return(list(
      matrix = log2(pmax(expression_matrix, 0) + 1),
      method = paste(
        "non-positive background estimates floored at zero;",
        "log2(x + 1)"
      ),
      n_genes_filtered = 0L
    ))
  }

  # Already normalised and logged by the submitter. Deliberately untouched: the
  # array series in this set are RMA/MAS5 outputs, and re-normalising them would
  # be a second, unjustified transformation. Downstream the benchmark applies
  # its own unsupervised variance filter, identically for every method.
  cat("  already log-normalised; used as-is\n")
  list(matrix = expression_matrix,
       method = "none (submitter-normalised log scale)",
       n_genes_filtered = 0L)
}


# ------------------------------------------------------------------------------
#  Metadata pipeline
# ------------------------------------------------------------------------------
#
#  derive -> author exclusions -> registry filters -> outcome mapping.
#  Each stage reports how many samples it removed, so a cohort that comes out
#  the wrong size can be traced to the stage that shrank it.

apply_sample_selection <- function(metadata, dataset) {

  n_start <- nrow(metadata)

  # --- Derived columns ------------------------------------------------------
  if (!is.null(dataset$derive)) {
    metadata <- dataset$derive(metadata)
    cat("  derived columns added\n")
  }

  # --- Author erratum exclusions -------------------------------------------
  # Declared but missing is a hard stop. See the registry entry for GSE57945:
  # silently running the uncorrected cohort while the log claims it is
  # corrected is the failure mode worth being noisy about.
  if (!is.null(dataset$excluded_samples_file)) {
    if (!file.exists(dataset$excluded_samples_file)) {
      stop(sprintf(paste0(
        "%s declares an author correction but the exclusion list is missing:\n",
        "    %s\n",
        "  Create it with one sample accession (GSM...) per line, taken from ",
        "the published erratum.\n",
        "  Running without it would produce the UNCORRECTED cohort while the ",
        "log claims otherwise."),
        dataset$accession, dataset$excluded_samples_file))
    }
    excluded <- trimws(readLines(dataset$excluded_samples_file, warn = FALSE))
    excluded <- excluded[nzchar(excluded) & !startsWith(excluded, "#")]

    n_before <- nrow(metadata)
    metadata <- metadata[!metadata$sample %in% excluded, , drop = FALSE]
    cat(sprintf("  author correction: dropped %d of %d listed exclusions (%d -> %d samples)\n",
                n_before - nrow(metadata), length(excluded),
                n_before, nrow(metadata)))
  }

  # --- Registry filters -----------------------------------------------------
  for (filter_spec in dataset$sample_filters) {
    if (!filter_spec$column %in% colnames(metadata)) {
      stop(sprintf("Filter column '%s' not found. Available: %s",
                   filter_spec$column,
                   paste(colnames(metadata), collapse = ", ")))
    }
    n_before <- nrow(metadata)
    metadata <- metadata[metadata[[filter_spec$column]] %in% filter_spec$keep, ,
                         drop = FALSE]
    cat(sprintf("  filter %s in {%s}: %d -> %d samples\n",
                filter_spec$column, paste(filter_spec$keep, collapse = ", "),
                n_before, nrow(metadata)))
  }

  # --- Outcome mapping ------------------------------------------------------
  outcome_spec <- dataset$outcome
  if (!outcome_spec$column %in% colnames(metadata)) {
    stop(sprintf("Outcome column '%s' not found. Available: %s",
                 outcome_spec$column,
                 paste(colnames(metadata), collapse = ", ")))
  }

  raw_outcome <- as.character(metadata[[outcome_spec$column]])
  outcome <- rep(NA_character_, length(raw_outcome))
  outcome[raw_outcome %in% outcome_spec$negative] <- outcome_spec$negative_label
  outcome[raw_outcome %in% outcome_spec$positive] <- outcome_spec$positive_label

  # Report what was dropped by VALUE, not just by count. An unmapped value is
  # either a class we deliberately excluded or a spelling we did not anticipate,
  # and those need different responses.
  unmapped <- table(raw_outcome[is.na(outcome)])
  if (length(unmapped) > 0) {
    cat("  unmapped outcome values (dropped):\n")
    for (value in names(unmapped)) {
      cat(sprintf("      %-40s %d\n", value, unmapped[[value]]))
    }
  }

  metadata <- metadata[!is.na(outcome), , drop = FALSE]
  outcome  <- outcome[!is.na(outcome)]

  # Negative level FIRST. Every metric downstream reads level 2 as the positive
  # class, so getting this backwards inverts every AUC in the run.
  metadata$outcome <- factor(outcome, levels = c(outcome_spec$negative_label,
                                                 outcome_spec$positive_label))

  # --- Grouping -------------------------------------------------------------
  if (!is.null(dataset$group_column)) {
    if (!dataset$group_column %in% colnames(metadata)) {
      stop(sprintf("Group column '%s' not found. Available: %s",
                   dataset$group_column,
                   paste(colnames(metadata), collapse = ", ")))
    }
    metadata$group_id <- as.character(metadata[[dataset$group_column]])
    if (any(is.na(metadata$group_id) | !nzchar(metadata$group_id))) {
      stop(sprintf(paste0(
        "%s: group column '%s' has missing values. Grouped CV would silently ",
        "degrade to ordinary CV for those samples."),
        dataset$accession, dataset$group_column))
    }
    group_sizes <- table(metadata$group_id)
    cat(sprintf("  grouping on %s: %d groups, sizes %d-%d\n",
                dataset$group_column, length(group_sizes),
                min(group_sizes), max(group_sizes)))
  } else {
    # One group per sample. Downstream code can then treat every dataset the
    # same way instead of branching on whether grouping applies.
    metadata$group_id <- metadata$sample
  }

  cat(sprintf("  sample selection: %d -> %d samples\n", n_start, nrow(metadata)))
  metadata
}


# Compares the realised class counts against the registry's pre-registered
# arithmetic. A warning rather than a stop: the data can legitimately differ
# from what was written down (GEO re-releases, a corrected exclusion list), but
# it must never differ silently, because the pre-registration is only worth
# something if a deviation from it is visible in the log.
check_expected_n <- function(metadata, dataset) {
  observed <- table(metadata$outcome)
  expected <- dataset$expected_n

  observed_negative <- as.integer(observed[[dataset$outcome$negative_label]])
  observed_positive <- as.integer(observed[[dataset$outcome$positive_label]])

  cat(sprintf("  class counts: %d %s / %d %s   (pre-registered: %d / %d)\n",
              observed_positive, dataset$outcome$positive_label,
              observed_negative, dataset$outcome$negative_label,
              expected[["positive"]], expected[["negative"]]))

  if (observed_negative != expected[["negative"]] ||
      observed_positive != expected[["positive"]]) {
    warning(sprintf(paste0(
      "%s: class counts differ from the pre-registered plan ",
      "(got %d/%d, expected %d/%d). This is not necessarily an error, but it ",
      "IS a deviation and belongs in the paper's methods section."),
      dataset$accession, observed_positive, observed_negative,
      expected[["positive"]], expected[["negative"]]), call. = FALSE)
  }

  if (min(observed_negative, observed_positive) < 20) {
    warning(sprintf(paste0(
      "%s: minority class has only %d samples. Stratified 5-fold CV leaves ",
      "~%d in each test fold, so per-fold AUC will be extremely noisy."),
      dataset$accession, min(observed_negative, observed_positive),
      round(min(observed_negative, observed_positive) / 5)), call. = FALSE)
  }

  invisible(NULL)
}


# ------------------------------------------------------------------------------
#  One dataset, end to end
# ------------------------------------------------------------------------------

prepare_dataset <- function(accession) {

  dataset <- get_validation_dataset(accession)

  locked_panels <- NULL
  if (dataset$analysis_role == "external_holdout") {
    if (!file.exists(dataset$locked_panel_file)) {
      stop(sprintf(paste0(
        "%s is an external holdout and its panel manifest is missing:\n",
        "    %s\n",
        "  Run lock_validation_panels.R using %s before preparing the holdout."),
        accession, dataset$locked_panel_file, dataset$source_accession))
    }
    locked_panels <- utils::read.csv(dataset$locked_panel_file,
                                     stringsAsFactors = FALSE)
    required_lock_fields <- c("source_accession", "holdout_accession",
                              "source_rankings_md5", "holdout_platform",
                              "holdout_annotation_md5", "locked_at_utc")
    if (length(setdiff(required_lock_fields, colnames(locked_panels))) > 0 ||
        !identical(unique(locked_panels$source_accession),
                   dataset$source_accession) ||
        !identical(unique(locked_panels$holdout_accession), accession)) {
      stop(accession, ": locked panel manifest has invalid provenance.")
    }
    platform <- unique(locked_panels$holdout_platform)
    annotation_path <- file.path(
      "data", accession, "raw", paste0(platform, ".annot.gz")
    )
    if (!file.exists(annotation_path)) {
      stop(accession, ": holdout platform annotation is missing: ",
           annotation_path)
    }
    observed_annotation_md5 <- unname(tools::md5sum(annotation_path))
    if (!identical(unique(locked_panels$holdout_annotation_md5),
                   observed_annotation_md5)) {
      stop(accession, ": platform annotation differs from the panel lock.")
    }
  }

  cat(sprintf("\n%s\n", paste(rep("=", 78), collapse = "")))
  print_validation_dataset(dataset)
  cat(sprintf("%s\n", paste(rep("-", 78), collapse = "")))

  data_dir        <- file.path("data", accession)
  expression_path <- file.path(data_dir, "expression.tsv")
  metadata_path   <- file.path(data_dir, "metadata.tsv")

  # The two RNA-seq cohorts ship expression in a supplementary file rather than
  # the series matrix, so a missing expression.tsv here is expected for them and
  # the message says what to do about it.
  if (!file.exists(expression_path)) {
    stop(sprintf(paste0(
      "%s: %s not found.\n",
      "  Fetch it first:   python fetch_geo_data.py %s\n",
      "  If that reports 'Expression is in a supplementary file', download the ",
      "counts matrix named in data/%s/manifest.json and write it to %s as a ",
      "tab-separated genes x samples table with gene SYMBOLS in column 1."),
      accession, expression_path, accession, accession, expression_path))
  }
  if (!file.exists(metadata_path)) {
    stop(sprintf("%s: %s not found. Run: python fetch_geo_data.py %s",
                 accession, metadata_path, accession))
  }

  # --- Load ----------------------------------------------------------------
  expression_table <- utils::read.delim(expression_path, row.names = 1,
                                        check.names = FALSE,
                                        stringsAsFactors = FALSE)
  expression_matrix <- as.matrix(expression_table)
  storage.mode(expression_matrix) <- "numeric"

  metadata <- utils::read.delim(metadata_path, check.names = FALSE,
                                stringsAsFactors = FALSE)
  if (!"sample" %in% colnames(metadata)) {
    stop(sprintf("%s: metadata has no 'sample' column.", accession))
  }

  cat(sprintf("  loaded %d genes x %d samples\n",
              nrow(expression_matrix), ncol(expression_matrix)))

  # --- Select samples -------------------------------------------------------
  metadata <- apply_sample_selection(metadata, dataset)

  common_samples <- intersect(colnames(expression_matrix), metadata$sample)
  if (length(common_samples) < nrow(metadata)) {
    cat(sprintf("  WARNING: %d selected samples have no expression column\n",
                nrow(metadata) - length(common_samples)))
  }
  if (length(common_samples) < 20) {
    stop(sprintf("%s: only %d samples survive selection and matching.",
                 accession, length(common_samples)))
  }

  metadata <- metadata[match(common_samples, metadata$sample), , drop = FALSE]
  expression_matrix <- expression_matrix[, common_samples, drop = FALSE]

  check_expected_n(metadata, dataset)

  # --- Scale check, then the platform branch -------------------------------
  # Verified on the ANALYSED subset: dropping healthy controls or a smoking arm
  # can change the matrix's character, and the subset is what gets normalised.
  scale_evidence <- detect_expression_scale(expression_matrix)
  verify_expression_scale(dataset$scale, scale_evidence, accession)

  normalization <- normalize_expression(expression_matrix, dataset$scale,
                                        accession)
  expression_matrix <- normalization$matrix

  # --- Collapse duplicate symbols ------------------------------------------
  # Counts for a shared gene symbol are additive. Normalized array probes retain
  # the highest-mean probe, matching the source preparation convention.
  gene_symbols <- rownames(expression_matrix)
  if (anyDuplicated(gene_symbols)) {
    n_duplicates <- sum(duplicated(gene_symbols))
    if (dataset$scale == "counts") {
      expression_matrix <- rowsum(expression_matrix, group = gene_symbols,
                                   reorder = FALSE)
      cat(sprintf("  collapsed %d duplicate gene symbols by summed counts\n",
                  n_duplicates))
    } else {
      row_means <- rowMeans(expression_matrix, na.rm = TRUE)
      keep_rows <- order(row_means, decreasing = TRUE)
      keep_rows <- keep_rows[!duplicated(gene_symbols[keep_rows])]
      cat(sprintf("  collapsed %d duplicate gene symbols (highest mean kept)\n",
                  nrow(expression_matrix) - length(keep_rows)))
      expression_matrix <- expression_matrix[sort(keep_rows), , drop = FALSE]
    }
  }

  # Genes that are constant across the analysed samples carry no information and
  # break standardisation (sd = 0). Removing them here keeps the benchmark's
  # variance filter honest -- otherwise "top 2000 by variance" silently starts
  # from a pool padded with zero-variance rows.
  gene_variances <- apply(expression_matrix, 1, stats::var, na.rm = TRUE)
  n_constant <- sum(!is.finite(gene_variances) | gene_variances == 0)
  if (n_constant > 0) {
    expression_matrix <- expression_matrix[is.finite(gene_variances) &
                                             gene_variances > 0, , drop = FALSE]
    cat(sprintf("  dropped %d zero-variance genes\n", n_constant))
  }

  # Complete-cohort mean imputation would expose held-out values to training
  # folds. Source matrices with non-finite values require an explicit imputation
  # protocol and are therefore rejected here.
  n_missing <- sum(!is.finite(expression_matrix))
  if (n_missing > 0) {
    stop(sprintf("%s contains %d non-finite expression values.",
                 accession, n_missing))
  }

  # --- Write ----------------------------------------------------------------
  # The benchmark reads exactly these columns and nothing else, so the confounder
  # columns named in the registry are carried through explicitly and the other
  # ~40 GEO boilerplate columns are dropped.
  output_metadata <- data.frame(
    sample_id = metadata$sample,
    outcome   = metadata$outcome,
    group_id  = metadata$group_id,
    stringsAsFactors = FALSE
  )
  for (confounder in dataset$confounders) {
    if (!confounder %in% colnames(metadata)) {
      stop(sprintf("Confounder column '%s' not found in %s metadata.",
                   confounder, accession))
    }
    values <- as.character(metadata[[confounder]])
    values[is.na(values) | !nzchar(values)] <- "NA_level"
    output_metadata[[confounder]] <- values
  }

  expression_filename <- if (dataset$scale == "counts") {
    "counts_prepared.csv"
  } else {
    "expression_prepared.csv"
  }
  expression_out <- file.path(data_dir, expression_filename)
  metadata_out   <- file.path(data_dir, "metadata_prepared.csv")

  utils::write.csv(as.data.frame(expression_matrix), expression_out,
                   row.names = TRUE)
  utils::write.csv(output_metadata, metadata_out, row.names = FALSE)

  manifest <- list(
    accession        = accession,
    label            = dataset$label,
    prepared_on      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    declared_scale   = dataset$scale,
    detected_scale   = scale_evidence$detected,
    normalization    = normalization$method,
    prepared_expression_file = expression_filename,
    normalization_scope = if (dataset$scale == "counts")
      "training_fold_fixed_reference_tmm" else "prepared_input",
    n_genes          = nrow(expression_matrix),
    n_samples        = ncol(expression_matrix),
    n_positive       = sum(output_metadata$outcome ==
                             dataset$outcome$positive_label),
    n_negative       = sum(output_metadata$outcome ==
                             dataset$outcome$negative_label),
    positive_label   = dataset$outcome$positive_label,
    negative_label   = dataset$outcome$negative_label,
    n_groups         = length(unique(output_metadata$group_id)),
    group_column     = dataset$group_column %||% NA_character_,
    confounders      = dataset$confounders,
    analysis_role    = dataset$analysis_role,
    source_accession = dataset$source_accession %||% NA_character_,
    locked_panel_file = dataset$locked_panel_file %||% NA_character_,
    locked_panel_md5 = if (is.null(dataset$locked_panel_file)) NA_character_
      else unname(tools::md5sum(dataset$locked_panel_file)),
    locked_at_utc = if (is.null(locked_panels)) NA_character_
      else unique(locked_panels$locked_at_utc),
    disease_term     = dataset$disease_term,
    target_go_terms  = dataset$target_go_terms,
    expected_n       = as.list(dataset$expected_n)
  )
  writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
             file.path(data_dir, "prep_manifest.json"))

  cat(sprintf("  wrote %s (%d genes x %d samples)\n",
              expression_out, nrow(expression_matrix), ncol(expression_matrix)))
  cat(sprintf("  wrote %s\n", metadata_out))

  invisible(manifest)
}


# ------------------------------------------------------------------------------
#  Entry point
# ------------------------------------------------------------------------------

command_args <- commandArgs(trailingOnly = TRUE)
if (length(command_args) == 0) {
  cat("Usage: Rscript benchmarks/validation_prepare.R <ACCESSION|all>\n\n")
  cat("Registered datasets:\n")
  for (name in names(validation_datasets)) {
    cat(sprintf("  %-12s %s\n", name, validation_datasets[[name]]$label))
  }
  quit(status = 1)
}

accessions <- if (identical(command_args[1], "all")) {
  names(validation_datasets)
} else {
  command_args
}

# One cohort failing must not stop the others: with five datasets and two of
# them needing a manual supplementary download, a hard stop on the first would
# hide the state of the remaining four. Failures are collected and reprinted at
# the end so the summary is a to-do list.
prep_results <- list()
for (accession in accessions) {
  outcome <- tryCatch(prepare_dataset(accession),
                      error = function(e) {
                        cat(sprintf("\n  !! %s FAILED: %s\n",
                                    accession, conditionMessage(e)))
                        structure(conditionMessage(e), class = "prep_failure")
                      })
  prep_results[[accession]] <- outcome
}

cat(sprintf("\n%s\n", paste(rep("=", 78), collapse = "")))
cat("Preparation summary\n")
cat(sprintf("%s\n", paste(rep("=", 78), collapse = "")))
for (accession in names(prep_results)) {
  result <- prep_results[[accession]]
  if (inherits(result, "prep_failure")) {
    cat(sprintf("  %-12s NOT READY  -- %s\n", accession,
                strsplit(as.character(result), "\n")[[1]][1]))
  } else {
    cat(sprintf("  %-12s ready      %d genes x %d samples (%d %s / %d %s)\n",
                accession, result$n_genes, result$n_samples,
                result$n_positive, result$positive_label,
                result$n_negative, result$negative_label))
  }
}
ready_results <- prep_results[!vapply(
  prep_results, inherits, logical(1), what = "prep_failure"
)]
ready_roles <- vapply(ready_results, function(result) result$analysis_role,
                      character(1))
if (any(ready_roles == "nested_cv")) {
  cat("\nRun a prepared nested-CV cohort with:\n")
  cat("  Rscript benchmarks/validation_benchmark.R <ACCESSION> kfold\n")
}
if (any(ready_roles == "external_holdout")) {
  cat("\nRun the external locked-panel evaluation with:\n")
  cat(paste(
    "  Rscript benchmarks/validation_holdout.R GSE19442",
    "locked_panels/GSE107994_to_GSE19442.csv\n"
  ))
}
