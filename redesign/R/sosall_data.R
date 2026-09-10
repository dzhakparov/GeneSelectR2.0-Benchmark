# ==============================================================================
#  SOS-ALL data loading and preprocessing, mirrored from the incumbent
#  benchmarks/sosall_benchmark_v2.0.R (log-CPM input, ENSG__SYMBOL collapse,
#  treatment parsing, location residualisation fitted on train only).
#
#  Provides:
#    load_sosall()                -> list(expression_matrix [samples x genes],
#                                    outcome [factor Healthy/AD], metadata)
#    sosall_residualise_split()   -> location residualisation, train-fitted
#
#  Pool selection and standardisation are NOT here: they happen per pool in
#  the benchmark runner (bio_prior.R), in the same order the incumbent uses
#  (residualise -> filter -> standardise).
# ==============================================================================

load_sosall <- function(expression_file = file.path("data",
                                                    "normalized_logcpm.csv"),
                        metadata_file = file.path("data", "metadata.csv")) {
  expression_long <- read.csv(expression_file, row.names = 1,
                              check.names = FALSE, stringsAsFactors = FALSE)
  gene_symbols <- sub(".*__", "", rownames(expression_long))
  expression_long$gene_symbol <- gene_symbols
  expression_collapsed <- aggregate(. ~ gene_symbol,
                                    data = expression_long, FUN = mean)
  rownames(expression_collapsed) <- expression_collapsed$gene_symbol
  expression_collapsed$gene_symbol <- NULL

  expression_matrix <- t(as.matrix(expression_collapsed))
  storage.mode(expression_matrix) <- "numeric"
  expression_matrix[is.na(expression_matrix)] <- 0

  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  treatment_parts <- strsplit(metadata[["treatment"]], "_", fixed = TRUE)
  metadata[["location"]]  <- vapply(treatment_parts, `[`, character(1), 1)
  metadata[["diagnosis"]] <- vapply(treatment_parts, function(parts) {
    if (length(parts) >= 2) paste(parts[-1], collapse = "_") else NA_character_
  }, character(1))

  common_samples <- intersect(rownames(expression_matrix), metadata[["X"]])
  expression_matrix <- expression_matrix[common_samples, , drop = FALSE]
  metadata <- metadata[match(common_samples, metadata[["X"]]), ,
                       drop = FALSE]

  #  Level order matters: level 2 is the event in every metric.
  outcome_factor <- factor(metadata[["diagnosis"]],
                           levels = c("Healthy", "AD"))
  if (any(is.na(outcome_factor))) {
    stop("SOS-ALL diagnosis contains values outside Healthy and AD.")
  }

  cat(sprintf("SOS-ALL: %d samples x %d genes (%s)\n",
              nrow(expression_matrix), ncol(expression_matrix),
              paste(table(outcome_factor), collapse = "/")))

  list(expression_matrix = expression_matrix, outcome = outcome_factor,
       metadata = metadata)
}

#  Residualise the location confounder out of train and test expression.
#  Coefficients are fitted on the training samples only (QR, as in the
#  incumbent); test is residualised with the training coefficients and the
#  training factor levels.
sosall_residualise_split <- function(expression_matrix, metadata,
                                     train_indices, test_indices) {
  build_design <- function(meta, reference_levels = NULL) {
    raw_values <- meta[["location"]]
    raw_values[is.na(raw_values) | raw_values == ""] <- "NA_level"
    if (is.null(reference_levels)) {
      levels_here <- levels(factor(raw_values))
    } else {
      levels_here <- reference_levels
    }
    #  Build the training contrast matrix manually. This retains the intercept
    #  and maps any test-only location to the training reference level without
    #  allowing model.matrix() to drop that test row.
    encoded <- matrix(1, nrow = length(raw_values), ncol = 1,
                      dimnames = list(NULL, "(Intercept)"))
    if (length(levels_here) > 1) {
      contrasts <- vapply(levels_here[-1],
                          function(level) as.numeric(raw_values == level),
                          numeric(length(raw_values)))
      colnames(contrasts) <- paste0("location", levels_here[-1])
      encoded <- cbind(encoded, contrasts)
    }
    list(design = encoded, levels = levels_here)
  }

  train_design <- build_design(metadata[train_indices, , drop = FALSE])
  coefs <- qr.coef(qr(train_design$design),
                   expression_matrix[train_indices, , drop = FALSE])
  #  qr.coef can leave NA rows for aliased columns; treat as zero effect.
  coefs[is.na(coefs)] <- 0

  train_resid <- expression_matrix[train_indices, , drop = FALSE] -
    train_design$design %*% coefs

  test_resid <- NULL
  if (!is.null(test_indices)) {
    test_design <- build_design(metadata[test_indices, , drop = FALSE],
                                reference_levels = train_design$levels)
    test_resid <- expression_matrix[test_indices, , drop = FALSE] -
      test_design$design %*% coefs
  }

  list(train = train_resid, test = test_resid)
}
