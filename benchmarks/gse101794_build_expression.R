# ==============================================================================
#  GSE101794 source builder: per-sample TPM files -> expression.tsv
# ==============================================================================
#
#  GEO distributes one two-column Gene/TPM file per participant. This script
#  assembles those files in metadata order and applies log2(TPM + 1). The result
#  belongs to the normalized-abundance branch of validation_prepare.R. TMM is
#  defined for count data and is not applied here.
#
#  Usage:
#      Rscript benchmarks/gse101794_build_expression.R
#      Rscript benchmarks/gse101794_build_expression.R --check-only
# ==============================================================================


# ------------------------------------------------------------------------------
#  Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
  accession       = "GSE101794",
  data_dir        = "data/GSE101794",
  metadata_file   = "data/GSE101794/metadata.tsv",
  raw_tar         = "data/GSE101794/raw/GSE101794_RAW.tar",
  extracted_dir   = "data/GSE101794/raw/tar",
  output_file     = "data/GSE101794/expression.tsv",
  expected_n      = 304L,
  expected_genes  = 13151L,
  expected_groups = c(CD = 254L, `Non-IBD` = 50L),
  comparison_tolerance = 5e-6
)


parse_arguments <- function(arguments) {
  unknown <- setdiff(arguments, "--check-only")
  if (length(unknown) > 0) {
    stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  }
  list(check_only = "--check-only" %in% arguments)
}


extract_raw_files <- function(config) {
  dir.create(config$extracted_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(config$extracted_dir,
                         pattern = "^GSM[0-9]+_.*\\.txt\\.gz$",
                         full.names = TRUE)
  if (length(existing) == config$expected_n) return(existing)

  if (!file.exists(config$raw_tar)) {
    stop("Raw GEO archive is missing: ", config$raw_tar)
  }

  cat(sprintf("Extracting %s\n", config$raw_tar))
  utils::untar(config$raw_tar, exdir = config$extracted_dir)

  list.files(config$extracted_dir,
             pattern = "^GSM[0-9]+_.*\\.txt\\.gz$",
             full.names = TRUE)
}


index_sample_files <- function(files, expected_n) {
  accessions <- sub("^(GSM[0-9]+)_.*$", "\\1", basename(files))
  malformed <- !grepl("^GSM[0-9]+$", accessions)
  if (any(malformed)) {
    stop("Could not parse GSM accession from: ",
         paste(basename(files[malformed]), collapse = ", "))
  }
  if (anyDuplicated(accessions)) {
    stop("More than one raw TPM file maps to the same GSM accession.")
  }
  if (length(accessions) != expected_n) {
    stop(sprintf("Expected %d raw files; found %d.",
                 expected_n, length(accessions)))
  }
  stats::setNames(files, accessions)
}


read_one_tpm <- function(path, expected_genes = NULL) {
  table <- utils::read.delim(path, check.names = FALSE,
                             stringsAsFactors = FALSE)
  if (!identical(colnames(table), c("Gene", "TPM"))) {
    stop(path, " must contain exactly the columns Gene and TPM.")
  }
  if (any(!is.finite(table$TPM)) || any(table$TPM < 0)) {
    stop(path, " contains non-finite or negative TPM values.")
  }
  # Kallisto output contains 55 symbols represented by more than one Gencode
  # feature. TPM is additive, so these rows are summed before log transformation.
  if (anyDuplicated(table$Gene)) {
    collapsed <- rowsum(table$TPM, group = table$Gene, reorder = FALSE)
    table <- data.frame(Gene = rownames(collapsed), TPM = collapsed[, 1],
                        row.names = NULL, stringsAsFactors = FALSE)
  }
  if (!is.null(expected_genes) && nrow(table) != expected_genes) {
    stop(sprintf("%s has %d genes; expected %d.",
                 path, nrow(table), expected_genes))
  }
  table
}


build_expression <- function(config) {
  if (!file.exists(config$metadata_file)) {
    stop("Metadata file is missing: ", config$metadata_file)
  }
  metadata <- utils::read.delim(config$metadata_file, check.names = FALSE,
                                stringsAsFactors = FALSE)
  required <- c("sample", "diagnosis")
  missing_columns <- setdiff(required, colnames(metadata))
  if (length(missing_columns) > 0) {
    stop("Metadata is missing: ", paste(missing_columns, collapse = ", "))
  }
  if (nrow(metadata) != config$expected_n || anyDuplicated(metadata$sample)) {
    stop("Metadata must contain 304 unique samples.")
  }

  observed_groups <- table(metadata$diagnosis)
  if (!identical(as.integer(observed_groups[names(config$expected_groups)]),
                 as.integer(config$expected_groups))) {
    stop("Diagnosis counts differ from 254 CD and 50 Non-IBD.")
  }

  files <- index_sample_files(extract_raw_files(config), config$expected_n)
  missing_files <- setdiff(metadata$sample, names(files))
  if (length(missing_files) > 0) {
    stop("Raw TPM files are missing for: ",
         paste(missing_files, collapse = ", "))
  }

  first <- read_one_tpm(files[[metadata$sample[1]]], config$expected_genes)
  genes <- first$Gene
  tpm <- matrix(NA_real_, nrow = length(genes), ncol = nrow(metadata),
                dimnames = list(genes, metadata$sample))
  tpm[, 1] <- first$TPM

  for (sample_index in seq.int(2L, nrow(metadata))) {
    sample_id <- metadata$sample[sample_index]
    current <- read_one_tpm(files[[sample_id]], config$expected_genes)
    if (!identical(current$Gene, genes)) {
      stop(files[[sample_id]], " has a different gene order or gene set.")
    }
    tpm[, sample_index] <- current$TPM
  }

  expression <- log2(tpm + 1)
  if (any(!is.finite(expression))) {
    stop("The assembled log2(TPM + 1) matrix contains non-finite values.")
  }
  expression
}


arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
expression <- build_expression(CONFIG)

cat(sprintf("Validated %d genes x %d samples; range [%.6f, %.6f]\n",
            nrow(expression), ncol(expression),
            min(expression), max(expression)))

if (arguments$check_only) {
  if (file.exists(CONFIG$output_file)) {
    existing <- as.matrix(utils::read.delim(
      CONFIG$output_file, row.names = 1, check.names = FALSE
    ))
    storage.mode(existing) <- "numeric"
    if (!identical(dim(existing), dim(expression)) ||
        !identical(rownames(existing), rownames(expression)) ||
        !identical(colnames(existing), colnames(expression))) {
      stop("Existing expression.tsv has different dimensions or identifiers.")
    }
    maximum_difference <- max(abs(existing - expression))
    if (maximum_difference > CONFIG$comparison_tolerance) {
      stop(sprintf(
        "Existing expression.tsv differs from the source build by %.8g.",
        maximum_difference
      ))
    }
    cat(sprintf("Existing expression.tsv reproduced within %.3g.\n",
                maximum_difference))
  }
  quit(status = 0)
}

utils::write.table(
  data.frame(gene = rownames(expression), expression, check.names = FALSE),
  file = CONFIG$output_file, sep = "\t", quote = FALSE,
  row.names = FALSE
)
cat(sprintf("Wrote %s\n", CONFIG$output_file))
