# Independent input and environment audit for the GeneSelectR benchmark suite.
# This script reads source data and code only. Historical result directories are
# intentionally outside its file inventory.

suppressPackageStartupMessages({
  library(jsonlite)
})

output_root <- Sys.getenv("GS_OUTPUT_ROOT", "independent_benchmark_runs")
audit_dir <- file.path(output_root, "preflight")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

source("benchmarks/validation_datasets.R")

rows <- list()
add_row <- function(dataset, status, n_samples = NA_integer_,
                    n_features = NA_integer_, negative = NA_integer_,
                    positive = NA_integer_, detail = "") {
  rows[[length(rows) + 1L]] <<- data.frame(
    dataset = dataset,
    status = status,
    n_samples = n_samples,
    n_features = n_features,
    negative = negative,
    positive = positive,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

# SOS-ALL: reproduce only the deterministic input parsing and alignment.
sos_expression <- read.csv("data/normalized_logcpm.csv", row.names = 1,
                           check.names = FALSE)
sos_metadata <- read.csv("data/metadata.csv", stringsAsFactors = FALSE)
sos_parts <- strsplit(sos_metadata$treatment, "_", fixed = TRUE)
sos_diagnosis <- vapply(sos_parts, function(x) paste(x[-1], collapse = "_"),
                        character(1))
sos_counts <- table(factor(sos_diagnosis, levels = c("Healthy", "AD")))
sos_ok <- !anyDuplicated(colnames(sos_expression)) &&
  !anyDuplicated(sos_metadata$X) &&
  all(sos_metadata$X %in% colnames(sos_expression)) &&
  all(is.finite(as.matrix(sos_expression)))
add_row("SOS-ALL", if (sos_ok) "ready" else "invalid",
        ncol(sos_expression), nrow(sos_expression),
        as.integer(sos_counts[1]), as.integer(sos_counts[2]),
        "Positive class fixed explicitly as AD.")
rm(sos_expression, sos_metadata)
gc()

# IMvigor210: inspect the primary Bioconductor resource from the local cache.
imvigor_status <- tryCatch({
  suppressPackageStartupMessages({
    library(ExperimentHub)
    library(SummarizedExperiment)
  })
  hub <- ExperimentHub(cache = "data/experimenthub_cache", localHub = TRUE)
  imvigor <- hub[["EH6677"]]
  response <- factor(as.character(colData(imvigor)$BOR),
                     levels = c("NR", "R"))
  response_counts <- table(response)
  counts_probe <- unlist(assay(imvigor, "counts")[seq_len(100),
                                                     seq_len(20)],
                         use.names = FALSE)
  valid_counts <- all(is.finite(counts_probe)) && all(counts_probe >= 0) &&
    all(abs(counts_probe - round(counts_probe)) < 1e-8)
  add_row("IMvigor210", if (valid_counts) "ready" else "invalid",
          ncol(imvigor), nrow(imvigor),
          as.integer(response_counts[1]), as.integer(response_counts[2]),
          "Bioconductor easierData resource EH6677; BOR NR/R.")
  rm(hub, imvigor)
  TRUE
}, error = function(e) {
  add_row("IMvigor210", "unavailable", detail = conditionMessage(e))
  FALSE
})
gc()

# Registered validation cohorts: verify exact class counts, sample alignment,
# finiteness, and unique feature/sample identifiers.
for (dataset in validation_datasets) {
  accession <- dataset$accession
  expression_filename <- if (dataset$scale == "counts") {
    "counts_prepared.csv"
  } else {
    "expression_prepared.csv"
  }
  expression_file <- file.path("data", accession, expression_filename)
  metadata_file <- file.path("data", accession, "metadata_prepared.csv")
  if (!file.exists(expression_file) || !file.exists(metadata_file)) {
    detail <- "Prepared expression or metadata is absent."
    if (identical(accession, "GSE57945")) {
      raw_probe <- read.delim(file.path("data", accession, "expression.tsv"),
                              nrows = 1000, check.names = FALSE)
      probe_values <- as.numeric(as.matrix(raw_probe[, -1, drop = FALSE]))
      fractional <- any(abs(probe_values - round(probe_values)) > 1e-8,
                        na.rm = TRUE)
      detail <- paste0(
        detail,
        " The supplied expression.tsv contains fractional values and is not ",
        "a raw-count matrix; the registered counts preprocessing is invalid ",
        "for this file. The published correction excludes 40 samples, and ",
        "their accession list is absent."
      )
      if (!fractional) {
        detail <- paste0(detail, " The limited scale probe was integer-valued.")
      }
    }
    add_row(accession, "unavailable", detail = detail)
    next
  }

  expression <- t(as.matrix(read.csv(expression_file, row.names = 1,
                                     check.names = FALSE)))
  storage.mode(expression) <- "numeric"
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE,
                       check.names = FALSE)
  outcome <- factor(metadata$outcome,
                    levels = c(dataset$outcome$negative_label,
                               dataset$outcome$positive_label))
  class_counts <- table(outcome)
  expected <- unname(as.integer(dataset$expected_n))
  observed <- unname(as.integer(class_counts))
  valid <- identical(rownames(expression), metadata$sample_id) &&
    !anyDuplicated(rownames(expression)) &&
    !anyDuplicated(colnames(expression)) &&
    all(is.finite(expression)) &&
    (dataset$scale != "counts" ||
       (all(expression >= 0) &&
        all(abs(expression - round(expression)) <= 1e-8))) &&
    !anyNA(outcome) && identical(observed, expected)
  add_row(accession, if (valid) "ready" else "invalid",
          nrow(expression), ncol(expression), observed[1], observed[2],
          if (valid) "Prepared data match the registered cohort." else
            "At least one alignment, finiteness, uniqueness, or cohort-size check failed.")
  rm(expression, metadata)
  gc()
}

audit_summary <- do.call(rbind, rows)
write.csv(audit_summary, file.path(audit_dir, "dataset_audit.csv"),
          row.names = FALSE)

# Hash all benchmark inputs and executable code used in the independent run.
input_files <- c(
  "data/normalized_logcpm.csv",
  "data/metadata.csv",
  "benchmarks/sosall_benchmark_v2.0.R",
  "imvigor210_benchmark_v2.0.R",
  "benchmarks/validation_benchmark.R",
  "benchmarks/count_preprocessing.R",
  "benchmarks/independent_preflight.R",
  "benchmarks/independent_analyse_results.R",
  "benchmarks/validation_datasets.R",
  "package/dist/GeneSelectR_0.99.0.tar.gz"
)

# External biological evidence is part of the benchmark input. Record the
# complete frozen GeneSelectR cache because it contains GO annotations,
# information-content values, target-term similarities, and disease-specific
# Open Targets seeds. Only the human MSigDB collections queried by the runner
# are included. STRING mapping and graph files are recorded separately.
external_input_files <- c(
  list.files("data/r_user_cache/R/GeneSelectR", full.names = TRUE,
             recursive = TRUE),
  file.path("data/r_user_cache/R/msigdbr", c(
    "msigdb.2026.1.Hs.C2.rds",
    "msigdb.2026.1.Hs.C7.rds",
    "msigdb.2026.1.Hs.H.rds"
  )),
  list.files("data/string_db_cache", full.names = TRUE, recursive = TRUE)
)
input_files <- c(input_files,
                 external_input_files[file.exists(external_input_files)])
for (dataset in validation_datasets) {
  candidates <- file.path("data", dataset$accession,
                          c("expression_prepared.csv", "counts_prepared.csv",
                            "metadata_prepared.csv",
                            "expression.tsv", "metadata.tsv"))
  input_files <- c(input_files, candidates[file.exists(candidates)])
}
input_files <- unique(input_files[file.exists(input_files)])
hashes <- data.frame(
  file = input_files,
  bytes = as.numeric(file.info(input_files)$size),
  md5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)
write.csv(hashes, file.path(audit_dir, "input_hashes.csv"), row.names = FALSE)

environment_record <- list(
  timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string,
  platform = R.version$platform,
  geneselectr_version = as.character(utils::packageVersion("GeneSelectR")),
  geneselectr_library = find.package("GeneSelectR"),
  r_user_cache_dir = Sys.getenv("R_USER_CACHE_DIR", ""),
  string_version = "12.0",
  string_confidence_threshold = 400L,
  output_root = normalizePath(output_root, mustWork = FALSE),
  global_worker_cap = Sys.getenv("GS_GLOBAL_WORKER_CAP", ""),
  top_variable_genes = Sys.getenv("GS_TOP_VARIABLE_GENES", "2000"),
  memory_fraction = Sys.getenv("GS_MEMORY_FRACTION", "0.55")
)
write_json(environment_record, file.path(audit_dir, "environment.json"),
           auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(audit_dir, "session_info.txt"))

print(audit_summary, row.names = FALSE)
cat("\nIndependent preflight written to:", normalizePath(audit_dir), "\n")
