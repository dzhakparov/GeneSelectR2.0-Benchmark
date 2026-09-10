# ==============================================================================
#  Lock source-cohort panels before external holdout preparation
# ==============================================================================
#
#  This script converts repeated source-cohort rankings into one fixed panel per
#  method and panel size. It reads no holdout data. A gene receives reciprocal
#  rank 1/r when present in a cached top-200 ranking and 0 when absent. Mean
#  reciprocal rank across the shared source resamples defines the consensus.
#  Ties are resolved by appearance frequency, mean observed rank, then symbol.
#
#  Usage:
#      Rscript benchmarks/lock_validation_panels.R \
#        independent_benchmark_runs/validation/GSE107994/\
#          2026-08-20_kfold_trainfilter \
#        locked_panels/GSE107994_to_GSE19442.csv \
#        data/GSE19442/raw/GPL6947.annot.gz
# ==============================================================================


# ------------------------------------------------------------------------------
#  Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
  expected_source_accession = "GSE107994",
  holdout_accession         = "GSE19442",
  holdout_platform          = "GPL6947",
  # Every method has at least 100 source-ranked genes represented on GPL6947.
  # LASSO has 189 of its cached top 200, so k=200 cannot be matched exactly.
  panel_sizes               = c(10L, 20L, 50L, 100L),
  aggregation_rule = paste(
    "mean reciprocal rank across cached source resamples;",
    "absent gene score 0; ties by frequency, mean observed rank, gene"
  )
)


arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3) {
  stop(paste(
    "Usage: Rscript benchmarks/lock_validation_panels.R",
    "<SOURCE_RESULT_DIR> <OUTPUT_CSV> <HOLDOUT_GPL_ANNOTATION>"
  ))
}

source_result_dir <- arguments[1]
output_file <- arguments[2]
annotation_file <- arguments[3]
data_dir <- file.path(source_result_dir, "data")
rankings_file <- file.path(data_dir, "cached_rankings.rds")
config_file <- file.path(data_dir, "config.rds")

if (!file.exists(rankings_file) || !file.exists(config_file)) {
  stop("Source result directory must contain data/cached_rankings.rds and data/config.rds.")
}
if (!file.exists(annotation_file)) {
  stop("Holdout platform annotation is missing: ", annotation_file)
}


read_platform_gene_universe <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection))
  line_number <- 0L
  table_marker <- NA_integer_
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(line) == 0) break
    line_number <- line_number + 1L
    if (identical(line, "!platform_table_begin")) {
      table_marker <- line_number
      break
    }
  }
  if (is.na(table_marker)) stop("Platform annotation has no table marker.")

  annotation <- utils::read.delim(
    path, skip = table_marker, check.names = FALSE,
    stringsAsFactors = FALSE, quote = ""
  )
  if (!"Gene symbol" %in% colnames(annotation)) {
    stop("Platform annotation has no 'Gene symbol' column.")
  }
  symbols <- trimws(annotation[["Gene symbol"]])
  sort(unique(symbols[
    nzchar(symbols) & symbols != "---" & !grepl("///", symbols, fixed = TRUE)
  ]))
}


holdout_gene_universe <- read_platform_gene_universe(annotation_file)
if (length(holdout_gene_universe) < max(CONFIG$panel_sizes)) {
  stop("Holdout platform has too few unambiguous gene symbols.")
}

source_config <- readRDS(config_file)
if (!identical(source_config$accession, CONFIG$expected_source_accession)) {
  stop(sprintf("Expected source accession %s; config contains %s.",
               CONFIG$expected_source_accession, source_config$accession))
}
if (!identical(source_config$subsample_scheme, "kfold")) {
  stop("Panels must be locked from the kfold source run.")
}
if (!identical(as.integer(source_config$top_variable_genes), 2000L)) {
  stop("Panels must be locked from the p=2000 source run.")
}

cached_rankings <- readRDS(rankings_file)
if (!is.list(cached_rankings) || length(cached_rankings) == 0) {
  stop("cached_rankings.rds contains no methods.")
}


aggregate_method <- function(method_name, ranking_entries) {
  if (!is.list(ranking_entries) || length(ranking_entries) == 0) {
    stop(method_name, " has no cached ranking entries.")
  }
  ranked_lists <- lapply(ranking_entries, function(entry) {
    genes <- as.character(entry$ranked)
    genes <- genes[nzchar(genes) & !duplicated(genes)]
    genes
  })
  gene_universe <- sort(unique(unlist(ranked_lists, use.names = FALSE)))
  if (length(gene_universe) < max(CONFIG$panel_sizes)) {
    stop(sprintf("%s has only %d distinct ranked genes.",
                 method_name, length(gene_universe)))
  }

  reciprocal_rank <- matrix(0, nrow = length(gene_universe),
                            ncol = length(ranked_lists),
                            dimnames = list(gene_universe, NULL))
  observed_rank <- matrix(NA_real_, nrow = length(gene_universe),
                          ncol = length(ranked_lists),
                          dimnames = list(gene_universe, NULL))

  for (entry_index in seq_along(ranked_lists)) {
    genes <- ranked_lists[[entry_index]]
    positions <- seq_along(genes)
    row_indices <- match(genes, gene_universe)
    reciprocal_rank[row_indices, entry_index] <- 1 / positions
    observed_rank[row_indices, entry_index] <- positions
  }

  consensus <- data.frame(
    gene = gene_universe,
    mean_reciprocal_rank = rowMeans(reciprocal_rank),
    appearance_frequency = rowMeans(reciprocal_rank > 0),
    mean_observed_rank = rowMeans(observed_rank, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  consensus <- consensus[order(
    -consensus$mean_reciprocal_rank,
    -consensus$appearance_frequency,
    consensus$mean_observed_rank,
    consensus$gene
  ), , drop = FALSE]
  consensus <- consensus[consensus$gene %in% holdout_gene_universe, , drop = FALSE]
  if (nrow(consensus) < max(CONFIG$panel_sizes)) {
    stop(sprintf(
      "%s has only %d source-ranked genes measurable on %s.",
      method_name, nrow(consensus), CONFIG$holdout_platform
    ))
  }

  do.call(rbind, lapply(CONFIG$panel_sizes, function(panel_size) {
    selected <- consensus[seq_len(panel_size), , drop = FALSE]
    data.frame(
      method = method_name,
      panel_size = panel_size,
      rank = seq_len(panel_size),
      selected,
      stringsAsFactors = FALSE
    )
  }))
}


locked <- do.call(rbind, lapply(names(cached_rankings), function(method_name) {
  aggregate_method(method_name, cached_rankings[[method_name]])
}))

locked$source_accession <- source_config$accession
locked$holdout_accession <- CONFIG$holdout_accession
locked$holdout_platform <- CONFIG$holdout_platform
locked$holdout_annotation_md5 <- unname(tools::md5sum(annotation_file))
locked$source_result_dir <- source_result_dir
locked$source_rankings_md5 <- unname(tools::md5sum(rankings_file))
locked$locked_at_utc <- format(Sys.time(), tz = "UTC",
                               format = "%Y-%m-%dT%H:%M:%SZ")
locked$aggregation_rule <- CONFIG$aggregation_rule

column_order <- c(
  "source_accession", "holdout_accession", "source_result_dir",
  "source_rankings_md5", "holdout_platform", "holdout_annotation_md5",
  "locked_at_utc", "aggregation_rule",
  "method", "panel_size", "rank", "gene", "mean_reciprocal_rank",
  "appearance_frequency", "mean_observed_rank"
)
locked <- locked[, column_order]

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(locked, output_file, row.names = FALSE)

cat(sprintf("Locked %d methods at panel sizes %s from %d source resamples.\n",
            length(unique(locked$method)),
            paste(CONFIG$panel_sizes, collapse = ", "),
            length(cached_rankings[[1]])))
cat(sprintf("Source rankings MD5: %s\n", unique(locked$source_rankings_md5)))
cat(sprintf("%s annotation MD5: %s\n", CONFIG$holdout_platform,
            unique(locked$holdout_annotation_md5)))
cat(sprintf("Wrote %s (%d rows)\n", output_file, nrow(locked)))
