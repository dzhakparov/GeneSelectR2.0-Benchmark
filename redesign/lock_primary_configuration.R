#!/usr/bin/env Rscript
# Lock the primary base GeneSelectR configuration after the deterministic
# development comparison. This panel is intended for a future external cohort.
# GSE19442 was accessed before this lock and is therefore ineligible as an
# independent test of this exact configuration.

results_root <- file.path("redesign", "results_corrected")
report_dir <- file.path(results_root, "full_benchmark_deterministic")
source_dir <- file.path(results_root, "validation_benchmark", "GSE107994")
output_file <- file.path(
  "locked_panels", "GSE107994_GS_full_ungrouped_primary_2026-08-29.csv")
source_manifest_file <- sub("[.]csv$", "_sources.csv", output_file)
panel_sizes <- c(10L, 20L, 50L, 100L)

comparison <- read.csv(file.path(report_dir, "cross_dataset_method_summary.csv"),
                       stringsAsFactors = FALSE)
base <- comparison[comparison$arm %in%
                     c("GS_full_grouped", "GS_full_ungrouped"), ]
primary <- base$arm[which.min(base$mean_dataset_rank)]
if (!identical(primary, "GS_full_ungrouped")) {
  stop("The deterministic comparison no longer selects GS_full_ungrouped.",
       call. = FALSE)
}

ranking_files <- file.path(source_dir, sprintf(
  "ranking_r%d_f%d_GS_full_ungrouped.csv",
  rep(1:3, each = 5), rep(1:5, times = 3)
))
if (!all(file.exists(ranking_files))) {
  stop("One or more GSE107994 source rankings are missing.", call. = FALSE)
}
rankings <- lapply(ranking_files, function(path) {
  genes <- read.csv(path, stringsAsFactors = FALSE)$gene
  genes[nzchar(genes) & !duplicated(genes)]
})
gene_universe <- sort(unique(unlist(rankings, use.names = FALSE)))
reciprocal_rank <- matrix(0, nrow = length(gene_universe),
                          ncol = length(rankings),
                          dimnames = list(gene_universe, NULL))
observed_rank <- matrix(NA_real_, nrow = length(gene_universe),
                        ncol = length(rankings),
                        dimnames = list(gene_universe, NULL))
for (index in seq_along(rankings)) {
  positions <- seq_along(rankings[[index]])
  rows <- match(rankings[[index]], gene_universe)
  reciprocal_rank[rows, index] <- 1 / positions
  observed_rank[rows, index] <- positions
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
), ]
stopifnot(nrow(consensus) >= max(panel_sizes))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
source_manifest <- data.frame(
  source_ranking = ranking_files,
  md5 = unname(tools::md5sum(ranking_files)),
  stringsAsFactors = FALSE
)
write.csv(source_manifest, source_manifest_file, row.names = FALSE)

locked_at <- format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
aggregation_rule <- paste(
  "mean reciprocal rank across 15 source outer-split rankings;",
  "absent gene score 0; ties by frequency, mean observed rank, gene"
)
locked <- do.call(rbind, lapply(panel_sizes, function(panel_size) {
  selected <- consensus[seq_len(panel_size), ]
  data.frame(
    source_accession = "GSE107994",
    method = primary,
    panel_size = panel_size,
    rank = seq_len(panel_size),
    gene = selected$gene,
    mean_reciprocal_rank = selected$mean_reciprocal_rank,
    appearance_frequency = selected$appearance_frequency,
    mean_observed_rank = selected$mean_observed_rank,
    source_rankings_manifest_md5 = unname(tools::md5sum(source_manifest_file)),
    locked_at_utc = locked_at,
    aggregation_rule = aggregation_rule,
    independent_test_eligibility = paste(
      "future external cohort only; GSE19442 was accessed before this lock"
    ),
    stringsAsFactors = FALSE
  )
}))
write.csv(locked, output_file, row.names = FALSE)
cat(sprintf("Locked %s at panel sizes %s to %s\n", primary,
            paste(panel_sizes, collapse = ", "), output_file))
