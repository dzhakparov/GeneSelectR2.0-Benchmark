# Build the score-filtered STRING edge cache used when the STRING version API
# is unavailable. The source is the frozen STRING v12 combined-score file.

command_args <- commandArgs(trailingOnly = TRUE)
score_threshold <- if (length(command_args) >= 1L) {
  as.integer(command_args[1])
} else {
  400L
}
if (length(score_threshold) != 1L || !is.finite(score_threshold) ||
    score_threshold < 1L) {
  stop("The STRING score threshold must be one positive integer.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Building the offline STRING cache requires data.table.")
}

string_version <- "12.0"
organism <- 9606L
cache_dir <- file.path("data", "string_db_cache")
source_file <- file.path(
  cache_dir, sprintf("%s.protein.links.v%s.txt.gz", organism, string_version)
)
output_file <- file.path(
  cache_dir,
  sprintf("%s.protein.links.score%d.v%s.rds", organism, score_threshold,
          string_version)
)
if (!file.exists(source_file)) {
  stop("STRING links source is absent: ", source_file)
}

# awk filters the decompressed stream before it enters R. Loading all 13.7
# million source rows before applying the threshold requires several gigabytes
# and provides no additional information for the configured graph.
filter_command <- sprintf(
  "gzip -dc %s | awk 'NR == 1 || $3 >= %d'",
  shQuote(source_file), score_threshold
)
edges <- data.table::fread(cmd = filter_command, showProgress = TRUE)
required_columns <- c("protein1", "protein2", "combined_score")
if (!all(required_columns %in% colnames(edges)) || nrow(edges) == 0L ||
    any(!is.finite(edges$combined_score)) ||
    any(edges$combined_score < score_threshold)) {
  stop("Filtered STRING edge table failed validation.")
}

temporary_file <- tempfile(pattern = basename(output_file),
                           tmpdir = dirname(output_file))
on.exit(unlink(temporary_file), add = TRUE)
saveRDS(as.data.frame(edges[, ..required_columns]), temporary_file,
        compress = "xz")
if (!file.rename(temporary_file, output_file)) {
  stop("Could not install the offline STRING edge cache atomically.")
}

cat(sprintf("Wrote %s: %d edges at combined score >= %d\n",
            output_file, nrow(edges), score_threshold))
