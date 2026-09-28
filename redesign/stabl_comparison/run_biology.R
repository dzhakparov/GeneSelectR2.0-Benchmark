#!/usr/bin/env Rscript

# Score saved Stabl rankings with the frozen seven-dataset biology references.
args <- commandArgs(trailingOnly = TRUE)
datasets <- c("GSE101794", "GSE107994", "GSE13355", "GSE65682",
              "GSE69683", "imvigor210", "sosall")
selected <- if (length(args)) args else datasets
stopifnot(all(selected %in% datasets))

suppressPackageStartupMessages({
  library(igraph)
  library(withr)
})
repo <- normalizePath(getwd(), mustWork = TRUE)
out <- Sys.getenv("STABL_BIOLOGY_OUT",
                  file.path(repo, "redesign/stabl_comparison/biology"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ref_dir <- file.path(out, "reference_inputs")
commit <- "a840a6bf6ab36c87ac69ed65a5d506b6f29de21d"
sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
n_null <- 1000L

read_frozen <- function(path, reader) {
  dest <- tempfile(fileext = if (identical(reader, readRDS)) ".rds" else ".csv")
  on.exit(unlink(dest))
  status <- system2("git", c("-C", shQuote(repo), "show",
                             shQuote(paste0(commit, ":", path))), stdout = dest)
  if (!identical(status, 0L)) stop("Missing frozen input: ", path)
  reader(dest)
}
ratio <- function(observed, expected) {
  if (expected > 0) observed / expected else NA_real_
}
edge_count <- function(graph, vertices) {
  vertices <- unique(intersect(vertices, igraph::V(graph)$name))
  if (length(vertices) < 2L) return(0L)
  igraph::ecount(igraph::induced_subgraph(graph, vertices))
}
add_missing <- function(graph, vertices) {
  missing <- setdiff(vertices, igraph::V(graph)$name)
  if (length(missing)) igraph::add_vertices(graph, length(missing), name = missing)
  else graph
}
build_hallmark_graph <- function(pool, sets) {
  edge_tables <- lapply(sets, function(members) {
    present <- intersect(unique(members), pool)
    if (length(present) < 2L) return(NULL)
    pairs <- utils::combn(sort(present), 2L)
    data.frame(from = pairs[1L, ], to = pairs[2L, ])
  })
  edge_tables <- Filter(Negate(is.null), edge_tables)
  graph <- if (length(edge_tables)) {
    igraph::graph_from_data_frame(unique(do.call(rbind, edge_tables)),
                                  directed = FALSE)
  } else igraph::make_empty_graph(directed = FALSE)
  add_missing(graph, pool)
}
build_string_graph <- function(pool_ids, edges) {
  in_pool <- edges$protein1 %in% pool_ids & edges$protein2 %in% pool_ids
  graph <- igraph::graph_from_data_frame(
    edges[in_pool, c("protein1", "protein2"), drop = FALSE], directed = FALSE)
  add_missing(graph, pool_ids)
}
one_expected <- function(values, label) {
  value <- unique(values)
  if (length(value) != 1L || !is.finite(value)) stop("Inconsistent ", label)
  value
}

source(file.path(repo, "redesign", "R", "bio_prior.R"))
hallmark_file <- file.path(ref_dir, "hallmark_sets.rds")
hallmark_sets <- if (file.exists(hallmark_file)) {
  readRDS(hallmark_file)
} else hallmark_sets_all()
hallmark_genes <- unique(unlist(hallmark_sets, use.names = FALSE))
string_info_path <- file.path(ref_dir, "9606.protein.info.v12.0.txt.gz")
string_edges_path <- file.path(ref_dir, "9606.protein.links.score400.v12.0.rds")
if (!file.exists(string_info_path))
  string_info_path <- file.path(repo, "data/string_db_cache/9606.protein.info.v12.0.txt.gz")
if (!file.exists(string_edges_path))
  string_edges_path <- file.path(repo, "data/string_db_cache/9606.protein.links.score400.v12.0.rds")
string_info <- read.delim(gzfile(string_info_path), stringsAsFactors = FALSE)
symbol_to_string <- stats::setNames(string_info$X.string_protein_id,
                                    string_info$preferred_name)
string_edges <- readRDS(string_edges_path)
ot_baseline_file <- file.path(ref_dir, "ot_dense_by_split.csv")
ot_baseline <- if (file.exists(ot_baseline_file)) {
  read.csv(ot_baseline_file)
} else read_frozen(
                 "redesign/results_corrected/older7_biology_ot_dense_2026-09-03/ot_dense_by_split.csv",
                 read.csv)
pool_archive_file <- file.path(ref_dir, "candidate_pools.rds")
pool_archive <- if (file.exists(pool_archive_file)) {
  readRDS(pool_archive_file)
} else NULL

for (dataset in selected) {
  cat("Scoring ", dataset, "\n", sep = "")
  subdir <- if (dataset %in% c("imvigor210", "sosall"))
    "full_recipe" else "validation_benchmark"
  frozen_dir <- paste("redesign/results_corrected", subdir, dataset, sep = "/")
  baseline_file <- file.path(ref_dir, paste0(dataset, "_biology.csv"))
  semantic_file <- file.path(ref_dir, paste0(dataset, "_semantic.rds"))
  baseline <- if (file.exists(baseline_file)) {
    read.csv(baseline_file)
  } else read_frozen(paste0(frozen_dir, "/older7_biology_multiaxis.csv"),
                               read.csv)
  semantic <- if (file.exists(semantic_file)) {
    readRDS(semantic_file)
  } else read_frozen(paste0(frozen_dir, "/older7_biology_semantic_reference.rds"),
                               readRDS)
  source(file.path(repo, "analysis", "config.R"), local = TRUE)
  config <- get_biology_config("older7")
  ontology_id <- config$ontology_ids[config$dataset == dataset]
  stopifnot(length(ontology_id) == 1L)
  ot_file <- file.path(repo, "data/r_user_cache/R/GeneSelectR",
    sprintf("ot_seeds_%s_n3000_s0.rds", gsub("[^A-Za-z0-9]", "", ontology_id)))
  archived_ot <- file.path(ref_dir, basename(ot_file))
  if (file.exists(archived_ot)) ot_file <- archived_ot
  ot_table <- readRDS(ot_file)
  stopifnot(!anyDuplicated(ot_table$symbol), all(is.finite(ot_table$score)))
  association <- stats::setNames(ot_table$score, ot_table$symbol)
  rows <- vector("list", 15L * length(sizes))
  row_index <- 0L

  for (repeat_idx in 1:3) for (fold_idx in 1:5) {
    pool_key <- sprintf("%s_r%d_f%d", dataset, repeat_idx, fold_idx)
    if (!is.null(pool_archive)) {
      pool <- pool_archive[[pool_key]]
    } else {
      split <- read_frozen(sprintf(
        "redesign/results_corrected/%s/%s/split_r%d_f%d.rds",
        if (subdir == "full_recipe") "grouped_benchmark" else subdir,
        dataset, repeat_idx, fold_idx), readRDS)
      pool <- split$pools$var2000
    }
    stopifnot(length(pool) == 2000L, !anyDuplicated(pool))
    ranking_file <- file.path(repo, "redesign/stabl_comparison", sprintf(
      "%s_r%d_f%d_Stabl.csv", dataset, repeat_idx, fold_idx))
    ranking <- read.csv(ranking_file)$gene
    stopifnot(length(ranking) == 2000L, !anyDuplicated(ranking),
              setequal(ranking, pool))
    base_split <- baseline[baseline$repeat_idx == repeat_idx &
                           baseline$fold_idx == fold_idx, , drop = FALSE]
    ot_split <- ot_baseline[ot_baseline$dataset == dataset &
                            ot_baseline$repeat_idx == repeat_idx &
                            ot_baseline$fold_idx == fold_idx &
                            ot_baseline$method == "GS_full_ungrouped", , drop = FALSE]
    stopifnot(nrow(base_split) == 7L * length(sizes),
              nrow(ot_split) == 2L * length(sizes))
    semantic_pool <- semantic$semantic_evidence_ratio[pool]
    depth_pool <- semantic$annotation_depth[pool]
    stopifnot(all(is.finite(semantic_pool)), !anyNA(depth_pool))
    association_pool <- association[pool]
    names(association_pool) <- pool
    association_pool[is.na(association_pool)] <- 0
    hallmark_graph <- build_hallmark_graph(pool, hallmark_sets)
    pool_string_ids <- unique(stats::na.omit(symbol_to_string[pool]))
    string_graph <- build_string_graph(pool_string_ids, string_edges)

    for (panel_size in sizes) {
      panel <- head(ranking, panel_size)
      base_size <- base_split[base_split$k == panel_size, , drop = FALSE]
      ot_size <- ot_split[ot_split$k == panel_size, , drop = FALSE]
      stopifnot(nrow(base_size) == 7L, nrow(ot_size) == 2L)
      null_seed <- 940000L + repeat_idx * 10000L + fold_idx * 100L +
        match(panel_size, sizes)
      stopifnot(all(base_size$null_seed == null_seed),
                all(ot_size$null_seed == null_seed),
                all(base_size$n_null == n_null), all(ot_size$n_null == n_null))

      go_observed <- mean(semantic_pool[panel])
      go_expected <- one_expected(base_size$GO_semantic_expected, "GO reference")
      hallmark_observed <- edge_count(hallmark_graph, panel)
      hallmark_expected <- one_expected(base_size$hallmark_expected,
                                         "Hallmark reference")
      panel_string_ids <- unique(stats::na.omit(symbol_to_string[panel]))
      panel_string_ids <- intersect(panel_string_ids, igraph::V(string_graph)$name)
      mapped <- length(panel_string_ids)
      string_observed <- edge_count(string_graph, panel_string_ids)
      matching <- base_size[base_size$string_mapped == mapped, , drop = FALSE]
      if (nrow(matching)) {
        string_expected <- one_expected(matching$string_expected,
                                        "STRING reference")
        string_reference <- "saved matched null"
      } else {
        null_string <- if (mapped < 2L) rep(0, n_null) else {
          withr::with_seed(null_seed + mapped * 1000L,
            vapply(seq_len(n_null), function(index) {
              edge_count(string_graph, sample(pool_string_ids, mapped))
            }, numeric(1)))
        }
        string_expected <- mean(null_string)
        string_reference <- "recomputed matched null"
      }
      ot_values <- lapply(c(0.05, 0.10), function(cutoff) {
        scores <- association_pool
        scores[scores < cutoff] <- 0
        reference <- ot_size[ot_size$cutoff == cutoff, , drop = FALSE]
        stopifnot(nrow(reference) == 1L)
        expected <- reference$open_targets_expected
        observed <- sum(scores[panel])
        c(observed = observed, expected = expected,
          ratio = ratio(observed, expected), annotated = sum(scores[panel] > 0))
      })
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        dataset, repeat_idx, fold_idx, method = "Stabl", k = panel_size,
        GO_observed = go_observed, GO_expected = go_expected,
        GO_enrichment = ratio(go_observed, go_expected),
        GO_annotated = sum(depth_pool[panel] > 0),
        hallmark_edges = hallmark_observed, hallmark_expected,
        hallmark_enrichment = ratio(hallmark_observed, hallmark_expected),
        hallmark_annotated = sum(panel %in% hallmark_genes),
        string_edges = string_observed, string_expected,
        string_enrichment = ratio(string_observed, string_expected),
        string_mapped = mapped, string_reference,
        ot_0.05_observed = ot_values[[1L]]["observed"],
        ot_0.05_expected = ot_values[[1L]]["expected"],
        ot_0.05_enrichment = ot_values[[1L]]["ratio"],
        ot_0.05_annotated = ot_values[[1L]]["annotated"],
        ot_0.10_observed = ot_values[[2L]]["observed"],
        ot_0.10_expected = ot_values[[2L]]["expected"],
        ot_0.10_enrichment = ot_values[[2L]]["ratio"],
        ot_0.10_annotated = ot_values[[2L]]["annotated"],
        pool_size = length(pool), n_null, null_seed)
    }
    cat(sprintf("  r%d f%d complete\n", repeat_idx, fold_idx))
    flush.console()
  }
  result <- do.call(rbind, rows)
  keys <- with(result, paste(repeat_idx, fold_idx, k))
  stopifnot(nrow(result) == 15L * length(sizes), !anyDuplicated(keys),
            all(result$pool_size == 2000L), all(result$n_null == 1000L))
  write.csv(result, file.path(out, paste0(dataset, "_stabl_biology.csv")),
            row.names = FALSE)
  cat("Saved ", dataset, "\n", sep = "")
}
