#!/usr/bin/env Rscript

# Four-axis post-selection biology assessment for the seven older benchmark
# datasets. Panels are taken from saved outer-fold rankings. Every panel is
# compared with 1,000 same-size panels from its split-specific candidate pool.

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1L) args[[1L]] else stop("dataset required")
stage <- if (length(args) >= 2L) args[[2L]] else "all"
budget <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 604800
validation_datasets_old <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
)
datasets <- c(validation_datasets_old, "imvigor210", "sosall")
stopifnot(dataset %in% datasets,
          stage %in% c("init", "prepare", "score", "assemble", "all"))
started_at <- proc.time()[["elapsed"]]
over_budget <- function() proc.time()[["elapsed"]] - started_at > budget

suppressPackageStartupMessages({
  library(igraph)
  library(withr)
})
options(warn = 1)
for (path in list.files("package/GeneSelectR/R", full.names = TRUE)) {
  source(path)
}
source(file.path("redesign", "R", "bio_prior.R"))
source(file.path("redesign", "R", "run_provenance.R"))

results_root <- redesign_results_root()
is_validation <- dataset %in% validation_datasets_old
split_dir <- if (is_validation) {
  file.path(results_root, "validation_benchmark", dataset)
} else {
  file.path(results_root, "grouped_benchmark", dataset)
}
dataset_dir <- if (is_validation) split_dir else {
  file.path(results_root, "full_recipe", dataset)
}
require_redesign_run(split_dir)
if (!identical(split_dir, dataset_dir)) require_redesign_run(dataset_dir)

methods <- c(
  "GS_full_ungrouped", "DGE", "LASSO", "ElasticNet", "mRMR", "Boruta",
  "RF_importance"
)
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
n_null <- 1000L
source(file.path("analysis", "config.R"))
config <- get_biology_config("older7")
config_row <- config[config$dataset == dataset, , drop = FALSE]
if (nrow(config_row) != 1L) {
  stop("Expected one biology configuration for ", dataset, call. = FALSE)
}
association_files <- strsplit(
  config_row$association_files[[1L]], ";", fixed = TRUE
)[[1L]]
target_terms <- strsplit(
  config_row$target_go_terms[[1L]], ";", fixed = TRUE
)[[1L]]
if (!all(file.exists(association_files))) {
  stop("Missing frozen Open Targets input for ", dataset, call. = FALSE)
}

string_info_path <- file.path(
  "data", "string_db_cache", "9606.protein.info.v12.0.txt.gz"
)
string_edges_path <- file.path(
  "data", "string_db_cache", "9606.protein.links.score400.v12.0.rds"
)
semantic_reference_path <- file.path(
  dataset_dir, "older7_biology_semantic_reference.rds"
)
biology_output <- file.path(dataset_dir, "older7_biology_multiaxis.csv")

ranking_path <- function(repeat_idx, fold_idx, method) {
  saved_method <- if (!is_validation && method == "GS_full_ungrouped") {
    "full_ungrouped"
  } else {
    method
  }
  file.path(dataset_dir, sprintf(
    "ranking_r%d_f%d_%s.csv", repeat_idx, fold_idx, saved_method
  ))
}
ranking_files <- unlist(lapply(1:3, function(repeat_idx) {
  unlist(lapply(1:5, function(fold_idx) {
    vapply(methods, function(method) {
      ranking_path(repeat_idx, fold_idx, method)
    }, character(1))
  }), use.names = FALSE)
}), use.names = FALSE)
split_files <- unlist(lapply(1:3, function(repeat_idx) {
  vapply(1:5, function(fold_idx) file.path(
    split_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
  ), character(1))
}), use.names = FALSE)
if (!all(file.exists(c(ranking_files, split_files)))) {
  stop("One or more saved rankings or splits are missing.", call. = FALSE)
}

extension <- "older7_biology_multiaxis_v1"
prepare_redesign_extension(
  dataset_dir,
  extension = extension,
  source_files = c(
    "redesign/run_older7_biology.R", config_path,
    list.files("package/GeneSelectR/R", full.names = TRUE),
    file.path("redesign", "R", c("bio_prior.R", "run_provenance.R")),
    association_files, string_info_path, string_edges_path,
    file.path(split_dir, "run_manifest.rds"), ranking_files, split_files
  ),
  config = list(
    dataset = dataset,
    methods = methods,
    panel_sizes = panel_sizes,
    n_null = n_null,
    string_version = "12.0",
    string_score_threshold = 400L,
    open_targets_ontology_ids = config_row$ontology_ids[[1L]],
    GO_terms = target_terms,
    GO_ontology = "BP",
    GO_similarity = "resnik",
    GO_calibration = "annotation_depth_matched",
    hallmark_collection = "H"
  ),
  output_files = c(semantic_reference_path, biology_output)
)
if (stage == "init") {
  cat(sprintf("%s older-seven biology extension initialized\n", dataset))
  quit(save = "no")
}

load_association_scores <- function(paths) {
  tables <- lapply(paths, function(path) {
    table <- readRDS(path)
    if (!all(c("symbol", "score") %in% names(table)) ||
        anyDuplicated(table$symbol) || any(!is.finite(table$score))) {
      stop("Invalid frozen Open Targets file: ", path, call. = FALSE)
    }
    stats::setNames(table$score, table$symbol)
  })
  symbols <- unique(unlist(lapply(tables, names), use.names = FALSE))
  score_matrix <- vapply(tables, function(scores) {
    output <- scores[symbols]
    output[is.na(output)] <- 0
    output
  }, numeric(length(symbols)))
  output <- if (is.null(dim(score_matrix))) score_matrix else {
    apply(score_matrix, 1L, max)
  }
  stats::setNames(as.numeric(output), symbols)
}

if (stage %in% c("prepare", "all") &&
    !file.exists(semantic_reference_path)) {
  pool_union <- unique(unlist(lapply(split_files, function(path) {
    readRDS(path)$pools$var2000
  }), use.names = FALSE))
  semantic_percentile <- score_semantic_layer(
    pool_union,
    target_terms = target_terms,
    ontology = "BP",
    sim_method = "resnik",
    use_cache = TRUE,
    verbose = TRUE
  )
  annotation_depth <- get_annotation_depth(
    pool_union, ontology = "BP", organism = "human"
  )
  semantic_ratio <- calibrate_by_depth(
    semantic_percentile, annotation_depth,
    n_bins = 10, min_bin_size = 20, epsilon = 0.01,
    winsorize_at = 4, verbose = FALSE
  )
  names(semantic_ratio) <- pool_union
  saveRDS(list(
    dataset = dataset,
    disease_label = config_row$disease_label[[1L]],
    target_go_terms = target_terms,
    pool_union = pool_union,
    semantic_percentile = semantic_percentile,
    annotation_depth = annotation_depth,
    semantic_evidence_ratio = semantic_ratio
  ), semantic_reference_path)
  cat(sprintf(
    "%s semantic reference prepared for %d genes\n",
    dataset, length(pool_union)
  ))
}

add_missing_vertices <- function(graph, vertices) {
  missing <- setdiff(vertices, igraph::V(graph)$name)
  if (length(missing) > 0L) {
    graph <- igraph::add_vertices(graph, length(missing), name = missing)
  }
  graph
}
edge_count <- function(graph, vertices) {
  vertices <- unique(vertices[vertices %in% igraph::V(graph)$name])
  if (length(vertices) < 2L) return(0)
  igraph::ecount(igraph::induced_subgraph(graph, vids = vertices))
}
empirical_summary <- function(observed, null_values) {
  if (length(null_values) != n_null || any(!is.finite(null_values))) {
    stop("Biology null contains invalid values.", call. = FALSE)
  }
  expected <- mean(null_values)
  null_sd <- stats::sd(null_values)
  c(
    expected = expected,
    ratio = if (expected > 0) observed / expected else NA_real_,
    z = if (is.finite(null_sd) && null_sd > 0) {
      (observed - expected) / null_sd
    } else {
      NA_real_
    },
    p = (1 + sum(null_values >= observed)) / (length(null_values) + 1)
  )
}

if (stage %in% c("score", "all")) {
  if (!file.exists(semantic_reference_path)) {
    stop("Run the biology prepare stage first.", call. = FALSE)
  }
  semantic_reference <- readRDS(semantic_reference_path)
  association_scores <- load_association_scores(association_files)
  hallmark_sets <- hallmark_sets_all()
  string_info <- read.delim(gzfile(string_info_path),
                            stringsAsFactors = FALSE)
  symbol_to_string <- stats::setNames(
    string_info$X.string_protein_id, string_info$preferred_name
  )
  string_edges <- readRDS(string_edges_path)

  build_string_graph <- function(pool) {
    ids <- unique(stats::na.omit(symbol_to_string[pool]))
    edge_table <- string_edges[
      string_edges$protein1 %in% ids & string_edges$protein2 %in% ids,
      c("protein1", "protein2"), drop = FALSE
    ]
    graph <- igraph::graph_from_data_frame(edge_table, directed = FALSE)
    add_missing_vertices(graph, ids)
  }
  build_hallmark_graph <- function(pool) {
    edge_tables <- lapply(hallmark_sets, function(members) {
      present <- intersect(unique(members), pool)
      if (length(present) < 2L) return(NULL)
      combinations <- utils::combn(sort(present), 2L)
      data.frame(from = combinations[1L, ], to = combinations[2L, ],
                 stringsAsFactors = FALSE)
    })
    edge_tables <- edge_tables[!vapply(edge_tables, is.null, logical(1))]
    graph <- if (length(edge_tables) == 0L) {
      igraph::make_empty_graph(directed = FALSE)
    } else {
      igraph::graph_from_data_frame(
        unique(do.call(rbind, edge_tables)), directed = FALSE
      )
    }
    add_missing_vertices(graph, pool)
  }

  for (repeat_idx in 1:3) {
    for (fold_idx in 1:5) {
      output_path <- file.path(dataset_dir, sprintf(
        "older7_biology_multiaxis_r%d_f%d.rds", repeat_idx, fold_idx
      ))
      if (file.exists(output_path)) next
      if (over_budget()) {
        cat("[budget] stop in biology scoring\n")
        quit(save = "no")
      }
      split <- readRDS(file.path(
        split_dir, sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
      ))
      pool <- split$pools$var2000
      rankings <- lapply(methods, function(method) {
        path <- ranking_path(repeat_idx, fold_idx, method)
        ranking <- read.csv(path, stringsAsFactors = FALSE)$gene
        if (length(ranking) != length(pool) || anyDuplicated(ranking) ||
            !setequal(ranking, pool)) {
          stop("Invalid ranking: ", path, call. = FALSE)
        }
        ranking
      })
      names(rankings) <- methods

      open_targets_pool <- stats::setNames(association_scores[pool], pool)
      open_targets_pool[is.na(open_targets_pool)] <- 0
      semantic_pool <- semantic_reference$semantic_evidence_ratio[pool]
      annotation_depth <- semantic_reference$annotation_depth[pool]
      stopifnot(all(is.finite(semantic_pool)), !anyNA(annotation_depth))
      string_graph <- build_string_graph(pool)
      hallmark_graph <- build_hallmark_graph(pool)
      string_ids_pool <- unique(stats::na.omit(symbol_to_string[pool]))
      hallmark_annotated <- unique(intersect(
        pool, unlist(hallmark_sets, use.names = FALSE)
      ))

      rows <- list()
      for (panel_size in panel_sizes) {
        null_seed <- 940000L + repeat_idx * 10000L + fold_idx * 100L +
          match(panel_size, panel_sizes)
        null_panels <- withr::with_seed(null_seed, lapply(
          seq_len(n_null), function(index) sample(pool, panel_size)
        ))
        null_open_targets <- vapply(null_panels, function(panel) {
          sum(open_targets_pool[panel])
        }, numeric(1))
        null_semantic <- vapply(null_panels, function(panel) {
          mean(semantic_pool[panel])
        }, numeric(1))
        null_hallmark <- vapply(null_panels, function(panel) {
          edge_count(hallmark_graph, panel)
        }, numeric(1))
        string_null_cache <- new.env(parent = emptyenv())

        for (method in methods) {
          panel <- head(rankings[[method]], panel_size)
          panel_string_ids <- unique(stats::na.omit(symbol_to_string[panel]))
          panel_string_ids <- intersect(
            panel_string_ids, igraph::V(string_graph)$name
          )
          n_string_mapped <- length(panel_string_ids)
          string_key <- as.character(n_string_mapped)
          if (!exists(string_key, envir = string_null_cache,
                      inherits = FALSE)) {
            null_string <- if (n_string_mapped < 2L) {
              rep(0, n_null)
            } else {
              withr::with_seed(
                null_seed + n_string_mapped * 1000L,
                vapply(seq_len(n_null), function(index) {
                  edge_count(
                    string_graph,
                    sample(string_ids_pool, n_string_mapped)
                  )
                }, numeric(1))
              )
            }
            assign(string_key, null_string, envir = string_null_cache)
          }
          null_string <- get(
            string_key, envir = string_null_cache, inherits = FALSE
          )
          summaries <- list(
            open_targets = empirical_summary(
              sum(open_targets_pool[panel]), null_open_targets
            ),
            semantic = empirical_summary(
              mean(semantic_pool[panel]), null_semantic
            ),
            hallmark = empirical_summary(
              edge_count(hallmark_graph, panel), null_hallmark
            ),
            string = empirical_summary(
              edge_count(string_graph, panel_string_ids), null_string
            )
          )
          rows[[length(rows) + 1L]] <- data.frame(
            repeat_idx = repeat_idx,
            fold_idx = fold_idx,
            method = method,
            k = panel_size,
            disease_label = config_row$disease_label[[1L]],
            ontology_ids = config_row$ontology_ids[[1L]],
            open_targets_sum = sum(open_targets_pool[panel]),
            open_targets_expected = summaries$open_targets["expected"],
            open_targets_enrichment = summaries$open_targets["ratio"],
            open_targets_z = summaries$open_targets["z"],
            open_targets_empirical_p = summaries$open_targets["p"],
            open_targets_overlap = sum(open_targets_pool[panel] > 0),
            open_targets_pool_targets = sum(open_targets_pool > 0),
            GO_semantic_mean_ratio = mean(semantic_pool[panel]),
            GO_semantic_expected = summaries$semantic["expected"],
            GO_semantic_enrichment = summaries$semantic["ratio"],
            GO_semantic_z = summaries$semantic["z"],
            GO_semantic_empirical_p = summaries$semantic["p"],
            GO_annotated = sum(annotation_depth[panel] > 0),
            hallmark_edges = edge_count(hallmark_graph, panel),
            hallmark_expected = summaries$hallmark["expected"],
            hallmark_enrichment = summaries$hallmark["ratio"],
            hallmark_z = summaries$hallmark["z"],
            hallmark_empirical_p = summaries$hallmark["p"],
            hallmark_annotated = sum(panel %in% hallmark_annotated),
            string_edges = edge_count(string_graph, panel_string_ids),
            string_expected = summaries$string["expected"],
            string_enrichment = summaries$string["ratio"],
            string_z = summaries$string["z"],
            string_empirical_p = summaries$string["p"],
            string_mapped = n_string_mapped,
            pool_size = length(pool),
            n_null = n_null,
            null_seed = null_seed,
            stringsAsFactors = FALSE
          )
        }
      }
      result <- do.call(rbind, rows)
      stopifnot(nrow(result) == length(methods) * length(panel_sizes))
      saveRDS(result, output_path)
      cat(sprintf(
        "[%s r%d f%d] multi-axis biology complete\n",
        dataset, repeat_idx, fold_idx
      ))
    }
  }
}

if (stage %in% c("assemble", "all")) {
  paths <- unlist(lapply(1:3, function(repeat_idx) {
    vapply(1:5, function(fold_idx) file.path(dataset_dir, sprintf(
      "older7_biology_multiaxis_r%d_f%d.rds", repeat_idx, fold_idx
    )), character(1))
  }), use.names = FALSE)
  if (!all(file.exists(paths))) {
    stop("Missing one or more split-level biology results.", call. = FALSE)
  }
  output <- do.call(rbind, lapply(paths, readRDS))
  keys <- with(output, paste(repeat_idx, fold_idx, method, k))
  stopifnot(
    nrow(output) == 15L * length(methods) * length(panel_sizes),
    !anyDuplicated(keys),
    all(output$pool_size == 2000L),
    all(output$n_null == n_null),
    all(output$open_targets_overlap <= output$k),
    all(output$GO_annotated <= output$k),
    all(output$hallmark_annotated <= output$k),
    all(output$string_mapped <= output$k)
  )
  write.csv(output, biology_output, row.names = FALSE)
  cat(sprintf("%s older-seven biology results assembled\n", dataset))
}

cat(sprintf("%s biology stage=%s complete\n", dataset, stage))
