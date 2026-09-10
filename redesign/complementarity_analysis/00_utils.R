# Shared functions for the GeneSelectR and differential-expression
# complementarity analyses.
#
# Run all scripts from the repository root. The functions below only read
# existing benchmark files and write to redesign/complementarity_analysis/results.

complementarity_root <- file.path("redesign", "complementarity_analysis")
complementarity_results_dir <- Sys.getenv(
  "GENESELECTR_COMPLEMENTARITY_RESULTS_DIR",
  unset = file.path(complementarity_root, "results")
)
complementarity_figure_dir <- file.path(complementarity_results_dir, "figures")

analysis_methods <- c(
  "GeneSelectR", "DGE", "Random forest", "Boruta", "mRMR", "LASSO",
  "Elastic net"
)

primary_gene_set_sizes <- c(10L, 20L, 50L)
group_gene_set_sizes <- c(20L, 50L)

split_grid <- expand.grid(
  repeat_idx = 1:3,
  fold_idx = 1:5,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
split_grid$split_id <- sprintf(
  "r%d_f%d", split_grid$repeat_idx, split_grid$fold_idx
)

assert_repository_root <- function() {
  required <- c(
    file.path("package", "GeneSelectR", "DESCRIPTION"),
    file.path("redesign", "R", "evaluator.R"),
    file.path("redesign", "results_corrected")
  )
  missing <- required[!file.exists(required)]
  if (length(missing) > 0L) {
    stop(
      "Run this script from the GeneSelectR2.0-Benchmark repository root. ",
      "Missing: ", paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

require_packages <- function(packages) {
  missing <- packages[!vapply(
    packages, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing) > 0L) {
    stop(
      "Required R package(s) are unavailable: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

dataset_registry <- function() {
  corrected <- file.path("redesign", "results_corrected")
  additional <- file.path(
    "redesign", "results_frozen_external_exact_2026-08-31",
    "validation_benchmark"
  )

  data.frame(
    dataset = c(
      "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683",
      "imvigor210", "sosall", "GSE16879", "GSE91061", "GSE92415",
      "GSE206285"
    ),
    dataset_group = c(
      rep("validation", 6L), "development", rep("additional", 4L)
    ),
    ranking_dir = c(
      file.path(corrected, "validation_benchmark", c(
        "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
      )),
      file.path(corrected, "full_recipe", c("imvigor210", "sosall")),
      file.path(additional, c(
        "GSE16879", "GSE91061", "GSE92415", "GSE206285"
      ))
    ),
    split_dir = c(
      file.path(corrected, "validation_benchmark", c(
        "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
      )),
      file.path(corrected, "grouped_benchmark", c("imvigor210", "sosall")),
      file.path(additional, c(
        "GSE16879", "GSE91061", "GSE92415", "GSE206285"
      ))
    ),
    ranking_layout = c(rep("development", 7L), rep("additional", 4L)),
    stringsAsFactors = FALSE
  )
}

dataset_record <- function(dataset) {
  registry <- dataset_registry()
  row <- registry[registry$dataset == dataset, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Unknown dataset: ", dataset, call. = FALSE)
  }
  row
}

normalise_dataset_argument <- function(arguments) {
  available <- dataset_registry()$dataset
  if (length(arguments) == 0L || identical(arguments[[1L]], "all")) {
    return(available)
  }
  requested <- unique(arguments)
  unknown <- setdiff(requested, available)
  if (length(unknown) > 0L) {
    stop("Unknown dataset(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  requested
}

saved_method_name <- function(dataset, method) {
  record <- dataset_record(dataset)
  if (method == "GeneSelectR") {
    if (dataset %in% c("imvigor210", "sosall")) {
      return("full_ungrouped")
    }
    return("GS_full_ungrouped")
  }
  switch(
    method,
    "DGE" = "DGE",
    "Random forest" = "RF_importance",
    "Boruta" = "Boruta",
    "mRMR" = "mRMR",
    "LASSO" = "LASSO",
    "Elastic net" = "ElasticNet",
    stop("Unknown method: ", method, call. = FALSE)
  )
}

ranking_path <- function(dataset, repeat_idx, fold_idx, method) {
  record <- dataset_record(dataset)
  saved <- saved_method_name(dataset, method)
  prefix <- if (
    record$dataset_group[[1L]] == "additional" && method != "GeneSelectR"
  ) {
    "competitor"
  } else {
    "ranking"
  }
  file.path(
    record$ranking_dir[[1L]],
    sprintf("%s_r%d_f%d_%s.csv", prefix, repeat_idx, fold_idx, saved)
  )
}

split_path <- function(dataset, repeat_idx, fold_idx) {
  record <- dataset_record(dataset)
  file.path(
    record$split_dir[[1L]],
    sprintf("split_r%d_f%d.rds", repeat_idx, fold_idx)
  )
}

read_ranking <- function(dataset, repeat_idx, fold_idx, method) {
  path <- ranking_path(dataset, repeat_idx, fold_idx, method)
  if (!file.exists(path)) {
    stop("Missing saved ranking: ", path, call. = FALSE)
  }
  table <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!"gene" %in% names(table)) {
    stop("Saved ranking has no gene column: ", path, call. = FALSE)
  }
  table$gene <- as.character(table$gene)
  if (nrow(table) == 0L || anyNA(table$gene) || any(table$gene == "") ||
      anyDuplicated(table$gene)) {
    stop("Saved ranking contains invalid or duplicate genes: ", path,
         call. = FALSE)
  }
  attr(table, "source_path") <- path
  table
}

read_split <- function(dataset, repeat_idx, fold_idx) {
  path <- split_path(dataset, repeat_idx, fold_idx)
  if (!file.exists(path)) {
    stop("Missing saved train/test division: ", path, call. = FALSE)
  }
  split <- readRDS(path)
  needed <- c("train_idx", "test_idx", "pools")
  if (!all(needed %in% names(split)) ||
      !"var2000" %in% names(split$pools)) {
    stop("Saved division has an unexpected structure: ", path,
         call. = FALSE)
  }
  if (anyDuplicated(split$train_idx) || anyDuplicated(split$test_idx) ||
      length(intersect(split$train_idx, split$test_idx)) > 0L) {
    stop("Training and test samples are not disjoint: ", path,
         call. = FALSE)
  }
  pool <- as.character(split$pools$var2000)
  if (length(pool) != 2000L || anyNA(pool) || anyDuplicated(pool)) {
    stop("Expected 2,000 unique candidate genes in: ", path,
         call. = FALSE)
  }
  split$pools$var2000 <- pool
  attr(split, "source_path") <- path
  split
}

validate_split_rankings <- function(dataset, repeat_idx, fold_idx) {
  split <- read_split(dataset, repeat_idx, fold_idx)
  pool <- split$pools$var2000
  rankings <- lapply(analysis_methods, function(method) {
    table <- read_ranking(dataset, repeat_idx, fold_idx, method)
    if (nrow(table) != length(pool) || !setequal(table$gene, pool)) {
      stop(
        "Ranking and candidate genes differ for ", dataset, " r", repeat_idx,
        " f", fold_idx, " ", method, call. = FALSE
      )
    }
    table
  })
  names(rankings) <- analysis_methods
  list(split = split, pool = pool, rankings = rankings)
}

expected_ranking_paths <- function(dataset) {
  unlist(lapply(seq_len(nrow(split_grid)), function(index) {
    row <- split_grid[index, ]
    vapply(analysis_methods, function(method) {
      ranking_path(dataset, row$repeat_idx, row$fold_idx, method)
    }, character(1))
  }), use.names = FALSE)
}

validate_dataset_inputs <- function(dataset) {
  ranking_paths <- expected_ranking_paths(dataset)
  split_paths <- vapply(seq_len(nrow(split_grid)), function(index) {
    row <- split_grid[index, ]
    split_path(dataset, row$repeat_idx, row$fold_idx)
  }, character(1))
  missing <- c(ranking_paths[!file.exists(ranking_paths)],
               split_paths[!file.exists(split_paths)])
  if (length(missing) > 0L) {
    stop(
      "Missing ", length(missing), " required saved file(s) for ", dataset,
      ": ", paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

top_genes <- function(ranking, k, allowed_genes = NULL) {
  genes <- if (is.data.frame(ranking)) ranking$gene else as.character(ranking)
  if (!is.null(allowed_genes)) {
    genes <- genes[genes %in% allowed_genes]
  }
  if (length(genes) < k) {
    stop("Ranking contains fewer than ", k, " eligible genes.",
         call. = FALSE)
  }
  selected <- genes[seq_len(k)]
  if (length(selected) != k || anyDuplicated(selected)) {
    stop("Top-gene-set construction did not return k unique genes.",
         call. = FALSE)
  }
  selected
}

selection_matrix <- function(gene_sets, universe) {
  if (length(gene_sets) < 2L || length(universe) == 0L ||
      anyDuplicated(universe)) {
    stop("Invalid gene sets or gene universe for stability calculation.",
         call. = FALSE)
  }
  matrix <- vapply(gene_sets, function(genes) universe %in% genes,
                   logical(length(universe)))
  rownames(matrix) <- universe
  matrix
}

pairwise_set_statistics <- function(gene_sets) {
  if (length(gene_sets) < 2L) {
    return(c(
      mean_jaccard = NA_real_, median_jaccard = NA_real_,
      mean_dice = NA_real_, mean_overlap_coefficient = NA_real_
    ))
  }
  pairs <- utils::combn(seq_along(gene_sets), 2L)
  values <- apply(pairs, 2L, function(indices) {
    first <- unique(gene_sets[[indices[[1L]]]])
    second <- unique(gene_sets[[indices[[2L]]]])
    intersection_size <- length(intersect(first, second))
    union_size <- length(union(first, second))
    c(
      jaccard = intersection_size / union_size,
      dice = 2 * intersection_size / (length(first) + length(second)),
      overlap = intersection_size / min(length(first), length(second))
    )
  })
  c(
    mean_jaccard = mean(values["jaccard", ]),
    median_jaccard = stats::median(values["jaccard", ]),
    mean_dice = mean(values["dice", ]),
    mean_overlap_coefficient = mean(values["overlap", ])
  )
}

compute_repository_nogueira <- local({
  function_cache <- NULL
  function(selection_matrix_value) {
    if (is.null(function_cache)) {
      require_packages("stabm")
      source_path <- file.path("package", "GeneSelectR", "R", "utils.R")
      source_environment <- new.env(parent = globalenv())
      sys.source(source_path, envir = source_environment)
      function_cache <<- get(
        "compute_nogueira_stability", envir = source_environment,
        inherits = FALSE
      )
    }
    result <- function_cache(
      selection_matrix_value, rownames(selection_matrix_value)
    )
    unname(result$nogueira_index)
  }
})

canonical_method_from_saved <- function(value) {
  output <- as.character(value)
  output[output %in% c("GS_full_ungrouped", "full_ungrouped")] <-
    "GeneSelectR"
  output[output == "RF_importance"] <- "Random forest"
  output[output == "ElasticNet"] <- "Elastic net"
  output[output == "Random"] <- "Random"
  output
}

load_predictive_results <- function(dataset) {
  record <- dataset_record(dataset)
  sizes <- primary_gene_set_sizes

  if (record$dataset_group[[1L]] != "additional") {
    path <- file.path(record$ranking_dir[[1L]], "eval_deterministic.csv")
    if (!file.exists(path)) {
      stop("Missing deterministic predictive results: ", path,
           call. = FALSE)
    }
    table <- utils::read.csv(path, stringsAsFactors = FALSE)
    needed <- c("repeat_idx", "fold_idx", "arm", "k", "AUC")
    if (!all(needed %in% names(table))) {
      stop("Unexpected predictive-results schema: ", path, call. = FALSE)
    }
    table$method <- canonical_method_from_saved(table$arm)
    table <- table[
      table$method %in% c(analysis_methods, "Random") & table$k %in% sizes,
      c("repeat_idx", "fold_idx", "method", "k", "AUC"), drop = FALSE
    ]
  } else {
    competitor_path <- file.path(
      record$ranking_dir[[1L]], "competitor_eval_results.csv"
    )
    gs_path <- file.path(record$ranking_dir[[1L]], "eval_results.csv")
    if (!all(file.exists(c(competitor_path, gs_path)))) {
      stop(
        "Missing additional-dataset predictive results: ",
        paste(c(competitor_path, gs_path)[
          !file.exists(c(competitor_path, gs_path))
        ], collapse = ", "), call. = FALSE
      )
    }
    competitors <- utils::read.csv(
      competitor_path, stringsAsFactors = FALSE
    )
    names(competitors)[names(competitors) == "method"] <- "stored_method"
    competitors$method <- canonical_method_from_saved(
      competitors$stored_method
    )
    competitors <- competitors[
      competitors$method %in% analysis_methods & competitors$k %in% sizes,
      c("repeat_idx", "fold_idx", "method", "k", "AUC"), drop = FALSE
    ]

    gs <- utils::read.csv(gs_path, stringsAsFactors = FALSE)
    gs$method <- canonical_method_from_saved(gs$arm)
    gs <- gs[
      gs$method %in% c("GeneSelectR", "Random") & gs$k %in% sizes,
      c("repeat_idx", "fold_idx", "method", "k", "AUC"), drop = FALSE
    ]
    table <- rbind(competitors, gs)
  }

  key <- with(table, paste(repeat_idx, fold_idx, method, k, sep = ":"))
  if (anyDuplicated(key)) {
    stop("Duplicate predictive result rows for ", dataset, call. = FALSE)
  }
  expected_methods <- c(analysis_methods, "Random")
  expected_rows <- nrow(split_grid) * length(expected_methods) * length(sizes)
  if (nrow(table) != expected_rows ||
      !setequal(unique(table$method), expected_methods)) {
    stop(
      "Predictive results are incomplete for ", dataset, ". Expected ",
      expected_rows, " rows and found ", nrow(table), ".", call. = FALSE
    )
  }

  random <- table[table$method == "Random", ]
  observed <- table[table$method != "Random", ]
  random_key <- with(random, paste(repeat_idx, fold_idx, k, sep = ":"))
  match_index <- match(
    with(observed, paste(repeat_idx, fold_idx, k, sep = ":")), random_key
  )
  if (anyNA(match_index)) {
    stop("Could not match random-gene-set AUC for ", dataset,
         call. = FALSE)
  }
  observed$random_AUC <- random$AUC[match_index]
  observed$AUC_minus_random <- observed$AUC - observed$random_AUC
  observed$dataset <- dataset
  observed$dataset_group <- record$dataset_group[[1L]]
  observed[, c(
    "dataset", "dataset_group", "repeat_idx", "fold_idx", "method", "k",
    "AUC", "random_AUC", "AUC_minus_random"
  )]
}

load_all_predictive_results <- function(datasets = dataset_registry()$dataset) {
  do.call(rbind, lapply(datasets, load_predictive_results))
}

load_outcome <- function(dataset) {
  record <- dataset_record(dataset)
  if (dataset %in% c("imvigor210", "sosall")) {
    path <- file.path(record$ranking_dir[[1L]], "base_outcome.rds")
    if (!file.exists(path)) stop("Missing outcome file: ", path,
                                 call. = FALSE)
    outcome <- readRDS(path)
  } else {
    path <- file.path(record$split_dir[[1L]], "base_data.rds")
    if (!file.exists(path)) stop("Missing dataset file: ", path,
                                 call. = FALSE)
    outcome <- readRDS(path)$outcome
  }
  if (!is.factor(outcome) || nlevels(droplevels(outcome)) != 2L) {
    stop("Expected a two-level factor outcome for ", dataset,
         call. = FALSE)
  }
  outcome
}

high_is_best_rank_percentile <- function(rank_value, total) {
  if (total <= 1L) return(rep(1, length(rank_value)))
  1 - (rank_value - 1) / (total - 1)
}

safe_spearman <- function(first, second) {
  keep <- is.finite(first) & is.finite(second)
  if (sum(keep) < 3L || length(unique(first[keep])) < 2L ||
      length(unique(second[keep])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(first[keep], second[keep], method = "spearman"))
}

calculate_benchmark_dge_metrics <- function(train_features, train_labels) {
  train_labels <- droplevels(train_labels)
  if (nlevels(train_labels) != 2L ||
      nrow(train_features) != length(train_labels) ||
      is.null(colnames(train_features)) || anyDuplicated(colnames(train_features))) {
    stop("Invalid training data for the benchmark DGE calculation.",
         call. = FALSE)
  }

  n_genes <- ncol(train_features)
  scores <- numeric(n_genes)
  t_statistics <- numeric(n_genes)
  p_values <- numeric(n_genes)

  for (gene_index in seq_len(n_genes)) {
    test <- tryCatch(
      stats::t.test(train_features[, gene_index] ~ train_labels),
      error = function(error) NULL
    )
    if (is.null(test)) {
      scores[[gene_index]] <- 0
      t_statistics[[gene_index]] <- 0
      p_values[[gene_index]] <- 1
    } else {
      t_statistic <- unname(test$statistic)
      p_value <- unname(test$p.value)
      scores[[gene_index]] <- -log10(max(p_value, 1e-300)) *
        abs(t_statistic)
      t_statistics[[gene_index]] <- t_statistic
      p_values[[gene_index]] <- p_value
    }
  }

  genes <- colnames(train_features)
  adjusted_p_values <- stats::p.adjust(p_values, method = "BH")
  table <- data.frame(
    gene = genes,
    dge_ranking_statistic = scores,
    dge_t_statistic = t_statistics,
    dge_p_value = p_values,
    dge_fdr = adjusted_p_values,
    dge_selected_bh_0_05 = adjusted_p_values < 0.05,
    stringsAsFactors = FALSE
  )
  table <- table[order(table$dge_ranking_statistic, decreasing = TRUE), ]
  rownames(table) <- NULL
  table
}

atomic_write_csv <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(path), "."), tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(table, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    stop("Could not install output file: ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(path), "."), tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not install output file: ", path, call. = FALSE)
  }
  invisible(path)
}

biology_config_record <- function(dataset) {
  source(file.path("analysis", "config.R"), local = TRUE)
  config <- if (dataset_record(dataset)$dataset_group[[1L]] !=
                "additional") {
    get_biology_config("older7")
  } else {
    get_biology_config("external")
  }
  row <- config[config$dataset == dataset, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Expected one biological configuration for ", dataset,
         call. = FALSE)
  }
  attr(row, "source_path") <- config_path
  row
}

semantic_reference_path <- function(dataset) {
  record <- dataset_record(dataset)
  filename <- if (record$dataset_group[[1L]] != "additional") {
    "older7_biology_semantic_reference.rds"
  } else {
    "biology_semantic_reference_v2.rds"
  }
  file.path(record$ranking_dir[[1L]], filename)
}

load_semantic_reference <- function(dataset) {
  path <- semantic_reference_path(dataset)
  if (!file.exists(path)) {
    stop("Missing saved GO semantic reference: ", path, call. = FALSE)
  }
  reference <- readRDS(path)
  needed <- c(
    "pool_union", "semantic_evidence_ratio", "annotation_depth"
  )
  if (!all(needed %in% names(reference)) ||
      is.null(names(reference$semantic_evidence_ratio)) ||
      is.null(names(reference$annotation_depth))) {
    stop("Saved GO semantic reference has an unexpected structure: ", path,
         call. = FALSE)
  }
  reference
}

load_open_targets_scores <- function(dataset) {
  config <- biology_config_record(dataset)
  paths <- strsplit(
    config$association_files[[1L]], ";", fixed = TRUE
  )[[1L]]
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing Open Targets input: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  tables <- lapply(paths, function(path) {
    table <- readRDS(path)
    if (!all(c("symbol", "score") %in% names(table)) ||
        anyNA(table$symbol) || anyDuplicated(table$symbol) ||
        any(!is.finite(table$score))) {
      stop("Invalid Open Targets input: ", path, call. = FALSE)
    }
    stats::setNames(table$score, table$symbol)
  })
  symbols <- unique(unlist(lapply(tables, names), use.names = FALSE))
  score_matrix <- vapply(tables, function(scores) {
    output <- scores[symbols]
    output[is.na(output)] <- 0
    as.numeric(output)
  }, numeric(length(symbols)))
  combined <- if (is.null(dim(score_matrix))) {
    score_matrix
  } else {
    apply(score_matrix, 1L, max)
  }
  stats::setNames(as.numeric(combined), symbols)
}

load_hallmark_sets <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      require_packages("msigdbr")
      source_environment <- new.env(parent = globalenv())
      sys.source(
        file.path("redesign", "R", "bio_prior.R"),
        envir = source_environment
      )
      cache <<- get(
        "hallmark_sets_all", envir = source_environment, inherits = FALSE
      )()
    }
    cache
  }
})

add_missing_graph_vertices <- function(graph, vertices) {
  missing <- setdiff(vertices, igraph::V(graph)$name)
  if (length(missing) > 0L) {
    graph <- igraph::add_vertices(graph, length(missing), name = missing)
  }
  graph
}

build_hallmark_graph <- function(pool, hallmark_sets) {
  require_packages("igraph")
  edges <- lapply(hallmark_sets, function(members) {
    present <- intersect(unique(members), pool)
    if (length(present) < 2L) return(NULL)
    pairs <- utils::combn(sort(present), 2L)
    data.frame(
      from = pairs[1L, ], to = pairs[2L, ], stringsAsFactors = FALSE
    )
  })
  edges <- edges[!vapply(edges, is.null, logical(1))]
  graph <- if (length(edges) == 0L) {
    igraph::make_empty_graph(directed = FALSE)
  } else {
    igraph::graph_from_data_frame(
      unique(do.call(rbind, edges)), directed = FALSE
    )
  }
  add_missing_graph_vertices(graph, pool)
}

load_string_mapping_and_edges <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      info_path <- file.path(
        "data", "string_db_cache", "9606.protein.info.v12.0.txt.gz"
      )
      edge_path <- file.path(
        "data", "string_db_cache", "9606.protein.links.score400.v12.0.rds"
      )
      if (!all(file.exists(c(info_path, edge_path)))) {
        stop(
          "Missing frozen STRING v12 input: ",
          paste(c(info_path, edge_path)[!file.exists(c(info_path, edge_path))],
                collapse = ", "), call. = FALSE
        )
      }
      info <- utils::read.delim(
        gzfile(info_path), stringsAsFactors = FALSE
      )
      cache <<- list(
        symbol_to_string = stats::setNames(
          info$X.string_protein_id, info$preferred_name
        ),
        edges = readRDS(edge_path),
        info_path = info_path,
        edge_path = edge_path
      )
    }
    cache
  }
})

build_string_graph <- function(pool, string_resources) {
  require_packages("igraph")
  identifiers <- unique(stats::na.omit(
    string_resources$symbol_to_string[pool]
  ))
  edge_table <- string_resources$edges[
    string_resources$edges$protein1 %in% identifiers &
      string_resources$edges$protein2 %in% identifiers,
    c("protein1", "protein2"), drop = FALSE
  ]
  graph <- igraph::graph_from_data_frame(edge_table, directed = FALSE)
  add_missing_graph_vertices(graph, identifiers)
}

count_graph_edges <- function(graph, vertices) {
  require_packages("igraph")
  vertices <- unique(vertices[vertices %in% igraph::V(graph)$name])
  if (length(vertices) < 2L) return(NA_real_)
  igraph::ecount(igraph::induced_subgraph(graph, vids = vertices))
}

matched_random_summary <- function(observed, null_values) {
  if (!is.finite(observed) || length(null_values) == 0L ||
      any(!is.finite(null_values))) {
    return(list(
      mean = NA_real_, median = NA_real_, ratio = NA_real_,
      reason = "non_finite_input"
    ))
  }
  mean_value <- mean(null_values)
  median_value <- stats::median(null_values)
  ratio <- if (mean_value > 0) observed / mean_value else NA_real_
  reason <- if (mean_value > 0) "" else "matched_random_mean_zero"
  list(mean = mean_value, median = median_value, ratio = ratio,
       reason = reason)
}

assert_repository_root()
