#!/usr/bin/env Rscript

# Measure biological properties of genes shared by DGE and GeneSelectR and of
# genes selected by only one method. Each observed value is compared with
# 1,000 same-size random gene sets from the same split-specific candidate genes.

source(file.path("redesign", "complementarity_analysis", "00_utils.R"))
require_packages(c("igraph", "msigdbr"))

membership_path <- file.path(
  complementarity_results_dir, "gene_group_membership.csv"
)
if (!file.exists(membership_path)) {
  stop(
    "Missing input: ", membership_path,
    ". Run 03_gene_group_characterization.R first.", call. = FALSE
  )
}

datasets <- normalise_dataset_argument(commandArgs(trailingOnly = TRUE))
membership <- utils::read.csv(
  membership_path, stringsAsFactors = FALSE, check.names = FALSE
)
membership <- membership[membership$dataset %in% datasets, , drop = FALSE]

needed <- c(
  "dataset", "dataset_group", "repeat_idx", "fold_idx", "split_id", "k",
  "group", "gene"
)
missing <- setdiff(needed, names(membership))
if (length(missing) > 0L) {
  stop("Membership table is missing columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

n_random <- 1000L
groups <- c("SHARED", "DGE_ONLY", "GS_ONLY")
hallmark_sets <- load_hallmark_sets()
string_resources <- load_string_mapping_and_edges()
result_rows <- list()

add_metric_row <- function(dataset, dataset_group, repeat_idx, fold_idx, split_id,
                           k, group, group_size, metric, n_metric_genes,
                           observed, random_values, random_seed,
                           reason = "") {
  if (nzchar(reason)) {
    comparison <- list(
      mean = NA_real_, median = NA_real_, ratio = NA_real_, reason = reason
    )
  } else {
    comparison <- matched_random_summary(observed, random_values)
  }
  data.frame(
    dataset = dataset,
    dataset_group = dataset_group,
    repeat_idx = repeat_idx,
    fold_idx = fold_idx,
    split_id = split_id,
    k = k,
    group = group,
    group_size = group_size,
    metric = metric,
    n_metric_genes = n_metric_genes,
    observed_value = observed,
    matched_random_mean = comparison$mean,
    matched_random_median = comparison$median,
    observed_to_random_ratio = comparison$ratio,
    number_of_random_draws = length(random_values),
    random_seed = random_seed,
    missing_reason = comparison$reason,
    stringsAsFactors = FALSE
  )
}

for (dataset_index in seq_along(datasets)) {
  dataset <- datasets[[dataset_index]]
  record <- dataset_record(dataset)
  validate_dataset_inputs(dataset)
  open_targets <- load_open_targets_scores(dataset)
  semantic_reference <- load_semantic_reference(dataset)

  for (split_index in seq_len(nrow(split_grid))) {
    split_row <- split_grid[split_index, ]
    split <- read_split(dataset, split_row$repeat_idx, split_row$fold_idx)
    pool <- split$pools$var2000

    open_targets_pool <- stats::setNames(open_targets[pool], pool)
    open_targets_pool[is.na(open_targets_pool)] <- 0
    semantic_pool <- semantic_reference$semantic_evidence_ratio[pool]
    annotation_depth <- semantic_reference$annotation_depth[pool]
    if (anyNA(semantic_pool) || any(!is.finite(semantic_pool)) ||
        anyNA(annotation_depth)) {
      stop("Incomplete saved GO scores for ", dataset, " ",
           split_row$split_id, call. = FALSE)
    }

    hallmark_graph <- build_hallmark_graph(pool, hallmark_sets)
    string_graph <- build_string_graph(pool, string_resources)
    string_ids_pool <- unique(stats::na.omit(
      string_resources$symbol_to_string[pool]
    ))
    raw_null_cache <- new.env(parent = emptyenv())
    string_null_cache <- new.env(parent = emptyenv())

    get_raw_null <- function(group_size, k_index) {
      key <- paste(k_index, group_size, sep = ":")
      if (!exists(key, envir = raw_null_cache, inherits = FALSE)) {
        seed <- 900000L + dataset_index * 100000L +
          split_row$repeat_idx * 10000L + split_row$fold_idx * 1000L +
          k_index * 100L + group_size
        set.seed(seed)
        random_sets <- replicate(
          n_random, sample(pool, group_size), simplify = FALSE
        )
        value <- list(
          seed = seed,
          open_targets = vapply(
            random_sets,
            function(genes) sum(open_targets_pool[genes]), numeric(1)
          ),
          semantic = vapply(
            random_sets,
            function(genes) mean(semantic_pool[genes]), numeric(1)
          ),
          hallmark = vapply(
            random_sets,
            function(genes) count_graph_edges(hallmark_graph, genes),
            numeric(1)
          )
        )
        assign(key, value, envir = raw_null_cache)
      }
      get(key, envir = raw_null_cache, inherits = FALSE)
    }

    get_string_null <- function(mapped_size, k_index) {
      key <- paste(k_index, mapped_size, sep = ":")
      if (!exists(key, envir = string_null_cache, inherits = FALSE)) {
        seed <- 5900000L + dataset_index * 100000L +
          split_row$repeat_idx * 10000L + split_row$fold_idx * 1000L +
          k_index * 100L + mapped_size
        if (mapped_size < 2L) {
          values <- rep(NA_real_, n_random)
        } else {
          if (length(string_ids_pool) < mapped_size) {
            stop("Too few STRING-mapped candidate genes for ", dataset,
                 " ", split_row$split_id, call. = FALSE)
          }
          set.seed(seed)
          values <- replicate(n_random, {
            identifiers <- sample(string_ids_pool, mapped_size)
            count_graph_edges(string_graph, identifiers)
          })
        }
        assign(
          key, list(seed = seed, values = values),
          envir = string_null_cache
        )
      }
      get(key, envir = string_null_cache, inherits = FALSE)
    }

    for (k_index in seq_along(group_gene_set_sizes)) {
      k <- group_gene_set_sizes[[k_index]]
      split_membership <- membership[
        membership$dataset == dataset &
          membership$repeat_idx == split_row$repeat_idx &
          membership$fold_idx == split_row$fold_idx & membership$k == k,
        , drop = FALSE
      ]

      for (group in groups) {
        genes <- unique(split_membership$gene[
          split_membership$group == group
        ])
        if (any(!genes %in% pool)) {
          stop("A grouped gene is absent from the candidate set for ", dataset,
               " ", split_row$split_id, " k=", k, call. = FALSE)
        }
        group_size <- length(genes)
        if (group_size == 0L) {
          for (metric in c(
            "Open Targets disease-association sum",
            "GO semantic mean", "Hallmark shared-membership edges",
            "STRING protein-association edges"
          )) {
            result_rows[[length(result_rows) + 1L]] <- add_metric_row(
              dataset, record$dataset_group[[1L]], split_row$repeat_idx,
              split_row$fold_idx, split_row$split_id, k, group, group_size,
              metric, 0L, NA_real_, numeric(0), NA_integer_, "empty_group"
            )
          }
          next
        }

        raw_null <- get_raw_null(group_size, k_index)
        result_rows[[length(result_rows) + 1L]] <- add_metric_row(
          dataset, record$dataset_group[[1L]], split_row$repeat_idx,
          split_row$fold_idx, split_row$split_id, k, group, group_size,
          "Open Targets disease-association sum",
          sum(open_targets_pool[genes] > 0),
          sum(open_targets_pool[genes]), raw_null$open_targets,
          raw_null$seed
        )
        result_rows[[length(result_rows) + 1L]] <- add_metric_row(
          dataset, record$dataset_group[[1L]], split_row$repeat_idx,
          split_row$fold_idx, split_row$split_id, k, group, group_size,
          "GO semantic mean", sum(annotation_depth[genes] > 0),
          mean(semantic_pool[genes]), raw_null$semantic, raw_null$seed
        )

        hallmark_reason <- if (group_size < 2L) {
          "fewer_than_two_genes"
        } else ""
        result_rows[[length(result_rows) + 1L]] <- add_metric_row(
          dataset, record$dataset_group[[1L]], split_row$repeat_idx,
          split_row$fold_idx, split_row$split_id, k, group, group_size,
          "Hallmark shared-membership edges",
          sum(genes %in% unique(unlist(hallmark_sets, use.names = FALSE))),
          if (nzchar(hallmark_reason)) NA_real_ else
            count_graph_edges(hallmark_graph, genes),
          if (nzchar(hallmark_reason)) numeric(0) else raw_null$hallmark,
          raw_null$seed, hallmark_reason
        )

        string_identifiers <- unique(stats::na.omit(
          string_resources$symbol_to_string[genes]
        ))
        string_identifiers <- intersect(
          string_identifiers, igraph::V(string_graph)$name
        )
        mapped_size <- length(string_identifiers)
        string_reason <- if (mapped_size < 2L) {
          "fewer_than_two_STRING_mapped_genes"
        } else ""
        string_null <- get_string_null(mapped_size, k_index)
        result_rows[[length(result_rows) + 1L]] <- add_metric_row(
          dataset, record$dataset_group[[1L]], split_row$repeat_idx,
          split_row$fold_idx, split_row$split_id, k, group, group_size,
          "STRING protein-association edges", mapped_size,
          if (nzchar(string_reason)) NA_real_ else
            count_graph_edges(string_graph, string_identifiers),
          if (nzchar(string_reason)) numeric(0) else string_null$values,
          string_null$seed, string_reason
        )
      }
    }
  }
  message("Completed grouped biological assessment for ", dataset)
}

results <- do.call(rbind, result_rows)
rownames(results) <- NULL
expected_rows <- length(datasets) * nrow(split_grid) *
  length(group_gene_set_sizes) * length(groups) * 4L
no_draw_reason <- results$missing_reason %in% c(
  "empty_group", "fewer_than_two_genes",
  "fewer_than_two_STRING_mapped_genes"
)
if (nrow(results) != expected_rows ||
    any(results$group_size < 0L) ||
    any(!results$number_of_random_draws %in% c(0L, n_random)) ||
    any(no_draw_reason & results$number_of_random_draws != 0L) ||
    any(!no_draw_reason & results$number_of_random_draws != n_random)) {
  stop("Grouped biological-assessment table failed validation.",
       call. = FALSE)
}

median_or_na <- function(value) {
  if (all(is.na(value))) NA_real_ else stats::median(value, na.rm = TRUE)
}
summary_table <- do.call(rbind, lapply(split(
  results,
  interaction(
    results$dataset, results$k, results$group, results$metric,
    drop = TRUE, lex.order = TRUE
  )
), function(table) {
  data.frame(
    dataset = table$dataset[[1L]],
    dataset_group = table$dataset_group[[1L]],
    k = table$k[[1L]],
    group = table$group[[1L]],
    metric = table$metric[[1L]],
    n_outer_splits = nrow(table),
    median_group_size = stats::median(table$group_size),
    median_observed_value = median_or_na(table$observed_value),
    median_matched_random_mean = median_or_na(table$matched_random_mean),
    median_observed_to_random_ratio = median_or_na(
      table$observed_to_random_ratio
    ),
    n_missing = sum(!is.finite(table$observed_to_random_ratio)),
    stringsAsFactors = FALSE
  )
}))
rownames(summary_table) <- NULL

atomic_write_csv(
  results,
  file.path(complementarity_results_dir, "gene_group_biology_by_split.csv")
)
atomic_write_csv(
  summary_table,
  file.path(complementarity_results_dir, "gene_group_biology_summary.csv")
)

message(
  "Wrote grouped biological results for ", length(datasets),
  " dataset(s) to ", complementarity_results_dir
)
