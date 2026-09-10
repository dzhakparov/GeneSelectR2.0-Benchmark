#!/usr/bin/env Rscript

# Compare the stability of the final top-k gene sets across the 15 saved
# train/test divisions. GeneSelectR's internal recurrence value is not used as
# a substitute for this across-division stability calculation.

source(file.path("redesign", "complementarity_analysis", "00_utils.R"))

datasets <- normalise_dataset_argument(commandArgs(trailingOnly = TRUE))
dir.create(complementarity_results_dir, recursive = TRUE, showWarnings = FALSE)

predictive_results <- load_all_predictive_results(datasets)

stability_rows <- list()
gene_recurrence_rows <- list()
top_gene_rows <- list()
pool_rows <- list()

for (dataset in datasets) {
  validate_dataset_inputs(dataset)
  record <- dataset_record(dataset)

  saved <- lapply(seq_len(nrow(split_grid)), function(index) {
    split_row <- split_grid[index, ]
    validate_split_rankings(
      dataset, split_row$repeat_idx, split_row$fold_idx
    )
  })
  names(saved) <- split_grid$split_id

  pools <- lapply(saved, `[[`, "pool")
  common_universe <- Reduce(intersect, pools)
  union_universe <- Reduce(union, pools)
  pool_similarity <- pairwise_set_statistics(pools)
  constant_pool <- length(common_universe) == length(union_universe)

  if (length(common_universe) < max(primary_gene_set_sizes)) {
    stop(
      "Fewer than 50 genes occur in every candidate set for ", dataset,
      ". The fixed-universe sensitivity calculation cannot be performed.",
      call. = FALSE
    )
  }

  pool_rows[[length(pool_rows) + 1L]] <- data.frame(
    dataset = dataset,
    dataset_group = record$dataset_group[[1L]],
    n_outer_splits = length(pools),
    candidate_genes_per_split = unique(lengths(pools)),
    candidate_union_size = length(union_universe),
    candidate_intersection_size = length(common_universe),
    candidate_pool_constant = constant_pool,
    mean_pairwise_candidate_jaccard = pool_similarity[["mean_jaccard"]],
    minimum_pairwise_candidate_jaccard = min(apply(
      utils::combn(seq_along(pools), 2L), 2L, function(indices) {
        length(intersect(pools[[indices[[1L]]]], pools[[indices[[2L]]]])) /
          length(union(pools[[indices[[1L]]]], pools[[indices[[2L]]]]))
      }
    )),
    stringsAsFactors = FALSE
  )

  availability_count <- table(unlist(pools, use.names = FALSE))
  availability_count <- as.integer(availability_count[union_universe])

  dataset_prediction <- predictive_results[
    predictive_results$dataset == dataset, , drop = FALSE
  ]

  for (method in analysis_methods) {
    method_rankings <- lapply(saved, function(item) {
      item$rankings[[method]]$gene
    })

    for (k in primary_gene_set_sizes) {
      exact_gene_sets <- lapply(method_rankings, top_genes, k = k)
      common_actual_gene_sets <- lapply(
        exact_gene_sets, intersect, y = common_universe
      )
      common_reranked_gene_sets <- lapply(
        method_rankings, top_genes, k = k, allowed_genes = common_universe
      )

      if (any(lengths(exact_gene_sets) != k) ||
          any(vapply(exact_gene_sets, anyDuplicated, integer(1)) > 0L)) {
        stop("Invalid top-k gene set for ", dataset, " ", method,
             " k=", k, call. = FALSE)
      }

      exact_matrix <- selection_matrix(exact_gene_sets, union_universe)
      common_actual_matrix <- selection_matrix(
        common_actual_gene_sets, common_universe
      )
      common_reranked_matrix <- selection_matrix(
        common_reranked_gene_sets, common_universe
      )
      set_statistics <- pairwise_set_statistics(exact_gene_sets)
      exact_nogueira <- compute_repository_nogueira(exact_matrix)
      common_actual_nogueira <- compute_repository_nogueira(
        common_actual_matrix
      )
      common_reranked_nogueira <- compute_repository_nogueira(
        common_reranked_matrix
      )

      prediction_subset <- dataset_prediction[
        dataset_prediction$method == method & dataset_prediction$k == k,
        , drop = FALSE
      ]
      if (nrow(prediction_subset) != nrow(split_grid)) {
        stop("Incomplete predictive results for ", dataset, " ", method,
             " k=", k, call. = FALSE)
      }

      stability_rows[[length(stability_rows) + 1L]] <- data.frame(
        dataset = dataset,
        dataset_group = record$dataset_group[[1L]],
        method = method,
        k = k,
        n_outer_splits = length(exact_gene_sets),
        nogueira_stability = common_actual_nogueira,
        nogueira_fixed_common_genes_actual_selection =
          common_actual_nogueira,
        nogueira_fixed_common_genes_reranked_topk =
          common_reranked_nogueira,
        nogueira_exact_topk_union = exact_nogueira,
        mean_number_actual_topk_in_common = mean(lengths(
          common_actual_gene_sets
        )),
        mean_fraction_actual_topk_in_common = mean(lengths(
          common_actual_gene_sets
        )) / k,
        mean_jaccard = set_statistics[["mean_jaccard"]],
        median_jaccard = set_statistics[["median_jaccard"]],
        mean_dice = set_statistics[["mean_dice"]],
        mean_overlap_coefficient = set_statistics[[
          "mean_overlap_coefficient"
        ]],
        mean_auc = mean(prediction_subset$AUC),
        mean_auc_minus_random = mean(prediction_subset$AUC_minus_random),
        candidate_union_size = length(union_universe),
        candidate_intersection_size = length(common_universe),
        candidate_pool_constant = constant_pool,
        stringsAsFactors = FALSE
      )

      selection_count <- rowSums(exact_matrix)
      gene_recurrence_rows[[length(gene_recurrence_rows) + 1L]] <- data.frame(
        dataset = dataset,
        dataset_group = record$dataset_group[[1L]],
        method = method,
        k = k,
        gene = union_universe,
        n_outer_splits = length(exact_gene_sets),
        n_candidate_splits = availability_count,
        n_selected = as.integer(selection_count),
        proportion_selected = selection_count / length(exact_gene_sets),
        proportion_selected_all_splits = selection_count /
          length(exact_gene_sets),
        proportion_selected_when_candidate = selection_count /
          availability_count,
        stringsAsFactors = FALSE
      )

      for (index in seq_along(exact_gene_sets)) {
        split_row <- split_grid[index, ]
        genes <- exact_gene_sets[[index]]
        top_gene_rows[[length(top_gene_rows) + 1L]] <- data.frame(
          dataset = dataset,
          dataset_group = record$dataset_group[[1L]],
          repeat_idx = split_row$repeat_idx,
          fold_idx = split_row$fold_idx,
          split_id = split_row$split_id,
          method = method,
          k = k,
          rank = seq_len(k),
          gene = genes,
          candidate_pool_size = length(pools[[index]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  message("Completed gene-set stability inputs for ", dataset)
}

stability <- do.call(rbind, stability_rows)
gene_recurrence <- do.call(rbind, gene_recurrence_rows)
top_gene_sets <- do.call(rbind, top_gene_rows)
candidate_diagnostics <- do.call(rbind, pool_rows)
rownames(stability) <- NULL
rownames(gene_recurrence) <- NULL
rownames(top_gene_sets) <- NULL
rownames(candidate_diagnostics) <- NULL

expected_stability_rows <- length(datasets) * length(analysis_methods) *
  length(primary_gene_set_sizes)
if (nrow(stability) != expected_stability_rows ||
    any(!is.finite(stability$mean_jaccard)) ||
    any(!is.finite(stability$mean_auc))) {
  stop("Final stability table failed validation.", call. = FALSE)
}

atomic_write_csv(
  stability,
  file.path(complementarity_results_dir, "gene_set_stability_summary.csv")
)
atomic_write_csv(
  gene_recurrence,
  file.path(complementarity_results_dir, "gene_cross_split_recurrence.csv")
)
atomic_write_csv(
  top_gene_sets,
  file.path(complementarity_results_dir, "top_gene_sets_by_split.csv")
)
atomic_write_csv(
  candidate_diagnostics,
  file.path(complementarity_results_dir, "candidate_gene_set_diagnostics.csv")
)

message(
  "Wrote stability results for ", length(datasets), " dataset(s) to ",
  complementarity_results_dir
)
