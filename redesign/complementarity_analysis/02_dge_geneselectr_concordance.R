#!/usr/bin/env Rscript

# Compare saved differential-expression and GeneSelectR rankings within each
# outer training division. The saved DGE files contain rank order only. The
# exact benchmark DGE calculation is repeated on the saved training data and
# its complete gene order must match the saved order.

source(file.path("redesign", "complementarity_analysis", "00_utils.R"))
source(file.path("redesign", "R", "bio_prior.R"))

datasets <- normalise_dataset_argument(commandArgs(trailingOnly = TRUE))
dir.create(complementarity_results_dir, recursive = TRUE, showWarnings = FALSE)

gene_rows <- list()
correlation_rows <- list()
overlap_rows <- list()

required_gs_columns <- c(
  "gene", "final_score", "raw_score", "pi_raw", "pi_scored", "u",
  "u_mi", "u_scored"
)

for (dataset in datasets) {
  validate_dataset_inputs(dataset)
  record <- dataset_record(dataset)
  outcome <- load_outcome(dataset)

  for (index in seq_len(nrow(split_grid))) {
    split_row <- split_grid[index, ]
    checked <- validate_split_rankings(
      dataset, split_row$repeat_idx, split_row$fold_idx
    )
    gs <- checked$rankings[["GeneSelectR"]]
    dge <- checked$rankings[["DGE"]]
    split <- checked$split

    missing_gs_columns <- setdiff(required_gs_columns, names(gs))
    if (length(missing_gs_columns) > 0L) {
      stop(
        "GeneSelectR ranking is missing required columns for ", dataset,
        " ", split_row$split_id, ": ",
        paste(missing_gs_columns, collapse = ", "), call. = FALSE
      )
    }
    if (!setequal(gs$gene, dge$gene) ||
        !setequal(gs$gene, checked$pool)) {
      stop(
        "DGE, GeneSelectR and candidate genes differ for ", dataset, " ",
        split_row$split_id, call. = FALSE
      )
    }

    if (!"train_raw" %in% names(split)) {
      stop("Saved division has no training matrix for ", dataset, " ",
           split_row$split_id, call. = FALSE)
    }
    y_train <- outcome[split$train_idx]
    if (length(y_train) != nrow(split$train_raw)) {
      stop("Training labels do not align with the saved division for ",
           dataset, " ", split_row$split_id, call. = FALSE)
    }
    standardised_train <- standardise_split(
      split$train_raw[, checked$pool, drop = FALSE], NULL
    )$train
    dge_metrics <- calculate_benchmark_dge_metrics(
      standardised_train, y_train
    )
    if (!identical(dge_metrics$gene, dge$gene)) {
      first_difference <- which(dge_metrics$gene != dge$gene)[[1L]]
      stop(
        "Recomputed DGE ranking does not match the saved ranking for ",
        dataset, " ", split_row$split_id, " at rank ", first_difference,
        call. = FALSE
      )
    }

    total <- nrow(gs)
    gs_rank <- seq_len(total)
    dge_rank <- match(gs$gene, dge$gene)
    dge_metric_index <- match(gs$gene, dge_metrics$gene)
    if (anyNA(dge_rank) || anyDuplicated(dge_rank)) {
      stop("DGE rank matching failed for ", dataset, " ",
           split_row$split_id, call. = FALSE)
    }
    gs_percentile <- high_is_best_rank_percentile(gs_rank, total)
    dge_percentile <- high_is_best_rank_percentile(dge_rank, total)

    gene_rows[[length(gene_rows) + 1L]] <- data.frame(
      dataset = dataset,
      dataset_group = record$dataset_group[[1L]],
      repeat_idx = split_row$repeat_idx,
      fold_idx = split_row$fold_idx,
      split_id = split_row$split_id,
      gene = gs$gene,
      candidate_pool_size = total,
      geneselectr_rank = gs_rank,
      geneselectr_rank_percentile = gs_percentile,
      geneselectr_final_score = gs$final_score,
      geneselectr_raw_combined_score = gs$raw_score,
      geneselectr_internal_recurrence = gs$pi_raw,
      geneselectr_calibrated_recurrence = gs$pi_scored,
      geneselectr_raw_utility = gs$u,
      geneselectr_mutual_information = gs$u_mi,
      geneselectr_calibrated_utility = gs$u_scored,
      dge_rank = dge_rank,
      dge_rank_percentile = dge_percentile,
      dge_ranking_statistic = dge_metrics$dge_ranking_statistic[
        dge_metric_index
      ],
      dge_t_statistic = dge_metrics$dge_t_statistic[dge_metric_index],
      dge_log_fold_change = NA_real_,
      dge_p_value = dge_metrics$dge_p_value[dge_metric_index],
      dge_fdr = dge_metrics$dge_fdr[dge_metric_index],
      dge_selected_bh_0_05 = dge_metrics$dge_selected_bh_0_05[
        dge_metric_index
      ],
      dge_saved_information = "rank_recomputed_statistic_p_value_and_fdr",
      stringsAsFactors = FALSE
    )

    correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
      dataset = dataset,
      dataset_group = record$dataset_group[[1L]],
      repeat_idx = split_row$repeat_idx,
      fold_idx = split_row$fold_idx,
      split_id = split_row$split_id,
      n_candidate_genes = total,
      spearman_dge_vs_geneselectr_rank = safe_spearman(
        dge_percentile, gs_percentile
      ),
      spearman_dge_vs_geneselectr_score = safe_spearman(
        dge_percentile, gs$final_score
      ),
      spearman_dge_vs_internal_recurrence = safe_spearman(
        dge_percentile, gs$pi_raw
      ),
      spearman_dge_vs_calibrated_recurrence = safe_spearman(
        dge_percentile, gs$pi_scored
      ),
      spearman_dge_vs_raw_utility = safe_spearman(
        dge_percentile, gs$u
      ),
      spearman_dge_vs_calibrated_utility = safe_spearman(
        dge_percentile, gs$u_scored
      ),
      stringsAsFactors = FALSE
    )

    for (k in primary_gene_set_sizes) {
      gs_top <- gs$gene[seq_len(k)]
      dge_top <- dge$gene[seq_len(k)]
      shared <- intersect(gs_top, dge_top)
      union_size <- length(union(gs_top, dge_top))
      overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(
        dataset = dataset,
        dataset_group = record$dataset_group[[1L]],
        repeat_idx = split_row$repeat_idx,
        fold_idx = split_row$fold_idx,
        split_id = split_row$split_id,
        k = k,
        n_shared = length(shared),
        jaccard = length(shared) / union_size,
        fraction_geneselectr_also_dge = length(shared) / k,
        fraction_dge_also_geneselectr = length(shared) / k,
        stringsAsFactors = FALSE
      )
    }
  }
  message("Completed DGE and GeneSelectR comparison for ", dataset)
}

gene_level <- do.call(rbind, gene_rows)
correlations <- do.call(rbind, correlation_rows)
overlap <- do.call(rbind, overlap_rows)
rownames(gene_level) <- NULL
rownames(correlations) <- NULL
rownames(overlap) <- NULL

expected_split_rows <- length(datasets) * nrow(split_grid)
if (nrow(correlations) != expected_split_rows ||
    nrow(overlap) != expected_split_rows * length(primary_gene_set_sizes) ||
    any(gene_level$candidate_pool_size != 2000L) ||
    any(overlap$n_shared < 0L | overlap$n_shared > overlap$k)) {
  stop("DGE and GeneSelectR output tables failed validation.",
       call. = FALSE)
}

correlation_summary <- do.call(rbind, lapply(
  split(correlations, correlations$dataset), function(table) {
    numeric_columns <- grep("^spearman_", names(table), value = TRUE)
    summarise_finite <- function(value, function_value) {
      value <- value[is.finite(value)]
      if (length(value) == 0L) NA_real_ else function_value(value)
    }
    output <- data.frame(
      dataset = table$dataset[[1L]],
      dataset_group = table$dataset_group[[1L]],
      n_outer_splits = nrow(table),
      stringsAsFactors = FALSE
    )
    for (column in numeric_columns) {
      output[[paste0("median_", column)]] <- summarise_finite(
        table[[column]], stats::median
      )
      output[[paste0("minimum_", column)]] <- summarise_finite(
        table[[column]], min
      )
      output[[paste0("maximum_", column)]] <- summarise_finite(
        table[[column]], max
      )
    }
    output
  }
))
rownames(correlation_summary) <- NULL

overlap_summary <- do.call(rbind, lapply(
  split(overlap, interaction(overlap$dataset, overlap$k, drop = TRUE)),
  function(table) {
    data.frame(
      dataset = table$dataset[[1L]],
      dataset_group = table$dataset_group[[1L]],
      k = table$k[[1L]],
      n_outer_splits = nrow(table),
      mean_shared = mean(table$n_shared),
      median_shared = stats::median(table$n_shared),
      mean_jaccard = mean(table$jaccard),
      median_jaccard = stats::median(table$jaccard),
      mean_fraction_geneselectr_also_dge = mean(
        table$fraction_geneselectr_also_dge
      ),
      mean_fraction_dge_also_geneselectr = mean(
        table$fraction_dge_also_geneselectr
      ),
      stringsAsFactors = FALSE
    )
  }
))
rownames(overlap_summary) <- NULL

atomic_write_csv(
  gene_level,
  file.path(complementarity_results_dir, "dge_geneselectr_gene_scores.csv")
)
atomic_write_csv(
  correlations,
  file.path(complementarity_results_dir, "dge_geneselectr_split_correlations.csv")
)
atomic_write_csv(
  correlation_summary,
  file.path(complementarity_results_dir, "dge_geneselectr_correlation_summary.csv")
)
atomic_write_csv(
  overlap,
  file.path(complementarity_results_dir, "dge_geneselectr_topk_overlap.csv")
)
atomic_write_csv(
  overlap_summary,
  file.path(complementarity_results_dir, "dge_geneselectr_overlap_summary.csv")
)

message(
  "Wrote DGE and GeneSelectR concordance results for ", length(datasets),
  " dataset(s) to ", complementarity_results_dir
)
