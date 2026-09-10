#!/usr/bin/env Rscript

# Exploratory analysis: restrict the candidate genes to the saved DGE top M,
# rank those genes by the saved GeneSelectR score, and evaluate the resulting
# top-k gene sets with the benchmark's unchanged predictive evaluator.
#
# This script fits prediction models and is computationally expensive. It does
# not refit DGE, GeneSelectR, or any other feature-selection method.

source(file.path("redesign", "complementarity_analysis", "00_utils.R"))

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) == 0L ||
    !arguments[[1L]] %in% c("run", "assemble")) {
  stop(
    "Usage:\n",
    "  Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run DATASET\n",
    "  Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R assemble",
    call. = FALSE
  )
}

action <- arguments[[1L]]
candidate_limits <- c(100L, 250L, 500L)
all_evaluator_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
checkpoint_dir <- file.path(
  complementarity_results_dir, "exploratory_dge_prefilter_checkpoints"
)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

if (action == "run") {
  if (length(arguments) != 2L) {
    stop("The run action requires exactly one dataset identifier.",
         call. = FALSE)
  }
  dataset <- normalise_dataset_argument(arguments[[2L]])
  if (length(dataset) != 1L) {
    stop("Run one dataset at a time.", call. = FALSE)
  }
  validate_dataset_inputs(dataset)
  require_packages(c("glmnet", "ranger", "withr", "xgboost", "pROC"))
  source(file.path("redesign", "R", "bio_prior.R"))
  source(file.path("redesign", "R", "evaluator.R"))

  outcome <- load_outcome(dataset)
  saved_prediction <- load_predictive_results(dataset)

  for (index in seq_len(nrow(split_grid))) {
    split_row <- split_grid[index, ]
    checkpoint_path <- file.path(
      checkpoint_dir,
      sprintf(
        "%s_r%d_f%d.rds", dataset,
        split_row$repeat_idx, split_row$fold_idx
      )
    )
    if (file.exists(checkpoint_path)) {
      message("Keeping existing checkpoint: ", checkpoint_path)
      next
    }

    checked <- validate_split_rankings(
      dataset, split_row$repeat_idx, split_row$fold_idx
    )
    split <- checked$split
    if (!all(c("train_raw", "test_raw") %in% names(split))) {
      stop("Saved division lacks train_raw or test_raw: ",
           split_path(dataset, split_row$repeat_idx, split_row$fold_idx),
           call. = FALSE)
    }
    pool <- checked$pool
    gs_ranking <- checked$rankings[["GeneSelectR"]]$gene
    dge_ranking <- checked$rankings[["DGE"]]$gene
    y_train <- outcome[split$train_idx]
    y_test <- outcome[split$test_idx]
    if (length(y_train) != nrow(split$train_raw) ||
        length(y_test) != nrow(split$test_raw)) {
      stop("Outcome and saved division rows do not align for ", dataset,
           " ", split_row$split_id, call. = FALSE)
    }

    standardised <- standardise_split(
      split$train_raw[, pool, drop = FALSE],
      split$test_raw[, pool, drop = FALSE]
    )
    result_rows <- list()
    gene_set_rows <- list()

    for (candidate_limit in candidate_limits) {
      dge_candidates <- dge_ranking[seq_len(candidate_limit)]
      gs_within_dge <- gs_ranking[gs_ranking %in% dge_candidates]
      if (length(gs_within_dge) != candidate_limit ||
          anyDuplicated(gs_within_dge)) {
        stop("DGE prefilter construction failed for ", dataset, " ",
             split_row$split_id, " M=", candidate_limit, call. = FALSE)
      }

      for (k in primary_gene_set_sizes) {
        genes <- top_genes(gs_within_dge, k)
        evaluation_seed <- 420000L + split_row$repeat_idx * 1000L +
          split_row$fold_idx * 100L + match(k, all_evaluator_sizes)
        predictions <- predict_with_ensemble(
          standardised$train[, genes, drop = FALSE], y_train,
          standardised$test[, genes, drop = FALSE],
          random_seed = evaluation_seed
        )
        auc <- bench_auc(y_test, predictions)

        baseline <- saved_prediction[
          saved_prediction$repeat_idx == split_row$repeat_idx &
            saved_prediction$fold_idx == split_row$fold_idx &
            saved_prediction$k == k &
            saved_prediction$method %in% c("GeneSelectR", "DGE"),
          , drop = FALSE
        ]
        if (nrow(baseline) != 2L ||
            length(unique(baseline$random_AUC)) != 1L) {
          stop("Could not match saved benchmark results for ", dataset,
               " ", split_row$split_id, " k=", k, call. = FALSE)
        }
        gs_auc <- baseline$AUC[baseline$method == "GeneSelectR"]
        dge_auc <- baseline$AUC[baseline$method == "DGE"]
        random_auc <- unique(baseline$random_AUC)
        ordinary_gs <- gs_ranking[seq_len(k)]
        ordinary_dge <- dge_ranking[seq_len(k)]

        result_rows[[length(result_rows) + 1L]] <- data.frame(
          dataset = dataset,
          dataset_group = dataset_record(dataset)$dataset_group[[1L]],
          repeat_idx = split_row$repeat_idx,
          fold_idx = split_row$fold_idx,
          split_id = split_row$split_id,
          candidate_limit_M = candidate_limit,
          k = k,
          AUC = auc,
          random_AUC = random_auc,
          AUC_minus_random = auc - random_auc,
          ordinary_geneselectr_AUC = gs_auc,
          ordinary_dge_AUC = dge_auc,
          AUC_minus_ordinary_geneselectr = auc - gs_auc,
          AUC_minus_ordinary_dge = auc - dge_auc,
          shared_with_ordinary_geneselectr = length(intersect(
            genes, ordinary_gs
          )),
          shared_with_ordinary_dge = length(intersect(
            genes, ordinary_dge
          )),
          evaluation_seed = evaluation_seed,
          evaluator_components = paste(
            attr(predictions, "components_used"), collapse = "+"
          ),
          ranking_source = "saved_training_only_rankings",
          analysis_status = "exploratory_dge_prefilter",
          stringsAsFactors = FALSE
        )
        gene_set_rows[[length(gene_set_rows) + 1L]] <- data.frame(
          dataset = dataset,
          dataset_group = dataset_record(dataset)$dataset_group[[1L]],
          repeat_idx = split_row$repeat_idx,
          fold_idx = split_row$fold_idx,
          split_id = split_row$split_id,
          candidate_limit_M = candidate_limit,
          k = k,
          rank = seq_len(k),
          gene = genes,
          stringsAsFactors = FALSE
        )
      }
    }

    checkpoint <- list(
      results = do.call(rbind, result_rows),
      gene_sets = do.call(rbind, gene_set_rows),
      metadata = list(
        dataset = dataset,
        repeat_idx = split_row$repeat_idx,
        fold_idx = split_row$fold_idx,
        candidate_limits = candidate_limits,
        gene_set_sizes = primary_gene_set_sizes,
        split_source = split_path(
          dataset, split_row$repeat_idx, split_row$fold_idx
        ),
        gs_ranking_source = ranking_path(
          dataset, split_row$repeat_idx, split_row$fold_idx, "GeneSelectR"
        ),
        dge_ranking_source = ranking_path(
          dataset, split_row$repeat_idx, split_row$fold_idx, "DGE"
        )
      )
    )
    if (nrow(checkpoint$results) !=
        length(candidate_limits) * length(primary_gene_set_sizes) ||
        any(table(
          checkpoint$gene_sets$candidate_limit_M,
          checkpoint$gene_sets$k
        ) != rep(primary_gene_set_sizes, each = length(candidate_limits)))) {
      stop("Exploratory checkpoint failed validation for ", dataset,
           " ", split_row$split_id, call. = FALSE)
    }
    atomic_save_rds(checkpoint, checkpoint_path)
    message("Wrote checkpoint: ", checkpoint_path)
  }
  message("Completed exploratory predictive evaluation for ", dataset)
}

if (action == "assemble") {
  datasets <- dataset_registry()$dataset
  expected_paths <- unlist(lapply(datasets, function(dataset) {
    file.path(
      checkpoint_dir,
      sprintf(
        "%s_r%d_f%d.rds", dataset,
        split_grid$repeat_idx, split_grid$fold_idx
      )
    )
  }), use.names = FALSE)
  missing <- expected_paths[!file.exists(expected_paths)]
  if (length(missing) > 0L) {
    stop(
      "Cannot assemble: ", length(missing), " checkpoint(s) are missing. ",
      "First missing file: ", missing[[1L]], call. = FALSE
    )
  }

  checkpoints <- lapply(expected_paths, readRDS)
  results <- do.call(rbind, lapply(checkpoints, `[[`, "results"))
  gene_sets <- do.call(rbind, lapply(checkpoints, `[[`, "gene_sets"))
  rownames(results) <- NULL
  rownames(gene_sets) <- NULL

  expected_result_rows <- length(datasets) * nrow(split_grid) *
    length(candidate_limits) * length(primary_gene_set_sizes)
  if (nrow(results) != expected_result_rows ||
      any(!is.finite(results$AUC)) ||
      any(!results$k %in% primary_gene_set_sizes)) {
    stop("Assembled exploratory results failed validation.",
         call. = FALSE)
  }

  stability_rows <- list()
  stability_inputs <- lapply(datasets, function(dataset) {
    lapply(seq_len(nrow(split_grid)), function(index) {
      split_row <- split_grid[index, ]
      list(
        pool = read_split(
          dataset, split_row$repeat_idx, split_row$fold_idx
        )$pools$var2000,
        geneselectr = read_ranking(
          dataset, split_row$repeat_idx, split_row$fold_idx,
          "GeneSelectR"
        )$gene,
        dge = read_ranking(
          dataset, split_row$repeat_idx, split_row$fold_idx, "DGE"
        )$gene
      )
    })
  })
  names(stability_inputs) <- datasets
  groups <- split(
    gene_sets,
    interaction(
      gene_sets$dataset, gene_sets$candidate_limit_M, gene_sets$k,
      drop = TRUE, lex.order = TRUE
    )
  )
  for (table in groups) {
    dataset <- table$dataset[[1L]]
    candidate_limit <- table$candidate_limit_M[[1L]]
    k <- table$k[[1L]]
    split_sets <- split(
      table$gene,
      interaction(
        table$repeat_idx, table$fold_idx,
        drop = TRUE, lex.order = TRUE
      )
    )
    if (length(split_sets) != nrow(split_grid) ||
        any(lengths(split_sets) != k)) {
      stop("Incomplete exploratory gene sets for ", dataset,
           " M=", candidate_limit, " k=", k, call. = FALSE)
    }
    saved_inputs <- stability_inputs[[dataset]]
    pools <- lapply(saved_inputs, `[[`, "pool")
    union_universe <- Reduce(union, pools)
    common_universe <- Reduce(intersect, pools)
    common_sets <- lapply(saved_inputs, function(item) {
      dge_candidates <- top_genes(
        item$dge, candidate_limit, allowed_genes = common_universe
      )
      geneselectr_within_dge <- item$geneselectr[
        item$geneselectr %in% dge_candidates
      ]
      top_genes(geneselectr_within_dge, k)
    })
    common_actual_sets <- lapply(
      split_sets, intersect, y = common_universe
    )
    exact_matrix <- selection_matrix(split_sets, union_universe)
    common_actual_matrix <- selection_matrix(
      common_actual_sets, common_universe
    )
    common_reranked_matrix <- selection_matrix(
      common_sets, common_universe
    )
    set_statistics <- pairwise_set_statistics(split_sets)
    common_actual_nogueira <- compute_repository_nogueira(
      common_actual_matrix
    )

    stability_rows[[length(stability_rows) + 1L]] <- data.frame(
      dataset = dataset,
      dataset_group = dataset_record(dataset)$dataset_group[[1L]],
      candidate_limit_M = candidate_limit,
      k = k,
      n_outer_splits = length(split_sets),
      nogueira_stability = common_actual_nogueira,
      nogueira_fixed_common_genes_actual_selection =
        common_actual_nogueira,
      nogueira_fixed_common_genes_reranked_topk =
        compute_repository_nogueira(common_reranked_matrix),
      nogueira_exact_topk_union = compute_repository_nogueira(exact_matrix),
      mean_number_actual_topk_in_common = mean(lengths(common_actual_sets)),
      mean_fraction_actual_topk_in_common =
        mean(lengths(common_actual_sets)) / k,
      mean_jaccard = set_statistics[["mean_jaccard"]],
      median_jaccard = set_statistics[["median_jaccard"]],
      mean_dice = set_statistics[["mean_dice"]],
      candidate_union_size = length(union_universe),
      candidate_intersection_size = length(common_universe),
      stringsAsFactors = FALSE
    )
  }
  stability <- do.call(rbind, stability_rows)
  rownames(stability) <- NULL

  atomic_write_csv(
    results,
    file.path(
      complementarity_results_dir,
      "exploratory_dge_prefilter_results.csv"
    )
  )
  atomic_write_csv(
    gene_sets,
    file.path(
      complementarity_results_dir,
      "exploratory_dge_prefilter_gene_sets.csv"
    )
  )
  atomic_write_csv(
    stability,
    file.path(
      complementarity_results_dir,
      "exploratory_dge_prefilter_stability.csv"
    )
  )
  message("Assembled exploratory DGE-prefilter results.")
}
