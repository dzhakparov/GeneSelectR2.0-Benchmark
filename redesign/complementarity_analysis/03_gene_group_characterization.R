#!/usr/bin/env Rscript

# Describe genes that are ranked in the top k by both methods, by DGE only, or
# by GeneSelectR only. These are ranking groups and do not imply statistical
# significance, causality, or biomarker validation.

source(file.path("redesign", "complementarity_analysis", "00_utils.R"))

input_path <- file.path(
  complementarity_results_dir, "dge_geneselectr_gene_scores.csv"
)
if (!file.exists(input_path)) {
  stop(
    "Missing input: ", input_path,
    ". Run 02_dge_geneselectr_concordance.R first.", call. = FALSE
  )
}

datasets <- normalise_dataset_argument(commandArgs(trailingOnly = TRUE))
gene_scores <- utils::read.csv(
  input_path, stringsAsFactors = FALSE, check.names = FALSE
)
gene_scores <- gene_scores[gene_scores$dataset %in% datasets, , drop = FALSE]

required <- c(
  "dataset", "dataset_group", "repeat_idx", "fold_idx", "split_id", "gene",
  "geneselectr_rank", "geneselectr_final_score",
  "geneselectr_internal_recurrence", "geneselectr_calibrated_recurrence",
  "geneselectr_raw_utility", "geneselectr_calibrated_utility", "dge_rank",
  "dge_rank_percentile", "dge_ranking_statistic", "dge_t_statistic",
  "dge_p_value", "dge_fdr", "dge_selected_bh_0_05",
  "dge_log_fold_change"
)
missing <- setdiff(required, names(gene_scores))
if (length(missing) > 0L) {
  stop("Input is missing columns: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

membership_rows <- list()

split_tables <- split(
  gene_scores,
  interaction(
    gene_scores$dataset, gene_scores$repeat_idx, gene_scores$fold_idx,
    drop = TRUE, lex.order = TRUE
  )
)

for (table in split_tables) {
  if (nrow(table) != 2000L || anyDuplicated(table$gene)) {
    stop(
      "Expected 2,000 unique candidate genes for ", table$dataset[[1L]],
      " ", table$split_id[[1L]], call. = FALSE
    )
  }
  for (k in group_gene_set_sizes) {
    gs_top <- table$geneselectr_rank <= k
    dge_top <- table$dge_rank <= k
    keep <- gs_top | dge_top
    group <- ifelse(
      gs_top & dge_top, "SHARED",
      ifelse(dge_top, "DGE_ONLY", "GS_ONLY")
    )
    output <- table[keep, required, drop = FALSE]
    output$k <- k
    output$group <- group[keep]
    output <- output[, c(
      "dataset", "dataset_group", "repeat_idx", "fold_idx", "split_id", "k",
      "group", "gene", "geneselectr_rank", "dge_rank",
      "geneselectr_final_score", "geneselectr_internal_recurrence",
      "geneselectr_calibrated_recurrence", "geneselectr_raw_utility",
      "geneselectr_calibrated_utility", "dge_rank_percentile",
      "dge_ranking_statistic", "dge_t_statistic", "dge_p_value", "dge_fdr",
      "dge_selected_bh_0_05",
      "dge_log_fold_change"
    )]

    counts <- table(output$group)
    shared_count <- if ("SHARED" %in% names(counts)) counts[["SHARED"]] else 0L
    dge_only_count <- if ("DGE_ONLY" %in% names(counts)) {
      counts[["DGE_ONLY"]]
    } else 0L
    gs_only_count <- if ("GS_ONLY" %in% names(counts)) {
      counts[["GS_ONLY"]]
    } else 0L
    if (shared_count + dge_only_count != k ||
        shared_count + gs_only_count != k ||
        nrow(output) != 2L * k - shared_count) {
      stop(
        "Gene-group counts failed validation for ", table$dataset[[1L]],
        " ", table$split_id[[1L]], " k=", k, call. = FALSE
      )
    }
    membership_rows[[length(membership_rows) + 1L]] <- output
  }
}

membership <- do.call(rbind, membership_rows)
rownames(membership) <- NULL

membership$shared_count <- as.integer(membership$group == "SHARED")
membership$dge_only_count <- as.integer(membership$group == "DGE_ONLY")
membership$gs_only_count <- as.integer(
  membership$group == "GS_ONLY"
)

count_summary <- stats::aggregate(
  membership[, c(
    "shared_count", "dge_only_count", "gs_only_count"
  )],
  by = membership[, c("dataset", "dataset_group", "k", "gene")],
  FUN = sum
)
names(count_summary)[names(count_summary) == "shared_count"] <-
  "n_splits_shared"
names(count_summary)[names(count_summary) == "dge_only_count"] <-
  "n_splits_dge_only"
names(count_summary)[names(count_summary) == "gs_only_count"] <-
  "n_splits_gs_only"

mean_columns <- c(
  "geneselectr_final_score", "geneselectr_internal_recurrence",
  "geneselectr_calibrated_recurrence", "geneselectr_raw_utility",
  "geneselectr_calibrated_utility", "dge_rank_percentile",
  "dge_ranking_statistic", "dge_t_statistic",
  "dge_log_fold_change"
)
gene_means <- stats::aggregate(
  gene_scores[, mean_columns],
  by = gene_scores[, c("dataset", "dataset_group", "gene")],
  FUN = function(value) {
    if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE)
  }
)
names(gene_means)[match(mean_columns, names(gene_means))] <- paste0(
  "mean_", mean_columns
)
candidate_counts <- stats::aggregate(
  data.frame(n_candidate_splits = rep.int(1L, nrow(gene_scores))),
  by = gene_scores[, c("dataset", "dataset_group", "gene")],
  FUN = sum
)
gene_means <- merge(
  candidate_counts, gene_means,
  by = c("dataset", "dataset_group", "gene"),
  all = TRUE, sort = FALSE
)

summary_table <- merge(
  count_summary, gene_means,
  by = c("dataset", "dataset_group", "gene"), all.x = TRUE, sort = FALSE
)
summary_table$n_outer_splits <- nrow(split_grid)
summary_table$n_splits_any_topk <- with(
  summary_table,
  n_splits_shared + n_splits_dge_only + n_splits_gs_only
)
summary_table$proportion_splits_shared <-
  summary_table$n_splits_shared / summary_table$n_outer_splits
summary_table$proportion_splits_dge_only <-
  summary_table$n_splits_dge_only / summary_table$n_outer_splits
summary_table$proportion_splits_gs_only <-
  summary_table$n_splits_gs_only / summary_table$n_outer_splits
summary_table$proportion_available_splits_shared <-
  summary_table$n_splits_shared / summary_table$n_candidate_splits
summary_table$proportion_available_splits_dge_only <-
  summary_table$n_splits_dge_only / summary_table$n_candidate_splits
summary_table$proportion_available_splits_gs_only <-
  summary_table$n_splits_gs_only / summary_table$n_candidate_splits
summary_table <- summary_table[order(
  summary_table$dataset, summary_table$k,
  -summary_table$n_splits_shared,
  -summary_table$n_splits_gs_only,
  -summary_table$n_splits_dge_only,
  summary_table$gene
), ]
rownames(summary_table) <- NULL

recurrent_long <- rbind(
  data.frame(
    summary_table[, c("dataset", "dataset_group", "k", "gene")],
    group = "SHARED", n_splits = summary_table$n_splits_shared
  ),
  data.frame(
    summary_table[, c("dataset", "dataset_group", "k", "gene")],
    group = "DGE_ONLY", n_splits = summary_table$n_splits_dge_only
  ),
  data.frame(
    summary_table[, c("dataset", "dataset_group", "k", "gene")],
    group = "GS_ONLY",
    n_splits = summary_table$n_splits_gs_only
  )
)
recurrent_long <- recurrent_long[recurrent_long$n_splits > 0L, ]
recurrent_long <- recurrent_long[order(
  recurrent_long$dataset, recurrent_long$k, recurrent_long$group,
  -recurrent_long$n_splits, recurrent_long$gene
), ]
rownames(recurrent_long) <- NULL

membership_output <- membership[, setdiff(
  names(membership),
  c("shared_count", "dge_only_count", "gs_only_count")
)]

atomic_write_csv(
  membership_output,
  file.path(complementarity_results_dir, "gene_group_membership.csv")
)
atomic_write_csv(
  summary_table,
  file.path(complementarity_results_dir, "gene_group_summary.csv")
)
atomic_write_csv(
  recurrent_long,
  file.path(complementarity_results_dir, "gene_group_recurrence_long.csv")
)

message(
  "Wrote shared, DGE-only and GeneSelectR-only gene summaries for ",
  length(datasets), " dataset(s) to ", complementarity_results_dir
)
