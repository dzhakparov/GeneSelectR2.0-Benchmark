#!/usr/bin/env Rscript

# Recompute ONLY the Open Targets axis of the older-seven biology assessment
# using the dense association files (ot_seeds_<id>_n3000_s0.rds), at two
# target definitions: association score >= 0.05 (primary, denser) and
# >= 0.1 (continuity with the frozen n100 runs).
#
# Everything else mirrors redesign/run_older7_biology.R EXACTLY: same pools
# (split_r*_f*.rds$pools$var2000), same rankings, same panel sizes, same
# n_null = 1000, same null-panel construction (sample from the split pool
# with withr::with_seed), the same null seed formula
# (940000 + r*10000 + f*100 + k_index), and the same empirical summary
# (expected / ratio / z / empirical p). The only change versus the frozen
# runs is the association data. Scores below the cutoff are set to 0, which
# reproduces the frozen semantics (targets = genes with positive score).
#
# Outputs (NEW directory, frozen results untouched):
#   redesign/results_corrected/older7_biology_ot_dense_2026-09-03/
#     ot_dense_by_split.csv              long: dataset x split x method x k x cutoff
#     ot_dense_axis_summary_k10_k20_k50.csv
#     ot_dense_overall_ranks.csv
#     ot_dense_coverage.csv

suppressPackageStartupMessages(library(withr))
options(warn = 1)

validation_datasets_old <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683"
)
datasets <- c(validation_datasets_old, "imvigor210", "sosall")
methods <- c(
  "GS_full_ungrouped", "DGE", "LASSO", "ElasticNet", "mRMR", "Boruta",
  "RF_importance"
)
panel_sizes <- c(10L, 20L, 50L, 100L, 200L, 500L)
n_null <- 1000L
cutoffs <- c(0.05, 0.1)

results_root <- file.path("redesign", "results_corrected")
output_dir <- file.path(results_root, "older7_biology_ot_dense_2026-09-03")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("analysis", "config.R"))
config <- get_biology_config("older7")

# Identical to run_older7_biology.R: per-symbol max across association files.
load_association_scores <- function(paths) {
  tables <- lapply(paths, function(path) {
    table <- readRDS(path)
    if (!all(c("symbol", "score") %in% names(table)) ||
        anyDuplicated(table$symbol) || any(!is.finite(table$score))) {
      stop("Invalid Open Targets file: ", path, call. = FALSE)
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

# Identical to run_older7_biology.R empirical_summary().
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

rows_all <- list()

for (dataset in datasets) {
  is_validation <- dataset %in% validation_datasets_old
  split_dir <- if (is_validation) {
    file.path(results_root, "validation_benchmark", dataset)
  } else {
    file.path(results_root, "grouped_benchmark", dataset)
  }
  dataset_dir <- if (is_validation) split_dir else {
    file.path(results_root, "full_recipe", dataset)
  }

  config_row <- config[config$dataset == dataset, , drop = FALSE]
  stopifnot(nrow(config_row) == 1L)
  dense_file <- file.path(
    "data", "r_user_cache", "R", "GeneSelectR",
    sprintf("ot_seeds_%s_n3000_s0.rds",
            gsub("[^A-Za-z0-9]", "", config_row$ontology_ids[[1L]]))
  )
  if (!file.exists(dense_file)) {
    stop("Missing dense Open Targets file: ", dense_file, call. = FALSE)
  }
  association_scores <- load_association_scores(dense_file)

  # Same ranking-path rule as run_older7_biology.R.
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

  for (repeat_idx in 1:3) {
    for (fold_idx in 1:5) {
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

      # Scores restricted to the pool, one vector per cutoff. Same rule as
      # the frozen run: absent from the association table -> 0.
      pool_scores_full <- stats::setNames(association_scores[pool], pool)
      pool_scores_full[is.na(pool_scores_full)] <- 0
      pool_scores_by_cutoff <- lapply(cutoffs, function(cutoff) {
        scores <- pool_scores_full
        scores[scores < cutoff] <- 0
        scores
      })
      names(pool_scores_by_cutoff) <- as.character(cutoffs)

      for (panel_size in panel_sizes) {
        # Identical null construction and seed formula.
        null_seed <- 940000L + repeat_idx * 10000L + fold_idx * 100L +
          match(panel_size, panel_sizes)
        null_panels <- withr::with_seed(null_seed, lapply(
          seq_len(n_null), function(index) sample(pool, panel_size)
        ))
        null_sums_by_cutoff <- lapply(pool_scores_by_cutoff, function(scores) {
          vapply(null_panels, function(panel) sum(scores[panel]), numeric(1))
        })

        for (method in methods) {
          panel <- head(rankings[[method]], panel_size)
          for (cutoff in cutoffs) {
            scores <- pool_scores_by_cutoff[[as.character(cutoff)]]
            summary <- empirical_summary(
              sum(scores[panel]), null_sums_by_cutoff[[as.character(cutoff)]]
            )
            rows_all[[length(rows_all) + 1L]] <- data.frame(
              dataset = dataset,
              repeat_idx = repeat_idx,
              fold_idx = fold_idx,
              method = method,
              k = panel_size,
              cutoff = cutoff,
              open_targets_sum = sum(scores[panel]),
              open_targets_expected = summary[["expected"]],
              open_targets_enrichment = summary[["ratio"]],
              open_targets_z = summary[["z"]],
              open_targets_empirical_p = summary[["p"]],
              open_targets_overlap = sum(scores[panel] > 0),
              open_targets_pool_targets = sum(scores > 0),
              pool_size = length(pool),
              n_null = n_null,
              null_seed = null_seed,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  cat(sprintf("[%s] dense Open Targets axis complete\n", dataset))
}

by_split <- do.call(rbind, rows_all)
rownames(by_split) <- NULL
stopifnot(
  nrow(by_split) ==
    length(datasets) * 15L * length(methods) * length(panel_sizes) *
      length(cutoffs),
  all(by_split$pool_size == 2000L),
  all(by_split$open_targets_overlap <= by_split$k),
  all(is.finite(by_split$open_targets_sum)),
  all(is.finite(by_split$open_targets_expected)),
  all(is.finite(by_split$open_targets_empirical_p))
)
if (any(!is.finite(by_split$open_targets_enrichment[
  by_split$open_targets_expected > 0
]))) {
  stop("Non-finite enrichment where the null expectation is positive.",
       call. = FALSE)
}
write.csv(by_split, file.path(output_dir, "ot_dense_by_split.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Summaries mirroring
#   integrated_interpretation_2026-09-03/dev7_biology_axis_summary_k10_k20_k50.csv
#   integrated_interpretation_2026-09-03/dev7_biology_overall_ranks.csv
# restricted to the Open Targets axis, with cutoff as an extra key.
# -----------------------------------------------------------------------------
headline_k <- c(10L, 20L, 50L)
hl <- by_split[by_split$k %in% headline_k, ]
axis_summary <- do.call(rbind, lapply(split(
  hl, interaction(hl$dataset, hl$method, hl$k, hl$cutoff, drop = TRUE)
), function(v) {
  data.frame(
    dataset = v$dataset[1L], method = v$method[1L], k = v$k[1L],
    cutoff = v$cutoff[1L], n_splits = nrow(v),
    open_targets_median_enrichment = stats::median(
      v$open_targets_enrichment, na.rm = TRUE
    ),
    open_targets_median_overlap = stats::median(v$open_targets_overlap),
    stringsAsFactors = FALSE
  )
}))
rownames(axis_summary) <- NULL
axis_summary$open_targets_rank <- ave(
  -axis_summary$open_targets_median_enrichment,
  interaction(axis_summary$dataset, axis_summary$k, axis_summary$cutoff,
              drop = TRUE),
  FUN = rank
)
axis_summary <- axis_summary[
  order(axis_summary$cutoff, axis_summary$dataset, axis_summary$k,
        axis_summary$method),
]
write.csv(
  axis_summary,
  file.path(output_dir, "ot_dense_axis_summary_k10_k20_k50.csv"),
  row.names = FALSE
)

overall_ranks <- do.call(rbind, lapply(cutoffs, function(cutoff) {
  v <- axis_summary[axis_summary$cutoff == cutoff, ]
  per_method <- stats::aggregate(
    v$open_targets_median_enrichment,
    by = list(method = v$method), FUN = mean
  )
  names(per_method)[2] <- "dataset_balanced_mean_median_enrichment"
  per_method$axis <- sprintf("open_targets_s%g", cutoff)
  per_method$rank <- rank(
    -per_method$dataset_balanced_mean_median_enrichment
  )
  per_method[, c("axis", "method",
                 "dataset_balanced_mean_median_enrichment", "rank")]
}))
write.csv(overall_ranks,
          file.path(output_dir, "ot_dense_overall_ranks.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# Coverage table mirroring dev7_opentargets_coverage.csv, plus the across-
# split range of pool targets requested for the run report.
# -----------------------------------------------------------------------------
coverage <- do.call(rbind, lapply(split(
  by_split, interaction(by_split$dataset, by_split$k, by_split$cutoff,
                        drop = TRUE)
), function(v) {
  data.frame(
    dataset = v$dataset[1L], k = v$k[1L], cutoff = v$cutoff[1L],
    median_pool_targets = stats::median(v$open_targets_pool_targets),
    min_pool_targets = min(v$open_targets_pool_targets),
    max_pool_targets = max(v$open_targets_pool_targets),
    median_panel_overlap = stats::median(v$open_targets_overlap),
    fraction_splits_zero_overlap = mean(v$open_targets_overlap == 0),
    stringsAsFactors = FALSE
  )
}))
rownames(coverage) <- NULL
coverage <- coverage[order(coverage$cutoff, coverage$dataset, coverage$k), ]
write.csv(coverage, file.path(output_dir, "ot_dense_coverage.csv"),
          row.names = FALSE)

cat("dense Open Targets recompute complete ->", output_dir, "\n")
