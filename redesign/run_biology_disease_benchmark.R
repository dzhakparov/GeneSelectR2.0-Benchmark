#!/usr/bin/env Rscript
# ==============================================================================
#  Disease-specific biology assessment for exact top-50 benchmark panels.
#
#  Each panel is scored against a frozen Open Targets association list fixed in
#  analysis/config.R. The null contains 1000 random top-50 gene sets from
#  the same split-specific var2000 pool. This controls for dataset-specific gene
#  availability and uses the identical panel that was evaluated for prediction.
#
#  This analysis was specified after predictive results had been examined. It is
#  exploratory and must not be presented as independent validation.
#
#  Usage: Rscript redesign/run_biology_disease_benchmark.R <dataset> [budget_sec]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
ds_arg <- if (length(args) >= 1) args[1] else stop("dataset required")
budget <- if (length(args) >= 2) as.numeric(args[2]) else 270

validation_ds <- c("GSE65682", "GSE69683", "GSE13355", "GSE107994",
                   "GSE101794")
stopifnot(ds_arg %in% c(validation_ds, "imvigor210", "sosall"))

t_start <- proc.time()[["elapsed"]]
over_budget <- function() (proc.time()[["elapsed"]] - t_start) > budget

source(file.path("redesign", "R", "run_provenance.R"))

is_validation <- ds_arg %in% validation_ds
results_root <- redesign_results_root()
split_dir <- if (is_validation) {
  file.path(results_root, "validation_benchmark", ds_arg)
} else {
  file.path(results_root, "grouped_benchmark", ds_arg)
}
out_dir <- if (is_validation) split_dir else
  file.path(results_root, "full_recipe", ds_arg)
require_redesign_run(split_dir)
if (!identical(out_dir, split_dir)) require_redesign_run(out_dir)

extension_drivers <- c(
  predfirst = "redesign/run_predfirst_benchmark.R",
  ensemble = "redesign/run_ensemble_eval.R",
  softprior = "redesign/run_softprior_benchmark.R",
  datadriven = "redesign/run_datadriven_benchmark.R",
  horseshoe = "redesign/run_horseshoe_benchmark.R",
  ensemble2 = "redesign/run_ensemble2_eval.R",
  wsweep = "redesign/run_wsweep.R",
  hsstab = "redesign/run_hsstab_benchmark.R",
  prune = "redesign/run_prune_benchmark.R",
  adapt = "redesign/run_adapt_benchmark.R"
)
for (extension in names(extension_drivers)) {
  require_redesign_extension(
    out_dir, extension,
    redesign_extension_sources(extension_drivers[[extension]]),
    config = list(dataset = ds_arg)
  )
}

source(file.path("analysis", "config.R"))
config <- get_biology_config("disease")
config_row <- config[config$dataset == ds_arg, , drop = FALSE]
if (nrow(config_row) != 1L) {
  stop("Expected one disease-biology configuration for ", ds_arg)
}
association_path <- config_row$association_file[[1]]
if (!file.exists(association_path)) {
  stop("Missing frozen Open Targets association file: ", association_path)
}
association_hash <- unname(tools::md5sum(association_path))

out_path <- file.path(out_dir, "biology_disease.csv")
prepare_redesign_extension(
  out_dir, "biology_disease",
  c(redesign_extension_sources("redesign/run_biology_disease_benchmark.R"),
    config_path, association_path),
  config = list(
    dataset = ds_arg,
    disease_label = config_row$disease_label[[1]],
    ontology_id = config_row$ontology_id[[1]],
    association_md5 = association_hash,
    max_associations = 100L,
    minimum_association_score = 0.1,
    panel_size = 50L,
    n_null = 1000L
  ),
  output_files = out_path
)

associations <- head(
  read.delim(association_path, stringsAsFactors = FALSE), 100L
)
required_association_columns <- c("analyte", "score")
if (!all(required_association_columns %in% names(associations)) ||
    anyDuplicated(associations$analyte) ||
    any(!is.finite(associations$score)) ||
    any(associations$score < 0)) {
  stop("Invalid frozen Open Targets association file: ", association_path)
}
association_score <- setNames(associations$score, associations$analyte)

evaluation_paths <- sort(list.files(
  out_dir, pattern = "^eval.*[.]csv$", full.names = TRUE
))
evaluation_arms <- sort(unique(unlist(lapply(evaluation_paths, function(path) {
  tab <- read.csv(path, stringsAsFactors = FALSE)
  if (!"arm" %in% names(tab)) return(character(0))
  as.character(tab$arm)
}))))
arms <- setdiff(evaluation_arms, "Random")
if (length(arms) != 28L) {
  stop("Expected 28 non-random evaluated arms, found ", length(arms), ": ",
       paste(arms, collapse = ", "))
}

nm_grouped <- if (is_validation) "GS_full_grouped" else "full_grouped"
nm_ungrouped <- if (is_validation) "GS_full_ungrouped" else "full_ungrouped"
split_files <- sort(list.files(split_dir, pattern = "^split_r.*[.]rds$",
                               full.names = TRUE))
stopifnot(length(split_files) == 15L)

done <- if (file.exists(out_path)) read.csv(out_path) else data.frame()
done_key <- if (nrow(done) > 0) {
  paste(done$repeat_idx, done$fold_idx, done$arm)
} else {
  character(0)
}

for (split_path in split_files) {
  matched <- regmatches(
    basename(split_path),
    regexec("split_r(\\d+)_f(\\d+)[.]rds", basename(split_path))
  )[[1]]
  repeat_idx <- as.integer(matched[2])
  fold_idx <- as.integer(matched[3])
  if (!any(!(paste(repeat_idx, fold_idx, arms) %in% done_key))) next
  if (over_budget()) {
    cat("[budget] stop\n")
    quit(save = "no")
  }

  split <- readRDS(split_path)
  pool <- unique(split$pools$var2000)
  panel_size <- min(50L, length(pool))
  pool_scores <- association_score[pool]
  pool_scores[is.na(pool_scores)] <- 0
  n_targets_in_pool <- sum(pool_scores > 0)

  set.seed(900000L + repeat_idx * 100L + fold_idx)
  null_sums <- vapply(seq_len(1000L), function(index) {
    sum(pool_scores[sample.int(length(pool), panel_size, replace = FALSE)])
  }, numeric(1))
  expected_sum <- mean(null_sums)
  null_sd <- stats::sd(null_sums)

  rankings <- list()
  for (arm in arms) {
    ranking_path <- file.path(
      out_dir, sprintf("ranking_r%d_f%d_%s.csv", repeat_idx, fold_idx, arm)
    )
    if (!file.exists(ranking_path)) {
      stop("Missing exact ranking for disease biology: ", ranking_path)
    }
    ranking <- read.csv(ranking_path, stringsAsFactors = FALSE)$gene
    if (length(ranking) == 0 && arm == nm_grouped) {
      ranking <- read.csv(file.path(
        out_dir,
        sprintf("ranking_r%d_f%d_%s.csv", repeat_idx, fold_idx, nm_ungrouped)
      ), stringsAsFactors = FALSE)$gene
    }
    if (anyDuplicated(ranking)) {
      stop("Duplicate genes in ranking: ", ranking_path)
    }
    rankings[[arm]] <- ranking
  }

  rows <- list()
  for (arm in arms) {
    key <- paste(repeat_idx, fold_idx, arm)
    if (key %in% done_key) next
    panel <- unique(head(rankings[[arm]], panel_size))
    panel_scores <- association_score[panel]
    panel_scores[is.na(panel_scores)] <- 0
    observed_sum <- sum(panel_scores)
    n_overlap <- sum(panel_scores > 0)
    rows[[length(rows) + 1L]] <- data.frame(
      repeat_idx = repeat_idx,
      fold_idx = fold_idx,
      arm = arm,
      disease_label = config_row$disease_label[[1]],
      ontology_id = config_row$ontology_id[[1]],
      association_md5 = association_hash,
      panel_size = length(panel),
      n_targets_in_pool = n_targets_in_pool,
      n_target_overlap = n_overlap,
      target_recall = if (n_targets_in_pool > 0) n_overlap / n_targets_in_pool
                      else NA_real_,
      association_sum = observed_sum,
      expected_sum = expected_sum,
      enrichment_ratio = if (expected_sum > 0) observed_sum / expected_sum
                         else NA_real_,
      enrichment_z = if (is.finite(null_sd) && null_sd > 0) {
        (observed_sum - expected_sum) / null_sd
      } else {
        NA_real_
      },
      empirical_p = (1 + sum(null_sums >= observed_sum)) /
                    (length(null_sums) + 1),
      stringsAsFactors = FALSE
    )
    done_key <- c(done_key, key)
  }
  if (length(rows) > 0) {
    write.table(
      do.call(rbind, rows), out_path,
      append = file.exists(out_path), sep = ",", row.names = FALSE,
      col.names = !file.exists(out_path)
    )
  }
  cat(sprintf("[%s r%d f%d] disease biology done\n", ds_arg, repeat_idx,
              fold_idx))
}

cat(sprintf("\nDone dataset=%s (%.0f s)\n", ds_arg,
            proc.time()[["elapsed"]] - t_start))
