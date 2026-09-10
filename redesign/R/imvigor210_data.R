# ==============================================================================
#  IMvigor210 data loading and preprocessing, shared by the redesign runners.
#
#  COPIED VERBATIM from redesign/run_gs_slim_imvigor210.R (which itself copied
#  the incumbent benchmark) so that every redesign variant sees byte-identical
#  data. Do not "improve" anything in here: if the incumbent changes, change
#  this file and every runner that uses it in the same commit.
#
#  Provides:
#    load_imvigor210()          -> list(X, y, clinical)  (log-CPM NOT applied;
#                                 returns raw counts + outcome, full cohort)
#    normalise_count_split()    -> train-fitted TMM + log-CPM
#    preprocess_split()         -> + variance filter to cfg$top_genes +
#                                  standardisation, all fitted on train only
#    make_stratified_folds()    -> the incumbent's fold construction
#
#  Constants match the incumbent: min_count_per_gene = 10,
#  min_samples_per_gene = 10, random_seed = 42.
# ==============================================================================

imv_min_count_per_gene   <- 10
imv_min_samples_per_gene <- 10
imv_random_seed          <- 42

load_imvigor210 <- function() {
  raw_counts <- NULL
  clinical <- NULL

  if (requireNamespace("easierData", quietly = TRUE)) {
    ok <- tryCatch({
      suppressPackageStartupMessages({
        library(ExperimentHub); library(easierData)
        library(SummarizedExperiment)
      })
      dat <- ExperimentHub(
        cache = file.path("data", "experimenthub_cache"))[["EH6677"]]
      raw_counts <- as.matrix(SummarizedExperiment::assay(dat, "counts"))
      clinical <- as.data.frame(SummarizedExperiment::colData(dat))
      TRUE
    }, error = function(e) {
      cat(sprintf("  easierData failed: %s\n", conditionMessage(e)))
      FALSE
    })
  }

  if (is.null(raw_counts)) {
    rds_path <- file.path("data", "IMvigor210.all.rds")
    if (!file.exists(rds_path)) {
      stop("No IMvigor210 source available (easierData failed, no local RDS).")
    }
    imv <- readRDS(rds_path)
    slot <- intersect(c("rawcounts", "counts"), names(imv))[1]
    raw_counts <- as.matrix(imv[[slot]])
    clinical <- as.data.frame(imv$clinical)
  }

  response_col <- intersect(c("binaryResponse", "BOR",
                              "Best.Confirmed.Overall.Response"),
                            names(clinical))[1]
  if (is.na(response_col)) stop("No response column found in clinical data.")

  response_values <- as.character(clinical[[response_col]])
  evaluable <- response_values %in% c("R", "NR")
  raw_counts <- raw_counts[, evaluable, drop = FALSE]
  clinical <- clinical[evaluable, , drop = FALSE]

  outcome_factor <- factor(
    ifelse(as.character(clinical[[response_col]]) == "R",
           "Responder", "NonResponder"),
    levels = c("NonResponder", "Responder")
  )
  cat(sprintf("Samples: %d (%s)\n", length(outcome_factor),
              paste(table(outcome_factor), collapse = "/")))

  if (any(!is.finite(raw_counts)) || any(raw_counts < 0) ||
      any(abs(raw_counts - round(raw_counts)) > .Machine$double.eps^0.5)) {
    stop("Count assay must contain finite non-negative integers.")
  }
  gene_ids <- rownames(raw_counts)
  if (any(duplicated(gene_ids))) {
    raw_counts <- rowsum(raw_counts, group = gene_ids, reorder = FALSE)
  }

  list(raw_counts = raw_counts, outcome = outcome_factor,
       clinical = clinical)
}

normalise_count_split <- function(raw_count_matrix, train_indices,
                                  test_indices = NULL) {
  train_counts <- raw_count_matrix[, train_indices, drop = FALSE]
  keep <- rowSums(train_counts >= imv_min_count_per_gene) >=
    min(imv_min_samples_per_gene, length(train_indices))
  train_counts <- train_counts[keep, , drop = FALSE]
  test_counts <- if (is.null(test_indices)) NULL else
    raw_count_matrix[keep, test_indices, drop = FALSE]

  library_sizes <- colSums(train_counts)
  upper_quartile_cpm <- apply(train_counts, 2, function(counts) {
    stats::quantile(counts / sum(counts) * 1e6, 0.75,
                    names = FALSE, type = 7)
  })
  reference_index <- which.min(
    abs(upper_quartile_cpm - mean(upper_quartile_cpm)))

  train_dge <- edgeR::DGEList(counts = train_counts)
  train_dge <- edgeR::calcNormFactors(train_dge,
                                      refColumn = reference_index)
  train_factors <- train_dge$samples$norm.factors
  train_logcpm <- edgeR::cpm(train_dge, log = TRUE, prior.count = 1)

  if (is.null(test_counts)) return(list(train = t(train_logcpm), test = NULL))

  test_logcpm <- vapply(seq_len(ncol(test_counts)), function(j) {
    pair_dge <- edgeR::DGEList(counts = cbind(
      reference = train_counts[, reference_index],
      held_out = test_counts[, j]))
    pair_dge <- edgeR::calcNormFactors(pair_dge, refColumn = 1L)
    pair_ratio <- pair_dge$samples$norm.factors[2] /
      pair_dge$samples$norm.factors[1]
    held_out_factor <- train_factors[reference_index] * pair_ratio
    pair_dge$samples$norm.factors <- c(train_factors[reference_index],
                                       held_out_factor)
    edgeR::cpm(pair_dge, log = TRUE, prior.count = 1)[, 2]
  }, numeric(nrow(train_counts)))
  rownames(test_logcpm) <- rownames(train_counts)

  list(train = t(train_logcpm), test = t(test_logcpm))
}

preprocess_split <- function(raw_count_matrix, train_indices, test_indices,
                             top_genes) {
  normalized <- normalise_count_split(raw_count_matrix, train_indices,
                                      test_indices)
  train_expression <- normalized$train
  test_expression <- normalized$test

  training_variances <- apply(train_expression, 2, var)
  eligible <- which(is.finite(training_variances) & training_variances > 0)
  if (length(eligible) > top_genes) {
    eligible <- eligible[
      order(training_variances[eligible], decreasing = TRUE)[
        seq_len(top_genes)]]
  }
  train_expression <- train_expression[, eligible, drop = FALSE]
  test_expression <- test_expression[, eligible, drop = FALSE]

  column_means <- colMeans(train_expression)
  column_sds <- apply(train_expression, 2, sd)
  column_sds[column_sds == 0 | is.na(column_sds)] <- 1
  train_expression <- sweep(sweep(train_expression, 2, column_means, "-"),
                            2, column_sds, "/")
  test_expression <- sweep(sweep(test_expression, 2, column_means, "-"),
                           2, column_sds, "/")
  list(train = train_expression, test = test_expression)
}

make_stratified_folds <- function(outcome_factor, k_folds, seed) {
  set.seed(seed)
  outcome_factor <- droplevels(outcome_factor)
  class1_indices <- which(outcome_factor == levels(outcome_factor)[1])
  class2_indices <- which(outcome_factor == levels(outcome_factor)[2])
  class1_folds <- split(sample(class1_indices),
                        rep(1:k_folds, length.out = length(class1_indices)))
  class2_folds <- split(sample(class2_indices),
                        rep(1:k_folds, length.out = length(class2_indices)))
  lapply(1:k_folds, function(fold) {
    sort(c(class1_folds[[fold]], class2_folds[[fold]]))
  })
}
