#!/usr/bin/env Rscript
# ==============================================================================
#  Greedy forward panel selection - IMvigor210 benchmark runner
# ==============================================================================
#
#  Same protocol as the GS_slim pilot and the incumbent benchmark: same outer
#  5-fold split (repeat 1), same train-fitted TMM + log-CPM + variance filter
#  + standardisation, same top-k panel sizes, same three-model ensemble
#  evaluator. Only the ranking method differs: the panel is picked directly by
#  greedy forward selection on inner-CV AUC.
#
#  Usage (from the project root):
#    Rscript redesign/run_greedy_imvigor210.R          # all 5 folds
#    Rscript redesign/run_greedy_imvigor210.R 2        # fold 2 only
#    Rscript redesign/run_greedy_imvigor210.R 2-4      # folds 2 to 4
#
#  Each fold writes its own CSV; finished folds are skipped on re-run.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
fold_range <- if (length(args) >= 1 && grepl("^[0-9]+(-[0-9]+)?$", args[1])) {
  parts <- as.integer(strsplit(args[1], "-")[[1]])
  if (length(parts) == 1) parts else seq(parts[1], parts[2])
} else 1:5

random_seed <- 42
panel_sizes <- c(10, 20, 50, 100, 200, 500)
glmnet_alpha_grid <- c(0.5, 1.0)
min_count_per_gene   <- 10
min_samples_per_gene <- 10
top_genes <- 2000

# Greedy settings.
n_candidates <- 200
k_max        <- 50
inner_k      <- 5

out_dir <- file.path("redesign", "results", "greedy_imvigor210")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
incumbent_csv <- file.path("results_imvigor210", "2026-08-08_kfold",
                           "data", "nested_results.csv")

n_cores <- max(1L, min(7L, parallel::detectCores(logical = FALSE) - 1L))
cat(sprintf("Greedy panel | folds %s | p=%d pool=%d k_max=%d | cores=%d\n\n",
            paste(fold_range, collapse = ","), top_genes, n_candidates,
            k_max, n_cores))

suppressPackageStartupMessages({library(glmnet); library(edgeR)})
for (f in c("gs_slim.R", "greedy_panel.R")) source(file.path("redesign", "R", f))


# ==============================================================================
#  Data loading + preprocessing (identical to the other runners)
# ==============================================================================

raw_counts <- NULL
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
  }, error = function(e) FALSE)
}
if (is.null(raw_counts)) {
  rds_path <- file.path("data", "IMvigor210.all.rds")
  if (!file.exists(rds_path)) stop("No IMvigor210 source available.")
  imv <- readRDS(rds_path)
  raw_counts <- as.matrix(imv[[intersect(c("rawcounts", "counts"),
                                         names(imv))[1]]])
  clinical <- as.data.frame(imv$clinical)
}

response_col <- intersect(c("binaryResponse", "BOR",
                            "Best.Confirmed.Overall.Response"),
                          names(clinical))[1]
response_values <- as.character(clinical[[response_col]])
evaluable <- response_values %in% c("R", "NR")
if (sum(evaluable) == 0) stop("No R/NR samples found.")
raw_counts <- raw_counts[, evaluable, drop = FALSE]
clinical <- clinical[evaluable, , drop = FALSE]
outcome_factor <- factor(
  ifelse(as.character(clinical[[response_col]]) == "R",
         "Responder", "NonResponder"),
  levels = c("NonResponder", "Responder"))
cat(sprintf("Samples: %d (%s)\n", length(outcome_factor),
            paste(table(outcome_factor), collapse = "/")))
gene_ids <- rownames(raw_counts)
if (any(duplicated(gene_ids))) {
  raw_counts <- rowsum(raw_counts, group = gene_ids, reorder = FALSE)
}
raw_count_matrix <- raw_counts

normalise_count_split <- function(train_indices, test_indices = NULL) {
  train_counts <- raw_count_matrix[, train_indices, drop = FALSE]
  keep <- rowSums(train_counts >= min_count_per_gene) >=
    min(min_samples_per_gene, length(train_indices))
  train_counts <- train_counts[keep, , drop = FALSE]
  test_counts <- if (is.null(test_indices)) NULL else
    raw_count_matrix[keep, test_indices, drop = FALSE]
  upper_quartile_cpm <- apply(train_counts, 2, function(counts) {
    stats::quantile(counts / sum(counts) * 1e6, 0.75, names = FALSE, type = 7)
  })
  reference_index <- which.min(
    abs(upper_quartile_cpm - mean(upper_quartile_cpm)))
  train_dge <- edgeR::DGEList(counts = train_counts)
  train_dge <- edgeR::calcNormFactors(train_dge, refColumn = reference_index)
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
    pair_dge$samples$norm.factors <- c(
      train_factors[reference_index],
      train_factors[reference_index] * pair_ratio)
    edgeR::cpm(pair_dge, log = TRUE, prior.count = 1)[, 2]
  }, numeric(nrow(train_counts)))
  rownames(test_logcpm) <- rownames(train_counts)
  list(train = t(train_logcpm), test = t(test_logcpm))
}

preprocess_split <- function(train_indices, test_indices) {
  normalized <- normalise_count_split(train_indices, test_indices)
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
  list(train = sweep(sweep(train_expression, 2, column_means, "-"),
                     2, column_sds, "/"),
       test  = sweep(sweep(test_expression, 2, column_means, "-"),
                     2, column_sds, "/"))
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
  lapply(1:k_folds, function(fold) sort(c(class1_folds[[fold]],
                                          class2_folds[[fold]])))
}


# ==============================================================================
#  Evaluator: the incumbent's three-model ensemble (identical)
# ==============================================================================

bench_auc <- function(true_labels, predicted_scores) {
  true_labels <- droplevels(true_labels)
  pos <- which(true_labels == levels(true_labels)[2])
  neg <- which(true_labels == levels(true_labels)[1])
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(predicted_scores, ties.method = "average")
  (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

compute_classification_metrics <- function(true_labels, probs, threshold = 0.5) {
  true_labels <- droplevels(true_labels)
  positive_class <- levels(true_labels)[2]
  negative_class <- levels(true_labels)[1]
  predicted_labels <- factor(ifelse(probs >= threshold,
                                    positive_class, negative_class),
                             levels = levels(true_labels))
  tp <- sum(predicted_labels == positive_class & true_labels == positive_class)
  tn <- sum(predicted_labels == negative_class & true_labels == negative_class)
  fp <- sum(predicted_labels == positive_class & true_labels == negative_class)
  fn <- sum(predicted_labels == negative_class & true_labels == positive_class)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else 0
  denom <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  list(balanced_accuracy = (sensitivity + specificity) / 2,
       mcc = if (denom > 0) (tp * tn - fp * fn) / denom else 0)
}

predict_with_ensemble <- function(train_features, train_labels, test_features) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  n_inner <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))
  best_fit <- NULL; best_auc <- -Inf
  for (alpha_value in glmnet_alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(train_features, numeric_labels,
                        family = "binomial", alpha = alpha_value,
                        nfolds = n_inner),
      error = function(e) NULL)
    if (!is.null(fit) && -min(fit$cvm) > best_auc) {
      best_auc <- -min(fit$cvm); best_fit <- fit
    }
  }
  p_glmnet <- if (is.null(best_fit)) rep(NA_real_, nrow(test_features)) else
    as.numeric(predict(best_fit, test_features,
                       s = "lambda.min", type = "response"))
  p_xgb <- tryCatch({
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_features),
                                   label = numeric_labels)
    fit <- xgboost::xgb.train(
      params = list(objective = "binary:logistic", eval_metric = "auc",
                    eta = 0.1, max_depth = 3, subsample = 0.8,
                    colsample_bytree = 0.8),
      data = dtrain, nrounds = 100, verbose = 0)
    predict(fit, as.matrix(test_features))
  }, error = function(e) rep(NA_real_, nrow(test_features)))
  p_rf <- tryCatch({
    fit <- ranger::ranger(x = train_features, y = droplevels(train_labels),
                          num.trees = 500, probability = TRUE,
                          seed = random_seed, num.threads = 1)
    predict(fit, data = test_features)$predictions[,
      levels(droplevels(train_labels))[2]]
  }, error = function(e) rep(NA_real_, nrow(test_features)))
  components <- cbind(p_glmnet, p_xgb, p_rf)
  if (any(!is.finite(components))) {
    stop("An ensemble component failed or returned non-finite predictions.",
         call. = FALSE)
  }
  ensemble <- rowMeans(components)
  list(glmnet = p_glmnet, xgboost = p_xgb, rf = p_rf, ensemble = ensemble)
}


# ==============================================================================
#  Main loop
# ==============================================================================

folds <- make_stratified_folds(outcome_factor, 5, seed = random_seed + 1000)

for (fold_idx in fold_range) {
  split_file <- file.path(out_dir, sprintf("split_f%d.csv", fold_idx))
  if (file.exists(split_file)) {
    cat(sprintf("[f%d] already done, skipping\n", fold_idx)); next
  }

  cat(sprintf("\n[f%d] preprocessing + greedy selection...\n", fold_idx))
  t0 <- proc.time()
  test_indices <- folds[[fold_idx]]
  train_indices <- setdiff(seq_along(outcome_factor), test_indices)
  train_labels <- outcome_factor[train_indices]
  test_labels <- outcome_factor[test_indices]
  split_data <- preprocess_split(train_indices, test_indices)

  greedy <- greedy_forward_panel(
    split_data$train, train_labels,
    n_candidates = n_candidates, k_max = k_max, inner_k = inner_k,
    n_cores = n_cores, random_seed = random_seed + fold_idx
  )

  # Full ranking: greedy path first, then the rest of the candidate pool in
  # univariate order, then everything else (for k > pool size).
  uni_rank <- univariate_auc_rank(split_data$train, train_labels)
  ranked <- c(greedy$path,
              setdiff(uni_rank, greedy$path),
              setdiff(colnames(split_data$train),
                      c(greedy$path, uni_rank)))

  write.csv(data.frame(rank = seq_along(ranked), gene = ranked,
                       in_greedy_path = ranked %in% greedy$path),
            file.path(out_dir, sprintf("ranking_f%d.csv", fold_idx)),
            row.names = FALSE)
  write.csv(data.frame(step = seq_along(greedy$step_auc),
                       gene = greedy$path, inner_cv_auc = greedy$step_auc),
            file.path(out_dir, sprintf("path_f%d.csv", fold_idx)),
            row.names = FALSE)

  elapsed <- (proc.time() - t0)["elapsed"]

  rows <- list()
  for (k in panel_sizes) {
    panel <- head(ranked, k)
    probs <- predict_with_ensemble(
      split_data$train[, panel, drop = FALSE], train_labels,
      split_data$test[, panel, drop = FALSE])
    for (evaluator in names(probs)) {
      p <- probs[[evaluator]]
      auc <- if (all(is.na(p))) NA_real_ else bench_auc(test_labels, p)
      met <- if (all(is.na(p))) list(balanced_accuracy = NA_real_,
                                     mcc = NA_real_)
      else compute_classification_metrics(test_labels, p)
      rows[[length(rows) + 1]] <- data.frame(
        Method = "Greedy", Evaluator = evaluator, Repeat = 1, Fold = fold_idx,
        k = k, AUC = auc, BalAcc = met$balanced_accuracy, MCC = met$mcc,
        Time = elapsed, stringsAsFactors = FALSE)
    }
  }
  write.csv(do.call(rbind, rows), split_file, row.names = FALSE)
  cat(sprintf("[f%d] done in %.1f min\n", fold_idx, elapsed / 60))
}


# ==============================================================================
#  Comparison (runs when all requested folds are present)
# ==============================================================================

split_files <- list.files(out_dir, pattern = "^split_.*\\.csv$",
                          full.names = TRUE)
greedy_rows <- do.call(rbind, lapply(split_files, read.csv))

if (!is.null(greedy_rows) && nrow(greedy_rows) > 0 &&
    file.exists(incumbent_csv)) {
  incumbent <- read.csv(incumbent_csv)
  incumbent <- incumbent[incumbent$Repeat == 1, ]
  slim_csv <- file.path("redesign", "results", "imvigor210_pilot",
                        "nested_results_gs_slim.csv")
  if (file.exists(slim_csv)) {
    incumbent <- rbind(incumbent, read.csv(slim_csv))
  }

  ens <- rbind(incumbent, greedy_rows)
  ens <- ens[ens$Evaluator == "ensemble", ]
  ks <- sort(unique(ens$k))
  random_auc <- tapply(ens$AUC[ens$Method == "Random"],
                       ens$k[ens$Method == "Random"], mean, na.rm = TRUE)
  methods <- c("Greedy", "GS_slim", "GS_harmonic", "GS_no_stability",
               "GS_semantic", "GS_utility_only", "DGE", "Boruta",
               "RF_importance", "LASSO", "ElasticNet", "Random")
  comparison <- do.call(rbind, lapply(methods, function(m) {
    data.frame(Method = m, k = ks,
               mean_AUC = sapply(ks, function(k) {
                 v <- ens$AUC[ens$Method == m & ens$k == k]
                 if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
               }))
  }))
  comparison$AUC_minus_Random <- comparison$mean_AUC -
    random_auc[as.character(comparison$k)]
  write.csv(comparison, file.path(out_dir, "comparison_vs_incumbents.csv"),
            row.names = FALSE)

  cat("\n=== AUC minus Random (ensemble, repeat 1) ===\n")
  print(comparison[comparison$Method %in% c("Greedy", "GS_harmonic",
                                            "GS_no_stability",
                                            "GS_utility_only", "DGE",
                                            "RF_importance", "Random"), ],
        row.names = FALSE, digits = 3)
}
