#!/usr/bin/env Rscript
# ==============================================================================
#  GeneSelectR Slim - IMvigor210 benchmark runner
# ==============================================================================
#
#  Evaluates GS_slim on IMvigor210 under the SAME protocol as the incumbent
#  benchmark (imvigor210_benchmark_v2.0.R): same outer CV, same train-fitted
#  TMM + log-CPM + variance filter + standardisation, same top-k panel sizes,
#  same three-model ensemble evaluator. The only thing that differs is the
#  ranking method. Incumbent numbers are read from the saved results of the
#  2026-08-08 kfold run, not re-run.
#
#  Usage (from the project root, in a terminal - see AGENTS.md on fork safety):
#
#    Rscript redesign/run_gs_slim_imvigor210.R smoke   # ~5 min:  1 fold,
#                                                      #   500 genes, tiny null
#    Rscript redesign/run_gs_slim_imvigor210.R pilot   # ~1 h:    1 repeat x
#                                                      #   5 folds, p=2000
#    Rscript redesign/run_gs_slim_imvigor210.R full    # ~3 h:    3 repeats x
#                                                      #   5 folds, M_null=100
#
#  Each split writes its own CSV, so an interrupted run RESUMES where it
#  stopped; re-running the same mode never redoes a finished split.
#
#  This script lives in redesign/ and touches nothing outside redesign/results/.
# ==============================================================================

mode <- {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[1] %in% c("smoke", "pilot", "full")) args[1]
  else "smoke"
}

cfg <- switch(mode,
  # smoke builds the same 5-fold split as the other modes but RUNS only fold
  # 1 (see the jobs loop). k_folds here is the split construction, not the
  # number of folds executed.
  smoke = list(n_repeats = 1, k_folds = 5, top_genes = 500,  B = 20,
               M_null = 30,  B_null = 10),
  pilot = list(n_repeats = 1, k_folds = 5, top_genes = 2000, B = 50,
               M_null = 50,  B_null = 20),
  full  = list(n_repeats = 3, k_folds = 5, top_genes = 2000, B = 50,
               M_null = 100, B_null = 20)
)

random_seed <- 42
panel_sizes <- c(10, 20, 50, 100, 200, 500)
glmnet_alpha_grid <- c(0.5, 1.0)     # evaluator's grid, as in the incumbent
alpha_grid <- c(0.5, 1.0)            # selector's grid, as in the incumbent
min_count_per_gene   <- 10
min_samples_per_gene <- 10

out_dir <- file.path("redesign", "results", paste0("imvigor210_", mode))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

incumbent_csv <- file.path("results_imvigor210", "2026-08-08_kfold",
                           "data", "nested_results.csv")

# Worker cap: 7 hard, per the project's worker budget. The budget lockfile
# machinery is not used here to keep redesign/ self-contained; do not run this
# concurrently with another benchmark on the same machine.
n_cores <- max(1L, min(7L, parallel::detectCores(logical = FALSE) - 1L))

cat(sprintf("GS_slim IMvigor210 | mode=%s | repeats=%d folds=%d p=%d | ",
            mode, cfg$n_repeats, cfg$k_folds, cfg$top_genes))
cat(sprintf("B=%d null=%dx%d | cores=%d\n\n", cfg$B, cfg$M_null, cfg$B_null,
            n_cores))

suppressPackageStartupMessages({
  library(glmnet)
  library(edgeR)
})

for (f in c("utility_deviance.R", "joint_null.R", "gs_slim.R")) {
  source(file.path("redesign", "R", f))
}


# ==============================================================================
#  Data loading (identical to the incumbent benchmark)
# ==============================================================================

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
raw_count_matrix <- raw_counts


# ==============================================================================
#  Preprocessing (identical to the incumbent benchmark)
# ==============================================================================
#
#  Count filter, TMM reference, variance filter and standardisation are all
#  fitted on the training fold only. This is verbatim the incumbent's
#  construction, so preprocessing cannot be the source of any difference in
#  results.

normalise_count_split <- function(train_indices, test_indices = NULL) {
  train_counts <- raw_count_matrix[, train_indices, drop = FALSE]
  keep <- rowSums(train_counts >= min_count_per_gene) >=
    min(min_samples_per_gene, length(train_indices))
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

preprocess_split <- function(train_indices, test_indices) {
  normalized <- normalise_count_split(train_indices, test_indices)
  train_expression <- normalized$train
  test_expression <- normalized$test

  training_variances <- apply(train_expression, 2, var)
  eligible <- which(is.finite(training_variances) & training_variances > 0)
  if (length(eligible) > cfg$top_genes) {
    eligible <- eligible[
      order(training_variances[eligible], decreasing = TRUE)[
        seq_len(cfg$top_genes)]]
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


# ==============================================================================
#  Evaluator: the incumbent's three-model soft-voting ensemble
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

  # Component 1: elastic net over the same alpha grid as the incumbent.
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

  # Component 2: XGBoost, incumbent defaults.
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

  # Component 3: random forest, incumbent defaults.
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
#  Main loop: one CSV per split, resumable
# ==============================================================================

jobs <- list()
for (repeat_idx in seq_len(cfg$n_repeats)) {
  folds <- make_stratified_folds(outcome_factor, cfg$k_folds,
                                 seed = random_seed + 1000 * repeat_idx)
  # smoke runs a single fold; pilot/full run all k_folds of each repeat.
  n_folds_here <- if (mode == "smoke") 1L else cfg$k_folds
  for (fold_idx in seq_len(n_folds_here)) {
    jobs[[length(jobs) + 1]] <- list(repeat_idx = repeat_idx,
                                     fold_idx = fold_idx,
                                     test_indices = folds[[fold_idx]])
  }
}

for (job in jobs) {
  split_file <- file.path(
    out_dir, sprintf("split_r%d_f%d.csv", job$repeat_idx, job$fold_idx))
  if (file.exists(split_file)) {
    cat(sprintf("[r%d f%d] already done, skipping\n",
                job$repeat_idx, job$fold_idx))
    next
  }

  cat(sprintf("\n[r%d f%d] preprocessing + GS_slim fit...\n",
              job$repeat_idx, job$fold_idx))
  t0 <- proc.time()
  train_indices <- setdiff(seq_along(outcome_factor), job$test_indices)
  train_labels <- outcome_factor[train_indices]
  test_labels <- outcome_factor[job$test_indices]
  split_data <- preprocess_split(train_indices, job$test_indices)

  fit <- gs_slim_fit(
    split_data$train, train_labels,
    alpha_grid = alpha_grid,
    B = cfg$B, k_folds = 5,
    M_null = cfg$M_null, B_null = cfg$B_null,
    n_cores = n_cores, random_seed = random_seed, verbose = TRUE
  )

  # Keep the per-gene table: this is the auditable output of the method.
  write.csv(fit$gene_table,
            file.path(out_dir, sprintf("genes_r%d_f%d.csv",
                                       job$repeat_idx, job$fold_idx)),
            row.names = FALSE)

  ranked <- fit$ranked
  elapsed <- (proc.time() - t0)["elapsed"]

  rows <- list()
  for (k in panel_sizes) {
    panel <- head(ranked[ranked %in% colnames(split_data$train)], k)
    if (length(panel) < 5L) next
    probs <- predict_with_ensemble(
      split_data$train[, panel, drop = FALSE], train_labels,
      split_data$test[, panel, drop = FALSE])
    for (evaluator in names(probs)) {
      p <- probs[[evaluator]]
      auc <- if (all(is.na(p))) NA_real_ else bench_auc(test_labels, p)
      met <- if (all(is.na(p))) list(balanced_accuracy = NA_real_, mcc = NA_real_)
      else compute_classification_metrics(test_labels, p)
      rows[[length(rows) + 1]] <- data.frame(
        Method = "GS_slim", Evaluator = evaluator,
        Repeat = job$repeat_idx, Fold = job$fold_idx, k = k,
        AUC = auc, BalAcc = met$balanced_accuracy, MCC = met$mcc,
        Time = elapsed, stringsAsFactors = FALSE)
    }
  }
  write.csv(do.call(rbind, rows), split_file, row.names = FALSE)
  cat(sprintf("[r%d f%d] done in %.1f min\n",
              job$repeat_idx, job$fold_idx, elapsed / 60))
}


# ==============================================================================
#  Comparison against incumbents (from their saved results, not re-run)
# ==============================================================================

slim_rows <- do.call(rbind, lapply(
  list.files(out_dir, pattern = "^split_.*\\.csv$", full.names = TRUE),
  read.csv))

if (nrow(slim_rows) == 0) {
  cat("\nNo split results found; nothing to compare yet.\n")
  quit(status = 0, save = "no")
}

write.csv(slim_rows, file.path(out_dir, "nested_results_gs_slim.csv"),
          row.names = FALSE)

# AUC-minus-Random at matched k is the primary metric (AGENTS.md: raw AUC is
# misleading because random signatures inherit dominant meta-genes).
summarise_auc <- function(df) {
  ens <- df[df$Evaluator == "ensemble", ]
  random_auc <- tapply(ens$AUC[ens$Method == "Random"],
                       ens$k[ens$Method == "Random"], mean, na.rm = TRUE)
  methods <- sort(unique(ens$Method))
  ks <- sort(unique(ens$k))
  out <- expand.grid(Method = methods, k = ks, stringsAsFactors = FALSE)
  out$mean_AUC <- mapply(function(m, k) {
    v <- ens$AUC[ens$Method == m & ens$k == k]
    if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
  }, out$Method, out$k)
  out$AUC_minus_Random <- out$mean_AUC -
    random_auc[as.character(out$k)]
  out
}

if (file.exists(incumbent_csv)) {
  incumbent <- read.csv(incumbent_csv)
  # Same-protocol subset: incumbents ran 3 repeats; when this run is pilot or
  # smoke, the fair comparison is repeat 1 only.
  if (cfg$n_repeats == 1) incumbent <- incumbent[incumbent$Repeat == 1, ]
  if (mode == "smoke") incumbent <- incumbent[incumbent$Fold == 1, ]

  comparison <- summarise_auc(rbind(incumbent, slim_rows))
  comparison <- comparison[order(comparison$k,
                                 -comparison$AUC_minus_Random), ]
  write.csv(comparison, file.path(out_dir, "comparison_vs_incumbents.csv"),
            row.names = FALSE)

  cat("\n=== AUC minus Random (ensemble evaluator) ===\n")
  show <- comparison[comparison$Method %in%
                       c("GS_slim", "GS_harmonic", "GS_no_stability",
                         "GS_semantic", "DGE", "Boruta", "RF_importance",
                         "Random"), ]
  print(show, row.names = FALSE, digits = 3)
  cat("\nResults written to", out_dir, "\n")
} else {
  cat(sprintf("\nIncumbent results not found at %s; GS_slim results saved ",
              incumbent_csv))
  cat("without comparison.\n")
}
