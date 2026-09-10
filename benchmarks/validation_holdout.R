# ==============================================================================
#  Locked-panel validation in an external cohort
# ==============================================================================
#
#  Feature panels are fixed in a source cohort by lock_validation_panels.R.
#  This script evaluates those fixed genes in repeated CV within the external
#  cohort. Feature selection is never repeated in the holdout. The result tests
#  panel portability across cohorts and platforms. It is distinct from direct
#  transport of a fitted source classifier, which is unsuitable here because
#  GSE107994 is RNA-seq and GSE19442 is an Illumina microarray.
#
#  Usage:
#      Rscript benchmarks/validation_holdout.R \
#        GSE19442 locked_panels/GSE107994_to_GSE19442.csv
# ==============================================================================


# ------------------------------------------------------------------------------
#  Configuration
# ------------------------------------------------------------------------------

CONFIG <- list(
  random_seed          = 42L,
  outer_folds          = 5L,
  outer_repeats        = 3L,
  glmnet_alpha_grid    = c(0.5, 1.0),
  minimum_gene_overlap = 0.80,
  output_root          = "results_validation"
)

# Warnings from small-fold model fits must appear beside the affected run.
# Delayed warnings lose the method and split context when hundreds of models run.
options(warn = 1)


`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
source("benchmarks/validation_datasets.R")


arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2) {
  stop(paste(
    "Usage: Rscript benchmarks/validation_holdout.R",
    "<HOLDOUT_ACCESSION> <LOCKED_PANELS_CSV>"
  ))
}

accession <- arguments[1]
panel_file <- arguments[2]
dataset <- get_validation_dataset(accession)

if (dataset$analysis_role != "external_holdout") {
  stop(accession, " is not registered as an external holdout.")
}
if (!file.exists(panel_file)) stop("Locked panel file is missing: ", panel_file)

required_packages <- c("glmnet", "ranger", "withr", "xgboost")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop("Required packages are missing: ", paste(missing_packages, collapse = ", "))
}


# ------------------------------------------------------------------------------
#  Inputs and lock checks
# ------------------------------------------------------------------------------

data_dir <- file.path("data", accession)
expression_file <- file.path(data_dir, "expression_prepared.csv")
metadata_file <- file.path(data_dir, "metadata_prepared.csv")

if (!file.exists(expression_file) || !file.exists(metadata_file)) {
  stop(sprintf(
    "%s is not prepared. Run Rscript benchmarks/validation_prepare.R %s",
    accession, accession
  ))
}

panels <- utils::read.csv(panel_file, stringsAsFactors = FALSE,
                          check.names = FALSE)
required_panel_columns <- c(
  "source_accession", "holdout_accession", "source_rankings_md5",
  "holdout_platform", "holdout_annotation_md5", "locked_at_utc",
  "aggregation_rule", "method", "panel_size", "rank", "gene"
)
missing_panel_columns <- setdiff(required_panel_columns, colnames(panels))
if (length(missing_panel_columns) > 0) {
  stop("Locked panel file is missing: ",
       paste(missing_panel_columns, collapse = ", "))
}
if (!identical(unique(panels$source_accession), dataset$source_accession)) {
  stop("Locked panels were not generated from ", dataset$source_accession, ".")
}
if (!identical(unique(panels$holdout_accession), accession)) {
  stop("Locked panel holdout accession does not match ", accession, ".")
}
if (anyDuplicated(panels[c("method", "panel_size", "gene")])) {
  stop("A locked method/panel contains duplicated genes.")
}

lock_times <- as.POSIXct(unique(panels$locked_at_utc),
                         format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
if (length(lock_times) != 1 || is.na(lock_times)) {
  stop("locked_at_utc must contain one valid UTC timestamp.")
}
if (lock_times > file.info(metadata_file)$mtime) {
  stop(paste(
    "Panels were locked after the prepared holdout metadata was created.",
    "Recreate the external-validation sequence with panels locked first."
  ))
}

panel_counts <- aggregate(gene ~ method + panel_size, panels, length)
if (any(panel_counts$gene != panel_counts$panel_size)) {
  stop("Each method/panel_size must contain exactly panel_size genes.")
}

expression <- as.matrix(utils::read.csv(
  expression_file, row.names = 1, check.names = FALSE
))
storage.mode(expression) <- "numeric"
metadata <- utils::read.csv(metadata_file, stringsAsFactors = FALSE,
                            check.names = FALSE)

required_metadata <- c("sample_id", "outcome", "group_id")
if (length(setdiff(required_metadata, colnames(metadata))) > 0) {
  stop("Prepared metadata has an invalid schema.")
}
if (!setequal(colnames(expression), metadata$sample_id)) {
  stop("Expression columns and metadata sample IDs differ.")
}
metadata <- metadata[match(colnames(expression), metadata$sample_id), , drop = FALSE]
labels <- factor(metadata$outcome,
                 levels = c(dataset$outcome$negative_label,
                            dataset$outcome$positive_label))
if (anyNA(labels)) stop("Prepared metadata contains an unexpected outcome label.")
features <- t(expression)


# ------------------------------------------------------------------------------
#  Stratified folds
# ------------------------------------------------------------------------------

make_stratified_folds <- function(labels, k, repeats, seed) {
  set.seed(seed)
  folds <- list()
  for (repeat_index in seq_len(repeats)) {
    assignment <- integer(length(labels))
    for (class_label in levels(labels)) {
      indices <- which(labels == class_label)
      indices <- sample(indices, length(indices))
      assignment[indices] <- rep(seq_len(k), length.out = length(indices))
    }
    for (fold_index in seq_len(k)) {
      folds[[length(folds) + 1L]] <- list(
        repeat_index = repeat_index,
        fold = fold_index,
        train = which(assignment != fold_index),
        test = which(assignment == fold_index)
      )
    }
  }
  folds
}


# ------------------------------------------------------------------------------
#  Evaluators
# ------------------------------------------------------------------------------

predict_glmnet <- function(train_x, train_y, test_x, random_seed) {
  numeric_y <- as.integer(train_y == levels(train_y)[2])
  # GSE19442 outer-training sets contain about 40 samples. Four inner folds
  # retain at least 10 observations per fold, which glmnet requires for AUC.
  inner_folds <- min(4L, max(3L, floor(nrow(train_x) / 10)))
  fold_id <- withr::with_seed(random_seed, {
    result <- integer(length(train_y))
    for (class_label in levels(train_y)) {
      indices <- which(train_y == class_label)
      result[indices] <- sample(rep(seq_len(inner_folds),
                                    length.out = length(indices)))
    }
    result
  })
  best_fit <- NULL
  best_auc <- -Inf

  for (alpha_value in CONFIG$glmnet_alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(
        train_x, numeric_y, family = "binomial", alpha = alpha_value,
        foldid = fold_id, type.measure = "auc"
      ),
      error = function(error) NULL
    )
    if (!is.null(fit) && max(fit$cvm, na.rm = TRUE) > best_auc) {
      best_auc <- max(fit$cvm, na.rm = TRUE)
      best_fit <- fit
    }
  }
  if (is.null(best_fit)) return(rep(NA_real_, nrow(test_x)))
  as.numeric(stats::predict(best_fit, test_x,
                            s = "lambda.min", type = "response"))
}


predict_xgboost <- function(train_x, train_y, test_x, random_seed) {
  numeric_y <- as.integer(train_y == levels(train_y)[2])
  train_matrix <- xgboost::xgb.DMatrix(as.matrix(train_x), label = numeric_y)
  fit <- tryCatch(
    withr::with_seed(random_seed + 1000L, xgboost::xgb.train(
      params = list(
        objective = "binary:logistic", eval_metric = "auc", eta = 0.1,
        max_depth = 3L, subsample = 0.8, colsample_bytree = 0.8,
        seed = random_seed + 1000L, nthread = 1L
      ),
      data = train_matrix, nrounds = 100L, verbose = 0
    )),
    error = function(error) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(test_x)))
  tryCatch(stats::predict(fit, as.matrix(test_x)),
           error = function(error) rep(NA_real_, nrow(test_x)))
}


predict_ranger <- function(train_x, train_y, test_x, random_seed) {
  positive <- levels(train_y)[2]
  fit <- tryCatch(
    ranger::ranger(
      x = data.frame(train_x, check.names = FALSE), y = train_y,
      num.trees = 500L, probability = TRUE,
      seed = random_seed + 2000L, num.threads = 1L
    ),
    error = function(error) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, nrow(test_x)))
  tryCatch(
    predict(fit, data = data.frame(test_x, check.names = FALSE))$predictions[, positive],
    error = function(error) rep(NA_real_, nrow(test_x))
  )
}


predict_ensemble <- function(train_x, train_y, test_x, random_seed) {
  components <- cbind(
    glmnet = predict_glmnet(train_x, train_y, test_x, random_seed),
    xgboost = predict_xgboost(train_x, train_y, test_x, random_seed),
    random_forest = predict_ranger(train_x, train_y, test_x, random_seed)
  )
  ensemble <- rowMeans(components, na.rm = TRUE)
  ensemble[is.nan(ensemble)] <- NA_real_
  cbind(ensemble = ensemble, components)
}


# ------------------------------------------------------------------------------
#  Metrics
# ------------------------------------------------------------------------------

compute_auc <- function(labels, probabilities) {
  keep <- is.finite(probabilities)
  y <- as.integer(labels[keep] == levels(labels)[2])
  scores <- probabilities[keep]
  positives <- sum(y == 1L)
  negatives <- sum(y == 0L)
  if (positives == 0L || negatives == 0L) return(NA_real_)
  (sum(rank(scores, ties.method = "average")[y == 1L]) -
     positives * (positives + 1) / 2) / (positives * negatives)
}


compute_average_precision <- function(labels, probabilities) {
  keep <- is.finite(probabilities)
  y <- as.integer(labels[keep] == levels(labels)[2])
  if (sum(y) == 0L) return(NA_real_)
  order_index <- order(probabilities[keep], decreasing = TRUE)
  ordered_y <- y[order_index]
  precision <- cumsum(ordered_y) / seq_along(ordered_y)
  sum(precision[ordered_y == 1L]) / sum(ordered_y)
}


compute_threshold_metrics <- function(labels, probabilities, threshold = 0.5) {
  keep <- is.finite(probabilities)
  truth <- labels[keep]
  predicted <- factor(
    ifelse(probabilities[keep] >= threshold, levels(labels)[2], levels(labels)[1]),
    levels = levels(labels)
  )
  positive <- levels(labels)[2]
  negative <- levels(labels)[1]
  tp <- sum(truth == positive & predicted == positive)
  tn <- sum(truth == negative & predicted == negative)
  fp <- sum(truth == negative & predicted == positive)
  fn <- sum(truth == positive & predicted == negative)
  sensitivity <- if (tp + fn == 0) NA_real_ else tp / (tp + fn)
  specificity <- if (tn + fp == 0) NA_real_ else tn / (tn + fp)
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (denominator == 0) NA_real_ else (tp * tn - fp * fn) / denominator
  c(balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    mcc = mcc)
}


# ------------------------------------------------------------------------------
#  Fixed-panel repeated CV
# ------------------------------------------------------------------------------

folds <- make_stratified_folds(labels, CONFIG$outer_folds,
                               CONFIG$outer_repeats, CONFIG$random_seed)
panel_keys <- unique(panels[c("method", "panel_size")])
results <- list()
predictions <- list()

for (panel_index in seq_len(nrow(panel_keys))) {
  method_name <- panel_keys$method[panel_index]
  panel_size <- panel_keys$panel_size[panel_index]
  panel_rows <- panels[
    panels$method == method_name & panels$panel_size == panel_size,
    , drop = FALSE
  ]
  panel_rows <- panel_rows[order(panel_rows$rank), , drop = FALSE]
  available_genes <- intersect(panel_rows$gene, colnames(features))
  coverage <- length(available_genes) / panel_size
  if (coverage < CONFIG$minimum_gene_overlap) {
    warning(sprintf("%s k=%d has %.1f%% gene overlap and was skipped.",
                    method_name, panel_size, 100 * coverage))
    next
  }

  for (split in folds) {
    train_x <- features[split$train, available_genes, drop = FALSE]
    test_x <- features[split$test, available_genes, drop = FALSE]
    train_y <- droplevels(labels[split$train])
    test_y <- labels[split$test]

    means <- colMeans(train_x)
    standard_deviations <- apply(train_x, 2, stats::sd)
    keep_genes <- is.finite(standard_deviations) & standard_deviations > 0
    train_x <- sweep(train_x[, keep_genes, drop = FALSE], 2, means[keep_genes], "-")
    train_x <- sweep(train_x, 2, standard_deviations[keep_genes], "/")
    test_x <- sweep(test_x[, keep_genes, drop = FALSE], 2, means[keep_genes], "-")
    test_x <- sweep(test_x, 2, standard_deviations[keep_genes], "/")

    panel_seed_index <- match(panel_size, sort(unique(panel_keys$panel_size)))
    evaluation_seed <- 420000L + split$repeat_index * 1000L +
      split$fold * 100L + panel_seed_index
    probability_matrix <- predict_ensemble(train_x, train_y, test_x,
                                           evaluation_seed)
    for (evaluator in colnames(probability_matrix)) {
      probabilities <- probability_matrix[, evaluator]
      threshold_metrics <- compute_threshold_metrics(test_y, probabilities)
      numeric_truth <- as.integer(test_y == levels(labels)[2])
      results[[length(results) + 1L]] <- data.frame(
        method = method_name,
        panel_size = panel_size,
        n_locked_genes_present = length(available_genes),
        repeat_index = split$repeat_index,
        fold = split$fold,
        evaluator = evaluator,
        evaluation_seed = evaluation_seed,
        n_test = length(split$test),
        n_positive = sum(test_y == levels(labels)[2]),
        auc = compute_auc(test_y, probabilities),
        pr_auc = compute_average_precision(test_y, probabilities),
        balanced_accuracy = threshold_metrics[["balanced_accuracy"]],
        mcc = threshold_metrics[["mcc"]],
        brier = mean((probabilities - numeric_truth)^2, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      predictions[[length(predictions) + 1L]] <- data.frame(
        method = method_name, panel_size = panel_size,
        repeat_index = split$repeat_index, fold = split$fold,
        evaluator = evaluator,
        sample_id = metadata$sample_id[split$test],
        outcome = as.character(test_y), probability = probabilities,
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, results)
predictions <- do.call(rbind, predictions)
if (is.null(results) || nrow(results) == 0) {
  stop("No locked panel passed the gene-overlap requirement.")
}

metric_columns <- c("auc", "pr_auc", "balanced_accuracy", "mcc", "brier")
summary_rows <- lapply(
  split(results, interaction(results$method, results$panel_size,
                             results$evaluator, drop = TRUE)),
  function(group) {
    output <- group[1, c("method", "panel_size", "n_locked_genes_present",
                         "evaluator"), drop = FALSE]
    for (metric in metric_columns) {
      output[[paste0(metric, "_mean")]] <- mean(group[[metric]], na.rm = TRUE)
      output[[paste0(metric, "_sd")]] <- stats::sd(group[[metric]], na.rm = TRUE)
    }
    output$n_splits <- nrow(group)
    output
  }
)
summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL

run_stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
output_dir <- file.path(CONFIG$output_root, accession,
                        paste0(run_stamp, "_locked_panel"))
dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)

utils::write.csv(results, file.path(output_dir, "data", "fold_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(predictions, file.path(output_dir, "data", "predictions.csv"),
                 row.names = FALSE)
utils::write.csv(summary_table, file.path(output_dir, "data", "summary.csv"),
                 row.names = FALSE)

run_manifest <- list(
  holdout_accession = accession,
  source_accession = dataset$source_accession,
  locked_panel_file = panel_file,
  locked_panel_md5 = unname(tools::md5sum(panel_file)),
  source_rankings_md5 = unique(panels$source_rankings_md5),
  holdout_platform = unique(panels$holdout_platform),
  holdout_annotation_md5 = unique(panels$holdout_annotation_md5),
  locked_at_utc = unique(panels$locked_at_utc),
  aggregation_rule = unique(panels$aggregation_rule),
  random_seed = CONFIG$random_seed,
  outer_folds = CONFIG$outer_folds,
  outer_repeats = CONFIG$outer_repeats,
  n_samples = nrow(metadata),
  class_counts = as.list(table(labels)),
  minimum_gene_overlap = CONFIG$minimum_gene_overlap,
  evaluator_seed_rule = paste(
    "420000 + repeat_index*1000 + fold*100 + panel_size_index;",
    "glmnet stratified foldid; xgboost seed +1000; ranger seed +2000"
  )
)
saveRDS(run_manifest, file.path(output_dir, "data", "config.rds"))

cat(sprintf("Wrote locked-panel validation to %s\n", output_dir))
