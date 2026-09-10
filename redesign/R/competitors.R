# ==============================================================================
#  Competitor rankers for the full benchmark.
#
#  Lifted verbatim (logic) from benchmarks/validation_benchmark.R lines
#  1724-1957 so the comparison against the incumbent roster is like-for-like.
#  The script-level globals `top_variable_genes` and `random_seed` are explicit
#  arguments. A failed method stops the affected split. Assigning a random
#  ranking to a named comparator would change its measured performance while
#  retaining the original method label.
#
#  Each ranker takes train features (samples x genes) and labels (2-level
#  factor, positive second) and returns list(ranked, selected, gs_object).
# ==============================================================================

competitor_note <- function(method, step, err) {
  cat(sprintf("    [%s] %s failed: %s\n", method, step,
              conditionMessage(err)))
  NULL
}

competitor_failure <- function(method, step, err = NULL) {
  detail <- if (is.null(err)) "no valid result was returned" else
    conditionMessage(err)
  stop(sprintf("%s failed during %s: %s", method, step, detail),
       call. = FALSE)
}

make_competitor_fold_id <- function(train_labels, n_folds) {
  labels <- droplevels(train_labels)
  if (length(levels(labels)) != 2L || min(table(labels)) < n_folds) {
    stop("Competitor CV requires at least one observation from each class ",
         "in every validation fold", call. = FALSE)
  }
  fold_id <- integer(length(labels))
  for (class_value in levels(labels)) {
    class_indices <- which(labels == class_value)
    fold_id[class_indices] <- sample(rep(
      seq_len(n_folds), length.out = length(class_indices)
    ))
  }
  fold_id
}

# Differential expression, scored as -log10(p) * |t|.
rank_by_differential_expression <- function(train_features, train_labels) {
  train_labels <- droplevels(train_labels)
  n_genes      <- ncol(train_features)
  scores       <- numeric(n_genes)
  raw_p_values <- numeric(n_genes)
  names(scores)       <- colnames(train_features)
  names(raw_p_values) <- colnames(train_features)

  for (gene_idx in seq_len(n_genes)) {
    t_test <- tryCatch(t.test(train_features[, gene_idx] ~ train_labels),
                       error = function(e) NULL)
    if (is.null(t_test)) {
      scores[gene_idx]       <- 0
      raw_p_values[gene_idx] <- 1
    } else {
      t_statistic            <- unname(t_test$statistic)
      scores[gene_idx]       <- -log10(max(t_test$p.value, 1e-300)) *
        abs(t_statistic)
      raw_p_values[gene_idx] <- unname(t_test$p.value)
    }
  }

  adjusted_p     <- p.adjust(raw_p_values, method = "BH")
  selected_genes <- names(raw_p_values)[adjusted_p < 0.05]

  list(ranked   = names(sort(scores, decreasing = TRUE)),
       selected = selected_genes,
       gs_object = NULL)
}

# LASSO with internal CV for lambda.
rank_by_lasso <- function(train_features, train_labels) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  n_inner_folds  <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))
  fold_id <- make_competitor_fold_id(train_labels, n_inner_folds)

  fit <- tryCatch(
    glmnet::cv.glmnet(train_features, numeric_labels,
                      family = "binomial", alpha = 1,
                      nfolds = n_inner_folds, foldid = fold_id),
    error = function(e) competitor_failure("LASSO", "cv.glmnet", e)
  )

  coefficients <- as.numeric(coef(fit, s = "lambda.min"))[-1]
  names(coefficients) <- colnames(train_features)

  list(ranked   = names(sort(abs(coefficients), decreasing = TRUE)),
       selected = names(coefficients)[coefficients != 0],
       gs_object = NULL)
}

# Elastic net with alpha tuning.
rank_by_elastic_net <- function(train_features, train_labels) {
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  n_inner_folds  <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))
  fold_id <- make_competitor_fold_id(train_labels, n_inner_folds)

  best_fit      <- NULL
  best_cv_error <- Inf

  for (alpha_value in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    fit <- tryCatch(
      glmnet::cv.glmnet(train_features, numeric_labels,
                        family = "binomial", alpha = alpha_value,
                        nfolds = n_inner_folds, foldid = fold_id),
      error = function(e) {
        competitor_note("ElasticNet", sprintf("cv.glmnet alpha=%.2f",
                                                alpha_value), e)
      }
    )
    if (!is.null(fit) && min(fit$cvm) < best_cv_error) {
      best_cv_error <- min(fit$cvm)
      best_fit      <- fit
    }
  }

  if (is.null(best_fit)) competitor_failure("ElasticNet", "alpha search")

  coefficients <- as.numeric(coef(best_fit, s = "lambda.min"))[-1]
  names(coefficients) <- colnames(train_features)

  list(ranked   = names(sort(abs(coefficients), decreasing = TRUE)),
       selected = names(coefficients)[coefficients != 0],
       gs_object = NULL)
}

# Minimum-redundancy maximum-relevance.
rank_by_mrmr <- function(train_features, train_labels,
                         top_variable_genes = 2000) {
  numeric_labels <- as.numeric(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])
  max_features <- min(ncol(train_features), top_variable_genes)

  feature_subset <- train_features[, 1:max_features, drop = FALSE]
  original_names <- colnames(feature_subset)

  # mRMRe mangles non-syntactic column names via make.names(); safe
  # placeholders are mapped back afterwards.
  safe_names <- paste0("V", seq_len(max_features))
  input_df <- data.frame(
    outcome = as.numeric(numeric_labels),
    matrix(as.numeric(feature_subset), nrow = nrow(feature_subset)),
    stringsAsFactors = FALSE
  )
  colnames(input_df) <- c("outcome", safe_names)

  mrmr_data <- tryCatch(mRMRe::mRMR.data(data = input_df),
                        error = function(e) {
                          competitor_failure("mRMR", "mRMR.data", e)
                        })

  mrmr_result <- tryCatch(
    mRMRe::mRMR.classic(data = mrmr_data,
                        target_indices = 1,
                        feature_count = min(200, max_features)),
    error = function(e) competitor_failure("mRMR", "mRMR.classic", e)
  )

  selected_indices <- tryCatch(
    as.integer(mRMRe::solutions(mrmr_result)[[1]]),
    error = function(e) competitor_failure("mRMR", "solutions", e)
  )
  if (is.null(selected_indices) || length(selected_indices) == 0) {
    competitor_failure("mRMR", "solutions")
  }

  selected_indices <- selected_indices[!is.na(selected_indices) &
                                         selected_indices >= 2 &
                                         selected_indices <= (max_features + 1)]
  if (length(selected_indices) == 0) {
    competitor_failure("mRMR", "solution index validation")
  }

  selected_genes  <- original_names[selected_indices - 1L]
  remaining_genes <- setdiff(colnames(train_features), selected_genes)

  list(ranked   = c(selected_genes, remaining_genes),
       selected = selected_genes,
       gs_object = NULL)
}

# Boruta wrapper around random forest.
rank_by_boruta <- function(train_features, train_labels,
                           top_variable_genes = 2000) {
  train_labels <- droplevels(train_labels)
  max_features <- min(ncol(train_features), top_variable_genes)

  if (ncol(train_features) > max_features) {
    variances      <- apply(train_features, 2, var)
    top_indices    <- order(variances, decreasing = TRUE)[1:max_features]
    train_features <- train_features[, top_indices, drop = FALSE]
  }

  boruta_result <- tryCatch(
    Boruta::Boruta(x = train_features, y = train_labels,
                   doTrace = 0, maxRuns = 50),
    error = function(e) competitor_failure("Boruta", "Boruta", e)
  )

  initial_stats <- Boruta::attStats(boruta_result)
  if (any(initial_stats$decision == "Tentative")) {
    boruta_result <- tryCatch(Boruta::TentativeRoughFix(boruta_result),
                              error = function(e) {
                                warning(sprintf(
                                  "Boruta tentative resolution failed: %s",
                                  conditionMessage(e)), call. = FALSE)
                                boruta_result
                              })
  }
  importance_stats <- Boruta::attStats(boruta_result)

  #  attStats() returns a data.frame whose ROWNAMES are check.names-mangled
  #  for non-syntactic gene symbols (about 3% of a real pool). Match by
  #  POSITION instead: Boruta preserves input column order in its stats.
  stopifnot(nrow(importance_stats) == ncol(train_features))
  feature_order <- order(-importance_stats$meanImp)
  ranked_genes  <- colnames(train_features)[feature_order]

  #  On small subsamples with many features Boruta often confirms NOTHING.
  #  That is a genuine property of the method here, not an error.
  selected_genes <- colnames(train_features)[
    importance_stats$decision == "Confirmed"
  ]

  list(ranked   = ranked_genes,
       selected = selected_genes,
       gs_object = NULL)
}

# Random forest variable importance.
rank_by_random_forest <- function(train_features, train_labels,
                                  top_variable_genes = 2000, random_seed = 42) {
  train_labels <- droplevels(train_labels)
  max_features <- min(ncol(train_features), top_variable_genes)

  if (ncol(train_features) > max_features) {
    variances      <- apply(train_features, 2, var)
    top_indices    <- order(variances, decreasing = TRUE)[1:max_features]
    train_features <- train_features[, top_indices, drop = FALSE]
  }

  rf_fit <- tryCatch(
    ranger::ranger(x = train_features, y = train_labels,
                   num.trees = 1000, importance = "impurity",
                   seed = random_seed),
    error = function(e) {
      competitor_failure("RF_importance", "ranger", e)
    }
  )

  sorted_importance <- sort(rf_fit$variable.importance, decreasing = TRUE)
  remaining_genes   <- setdiff(colnames(train_features),
                               names(sorted_importance))
  list(ranked   = c(names(sorted_importance), remaining_genes),
       selected = NULL,
       gs_object = NULL)
}

# Random baseline -- the most important comparator.
rank_at_random <- function(train_features, train_labels) {
  list(ranked = sample(colnames(train_features)),
       selected = NULL,
       gs_object = NULL)
}
