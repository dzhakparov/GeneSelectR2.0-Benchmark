# ==============================================================================
#  Evaluator: the incumbent benchmark's three-model soft-voting ensemble,
#  copied verbatim from redesign/run_gs_slim_imvigor210.R (which copied it
#  from the incumbent benchmark). Do not tune anything in here: every arm in
#  the grouped benchmark must be scored by the same evaluator as the
#  incumbents it is compared against.
#
#  Provides: bench_auc, predict_with_ensemble, random_panel_aucs.
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

predict_with_ensemble <- function(train_features, train_labels, test_features,
                                  glmnet_alpha_grid = c(0.5, 1.0),
                                  random_seed = 42) {
  random_seed <- as.integer(random_seed)
  if (length(random_seed) != 1L || is.na(random_seed)) {
    stop("random_seed must be one non-missing integer", call. = FALSE)
  }
  numeric_labels <- as.integer(droplevels(train_labels) ==
                                 levels(droplevels(train_labels))[2])

  # Component 1: elastic net over the same alpha grid as the incumbent.
  # Use one explicit stratified fold assignment for every alpha value. This
  # makes comparisons paired at the model-fitting level and prevents method
  # order or prior RNG use from changing the reported AUC.
  n_inner <- min(5, max(3, floor(min(table(train_labels)) * 0.8)))
  fold_id <- withr::with_seed(random_seed, {
    result <- integer(length(numeric_labels))
    for (class_value in sort(unique(numeric_labels))) {
      indices <- which(numeric_labels == class_value)
      result[indices] <- sample(rep(seq_len(n_inner),
                                    length.out = length(indices)))
    }
    result
  })
  best_fit <- NULL; best_auc <- -Inf
  glmnet_errors <- character(0)
  for (alpha_value in glmnet_alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(train_features, numeric_labels,
                        family = "binomial", alpha = alpha_value,
                        foldid = fold_id, parallel = FALSE),
      error = function(e) {
        glmnet_errors <<- c(glmnet_errors, conditionMessage(e))
        NULL
      })
    if (!is.null(fit) && -min(fit$cvm) > best_auc) {
      best_auc <- -min(fit$cvm); best_fit <- fit
    }
  }
  if (is.null(best_fit)) {
    stop(sprintf("ensemble glmnet component failed: %s",
                 paste(unique(glmnet_errors), collapse = "; ")),
         call. = FALSE)
  }
  p_glmnet <- as.numeric(predict(best_fit, test_features,
                                 s = "lambda.min", type = "response"))

  # Component 2: XGBoost, incumbent defaults.
  p_xgb <- tryCatch({
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_features),
                                   label = numeric_labels)
    fit <- withr::with_seed(random_seed + 1000L, xgboost::xgb.train(
      params = list(objective = "binary:logistic", eval_metric = "auc",
                    eta = 0.1, max_depth = 3, subsample = 0.8,
                    colsample_bytree = 0.8, seed = random_seed + 1000L,
                    nthread = 1),
      data = dtrain, nrounds = 100, verbose = 0))
    predict(fit, as.matrix(test_features))
  }, error = function(e) {
    stop(sprintf("ensemble XGBoost component failed: %s",
                 conditionMessage(e)), call. = FALSE)
  })

  # Component 3: random forest, incumbent defaults.
  p_rf <- tryCatch({
    fit <- ranger::ranger(x = train_features, y = droplevels(train_labels),
                          num.trees = 500, probability = TRUE,
                          seed = random_seed + 2000L, num.threads = 1)
    predict(fit, data = test_features)$predictions[,
      levels(droplevels(train_labels))[2]]
  }, error = function(e) {
    stop(sprintf("ensemble random-forest component failed: %s",
                 conditionMessage(e)), call. = FALSE)
  })

  components <- cbind(p_glmnet, p_xgb, p_rf)
  if (any(!is.finite(components))) {
    stop("ensemble component returned non-finite predictions", call. = FALSE)
  }
  scores <- rowMeans(components)
  attr(scores, "components_used") <- c("glmnet", "xgboost", "ranger")
  scores
}

#  Random-panel baseline: mean AUC of `n_draws` random k-gene panels drawn
#  from the SAME pool the arm under test selected from. AUC minus this
#  baseline is the project's primary metric (random signatures inherit
#  dominant meta-genes, so raw AUC overstates every method).
random_panel_aucs <- function(train_features, test_features, train_labels,
                              test_labels, pool_genes, k, n_draws = 10,
                              seed = 99, model_seed = 42) {
  k_eff <- min(k, length(pool_genes))
  set.seed(seed)
  vapply(seq_len(n_draws), function(draw) {
    panel <- sample(pool_genes, k_eff)
    scores <- predict_with_ensemble(
      train_features[, panel, drop = FALSE], train_labels,
      test_features[, panel, drop = FALSE], random_seed = model_seed)
    bench_auc(test_labels, scores)
  }, numeric(1))
}
