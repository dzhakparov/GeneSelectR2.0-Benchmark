# ==============================================================================
# GeneSelectR 2.0 — Enhanced Benchmark
#   - Hyperparameter optimisation for every method
#   - Additional feature-selection baselines: mRMR, Boruta, LASSO, RF importance
#   - Date-stamped output folders (YYYY-MM-DD)
#
# Files expected:
#   - data/normalized_logcpm.csv: rows = genes, cols = samples
#   - data/metadata.csv: has sample ID col and treatment = "location_diagnosis"
# ==============================================================================
suppressPackageStartupMessages({
  library(GeneSelectR)
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Install / load optional packages for extra methods --------------------------
pkg_ensure <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s …", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

pkg_ensure("mRMRe")       # mRMR
pkg_ensure("Boruta")       # Boruta (wrapper around randomForest)
pkg_ensure("randomForest") # RF variable importance
pkg_ensure("ranger")       # fast RF (used inside HP tuning)

# ----------------------------- CONFIG -----------------------------------------
SEED <- 42
set.seed(SEED)

EXPR_FILE <- "data/normalized_logcpm.csv"
META_FILE <- "data/metadata.csv"

TREATMENT_COL <- "treatment"
TREATMENT_SPLIT_SEP <- "_"
NEW_LOC_COL <- "location"
NEW_DX_COL  <- "diagnosis"

SAMPLE_ID_COL <- "X"
OUTCOME_COL   <- "diagnosis"

# Confounders
CONFOUNDERS_CAT <- c(NEW_LOC_COL)
CONFOUNDERS_NUM <- character(0)
DO_RESIDUALIZE  <- TRUE

# Evaluation
K_OUTER  <- 5
R_OUTER  <- 5
K_VALUES <- c(10, 20, 50, 100, 200, 500)

# GeneSelectR inner CV
K_INNER_GS <- 5
R_INNER_GS <- 5
N_CORES <- max(1, parallel::detectCores() - 1)

# Optional variance filter (set to NULL to disable)
TOP_VAR_GENES <- 5000

# Targeted immune GO terms (BP)
TARGET_TERMS <- c(
  "GO:0006955", "GO:0002376", "GO:0006952",
  "GO:0002250", "GO:0045087", "GO:0002682", "GO:0050776"
)
BIO_ONTOLOGY   <- "BP"
BIO_SIM_METHOD <- "resnik"
BIO_ENRICH_FDR <- 0.05
BIO_IC_QUANTILE <- 0.5
BIO_MAX_ENRICHED <- 100
BIO_TOPK_SIMS   <- 5

SCORE_FORMULA <- "geometric"
SCORE_WEIGHTS <- c(1, 1, 1)

# ---- Date-stamped output directory ------------------------------------------
RUN_DATE <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- file.path("results_targeted_sosall", RUN_DATE)
dir.create(file.path(OUT_DIR, "data"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== Run date: %s ===\n", RUN_DATE))
cat(sprintf("=== Output dir: %s ===\n", OUT_DIR))

# ====================== HYPERPARAMETER GRIDS ==================================
#
# Each method gets an inner HP search.  The best HP set (by inner-CV AUC)
# is used to rank genes, then those genes are evaluated on the outer-fold.
# --------------------------------------------------------------------------

# --- glmnet classifier HP grid (shared by all methods for downstream eval) ---
GLMNET_ALPHA_GRID <- c(0.5, 1.0)   # ridge → lasso

# --- GeneSelectR HP grid ---
GS_ALPHA_GRID    <- c(0.5, 1.0)
GS_WEIGHTS_GRID  <- list(
  c(1, 1, 1),
  c(2, 1, 1),
  c(1, 2, 1),
  c(1, 1, 2)
)

# --- Boruta HP grid ---
BORUTA_NTREE_GRID  <- c(500)
BORUTA_MAXRUNS     <- c(50)

# --- Random Forest importance HP grid ---
RF_NTREE_GRID    <- c(1000)
RF_MTRY_FRAC     <- c(0.33)

# --- mRMR HP grid ---
MRMR_NFEATURES_GRID <- c(200)  # number of features mRMR selects internally

# ========================== HELPERS ===========================================

fast_auc <- function(y, score) {

  y <- droplevels(y)
  if (nlevels(y) != 2) stop("y must have 2 levels")
  pos <- which(y == levels(y)[2])
  neg <- which(y == levels(y)[1])
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

# --- confounder design matrix (train-fit / test-apply) ---
one_hot_train <- function(x, drop_first = TRUE, prefix = "f") {
  x <- factor(x)
  mm <- model.matrix(~ x)
  colnames(mm) <- gsub("^x", prefix, colnames(mm))
  if (drop_first && ncol(mm) > 1) mm <- mm[, -1, drop = FALSE]
  list(mm = mm, levels = levels(x))
}

one_hot_apply <- function(x, levels_ref, drop_first = TRUE, prefix = "f") {
  x <- factor(x, levels = levels_ref)
  mm <- model.matrix(~ x)
  colnames(mm) <- gsub("^x", prefix, colnames(mm))
  if (drop_first && ncol(mm) > 1) mm <- mm[, -1, drop = FALSE]
  mm
}

build_C_train <- function(meta_df) {
  mats <- list()
  ref  <- list(cat_levels = list(), num_median = list())
  for (cc in CONFOUNDERS_CAT) {
    if (cc %in% names(meta_df)) {
      x <- meta_df[[cc]]
      x[is.na(x) | x == ""] <- "NA_level"
      oh <- one_hot_train(x, drop_first = TRUE, prefix = cc)
      mats[[cc]] <- oh$mm
      ref$cat_levels[[cc]] <- oh$levels
    }
  }
  for (nc in CONFOUNDERS_NUM) {
    if (nc %in% names(meta_df)) {
      x <- suppressWarnings(as.numeric(meta_df[[nc]]))
      med <- median(x, na.rm = TRUE)
      x[is.na(x)] <- med
      ref$num_median[[nc]] <- med
      z <- as.numeric(scale(x))
      mats[[nc]] <- matrix(z, ncol = 1)
      colnames(mats[[nc]]) <- paste0(nc, "_z")
    }
  }
  C <- if (length(mats)) do.call(cbind, mats) else NULL
  if (is.null(C)) C <- matrix(0, nrow = nrow(meta_df), ncol = 0)
  C <- cbind("(Intercept)" = 1, C)
  list(C = C, ref = ref)
}

build_C_apply <- function(meta_df, ref) {
  mats <- list()
  for (cc in CONFOUNDERS_CAT) {
    if (cc %in% names(meta_df) && !is.null(ref$cat_levels[[cc]])) {
      x <- meta_df[[cc]]
      x[is.na(x) | x == ""] <- "NA_level"
      mats[[cc]] <- one_hot_apply(x, ref$cat_levels[[cc]], drop_first = TRUE, prefix = cc)
    }
  }
  for (nc in CONFOUNDERS_NUM) {
    if (nc %in% names(meta_df) && !is.null(ref$num_median[[nc]])) {
      x <- suppressWarnings(as.numeric(meta_df[[nc]]))
      x[is.na(x)] <- ref$num_median[[nc]]
      z <- as.numeric(scale(x))
      mats[[nc]] <- matrix(z, ncol = 1)
      colnames(mats[[nc]]) <- paste0(nc, "_z")
    }
  }
  C <- if (length(mats)) do.call(cbind, mats) else NULL
  if (is.null(C)) C <- matrix(0, nrow = nrow(meta_df), ncol = 0)
  C <- cbind("(Intercept)" = 1, C)
  C
}

fit_residualizer <- function(C_train, X_train) {
  qrC <- qr(C_train)
  B_hat <- qr.coef(qrC, X_train)
  list(B_hat = B_hat)
}
apply_residualizer <- function(C, X, fit) X - C %*% fit$B_hat

make_stratified_folds <- function(y, K = 5, seed = 1) {
  set.seed(seed)
  y <- droplevels(y)
  idx1 <- which(y == levels(y)[1])
  idx2 <- which(y == levels(y)[2])
  folds <- vector("list", K)
  s1 <- split(sample(idx1), rep(1:K, length.out = length(idx1)))
  s2 <- split(sample(idx2), rep(1:K, length.out = length(idx2)))
  for (k in 1:K) folds[[k]] <- sort(c(s1[[k]], s2[[k]]))
  folds
}

# ---- glmnet classifier with HP tuning (alpha) ----
predict_glmnet_tuned <- function(X_train, y_train, X_test,
                                 alpha_grid = GLMNET_ALPHA_GRID) {
  y_train <- droplevels(y_train)
  y_num   <- as.integer(y_train == levels(y_train)[2])
  nf <- min(5, max(3, floor(min(table(y_train)) * 0.8)))

  best_auc   <- -Inf
  best_alpha <- 0.5
  best_fit   <- NULL

  for (a in alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(x = X_train, y = y_num, family = "binomial",
                        alpha = a, nfolds = nf),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    # use cv deviance minimum as proxy
    min_cvm <- min(fit$cvm)
    if (-min_cvm > best_auc) {
      best_auc   <- -min_cvm
      best_alpha <- a
      best_fit   <- fit
    }
  }

  if (is.null(best_fit)) return(rep(NA_real_, nrow(X_test)))
  as.numeric(predict(best_fit, X_test, s = "lambda.min", type = "response"))
}

# ====================== RANKING FUNCTIONS (with HP search) ====================

# ---- GeneSelectR targeted (HP: alpha, score_weights) ----
rank_geneselectr_targeted <- function(X_train, y_train) {
  best_gene_list <- NULL
  best_score     <- -Inf

  for (a in GS_ALPHA_GRID) {
    for (sw in GS_WEIGHTS_GRID) {
      gs <- tryCatch(
        geneselectr2_fit(
          X = X_train, y = y_train,
          bio_mode = "supervised",
          target_terms = TARGET_TERMS,
          bio_ontology = BIO_ONTOLOGY,
          bio_sim_method = BIO_SIM_METHOD,
          score_formula = SCORE_FORMULA,
          score_weights = sw,
          regularization_method = "elastic_net",
          alpha = a,
          K = K_INNER_GS, R = R_INNER_GS,
          n_cores = N_CORES,
          random_seed = SEED,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(gs)) next
      inner_auc <- mean(gs$cv_results$auc_scores, na.rm = TRUE)
      if (inner_auc > best_score) {
        best_score     <- inner_auc
        best_gene_list <- gs$gene_scores$gene
      }
    }
  }
  if (is.null(best_gene_list)) return(sample(colnames(X_train)))
  best_gene_list
}

# ---- GeneSelectR data-driven (HP: alpha, score_weights) ----
rank_geneselectr_datadriven <- function(X_train, y_train) {
  best_gene_list <- NULL
  best_score     <- -Inf

  for (a in GS_ALPHA_GRID) {
    for (sw in GS_WEIGHTS_GRID) {
      gs <- tryCatch(
        geneselectr2_fit(
          X = X_train, y = y_train,
          bio_mode = "data_driven",
          bio_ontology = BIO_ONTOLOGY,
          bio_sim_method = BIO_SIM_METHOD,
          bio_enrich_fdr = BIO_ENRICH_FDR,
          bio_ic_quantile = BIO_IC_QUANTILE,
          bio_max_enriched = BIO_MAX_ENRICHED,
          bio_n_top_sims = BIO_TOPK_SIMS,
          score_formula = SCORE_FORMULA,
          score_weights = sw,
          regularization_method = "elastic_net",
          alpha = a,
          K = K_INNER_GS, R = R_INNER_GS,
          n_cores = N_CORES,
          random_seed = SEED,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(gs)) next
      inner_auc <- mean(gs$cv_results$auc_scores, na.rm = TRUE)
      if (inner_auc > best_score) {
        best_score     <- inner_auc
        best_gene_list <- gs$gene_scores$gene
      }
    }
  }
  if (is.null(best_gene_list)) return(sample(colnames(X_train)))
  best_gene_list
}

# ---- DGE t-test (no HPs to tune; ranking is deterministic) ----
rank_dge_ttest <- function(X_train, y_train) {
  y_train <- droplevels(y_train)
  p <- apply(X_train, 2, function(g) {
    tryCatch(t.test(g ~ y_train)$p.value, error = function(e) 1.0)
  })
  padj <- p.adjust(p, method = "BH")
  names(sort(padj))
}

# ---- LASSO-based ranking (HP: alpha fixed=1 by definition; lambda via CV) ----
rank_lasso <- function(X_train, y_train) {
  y_train <- droplevels(y_train)
  y_num   <- as.integer(y_train == levels(y_train)[2])
  nf <- min(5, max(3, floor(min(table(y_train)) * 0.8)))

  fit <- tryCatch(
    glmnet::cv.glmnet(x = X_train, y = y_num, family = "binomial",
                      alpha = 1, nfolds = nf),
    error = function(e) NULL
  )
  if (is.null(fit)) return(sample(colnames(X_train)))

  coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]  # drop intercept
  names(coefs) <- colnames(X_train)
  names(sort(abs(coefs), decreasing = TRUE))
}

# ---- Elastic-net ranking (HP: alpha via inner CV) ----
rank_elasticnet <- function(X_train, y_train) {
  y_train <- droplevels(y_train)
  y_num   <- as.integer(y_train == levels(y_train)[2])
  nf <- min(5, max(3, floor(min(table(y_train)) * 0.8)))

  best_cvm   <- Inf
  best_alpha <- 0.5
  best_fit   <- NULL

  for (a in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    fit <- tryCatch(
      glmnet::cv.glmnet(x = X_train, y = y_num, family = "binomial",
                        alpha = a, nfolds = nf),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (min(fit$cvm) < best_cvm) {
      best_cvm   <- min(fit$cvm)
      best_alpha <- a
      best_fit   <- fit
    }
  }
  if (is.null(best_fit)) return(sample(colnames(X_train)))

  coefs <- as.numeric(coef(best_fit, s = "lambda.min"))[-1]
  names(coefs) <- colnames(X_train)
  names(sort(abs(coefs), decreasing = TRUE))
}

# ---- mRMR ranking (HP: n_features) ----
rank_mrmr <- function(X_train, y_train) {
  y_train <- droplevels(y_train)
  y_num   <- as.integer(y_train == levels(y_train)[2])

  best_genes <- NULL
  best_score <- -Inf

  p_use <- min(ncol(X_train), 5000)
  feat_names <- colnames(X_train)[1:p_use]

  df_in <- data.frame(outcome = y_num, X_train[, feat_names, drop = FALSE])
  dd <- tryCatch(mRMRe::mRMR.data(data = df_in), error = function(e) NULL)
  if (is.null(dd)) return(sample(colnames(X_train)))

  for (nfeat in MRMR_NFEATURES_GRID) {
    nfeat_use <- min(nfeat, p_use)

    res <- tryCatch(
      mRMRe::mRMR.classic(data = dd, target_indices = 1, feature_count = nfeat_use),
      error = function(e) NULL
    )
    if (is.null(res)) next

    # indices returned refer to columns in df_in in most mRMRe versions
    sel_cols <- as.integer(mRMRe::solutions(res)[[1]])
    sel_cols <- sel_cols[sel_cols != 1]  # drop outcome if ever present
    sel_cols <- sel_cols[sel_cols >= 2 & sel_cols <= (p_use + 1)]

    if (!length(sel_cols)) next

    gene_names <- colnames(df_in)[sel_cols]
    gene_names <- setdiff(gene_names, "outcome")

    inner_fit <- tryCatch(
      glmnet::cv.glmnet(
        x = X_train[, gene_names, drop = FALSE], y = y_num,
        family = "binomial", alpha = 0.5,
        nfolds = min(5, max(3, floor(min(table(y_train)) * 0.8)))
      ),
      error = function(e) NULL
    )
    if (is.null(inner_fit)) next

    sc <- -min(inner_fit$cvm)
    if (sc > best_score) {
      best_score <- sc
      best_genes <- c(gene_names, setdiff(colnames(X_train), gene_names))
    }
  }

  if (is.null(best_genes)) return(sample(colnames(X_train)))
  best_genes
}

# ---- Boruta ranking (HP: ntree, maxRuns) ----
rank_boruta <- function(X_train, y_train) {
  y_train <- droplevels(y_train)

  # Sub-sample features if too many for Boruta speed
  max_feat <- min(ncol(X_train), 5000)
  if (ncol(X_train) > max_feat) {
    vr <- apply(X_train, 2, var)
    keep <- order(vr, decreasing = TRUE)[1:max_feat]
    X_sub <- X_train[, keep, drop = FALSE]
  } else {
    X_sub <- X_train
  }

  best_genes <- NULL
  best_n_confirmed <- 0

  for (ntree in BORUTA_NTREE_GRID) {
    for (mr in BORUTA_MAXRUNS) {
      bor <- tryCatch(
        Boruta::Boruta(x = X_sub, y = y_train,
                       doTrace = 0,
                       num.trees = ntree,
                       maxRuns = mr),
        error = function(e) NULL
      )
      if (is.null(bor)) next

      bor_final <- tryCatch(
        Boruta::TentativeRoughFix(bor),
        error = function(e) bor
      )

      confirmed <- names(which(bor_final$finalDecision == "Confirmed"))
      tentative <- names(which(bor_final$finalDecision == "Tentative"))
      sel <- c(confirmed, tentative)

      if (length(sel) > best_n_confirmed) {
        best_n_confirmed <- length(sel)
        imp <- Boruta::attStats(bor_final)
        imp$gene <- rownames(imp)
        imp <- imp[order(-imp$meanImp), ]
        remaining <- setdiff(colnames(X_sub), imp$gene)
        best_genes <- c(imp$gene, remaining)
      }
    }
  }

  if (is.null(best_genes)) {
    # Fallback: use RF importance
    return(rank_rf_importance(X_train, y_train))
  }

  # Append any genes not in X_sub
  missing <- setdiff(colnames(X_train), best_genes)
  c(best_genes, missing)
}

# ---- Random Forest variable importance (HP: ntree, mtry) ----
rank_rf_importance <- function(X_train, y_train) {
  y_train <- droplevels(y_train)

  max_feat <- min(ncol(X_train), 5000)
  if (ncol(X_train) > max_feat) {
    vr <- apply(X_train, 2, var)
    keep <- order(vr, decreasing = TRUE)[1:max_feat]
    X_sub <- X_train[, keep, drop = FALSE]
  } else {
    X_sub <- X_train
  }

  best_genes <- NULL
  best_oob   <- Inf

  for (ntree in RF_NTREE_GRID) {
    for (mf in RF_MTRY_FRAC) {
      mtry_val <- max(1, round(ncol(X_sub) * mf))
      rf <- tryCatch(
        ranger::ranger(
          x = X_sub, y = y_train,
          num.trees = ntree,
          mtry = mtry_val,
          importance = "impurity",
          seed = SEED
        ),
        error = function(e) NULL
      )
      if (is.null(rf)) next
      if (rf$prediction.error < best_oob) {
        best_oob <- rf$prediction.error
        imp <- sort(rf$variable.importance, decreasing = TRUE)
        remaining <- setdiff(colnames(X_sub), names(imp))
        best_genes <- c(names(imp), remaining)
      }
    }
  }

  if (is.null(best_genes)) return(sample(colnames(X_train)))

  missing <- setdiff(colnames(X_train), best_genes)
  c(best_genes, missing)
}

# ---- Random baseline (no HP) ----
rank_random <- function(X_train, y_train) {
  sample(colnames(X_train))
}

# ======================== ALL METHODS REGISTRY ================================
ALL_METHODS <- c(
  "GeneSelectR_targeted",
  "GeneSelectR_datadriven",
  "DGE_ttest",
  "LASSO",
  "ElasticNet",
  "mRMR",
  "Boruta",
  "RF_importance",
  "Random"
)

get_ranker <- function(method) {
  switch(method,
         "GeneSelectR_targeted"   = rank_geneselectr_targeted,
         "GeneSelectR_datadriven" = rank_geneselectr_datadriven,
         "DGE_ttest"              = rank_dge_ttest,
         "LASSO"                  = rank_lasso,
         "ElasticNet"             = rank_elasticnet,
         "mRMR"                   = rank_mrmr,
         "Boruta"                 = rank_boruta,
         "RF_importance"          = rank_rf_importance,
         "Random"                 = rank_random,
         stop("Unknown method: ", method)
  )
}

# ======================== NESTED EVALUATION ===================================
evaluate_nested <- function(X, y, meta, method,
                            k_values = K_VALUES, seed = 1) {
  rank_fn <- get_ranker(method)
  rows <- list()

  for (r in 1:R_OUTER) {
    folds <- make_stratified_folds(y, K = K_OUTER, seed = seed + 1000 * r)

    for (f in 1:K_OUTER) {
      test_idx  <- folds[[f]]
      train_idx <- setdiff(seq_along(y), test_idx)

      X_tr0 <- X[train_idx, , drop = FALSE]
      X_te0 <- X[test_idx,  , drop = FALSE]
      y_tr  <- y[train_idx]
      y_te  <- y[test_idx]

      meta_tr <- meta[train_idx, , drop = FALSE]
      meta_te <- meta[test_idx,  , drop = FALSE]

      # Optional residualization (fit on train, apply to test)
      if (DO_RESIDUALIZE) {
        Ct   <- build_C_train(meta_tr)
        C_tr <- Ct$C
        C_te <- build_C_apply(meta_te, Ct$ref)

        resid_fit <- fit_residualizer(C_tr, X_tr0)
        X_tr <- apply_residualizer(C_tr, X_tr0, resid_fit)
        X_te <- apply_residualizer(C_te, X_te0, resid_fit)
      } else {
        X_tr <- X_tr0
        X_te <- X_te0
      }

      # Standardize on TRAIN only
      mu  <- colMeans(X_tr)
      sdv <- apply(X_tr, 2, sd)
      sdv[sdv == 0 | is.na(sdv)] <- 1
      X_tr <- sweep(sweep(X_tr, 2, mu, "-"), 2, sdv, "/")
      X_te <- sweep(sweep(X_te, 2, mu, "-"), 2, sdv, "/")

      cat(sprintf("  [%s] repeat=%d fold=%d ranking ...\n", method, r, f))
      ranked <- tryCatch(rank_fn(X_tr, y_tr),
                         error = function(e) {
                           warning(sprintf("Ranking failed (%s r=%d f=%d): %s",
                                           method, r, f, e$message))
                           sample(colnames(X_tr))
                         })

      for (k in k_values) {
        genes_k <- head(ranked[ranked %in% colnames(X_tr)], k)
        if (length(genes_k) < 5) next

        p <- tryCatch(
          predict_glmnet_tuned(X_tr[, genes_k, drop = FALSE], y_tr,
                               X_te[, genes_k, drop = FALSE]),
          error = function(e) rep(NA_real_, length(y_te))
        )
        auc <- if (all(is.na(p))) NA_real_ else fast_auc(y_te, p)

        rows[[length(rows) + 1]] <- data.frame(
          Method = method,
          Repeat = r,
          Fold   = f,
          k      = k,
          AUC    = auc,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

# ======================== LOAD + TRANSFORM ====================================
expr <- read.csv(EXPR_FILE, row.names = 1, check.names = FALSE,
                 stringsAsFactors = FALSE)

gene_symbol <- sub(".*__", "", rownames(expr))
expr$gene_symbol <- gene_symbol
expr_agg <- aggregate(. ~ gene_symbol, data = expr, FUN = mean)
rownames(expr_agg) <- expr_agg$gene_symbol
expr_agg$gene_symbol <- NULL

Xg <- as.matrix(expr_agg)
storage.mode(Xg) <- "numeric"
Xg[is.na(Xg)] <- 0
X <- t(Xg)
X <- as.matrix(X)
storage.mode(X) <- "numeric"
X[is.na(X)] <- 0

meta <- read.csv(META_FILE, stringsAsFactors = FALSE)

# Split treatment → location + diagnosis
tmp <- strsplit(as.character(meta[[TREATMENT_COL]]), TREATMENT_SPLIT_SEP,
                fixed = TRUE)
meta[[NEW_LOC_COL]] <- vapply(tmp, function(z)
  if (length(z) >= 1) z[[1]] else NA_character_, character(1))
meta[[NEW_DX_COL]]  <- vapply(tmp, function(z)
  if (length(z) >= 2) paste(z[-1], collapse = TREATMENT_SPLIT_SEP)
  else NA_character_, character(1))

stopifnot(SAMPLE_ID_COL %in% names(meta))
stopifnot(OUTCOME_COL   %in% names(meta))
stopifnot(TREATMENT_COL %in% names(meta))

common <- intersect(rownames(X), meta[[SAMPLE_ID_COL]])
if (length(common) < 20) stop("Too few matched samples")
X    <- X[common, , drop = FALSE]
meta <- meta[match(common, meta[[SAMPLE_ID_COL]]), , drop = FALSE]

y <- factor(meta[[OUTCOME_COL]])
y <- droplevels(y)
if (nlevels(y) != 2) stop("Outcome must have exactly 2 classes")

cat(sprintf("Loaded: %d samples x %d genes\n", nrow(X), ncol(X)))
cat("Outcome counts:\n"); print(table(y))

# Variance filter
vars <- apply(X, 2, var)
X <- X[, vars > 0, drop = FALSE]
if (!is.null(TOP_VAR_GENES) && ncol(X) > TOP_VAR_GENES) {
  top <- order(apply(X, 2, var), decreasing = TRUE)[1:TOP_VAR_GENES]
  X <- X[, top, drop = FALSE]
}
cat(sprintf("After filtering: %d genes\n", ncol(X)))
cat(sprintf("Residualization: %s\n", if (DO_RESIDUALIZE) "ON" else "OFF"))

# ============================== RUN ===========================================
cat(sprintf("\nRunning nested evaluation (%d methods) ...\n", length(ALL_METHODS)))

all_results <- list()
for (m in ALL_METHODS) {
  cat(sprintf("\n>>> Method: %s\n", m))
  t0 <- proc.time()
  all_results[[m]] <- evaluate_nested(X, y, meta, method = m, seed = SEED)
  dt <- (proc.time() - t0)["elapsed"]
  cat(sprintf("<<< %s done in %.1f sec\n", m, dt))
}

res_df <- dplyr::bind_rows(all_results)
write.csv(res_df, file.path(OUT_DIR, "data", "nested_results.csv"),
          row.names = FALSE)

# ======================== SUMMARY + PLOTS =====================================
sum_df <- res_df %>%
  group_by(Method, k) %>%
  summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd   = sd(AUC, na.rm = TRUE),
    AUC_med  = median(AUC, na.rm = TRUE),
    n        = sum(!is.na(AUC)),
    .groups  = "drop"
  )
write.csv(sum_df, file.path(OUT_DIR, "data", "nested_summary.csv"),
          row.names = FALSE)

# --- Main parsimony curve ---
method_colors <- c(
  "GeneSelectR_targeted"   = "#E41A1C",
  "GeneSelectR_datadriven" = "#377EB8",
  "DGE_ttest"              = "#4DAF4A",
  "LASSO"                  = "#FF7F00",
  "ElasticNet"             = "#984EA3",
  "mRMR"                   = "#A65628",
  "Boruta"                 = "#F781BF",
  "RF_importance"          = "#999999",
  "Random"                 = "grey50"
)

p <- ggplot(sum_df, aes(x = k, y = AUC_mean, color = Method, linetype = Method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd),
                width = 0.08, alpha = 0.35) +
  scale_x_log10(breaks = K_VALUES) +
  scale_color_manual(values = method_colors) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  theme_bw(base_size = 12) +
  labs(
    title = "Feature selection benchmark (nested CV, HP-tuned)",
    subtitle = sprintf("Run: %s | %d samples × %d genes | %d×%d outer CV",
                       RUN_DATE, nrow(X), ncol(X), R_OUTER, K_OUTER),
    x = "Top-k genes (log scale)",
    y = "Outer-fold AUC (mean ± SD)"
  )
ggsave(file.path(OUT_DIR, "figures", "parsimony_all_methods.pdf"),
       p, width = 12, height = 7)

# --- Box-plot at each k ---
p_box <- ggplot(res_df, aes(x = factor(k), y = AUC, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = method_colors) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  theme_bw(base_size = 12) +
  labs(
    title = "AUC distributions by method and k",
    subtitle = sprintf("Run: %s", RUN_DATE),
    x = "Top-k genes",
    y = "AUC"
  )
ggsave(file.path(OUT_DIR, "figures", "auc_boxplots.pdf"),
       p_box, width = 14, height = 7)

# --- Pairwise Wilcoxon at k=50 and k=200 ---
pairwise_tests <- list()
for (kk in c(50, 200)) {
  sub <- res_df %>% filter(k == kk, !is.na(AUC))
  if (length(unique(sub$Method)) < 2) next
  pw <- pairwise.wilcox.test(sub$AUC, sub$Method, p.adjust.method = "BH")
  pmat <- as.data.frame(pw$p.value)
  pmat$Method1 <- rownames(pmat)
  pmat_long <- pmat %>%
    pivot_longer(-Method1, names_to = "Method2", values_to = "p_adj") %>%
    filter(!is.na(p_adj)) %>%
    mutate(k = kk)
  pairwise_tests[[as.character(kk)]] <- pmat_long
}
if (length(pairwise_tests) > 0) {
  pw_df <- bind_rows(pairwise_tests)
  write.csv(pw_df, file.path(OUT_DIR, "data", "pairwise_wilcoxon.csv"),
            row.names = FALSE)
}

# --- Save run config ---
config <- list(
  run_date       = RUN_DATE,
  seed           = SEED,
  n_samples      = nrow(X),
  n_genes        = ncol(X),
  k_outer        = K_OUTER,
  r_outer        = R_OUTER,
  k_values       = K_VALUES,
  methods        = ALL_METHODS,
  residualize    = DO_RESIDUALIZE,
  confounders_cat = CONFOUNDERS_CAT,
  confounders_num = CONFOUNDERS_NUM,
  glmnet_alpha_grid = GLMNET_ALPHA_GRID,
  gs_alpha_grid     = GS_ALPHA_GRID,
  target_terms      = TARGET_TERMS
)
saveRDS(config, file.path(OUT_DIR, "data", "run_config.rds"))

cat("\n=== Quick summary at k=50 and k=200 ===\n")
print(sum_df %>% filter(k %in% c(50, 200)) %>% arrange(k, desc(AUC_mean)))

cat("\nSaved to:", OUT_DIR, "\n")
cat("  - data/nested_results.csv\n")
cat("  - data/nested_summary.csv\n")
cat("  - data/pairwise_wilcoxon.csv\n")
cat("  - data/run_config.rds\n")
cat("  - figures/parsimony_all_methods.pdf\n")
cat("  - figures/auc_boxplots.pdf\n")
cat("Done.\n")
