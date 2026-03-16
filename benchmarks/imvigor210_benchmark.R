# ==============================================================================
# GeneSelectR 2.0 — Benchmark on Mariathasan et al. 2018 (IMvigor210)
#
# Dataset: IMvigor210 — anti-PD-L1 (atezolizumab) in metastatic urothelial CA
#   Binary outcome: Responder (CR/PR) vs Non-responder (SD/PD)
#
# DATA ACCESS (two options, tried in order):
#
#   Option A — easierData (Bioconductor, maintained, 192 patients):
#     BiocManager::install("easierData")
#
#   Option B — Full cohort RDS from GitHub (~348 patients):
#     download.file(
#       "https://github.com/snijeshvp/IMvigor210/raw/main/IMvigor210.all.rds",
#       destfile = "data/IMvigor210.all.rds"
#     )
#
# The script auto-detects which source is available.
# ==============================================================================

suppressPackageStartupMessages({
  library(GeneSelectR)
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

pkg_ensure <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s ...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

pkg_ensure("mRMRe")
pkg_ensure("Boruta")
pkg_ensure("randomForest")
pkg_ensure("ranger")
pkg_ensure("edgeR")

# ----------------------------- CONFIG -----------------------------------------
SEED <- 42
set.seed(SEED)

RESPONSE_COL   <- "BOR"             # Best Overall Response
RESPONDER_VALS <- c("R")
NONRESP_VALS   <- c("NR")

CONFOUNDERS_CAT <- character(0)      # adjust if metadata available
CONFOUNDERS_NUM <- character(0)
DO_RESIDUALIZE  <- FALSE             # set TRUE once confounders are confirmed

K_OUTER  <- 5
R_OUTER  <- 5
K_VALUES <- c(10, 20, 50, 100, 200, 500)

K_INNER_GS <- 5
R_INNER_GS <- 5
N_CORES <- max(1, parallel::detectCores() - 1)

TOP_VAR_GENES <- 5000

# Targeted GO terms: immune + TGFbeta (Mariathasan et al. biology)
TARGET_TERMS <- c(
  "GO:0006955",   # immune response
  "GO:0002376",   # immune system process
  "GO:0045087",   # innate immune response
  "GO:0002250",   # adaptive immune response
  "GO:0042110",   # T cell activation
  "GO:0050863",   # regulation of T cell activation
  "GO:0007179",   # TGFbeta receptor signaling pathway
  "GO:0071559",   # response to TGFbeta
  "GO:0002682",   # regulation of immune system process
  "GO:0050776"    # regulation of immune response
)
BIO_ONTOLOGY    <- "BP"
BIO_SIM_METHOD  <- "resnik"
BIO_ENRICH_FDR  <- 0.05
BIO_IC_QUANTILE <- 0.5
BIO_MAX_ENRICHED <- 100
BIO_TOPK_SIMS   <- 5
SCORE_FORMULA <- "geometric"
SCORE_WEIGHTS <- c(1, 1, 1)

GS_ALPHA_GRID    <- c(0.5, 1.0)
GS_WEIGHTS_GRID  <- list(c(1, 1, 1))
GLMNET_ALPHA_GRID <- c(0.0, 0.25, 0.5, 0.75, 1.0)
BORUTA_NTREE_GRID <- c(500)
BORUTA_MAXRUNS    <- c(100)
RF_NTREE_GRID     <- c(1000)
RF_MTRY_FRAC      <- c(0.33)
MRMR_NFEATURES_GRID <- c(200)

RUN_DATE <- format(Sys.Date(), "%Y-%m-%d")
OUT_DIR  <- file.path("results_cancer", RUN_DATE)
dir.create(file.path(OUT_DIR, "data"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== IMvigor210 benchmark | Run date: %s ===\n", RUN_DATE))

# ========================== DATA LOADING ======================================

raw_counts <- NULL
pheno      <- NULL
data_source <- "none"

# ---------- Option A: easierData (Bioconductor) ----------
if (requireNamespace("easierData", quietly = TRUE)) {
  cat("Loading IMvigor210 via easierData (Bioconductor)...\n")
  suppressPackageStartupMessages({
    library(ExperimentHub)
    library(easierData)
    library(SummarizedExperiment)
  })

  eh <- ExperimentHub()
  dat <- eh[["EH6677"]]  # Mariathasan2018_PDL1_treatment SummarizedExperiment

  # Extract counts and metadata
  raw_counts <- assay(dat, "counts")
  pheno      <- as.data.frame(colData(dat))

  # Map the BOR column
  if (!"BOR" %in% names(pheno)) {
    # easierData may use different column names; check alternatives
    bor_candidates <- c("BOR", "Best.Confirmed.Overall.Response",
                        "best_confirmed_overall_response", "binaryResponse")
    found <- intersect(bor_candidates, names(pheno))
    if (length(found) > 0) {
      RESPONSE_COL <- found[1]
      cat(sprintf("  Response column: %s\n", RESPONSE_COL))
    } else {
      cat("  Available columns:", paste(names(pheno), collapse = ", "), "\n")
      stop("Cannot find response column in easierData. Check colData(dat).")
    }
  }
  data_source <- "easierData"
  cat(sprintf("  Loaded: %d genes x %d samples\n", nrow(raw_counts), ncol(raw_counts)))
}

# ---------- Option B: GitHub RDS (full cohort) ----------
if (is.null(raw_counts)) {
  rds_path <- "data/IMvigor210.all.rds"
  if (!file.exists(rds_path)) {
    cat("easierData not installed. Attempting GitHub download...\n")
    dir.create("data", showWarnings = FALSE)
    tryCatch({
      download.file(
        "https://github.com/snijeshvp/IMvigor210/raw/main/IMvigor210.all.rds",
        destfile = rds_path, mode = "wb", quiet = FALSE
      )
    }, error = function(e) {
      stop(
        "Could not download IMvigor210 data. Please install one of:\n",
        "  BiocManager::install('easierData')          # 192 patients, Bioconductor\n",
        "  # OR manually download the RDS to data/IMvigor210.all.rds from:\n",
        "  # https://github.com/snijeshvp/IMvigor210\n",
        call. = FALSE
      )
    })
  }

  cat("Loading IMvigor210 from GitHub RDS...\n")
  imv <- readRDS(rds_path)

  # The RDS contains a list with: rawcounts, clinical, gene_info, etc.
  if (is.list(imv) && "rawcounts" %in% names(imv)) {
    raw_counts <- as.matrix(imv$rawcounts)
    pheno      <- as.data.frame(imv$clinical)
    data_source <- "github_rds"
  } else if (is.list(imv) && "counts" %in% names(imv)) {
    raw_counts <- as.matrix(imv$counts)
    pheno      <- as.data.frame(imv$clinical)
    data_source <- "github_rds"
  } else {
    stop("Unexpected RDS structure. Expected list with 'rawcounts' and 'clinical'.")
  }

  # Map response column
  bor_candidates <- c("Best.Confirmed.Overall.Response", "BOR",
                      "best_confirmed_overall_response", "binaryResponse")
  found <- intersect(bor_candidates, names(pheno))
  if (length(found) > 0) {
    RESPONSE_COL <- found[1]
  } else {
    cat("  Available columns:", paste(head(names(pheno), 20), collapse = ", "), "\n")
    stop("Cannot find response column in GitHub RDS. Check names(imv$clinical).")
  }
  cat(sprintf("  Loaded: %d genes x %d samples | Response col: %s\n",
              nrow(raw_counts), ncol(raw_counts), RESPONSE_COL))
}

cat(sprintf("Data source: %s\n", data_source))
cat(sprintf("Raw data: %d genes x %d samples\n", nrow(raw_counts), ncol(raw_counts)))

# ========================== PREPROCESSING =====================================

# ---- Filter to evaluable response ----
resp_vals <- as.character(pheno[[RESPONSE_COL]])
valid <- resp_vals %in% c(RESPONDER_VALS, NONRESP_VALS)
# vslifcat(sprintf("Response distribution (all):\n"))
print(table(resp_vals, useNA = "ifany"))

if (sum(valid) < 30) stop("Too few evaluable samples")

raw_counts <- raw_counts[, valid, drop = FALSE]
pheno      <- pheno[valid, , drop = FALSE]

cat(sprintf("After response filter: %d samples\n", ncol(raw_counts)))

# ---- Binary outcome ----
pheno$response_binary <- ifelse(
  as.character(pheno[[RESPONSE_COL]]) %in% RESPONDER_VALS,
  "Responder", "NonResponder"
)
cat("Binary outcome:\n"); print(table(pheno$response_binary))
OUTCOME_COL <- "response_binary"

# ---- edgeR normalization → log2 CPM ----
dge <- edgeR::DGEList(counts = raw_counts)
min_samples <- max(5, round(ncol(raw_counts) * 0.10))
keep <- rowSums(edgeR::cpm(dge) > 1) >= min_samples
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

cat(sprintf("After expression filter: %d genes\n", nrow(logcpm)))

# ---- Gene symbol mapping (if Ensembl IDs) ----
gene_ids <- rownames(logcpm)

if (any(grepl("^ENSG", gene_ids))) {
  cat("Detected Ensembl IDs — mapping to gene symbols...\n")
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    mapping <- tryCatch({
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = sub("\\..*", "", gene_ids),  # strip version suffix
        column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
      )
    }, error = function(e) NULL)
    if (!is.null(mapping)) {
      mapped <- mapping[sub("\\..*", "", gene_ids)]
      has_symbol <- !is.na(mapped)
      logcpm <- logcpm[has_symbol, ]
      rownames(logcpm) <- mapped[has_symbol]
      cat(sprintf("  Mapped %d / %d to symbols\n", sum(has_symbol), length(has_symbol)))
    }
  } else {
    cat("  org.Hs.eg.db not installed; using raw IDs.\n")
  }
}

# ---- Aggregate duplicate gene symbols ----
if (any(duplicated(rownames(logcpm)))) {
  cat("Aggregating duplicated gene symbols...\n")
  expr_df <- as.data.frame(logcpm)
  expr_df$gene <- rownames(logcpm)
  expr_agg <- aggregate(. ~ gene, data = expr_df, FUN = mean)
  rownames(expr_agg) <- expr_agg$gene
  expr_agg$gene <- NULL
  logcpm <- as.matrix(expr_agg)
}

# ---- Build X (samples x genes) ----
X <- t(logcpm)
storage.mode(X) <- "numeric"
X[is.na(X)] <- 0

common <- intersect(rownames(X), rownames(pheno))
if (length(common) < 30) stop("Too few matched samples")
X    <- X[common, , drop = FALSE]
meta <- pheno[common, , drop = FALSE]

y <- factor(meta[[OUTCOME_COL]])
y <- droplevels(y)
if (nlevels(y) != 2) stop("Outcome must have exactly 2 classes")

cat(sprintf("\nFinal dataset: %d samples x %d genes\n", nrow(X), ncol(X)))
cat("Outcome:\n"); print(table(y))

# ---- Variance filter ----
vars <- apply(X, 2, var)
X <- X[, vars > 0, drop = FALSE]
if (!is.null(TOP_VAR_GENES) && ncol(X) > TOP_VAR_GENES) {
  top_idx <- order(apply(X, 2, var), decreasing = TRUE)[1:TOP_VAR_GENES]
  X <- X[, top_idx, drop = FALSE]
}
cat(sprintf("After variance filter: %d genes\n", ncol(X)))

# ========================== HELPER FUNCTIONS ==================================

fast_auc <- function(y, score) {
  y <- droplevels(y)
  if (nlevels(y) != 2) return(NA_real_)
  pos <- which(y == levels(y)[2])
  neg <- which(y == levels(y)[1])
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[pos]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

make_stratified_folds <- function(y, K = 5, seed = 1) {
  set.seed(seed)
  y <- droplevels(y)
  idx1 <- which(y == levels(y)[1])
  idx2 <- which(y == levels(y)[2])
  s1 <- split(sample(idx1), rep(1:K, length.out = length(idx1)))
  s2 <- split(sample(idx2), rep(1:K, length.out = length(idx2)))
  folds <- vector("list", K)
  for (k in 1:K) folds[[k]] <- sort(c(s1[[k]], s2[[k]]))
  folds
}

predict_glmnet_tuned <- function(X_train, y_train, X_test,
                                 alpha_grid = GLMNET_ALPHA_GRID) {
  y_num <- as.integer(droplevels(y_train) == levels(droplevels(y_train))[2])
  nf <- min(5, max(3, floor(min(table(y_train)) * 0.8)))
  best_cvm <- Inf; best_fit <- NULL
  for (a in alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(x = X_train, y = y_num, family = "binomial",
                        alpha = a, nfolds = nf),
      error = function(e) NULL)
    if (!is.null(fit) && min(fit$cvm) < best_cvm) {
      best_cvm <- min(fit$cvm); best_fit <- fit
    }
  }
  if (is.null(best_fit)) return(rep(NA_real_, nrow(X_test)))
  as.numeric(predict(best_fit, X_test, s = "lambda.min", type = "response"))
}

# ========================== RANKING FUNCTIONS =================================

rank_geneselectr_targeted <- function(X_train, y_train) {
  best_genes <- NULL; best_score <- -Inf
  for (a in GS_ALPHA_GRID) {
    for (sw in GS_WEIGHTS_GRID) {
      gs <- tryCatch(
        geneselectr2_fit(
          X = X_train, y = y_train, bio_mode = "supervised",
          target_terms = TARGET_TERMS, bio_ontology = BIO_ONTOLOGY,
          bio_sim_method = BIO_SIM_METHOD, score_formula = SCORE_FORMULA,
          score_weights = sw, regularization_method = "elastic_net",
          alpha = a, K = K_INNER_GS, R = R_INNER_GS,
          n_cores = N_CORES, random_seed = SEED, verbose = FALSE
        ), error = function(e) NULL)
      if (is.null(gs)) next
      sc <- mean(gs$cv_results$auc_scores, na.rm = TRUE)
      if (sc > best_score) { best_score <- sc; best_genes <- gs$gene_scores$gene }
    }
  }
  if (is.null(best_genes)) return(sample(colnames(X_train)))
  best_genes
}

rank_geneselectr_datadriven <- function(X_train, y_train) {
  best_genes <- NULL; best_score <- -Inf
  for (a in GS_ALPHA_GRID) {
    for (sw in GS_WEIGHTS_GRID) {
      gs <- tryCatch(
        geneselectr2_fit(
          X = X_train, y = y_train, bio_mode = "data_driven",
          bio_ontology = BIO_ONTOLOGY, bio_sim_method = BIO_SIM_METHOD,
          bio_enrich_fdr = BIO_ENRICH_FDR, bio_ic_quantile = BIO_IC_QUANTILE,
          bio_max_enriched = BIO_MAX_ENRICHED, bio_n_top_sims = BIO_TOPK_SIMS,
          score_formula = SCORE_FORMULA, score_weights = sw,
          regularization_method = "elastic_net", alpha = a,
          K = K_INNER_GS, R = R_INNER_GS, n_cores = N_CORES,
          random_seed = SEED, verbose = FALSE
        ), error = function(e) NULL)
      if (is.null(gs)) next
      sc <- mean(gs$cv_results$auc_scores, na.rm = TRUE)
      if (sc > best_score) { best_score <- sc; best_genes <- gs$gene_scores$gene }
    }
  }
  if (is.null(best_genes)) return(sample(colnames(X_train)))
  best_genes
}

rank_dge_ttest <- function(X_train, y_train) {
  y_train <- droplevels(y_train)
  p <- apply(X_train, 2, function(g) {
    tryCatch(t.test(g ~ y_train)$p.value, error = function(e) 1.0)
  })
  names(sort(p.adjust(p, method = "BH")))
}

rank_lasso <- function(X_train, y_train) {
  y_num <- as.integer(droplevels(y_train) == levels(droplevels(y_train))[2])
  nf <- min(5, max(3, floor(min(table(y_train)) * 0.8)))
  fit <- tryCatch(
    glmnet::cv.glmnet(x = X_train, y = y_num, family = "binomial",
                      alpha = 1, nfolds = nf),
    error = function(e) NULL)
  if (is.null(fit)) return(sample(colnames(X_train)))
  coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
  names(coefs) <- colnames(X_train)
  names(sort(abs(coefs), decreasing = TRUE))
}

rank_mrmr <- function(X_train, y_train) {
  y_num <- as.integer(droplevels(y_train) == levels(droplevels(y_train))[2])
  best_genes <- NULL; best_score <- -Inf
  for (nfeat in MRMR_NFEATURES_GRID) {
    nfeat_use <- min(nfeat, ncol(X_train))
    dd <- tryCatch({
      df_in <- data.frame(outcome = y_num, X_train[, 1:min(ncol(X_train), 5000)])
      mRMR.data(data = df_in)
    }, error = function(e) NULL)
    if (is.null(dd)) next
    res <- tryCatch(
      mRMR.classic(data = dd, target_indices = 1, feature_count = nfeat_use),
      error = function(e) NULL)
    if (is.null(res)) next
    sel_idx <- as.integer(solutions(res)[[1]])
    sel_idx <- sel_idx[sel_idx > 0 & sel_idx <= ncol(X_train)]
    if (length(sel_idx) == 0) next
    gene_names <- colnames(X_train)[sel_idx]
    inner_fit <- tryCatch(
      glmnet::cv.glmnet(x = X_train[, gene_names, drop = FALSE], y = y_num,
                        family = "binomial", alpha = 0.5,
                        nfolds = min(5, max(3, floor(min(table(y_train)) * 0.8)))),
      error = function(e) NULL)
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


rank_rf_importance <- function(X_train, y_train) {
  y_train <- droplevels(y_train)
  max_feat <- min(ncol(X_train), 5000)
  X_sub <- if (ncol(X_train) > max_feat) {
    vr <- apply(X_train, 2, var)
    X_train[, order(vr, decreasing = TRUE)[1:max_feat], drop = FALSE]
  } else X_train
  best_genes <- NULL; best_oob <- Inf
  for (ntree in RF_NTREE_GRID) {
    for (mf in RF_MTRY_FRAC) {
      mtry_val <- max(1, round(ncol(X_sub) * mf))
      rf <- tryCatch(
        ranger::ranger(x = X_sub, y = y_train, num.trees = ntree,
                       mtry = mtry_val, importance = "impurity", seed = SEED),
        error = function(e) NULL)
      if (is.null(rf)) next
      if (rf$prediction.error < best_oob) {
        best_oob <- rf$prediction.error
        imp <- sort(rf$variable.importance, decreasing = TRUE)
        best_genes <- c(names(imp), setdiff(colnames(X_sub), names(imp)))
      }
    }
  }
  if (is.null(best_genes)) return(sample(colnames(X_train)))
  c(best_genes, setdiff(colnames(X_train), best_genes))
}

rank_random <- function(X_train, y_train) sample(colnames(X_train))

# ======================== METHOD REGISTRY =====================================
ALL_METHODS <- c(
  "GeneSelectR_targeted", "GeneSelectR_datadriven",
  "DGE_ttest", "LASSO", "mRMR", "Boruta", "RF_importance", "Random"
)

get_ranker <- function(method) {
  switch(method,
         "GeneSelectR_targeted"   = rank_geneselectr_targeted,
         "GeneSelectR_datadriven" = rank_geneselectr_datadriven,
         "DGE_ttest"   = rank_dge_ttest,
         "LASSO"       = rank_lasso,
         "mRMR"        = rank_mrmr,
         "Boruta"      = rank_boruta,
         "RF_importance" = rank_rf_importance,
         "Random"      = rank_random,
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
      X_tr <- X[train_idx, , drop = FALSE]
      X_te <- X[test_idx,  , drop = FALSE]
      y_tr <- y[train_idx]; y_te <- y[test_idx]

      # Standardize (fit on train, apply to test)
      mu <- colMeans(X_tr); sdv <- apply(X_tr, 2, sd)
      sdv[sdv == 0 | is.na(sdv)] <- 1
      X_tr <- sweep(sweep(X_tr, 2, mu, "-"), 2, sdv, "/")
      X_te <- sweep(sweep(X_te, 2, mu, "-"), 2, sdv, "/")

      cat(sprintf("  [%s] repeat=%d fold=%d ...\n", method, r, f))
      ranked <- tryCatch(rank_fn(X_tr, y_tr), error = function(e) {
        warning(sprintf("Ranking failed (%s r=%d f=%d): %s", method, r, f, e$message))
        sample(colnames(X_tr))
      })

      for (k in k_values) {
        genes_k <- head(ranked[ranked %in% colnames(X_tr)], k)
        if (length(genes_k) < 5) next
        p <- tryCatch(
          predict_glmnet_tuned(X_tr[, genes_k, drop = FALSE], y_tr,
                               X_te[, genes_k, drop = FALSE]),
          error = function(e) rep(NA_real_, length(y_te)))
        auc <- if (all(is.na(p))) NA_real_ else fast_auc(y_te, p)
        rows[[length(rows) + 1]] <- data.frame(
          Method = method, Repeat = r, Fold = f, k = k, AUC = auc,
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

# ============================== RUN ===========================================
cat(sprintf("\nRunning nested evaluation (%d methods, %dx%d CV) ...\n",
            length(ALL_METHODS), R_OUTER, K_OUTER))

all_results <- list()
for (m in ALL_METHODS) {
  cat(sprintf("\n>>> Method: %s\n", m))
  t0 <- proc.time()
  all_results[[m]] <- evaluate_nested(X, y, meta, method = m, seed = SEED)
  dt <- (proc.time() - t0)["elapsed"]
  cat(sprintf("<<< %s done in %.1f sec\n", m, dt))
}

res_df <- dplyr::bind_rows(all_results)
write.csv(res_df, file.path(OUT_DIR, "data", "nested_results.csv"), row.names = FALSE)

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
write.csv(sum_df, file.path(OUT_DIR, "data", "nested_summary.csv"), row.names = FALSE)

perf_table <- sum_df %>%
  mutate(AUC_str = sprintf("%.3f +/- %.3f", AUC_mean, AUC_sd)) %>%
  select(Method, k, AUC_str) %>%
  pivot_wider(names_from = k, values_from = AUC_str, names_prefix = "k=")
write.csv(perf_table, file.path(OUT_DIR, "data", "performance_table.csv"),
          row.names = FALSE)

method_colors <- c(
  "GeneSelectR_targeted"   = "#E41A1C",
  "GeneSelectR_datadriven" = "#377EB8",
  "DGE_ttest"              = "#4DAF4A",
  "LASSO"                  = "#FF7F00",
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
    title = "IMvigor210: Feature selection benchmark (nested CV, HP-tuned)",
    subtitle = sprintf("Run: %s | %d samples x %d genes | %dx%d outer CV | CR/PR vs SD/PD",
                       RUN_DATE, nrow(X), ncol(X), R_OUTER, K_OUTER),
    x = "Top-k genes (log scale)",
    y = "Outer-fold AUC (mean +/- SD)"
  )
ggsave(file.path(OUT_DIR, "figures", "parsimony_all_methods.pdf"),
       p, width = 12, height = 7)

p_box <- ggplot(res_df, aes(x = factor(k), y = AUC, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = method_colors) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  theme_bw(base_size = 12) +
  labs(
    title = "IMvigor210: AUC distributions by method and k",
    subtitle = sprintf("Run: %s", RUN_DATE),
    x = "Top-k genes", y = "AUC"
  )
ggsave(file.path(OUT_DIR, "figures", "auc_boxplots.pdf"),
       p_box, width = 14, height = 7)

# ---- Pairwise Wilcoxon ----
pairwise_tests <- list()
for (kk in c(50, 100, 200)) {
  sub <- res_df %>% filter(k == kk, !is.na(AUC))
  if (length(unique(sub$Method)) < 2) next
  pw <- pairwise.wilcox.test(sub$AUC, sub$Method, p.adjust.method = "BH")
  pmat <- as.data.frame(pw$p.value); pmat$Method1 <- rownames(pmat)
  pmat_long <- pmat %>%
    pivot_longer(-Method1, names_to = "Method2", values_to = "p_adj") %>%
    filter(!is.na(p_adj)) %>% mutate(k = kk)
  pairwise_tests[[as.character(kk)]] <- pmat_long
}
if (length(pairwise_tests) > 0) {
  pw_df <- bind_rows(pairwise_tests)
  write.csv(pw_df, file.path(OUT_DIR, "data", "pairwise_wilcoxon.csv"),
            row.names = FALSE)
}

# ---- Save config ----
config <- list(
  dataset     = "IMvigor210 (Mariathasan et al. 2018, Nature)",
  data_source = data_source,
  run_date    = RUN_DATE,
  seed        = SEED,
  n_samples   = nrow(X),
  n_genes     = ncol(X),
  outcome     = "Responder (CR/PR) vs Non-responder (SD/PD)",
  k_outer     = K_OUTER,
  r_outer     = R_OUTER,
  k_values    = K_VALUES,
  methods     = ALL_METHODS,
  target_terms = TARGET_TERMS
)
saveRDS(config, file.path(OUT_DIR, "data", "run_config.rds"))

cat("\n=== Summary at k=50 and k=100 ===\n")
print(sum_df %>% filter(k %in% c(50, 100)) %>% arrange(k, desc(AUC_mean)))
cat("\nAll outputs saved to:", OUT_DIR, "\n")
cat("Done.\n")
