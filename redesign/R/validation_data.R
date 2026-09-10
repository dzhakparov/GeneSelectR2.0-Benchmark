# ==============================================================================
#  Validation-dataset specs, loaders, fold builders, residualisation.
#
#  Five datasets. GSE101794 replaces GSE57945 because it provides gene-level
#  TPM for the RISK ileal cohort. The local files do not establish that the
#  labels implement the GSE57945 corrigendum, so no correction claim is made.
#
#  Scale branches:
#    "counts"     -> train-fitted TMM + logCPM (normalise_count_split)
#    "normalized" -> used as-is (log-scale array / log2(TPM+1))
#
#  The counts branch STOPS on non-integer values: TMM on RPKM/TPM is wrong,
#  and nothing downstream would complain about it. That loud stop is the
#  TPM-in-counts-clothing trap.
# ==============================================================================

validation_specs <- list(
  GSE65682 = list(
    accession   = "GSE65682",
    label       = "Sepsis, MARS (blood) - 28-day mortality",
    scale       = "normalized",
    expr_file   = "data/GSE65682/expression_prepared.csv",
    meta_file   = "data/GSE65682/metadata_prepared.csv",
    positive    = "non_survivor",
    confounders = c("ch_gender"),
    paired      = FALSE
  ),
  GSE69683 = list(
    accession   = "GSE69683",
    label       = "Asthma, U-BIOPRED (blood) - severe vs moderate",
    scale       = "normalized",
    expr_file   = "data/GSE69683/expression_prepared.csv",
    meta_file   = "data/GSE69683/metadata_prepared.csv",
    positive    = "severe",
    confounders = c("ch_gender"),
    paired      = FALSE
  ),
  GSE13355 = list(
    accession   = "GSE13355",
    label       = "Psoriasis (skin) - lesional vs uninvolved [PAIRED, control]",
    scale       = "normalized",
    expr_file   = "data/GSE13355/expression_prepared.csv",
    meta_file   = "data/GSE13355/metadata_prepared.csv",
    positive    = "lesional",
    confounders = character(0),
    paired      = TRUE
  ),
  GSE107994 = list(
    accession   = "GSE107994",
    label       = "Tuberculosis, Leicester (blood) - active vs latent",
    scale       = "counts",
    expr_file   = "data/GSE107994/counts_prepared.csv",
    meta_file   = "data/GSE107994/metadata_prepared.csv",
    positive    = "active",
    confounders = c("ch_gender", "ch_ethnicity"),
    paired      = FALSE
  ),
  GSE101794 = list(
    accession   = "GSE101794",
    label       = "Crohn's, RISK ileum - CD vs non-IBD",
    scale       = "normalized",
    expr_file   = "data/GSE101794/expression.tsv",
    meta_file   = "data/GSE101794/metadata.tsv",
    positive    = "CD",
    confounders = character(0),
    paired      = FALSE
  ),
  GSE16879 = list(
    accession   = "GSE16879",
    label       = "IBD, infliximab response (pre-treatment mucosa)",
    scale       = "normalized",
    expr_file   = "data/GSE16879/expression_prepared.csv",
    meta_file   = "data/GSE16879/metadata_prepared.csv",
    positive    = "responder",
    confounders = c("ch_tissue"),
    paired      = FALSE
  ),
  GSE91061 = list(
    accession   = "GSE91061",
    label       = "Melanoma, anti-PD-1 - responder vs non-responder",
    scale       = "counts",
    expr_file   = "data/GSE91061/counts_prepared.csv",
    meta_file   = "data/GSE91061/metadata_prepared.csv",
    positive    = "responder",
    confounders = character(0),
    paired      = FALSE
  ),
  GSE92415 = list(
    accession   = "GSE92415",
    label       = "Ulcerative colitis, golimumab - week-6 response",
    scale       = "normalized",
    expr_file   = "data/GSE92415/expression_prepared.csv",
    meta_file   = "data/GSE92415/metadata_prepared.csv",
    positive    = "responder",
    confounders = character(0),
    paired      = FALSE
  ),
  GSE206285 = list(
    accession   = "GSE206285",
    label       = "Ulcerative colitis, ustekinumab - week-8 remission",
    scale       = "normalized",
    expr_file   = "data/GSE206285/expression_prepared.csv",
    meta_file   = "data/GSE206285/metadata_prepared.csv",
    positive    = "remission",
    confounders = character(0),
    paired      = FALSE
  )
)

#  Stop loudly if a matrix declared as counts is not non-negative integers.
assert_counts_matrix <- function(mat, accession) {
  if (any(!is.finite(mat)) || any(mat < 0)) {
    stop(sprintf(paste0(
      "%s: count matrix has negative or non-finite values. This is the ",
      "TPM-in-counts-clothing trap; TMM would be wrong, not merely ",
      "suboptimal."), accession))
  }
  frac <- sum(mat != round(mat))
  if (frac > 0) {
    stop(sprintf(paste0(
      "%s: %d count values are fractional. The file is not raw counts ",
      "(likely RPKM/TPM). Route it through the normalized branch instead; ",
      "TMM-normalising it would be wrong."), accession, frac))
  }
  invisible(TRUE)
}

#  Load one dataset -> list(expr samples x genes, outcome factor (negative,
#  positive), groups, confounder data.frame, spec).
load_validation_dataset <- function(accession, specs = validation_specs) {
  stopifnot(accession %in% names(specs))
  spec <- specs[[accession]]

  sep <- if (grepl("\\.tsv$", spec$expr_file)) "\t" else ","
  expr_raw <- read.delim(spec$expr_file, sep = sep, row.names = 1,
                         check.names = FALSE)
  expr <- t(as.matrix(expr_raw))                 # samples x genes
  storage.mode(expr) <- "numeric"

  msep <- if (grepl("\\.tsv$", spec$meta_file)) "\t" else ","
  meta <- read.delim(spec$meta_file, sep = msep, check.names = FALSE,
                     stringsAsFactors = FALSE)
  if (!file.exists(spec$meta_file) || nrow(meta) == 0)
    stop(accession, ": metadata file missing or empty")

  #  Sample id column: prepared files use sample_id; GSE101794 uses sample.
  id_col <- if ("sample_id" %in% colnames(meta)) "sample_id" else "sample"
  meta <- meta[match(rownames(expr), meta[[id_col]]), , drop = FALSE]
  if (any(is.na(meta[[id_col]])))
    stop(accession, ": expression samples do not all appear in metadata")

  #  Outcome column: prepared files use `outcome`; GSE101794 uses diagnosis.
  if ("outcome" %in% colnames(meta)) {
    y_raw <- meta$outcome
  } else {
    y_raw <- meta$diagnosis
  }
  negative <- setdiff(unique(y_raw), spec$positive)
  if (length(negative) != 1)
    stop(accession, ": expected exactly one negative level, got: ",
         paste(negative, collapse = ", "))
  outcome <- factor(y_raw, levels = c(negative, spec$positive))

  groups <- if ("group_id" %in% colnames(meta)) meta$group_id else
    rownames(expr)

  conf <- if (length(spec$confounders) > 0)
    meta[, spec$confounders, drop = FALSE] else
    data.frame(row.names = rownames(expr))

  if (spec$scale == "counts") assert_counts_matrix(expr, accession)

  list(expr = expr, outcome = outcome, groups = groups,
       confounders = conf, spec = spec)
}

#  Residualise categorical confounders, FIT ON TRAIN ONLY, applied to test.
#  Same construction as sosall_residualise_split but column-agnostic and with
#  one combined design matrix across confounders.
residualise_split_generic <- function(train_expr, test_expr,
                                      conf_train, conf_test) {
  if (ncol(conf_train) == 0) return(list(train = train_expr, test = test_expr))

  build_design <- function(conf, reference_levels = NULL) {
    #  The intercept is part of the fitted confounder model. Omitting it makes
    #  the reference group prediction zero and leaves the expression baseline
    #  in every residual. With non-zero expression scales this can increase,
    #  rather than remove, apparent between-group differences.
    blocks <- list("(Intercept)" = matrix(1, nrow = nrow(conf), ncol = 1,
                                           dimnames = list(NULL,
                                                           "(Intercept)")))
    levels_out <- list()
    for (cn in colnames(conf)) {
      vals <- as.character(conf[[cn]])
      vals[is.na(vals) | vals == ""] <- "NA_level"
      lvls <- if (is.null(reference_levels)) levels(factor(vals)) else
        reference_levels[[cn]]
      #  Single-level subsets contribute no contrast columns.
      if (length(lvls) < 2) {
        levels_out[[cn]] <- lvls
        next
      }
      #  Manual dummy matrix: rows whose level is absent from lvls (test-only
      #  levels) get all zeros = reference level. model.matrix() is avoided
      #  because its na.action silently drops those rows.
      mm <- vapply(lvls[-1], function(l) as.numeric(vals == l),
                   numeric(length(vals)))
      colnames(mm) <- paste0(cn, "_", lvls[-1])
      blocks[[cn]] <- mm
      levels_out[[cn]] <- lvls
    }
    design <- do.call(cbind, blocks)
    list(design = design, levels = levels_out)
  }

  tr <- build_design(conf_train)
  coefs <- qr.coef(qr(tr$design), train_expr)
  coefs[is.na(coefs)] <- 0                     # aliased columns -> no effect
  train_resid <- train_expr - tr$design %*% coefs

  test_resid <- NULL
  if (!is.null(test_expr)) {
    te <- build_design(conf_test, reference_levels = tr$levels)
    test_resid <- test_expr - te$design %*% coefs
  }
  list(train = train_resid, test = test_resid)
}

#  Paired folds: whole groups (patients) move together. Class balance is
#  automatic when every group spans both classes (psoriasis PP/PN pairs).
make_paired_folds <- function(groups, k_folds, seed) {
  set.seed(seed)
  ug <- sample(unique(groups))
  fold_id <- rep(seq_len(k_folds), length.out = length(ug))
  names(fold_id) <- ug
  folds <- lapply(seq_len(k_folds), function(f)
    unname(which(fold_id[as.character(groups)] == f)))
  folds
}

#  Per-split fold job list. Unpaired datasets use the same stratified folds
#  and seeds as the incumbents (seed = 42 + 1000 * repeat_idx); paired ones
#  use make_paired_folds with the same seeds.
validation_split_jobs <- function(dat, k_folds = 5, n_repeats = 3,
                                  random_seed = 42) {
  jobs <- list()
  for (rep in seq_len(n_repeats)) {
    if (dat$spec$paired) {
      folds <- make_paired_folds(dat$groups, k_folds,
                                 seed = random_seed + 1000 * rep)
    } else {
      folds <- make_stratified_folds(dat$outcome, k_folds,
                                     seed = random_seed + 1000 * rep)
    }
    for (f in seq_len(k_folds)) {
      jobs[[length(jobs) + 1]] <- list(repeat_idx = rep, fold_idx = f,
                                       test_indices = folds[[f]])
    }
  }
  jobs
}
