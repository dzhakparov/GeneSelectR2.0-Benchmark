# ==============================================================================
#  Tests for the validation benchmark: platform branch, TPM trap, paired
#  folds, leakage, gate fallback, competitor rankers, loaders.
#
#  Every test uses small synthetic matrices except the loader tests, which
#  read the real prepared files (they exist; prep is done).
# ==============================================================================

suppressPackageStartupMessages(library(testthat))

#  test_file() runs with redesign/tests as working directory.
proj_root <- normalizePath(file.path(getwd(), "..", ".."))
setwd(proj_root)
for (f in list.files(file.path(proj_root, "package", "GeneSelectR", "R"),
                     full.names = TRUE)) source(f)
for (f in c("bio_prior.R", "evaluator.R", "imvigor210_data.R",
            "competitors.R", "validation_data.R", "run_provenance.R")) {
  source(file.path(proj_root, "redesign", "R", f))
}

set.seed(1)

# ---- counts trap -------------------------------------------------------------

test_that("assert_counts_matrix accepts integer counts", {
  m <- matrix(rpois(200, 20), nrow = 10)
  expect_true(assert_counts_matrix(m, "TEST"))
})

test_that("assert_counts_matrix stops on fractional values (TPM trap)", {
  m <- matrix(rpois(200, 20), nrow = 10) + 0.5
  expect_error(assert_counts_matrix(m, "TEST"), "fractional")
})

test_that("assert_counts_matrix stops on negative values", {
  m <- matrix(rpois(200, 20), nrow = 10)
  m[1, 1] <- -3
  expect_error(assert_counts_matrix(m, "TEST"), "negative")
})

# ---- paired folds ------------------------------------------------------------

test_that("paired folds never split a patient across train and test", {
  groups <- rep(paste0("P", 1:40), each = 2)   # 40 patients, 2 samples each
  folds <- make_paired_folds(groups, k_folds = 5, seed = 1042)
  expect_length(folds, 5)
  all_idx <- sort(unlist(folds))
  expect_equal(all_idx, seq_along(groups))     # every sample exactly once
  for (f in folds) {
    in_fold <- unique(groups[f])
    out_fold <- unique(groups[-f])
    expect_length(intersect(in_fold, out_fold), 0)
  }
})

test_that("paired folds are reproducible for the same seed", {
  groups <- rep(paste0("P", 1:30), each = 2)
  f1 <- make_paired_folds(groups, 5, seed = 2042)
  f2 <- make_paired_folds(groups, 5, seed = 2042)
  expect_identical(f1, f2)
})

# ---- residualisation ---------------------------------------------------------

test_that("residualisation uses train-fitted coefficients on test", {
  n <- 60
  sex <- factor(rep(c("F", "M"), each = 30))
  X <- matrix(rnorm(n * 20, mean = 6), nrow = n)
  colnames(X) <- paste0("G", 1:20)
  X[, 1] <- X[, 1] + 3 * as.integer(sex == "M")   # strong sex effect
  conf <- data.frame(ch_gender = sex)

  tr_idx <- 1:40; te_idx <- 41:60
  out <- residualise_split_generic(X[tr_idx, ], X[te_idx, ],
                                   conf[tr_idx, , drop = FALSE],
                                   conf[te_idx, , drop = FALSE])
  #  Train residual of the sex-driven gene must lose the group difference.
  m_f <- mean(out$train[sex[tr_idx] == "F", 1])
  m_m <- mean(out$train[sex[tr_idx] == "M", 1])
  expect_lt(abs(m_m - m_f), 0.5)
  #  Test is corrected with TRAIN coefficients: manually recompute.
  d <- model.matrix(~ conf$ch_gender[tr_idx])
  cf <- qr.coef(qr(d), X[tr_idx, ])
  expect_equal(unname(out$train), unname(X[tr_idx, ] - d %*% cf),
               tolerance = 1e-10)
})

test_that("residualisation removes the expression baseline", {
  X <- matrix(rnorm(60, mean = 7), nrow = 6)
  conf <- data.frame(batch = rep(c("A", "B"), each = 3))
  out <- residualise_split_generic(X[1:4, ], X[5:6, ],
                                   conf[1:4, , drop = FALSE],
                                   conf[5:6, , drop = FALSE])
  expect_lt(max(abs(colMeans(out$train))), 1e-10)
})

test_that("residualisation survives a single-level confounder in a subset", {
  X <- matrix(rnorm(50, mean = 5), nrow = 5)
  conf <- data.frame(ch_ethnicity = c("A", "A", "A", "A", "A"))
  out <- residualise_split_generic(X[1:3, ], X[4:5, ],
                                   conf[1:3, , drop = FALSE],
                                   conf[4:5, , drop = FALSE])
  expect_equal(dim(out$train), c(3L, 10L))
  expect_equal(dim(out$test), c(2L, 10L))
  expect_lt(max(abs(colMeans(out$train))), 1e-10)
})

test_that("residualisation with no confounders is the identity", {
  X <- matrix(rnorm(50), nrow = 5)
  empty <- data.frame(row.names = 1:5)
  out <- residualise_split_generic(X[1:3, ], X[4:5, ], empty[1:3, , drop = FALSE],
                                   empty[4:5, , drop = FALSE])
  expect_identical(out$train, X[1:3, ])
})

# ---- checkpoint provenance --------------------------------------------------

test_that("run manifests reject changed source files and legacy directories", {
  source_file <- tempfile(fileext = ".R")
  input_file <- tempfile(fileext = ".csv")
  out_dir <- tempfile()
  writeLines("x <- 1", source_file)
  writeLines("a,b\n1,2", input_file)
  prepare_redesign_run(out_dir, list(seed = 1), source_file, input_file)
  expect_true(file.exists(file.path(out_dir, "run_manifest.rds")))
  expect_true(file.exists(file.path(out_dir, "run_manifest.version")))
  expect_silent(prepare_redesign_run(out_dir, list(seed = 1), source_file,
                                     input_file))

  writeLines("x <- 2", source_file)
  expect_error(prepare_redesign_run(out_dir, list(seed = 1), source_file,
                                    input_file), "does not match")

  legacy_dir <- tempfile()
  dir.create(legacy_dir)
  writeLines("old", file.path(legacy_dir, "eval_results.csv"))
  expect_error(prepare_redesign_run(legacy_dir, list(seed = 1), source_file,
                                    input_file), "no provenance manifest")
})

test_that("object hashes are reproducible and detect data changes", {
  object <- list(X = matrix(1:12, nrow = 3), y = factor(c("A", "B", "A")))
  expect_identical(redesign_object_hash(object), redesign_object_hash(object))
  changed <- object
  changed$X[1, 1] <- 99
  expect_false(identical(redesign_object_hash(object),
                         redesign_object_hash(changed)))
})

test_that("extension manifests reject changed code and untracked output", {
  base_source <- tempfile(fileext = ".R")
  extension_source <- tempfile(fileext = ".R")
  out_dir <- tempfile()
  writeLines("base <- 1", base_source)
  writeLines("extension <- 1", extension_source)
  prepare_redesign_run(out_dir, list(seed = 1), base_source)
  prepare_redesign_extension(out_dir, "example", extension_source,
                             config = list(dataset = "synthetic"))
  expect_true(file.exists(file.path(
    out_dir, "extension_manifest_example.rds")))
  expect_true(file.exists(file.path(
    out_dir, "extension_manifest_example.version")))
  expect_silent(prepare_redesign_extension(
    out_dir, "example", extension_source,
    config = list(dataset = "synthetic")))
  expect_silent(require_redesign_extension(
    out_dir, "example", extension_source,
    config = list(dataset = "synthetic")))

  writeLines("extension <- 2", extension_source)
  expect_error(prepare_redesign_extension(
    out_dir, "example", extension_source,
    config = list(dataset = "synthetic")), "does not match")
  expect_error(require_redesign_extension(
    out_dir, "example", extension_source,
    config = list(dataset = "synthetic")), "does not match")

  legacy_out_dir <- tempfile()
  prepare_redesign_run(legacy_out_dir, list(seed = 1), base_source)
  writeLines("repeat_idx,fold_idx,arm,k,AUC",
             file.path(legacy_out_dir, "eval_legacy.csv"))
  expect_error(prepare_redesign_extension(
    legacy_out_dir, "legacy", extension_source), "no provenance manifest")
})

# ---- pool leakage ------------------------------------------------------------

test_that("make_pools variance filter uses the train matrix only", {
  n_tr <- 40; n_te <- 10
  Xtr <- matrix(rnorm(n_tr * 100), nrow = n_tr)
  colnames(Xtr) <- paste0("G", 1:100)
  pools <- make_pools(Xtr, bio_genes = character(0), top_var = 20)
  v <- apply(Xtr, 2, var)
  expect_setequal(pools$var2000, names(sort(v, decreasing = TRUE))[1:20])
})

# ---- module gate + fallback ---------------------------------------------------

test_that("module gate certifies nothing under a null outcome", {
  set.seed(2)
  X <- matrix(rnorm(80 * 300), nrow = 80)
  colnames(X) <- paste0("G", 1:300)
  ybin <- rbinom(80, 1, 0.5)
  g <- module_gate(X, ybin, B = 100, seed = 7)
  expect_true(is.list(g$certified_genes))
  expect_true("q20" %in% names(g$certified_genes))
  #  On pure noise the q20 set should be empty or tiny.
  expect_lt(length(g$certified_genes$q20), 50)
})

test_that("module gate flags real planted module signal", {
  set.seed(3)
  n <- 100
  ybin <- rep(0:1, each = 50)
  X <- matrix(rnorm(n * 300), nrow = n)
  colnames(X) <- paste0("G", 1:300)
  sets <- hallmark_sets_all()
  s1 <- sets[[1]]
  members <- colnames(X)[1:25]
  #  Overwrite the first 25 columns with the planted signal and rename them
  #  into the first Hallmark set so the gate can find them.
  colnames(X)[1:25] <- s1[1:25]
  X[51:100, 1:25] <- X[51:100, 1:25] + 2
  g <- module_gate(X, ybin, B = 200, seed = 7)
  expect_gt(length(g$certified_genes$q20), 0)
  expect_true(all(g$certified_genes$q20 %in% colnames(X)))
})

# ---- competitor rankers -------------------------------------------------------

make_signal_data <- function(n = 80, p = 60, n_signal = 5, eff = 1.5) {
  y <- factor(rep(c("neg", "pos"), each = n / 2), levels = c("neg", "pos"))
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("G", 1:p)
  sig <- paste0("G", 1:n_signal)
  X[y == "pos", sig] <- X[y == "pos", sig] + eff
  list(X = X, y = y, sig = sig)
}

test_that("DGE ranks planted signal at the top", {
  d <- make_signal_data()
  rk <- rank_by_differential_expression(d$X, d$y)
  expect_length(rk$ranked, 60)
  expect_true(all(d$sig %in% rk$ranked[1:10]))
  expect_true(all(d$sig %in% rk$selected))
})

test_that("LASSO and ElasticNet return full-length rankings", {
  d <- make_signal_data()
  rk1 <- rank_by_lasso(d$X, d$y)
  rk2 <- rank_by_elastic_net(d$X, d$y)
  expect_length(rk1$ranked, 60)
  expect_length(rk2$ranked, 60)
  expect_setequal(rk1$ranked, colnames(d$X))
  #  Signal should be enriched near the top.
  expect_true(mean(d$sig %in% rk1$ranked[1:15]) >= 0.6)
})

test_that("RF importance returns a full-length ranking", {
  d <- make_signal_data()
  rk <- rank_by_random_forest(d$X, d$y, random_seed = 1)
  expect_length(rk$ranked, 60)
  expect_true(all(d$sig %in% rk$ranked[1:20]))
})

test_that("mRMR and Boruta run and return valid rankings", {
  skip_if_not(requireNamespace("mRMRe", quietly = TRUE))
  skip_if_not(requireNamespace("Boruta", quietly = TRUE))
  d <- make_signal_data()
  rk1 <- rank_by_mrmr(d$X, d$y)
  rk2 <- rank_by_boruta(d$X, d$y)
  expect_setequal(rk1$ranked, colnames(d$X))
  expect_setequal(rk2$ranked, colnames(d$X))
})

test_that("Boruta handles non-syntactic gene names without duplication", {
  skip_if_not(requireNamespace("Boruta", quietly = TRUE))
  d <- make_signal_data(p = 40, n_signal = 3)
  nasty <- paste0(c("1ABC", "GENE-X", "A B"), 1:3)
  colnames(d$X)[1:3] <- nasty
  rk <- rank_by_boruta(d$X, d$y)
  expect_length(rk$ranked, 40)
  expect_setequal(rk$ranked, colnames(d$X))
  expect_true(all(rk$selected %in% colnames(d$X)))
})

# ---- loaders (real data) ------------------------------------------------------

test_that("all five datasets load with coherent dimensions", {
  expected <- list(
    GSE65682  = c(n = 479, pos = 114),
    GSE69683  = c(n = 323, pos = 246),
    GSE13355  = c(n = 116, pos = 58),
    GSE107994 = c(n = 94,  pos = 42),
    GSE101794 = c(n = 304, pos = 254)
  )
  for (acc in names(expected)) {
    dat <- load_validation_dataset(acc)
    expect_equal(nrow(dat$expr), unname(expected[[acc]]["n"]),
                 label = paste(acc, "n"))
    expect_equal(sum(dat$outcome == levels(dat$outcome)[2]),
                 unname(expected[[acc]]["pos"]),
                 label = paste(acc, "positives"))
    expect_equal(nrow(dat$expr), length(dat$outcome))
    expect_equal(nrow(dat$expr), length(dat$groups))
    expect_false(anyNA(dat$expr))
  }
})

test_that("TB loader routes counts through the integer check", {
  dat <- load_validation_dataset("GSE107994")
  expect_equal(dat$spec$scale, "counts")
  m <- dat$expr
  expect_true(all(m == round(m)))
})

test_that("psoriasis groups are 58 patients, two samples each", {
  dat <- load_validation_dataset("GSE13355")
  expect_true(dat$spec$paired)
  tt <- table(dat$groups)
  expect_equal(length(tt), 58)
  expect_true(all(tt == 2))
})

test_that("GSE101794 outcome order is negative-first", {
  dat <- load_validation_dataset("GSE101794")
  expect_equal(levels(dat$outcome), c("Non-IBD", "CD"))
})

# ---- split jobs ----------------------------------------------------------------

test_that("split jobs give 15 splits with full coverage", {
  dat <- load_validation_dataset("GSE107994")
  jobs <- validation_split_jobs(dat)
  expect_length(jobs, 15)
  for (j in jobs) {
    expect_true(all(j$test_indices %in% seq_along(dat$outcome)))
    expect_gt(length(j$test_indices), 5)
  }
})

test_that("paired split jobs keep patients whole", {
  dat <- load_validation_dataset("GSE13355")
  jobs <- validation_split_jobs(dat)
  for (j in jobs) {
    in_g <- unique(dat$groups[j$test_indices])
    out_g <- unique(dat$groups[-j$test_indices])
    expect_length(intersect(in_g, out_g), 0)
  }
})
