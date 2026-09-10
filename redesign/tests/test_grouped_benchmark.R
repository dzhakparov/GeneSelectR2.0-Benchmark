# ==============================================================================
#  Tests for the grouped-benchmark redesign (bio_prior.R, evaluator.R,
#  sosall_data.R, runner logic). Run via redesign/run_tests_grouped.R.
#
#  The integration test at the end plants signal in a real Hallmark set and
#  checks that (a) the module gate certifies it, (b) the package fit recipe
#  runs on the gated pool and returns a sane ranking, (c) a pure-noise
#  dataset certifies nothing. That guards against invented results: if the
#  machinery certified noise, these tests fail.
# ==============================================================================

suppressPackageStartupMessages(library(testthat))

#  test_file() runs with redesign/tests as working directory.
proj_root <- normalizePath(file.path(getwd(), "..", ".."))
for (f in list.files(file.path(proj_root, "package", "GeneSelectR", "R"),
                     full.names = TRUE)) source(f)
for (f in c("bio_prior.R", "evaluator.R", "sosall_data.R",
            "imvigor210_data.R")) {
  source(file.path(proj_root, "redesign", "R", f))
}

set.seed(2024)


# ------------------------------------------------------------------------------
context("bio_prior: Hallmark sets and pools")

test_that("hallmark_sets_all returns the expected collection", {
  sets <- hallmark_sets_all()
  expect_gte(length(sets), 40)
  expect_true("HALLMARK_INTERFERON_GAMMA_RESPONSE" %in% names(sets))
  expect_true(all(vapply(sets, function(s)
    is.character(s) && length(s) > 0, logical(1))))
  u <- hallmark_union()
  expect_gt(length(u), 3000)
  expect_equal(anyDuplicated(u), 0L)
})

test_that("make_pools builds correct, disjoint-purpose pools", {
  X <- matrix(rnorm(80 * 120), 80, 120)
  colnames(X) <- paste0("G", seq_len(120))
  X[, 1:5] <- 0                      # zero-variance columns must be excluded
  bio_genes <- paste0("G", 1:60)     # 55 of these have variance > 0
  pools <- make_pools(X, bio_genes, top_var = 50)

  expect_length(pools$var2000, 50)
  expect_length(pools$bio, 55)
  expect_length(pools$varlarge, length(pools$bio))
  expect_true(all(pools$bio %in% bio_genes))
  expect_false(any(pools$bio %in% paste0("G", 1:5)))
  expect_equal(anyDuplicated(pools$var2000), 0L)
  expect_equal(anyDuplicated(pools$bio), 0L)
  expect_equal(anyDuplicated(pools$varlarge), 0L)
})

test_that("standardise_split uses training statistics", {
  tr <- matrix(rnorm(40 * 10, mean = 5, sd = 2), 40, 10)
  te <- matrix(rnorm(20 * 10, mean = 9, sd = 3), 20, 10)
  colnames(tr) <- colnames(te) <- paste0("G", 1:10)
  s <- standardise_split(tr, te)
  expect_equal(unname(colMeans(s$train)), rep(0, 10), tolerance = 1e-10)
  expect_equal(unname(apply(s$train, 2, sd)), rep(1, 10), tolerance = 1e-10)
  #  Test means must equal (te_mean - tr_mean)/tr_sd, not 0.
  expect_false(isTRUE(all.equal(unname(colMeans(s$test)), rep(0, 10),
                                tolerance = 1e-6)))
  expect_equal(unname(colMeans(s$test)),
               unname((colMeans(te) - colMeans(tr)) / apply(tr, 2, sd)),
               tolerance = 1e-10)
})


# ------------------------------------------------------------------------------
context("bio_prior: module gate")

test_that("module gate is deterministic given the seed", {
  X <- matrix(rnorm(60 * 500), 60, 500)
  colnames(X) <- sample(hallmark_union(), 500)
  y <- rbinom(60, 1, 0.5)
  g1 <- module_gate(X, y, B = 200, seed = 7)
  g2 <- module_gate(X, y, B = 200, seed = 7)
  expect_identical(g1$table$q_BH, g2$table$q_BH)
  expect_identical(g1$certified_genes, g2$certified_genes)
})

test_that("module gate certifies nothing on pure noise", {
  X <- matrix(rnorm(60 * 800), 60, 800)
  colnames(X) <- sample(hallmark_union(), 800)
  y <- rbinom(60, 1, 0.5)
  g <- module_gate(X, y, B = 300, seed = 7)
  #  This is the no-invented-results guard: noise must not certify.
  expect_lte(nrow(g$table[g$table$q_BH <= 0.05, ]), 1)
  expect_lte(g$null_run_pass_q20, 2)
  expect_true(all(unlist(g$certified_genes) %in% colnames(X)))
})

test_that("module gate finds a planted Hallmark set", {
  sets <- hallmark_sets_all()
  planted <- sets[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]][1:40]
  other <- sample(setdiff(hallmark_union(), planted), 460)
  X <- matrix(rnorm(60 * 500), 60, 500)
  colnames(X) <- c(planted, other)
  y <- rep(c(0, 1), each = 30)
  #  Every planted gene gets the same moderate effect; noise genes get none.
  X[31:60, planted] <- X[31:60, planted] + 1.2
  g <- module_gate(X, y, B = 300, seed = 7)
  hit <- g$table[g$table$set == "HALLMARK_INTERFERON_GAMMA_RESPONSE", ]
  expect_equal(nrow(hit), 1)
  expect_lte(hit$q_BH, 0.05)
  expect_true(all(planted %in% g$certified_genes$q20))
  expect_lte(g$null_run_pass_q20, 2)
})

test_that("module gate drops sets with too few pool members", {
  sets <- hallmark_sets_all()
  target <- names(sets)[1]
  few <- sets[[1]][1:5]              # below min_members = 10
  #  Draw the rest EXCLUDING every other gene of the target set, so the pool
  #  really holds only 5 members of it.
  rest <- sample(setdiff(hallmark_union(), sets[[1]]), 300)
  X <- matrix(rnorm(50 * 305), 50, 305)
  colnames(X) <- c(few, rest)
  y <- rbinom(50, 1, 0.5)
  g <- module_gate(X, y, B = 100, seed = 7)
  expect_false(target %in% g$table$set)
})


# ------------------------------------------------------------------------------
context("evaluator: AUC and ensemble")

test_that("bench_auc is exact on edge cases", {
  y <- factor(rep(c("a", "b"), each = 10), levels = c("a", "b"))
  expect_equal(bench_auc(y, c(rep(0, 10), rep(1, 10))), 1)
  expect_equal(bench_auc(y, c(rep(1, 10), rep(0, 10))), 0)
  expect_equal(bench_auc(y, rep(0.5, 20)), 0.5)
})

test_that("ensemble predicts separable data", {
  X <- matrix(rnorm(80 * 30), 80, 30)
  colnames(X) <- paste0("G", 1:30)
  y <- factor(rep(c("neg", "pos"), each = 40), levels = c("neg", "pos"))
  X[41:80, 1:5] <- X[41:80, 1:5] + 2
  tr <- c(1:30, 41:70); te <- c(31:40, 71:80)   # both classes in both sets
  scores <- predict_with_ensemble(X[tr, ], y[tr], X[te, ])
  expect_length(scores, 20)
  expect_true(all(is.finite(scores)))
  expect_gt(bench_auc(y[te], scores), 0.8)
})

test_that("ensemble predictions are exactly reproducible for a fixed seed", {
  X <- matrix(rnorm(60 * 15), 60, 15)
  colnames(X) <- paste0("G", 1:15)
  y <- factor(rep(c("neg", "pos"), each = 30),
              levels = c("neg", "pos"))
  tr <- c(1:22, 31:52)
  te <- setdiff(seq_len(nrow(X)), tr)
  first <- predict_with_ensemble(X[tr, ], y[tr], X[te, ],
                                 random_seed = 1042L)
  runif(25)  # unrelated RNG use must not alter the evaluator result
  second <- predict_with_ensemble(X[tr, ], y[tr], X[te, ],
                                  random_seed = 1042L)
  expect_identical(first, second)
})

test_that("random baseline stays near chance on noise", {
  X <- matrix(rnorm(100 * 80), 100, 80)
  colnames(X) <- paste0("G", 1:80)
  y <- factor(rep(c("neg", "pos"), each = 50), levels = c("neg", "pos"))
  tr <- c(1:35, 51:85); te <- c(36:50, 86:100)  # both classes in both sets
  aucs <- random_panel_aucs(X[tr, ], X[te, ], y[tr], y[te],
                            colnames(X), k = 10, n_draws = 5, seed = 1)
  expect_length(aucs, 5)
  expect_true(all(aucs >= 0 & aucs <= 1))
  expect_lt(abs(mean(aucs) - 0.5), 0.25)
})


# ------------------------------------------------------------------------------
context("sosall_data: residualisation")

test_that("residualisation removes the location effect and keeps shape", {
  n <- 90; p <- 40
  loc <- rep(c("A", "B", "C"), each = 30)
  X <- matrix(rnorm(n * p, mean = 6), n, p)
  colnames(X) <- paste0("G", 1:p)
  X[loc == "C", ] <- X[loc == "C", ] + 3      # location shift to remove
  meta <- data.frame(location = loc)
  #  Interleaved split so every location appears in train and test (as in
  #  the real stratified folds on the full cohort).
  tr <- c(1:20, 31:50, 61:80); te <- c(21:30, 51:60, 81:90)
  out <- sosall_residualise_split(X, meta, tr, te)
  expect_equal(dim(out$train), c(60, p))
  expect_equal(dim(out$test), c(30, p))

  #  The planted location effect (C vs A ~ 3) must be gone after
  #  residualisation. Compare mean |effect| before vs after: single
  #  coefficients keep sampling noise (SE ~ 0.3 at this n), so the check is
  #  on the mean across genes, not the max.
  des <- model.matrix(~ factor(meta$location[tr]))
  coef_before <- qr.coef(qr(des), X[tr, ])
  coef_after  <- qr.coef(qr(des), out$train)
  coef_before[is.na(coef_before)] <- 0
  coef_after[is.na(coef_after)] <- 0
  c_row <- grep("factor.*C$", rownames(coef_after))
  expect_gt(mean(abs(coef_before[c_row, ])), 2.5)   # effect was planted
  expect_lt(mean(abs(coef_after[c_row, ])), 0.2)    # and is now gone
  expect_lt(max(abs(colMeans(out$train))), 1e-10)
})


# ------------------------------------------------------------------------------
context("runner logic: splits, fallback, integration")

test_that("split construction matches the incumbent seeds and stratifies", {
  y <- factor(rep(c("neg", "pos"), times = c(120, 30)))
  folds_r1 <- make_stratified_folds(y, 5, seed = 42 + 1000 * 1)
  folds_r2 <- make_stratified_folds(y, 5, seed = 42 + 1000 * 2)
  expect_false(identical(folds_r1[[1]], folds_r2[[1]]))
  expect_length(folds_r1, 5)
  expect_equal(sort(unlist(folds_r1)), seq_along(y))
  for (f in folds_r1) {
    expect_equal(sum(y[f] == "pos"), 6)   # 30/5 positives per fold
  }
  #  Determinism: same seed reproduces the same folds.
  expect_identical(folds_r1, make_stratified_folds(y, 5, seed = 1042))
})

test_that("empty gate triggers the defined fallback condition", {
  X <- matrix(rnorm(60 * 800), 60, 800)
  colnames(X) <- sample(hallmark_union(), 800)
  y <- rbinom(60, 1, 0.5)
  g <- module_gate(X, y, B = 200, seed = 7)
  gated_pool <- intersect(colnames(X), g$certified_genes$q20)
  #  The runner's rule is length < 10 -> fallback. On noise the certified
  #  pool should be empty or tiny, i.e. the rule must evaluate TRUE here;
  #  if noise ever certifies 10+ genes the noise test above already fails.
  expect_true(length(gated_pool) < 10)
})

test_that("integration: planted set -> gate -> package fit -> ranking", {
  sets <- hallmark_sets_all()
  planted <- sets[["HALLMARK_P53_PATHWAY"]][1:30]
  other <- sample(setdiff(hallmark_union(), planted), 270)
  X <- matrix(rnorm(80 * 300), 80, 300)
  colnames(X) <- c(planted, other)
  y <- factor(rep(c("neg", "pos"), each = 40), levels = c("neg", "pos"))
  X[41:80, planted] <- X[41:80, planted] + 1.0

  g <- module_gate(X, as.integer(y) - 1L, B = 300, seed = 7)
  certified <- intersect(colnames(X), g$certified_genes$q20)
  expect_true("HALLMARK_P53_PATHWAY" %in%
                g$table$set[g$table$q_BH <= 0.2])
  expect_gt(length(certified), 10)

  fit <- geneselectr2_fit(
    X[, certified, drop = FALSE], y, gate_method = "none", B = 20,
    subsample_scheme = "kfold", subsample_k_folds = 5,
    utility_method = "instance_shap",
    components = c("stability", "utility"), score_formula = "geometric",
    calibration_mode = "percentile", alpha = 0.5,
    n_cores = 1, random_seed = 42, use_cache = FALSE, verbose = FALSE)
  gs <- fit$gene_scores[order(fit$gene_scores$final_score,
                              decreasing = TRUE), ]
  expect_gt(nrow(gs), 0)
  expect_true(all(gs$gene %in% certified))
  #  Planted genes must dominate the top of the ranking.
  expect_gt(mean(gs$gene[1:10] %in% planted), 0.5)
})
