#!/usr/bin/env Rscript

# Diagnose why the utility-pillar permutation null is LOPO-unstable on
# GSE107994 and GSE13355 but stable on the other five older-seven datasets.
# Uses ONLY saved artifacts: the rebuilt r1f1 null checkpoints
# (older7_calibration_null_r1f1.rds, 20 permutations x 20 subsamples, written
# 2026-09-03 by redesign/run_older7_calibration_diagnostics.R) plus the saved
# r1f1 fits, splits, and outcomes needed to re-derive the observed raw utility
# exactly. No model fitting happens here.
#
# Three hypotheses are discriminated:
#   H1 epsilon-floor interaction: instability concentrates in genes with
#      observed utility ~ 0, where the calibrated ratio is set by
#      eps / null_mean and is hypersensitive to dropping one permutation.
#   H2 small null: 20 permutations give noisy per-gene null estimates;
#      instability is spread across the utility distribution and shrinks when
#      more permutations are added.
#   H3 dataset property: the utility distribution itself is degenerate on
#      these datasets (most observed utility ~ 0 and/or the per-permutation
#      null values are intrinsically dispersed or zero-inflated), so neither
#      more permutations nor eps changes fix it.
#
# Key structural fact the analysis exploits: the LOPO diagnostic uses ONLY the
# null matrix (each permutation's utility vector is calibrated against the
# mean of the other 19). The observed outcome utility enters only through
# epsilon. So LOPO dispersion decomposes into (a) the intrinsic spread of a
# gene's 20 null values around their own mean -- which more permutations do
# NOT shrink -- and (b) noise in the leave-out reference mean, which shrinks
# like 1/(n_perm - 1). The variance decomposition below separates (a) from (b)
# directly, and the simulation quantifies what pure resampling noise (H2)
# would produce at n = 20 and n = 100 permutations.
#
# Usage: Rscript redesign/summarise_older7_null_instability.R [dataset]
#   dataset: one of GSE107994 GSE13355 GSE65682 GSE101794, or "all" (default).
#   Per-dataset CSVs are written per run; run "assemble" after all four to
#   write the combined files.

args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- if (length(args) >= 1L) args[[1L]] else "all"

unstable <- c("GSE107994", "GSE13355")
stable_refs <- c("GSE65682", "GSE101794")
all_datasets <- c(unstable, stable_refs)
stopifnot(dataset_arg %in% c(all_datasets, "all", "assemble"))

# Only the function definitions are needed; no model fitting, so glmnet is
# never called. utils.R: percentile01, compute_mi_vectorized, create_subsamples.
# calibration.R: calibrate_by_null. bio_prior.R: standardise_split.
source(file.path("package", "GeneSelectR", "R", "utils.R"))
source(file.path("package", "GeneSelectR", "R", "calibration.R"))
source(file.path("redesign", "R", "bio_prior.R"))

results_root <- file.path("redesign", "results_corrected")
out_dir <- file.path(results_root, "older7_null_instability_2026-09-03")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- function(name) file.path(out_dir, name)

# ---------------------------------------------------------------------------
# Helpers copied verbatim from redesign/run_older7_calibration_diagnostics.R
# so the observed raw utility is re-derived on the identical code path that
# produced the saved diagnostics. WHY copy instead of source: the driver is a
# script that executes on source; its functions are small and stable.
# ---------------------------------------------------------------------------

reconstruct_instance_frequency <- function(fit, genes) {
  shap_matrix <- fit$instance_importance
  if (is.null(shap_matrix) || is.null(colnames(shap_matrix)) ||
      anyDuplicated(colnames(shap_matrix)) ||
      !all(colnames(shap_matrix) %in% genes) ||
      any(!is.finite(shap_matrix)) || any(shap_matrix < 0)) {
    stop("Saved instance-level SHAP matrix is invalid.", call. = FALSE)
  }
  important <- matrix(
    FALSE, nrow = nrow(shap_matrix), ncol = ncol(shap_matrix),
    dimnames = dimnames(shap_matrix)
  )
  for (sample_idx in seq_len(nrow(shap_matrix))) {
    values <- shap_matrix[sample_idx, ]
    positive <- values[values > 0]
    if (length(positive) == 0L) next
    threshold <- as.numeric(stats::quantile(
      positive, 0.75, na.rm = TRUE
    ))
    important[sample_idx, ] <- values >= threshold
  }
  frequency <- stats::setNames(rep(0, length(genes)), genes)
  frequency[colnames(shap_matrix)] <- colMeans(important)
  frequency
}

reconstruct_raw_mi <- function(X, y, B = 50L, seed = 42L) {
  subsamples <- create_subsamples(
    y, B = B, random_seed = seed, scheme = "kfold", k_folds = 5L
  )
  mi_matrix <- vapply(subsamples, function(subsample) {
    compute_mi_vectorized(
      X[subsample$train, , drop = FALSE], y[subsample$train],
      method = "discrete", n_bins = 5L
    )
  }, numeric(ncol(X)))
  rowMeans(mi_matrix)
}

leave_one_permutation_out <- function(null_matrix, epsilon) {
  n_permutations <- nrow(null_matrix)
  t(vapply(seq_len(n_permutations), function(permutation_idx) {
    calibrate_by_null(
      null_matrix[permutation_idx, ],
      null_matrix[-permutation_idx, , drop = FALSE],
      epsilon = epsilon, winsorize_at = 4
    )
  }, numeric(ncol(null_matrix))))
}

ratio_metrics <- function(ratios) {
  c(
    median_ratio = stats::median(ratios),
    sd_log2_ratio = stats::sd(log2(ratios)),
    fraction_between_half_and_two = mean(ratios >= 0.5 & ratios <= 2),
    fraction_at_winsor_floor = mean(ratios <= 2^-4 + 1e-12),
    fraction_at_winsor_ceiling = mean(ratios >= 2^4 - 1e-12)
  )
}

# ---------------------------------------------------------------------------
# Per-dataset data loading + exact re-derivation of observed raw utility
# ---------------------------------------------------------------------------

load_dataset <- function(dataset) {
  dataset_dir <- file.path(results_root, "validation_benchmark", dataset)

  null_check <- readRDS(file.path(
    dataset_dir, "older7_calibration_null_r1f1.rds"
  ))
  null_utility <- null_check$null_utility
  utility_epsilon <- null_check$qa$utility_epsilon
  stopifnot(nrow(null_utility) == 20L)

  split <- readRDS(file.path(dataset_dir, "split_r1_f1.rds"))
  pool <- split$pools$var2000
  standardized <- standardise_split(
    split$train_raw[, pool, drop = FALSE],
    split$test_raw[, pool, drop = FALSE]
  )
  X <- standardized$train
  outcome <- droplevels(as.factor(
    readRDS(file.path(dataset_dir, "base_data.rds"))$outcome
  ))
  y <- droplevels(outcome[split$train_idx])

  metadata <- read.csv(file.path(
    dataset_dir, "ranking_r1_f1_GS_full_ungrouped_meta.csv"
  ), stringsAsFactors = FALSE)
  alpha <- as.numeric(metadata$value[metadata$key == "alpha"])
  alpha_tag <- if (identical(alpha, 0.5)) "0p5" else "1"
  fit <- readRDS(file.path(dataset_dir, sprintf(
    "fit_r1_f1_GS_full_ungrouped_a%s.rds", alpha_tag
  )))$fit

  scores <- fit$gene_scores[match(colnames(X), fit$gene_scores$gene), ,
                            drop = FALSE]
  stopifnot(identical(as.character(scores$gene), colnames(X)))

  shap_frequency <- unname(
    reconstruct_instance_frequency(fit, colnames(X))[colnames(X)]
  )
  raw_mi <- reconstruct_raw_mi(X, y, B = 50L, seed = 42L)
  observed_utility <- sqrt(shap_frequency * raw_mi + 1e-10)

  # Exactness gate: the re-derived observed utility, calibrated against the
  # saved null with the saved epsilon, must reproduce the saved u_scored
  # to 1e-12 (the same tolerance the driver's own QA used). If this fails the
  # re-derivation drifted from the pipeline and nothing downstream is valid.
  recomputed_mi_pct <- percentile01(raw_mi)
  recomputed_u <- calibrate_by_null(
    observed_utility, null_utility,
    epsilon = utility_epsilon, winsorize_at = 4
  )
  mi_diff <- max(abs(recomputed_mi_pct - scores$u_mi))
  u_diff <- max(abs(recomputed_u - scores$u_scored))
  if (mi_diff > 1e-12 || u_diff > 1e-12) {
    stop(sprintf(
      paste0("%s: observed-utility re-derivation does not reproduce the saved ",
             "fit (MI diff %.3g, utility diff %.3g)."),
      dataset, mi_diff, u_diff
    ), call. = FALSE)
  }

  list(
    dataset = dataset,
    genes = colnames(X),
    null_utility = null_utility,
    epsilon = utility_epsilon,
    observed_utility = observed_utility,
    saved_u_scored = scores$u_scored
  )
}

# ---------------------------------------------------------------------------
# Per-dataset analysis
# ---------------------------------------------------------------------------

analyse_dataset <- function(data) {
  dataset <- data$dataset
  null_utility <- data$null_utility
  eps <- data$epsilon
  obs <- data$observed_utility
  genes <- data$genes
  n_genes <- length(genes)
  n_perm <- nrow(null_utility)

  null_mean <- colMeans(null_utility)
  null_sd <- apply(null_utility, 2L, stats::sd)
  null_cv <- ifelse(null_mean > 0, null_sd / null_mean, NA_real_)
  null_n_positive <- colSums(null_utility > 0)

  # LOPO ratio matrix (20 x n_genes), the exact diagnostic from the driver.
  lopo <- leave_one_permutation_out(null_utility, eps)
  lopo_log2 <- log2(lopo)

  # Leave-one-permutation-out applied to the OBSERVED ratio: how much each
  # gene's actual evidence ratio u_scored moves when one permutation is
  # dropped from the reference mean. This is the decision-relevant version of
  # H1 -- LOPO-on-null is only a proxy for it.
  obs_ratio_lopo <- vapply(seq_len(n_perm), function(i) {
    ref_mean <- colMeans(null_utility[-i, , drop = FALSE])
    ratio <- (obs + eps) / (ref_mean + eps)
    2^pmax(pmin(log2(ratio), 4), -4)
  }, numeric(n_genes))

  # Observed-utility bins: 0 = at the raw-utility floor, 1..10 = deciles of
  # the above-floor values. The observed raw utility is
  # sqrt(shap_freq * raw_mi + 1e-10), so a gene with zero SHAP frequency or
  # zero MI does not get exactly 0, it gets sqrt(1e-10) = 1e-5. "At the
  # floor" means at that additive constant, not at literal zero -- on the
  # unstable datasets >= 95% of genes sit there, and rank-based deciles over
  # all genes collapse those ties into a single bin (verified on GSE107994:
  # 1974 of 2000 genes tied at 1e-5). H1 is precisely about that mass, so it
  # must be its own bin.
  utility_floor <- sqrt(1e-10)
  at_floor <- obs <= utility_floor * 1.0001
  obs_bin <- integer(n_genes)
  above_floor <- !at_floor
  if (any(above_floor)) {
    obs_bin[above_floor] <- ceiling(
      10 * rank(obs[above_floor], ties.method = "average") / sum(above_floor)
    )
  }

  # Epsilon contribution per gene, in log2 units: log2 ratio with the saved
  # eps minus log2 ratio with a near-zero eps (1e-12 guard against 0/0; both
  # winsorised the same way). Positive = eps inflated the ratio.
  ratio_epsfree <- (obs + 1e-12) / (null_mean + 1e-12)
  ratio_epsfree <- 2^pmax(pmin(log2(ratio_epsfree), 4), -4)
  eps_contribution_log2 <- log2(data$saved_u_scored) - log2(ratio_epsfree)

  gene_level <- data.frame(
    dataset = dataset,
    gene = genes,
    obs_utility = obs,
    obs_bin = obs_bin,
    null_mean = null_mean,
    null_sd = null_sd,
    null_cv = null_cv,
    null_n_positive = null_n_positive,
    null_max = apply(null_utility, 2L, max),
    ratio_observed = data$saved_u_scored,
    eps_contribution_log2 = eps_contribution_log2,
    lopo_min = apply(lopo, 2L, min),
    lopo_max = apply(lopo, 2L, max),
    lopo_log2_range = apply(lopo_log2, 2L, max) - apply(lopo_log2, 2L, min),
    obsratio_lopo_min = apply(obs_ratio_lopo, 2L, min),
    obsratio_lopo_max = apply(obs_ratio_lopo, 2L, max),
    obsratio_lopo_log2_range = apply(log2(obs_ratio_lopo), 2L, max) -
      apply(log2(obs_ratio_lopo), 2L, min),
    stringsAsFactors = FALSE
  )

  # Per-bin summary: LOPO range and null dispersion against the observed-
  # utility bins. If H1 holds, instability (lopo_log2_range, fraction of
  # genes with any LOPO value outside [0.5, 2]) concentrates in bin 0.
  decile_summary <- do.call(rbind, lapply(sort(unique(obs_bin)), function(b) {
    keep <- obs_bin == b
    lopo_keep <- lopo_log2[, keep, drop = FALSE]
    data.frame(
      dataset = dataset,
      obs_bin = b,
      n_genes = sum(keep),
      median_obs_utility = stats::median(obs[keep]),
      frac_obs_at_floor = mean(at_floor[keep]),
      median_null_mean = stats::median(null_mean[keep]),
      median_null_cv = stats::median(null_cv[keep], na.rm = TRUE),
      frac_null_all_zero = mean(null_mean[keep] == 0),
      median_null_n_positive = stats::median(null_n_positive[keep]),
      median_lopo_log2_range = stats::median(
        gene_level$lopo_log2_range[keep]
      ),
      frac_genes_any_lopo_outside_half_two = mean(vapply(
        which(keep),
        function(j) any(lopo[, j] < 0.5 | lopo[, j] > 2),
        logical(1)
      )),
      median_obsratio_lopo_log2_range = stats::median(
        gene_level$obsratio_lopo_log2_range[keep]
      ),
      median_ratio_observed = stats::median(data$saved_u_scored[keep]),
      median_eps_contribution_log2 = stats::median(
        eps_contribution_log2[keep]
      ),
      stringsAsFactors = FALSE
    )
  }))

  # Epsilon-floor summary: dataset-level scalars describing how degenerate
  # the observed and null utility distributions are (H3) and how much of the
  # ratio is set by eps rather than data (H1).
  total_utility <- sum(obs)
  sorted_obs <- sort(obs, decreasing = TRUE)
  eps_floor_summary <- data.frame(
    dataset = dataset,
    epsilon = eps,
    utility_floor = utility_floor,
    frac_obs_at_floor = mean(at_floor),
    n_obs_above_floor = sum(above_floor),
    top10_utility_share = sum(sorted_obs[seq_len(min(10L, n_genes))]) /
      total_utility,
    top50_utility_share = sum(sorted_obs[seq_len(min(50L, n_genes))]) /
      total_utility,
    frac_null_mean_zero = mean(null_mean == 0),
    median_null_mean = stats::median(null_mean),
    median_null_mean_over_eps = stats::median(null_mean / eps),
    median_null_cv = stats::median(null_cv, na.rm = TRUE),
    frac_null_cv_above_1 = mean(null_cv > 1, na.rm = TRUE),
    median_null_n_positive = stats::median(null_n_positive),
    # Numerator eps-dominated: obs contributes < 1% of (obs + eps).
    frac_genes_numerator_eps_dominated = mean(obs < 0.01 * eps),
    # Denominator eps-dominated: null_mean contributes < 1% of
    # (null_mean + eps); for these genes the ratio is pinned near 1 no
    # matter what the null does.
    frac_genes_denominator_eps_dominated = mean(null_mean < 0.01 * eps),
    # Among floor genes: where does their null mean sit relative to eps?
    # ratio ~ eps / (null_mean + eps) is only sensitive to dropping a
    # permutation when null_mean is comparable to or larger than eps.
    median_null_mean_over_eps_floor_genes = stats::median(
      (null_mean / eps)[at_floor]
    ),
    stringsAsFactors = FALSE
  )

  # Epsilon sensitivity: rerun LOPO and the observed-outcome calibration with
  # the saved eps, eps / 100, and a fixed 1e-6. If H1 (eps-floor interaction)
  # is the whole story, shrinking eps collapses the LOPO dispersion. If the
  # dispersion survives, the null values themselves are dispersed (H2/H3).
  eps_variants <- list(saved = eps, saved_over_100 = eps / 100, fixed_1e_6 = 1e-6)
  eps_sensitivity <- do.call(rbind, lapply(names(eps_variants), function(v) {
    e <- eps_variants[[v]]
    lopo_v <- leave_one_permutation_out(null_utility, e)
    obs_v <- calibrate_by_null(obs, null_utility, epsilon = e, winsorize_at = 4)
    m_lopo <- ratio_metrics(lopo_v)
    m_obs <- ratio_metrics(obs_v)
    data.frame(
      dataset = dataset,
      eps_variant = v,
      eps_value = e,
      lopo_median_ratio = m_lopo[["median_ratio"]],
      lopo_sd_log2 = m_lopo[["sd_log2_ratio"]],
      lopo_frac_within_half_two = m_lopo[["fraction_between_half_and_two"]],
      lopo_frac_at_floor = m_lopo[["fraction_at_winsor_floor"]],
      observed_median_ratio = m_obs[["median_ratio"]],
      observed_sd_log2 = m_obs[["sd_log2_ratio"]],
      observed_frac_within_half_two = m_obs[["fraction_between_half_and_two"]],
      stringsAsFactors = FALSE
    )
  }))

  # Variance decomposition of the LOPO log2 ratios (law of total variance
  # over the balanced 20 x n_genes layout):
  #   total = between-gene + within-gene
  # Between-gene: genes differ in their mean LOPO log-ratio -- driven by
  #   per-gene null shape (zero-inflated genes sit far below 1 in log space).
  # Within-gene: the 20 leave-out answers for ONE gene disagree. This is the
  #   part that measures null instability proper. It decomposes further into
  #   the intrinsic spread of the gene's null values around their own mean
  #   (NOT shrinkable by adding permutations) and the leave-out reference-
  #   mean noise (shrinks ~1/(n_perm - 1), the only part H2 can claim).
  between_var <- stats::var(colMeans(lopo_log2))
  within_var_per_gene <- apply(lopo_log2, 2L, stats::var)
  within_var <- mean(within_var_per_gene)
  total_var <- stats::var(as.vector(lopo_log2))
  xspread_var_per_gene <- vapply(seq_len(n_genes), function(j) {
    stats::var(log2(
      (null_utility[, j] + eps) / (null_mean[[j]] + eps)
    ))
  }, numeric(1))
  variance_decomposition <- data.frame(
    dataset = dataset,
    total_sd_log2 = sqrt(total_var),
    between_gene_sd_log2 = sqrt(between_var),
    within_gene_sd_log2 = sqrt(within_var),
    within_gene_xspread_sd_log2 = sqrt(mean(xspread_var_per_gene)),
    # Residual: leave-out reference-mean noise (approximate; the two
    # components are correlated, so this is within minus x-spread).
    within_gene_mean_noise_sd_log2 = sqrt(pmax(
      within_var - mean(xspread_var_per_gene), 0
    )),
    stringsAsFactors = FALSE
  )

  # H2 simulation: would pure resampling noise reproduce the observed LOPO
  # dispersion? Per gene, the null-generating process is modelled two ways:
  #   empirical_bootstrap: draw with replacement from the gene's own 20 null
  #     values (keeps the observed zero-inflation and tail shape exactly);
  #   gamma_matched: zero-inflated gamma -- empirical zero probability, gamma
  #     by moments on the positive values (the task's parametric version).
  # Each is simulated at n_perm = 20 (does it reproduce 2.2-2.3?) and at
  # n_perm = 100 (the H2 prediction for the confirmation run: LOPO with a
  # 99-permutation reference). Genes whose null is all zeros are degenerate
  # (ratio eps/eps = 1) and stay all zero.
  set.seed(20260903L)
  n_reps <- 25L
  simulate_lopo <- function(n_draw, method) {
    rep_metrics <- vapply(seq_len(n_reps), function(rep_idx) {
      sim_null <- matrix(0, nrow = n_draw, ncol = n_genes)
      for (j in which(null_mean > 0)) {
        values <- null_utility[, j]
        if (method == "empirical_bootstrap") {
          sim_null[, j] <- sample(values, n_draw, replace = TRUE)
        } else {
          positives <- values[values > 0]
          zero_prob <- 1 - length(positives) / length(values)
          is_zero <- stats::runif(n_draw) < zero_prob
          drawn <- numeric(n_draw)
          n_pos <- sum(!is_zero)
          if (n_pos > 0L) {
            if (length(positives) == 1L) {
              drawn[!is_zero] <- positives
            } else {
              m <- mean(positives)
              s <- stats::sd(positives)
              # Moment match; a gene with sd ~ 0 is near-degenerate, floor
              # the shape to avoid a divide-by-zero rate.
              shape <- max((m / max(s, 1e-12))^2, 1e-3)
              drawn[!is_zero] <- stats::rgamma(
                n_pos, shape = shape, rate = shape / m
              )
            }
          }
          sim_null[, j] <- drawn
        }
      }
      ratio_metrics(leave_one_permutation_out(sim_null, eps))
    }, numeric(5L))
    data.frame(
      dataset = dataset,
      method = method,
      n_permutations = n_draw,
      median_ratio_mean = mean(rep_metrics["median_ratio", ]),
      sd_log2_mean = mean(rep_metrics["sd_log2_ratio", ]),
      sd_log2_sd = stats::sd(rep_metrics["sd_log2_ratio", ]),
      frac_within_half_two_mean = mean(
        rep_metrics["fraction_between_half_and_two", ]
      ),
      stringsAsFactors = FALSE
    )
  }
  observed_row <- data.frame(
    dataset = dataset,
    method = "observed",
    n_permutations = n_perm,
    median_ratio_mean = unname(ratio_metrics(lopo)[["median_ratio"]]),
    sd_log2_mean = unname(ratio_metrics(lopo)[["sd_log2_ratio"]]),
    sd_log2_sd = NA_real_,
    frac_within_half_two_mean = unname(
      ratio_metrics(lopo)[["fraction_between_half_and_two"]]
    ),
    stringsAsFactors = FALSE
  )
  simulation <- rbind(
    observed_row,
    simulate_lopo(20L, "empirical_bootstrap"),
    simulate_lopo(100L, "empirical_bootstrap"),
    simulate_lopo(20L, "gamma_matched"),
    simulate_lopo(100L, "gamma_matched")
  )

  list(
    gene_level = gene_level,
    decile_summary = decile_summary,
    eps_floor_summary = eps_floor_summary,
    eps_sensitivity = eps_sensitivity,
    variance_decomposition = variance_decomposition,
    simulation = simulation
  )
}

write_dataset_outputs <- function(result, dataset) {
  for (name in names(result)) {
    write.csv(
      result[[name]],
      out_path(sprintf("%s_%s.csv", name, dataset)),
      row.names = FALSE
    )
  }
  cat(sprintf(
    paste0("[%s] LOPO sd_log2 %.3f | within-gene sd %.3f (x-spread %.3f) | ",
           "obs at floor %.1f%% | null CV>1 %.1f%%\n"),
    dataset,
    result$variance_decomposition$total_sd_log2,
    result$variance_decomposition$within_gene_sd_log2,
    result$variance_decomposition$within_gene_xspread_sd_log2,
    100 * result$eps_floor_summary$frac_obs_at_floor,
    100 * result$eps_floor_summary$frac_null_cv_above_1
  ))
}

if (dataset_arg %in% c(all_datasets, "all")) {
  targets <- if (dataset_arg == "all") all_datasets else dataset_arg
  for (dataset in targets) {
    cat(sprintf("[start] %s\n", dataset))
    data <- load_dataset(dataset)
    result <- analyse_dataset(data)
    write_dataset_outputs(result, dataset)
  }
}

if (dataset_arg %in% c("all", "assemble")) {
  # Combine the per-dataset CSVs into the tidy deliverables. Reads what the
  # per-dataset runs wrote, so "assemble" alone refreshes the combined files.
  names_to_combine <- c(
    "gene_level", "decile_summary", "eps_floor_summary",
    "eps_sensitivity", "variance_decomposition", "simulation"
  )
  for (name in names_to_combine) {
    pieces <- lapply(all_datasets, function(dataset) {
      path <- out_path(sprintf("%s_%s.csv", name, dataset))
      if (!file.exists(path)) {
        stop(sprintf("Missing per-dataset output: %s", path), call. = FALSE)
      }
      read.csv(path, stringsAsFactors = FALSE)
    })
    write.csv(
      do.call(rbind, pieces),
      out_path(sprintf("combined_%s.csv", name)),
      row.names = FALSE
    )
  }
  cat(sprintf("[assemble] combined CSVs written to %s\n", out_dir))
}
