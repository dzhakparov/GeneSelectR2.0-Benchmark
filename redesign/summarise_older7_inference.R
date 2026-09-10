# summarise_older7_inference.R
#
# Uncertainty statements for the benchmark's key paired comparisons, replacing
# the paired Wilcoxon over the 15 overlapping CV splits (invalid: folds within
# a repeat share observations, so the 15 differences are not independent).
#
# Structure per dataset: 15 outer splits = 3 repeats x 5 folds. Folds within a
# repeat share observations; repeats are re-partitions. The unit that is
# approximately independent is the repeat, not the split.
#
# Three intervals are reported for every contrast:
#
#   1. Repeat-cluster bootstrap percentile interval (primary). The 3 repeats
#      are resampled with replacement within each dataset, keeping all folds
#      of a resampled repeat. B = 10000, fixed seed. WARNING: with only 3
#      clusters the bootstrap distribution is coarse (the resampled mean can
#      take few distinct values), so these intervals are APPROXIMATE and tend
#      to be optimistic. Treat them as a lower bound on the true width.
#   2. t interval over the 3 repeat-level means (df = 2). This is the
#      conservative companion: 3 independent observations, heavy tails.
#      For dataset-balanced estimates the t interval is taken over the 7
#      per-dataset means (df = 6).
#   3. Nadeau-Bengio corrected resampled t interval over all 15 splits
#      (SE = sqrt((1/n + n_test/n_train) * var), 5-fold so ratio 1/4, df = 14).
#      Mirrored from redesign/summarise_frozen_external_competitors.R
#      (lines ~112-146) for consistency with the frozen external runs. It
#      corrects for overlap via the train/test ratio but does not model the
#      repeat clusters; reported as a third, comparability interval.
#
# Inputs:
#   - redesign/results_corrected/integrated_interpretation_2026-09-03/
#       dev7_fixed_panel_competitor_by_split.csv
#       (per dataset x repeat x fold x method x k AUC and delta vs 30-draw
#        Random; competitor comparisons use delta_random_30draw because the
#        primary metric is AUC-minus-Random at matched k)
#   - redesign/results_corrected/older7_remaining_analyses_2026-09-02/
#       combined_ablation_vs_random30_by_split.csv
#       (per-split per-variant delta_random_30draw, including GS_current;
#        ablation differences are computed per split from this file)
#   - combined_ablation_vs_current.csv and combined_component_contributions.csv
#       are used only to reconcile means.
#
# Outputs:
#   - redesign/results_corrected/integrated_interpretation_2026-09-03/
#       dev7_paired_inference_2026-09-03.csv
#   - redesign/results_corrected/integrated_interpretation_2026-09-03/
#       inference_notes.md

# The script lives in redesign/; resolve the repo root from its own path so it
# can be run from anywhere (Rscript redesign/summarise_older7_inference.R).
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- dirname(script_dir)

competitor_file <- file.path(
  repo_root, "redesign/results_corrected/integrated_interpretation_2026-09-03",
  "dev7_fixed_panel_competitor_by_split.csv"
)
ablation_split_file <- file.path(
  repo_root, "redesign/results_corrected/older7_remaining_analyses_2026-09-02",
  "combined_ablation_vs_random30_by_split.csv"
)
ablation_summary_file <- file.path(
  repo_root, "redesign/results_corrected/older7_remaining_analyses_2026-09-02",
  "combined_ablation_vs_current.csv"
)
out_csv <- file.path(
  repo_root, "redesign/results_corrected/integrated_interpretation_2026-09-03",
  "dev7_paired_inference_2026-09-03.csv"
)
out_notes <- file.path(
  repo_root, "redesign/results_corrected/integrated_interpretation_2026-09-03",
  "inference_notes.md"
)

K_PRIMARY <- c(10, 20, 50)
COMPETITORS <- c("RF_importance", "DGE", "Boruta", "mRMR", "LASSO", "ElasticNet")
GS_ARM <- "GS_full_ungrouped"
ABLATION_VARIANTS <- c("GS_recurrence_only", "GS_SHAP_only", "GS_MI_only",
                       "GS_SHAPxMI")
B_BOOT <- 10000L
BOOT_SEED <- 20260903L
ALPHA <- 0.05

# ---------------------------------------------------------------------------
# Interval helpers
# ---------------------------------------------------------------------------

# Repeat-cluster bootstrap percentile interval for one dataset's split-level
# paired differences. diffs_by_repeat is a list of length 3, each element the
# vector of per-split differences for that repeat (length = n_folds, or
# n_folds when the endpoint already averages over k inside the split).
boot_repeat_cluster <- function(diffs_by_repeat, B = B_BOOT) {
  n_repeats <- length(diffs_by_repeat)
  reps <- seq_len(n_repeats)
  boot_means <- replicate(B, {
    picked <- sample(reps, n_repeats, replace = TRUE)
    mean(unlist(diffs_by_repeat[picked], use.names = FALSE))
  })
  stats::quantile(boot_means, c(ALPHA / 2, 1 - ALPHA / 2), names = FALSE)
}

# t interval over repeat-level means, df = n_repeats - 1 (= 2 here).
t_over_repeat_means <- function(diffs_by_repeat) {
  repeat_means <- vapply(diffs_by_repeat, mean, numeric(1L))
  m <- mean(repeat_means)
  if (length(repeat_means) < 2L) return(c(m, m))
  se <- stats::sd(repeat_means) / sqrt(length(repeat_means))
  crit <- stats::qt(1 - ALPHA / 2, df = length(repeat_means) - 1L)
  c(m - crit * se, m + crit * se)
}

# Nadeau-Bengio corrected resampled t interval over all splits, mirroring
# redesign/summarise_frozen_external_competitors.R: 5-fold CV, so the
# test/train ratio is 1/4.
nb_corrected_t <- function(diffs) {
  n <- length(diffs)
  m <- mean(diffs)
  se <- sqrt((1 / n + 1 / 4) * stats::var(diffs))
  crit <- stats::qt(1 - ALPHA / 2, df = n - 1L)
  c(m - crit * se, m + crit * se)
}

# All three intervals for a vector of per-split differences tagged by repeat.
# diffs: numeric; repeat_id: same length, values in 1..3.
all_intervals_one_dataset <- function(diffs, repeat_id) {
  by_repeat <- split(diffs, repeat_id)
  by_repeat <- by_repeat[order(as.integer(names(by_repeat)))]
  boot <- boot_repeat_cluster(by_repeat)
  tr <- t_over_repeat_means(by_repeat)
  nb <- nb_corrected_t(diffs)
  list(mean_diff = mean(diffs), split_sd = stats::sd(diffs),
       boot_lo = boot[1L], boot_hi = boot[2L],
       t_lo = tr[1L], t_hi = tr[2L],
       nb_lo = nb[1L], nb_hi = nb[2L])
}

# Dataset-balanced estimate: average the per-dataset means first, then
# intervals. Bootstrap: resample repeats within each dataset, average the
# per-dataset resample means. t interval: over the 7 per-dataset means (df=6).
# The Nadeau-Bengio interval is not defined across datasets (its correction
# assumes one dataset's CV splits), so it is NA at balanced scope.
all_intervals_balanced <- function(per_dataset) {
  # per_dataset: named list; each element a list with $by_repeat (list of
  # split-diff vectors) and $mean_diff.
  ds_names <- names(per_dataset)
  dataset_means <- vapply(per_dataset, function(x) x$mean_diff, numeric(1L))
  balanced_mean <- mean(dataset_means)

  n_repeats <- length(per_dataset[[1L]]$by_repeat)
  boot_means <- replicate(B_BOOT, {
    per_ds <- vapply(ds_names, function(ds) {
      by_repeat <- per_dataset[[ds]]$by_repeat
      picked <- sample(seq_len(n_repeats), n_repeats, replace = TRUE)
      mean(unlist(by_repeat[picked], use.names = FALSE))
    }, numeric(1L))
    mean(per_ds)
  })
  boot <- stats::quantile(boot_means, c(ALPHA / 2, 1 - ALPHA / 2),
                          names = FALSE)

  se <- stats::sd(dataset_means) / sqrt(length(dataset_means))
  crit <- stats::qt(1 - ALPHA / 2, df = length(dataset_means) - 1L)
  tr <- c(balanced_mean - crit * se, balanced_mean + crit * se)

  list(mean_diff = balanced_mean,
       split_sd = NA_real_,  # split-level sd has no meaning across datasets
       boot_lo = boot[1L], boot_hi = boot[2L],
       t_lo = tr[1L], t_hi = tr[2L],
       nb_lo = NA_real_, nb_hi = NA_real_)
}

pack_row <- function(comparison, scope, k_label, est, n_splits, n_repeats,
                     note) {
  data.frame(
    comparison = comparison,
    scope = scope,
    k = k_label,
    mean_diff = est$mean_diff,
    split_sd = est$split_sd,
    boot_ci_lo = est$boot_lo,
    boot_ci_hi = est$boot_hi,
    t_repeat_ci_lo = est$t_lo,
    t_repeat_ci_hi = est$t_hi,
    nadeau_bengio_ci_lo = est$nb_lo,
    nadeau_bengio_ci_hi = est$nb_hi,
    n_splits = n_splits,
    n_repeats = n_repeats,
    note = note,
    stringsAsFactors = FALSE
  )
}

BOOT_NOTE <- paste(
  "boot_* = repeat-cluster bootstrap percentile interval (resample the 3",
  "repeats within each dataset, keep all folds, B=10000, seed 20260903);",
  "APPROXIMATE with only 3 clusters and likely optimistic.",
  "t_repeat_* = t interval over repeat-level means (df=2 per dataset; df=6",
  "over per-dataset means at balanced scope).",
  "nadeau_bengio_* = corrected resampled t over all 15 splits (df=14),",
  "mirrored from redesign/summarise_frozen_external_competitors.R for",
  "consistency with the frozen external runs; NA at balanced scope."
)

set.seed(BOOT_SEED)

rows <- list()

# ---------------------------------------------------------------------------
# Competitor comparisons: GS_full_ungrouped minus each competitor, on
# delta_random_30draw (AUC minus 30-draw Random, the primary metric).
# ---------------------------------------------------------------------------

comp <- read.csv(competitor_file, stringsAsFactors = FALSE)
comp <- comp[comp$k %in% K_PRIMARY, ]
comp <- comp[comp$arm %in% c(GS_ARM, COMPETITORS), ]

gs_splits <- comp[comp$arm == GS_ARM,
                  c("dataset", "repeat_idx", "fold_idx", "k",
                    "delta_random_30draw")]
names(gs_splits)[5L] <- "gs_delta"
paired_comp <- merge(
  comp[comp$arm != GS_ARM,
       c("dataset", "repeat_idx", "fold_idx", "k", "arm",
         "delta_random_30draw")],
  gs_splits,
  by = c("dataset", "repeat_idx", "fold_idx", "k")
)
paired_comp$diff <- paired_comp$gs_delta - paired_comp$delta_random_30draw
stopifnot(nrow(paired_comp) ==
            7L * 3L * 5L * length(K_PRIMARY) * length(COMPETITORS))

for (method in COMPETITORS) {
  sub <- paired_comp[paired_comp$arm == method, ]

  # (a) primary endpoint: mean of the per-k differences over k in {10,20,50}
  # inside each split.
  primary <- aggregate(diff ~ dataset + repeat_idx + fold_idx, sub, mean)

  # (b) each k separately.
  endpoints <- c(list(primary = primary),
                 lapply(K_PRIMARY, function(kk) sub[sub$k == kk, ]))
  endpoint_labels <- c("primary", as.character(K_PRIMARY))

  per_dataset_est <- list()
  for (ds in sort(unique(sub$dataset))) {
    for (ei in seq_along(endpoints)) {
      dd <- endpoints[[ei]][endpoints[[ei]]$dataset == ds, ]
      dd <- dd[order(dd$repeat_idx, dd$fold_idx), ]
      est <- all_intervals_one_dataset(dd$diff, dd$repeat_idx)
      rows[[length(rows) + 1L]] <- pack_row(
        paste0("GS_full_ungrouped_minus_", method), ds, endpoint_labels[ei],
        est, nrow(dd), length(unique(dd$repeat_idx)), BOOT_NOTE
      )
      if (endpoint_labels[ei] == "primary") {
        per_dataset_est[[ds]] <- list(
          by_repeat = split(dd$diff, dd$repeat_idx)[order(unique(dd$repeat_idx))],
          mean_diff = est$mean_diff
        )
      }
    }
  }

  # Dataset-balanced scope, primary endpoint and per k.
  for (ei in seq_along(endpoints)) {
    per_ds_ei <- lapply(sort(unique(sub$dataset)), function(ds) {
      dd <- endpoints[[ei]][endpoints[[ei]]$dataset == ds, ]
      dd <- dd[order(dd$repeat_idx, dd$fold_idx), ]
      list(by_repeat = split(dd$diff, dd$repeat_idx)[order(unique(dd$repeat_idx))],
           mean_diff = mean(dd$diff))
    })
    names(per_ds_ei) <- sort(unique(sub$dataset))
    est <- all_intervals_balanced(per_ds_ei)
    rows[[length(rows) + 1L]] <- pack_row(
      paste0("GS_full_ungrouped_minus_", method), "balanced",
      endpoint_labels[ei], est,
      nrow(endpoints[[ei]]), 3L, BOOT_NOTE
    )
  }
}

# ---------------------------------------------------------------------------
# Ablation comparisons: variant minus GS_current, paired per split, per k.
# Differences are computed from the per-split file so the repeat structure is
# available for the cluster bootstrap.
# ---------------------------------------------------------------------------

abl <- read.csv(ablation_split_file, stringsAsFactors = FALSE)
abl <- abl[abl$k %in% K_PRIMARY, ]
abl_cur <- abl[abl$variant == "GS_current",
               c("dataset", "repeat_idx", "fold_idx", "k",
                 "delta_random_30draw")]
names(abl_cur)[5L] <- "current_delta"
paired_abl <- merge(
  abl[abl$variant %in% ABLATION_VARIANTS,
      c("dataset", "repeat_idx", "fold_idx", "k", "variant",
        "delta_random_30draw")],
  abl_cur,
  by = c("dataset", "repeat_idx", "fold_idx", "k")
)
paired_abl$diff <- paired_abl$delta_random_30draw - paired_abl$current_delta

for (variant in ABLATION_VARIANTS) {
  sub <- paired_abl[paired_abl$variant == variant, ]
  for (kk in K_PRIMARY) {
    dk <- sub[sub$k == kk, ]
    per_ds_ei <- list()
    for (ds in sort(unique(dk$dataset))) {
      dd <- dk[dk$dataset == ds, ]
      dd <- dd[order(dd$repeat_idx, dd$fold_idx), ]
      est <- all_intervals_one_dataset(dd$diff, dd$repeat_idx)
      rows[[length(rows) + 1L]] <- pack_row(
        paste0(variant, "_minus_GS_current"), ds, as.character(kk),
        est, nrow(dd), length(unique(dd$repeat_idx)), BOOT_NOTE
      )
      per_ds_ei[[ds]] <- list(
        by_repeat = split(dd$diff, dd$repeat_idx)[order(unique(dd$repeat_idx))],
        mean_diff = est$mean_diff
      )
    }
    est <- all_intervals_balanced(per_ds_ei)
    rows[[length(rows) + 1L]] <- pack_row(
      paste0(variant, "_minus_GS_current"), "balanced", as.character(kk),
      est, nrow(dk), 3L, BOOT_NOTE
    )
  }
}

result <- do.call(rbind, rows)
rownames(result) <- NULL

# ---------------------------------------------------------------------------
# Reconciliation checks: reported means must equal plain averages of inputs.
# ---------------------------------------------------------------------------

# 1. Competitor primary means vs plain average of per-split differences.
check_comp <- aggregate(diff ~ dataset + arm, paired_comp, mean)
our_comp <- result[result$k == "primary" & result$scope != "balanced",
                   c("scope", "comparison", "mean_diff")]
m1 <- merge(check_comp, our_comp,
            by.x = c("dataset", "arm"),
            by.y = c("scope", "comparison"))
m1$match <- abs(m1$diff - m1$mean_diff) < 1e-12
# comparison column carries the "GS_full_ungrouped_minus_" prefix; strip it.
our_comp$arm_stripped <- sub("^GS_full_ungrouped_minus_", "",
                             our_comp$comparison)
m1 <- merge(check_comp, our_comp, by.x = c("dataset", "arm"),
            by.y = c("scope", "arm_stripped"))
stopifnot(all(abs(m1$diff - m1$mean_diff) < 1e-12))

# 2. Ablation per-dataset per-k means vs combined_ablation_vs_current.csv.
ref_abl <- read.csv(ablation_summary_file, stringsAsFactors = FALSE)
ref_abl <- ref_abl[ref_abl$k %in% K_PRIMARY, ]
our_abl <- result[grepl("_minus_GS_current$", result$comparison) &
                    result$scope != "balanced",
                  c("scope", "comparison", "k", "mean_diff")]
our_abl$variant <- sub("_minus_GS_current$", "", our_abl$comparison)
our_abl$k <- as.integer(our_abl$k)
m2 <- merge(ref_abl, our_abl,
            by.x = c("dataset", "k", "variant"),
            by.y = c("scope", "k", "variant"))
stopifnot(nrow(m2) == nrow(ref_abl))
stopifnot(all(abs(m2$mean_difference_vs_current - m2$mean_diff) < 1e-9))
cat(sprintf("Reconciliation OK: %d competitor means and %d ablation means",
            nrow(m1), nrow(m2)), "match plain averages of the inputs.\n")

write.csv(result, out_csv, row.names = FALSE, quote = TRUE)
cat("Wrote", out_csv, "with", nrow(result), "rows.\n")

# ---------------------------------------------------------------------------
# inference_notes.md: which intervals exclude zero.
# ---------------------------------------------------------------------------

excl <- function(lo, hi) ifelse(is.na(lo), NA, lo > 0 | hi < 0)

fmt <- function(x) sprintf("%+.4f", x)

notes <- c(
  "# Paired inference notes (2026-09-03)",
  "",
  "Replaces the paired Wilcoxon over the 15 overlapping CV splits. Structure:",
  "15 splits per dataset = 3 repeats x 5 folds; folds within a repeat share",
  "observations, so the repeat is the independent unit.",
  "",
  "Intervals reported per contrast (all in `dev7_paired_inference_2026-09-03.csv`):",
  "",
  "- `boot_*`: repeat-cluster bootstrap percentile interval (resample the 3",
  "  repeats within each dataset, keep all folds, B = 10000, seed 20260903).",
  "  **With only 3 clusters these intervals are approximate and likely",
  "  optimistic**; do not read them as exact 95% coverage.",
  "- `t_repeat_*`: t interval over the 3 repeat-level means (df = 2; at",
  "  balanced scope over the 7 per-dataset means, df = 6). Conservative",
  "  companion to the bootstrap.",
  "- `nadeau_bengio_*`: corrected resampled t over all 15 splits (df = 14),",
  "  same method as the frozen external runs",
  "  (`redesign/summarise_frozen_external_competitors.R`); per dataset only.",
  "",
  "A difference is counted as excluding zero below only when BOTH the",
  "bootstrap and the repeat-level t interval exclude zero (conservative",
  "rule). Primary endpoint = mean delta vs 30-draw Random over k in",
  "{10, 20, 50} per split. Balanced = per-dataset means averaged with equal",
  "weight over the 7 datasets.",
  ""
)

# GS vs competitors, primary endpoint.
notes <- c(notes, "## GS_full_ungrouped vs competitors, primary endpoint", "")
prim <- result[result$k == "primary" &
                 grepl("^GS_full_ungrouped_minus_", result$comparison), ]
prim$both_excl <- excl(prim$boot_ci_lo, prim$boot_ci_hi) &
  excl(prim$t_repeat_ci_lo, prim$t_repeat_ci_hi)

bal <- prim[prim$scope == "balanced", ]
for (i in seq_len(nrow(bal))) {
  r <- bal[i, ]
  cmp <- sub("^GS_full_ungrouped_minus_", "", r$comparison)
  verdict <- if (r$both_excl) "excludes zero" else "includes zero"
  notes <- c(notes, sprintf(
    "- balanced, GS - %s: mean %s, boot [%s, %s], t-repeat [%s, %s] -> %s",
    cmp, fmt(r$mean_diff), fmt(r$boot_ci_lo), fmt(r$boot_ci_hi),
    fmt(r$t_repeat_ci_lo), fmt(r$t_repeat_ci_hi), verdict))
}
n_ds_excl <- sum(prim$both_excl[prim$scope != "balanced"], na.rm = TRUE)
n_ds_tot <- sum(prim$scope != "balanced")
notes <- c(notes, "",
           sprintf("Per dataset: %d of %d dataset-level competitor contrasts",
                   n_ds_excl, n_ds_tot),
           "have both intervals excluding zero.", "")

# Ablation, balanced scope per k.
notes <- c(notes, "## Ablation: variant - GS_current, dataset-balanced", "")
abl_rows <- result[grepl("_minus_GS_current$", result$comparison) &
                     result$scope == "balanced", ]
abl_rows$both_excl <- excl(abl_rows$boot_ci_lo, abl_rows$boot_ci_hi) &
  excl(abl_rows$t_repeat_ci_lo, abl_rows$t_repeat_ci_hi)
for (i in seq_len(nrow(abl_rows))) {
  r <- abl_rows[i, ]
  verdict <- if (r$both_excl) "excludes zero" else "includes zero"
  notes <- c(notes, sprintf(
    "- k=%s, %s: mean %s, boot [%s, %s], t-repeat [%s, %s] -> %s",
    r$k, sub("_minus_GS_current$", "", r$comparison), fmt(r$mean_diff),
    fmt(r$boot_ci_lo), fmt(r$boot_ci_hi), fmt(r$t_repeat_ci_lo),
    fmt(r$t_repeat_ci_hi), verdict))
}
notes <- c(notes, "",
           "Per-dataset ablation results are in the CSV; intervals wider than",
           "the competitor contrasts because the differences are smaller.", "")

writeLines(notes, out_notes)
cat("Wrote", out_notes, "\n")
