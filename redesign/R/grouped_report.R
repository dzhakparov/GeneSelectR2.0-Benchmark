# ==============================================================================
#  Report builder for the grouped benchmark. Reads the eval CSVs of both
#  cohorts plus the saved incumbent results and writes:
#    - <out_root>/REPORT.md            human-readable comparison
#    - <out_root>/summary_by_arm.csv   mean AUC / delta per arm x k x cohort
#    - <out_root>/gate_summary.csv     module-gate statistics per split
#
#  Incumbent anchors (read, never re-run): GS_stab_util, GS_uncalibrated,
#  GS_utility_only, GS_harmonic, RF_importance, Boruta, Random from
#  results_imvigor210/2026-08-08_kfold and results_sosall/2026-08-13_kfold.
#  Anchor caveat stated in the report: incumbent GS variants ran the
#  evidence-ratio calibration null; the new arms ran calibration_mode =
#  "percentile". GS_uncalibrated is the closest anchor to the new recipe.
# ==============================================================================

build_grouped_report <- function(out_root) {
  inc_paths <- c(
    imvigor210 = file.path("results_imvigor210", "2026-08-08_kfold", "data"),
    sosall     = file.path("results_sosall", "2026-08-13_kfold", "data"))

  anchor_methods <- c("GS_stab_util", "GS_uncalibrated", "GS_utility_only",
                      "GS_harmonic", "RF_importance", "Boruta")

  summaries <- list()
  gate_rows <- list()
  report_lines <- c(
    "# Grouped benchmark report",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    "New arms: ranking = stability x utility (geometric, instance_shap,",
    "B=50 kfold, alpha=0.5, calibration_mode=percentile). Biology enters",
    "only as a prior (Hallmark pool and/or pre-selection module gate).",
    "Delta = AUC minus the matched-pool random baseline (3 draws/split).",
    "")

  for (cohort in c("imvigor210", "sosall")) {
    eval_path <- file.path(out_root, cohort, "eval_results.csv")
    if (!file.exists(eval_path)) next
    ev <- read.csv(eval_path)

    #  Mean AUC per arm x k, and delta vs the matched-pool Random arm.
    pool_of_arm <- function(arm) sub("_ungrouped$|_grouped$", "", arm)
    arms <- unique(ev$arm[!grepl("^Random_", ev$arm)])
    for (arm in arms) {
      pool <- pool_of_arm(arm)
      rnd <- ev[ev$arm == paste0("Random_", pool), ]
      for (k in sort(unique(ev$k))) {
        a <- ev[ev$arm == arm & ev$k == k, ]
        r <- rnd[rnd$k == k, ]
        if (nrow(a) == 0) next
        summaries[[length(summaries) + 1]] <- data.frame(
          cohort = cohort, arm = arm, k = k,
          n_splits = nrow(a),
          mean_AUC = mean(a$AUC, na.rm = TRUE),
          sd_AUC = sd(a$AUC, na.rm = TRUE),
          mean_Random = ifelse(nrow(r) > 0, mean(r$AUC, na.rm = TRUE),
                               NA_real_),
          delta = mean(a$AUC, na.rm = TRUE) -
            ifelse(nrow(r) > 0, mean(r$AUC, na.rm = TRUE), NA_real_),
          n_fallback = sum(a$gate_fallback))
      }
    }

    #  Gate statistics from the split RDS files.
    split_files <- list.files(file.path(out_root, cohort),
                              pattern = "^split_r.*\\.rds$", full.names = TRUE)
    for (sf in split_files) {
      sp <- readRDS(sf)
      m <- regmatches(basename(sf),
                      regexec("split_r(\\d+)_f(\\d+)\\.rds", basename(sf)))[[1]]
      for (pn in c("var2000", "bio")) {
        g <- sp$gates[[pn]]
        gate_rows[[length(gate_rows) + 1]] <- data.frame(
          cohort = cohort, repeat_idx = as.integer(m[2]),
          fold_idx = as.integer(m[3]), pool = pn,
          pool_size = length(sp$pools[[pn]]),
          sets_tested = nrow(g$table),
          sets_q05 = sum(g$table$q_BH <= 0.05),
          sets_q20 = sum(g$table$q_BH <= 0.20),
          genes_q20 = length(g$certified_genes$q20),
          null_selftest_pass = g$null_run_pass_q20)
      }
    }

    #  Incumbent anchors for this cohort.
    inc_delta <- read.csv(file.path(inc_paths[[cohort]],
                                    "delta_vs_random.csv"))
    inc_delta <- inc_delta[inc_delta$Method %in% anchor_methods, ]

    report_lines <- c(report_lines,
      sprintf("## %s", cohort), "",
      "New arms (delta vs matched-pool Random):", "",
      "| arm | k | mean AUC | Random | delta | splits | gate fallbacks |",
      "|---|---|---|---|---|---|---|")
    s <- do.call(rbind, summaries)
    s <- s[s$cohort == cohort, ]
    for (i in seq_len(nrow(s))) {
      report_lines <- c(report_lines, sprintf(
        "| %s | %d | %.3f | %.3f | %+.3f | %d | %d |",
        s$arm[i], s$k[i], s$mean_AUC[i], s$mean_Random[i], s$delta[i],
        s$n_splits[i], s$n_fallback[i]))
    }
    report_lines <- c(report_lines, "",
      "Incumbent anchors (delta vs Random, saved runs):", "",
      "| method | k | delta |", "|---|---|---|")
    for (i in seq_len(nrow(inc_delta))) {
      report_lines <- c(report_lines, sprintf(
        "| %s | %d | %+.3f |", inc_delta$Method[i], inc_delta$k[i],
        inc_delta$AUC_delta[i]))
    }
    report_lines <- c(report_lines, "")
  }

  summary_df <- if (length(summaries) > 0) do.call(rbind, summaries) else
    data.frame()
  gate_df <- if (length(gate_rows) > 0) do.call(rbind, gate_rows) else
    data.frame()
  write.csv(summary_df, file.path(out_root, "summary_by_arm.csv"),
            row.names = FALSE)
  write.csv(gate_df, file.path(out_root, "gate_summary.csv"),
            row.names = FALSE)

  if (nrow(gate_df) > 0) {
    report_lines <- c(report_lines, "## Module gate statistics", "",
      "| cohort | pool | splits | mean sets tested | mean sets q<=0.05 | mean sets q<=0.20 | mean genes certified | splits with empty gate | null self-test violations |",
      "|---|---|---|---|---|---|---|---|---|")
    for (coh in unique(gate_df$cohort)) {
      for (pn in c("var2000", "bio")) {
        g <- gate_df[gate_df$cohort == coh & gate_df$pool == pn, ]
        if (nrow(g) == 0) next
        report_lines <- c(report_lines, sprintf(
          "| %s | %s | %d | %.1f | %.1f | %.1f | %.1f | %d | %d |",
          coh, pn, nrow(g), mean(g$sets_tested), mean(g$sets_q05),
          mean(g$sets_q20), mean(g$genes_q20), sum(g$genes_q20 == 0),
          sum(g$null_selftest_pass > 2)))
      }
    }
    report_lines <- c(report_lines, "",
      "Null self-test: one fixed label permutation analysed as observed;",
      "a correct calibration certifies ~0-2 sets. Violations = splits",
      "where the null run itself certified >2 sets.", "")
  }

  report_lines <- c(report_lines,
    "## Caveats",
    "",
    "- Incumbent GS anchors used the evidence-ratio calibration null; the",
    "  new arms used calibration_mode=percentile (the null is 97% of fit",
    "  cost). GS_uncalibrated is the closest recipe anchor.",
    "- Random baselines: 3 draws per split/pool/k (new arms); incumbent",
    "  deltas use their own saved Random method.",
    "- Gate fallbacks counted per arm-k row; fallback means the gate was",
    "  empty and the ungrouped ranking was used (defined behaviour).")

  writeLines(report_lines, file.path(out_root, "REPORT.md"))
  cat(sprintf("[report] wrote %s\n", file.path(out_root, "REPORT.md")))
}
