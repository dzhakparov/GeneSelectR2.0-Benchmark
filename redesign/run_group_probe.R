#!/usr/bin/env Rscript
# ==============================================================================
#  Group-based selection probe on IMvigor210 fold 1 (cached matrix).
#
#  Two approaches, run separately, then one combined report:
#
#    Rscript redesign/run_group_probe.R gvs       # approach 1: data-driven
#                                                 #   groups (trex+GVS)
#    Rscript redesign/run_group_probe.R hallmark  # approach 2: pathway groups
#                                                 #   (Hallmark + perm null)
#    Rscript redesign/run_group_probe.R combined  # combined report, reads the
#                                                 #   two results above
#
#  Approach 1 (gvs): trex+GVS at tFDR=0.2, corr_max=0.5, K=40. This is the
#  fair-test configuration: 2x the experiments of the first failed run, group
#  level. If it still certifies nothing, data-driven group dummy-competition
#  is dead on this data.
#
#  Approach 2 (hallmark): module-level association test. Statistic per set:
#  mean |two-sample t| over member genes present in the 2000-gene candidate
#  pool. Null: 1000 label permutations. BH at 0.05/0.1/0.2. The null is the
#  LABEL permutation (tests association with outcome), which is the null an
#  FDR claim needs. Built-in self-test: one fixed permutation is analysed as
#  if it were the observed data; if that "null run" certifies many sets, the
#  calibration is broken and the real result is not to be trusted.
#
#  Everything reads redesign/results/trex_sanity/fold1_train.rds (n=153,
#  p=2000, 20 responders) and writes to redesign/results/group_probe/.
# ==============================================================================

mode <- {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[1] %in% c("gvs", "hallmark", "combined"))
    args[1] else stop("usage: run_group_probe.R [gvs|hallmark|combined]")
}

out_dir <- file.path("redesign", "results", "group_probe")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cached <- readRDS(file.path("redesign", "results", "trex_sanity",
                            "fold1_train.rds"))
X <- cached$X
y <- cached$y
cat(sprintf("[%s] data: n=%d p=%d responders=%d\n", mode,
            nrow(X), ncol(X), sum(y)))


# ==============================================================================
#  Approach 1: trex+GVS (data-driven groups, dummy competition)
# ==============================================================================

run_gvs <- function() {
  suppressPackageStartupMessages(library(TRexSelector))
  t0 <- proc.time()[["elapsed"]]
  fit <- trex(X, y, tFDR = 0.2, K = 30, method = "trex+GVS",
              corr_max = 0.5, seed = 1,
              parallel_process = TRUE, parallel_max_cores = 7,
              verbose = FALSE)
  elapsed <- proc.time()[["elapsed"]] - t0
  sel <- which(fit$selected_var > 0)
  top_phi <- sort(fit$Phi_prime, decreasing = TRUE)[1:15]

  write.csv(data.frame(gene = colnames(X)[sel], index = sel),
            file.path(out_dir, "gvs_selected.csv"), row.names = FALSE)
  diag <- data.frame(
    setting = c("tFDR", "K", "corr_max", "T_stop", "v_thresh",
                "n_selected", "elapsed_sec"),
    value = c(0.2, 30, 0.5, fit$T_stop, fit$v_thresh,
              length(sel), round(elapsed, 1)))
  write.csv(diag, file.path(out_dir, "gvs_diagnostics.csv"),
            row.names = FALSE)
  phi_sorted <- sort(fit$Phi_prime, decreasing = TRUE)
  n_show <- min(100, length(phi_sorted))
  write.csv(data.frame(
    gene = colnames(X)[order(fit$Phi_prime, decreasing = TRUE)[
      seq_len(n_show)]],
    phi = round(phi_sorted[seq_len(n_show)], 4)),
    file.path(out_dir, "gvs_top_phi.csv"), row.names = FALSE)

  cat(sprintf("[gvs] %.0f s | T_stop=%s v_thresh=%.3f | selected %d genes\n",
              elapsed, fit$T_stop, fit$v_thresh, length(sel)))
  cat(sprintf("[gvs] top Phi: %s\n",
              paste(round(top_phi, 3), collapse = " ")))
  if (length(sel) > 0) {
    cat(sprintf("[gvs] genes: %s\n",
                paste(colnames(X)[sel][1:min(20, length(sel))],
                      collapse = ", ")))
  }
}


# ==============================================================================
#  Approach 2: Hallmark module test with label-permutation null
# ==============================================================================

#  Two-sample t-statistic for every gene at once (equal variance).
t_stats <- function(Xmat, ybin) {
  n1 <- sum(ybin == 1); n0 <- sum(ybin == 0)
  m1 <- colMeans(Xmat[ybin == 1, , drop = FALSE])
  m0 <- colMeans(Xmat[ybin == 0, , drop = FALSE])
  v1 <- apply(Xmat[ybin == 1, , drop = FALSE], 2, var)
  v0 <- apply(Xmat[ybin == 0, , drop = FALSE], 2, var)
  sp <- sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
  sp[sp == 0 | !is.finite(sp)] <- NA_real_
  (m1 - m0) / (sp * sqrt(1 / n1 + 1 / n0))
}

hallmark_sets <- function(pool_genes) {
  ms <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e)
      msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  name_col <- intersect(c("gs_name", "gs_cat"), names(ms))[1]
  sym_col <- intersect(c("gene_symbol", "human_gene_symbol"), names(ms))[1]
  sets <- split(ms[[sym_col]], ms[[name_col]])
  sets <- lapply(sets, function(g) intersect(g, pool_genes))
  sets[vapply(sets, length, integer(1)) >= 10]
}

run_hallmark <- function(B = 1000, seed = 7) {
  sets <- hallmark_sets(colnames(X))
  cat(sprintf("[hallmark] %d sets with >=10 members in pool\n",
              length(sets)))
  stopifnot(length(sets) >= 20)   # sanity: expect ~40-50 Hallmark sets

  t_obs <- abs(t_stats(X, y))
  set_index <- lapply(sets, function(g) match(g, colnames(X)))
  stat_obs <- vapply(set_index, function(idx)
    mean(t_obs[idx], na.rm = TRUE), numeric(1))

  set.seed(seed)
  null_stats <- matrix(NA_real_, nrow = B, ncol = length(sets))
  colnames(null_stats) <- names(sets)
  for (b in seq_len(B)) {
    t_b <- abs(t_stats(X, sample(y)))
    null_stats[b, ] <- vapply(set_index, function(idx)
      mean(t_b[idx], na.rm = TRUE), numeric(1))
  }

  p_perm <- (1 + colSums(sweep(null_stats, 2, stat_obs, `>=`))) / (B + 1)
  q_bh <- p.adjust(p_perm, method = "BH")

  res <- data.frame(set = names(sets), n_in_pool = lengths(sets),
                    stat_obs = round(stat_obs, 4),
                    null_mean = round(colMeans(null_stats), 4),
                    p_perm = p_perm, q_BH = q_bh,
                    pass_0.05 = q_bh <= 0.05,
                    pass_0.10 = q_bh <= 0.10,
                    pass_0.20 = q_bh <= 0.20,
                    row.names = NULL)
  res <- res[order(res$q_BH), ]
  write.csv(res, file.path(out_dir, "hallmark_sets.csv"), row.names = FALSE)

  #  Self-test: one fixed permutation treated as observed. Under a correct
  #  calibration essentially nothing should pass BH 0.2 here.
  set.seed(123)
  t_null_run <- abs(t_stats(X, sample(y)))
  stat_null_run <- vapply(set_index, function(idx)
    mean(t_null_run[idx], na.rm = TRUE), numeric(1))
  p_null_run <- (1 + colSums(sweep(null_stats, 2, stat_null_run,
                                   `>=`))) / (B + 1)
  q_null_run <- p.adjust(p_null_run, method = "BH")
  n_null_pass <- sum(q_null_run <= 0.2)
  cat(sprintf("[hallmark] SELF-TEST: null-run sets passing BH 0.2: %d ",
              n_null_pass))
  cat(sprintf("(expect ~0-2; if large, calibration is broken)\n"))

  cat(sprintf("[hallmark] sets passing BH 0.05/0.10/0.20: %d / %d / %d\n",
              sum(res$pass_0.05), sum(res$pass_0.10), sum(res$pass_0.20)))
  cat("[hallmark] top 10 sets:\n")
  print(head(res[, c("set", "n_in_pool", "stat_obs", "null_mean",
                     "p_perm", "q_BH")], 10), row.names = FALSE)

  meta <- data.frame(setting = c("B_permutations", "seed", "n_sets",
                                 "null_run_pass_BH0.2"),
                     value = c(B, seed, length(sets), n_null_pass))
  write.csv(meta, file.path(out_dir, "hallmark_meta.csv"),
            row.names = FALSE)
}


# ==============================================================================
#  Combined report: reads both approaches' saved results
# ==============================================================================

run_combined <- function() {
  gvs_sel_path <- file.path(out_dir, "gvs_selected.csv")
  gvs_diag_path <- file.path(out_dir, "gvs_diagnostics.csv")
  hallmark_path <- file.path(out_dir, "hallmark_sets.csv")
  meta_path <- file.path(out_dir, "hallmark_meta.csv")
  missing <- c(gvs_diag_path, hallmark_path, meta_path)[
    !file.exists(c(gvs_diag_path, hallmark_path, meta_path))]
  if (length(missing) > 0) {
    stop("run the separate approaches first; missing: ",
         paste(missing, collapse = ", "))
  }

  gvs_diag <- read.csv(gvs_diag_path)
  gvs_sel <- if (file.exists(gvs_sel_path)) read.csv(gvs_sel_path)$gene
             else character(0)
  hm <- read.csv(hallmark_path)
  meta <- read.csv(meta_path)

  lines <- c(
    "GROUP-BASED SELECTION PROBE — COMBINED REPORT",
    sprintf("data: fold1 train, n=%d p=%d responders=%d",
            nrow(X), ncol(X), sum(y)),
    "",
    "Approach 1 — trex+GVS (data-driven groups, dummy competition):",
    sprintf("  tFDR=0.2 K=30 corr_max=0.5 | v_thresh=%s | selected %s genes",
            gvs_diag$value[gvs_diag$setting == "v_thresh"],
            gvs_diag$value[gvs_diag$setting == "n_selected"]),
    sprintf("  genes: %s", ifelse(length(gvs_sel) > 0,
                                  paste(gvs_sel, collapse = ", "),
                                  "(none)")),
    "",
    "Approach 2 — Hallmark module test (1000 label permutations):",
    sprintf("  sets tested: %s | null-run self-test passes (BH 0.2): %s",
            meta$value[meta$setting == "n_sets"],
            meta$value[meta$setting == "null_run_pass_BH0.2"]),
    sprintf("  sets passing BH 0.05 / 0.10 / 0.20: %d / %d / %d",
            sum(hm$pass_0.05), sum(hm$pass_0.10), sum(hm$pass_0.20)),
    "",
    "  top 10 Hallmark sets by q-value:")
  top10 <- head(hm[, c("set", "n_in_pool", "stat_obs", "null_mean",
                       "p_perm", "q_BH")], 10)
  lines <- c(lines, capture.output(print(top10, row.names = FALSE)))

  #  Cross-check: do the GVS-selected genes (if any) sit inside the
  #  significant Hallmark sets?
  sig_sets <- hm$set[hm$pass_0.20]
  lines <- c(lines, "", sprintf(
    "Cross-check: %d GVS genes; %d Hallmark sets pass BH 0.2.",
    length(gvs_sel), length(sig_sets)))

  report <- paste(lines, collapse = "\n")
  cat(report, "\n")
  writeLines(report, file.path(out_dir, "combined_report.txt"))
  cat(sprintf("\n[combined] report written to %s\n",
              file.path(out_dir, "combined_report.txt")))
}


switch(mode,
       gvs = run_gvs(),
       hallmark = run_hallmark(),
       combined = run_combined())
