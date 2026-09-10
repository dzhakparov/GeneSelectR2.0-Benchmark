# ==============================================================================
#  Biological prior machinery: Hallmark sets, module gate, candidate pools.
#
#  Design rule for this redesign (user requirement): biology enters ONLY as a
#  prior - as a candidate-pool constraint and as a pre-selection module gate.
#  There is NO post-hoc biological scoring anywhere in the new arms, so the
#  annotation used here can never re-rank genes after the fact. That removes
#  the circularity the old three-pillar design was criticised for.
#
#  Two biological priors, both Hallmark-based (MSigDB "H" collection):
#    1. bio pool: candidate pool = genes annotated in ANY Hallmark set.
#    2. module gate: certify individual Hallmark sets on the training data
#       (mean |t| vs label-permutation null, BH), restrict the pool to genes
#       in certified sets BEFORE the stability+utility ranking runs.
#
#  The gate's null is label permutation, which is the null an association
#  claim needs; pool composition does not affect it.
# ==============================================================================


#  All Hallmark sets as a named list of gene-symbol vectors. Handles both the
#  old msigdbr API (category = "H") and the new one (collection = "H").
hallmark_sets_all <- function() {
  ms <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e)
      msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  name_col <- intersect(c("gs_name", "gs_cat"), names(ms))[1]
  sym_col  <- intersect(c("gene_symbol", "human_gene_symbol"), names(ms))[1]
  if (is.na(name_col) || is.na(sym_col)) {
    stop("msigdbr output has no recognisable set-name / symbol columns.")
  }
  split(ms[[sym_col]], ms[[name_col]])
}

#  Union of all Hallmark-annotated genes: the biological candidate-pool prior.
hallmark_union <- function() {
  unique(unlist(hallmark_sets_all(), use.names = FALSE))
}


#  Two-sample t-statistic for every gene at once (equal variance). X is
#  samples x genes, ybin is 0/1. Genes with zero pooled variance get NA and
#  are ignored by the module statistic (na.rm = TRUE there).
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


#  Module gate: certify Hallmark sets on a training matrix.
#
#  Statistic per set: mean |t| over member genes present in `colnames(X)`.
#  Null: B label permutations, recomputing the whole statistic each time.
#  p = (1 + #{null >= obs}) / (B + 1); BH across tested sets.
#
#  Self-test: one fixed permutation (seed 123) is analysed as if it were the
#  observed data and passed through the same null; the number of sets it
#  "certifies" at q <= 0.2 is returned as null_run_pass_q20. Under a correct
#  calibration this is ~0-2. If it is large the result is not to be trusted.
module_gate <- function(X, ybin, B = 1000, seed = 7, min_members = 10,
                        sets = NULL) {
  if (is.null(sets)) sets <- hallmark_sets_all()
  pool_genes <- colnames(X)
  sets <- lapply(sets, function(g) intersect(g, pool_genes))
  sets <- sets[vapply(sets, length, integer(1)) >= min_members]
  if (length(sets) == 0) {
    return(list(table = data.frame(), certified_genes = list(
      q05 = character(0), q10 = character(0), q20 = character(0)),
      null_run_pass_q20 = NA_integer_))
  }

  set_index <- lapply(sets, function(g) match(g, pool_genes))
  t_obs <- abs(t_stats(X, ybin))
  stat_obs <- vapply(set_index, function(idx)
    mean(t_obs[idx], na.rm = TRUE), numeric(1))

  set.seed(seed)
  null_stats <- matrix(NA_real_, nrow = B, ncol = length(sets))
  colnames(null_stats) <- names(sets)
  for (b in seq_len(B)) {
    t_b <- abs(t_stats(X, sample(ybin)))
    null_stats[b, ] <- vapply(set_index, function(idx)
      mean(t_b[idx], na.rm = TRUE), numeric(1))
  }

  p_perm <- (1 + colSums(sweep(null_stats, 2, stat_obs, `>=`))) / (B + 1)
  q_bh <- p.adjust(p_perm, method = "BH")

  tab <- data.frame(set = names(sets), n_in_pool = lengths(sets),
                    stat_obs = stat_obs, null_mean = colMeans(null_stats),
                    p_perm = p_perm, q_BH = q_bh, row.names = NULL)
  tab <- tab[order(tab$q_BH), ]

  genes_in <- function(level) {
    ok <- names(sets)[q_bh <= level]
    unique(unlist(sets[ok], use.names = FALSE))
  }

  #  Self-test on one fixed null realisation.
  set.seed(123)
  t_null_run <- abs(t_stats(X, sample(ybin)))
  stat_null_run <- vapply(set_index, function(idx)
    mean(t_null_run[idx], na.rm = TRUE), numeric(1))
  p_null_run <- (1 + colSums(sweep(null_stats, 2, stat_null_run,
                                   `>=`))) / (B + 1)
  null_run_pass <- sum(p.adjust(p_null_run, method = "BH") <= 0.2)

  list(table = tab,
       certified_genes = list(q05 = genes_in(0.05), q10 = genes_in(0.10),
                              q20 = genes_in(0.20)),
       null_run_pass_q20 = null_run_pass)
}


#  Candidate pools for one training matrix. All selection is fitted on the
#  training columns only; the returned column name vectors are then applied
#  unchanged to the test matrix.
#
#    var2000  : top `top_var` genes by training variance (the incumbent pool)
#    bio      : Hallmark-annotated genes with positive training variance,
#               NO variance cap (the biological-prior pool)
#    varlarge : top length(bio) genes by training variance - the size control
#               that separates "biological" from "bigger pool"
make_pools <- function(train_expr, bio_genes, top_var = 2000) {
  v <- apply(train_expr, 2, var)
  eligible <- which(is.finite(v) & v > 0)
  by_var <- eligible[order(v[eligible], decreasing = TRUE)]

  var2000 <- colnames(train_expr)[by_var[seq_len(min(top_var,
                                                     length(by_var)))]]

  bio <- intersect(colnames(train_expr)[eligible], bio_genes)

  varlarge <- colnames(train_expr)[by_var[seq_len(min(length(bio),
                                                      length(by_var)))]]

  list(var2000 = var2000, bio = bio, varlarge = varlarge)
}


#  Standardise train and test using TRAINING statistics only (same as the
#  incumbent benchmark).
standardise_split <- function(train_expr, test_expr) {
  cm <- colMeans(train_expr)
  cs <- apply(train_expr, 2, sd)
  cs[cs == 0 | is.na(cs)] <- 1
  tr <- sweep(sweep(train_expr, 2, cm, "-"), 2, cs, "/")
  te <- if (is.null(test_expr)) NULL else
    sweep(sweep(test_expr, 2, cm, "-"), 2, cs, "/")
  list(train = tr, test = te)
}


# ==============================================================================
#  Soft module prior (added 2026-08-22). The all-or-nothing gate fired on
#  too few splits to matter. This is the continuous replacement: every set
#  gets a z-score against the same label-permutation null, and gene scores
#  are multiplied by (1 + w * max(0, z)) for the best set the gene belongs
#  to. No certification threshold, no fallback branch, works at any n.
# ==============================================================================

#  Continuous module scores. Same statistic and null as module_gate
#  (mean |t| over members, B label permutations), but returns
#  z = (observed - null_mean) / null_sd per set instead of a q-value.
#  sd == 0 (degenerate null) -> z = 0, no boost.
module_scores <- function(X, ybin, B = 1000, seed = 7, min_members = 10,
                          sets = NULL) {
  if (is.null(sets)) sets <- hallmark_sets_all()
  pool_genes <- colnames(X)
  sets <- lapply(sets, function(g) intersect(g, pool_genes))
  sets <- sets[vapply(sets, length, integer(1)) >= min_members]
  if (length(sets) == 0) {
    return(list(z = setNames(numeric(0), character(0)), sets = sets,
                table = data.frame()))
  }

  set_index <- lapply(sets, function(g) match(g, pool_genes))
  stat_obs <- vapply(set_index, function(idx)
    mean(abs(t_stats(X, ybin))[idx], na.rm = TRUE), numeric(1))

  set.seed(seed)
  null_stats <- matrix(NA_real_, nrow = B, ncol = length(sets))
  colnames(null_stats) <- names(sets)
  for (b in seq_len(B)) {
    t_b <- abs(t_stats(X, sample(ybin)))
    null_stats[b, ] <- vapply(set_index, function(idx)
      mean(t_b[idx], na.rm = TRUE), numeric(1))
  }

  mu <- colMeans(null_stats)
  sd_ <- apply(null_stats, 2, sd)
  z <- (stat_obs - mu) / sd_
  z[!is.finite(z)] <- 0

  list(z = z, sets = sets,
       table = data.frame(set = names(sets), z = z, stat_obs = stat_obs,
                          null_mean = mu, row.names = NULL))
}

#  Apply the soft prior to a named score vector. multiplier per gene =
#  1 + w * max(0, best z among the gene's sets); genes in no set get 1.
#  Returns the re-scored ranking (names, best first).
soft_prior_rescore <- function(scores, sets, z, w = 0.5) {
  boost <- setNames(rep(1, length(scores)), names(scores))
  for (s in names(z)) {
    if (z[s] <= 0) next
    members <- intersect(sets[[s]], names(scores))
    if (length(members) == 0) next
    boost[members] <- pmax(boost[members], 1 + w * z[s])
  }
  boosted <- scores * boost[names(scores)]
  names(sort(boosted, decreasing = TRUE))
}
