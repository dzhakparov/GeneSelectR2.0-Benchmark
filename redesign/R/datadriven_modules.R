# ==============================================================================
#  Data-driven modules: correlation clusters from the TRAINING matrix only
#  (added 2026-08-22). Hallmark sets are fixed a priori; these are learned
#  per split. Same z-score + soft-boost machinery as the Hallmark prior,
#  so the two sources can be combined by taking the max boost per gene.
#
#  Clustering: 1 - |cor| distance, average linkage, cutree at fixed k.
#  dynamicTreeCut was tried first and degenerates on weakly correlated
#  data (1 module on GSE69683 at deepSplit 4); fixed k = min(200, p/10)
#  gives ~10 usable modules on every dataset after size filtering.
#  Modules smaller than min_size or larger than max_size are dropped:
#  tiny modules give noisy z-scores, huge ones are the meta-gene and
#  carry no specific signal. Clustering is unsupervised -- labels enter
#  only through module_scores, so there is no leak into module structure.
# ==============================================================================

datadriven_modules <- function(X, min_size = 10, max_size = 200,
                               k = min(200L, max(2L, ncol(X) %/% 10L))) {
  cc <- cor(X)
  cc[!is.finite(cc)] <- 0
  hc <- hclust(as.dist(1 - abs(cc)), method = "average")
  lab <- cutree(hc, k = k)
  sets <- split(names(lab), paste0("DC", lab))
  sizes <- vapply(sets, length, integer(1))
  sets[sizes >= min_size & sizes <= max_size]
}
