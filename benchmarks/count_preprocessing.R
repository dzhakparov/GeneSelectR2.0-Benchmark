# Fold-local normalization for integer RNA-seq count matrices.
#
# The count matrix is genes by samples. Gene eligibility and the TMM reference
# are estimated from the training samples. Each held-out sample is normalized
# independently against the fixed training reference. This prevents one test
# sample, or the composition of the complete test fold, from changing either
# the training matrix or another held-out sample.

normalise_count_split_fixed_reference <- function(
    raw_count_matrix, train_indices, test_indices = NULL,
    min_count_per_gene = 10L, min_samples_per_gene = 10L,
    prior_count = 1) {

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Fold-local count normalization requires the edgeR package.")
  }
  if (!is.matrix(raw_count_matrix) || nrow(raw_count_matrix) == 0L ||
      ncol(raw_count_matrix) == 0L) {
    stop("raw_count_matrix must be a non-empty genes-by-samples matrix.")
  }
  if (any(!is.finite(raw_count_matrix)) || any(raw_count_matrix < 0) ||
      any(abs(raw_count_matrix - round(raw_count_matrix)) > 1e-8)) {
    stop("Count input must contain finite non-negative integers.")
  }
  if (length(train_indices) == 0L || anyDuplicated(train_indices) ||
      any(train_indices < 1L | train_indices > ncol(raw_count_matrix))) {
    stop("train_indices must identify distinct samples in raw_count_matrix.")
  }
  if (!is.null(test_indices) &&
      (anyDuplicated(test_indices) ||
       any(test_indices < 1L | test_indices > ncol(raw_count_matrix)) ||
       length(intersect(train_indices, test_indices)) > 0L)) {
    stop("test_indices must be distinct from the training samples.")
  }

  train_counts <- raw_count_matrix[, train_indices, drop = FALSE]

  # The fixed abundance rule is outcome-independent and is fitted within the
  # training fold. The minimum sample requirement is reduced only when the
  # complete training fold contains fewer samples than the configured value.
  required_samples <- min(as.integer(min_samples_per_gene), ncol(train_counts))
  keep <- rowSums(train_counts >= min_count_per_gene) >= required_samples
  if (!any(keep)) {
    stop("No genes passed the training-fold count filter.")
  }
  train_counts <- train_counts[keep, , drop = FALSE]
  test_counts <- if (is.null(test_indices)) NULL else
    raw_count_matrix[keep, test_indices, drop = FALSE]

  library_sizes <- colSums(train_counts)
  if (any(library_sizes <= 0)) {
    stop("A training sample has a zero count-library size after filtering.")
  }

  # This reproduces edgeR's upper-quartile reference choice using training
  # samples only. The selected reference remains fixed for all held-out samples.
  upper_quartile_cpm <- apply(train_counts, 2, function(counts) {
    stats::quantile(counts / sum(counts) * 1e6, 0.75,
                    names = FALSE, type = 7)
  })
  reference_index <- which.min(
    abs(upper_quartile_cpm - mean(upper_quartile_cpm))
  )

  train_dge <- edgeR::DGEList(counts = train_counts)
  train_dge <- edgeR::calcNormFactors(train_dge,
                                      refColumn = reference_index)
  train_factors <- train_dge$samples$norm.factors
  train_logcpm <- edgeR::cpm(train_dge, log = TRUE,
                             prior.count = prior_count)

  if (is.null(test_counts)) {
    return(list(
      train = t(train_logcpm),
      test = NULL,
      kept_genes = rownames(train_counts),
      reference_sample = colnames(train_counts)[reference_index]
    ))
  }

  # Pairwise TMM factors have product one. Their ratio is invariant to that
  # rescaling and places the held-out sample on the fixed training-reference
  # scale. No other held-out sample enters this calculation.
  test_logcpm <- vapply(seq_len(ncol(test_counts)), function(test_index) {
    pair_dge <- edgeR::DGEList(counts = cbind(
      reference = train_counts[, reference_index],
      held_out = test_counts[, test_index]
    ))
    pair_dge <- edgeR::calcNormFactors(pair_dge, refColumn = 1L)
    pair_ratio <- pair_dge$samples$norm.factors[2] /
      pair_dge$samples$norm.factors[1]
    held_out_factor <- train_factors[reference_index] * pair_ratio
    pair_dge$samples$norm.factors <- c(
      train_factors[reference_index], held_out_factor
    )
    edgeR::cpm(pair_dge, log = TRUE,
               prior.count = prior_count)[, 2]
  }, numeric(nrow(train_counts)))
  rownames(test_logcpm) <- rownames(train_counts)
  colnames(test_logcpm) <- colnames(test_counts)

  list(
    train = t(train_logcpm),
    test = t(test_logcpm),
    kept_genes = rownames(train_counts),
    reference_sample = colnames(train_counts)[reference_index]
  )
}
