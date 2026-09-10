# ==============================================================================
#  Rank-ensemble combinators for the post-benchmark ensemble evaluation.
#
#  rank_average(rankings, pool): mean rank across rankings; a gene absent
#    from a ranking (e.g. a gated ranking that covers only certified genes)
#    gets length(that ranking) + its position in the first ranking that
#    contains it -- penalised but not dropped.
#  kswitch(k, k_cut): pick the ranking per panel size (predfirst_raw for
#    small panels, GS_full_grouped above the cut). This rule was constructed
#    after the development results were examined.
# ==============================================================================

rank_average <- function(rankings, pool) {
  n <- length(pool)
  rank_sum <- setNames(rep(0, n), pool)
  for (rk in rankings) {
    r <- setNames(seq_along(rk), rk)
    missing <- setdiff(pool, rk)
    if (length(missing) > 0) {
      #  Penalised rank: after every gene the ranking did cover.
      r[missing] <- length(rk) + match(missing, pool) / (n + 1)
    }
    rank_sum <- rank_sum + r[pool]
  }
  names(sort(rank_sum))
}

kswitch_ranking <- function(k, k_cut, ranking_small, ranking_large) {
  if (k <= k_cut) ranking_small else ranking_large
}
