# ==============================================================================
#  Horseshoe probit regression (added 2026-08-23).
#
#  One model fit instead of B subsample fits: every gene gets a posterior
#  inclusion probability and posterior mean effect from the SAME joint
#  posterior, so the stability score falls out of the math.
#
#  Model: Albert-Chib probit (latent z, sigma = 1) with a Carvalho-Polson-
#  Scott horseshoe on beta:
#    beta_j | lambda_j, tau ~ N(0, tau^2 lambda_j^2)
#    lambda_j, tau        ~ C+(0, 1)   (via InvGamma parameter expansion)
#
#  beta sampling uses Bhattacharya-Chakraborty-Mallick (2016) fast
#  horseshoe: O(n^2 p) per iteration instead of O(p^3).
#
#  Score for ranking: posterior mean |beta_j| x PIP_j, where
#  PIP_j = posterior fraction of draws with |beta_j| > eps. Mean |beta|
#  alone cannot separate "always medium" from "rarely huge"; PIP alone
#  ignores effect size. eps fixed at 0.05 on standardised-X scale.
# ==============================================================================

hs_probit <- function(X, ybin, n_iter = 4000, burn = 1000, eps = 0.05,
                      seed = 42, thin = 1) {
  stopifnot(all(ybin %in% c(0, 1)), nrow(X) == length(ybin))
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  Xt <- t(X)

  beta  <- rep(0, p)
  lam2  <- rep(1, p)     # lambda_j^2
  nu    <- rep(1, p)
  tau2  <- 1
  xi    <- 1
  z     <- ifelse(ybin == 1, 0.5, -0.5)

  keep <- seq(burn + 1, n_iter, by = thin)
  n_keep <- length(keep)
  draws <- matrix(NA_real_, n_keep, p, dimnames = list(NULL, colnames(X)))
  ki <- 0L

  for (it in seq_len(n_iter)) {
    #  1. latent z | beta  (truncated normal, vectorised via inverse CDF)
    eta <- as.vector(X %*% beta)
    plo <- pnorm(0, mean = eta)
    u <- runif(n)
    z <- ifelse(ybin == 1,
                qnorm(plo + u * (1 - plo), mean = eta),
                qnorm(u * plo, mean = eta))
    z[!is.finite(z)] <- sign(z[!is.finite(z)]) * 8  # clamp extreme tails

    #  2. beta | rest  (Bhattacharya fast horseshoe, D = tau2 * lam2)
    Dv <- tau2 * lam2
    u_s <- rnorm(p) * sqrt(Dv)
    delta <- rnorm(n)
    v <- as.vector(X %*% u_s) + delta
    #  Solve (X D X' + I) w = z - v via n x n system
    DXt <- Dv * Xt                      # p x n, each row scaled
    S <- X %*% DXt                      # n x n
    diag(S) <- diag(S) + 1
    w <- solve(S, z - v)
    beta <- u_s + as.vector(DXt %*% w)

    #  3. horseshoe hyperparameters (parameter-expanded Gibbs)
    lam2 <- 1 / rgamma(p, 1, rate = 1 / nu + beta^2 / (2 * tau2))
    nu   <- 1 / rgamma(p, 1, rate = 1 + 1 / lam2)
    tau2 <- 1 / rgamma(1, (p + 1) / 2, rate = 1 / xi + sum(beta^2 / lam2) / 2)
    xi   <- 1 / rgamma(1, 1, rate = 1 + 1 / tau2)

    if (it %in% keep) { ki <- ki + 1L; draws[ki, ] <- beta }
  }

  mean_abs <- colMeans(abs(draws))
  pip <- colMeans(abs(draws) > eps)
  score <- mean_abs * pip
  names(score) <- colnames(X)

  list(ranking = names(sort(score, decreasing = TRUE)),
       score = sort(score, decreasing = TRUE),
       mean_abs = sort(mean_abs, decreasing = TRUE),
       pip = sort(pip, decreasing = TRUE),
       draws = draws)
}
