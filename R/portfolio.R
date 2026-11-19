library(quadprog)
library(Rglpk)

W_MAX <- 0.5

min_var_weights <- function(S, wmax = W_MAX) {
  d <- ncol(S)
  S <- S + diag(1e-10, d)
  Amat <- cbind(rep(1, d), diag(d), -diag(d))
  bvec <- c(1, rep(0, d), rep(-wmax, d))
  sol <- solve.QP(Dmat = 2 * S, dvec = rep(0, d), Amat = Amat, bvec = bvec, meq = 1)
  w <- pmax(sol$solution, 0); w / sum(w)
}

# Rockafellar-Uryasev linear program. alpha = 0.05 targets 95% CVaR.
min_cvar_weights <- function(R, alpha = 0.05, wmax = W_MAX) {
  n <- nrow(R); d <- ncol(R)
  # Minimize VaR + mean(u) / alpha, where u is the excess loss.
  obj <- c(rep(0, d), 1, rep(1 / (alpha * n), n))
  # Each scenario imposes R_i w + VaR + u_i >= 0. The last row sets the budget.
  ii <- c(rep(seq_len(n), d), seq_len(n), seq_len(n), rep(n + 1, d))
  jj <- c(rep(seq_len(d), each = n), rep(d + 1, n), d + 1 + seq_len(n), seq_len(d))
  vv <- c(as.vector(R), rep(1, n), rep(1, n), rep(1, d))
  mat <- slam::simple_triplet_matrix(ii, jj, vv, nrow = n + 1, ncol = d + 1 + n)
  dir <- c(rep(">=", n), "==")
  rhs <- c(rep(0, n), 1)
  bounds <- list(lower = list(ind = c(seq_len(d), d + 1), val = c(rep(0, d), -Inf)),
                 upper = list(ind = seq_len(d), val = rep(wmax, d)))
  sol <- Rglpk_solve_LP(obj, mat, dir, rhs, bounds = bounds, max = FALSE)
  w <- pmax(sol$solution[seq_len(d)], 0); w / sum(w)
}

scenarios_from_uniforms <- function(U_sim, fits, pit = "correct") {
  sapply(seq_along(fits), function(j) {
    f <- fits[[j]]; nu <- f$coef$nu
    q <- if (pit == "correct") qstd(U_sim[, j], nu, fits[[j]]$backend) else qt(U_sim[, j], df = nu)
    f$mean_next + f$sigma_next * q
  })
}

perf_metrics <- function(r, rf = 0.02, turnover = NULL) {
  r <- r[is.finite(r)]
  ann <- 252
  wealth <- cumprod(1 + r)
  dd <- wealth / cummax(wealth) - 1
  ar <- prod(1 + r)^(ann / length(r)) - 1
  vol <- sd(r) * sqrt(ann)
  rf_d <- (1 + rf)^(1 / ann) - 1
  down <- sqrt(mean(pmin(r - rf_d, 0)^2)) * sqrt(ann)
  var95 <- -quantile(r, 0.05)
  out <- c(`Total return (%)` = 100 * (wealth[length(wealth)] - 1),
           `Annualized return (%)` = 100 * ar,
           `Annualized volatility (%)` = 100 * vol,
           Sharpe = (ar - rf) / vol,
           Sortino = (ar - rf) / down,
           `Max drawdown (%)` = 100 * min(dd),
           Calmar = ar / abs(min(dd)),
           `Daily VaR 95% (%)` = 100 * unname(var95),
           `Daily CVaR 95% (%)` = 100 * -mean(r[r <= -var95]))
  if (!is.null(turnover)) out <- c(out, `Annual turnover (x)` = turnover)
  out
}

max_drawdown <- function(r) { w <- cumprod(1 + r); min(w / cummax(w) - 1) }
sharpe_ann <- function(r, rf = 0.02) {
  ar <- prod(1 + r)^(252 / length(r)) - 1
  (ar - rf) / (sd(r) * sqrt(252))
}

paired_block_bootstrap <- function(r1, r2, B = 2000, block = 20, seed = 1) {
  set.seed(seed)
  n <- length(r1)
  stat <- function(a, b) c(
    dSharpe = sharpe_ann(a) - sharpe_ann(b),
    dCVaR = (-mean(a[a <= quantile(a, 0.05)])) - (-mean(b[b <= quantile(b, 0.05)])),
    dMDD = max_drawdown(a) - max_drawdown(b))
  obs <- stat(r1, r2)
  nb <- ceiling(n / block)
  boot <- replicate(B, {
    starts <- sample.int(n, nb, replace = TRUE)
    idx <- as.vector(sapply(starts, function(s) ((s - 1 + 0:(block - 1)) %% n) + 1))[1:n]
    stat(r1[idx], r2[idx])
  })
  # Center the bootstrap distribution under the null of equality.
  p <- sapply(seq_along(obs), function(k) {
    cen <- boot[k, ] - mean(boot[k, ])
    mean(abs(cen) >= abs(obs[k]))
  })
  data.frame(stat = names(obs), diff = unname(obs), p_value = p)
}
