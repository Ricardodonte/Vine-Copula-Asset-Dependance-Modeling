pstd_manual <- function(z, nu) pt(z * sqrt(nu / (nu - 2)), df = nu)
qstd_manual <- function(u, nu) qt(u, df = nu) * sqrt((nu - 2) / nu)
dstd_log <- function(z, nu) {
  lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(pi * (nu - 2)) -
    (nu + 1) / 2 * log1p(z^2 / (nu - 2))
}
std_abs_mean <- function(nu) {
  2 * sqrt(nu - 2) * gamma((nu + 1) / 2) / (sqrt(pi) * (nu - 1) * gamma(nu / 2))
}

egarch_filter_manual <- function(x, mu, phi, theta, omega, alpha, gamma, beta, nu, ln_s2_0) {
  n <- length(x)
  p <- length(phi)
  q <- length(theta)
  eps <- numeric(n)
  ln_s2 <- numeric(n)
  z <- numeric(n)
  Ez <- std_abs_mean(nu)
  y <- x - mu
  for (t in seq_len(n)) {
    m <- 0
    if (p > 0) for (i in seq_len(p)) if (t - i >= 1) m <- m + phi[i] * y[t - i]
    if (q > 0) for (j in seq_len(q)) if (t - j >= 1) m <- m + theta[j] * eps[t - j]
    eps[t] <- y[t] - m
    if (t == 1) {
      ln_s2[t] <- ln_s2_0
    } else {
      ln_s2[t] <- omega + alpha * z[t - 1] + gamma * (abs(z[t - 1]) - Ez) + beta * ln_s2[t - 1]
    }
    z[t] <- eps[t] / exp(0.5 * ln_s2[t])
  }
  list(eps = eps, ln_s2 = ln_s2, z = z)
}

.unpack <- function(par, p, q) {
  list(mu = par[1],
       phi = if (p > 0) par[1 + seq_len(p)] else numeric(0),
       theta = if (q > 0) par[1 + p + seq_len(q)] else numeric(0),
       omega = par[2 + p + q], alpha = par[3 + p + q], gamma = par[4 + p + q],
       beta = par[5 + p + q], nu = par[6 + p + q])
}

egarch_nll <- function(par, x, p, q, ln_s2_0) {
  k <- .unpack(par, p, q)
  f <- egarch_filter_manual(x, k$mu, k$phi, k$theta, k$omega, k$alpha, k$gamma, k$beta, k$nu, ln_s2_0)
  ll <- sum(dstd_log(f$z, k$nu) - 0.5 * f$ln_s2)
  if (!is.finite(ll)) return(1e10)
  -ll
}

fit_egarch_manual <- function(x, p = 0, q = 0, start = NULL) {
  x <- as.numeric(x)
  ln_s2_0 <- log(mean((x - mean(x))^2))
  if (is.null(start)) {
    start <- c(mean(x), rep(0, p), rep(0, q),
               ln_s2_0 * (1 - 0.95), 0, 0.1, 0.95, 6)
  }
  lower <- c(-Inf, rep(-0.99, p), rep(-0.99, q), -Inf, -Inf, -Inf, -0.9999, 2.05)
  upper <- c(Inf, rep(0.99, p), rep(0.99, q), Inf, Inf, Inf, 0.9999, 100)
  opt <- nlminb(start, egarch_nll, x = x, p = p, q = q, ln_s2_0 = ln_s2_0,
                lower = lower, upper = upper, control = list(eval.max = 2000, iter.max = 1000))
  if (opt$convergence != 0 || opt$objective > 1e9) {
    start2 <- start
    start2[length(start2)] <- 4
    start2[length(start2) - 1] <- 0.9
    opt2 <- nlminb(start2, egarch_nll, x = x, p = p, q = q, ln_s2_0 = ln_s2_0,
                   lower = lower, upper = upper, control = list(eval.max = 2000, iter.max = 1000))
    if (opt2$objective < opt$objective) opt <- opt2
  }
  k <- .unpack(opt$par, p, q)
  f <- egarch_filter_manual(x, k$mu, k$phi, k$theta, k$omega, k$alpha, k$gamma, k$beta, k$nu, ln_s2_0)
  n <- length(x)
  y <- x - k$mu
  m_next <- k$mu
  if (p > 0) for (i in seq_len(p)) m_next <- m_next + k$phi[i] * y[n + 1 - i]
  if (q > 0) for (j in seq_len(q)) m_next <- m_next + k$theta[j] * f$eps[n + 1 - j]
  ln_s2_next <- k$omega + k$alpha * f$z[n] + k$gamma * (abs(f$z[n]) - std_abs_mean(k$nu)) +
    k$beta * f$ln_s2[n]
  structure(list(par = opt$par, coef = k, p = p, q = q, loglik = -opt$objective,
                 z = f$z, sigma = exp(0.5 * f$ln_s2), eps = f$eps, x = x, ln_s2_0 = ln_s2_0,
                 mean_next = m_next, sigma_next = exp(0.5 * ln_s2_next),
                 convergence = opt$convergence), class = "egarch_fit")
}

egarch_update_manual <- function(fit, x_new) {
  x_all <- c(fit$x, as.numeric(x_new))
  k <- fit$coef
  f <- egarch_filter_manual(x_all, k$mu, k$phi, k$theta, k$omega, k$alpha, k$gamma, k$beta, k$nu, fit$ln_s2_0)
  n <- length(x_all)
  y <- x_all - k$mu
  m_next <- k$mu
  if (fit$p > 0) for (i in seq_len(fit$p)) m_next <- m_next + k$phi[i] * y[n + 1 - i]
  if (fit$q > 0) for (j in seq_len(fit$q)) m_next <- m_next + k$theta[j] * f$eps[n + 1 - j]
  ln_s2_next <- k$omega + k$alpha * f$z[n] + k$gamma * (abs(f$z[n]) - std_abs_mean(k$nu)) +
    k$beta * f$ln_s2[n]
  list(mean_next = m_next, sigma_next = exp(0.5 * ln_s2_next))
}

egarch_se_manual <- function(fit) {
  if (fit$convergence != 0) {
    warning("The EGARCH fit did not converge. Standard errors are unavailable.", call. = FALSE)
    return(rep(NA_real_, length(fit$par)))
  }
  nll <- function(th) egarch_nll(th, fit$x, fit$p, fit$q, fit$ln_s2_0)
  # The default step can push beta above 1 and destabilize the filter.
  H <- numDeriv::hessian(nll, fit$par, method.args = list(d = 1e-4, eps = 1e-5))
  tryCatch(sqrt(diag(chol2inv(chol(H)))), error = function(e) rep(NA_real_, length(fit$par)))
}
