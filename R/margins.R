egarch_backend <- function() {
  if (isTRUE(getOption("vine.use_rugarch", TRUE))) "rugarch" else "manual"
}

pstd <- function(z, nu, backend = egarch_backend()) {
  if (backend == "manual") return(pstd_manual(z, nu))
  rugarch::pdist("std", q = z, shape = nu)
}

qstd <- function(u, nu, backend = egarch_backend()) {
  if (backend == "manual") return(qstd_manual(u, nu))
  rugarch::qdist("std", p = u, shape = nu)
}

fit_egarch <- function(x, p = 0, q = 0) {
  backend <- egarch_backend()
  fit <- if (backend == "manual") fit_egarch_manual(x, p, q) else fit_egarch_rugarch(x, p, q)
  fit$backend <- backend
  fit
}

egarch_update <- function(fit, x_new) {
  if (fit$backend == "manual") return(egarch_update_manual(fit, x_new))
  egarch_update_rugarch(fit, x_new)
}

egarch_se <- function(fit) {
  if (fit$backend == "manual") return(egarch_se_manual(fit))
  fit$se
}

egarch_spec <- function(p, q, fixed = list()) {
  rugarch::ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(p, q), include.mean = TRUE),
    distribution.model = "std",
    fixed.pars = fixed
  )
}

fit_egarch_rugarch <- function(x, p = 0, q = 0) {
  fit <- rugarch::ugarchfit(egarch_spec(p, q), data = x, solver = "hybrid",
                           solver.control = list(rseed = 42L))
  if (rugarch::convergence(fit) != 0L) {
    stop("rugarch : l'ajustement EGARCH n'a pas convergé.", call. = FALSE)
  }

  par <- rugarch::coef(fit)
  forecast <- rugarch::ugarchforecast(fit, n.ahead = 1)
  covariance <- rugarch::vcov(fit)
  se <- rep(NA_real_, length(par))
  if (!is.null(covariance)) {
    variance <- diag(covariance)
    valid <- is.finite(variance) & variance > 0
    se[valid] <- sqrt(variance[valid])
  }

  order <- match(param_names(p, q, labels = FALSE), names(par))
  list(par = par[order], coef = list(nu = unname(par[["shape"]])), p = p, q = q,
       z = as.numeric(rugarch::residuals(fit, standardize = TRUE)),
       sigma = as.numeric(rugarch::sigma(fit)), loglik = rugarch::likelihood(fit),
       mean_next = as.numeric(rugarch::fitted(forecast))[1],
       sigma_next = as.numeric(rugarch::sigma(forecast))[1], x = as.numeric(x), se = se[order])
}

egarch_update_rugarch <- function(fit, x_new) {
  n <- length(x_new)
  if (n == 0L) return(list(mean_next = fit$mean_next, sigma_next = fit$sigma_next))

  # Même initialisation que sur la fenêtre d'ajustement.
  spec <- egarch_spec(fit$p, fit$q, fixed = as.list(fit$par))
  forecast <- rugarch::ugarchforecast(spec, data = c(fit$x, x_new),
                                     n.ahead = 1, n.roll = n, out.sample = n)
  list(mean_next = as.numeric(rugarch::fitted(forecast))[n + 1L],
       sigma_next = as.numeric(rugarch::sigma(forecast))[n + 1L])
}

param_names <- function(p, q, labels = TRUE) {
  parameters <- c("mu", if (p > 0) paste0("ar", seq_len(p)), if (q > 0) paste0("ma", seq_len(q)),
                  "omega", "alpha1", "gamma1", "beta1", "shape")
  if (labels) {
    parameters[match(c("alpha1", "gamma1", "shape"), parameters)] <-
      c("alpha1 (signe)", "gamma1 (taille)", "shape (nu)")
  }
  parameters
}

fit_margins <- function(R, orders, pit = c("correct", "original")) {
  pit <- match.arg(pit)
  R <- as.matrix(R)
  fits <- lapply(seq_len(ncol(R)), function(j) {
    fit_egarch(R[, j], p = orders[j, 1], q = orders[j, 2])
  })
  names(fits) <- colnames(R)
  U <- sapply(fits, function(f) {
    if (pit == "correct") pstd(f$z, f$coef$nu, f$backend) else pt(f$z, df = f$coef$nu)
  })
  colnames(U) <- colnames(R)
  list(fits = fits, U = U)
}
