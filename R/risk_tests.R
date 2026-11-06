# Kupiec (1995) tests VaR exceedance frequency.
kupiec_test <- function(hit, alpha) {
  n <- length(hit); x <- sum(hit); pi_hat <- x / n
  ll0 <- (n - x) * log(1 - alpha) + x * log(alpha)
  ll1 <- if (x == 0) n * log(1) else if (x == n) 0 else (n - x) * log(1 - pi_hat) + x * log(pi_hat)
  lr <- -2 * (ll0 - ll1)
  c(rate = pi_hat, LR = lr, p = 1 - pchisq(lr, 1))
}

# Christoffersen (1998) tests exceedance independence and conditional coverage.
christoffersen_test <- function(hit, alpha) {
  h0 <- hit[-length(hit)]; h1 <- hit[-1]
  n00 <- sum(h0 == 0 & h1 == 0); n01 <- sum(h0 == 0 & h1 == 1)
  n10 <- sum(h0 == 1 & h1 == 0); n11 <- sum(h0 == 1 & h1 == 1)
  p01 <- n01 / max(n00 + n01, 1); p11 <- n11 / max(n10 + n11, 1)
  p <- (n01 + n11) / (n00 + n01 + n10 + n11)
  safe <- function(k, q) if (k == 0) 0 else k * log(q)
  ll_ind0 <- safe(n00 + n10, 1 - p) + safe(n01 + n11, p)
  ll_ind1 <- safe(n00, 1 - p01) + safe(n01, p01) + safe(n10, 1 - p11) + safe(n11, p11)
  lr_ind <- -2 * (ll_ind0 - ll_ind1)
  lr_uc <- kupiec_test(hit, alpha)[["LR"]]
  c(LR_ind = lr_ind, p_ind = 1 - pchisq(lr_ind, 1),
    LR_cc = lr_uc + lr_ind, p_cc = 1 - pchisq(lr_uc + lr_ind, 2))
}

# McNeil and Frey (2000) test whether losses beyond VaR exceed ES on average.
# Scale residuals by ES and estimate the p-value by bootstrap.
es_exceedance_test <- function(loss, var, es, B = 5000, seed = 1) {
  hit <- loss > var
  if (sum(hit) < 3) return(c(n = sum(hit), mean = NA, p = NA))
  e <- (loss[hit] - es[hit]) / es[hit]
  set.seed(seed)
  e0 <- e - mean(e)
  boot <- replicate(B, mean(sample(e0, replace = TRUE)))
  c(n = sum(hit), mean = mean(e), p = mean(boot >= mean(e)))
}

# Acerbi and Szekely (2014) Z2. A strongly negative value suggests underestimated ES.
# For 97.5% ES, the indicative 5% threshold is -0.70.
acerbi_szekely_z2 <- function(loss, var, es, alpha) {
  hit <- loss > var
  sum(loss * hit / (length(loss) * alpha * es)) * -1 + 1
}

risk_backtest_table <- function(loss, var95, var99, es975, var975) {
  h95 <- as.numeric(loss > var95); h99 <- as.numeric(loss > var99)
  k95 <- kupiec_test(h95, 0.05); k99 <- kupiec_test(h99, 0.01)
  c95 <- christoffersen_test(h95, 0.05); c99 <- christoffersen_test(h99, 0.01)
  es <- es_exceedance_test(loss, var975, es975)
  z2 <- acerbi_szekely_z2(loss, var975, es975, 0.025)
  c(`VaR95 exceedances (%)` = 100 * k95[["rate"]], `p Kupiec 95` = k95[["p"]],
    `p Christoffersen cc 95` = c95[["p_cc"]],
    `VaR99 exceedances (%)` = 100 * k99[["rate"]], `p Kupiec 99` = k99[["p"]],
    `p Christoffersen cc 99` = c99[["p_cc"]],
    `ES97.5 : mean residual` = es[["mean"]], `ES97.5 : p McNeil-Frey` = es[["p"]],
    `ES97.5 : Z2 Acerbi-Szekely` = z2)
}
