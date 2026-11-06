source("scripts/00_setup.R")
library(parallel)

R <- DATA$R
dates <- DATA$dates
OOS_START <- as.Date("2019-11-04")
REBAL <- 5; W_MARG <- 1000; W_COP <- 504; W_HIST <- 252
N_SIM <- 10000; N_SIM_SMALL <- 500; SEEDS <- 1:5
COST_BPS <- 10
d <- ncol(R); ew <- rep(1 / d, d)

oos <- which(dates >= OOS_START)
rebal_idx <- oos[seq(1, length(oos), by = REBAL)]

simulate_gauss <- function(n, Rho) {
  Z <- matrix(rnorm(n * ncol(Rho)), n) %*% chol(Rho)
  pnorm(Z)
}

risk_from_sims <- function(port) {
  loss <- -port
  v95 <- unname(quantile(loss, 0.95)); v99 <- unname(quantile(loss, 0.99))
  v975 <- unname(quantile(loss, 0.975))
  c(var95 = v95, var99 = v99, var975 = v975, es975 = mean(loss[loss >= v975]))
}

egarch_update_series <- function(fits, t0, t) {
  lapply(names(fits), function(a) {
    f <- fits[[a]]
    if (t > t0) egarch_update(f, R[t0:(t - 1), a]) else list(mean_next = f$mean_next, sigma_next = f$sigma_next)
  })
}

one_week <- function(k) {
  t0 <- rebal_idx[k]
  t_end <- if (k < length(rebal_idx)) rebal_idx[k + 1] - 1 else nrow(R)
  win <- max(1, t0 - W_MARG):(t0 - 1)
  m <- fit_margins(R[win, ], ORDERS, pit = "correct")
  nw <- length(win); cw <- max(1, nw - W_COP + 1):nw
  U_ok <- m$U[cw, ]

  set.seed(1000 + k)
  v_corr <- fit_vine(U_ok)
  Rho <- cor(qnorm(U_ok))
  Rho <- (Rho + t(Rho)) / 2

  H <- R[max(1, t0 - W_HIST):(t0 - 1), ]
  w <- list(
    `Equal weight (1/N)` = ew,
    `Historical MinVar` = min_var_weights(cov(H)),
    `Historical MinCVaR` = min_cvar_weights(H))
  for (s in SEEDS) {
    set.seed(s * 100000 + k)
    Uc <- RVineSim(N_SIM, v_corr)
    Sc <- scenarios_from_uniforms(Uc, m$fits, "correct")
    w[[paste0("MinCVaR R-vine ARMA-EGARCH #", s)]] <- min_cvar_weights(Sc)
    if (s == 1) w[["MinVar R-vine ARMA-EGARCH"]] <- min_var_weights(cov(Sc))
    set.seed(s * 200000 + k)
    Ug <- simulate_gauss(N_SIM, Rho)
    w[[paste0("Gaussian copula MinCVaR #", s)]] <- min_cvar_weights(scenarios_from_uniforms(Ug, m$fits, "correct"))
    # Use the first 500 scenarios to measure Monte Carlo variation.
    w[[paste0("MinCVaR R-vine ARMA-EGARCH 500 sims #", s)]] <- min_cvar_weights(Sc[seq_len(N_SIM_SMALL), ])
  }

  risk <- list()
  for (t in t0:t_end) {
    upd <- egarch_update_series(m$fits, t0, t)
    fits_t <- m$fits
    for (j in seq_along(fits_t)) {
      fits_t[[j]]$mean_next <- upd[[j]]$mean_next; fits_t[[j]]$sigma_next <- upd[[j]]$sigma_next
    }
    set.seed(7 * t)
    pc <- scenarios_from_uniforms(RVineSim(N_SIM, v_corr), fits_t, "correct") %*% ew
    pg <- scenarios_from_uniforms(simulate_gauss(N_SIM, Rho), fits_t, "correct") %*% ew
    pi <- scenarios_from_uniforms(matrix(runif(N_SIM * d), N_SIM), fits_t, "correct") %*% ew
    hs <- R[max(1, t - W_HIST):(t - 1), ] %*% ew
    mu <- colMeans(R[max(1, t - W_HIST):(t - 1), ]); S <- cov(R[max(1, t - W_HIST):(t - 1), ])
    m_p <- sum(ew * mu); s_p <- sqrt(drop(t(ew) %*% S %*% ew))
    gaus_vc <- c(var95 = -m_p + s_p * qnorm(0.95), var99 = -m_p + s_p * qnorm(0.99),
                 var975 = -m_p + s_p * qnorm(0.975),
                 es975 = -m_p + s_p * dnorm(qnorm(0.975)) / 0.025)
    risk[[length(risk) + 1]] <- rbind(
      data.frame(t = t, model = "R-vine ARMA-EGARCH", t(risk_from_sims(pc))),
      data.frame(t = t, model = "EGARCH + Gaussian copula", t(risk_from_sims(pg))),
      data.frame(t = t, model = "EGARCH + independence", t(risk_from_sims(pi))),
      data.frame(t = t, model = "Historical simulation (252 days)", t(risk_from_sims(hs))),
      data.frame(t = t, model = "Normal variance-covariance (252 days)", t(gaus_vc)))
  }
  fam <- v_corr$family[lower.tri(v_corr$family)]
  list(k = k, t0 = t0, weights = w, risk = do.call(rbind, risk), families = family_label(fam),
       n_rot = sum(fam > 10))
}

cache_name <- if (BACKEND == "rugarch") {
  loadNamespace("rugarch")
  paste0("backtest_weeks_rvine_", utils::packageVersion("rugarch"), ".rds")
} else {
  "backtest_weeks_rvine.rds"
}
cache_file <- file.path(OUT_CACHE, cache_name)
use_cache <- file.exists(cache_file)
if (use_cache) {
  weeks <- readRDS(cache_file)
} else {
  t_start <- Sys.time()
  weeks <- mclapply(seq_along(rebal_idx), function(k) {
    tryCatch(one_week(k), error = function(e) list(k = k, error = conditionMessage(e)))
  }, mc.cores = if (.Platform$OS.type == "windows") 1 else max(1, detectCores()), mc.preschedule = FALSE)
  cat("Elapsed:", format(Sys.time() - t_start), "\n")
}
errs <- Filter(function(x) !is.null(x$error), weeks)
if (length(errs)) { cat(length(errs), "weekly fits failed\n"); print(head(sapply(errs, `[[`, "error"))) }
stopifnot(length(errs) == 0)
if (!use_cache) saveRDS(weeks, cache_file)

simple <- exp(R) - 1
strategies <- names(weeks[[1]]$weights)
run_strategy <- function(name) {
  ret <- numeric(0); turn <- 0; w_cur <- rep(0, d); tt <- integer(0)
  for (k in seq_along(weeks)) {
    t0 <- weeks[[k]]$t0
    t_end <- if (k < length(weeks)) weeks[[k + 1]]$t0 - 1 else nrow(R)
    w_new <- weeks[[k]]$weights[[name]]
    tr <- sum(abs(w_new - w_cur)); turn <- turn + tr
    cost <- tr * COST_BPS / 1e4
    w <- w_new
    for (t in t0:t_end) {
      rp <- sum(w * simple[t, ])
      if (t == t0) rp <- rp - cost
      ret <- c(ret, rp); tt <- c(tt, t)
      w <- w * (1 + simple[t, ]); w <- w / sum(w)
    }
    w_cur <- w
  }
  years <- length(ret) / 252
  # Annual turnover excludes the initial purchase.
  list(ret = ret, t = tt, turnover = (turn - 1) / years)
}
paths <- lapply(setNames(strategies, strategies), run_strategy)
saveRDS(paths, file.path(OUT_CACHE, "backtest_paths.rds"))
oos_dates <- dates[paths[[1]]$t]

perf <- t(sapply(paths, function(p) perf_metrics(p$ret, turnover = p$turnover)))

group_of <- function(n) sub(" #[0-9]+$", "", n)
grp <- group_of(rownames(perf))
main_tab <- do.call(rbind, lapply(unique(grp), function(g) {
  rows <- perf[grp == g, , drop = FALSE]
  colMeans(rows)
}))
rownames(main_tab) <- unique(grp)
write_table(round(main_tab, 3), "07_strategy_performance", digits = 3)

mc_tab <- do.call(rbind, lapply(unique(grp[grepl("#", rownames(perf))]), function(g) {
  rows <- perf[grp == g, , drop = FALSE]
  data.frame(Strategy = g,
             `Sharpe min` = min(rows[, "Sharpe"]), `Sharpe max` = max(rows[, "Sharpe"]),
             `Min total return (%)` = min(rows[, "Total return (%)"]),
             `Max total return (%)` = max(rows[, "Total return (%)"]),
             `Max DD min (%)` = min(rows[, "Max drawdown (%)"]), `Max DD max (%)` = max(rows[, "Max drawdown (%)"]),
             check.names = FALSE)
}))
write_table(mc_tab, "08_monte_carlo_variation", digits = 3)

# Bootstrap comparisons use the first random seed.
ref <- "MinCVaR R-vine ARMA-EGARCH #1"
comparisons <- c("Equal weight (1/N)", "Historical MinVar", "Historical MinCVaR", "MinVar R-vine ARMA-EGARCH",
                 "Gaussian copula MinCVaR #1")
tests <- do.call(rbind, lapply(comparisons, function(cn) {
  bt <- paired_block_bootstrap(paths[[ref]]$ret, paths[[cn]]$ret)
  data.frame(Comparison = paste(sub(" #1", "", ref), "vs", sub(" #1", "", cn)),
             `Sharpe difference` = bt$diff[1], `p Sharpe` = bt$p_value[1],
             `Daily 95% CVaR difference` = 100 * bt$diff[2], `p CVaR` = bt$p_value[2],
             `Max drawdown difference (pp)` = 100 * bt$diff[3], `p Max DD` = bt$p_value[3], check.names = FALSE)
}))
write_table(tests, "09_significance_tests", digits = 3)

show <- c("Equal weight (1/N)", "Historical MinVar", "Historical MinCVaR", "MinVar R-vine ARMA-EGARCH",
          "MinCVaR R-vine ARMA-EGARCH #1", "Gaussian copula MinCVaR #1")
cum <- do.call(rbind, lapply(show, function(s) data.frame(date = oos_dates, strategy = sub(" #1", "", s),
                                                         wealth = cumprod(1 + paths[[s]]$ret))))
g <- ggplot(cum, aes(date, wealth, colour = strategy)) + geom_line(linewidth = 0.5) +
  labs(title = "Out-of-sample cumulative wealth (10 bps trading cost)", x = NULL, y = "Value of 1 invested",
       colour = NULL) + theme_minimal(base_size = 11) + theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))
ggsave(file.path(OUT_FIG, "03_cumulative_wealth.png"), g, width = 10, height = 6, dpi = 150)

base <- cumprod(1 + paths[["Historical MinVar"]]$ret)
rel <- do.call(rbind, lapply(setdiff(show, c("Historical MinVar", "Equal weight (1/N)")), function(s)
  data.frame(date = oos_dates, strategy = sub(" #1", "", s),
             difference = 100 * (cumprod(1 + paths[[s]]$ret) / base - 1))))
g <- ggplot(rel, aes(date, difference, colour = strategy)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.5) +
  labs(title = "Performance relative to historical MinVar (%), excluding 1/N", x = NULL, y = "%", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 2))
ggsave(file.path(OUT_FIG, "04_relative_performance.png"), g, width = 10, height = 6, dpi = 150)

wdf <- do.call(rbind, lapply(c("MinCVaR R-vine ARMA-EGARCH #1", "Historical MinCVaR"),
  function(s) do.call(rbind, lapply(weeks, function(wk) data.frame(
    date = dates[wk$t0], strategy = sub(" #1", "", s), asset = colnames(R), weight = wk$weights[[s]])))))
g <- ggplot(wdf, aes(date, weight, fill = asset)) + geom_area() + facet_wrap(~strategy, ncol = 1) +
  labs(title = "Weekly portfolio weights", x = NULL, y = NULL, fill = NULL) + theme_minimal(base_size = 11)
ggsave(file.path(OUT_FIG, "05_weights.png"), g, width = 10, height = 8, dpi = 150)

avg_w <- do.call(rbind, lapply(c("Historical MinVar", "Historical MinCVaR", "MinVar R-vine ARMA-EGARCH",
                                 "MinCVaR R-vine ARMA-EGARCH #1", "Gaussian copula MinCVaR #1"), function(s)
  colMeans(do.call(rbind, lapply(weeks, function(wk) wk$weights[[s]])))))
rownames(avg_w) <- sub(" #1", "", c("Historical MinVar", "Historical MinCVaR", "MinVar R-vine ARMA-EGARCH",
                                    "MinCVaR R-vine ARMA-EGARCH #1", "Gaussian copula MinCVaR #1"))
colnames(avg_w) <- colnames(R)
write_table(round(100 * avg_w, 1), "10_average_weights_pct", digits = 1)

fam_share <- table(unlist(lapply(weeks, `[[`, "families")))
write_table(as.data.frame(sort(fam_share, decreasing = TRUE)), "11_backtest_copula_families", digits = 0)
cat("Share of weeks with at least one rotated copula:",
    mean(sapply(weeks, `[[`, "n_rot") > 0), "\n")

saveRDS(do.call(rbind, lapply(weeks, `[[`, "risk")), file.path(OUT_CACHE, "risk_forecasts.rds"))
