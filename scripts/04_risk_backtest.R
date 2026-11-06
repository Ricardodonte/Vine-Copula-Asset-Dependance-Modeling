source("scripts/00_setup.R")

R <- DATA$R
dates <- DATA$dates
rf <- readRDS(file.path(OUT_CACHE, "risk_forecasts.rds"))
ew <- rep(1 / ncol(R), ncol(R))
rf$loss <- -drop(R[rf$t, ] %*% ew)
rf$date <- dates[rf$t]

tab <- do.call(rbind, lapply(split(rf, rf$model), function(x) {
  x <- x[order(x$t), ]
  risk_backtest_table(x$loss, x$var95, x$var99, x$es975, x$var975)
}))
write_table(round(tab, 3), "12_var_es_backtest", digits = 3)

rf$regime <- cut(rf$date, as.Date(c("2019-11-01", "2022-01-01", "2026-01-01")),
                 labels = c("COVID (nov. 2019 → 2021)", "Post-COVID (2022 → 2025)"), right = FALSE)
tab_reg <- do.call(rbind, lapply(split(rf, list(rf$model, rf$regime)), function(x) {
  x <- x[order(x$t), ]
  data.frame(Model = x$model[1], Regime = x$regime[1],
             t(risk_backtest_table(x$loss, x$var95, x$var99, x$es975, x$var975)), check.names = FALSE)
}))
rownames(tab_reg) <- NULL
write_table(tab_reg, "13_var_es_backtest_by_regime", digits = 3)

sub_rf <- subset(rf, model %in% c("R-vine ARMA-EGARCH", "EGARCH + Gaussian copula",
                                   "Historical simulation (252 days)"))
g <- ggplot(sub_rf, aes(date)) +
  geom_point(data = subset(rf, model == "R-vine ARMA-EGARCH"), aes(y = 100 * loss),
             size = 0.3, colour = "grey60") +
  geom_line(aes(y = 100 * var99, colour = model), linewidth = 0.45) +
  labs(title = "Daily losses and 99% VaR forecasts for the 1/N portfolio", x = NULL, y = "Loss (%)",
       colour = NULL) + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT_FIG, "06_var99.png"), g, width = 10, height = 5, dpi = 150)
