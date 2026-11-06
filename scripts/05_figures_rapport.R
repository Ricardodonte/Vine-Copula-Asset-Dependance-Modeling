source("scripts/00_setup.R")

roll <- readRDS(file.path(OUT_CACHE, "rolling.rds"))
roll <- roll[roll$version == "Corrigée", ]
lab <- c(`BND-GOLD` = "GOLD-BND", `SPY-BND` = "SPY-BND", `BTC-GOLD` = "BTC-GOLD", `SPY-BTC` = "SPY-BTC")

periods <- list(
  pre  = list(range = c("2017-01-01", "2019-12-31"), pairs = c("BND-GOLD", "SPY-BND", "BTC-GOLD"),
              title = "pré-COVID"),
  cov  = list(range = c("2020-01-01", "2021-12-31"), pairs = c("SPY-BTC", "BTC-GOLD", "BND-GOLD"),
              title = "pendant la crise COVID-19"),
  post = list(range = c("2022-01-01", "2025-11-30"), pairs = c("SPY-BTC", "BND-GOLD", "SPY-BND"),
              title = "post-COVID-19"))
measures <- c(tau = "Évolution de la dépendance (tau de Kendall)", LTD = "Évolution de la LTD",
              UTD = "Évolution de l'UTD")
cols <- c(`GOLD-BND` = "#1b9e3e", `SPY-BND` = "#3b6fd8", `BTC-GOLD` = "#e8524a", `SPY-BTC` = "#9b4dca")

for (p in names(periods)) {
  P <- periods[[p]]
  x <- roll[roll$date >= as.Date(P$range[1]) & roll$date <= as.Date(P$range[2]) & roll$paire %in% P$pairs, ]
  x$Série <- lab[x$paire]
  for (m in names(measures)) {
    g <- ggplot(x, aes(date, .data[[m]], colour = Série)) + geom_line(linewidth = 0.55) +
      scale_colour_manual(values = cols) +
      labs(title = paste(measures[[m]], P$title), x = "Date", y = "Valeur") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(OUT_FIG, sprintf("rapport_%s_%s.png", p, m)), g, width = 8, height = 4.2, dpi = 170)
  }
}
