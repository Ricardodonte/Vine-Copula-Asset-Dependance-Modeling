# On ajuste les marges une fois sur tout l’échantillon. Sur 252 jours, l’EGARCH-t est trop instable.
source("scripts/00_setup.R")
library(parallel)

R <- DATA$R; dates <- DATA$dates
m <- fit_margins(R, ORDERS, pit = "correct")
U_ok <- m$U
U_bug <- sapply(m$fits, function(f) pt(f$z, df = f$coef$nu))

PAIRS <- list(c("SPY", "BTC"), c("SPY", "BND"), c("SPY", "GOLD"),
              c("BTC", "GOLD"), c("BND", "GOLD"), c("BTC", "BND"))
WIN <- 252; STEP <- 5
ends <- seq(WIN, nrow(R), by = STEP)

roll <- mclapply(ends, function(e) {
  w <- (e - WIN + 1):e
  do.call(rbind, lapply(PAIRS, function(pr) {
    a <- pair_dependence(U_ok[w, pr[1]], U_ok[w, pr[2]], rotations = TRUE)
    b <- pair_dependence(U_bug[w, pr[1]], U_bug[w, pr[2]], rotations = FALSE)
    data.frame(date = dates[e], paire = paste(pr, collapse = "-"),
               version = c("Corrigée", "Origine"),
               tau = c(a[["tau"]], b[["tau"]]), LTD = c(a[["LTD"]], b[["LTD"]]),
               UTD = c(a[["UTD"]], b[["UTD"]]), famille = c(a[["family"]], b[["family"]]))
  }))
}, mc.cores = if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1))
roll <- do.call(rbind, roll)
saveRDS(roll, file.path(OUT_CACHE, "rolling.rds"))
write.csv(roll, file.path(OUT_TAB, "05_dependance_glissante.csv"), row.names = FALSE)

roll$regime <- cut(roll$date, as.Date(c("2016-01-01", "2020-01-01", "2022-01-01", "2026-01-01")),
                   labels = c("Pré-COVID", "COVID", "Post-COVID"), right = FALSE)
agg <- aggregate(cbind(tau, LTD, UTD) ~ paire + version + regime, data = roll, FUN = mean)
agg <- agg[order(agg$paire, agg$regime, agg$version), ]
rownames(agg) <- NULL
write_table(agg, "06_dependance_glissante_moyennes", digits = 3)

long <- rbind(
  data.frame(roll[, c("date", "paire", "version")], mesure = "tau de Kendall", valeur = roll$tau),
  data.frame(roll[, c("date", "paire", "version")], mesure = "LTD (co-crash)", valeur = roll$LTD),
  data.frame(roll[, c("date", "paire", "version")], mesure = "UTD (co-hausse)", valeur = roll$UTD))
g <- ggplot(subset(long, version == "Corrigée"), aes(date, valeur, colour = paire)) +
  geom_line(linewidth = 0.4) +
  geom_vline(xintercept = as.Date(c("2020-01-01", "2022-01-01")), linetype = 2, colour = "grey50") +
  facet_wrap(~mesure, ncol = 1, scales = "free_y") +
  labs(title = "Dépendance glissante sur 252 jours (PIT corrigée, copules tournées)",
       x = NULL, y = NULL, colour = "Paire") + theme_minimal(base_size = 11)
ggsave(file.path(OUT_FIG, "01_dependance_glissante.png"), g, width = 10, height = 9, dpi = 150)

g2 <- ggplot(subset(long, paire %in% c("SPY-BTC", "SPY-BND", "BND-GOLD") & mesure != "tau de Kendall"),
             aes(date, valeur, colour = version)) +
  geom_line(linewidth = 0.4) + facet_grid(mesure ~ paire) +
  labs(title = "Effet de la correction de la PIT sur les dépendances de queue",
       x = NULL, y = NULL, colour = NULL) + theme_minimal(base_size = 11)
ggsave(file.path(OUT_FIG, "02_queues_origine_vs_corrigee.png"), g2, width = 11, height = 6, dpi = 150)
