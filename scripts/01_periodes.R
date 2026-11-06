source("scripts/00_setup.R")

R <- DATA$R
dates <- DATA$dates
write_table(describe_returns(R), "01_statistiques_descriptives")

sel <- function(p) which(dates >= as.Date(p[1]) & dates <= as.Date(p[2]))

egarch_tables <- list(); vine_tabs <- list(); pit_check <- list(); fit_summary <- list()

for (pn in names(PERIODS)) {
  idx <- sel(PERIODS[[pn]])
  cat("\n====", pn, ":", length(idx), "observations ====\n")
  m <- fit_margins(R[idx, ], ORDERS, pit = "correct")
  U_ok <- m$U
  U_bug <- sapply(m$fits, function(f) pt(f$z, df = f$coef$nu))

  for (a in names(m$fits)) {
    f <- m$fits[[a]]
    se <- egarch_se(f)
    egarch_tables[[paste(pn, a)]] <- data.frame(
      Période = pn, Actif = a, Paramètre = param_names(f$p, f$q),
      Estimation = round(unname(f$par), 6), `Écart-type` = round(se, 6), check.names = FALSE)
  }

  for (a in colnames(U_ok)) {
    pit_check[[paste(pn, a)]] <- data.frame(
      Période = pn, Actif = a, nu = round(m$fits[[a]]$coef$nu, 2),
      `Queue basse, PIT d'origine (%)` = round(100 * mean(U_bug[, a] < 0.05), 2),
      `Queue basse, PIT corrigée (%)` = round(100 * mean(U_ok[, a] < 0.05), 2),
      `Queue haute, PIT d'origine (%)` = round(100 * mean(U_bug[, a] > 0.95), 2),
      `Queue haute, PIT corrigée (%)` = round(100 * mean(U_ok[, a] > 0.95), 2),
      `p KS, PIT d'origine` = signif(suppressWarnings(ks.test(U_bug[, a], "punif")$p.value), 3),
      `p KS, PIT corrigée` = signif(suppressWarnings(ks.test(U_ok[, a], "punif")$p.value), 3),
      check.names = FALSE)
  }

  set.seed(42)
  v_orig <- fit_vine(U_bug, original = TRUE)
  set.seed(42)
  v_corr <- fit_vine(U_ok)
  t_orig <- vine_table(v_orig); t_corr <- vine_table(v_corr)
  # Le tau empirique repose sur les rangs et reste identique après la PIT.
  emp <- cor(U_ok, method = "kendall")
  t_orig$`tau empirique` <- t_corr$`tau empirique` <- NA
  for (tb in c("t_orig", "t_corr")) {
    x <- get(tb)
    for (r in which(x$Arbre == 1)) {
      pr <- strsplit(x$Paire[r], "-")[[1]]
      x$`tau empirique`[r] <- round(emp[pr[1], pr[2]], 3)
    }
    attr(x, "fit") <- attr(get(tb), "fit")
    assign(tb, x)
  }
  vine_tabs[[pn]] <- list(orig = t_orig, corr = t_corr)
  fit_summary[[pn]] <- rbind(
    `Origine (pt non réduite, sans rotations)` = attr(t_orig, "fit"),
    `Corrigée (loi std, avec rotations)` = attr(t_corr, "fit"))
  print(t_orig); print(t_corr)
  saveRDS(list(margins = m, v_orig = v_orig, v_corr = v_corr),
          file.path(OUT_CACHE, paste0("periode_", which(names(PERIODS) == pn), ".rds")))
}

write_table(do.call(rbind, egarch_tables), "02_parametres_egarch", digits = 6)
write_table(do.call(rbind, pit_check), "03_diagnostic_pit", digits = 2)

lines <- c()
for (pn in names(PERIODS)) {
  lines <- c(lines, paste("##", pn), "", "**Version d'origine (PIT pt(z, nu), sans copules tournées)**", "",
             knitr_like(vine_tabs[[pn]]$orig, 3), "",
             sprintf("logLik = %.1f, AIC = %.1f, BIC = %.1f", attr(vine_tabs[[pn]]$orig, "fit")[1],
                     attr(vine_tabs[[pn]]$orig, "fit")[2], attr(vine_tabs[[pn]]$orig, "fit")[3]), "",
             "**Version corrigée (PIT loi std, copules tournées autorisées)**", "",
             knitr_like(vine_tabs[[pn]]$corr, 3), "",
             sprintf("logLik = %.1f, AIC = %.1f, BIC = %.1f", attr(vine_tabs[[pn]]$corr, "fit")[1],
                     attr(vine_tabs[[pn]]$corr, "fit")[2], attr(vine_tabs[[pn]]$corr, "fit")[3]), "")
}
writeLines(lines, file.path(OUT_TAB, "04_vines_par_periode.md"))
saveRDS(vine_tabs, file.path(OUT_CACHE, "vine_tabs.rds"))
