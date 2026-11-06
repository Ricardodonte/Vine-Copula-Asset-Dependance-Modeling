source("scripts/00_setup.R")

R <- DATA$R
dates <- DATA$dates
write_table(describe_returns(R), "01_descriptive_statistics")

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
      Period = pn, Asset = a, Parameter = param_names(f$p, f$q),
      Estimate = round(unname(f$par), 6), `Standard error` = round(se, 6), check.names = FALSE)
  }

  for (a in colnames(U_ok)) {
    pit_check[[paste(pn, a)]] <- data.frame(
      Period = pn, Asset = a, nu = round(m$fits[[a]]$coef$nu, 2),
      `Lower tail, original PIT (%)` = round(100 * mean(U_bug[, a] < 0.05), 2),
      `Lower tail, corrected PIT (%)` = round(100 * mean(U_ok[, a] < 0.05), 2),
      `Upper tail, original PIT (%)` = round(100 * mean(U_bug[, a] > 0.95), 2),
      `Upper tail, corrected PIT (%)` = round(100 * mean(U_ok[, a] > 0.95), 2),
      `KS p-value, original PIT` = signif(suppressWarnings(ks.test(U_bug[, a], "punif")$p.value), 3),
      `KS p-value, corrected PIT` = signif(suppressWarnings(ks.test(U_ok[, a], "punif")$p.value), 3),
      check.names = FALSE)
  }

  set.seed(42)
  v_orig <- fit_vine(U_bug, original = TRUE)
  set.seed(42)
  v_corr <- fit_vine(U_ok)
  t_orig <- vine_table(v_orig); t_corr <- vine_table(v_corr)
  # Empirical Kendall's tau is rank-based and is preserved by the PIT.
  emp <- cor(U_ok, method = "kendall")
  t_orig$`empirical tau` <- t_corr$`empirical tau` <- NA
  for (tb in c("t_orig", "t_corr")) {
    x <- get(tb)
    for (r in which(x$Tree == 1)) {
      pr <- strsplit(x$Pair[r], "-")[[1]]
      x$`empirical tau`[r] <- round(emp[pr[1], pr[2]], 3)
    }
    attr(x, "fit") <- attr(get(tb), "fit")
    assign(tb, x)
  }
  vine_tabs[[pn]] <- list(orig = t_orig, corr = t_corr)
  fit_summary[[pn]] <- rbind(
    `Original (unscaled pt, no rotations)` = attr(t_orig, "fit"),
    `Corrected (standardized Student-t, rotations)` = attr(t_corr, "fit"))
  print(t_orig); print(t_corr)
  saveRDS(list(margins = m, v_orig = v_orig, v_corr = v_corr),
          file.path(OUT_CACHE, paste0("period_", which(names(PERIODS) == pn), ".rds")))
}

write_table(do.call(rbind, egarch_tables), "02_egarch_parameters", digits = 6)
write_table(do.call(rbind, pit_check), "03_pit_diagnostics", digits = 2)

lines <- c()
for (pn in names(PERIODS)) {
  lines <- c(lines, paste("##", pn), "", "**Original version (PIT pt(z, nu), no rotated copulas)**", "",
             knitr_like(vine_tabs[[pn]]$orig, 3), "",
             sprintf("logLik = %.1f, AIC = %.1f, BIC = %.1f", attr(vine_tabs[[pn]]$orig, "fit")[1],
                     attr(vine_tabs[[pn]]$orig, "fit")[2], attr(vine_tabs[[pn]]$orig, "fit")[3]), "",
             "**Corrected version (standardized Student-t PIT, rotations allowed)**", "",
             knitr_like(vine_tabs[[pn]]$corr, 3), "",
             sprintf("logLik = %.1f, AIC = %.1f, BIC = %.1f", attr(vine_tabs[[pn]]$corr, "fit")[1],
                     attr(vine_tabs[[pn]]$corr, "fit")[2], attr(vine_tabs[[pn]]$corr, "fit")[3]), "")
}
writeLines(lines, file.path(OUT_TAB, "04_vines_by_period.md"))
saveRDS(vine_tabs, file.path(OUT_CACHE, "vine_tabs.rds"))
