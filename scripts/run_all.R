args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) == 1L && !args[1] %in% c("rugarch", "manual"))) {
  stop("Usage : Rscript scripts/run_all.R [rugarch|manual]", call. = FALSE)
}
if (length(args) == 1L) options(vine.use_rugarch = args[1] == "rugarch")

for (s in c("01_periods.R", "02_rolling.R", "03_backtest.R", "04_risk_backtest.R", "05_report_figures.R")) {
  cat("\n########", s, "########\n")
  local(source(file.path("scripts", s), encoding = "UTF-8"))
}
