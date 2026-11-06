suppressPackageStartupMessages({
  library(VineCopula)
  library(ggplot2)
})

for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f, encoding = "UTF-8")

PRICES_CSV <- "data/prix_spy_btc_bnd.csv"   # exporté avec export_prix.R (quantmod / Yahoo)
GOLD_CSV   <- "data/XAU_USD.csv"            # TradingView, comme dans le projet d'origine
BACKEND <- egarch_backend()
OUT_DIR <- file.path("results", BACKEND)
OUT_TAB <- file.path(OUT_DIR, "tables")
OUT_FIG <- file.path(OUT_DIR, "figures")
OUT_CACHE <- file.path(OUT_DIR, "cache")
for (d in c(OUT_TAB, OUT_FIG, OUT_CACHE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

ORDERS <- rbind(SPY = c(1, 1), BTC = c(0, 0), BND = c(0, 0), GOLD = c(0, 0))

PERIODS <- list(
  `Pré-COVID (2016-01 → 2019-11)`  = c("2016-01-04", "2019-11-01"),
  `COVID (2020-01 → 2021-12)`      = c("2020-01-01", "2022-01-01"),
  `Post-COVID (2022-01 → 2025-11)` = c("2022-01-01", "2025-11-03"))

DATA <- load_returns(PRICES_CSV, GOLD_CSV)
cat(sprintf("Données : %d rendements quotidiens du %s au %s\n",
            nrow(DATA$R), min(DATA$dates), max(DATA$dates)))

write_table <- function(df, name, digits = 4) {
  write.csv(df, file.path(OUT_TAB, paste0(name, ".csv")), row.names = TRUE)
  md <- knitr_like(df, digits)
  writeLines(md, file.path(OUT_TAB, paste0(name, ".md")))
}

knitr_like <- function(df, digits = 4) {
  df <- as.data.frame(df, check.names = FALSE)
  rn <- rownames(df)
  has_rn <- !is.null(rn) && !all(rn == as.character(seq_len(nrow(df))))
  fmt <- function(x) if (is.numeric(x)) formatC(x, format = "f", digits = digits) else as.character(x)
  cells <- as.data.frame(lapply(df, fmt), check.names = FALSE)
  if (has_rn) cells <- cbind(` ` = rn, cells)
  hdr <- paste0("| ", paste(names(cells), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(cells)), collapse = "|"), "|")
  body <- apply(cells, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(hdr, sep, body)
}
