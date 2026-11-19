read_gold <- function(path) {
  g <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")[, 1:2]
  names(g) <- c("date", "GOLD")
  g$date <- as.Date(g$date, format = "%m/%d/%Y")
  g$GOLD <- as.numeric(gsub(",", "", g$GOLD))
  g[order(g$date), ]
}

read_prices <- function(path) {
  p <- read.csv(path, check.names = FALSE)
  names(p)[1] <- "date"
  p$date <- as.Date(p$date)
  names(p) <- sub("^BTC.*", "BTC", names(p))
  p[, c("date", "SPY", "BTC", "BND")]
}

load_returns <- function(prices_path, gold_path, start = "2016-01-01", end = "2025-11-03") {
  px <- merge(read_prices(prices_path), read_gold(gold_path), by = "date")
  px <- na.omit(px)
  px <- px[px$date >= as.Date(start) & px$date <= as.Date(end), ]
  px <- px[order(px$date), ]
  lr <- apply(log(as.matrix(px[, -1])), 2, diff)
  list(dates = px$date[-1], R = lr, prices = px)
}

describe_returns <- function(R) {
  f <- function(x) c(Mean = mean(x), `Std. dev.` = sd(x), Skewness = moments::skewness(x),
                     Kurtosis = moments::kurtosis(x), Min = min(x), Q1 = unname(quantile(x, 0.25)),
                     Median = median(x), Q3 = unname(quantile(x, 0.75)), Max = max(x))
  round(sapply(as.data.frame(R), f), 4)
}
