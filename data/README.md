# Local data

Place both input files in this directory before running the analysis.

| File | Contents |
| --- | --- |
| `prix_spy_btc_bnd.csv` | Adjusted prices for SPY, BTC-USD, and BND. The first column contains dates in `YYYY-MM-DD` format. The price columns are `SPY`, `BTC`, and `BND`. |
| `XAU_USD.csv` | TradingView gold export. The reader uses the first two columns for dates in `MM/DD/YYYY` format and closing prices. |

To regenerate the first file from Yahoo Finance, install `quantmod` and run `Rscript scripts/export_prices.R` from the repository root. Place the gold export in this directory manually. Use the same price exports as the study to compare the numbers.

The analysis aligns prices on their common dates and calculates log returns for 2016–2025. Both CSV files stay local.
