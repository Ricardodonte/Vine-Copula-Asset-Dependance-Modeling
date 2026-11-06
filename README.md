# Vine Copula Asset Dependence Modeling

This study models tail dependence among SPY, Bitcoin, BND, and gold with ARMA-EGARCH margins and an R-vine copula. The out-of-sample backtests compare portfolio allocation strategies and evaluate Value at Risk (VaR) and Expected Shortfall (ES) forecasts.

[Read the report](report/report.pdf)

## Repository structure

| Directory | Contents |
| --- | --- |
| `R/` | Models, portfolio optimization, and risk tests |
| `scripts/` | Analysis entry point and individual stages |
| `data/` | Instructions for preparing the price files |
| `report/` | LaTeX sources, figures, and the PDF |
| `results/` | Tables, figures, and caches generated locally |

## Run the analysis

Install the packages in R and prepare the inputs described in [data/README.md](data/README.md).

```r
install.packages(c("rugarch", "VineCopula", "quadprog", "Rglpk", "slam",
                   "ggplot2", "moments", "numDeriv"))
```

From the repository root, run:

```sh
Rscript scripts/run_all.R
Rscript scripts/run_all.R manual
```

The first command uses `rugarch`. The second uses the manual implementation of the EGARCH margins and standardized Student-t distribution. The two modes write to `results/rugarch/` and `results/manual/`, respectively. The full backtest can take a while, particularly on Windows.

## Build the report

From the `report/` directory, run:

```sh
latexmk -pdf -jobname=report main.tex
```

The figures required to compile the report are in `report/figures/`. Git ignores raw data, generated results, and LaTeX build files.
