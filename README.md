# Vine copulas et optimisation de portefeuille

Étude des dépendances extrêmes entre le SPY, le Bitcoin, le BND et l’or. Le projet estime des marges ARMA-EGARCH et une R-vine, puis compare les stratégies d’allocation et les prévisions de VaR et d’Expected Shortfall hors échantillon.

[Lire le rapport](report/rapport.pdf)

## Structure

| Dossier | Contenu |
| --- | --- |
| `R/` | Fonctions de modélisation, d’optimisation et de backtest de risque |
| `scripts/` | Lancement de l’analyse et scripts par étape |
| `data/` | Instructions pour préparer les deux fichiers de prix |
| `report/` | Sources LaTeX, figures intégrées et PDF |
| `results/` | Tableaux, figures et caches produits localement |

## Relancer l’analyse

Installer les packages dans R, puis préparer les fichiers décrits dans [data/README.md](data/README.md).

```r
install.packages(c("rugarch", "VineCopula", "quadprog", "Rglpk", "slam",
                   "ggplot2", "moments", "numDeriv"))
```

Depuis la racine du dépôt :

```sh
Rscript scripts/run_all.R
Rscript scripts/run_all.R manual
```

La première commande utilise `rugarch`. La seconde utilise l’implémentation manuelle des marges EGARCH et de la Student standardisée. Les deux modes écrivent séparément dans `results/rugarch/` et `results/manual/`. Le backtest complet peut prendre du temps, en particulier sous Windows.

## Recompiler le rapport

Depuis le dossier `report/` :

```sh
latexmk -pdf -jobname=rapport main.tex
```

Les figures nécessaires à la compilation sont dans `report/figures/`. Les fichiers intermédiaires de LaTeX, les données brutes et les résultats générés sont ignorés par Git.
