# Données locales

Placer les deux fichiers dans ce dossier avant de lancer l’analyse.

| Fichier | Contenu |
| --- | --- |
| `prix_spy_btc_bnd.csv` | Prix ajustés de SPY, BTC-USD et BND. La première colonne contient les dates au format `AAAA-MM-JJ`. Les colonnes de prix sont `SPY`, `BTC` et `BND`. |
| `XAU_USD.csv` | Export TradingView de l’or. Le lecteur utilise les deux premières colonnes, avec les dates au format `MM/JJ/AAAA` et les prix de clôture. |

Pour régénérer le premier fichier depuis Yahoo Finance, installer `quantmod` puis lancer `Rscript scripts/export_prix.R` depuis la racine. Le fichier d’origine de l’or est à placer manuellement. Utiliser les mêmes exports de prix que le rapport pour comparer les chiffres.

L’analyse aligne les prix sur les dates communes et calcule les log-rendements sur la période 2016–2025. Les deux CSV restent locaux.
