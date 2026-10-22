## Tableau principal (budget 1000 lectures, 10 fonctions × 100 départs)

| méthode | bits / lecture | lectures / pas | succès (IC 95 %) | lectures médianes (succès) |
|---|---|---|---|---|
| sonde-complexe (dépôt d'origine) (8 bits) | 24 | 2 | 98 % (97–99) | 27 |
| sonde-complexe (dépôt d'origine) (12 bits) | 36 | 2 | 100 % (99–100) | 29 |
| sonde-couleur-sans-axe-signée (16 bits) | 16 | 2 | 99 % (99–100) | 30 |
| sonde-couleur-sans-axe-signée (24 bits) | 24 | 2 | 100 % (99–100) | 17 |
| sonde-couleur-sans-axe (24 bits) | 24 | 2 | 59 % (56–62) | 48 |
| sonde-couleur (4D) (24 bits) | 24 | 5 | 76 % (73–78) | 73 |
| Newton-complexe-sans-couleur | exacte | 2 (f + f') | 80 % (78–82) | 11 |
| Newton-sans-couleur | exacte | 5 | 92 % (91–94) | 26 |
| Newton-jacobienne-exacte | exacte | 2 (F + J) | 92 % (91–94) | 11 |
| Broyden-sans-couleur | exacte | 1 (+5 init.) | 92 % (90–93) | 18 |
| BFGS-sans-couleur | exacte | 5 | 99 % (99–100) | 61 |
| Nelder-Mead-sans-couleur | exacte | 1 | 100 % (99–100) | 98 |
| CMA-ES-sans-couleur | exacte | 1 | 100 % (100–100) | 230 |

## Comparaison à fonctions identiques avec le repo complexe

| méthode (6 fonctions identiques au repo complexe) | succès (IC 95 %) | lectures médianes |
|---|---|---|
| sonde-complexe (dépôt d'origine) (8 bits) | 97 % (95–98) | 28 |
| sonde-complexe (dépôt d'origine) (12 bits) | 99 % (98–100) | 31 |
| sonde-couleur-sans-axe-signée (16 bits) | 99 % (98–100) | 34 |
| sonde-couleur-sans-axe-signée (24 bits) | 100 % (99–100) | 19 |

## Succès (%) par fonction — type de cible entre parenthèses

| méthode | q2+1<br>(sphère) | q2-1<br>(2 points réels) | z3-1<br>(point + sphère) | quad<br>(sphère) | poly5<br>(point + 2 sphères) | tanh<br>(point + sphères) | sinz-z<br>(point + sphères) | exp-z-1<br>(point + sphères) | rationnelle<br>(sphère) | sigmoide-0.5<br>(point + sphères) |
|---|---|---|---|---|---|---|---|---|---|---|
| sonde-complexe (dépôt d'origine) (8 bits) | 100 | 100 | 100 | 100 | 100 | 91 | 100 | 100 | 92 | 97 |
| sonde-complexe (dépôt d'origine) (12 bits) | 99 | 100 | 100 | 100 | 100 | 96 | 100 | 100 | 100 | 100 |
| sonde-couleur-sans-axe-signée (16 bits) | 98 | 100 | 100 | 100 | 100 | 98 | 100 | 100 | 100 | 98 |
| sonde-couleur-sans-axe-signée (24 bits) | 100 | 100 | 100 | 99 | 100 | 97 | 100 | 100 | 100 | 100 |
| sonde-couleur-sans-axe (24 bits) | 95 | 82 | 65 | 81 | 68 | 10 | 59 | 60 | 47 | 19 |
| sonde-couleur (4D) (24 bits) | 91 | 95 | 83 | 93 | 85 | 28 | 91 | 93 | 63 | 37 |
| Newton-complexe-sans-couleur | 100 | 100 | 100 | 100 | 100 | 25 | 100 | 100 | 40 | 36 |
| Newton-sans-couleur | 100 | 100 | 100 | 99 | 100 | 65 | 100 | 100 | 90 | 70 |
| Newton-jacobienne-exacte | 100 | 100 | 100 | 99 | 100 | 65 | 100 | 100 | 90 | 70 |
| Broyden-sans-couleur | 100 | 100 | 100 | 99 | 100 | 66 | 100 | 99 | 87 | 68 |
| BFGS-sans-couleur | 100 | 100 | 100 | 100 | 100 | 99 | 100 | 100 | 96 | 99 |
| Nelder-Mead-sans-couleur | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 98 | 100 |
| CMA-ES-sans-couleur | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
