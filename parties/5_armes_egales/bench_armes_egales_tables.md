# Partie 5 — tableaux

| méthode | bits | σ | succès % | lectures médianes |
|---|---|---|---|---|
| sonde-couleur-mémoire | 24 | 0.0 | 100 | 19 |
| sonde-couleur | 24 | 0.0 | 95 | 36 |
| Broyden-mêmes-mesures | 24 | 0.0 | 98 | 13 |
| Newton-mêmes-mesures | 24 | 0.0 | 97 | 26 |
| BFGS-mêmes-mesures | 24 | 0.0 | 0 | nan |
| Nelder-Mead-mêmes-mesures | 24 | 0.0 | 82 | 141 |
| sonde-couleur-mémoire | 24 | 0.005 | 99 | 19 |
| sonde-couleur | 24 | 0.005 | 98 | 41 |
| Broyden-mêmes-mesures | 24 | 0.005 | 99 | 13 |
| Newton-mêmes-mesures | 24 | 0.005 | 96 | 26 |
| BFGS-mêmes-mesures | 24 | 0.005 | 90 | 71 |
| Nelder-Mead-mêmes-mesures | 24 | 0.005 | 89 | 136 |
| sonde-couleur-mémoire | 24 | 0.02 | 100 | 16 |
| sonde-couleur | 24 | 0.02 | 99 | 46 |
| Broyden-mêmes-mesures | 24 | 0.02 | 99 | 26 |
| Newton-mêmes-mesures | 24 | 0.02 | 77 | 36 |
| BFGS-mêmes-mesures | 24 | 0.02 | 0 | nan |
| Nelder-Mead-mêmes-mesures | 24 | 0.02 | 82 | 134 |
| sonde-couleur-mémoire | 36 | 0.0 | 92 | 13 |
| sonde-couleur | 36 | 0.0 | 99 | 65 |
| Broyden-mêmes-mesures | 36 | 0.0 | 100 | 12 |
| Newton-mêmes-mesures | 36 | 0.0 | 98 | 21 |
| BFGS-mêmes-mesures | 36 | 0.0 | 100 | 51 |
| Nelder-Mead-mêmes-mesures | 36 | 0.0 | 99 | 139 |
| sonde-couleur-mémoire | 36 | 0.005 | 100 | 14 |
| sonde-couleur | 36 | 0.005 | 100 | 46 |
| Broyden-mêmes-mesures | 36 | 0.005 | 100 | 12 |
| Newton-mêmes-mesures | 36 | 0.005 | 95 | 21 |
| BFGS-mêmes-mesures | 36 | 0.005 | 98 | 56 |
| Nelder-Mead-mêmes-mesures | 36 | 0.005 | 99 | 137 |
| sonde-couleur-mémoire | 36 | 0.02 | 99 | 19 |
| sonde-couleur | 36 | 0.02 | 100 | 36 |
| Broyden-mêmes-mesures | 36 | 0.02 | 100 | 14 |
| Newton-mêmes-mesures | 36 | 0.02 | 82 | 31 |
| BFGS-mêmes-mesures | 36 | 0.02 | 71 | 81 |
| Nelder-Mead-mêmes-mesures | 36 | 0.02 | 80 | 140 |


## Apparié : sonde-rank1 contre Broyden (gain en % de lectures, + = la sonde est moins chère ; IC95 % bootstrap)

| bits | σ | n | gain % | IC95 |
|---|---|---|---|---|
| 24 | 0.0 | 350 | -43 | [-46, -39] |
| 24 | 0.005 | 355 | -42 | [-46, -38] |
| 24 | 0.02 | 356 | 30 | [26, 33] |
| 36 | 0.0 | 332 | -17 | [-20, -15] |
| 36 | 0.005 | 360 | -24 | [-27, -21] |
| 36 | 0.02 | 357 | -34 | [-38, -30] |
