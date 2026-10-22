# Partie 1.3 — tableaux (générés par bench_mobius_figs.py)

Budget 600 lectures ; 5 instances × 30 départs par ε. Défaut = moyenne du défaut de conformité au zéro.


## Succès % / lectures médianes, 24 bits

| méthode | ε=0.0 (déf 0.00) | ε=0.03 (déf 0.02) | ε=0.1 (déf 0.07) | ε=0.3 (déf 0.20) | ε=1.0 (déf 0.50) | ε=3.0 (déf 0.61) |
|---|---|---|---|---|---|---|
| sonde-couleur (24 bits) | 100 / 45 | 100 / 45 | 99 / 45 | 97 / 49 | 49 / 76 | 37 / 143 |
| sonde-couleur-conforme-3-sondes (24 bits) | 100 / 28 | 100 / 28 | 100 / 28 | 100 / 40 | 82 / 233 | 36 / 404 |
| sonde-couleur-conforme-projetée (24 bits) | 100 / 36 | 100 / 36 | 100 / 37 | 100 / 44 | 100 / 203 | 87 / 274 |
| sonde-couleur-mémoire (24 bits) | 100 / 18 | 100 / 18 | 100 / 18 | 100 / 18 | 83 / 24 | 99 / 32 |
| sonde-couleur-mémoire-multiple (24 bits) | 100 / 19 | 100 / 18 | 100 / 19 | 100 / 19 | 81 / 26 | 99 / 31 |
| sonde-couleur-mémoire-conforme (24 bits) | 100 / 17 | 100 / 17 | 100 / 18 | 100 / 18 | 82 / 49 | 47 / 130 |
| Broyden-sans-couleur | 100 / 13 | 100 / 14 | 100 / 14 | 99 / 15 | 74 / 20 | 100 / 19 |
| Newton-sans-couleur | 100 / 21 | 100 / 21 | 100 / 21 | 99 / 21 | 67 / 26 | 99 / 26 |
| Newton-jacobienne-exacte | 100 / 9 | 100 / 9 | 100 / 9 | 99 / 9 | 67 / 11 | 99 / 12 |
| BFGS-sans-couleur | 100 / 36 | 100 / 36 | 100 / 41 | 100 / 46 | 81 / 56 | 100 / 81 |
| Nelder-Mead-sans-couleur | 100 / 110 | 100 / 109 | 100 / 109 | 99 / 110 | 85 / 128 | 100 / 163 |


## Succès % / lectures médianes, 36 bits

| méthode | ε=0.0 (déf 0.00) | ε=0.03 (déf 0.02) | ε=0.1 (déf 0.07) | ε=0.3 (déf 0.20) | ε=1.0 (déf 0.50) | ε=3.0 (déf 0.61) |
|---|---|---|---|---|---|---|
| sonde-couleur (36 bits) | 100 / 27 | 100 / 27 | 100 / 27 | 99 / 27 | 79 / 32 | 96 / 46 |
| sonde-couleur-conforme-3-sondes (36 bits) | 100 / 45 | 100 / 45 | 100 / 45 | 100 / 57 | 51 / 230 | 3 / 274 |
| sonde-couleur-conforme-projetée (36 bits) | 100 / 26 | 100 / 26 | 100 / 28 | 100 / 36 | 98 / 214 | 77 / 250 |
| sonde-couleur-mémoire (36 bits) | 100 / 15 | 100 / 15 | 100 / 15 | 100 / 15 | 80 / 17 | 100 / 18 |
| sonde-couleur-mémoire-multiple (36 bits) | 100 / 18 | 100 / 17 | 100 / 17 | 100 / 18 | 81 / 21 | 100 / 22 |
| sonde-couleur-mémoire-conforme (36 bits) | 100 / 15 | 100 / 15 | 100 / 15 | 100 / 16 | 76 / 34 | 48 / 86 |
| Broyden-sans-couleur | 100 / 13 | 100 / 14 | 100 / 14 | 99 / 15 | 74 / 20 | 100 / 19 |
| Newton-sans-couleur | 100 / 21 | 100 / 21 | 100 / 21 | 99 / 21 | 67 / 26 | 99 / 26 |
| Newton-jacobienne-exacte | 100 / 9 | 100 / 9 | 100 / 9 | 99 / 9 | 67 / 11 | 99 / 12 |
| BFGS-sans-couleur | 100 / 36 | 100 / 36 | 100 / 41 | 100 / 46 | 81 / 56 | 100 / 81 |
| Nelder-Mead-sans-couleur | 100 / 110 | 100 / 109 | 100 / 109 | 99 / 110 | 85 / 128 | 100 / 163 |


## Gain apparié (en % de lectures, + = la 2ᵉ méthode est moins chère), IC95 % bootstrap, n = départs réussis par les deux

| bits | comparaison | ε | défaut | gain % | IC95 | n | Δ succès (pts) |
|---|---|---|---|---|---|---|---|
| 24 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.0 | 0.00 | 41 | [40, 42] | 150 | -0 |
| 24 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.03 | 0.02 | 42 | [40, 43] | 150 | -0 |
| 24 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.1 | 0.07 | 39 | [38, 41] | 148 | +1 |
| 24 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.3 | 0.20 | 12 | [7, 17] | 145 | +3 |
| 24 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 1.0 | 0.50 | -202 | [-247, -163] | 62 | +33 |
| 24 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 3.0 | 0.61 | -145 | [-204, -93] | 19 | -1 |
| 24 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.0 | 0.00 | 21 | [20, 23] | 150 | -0 |
| 24 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.03 | 0.02 | 22 | [20, 24] | 150 | -0 |
| 24 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.1 | 0.07 | 20 | [18, 22] | 148 | +1 |
| 24 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.3 | 0.20 | 11 | [8, 14] | 145 | +3 |
| 24 | sonde-couleur-conforme-projetée vs sonde-couleur | 1.0 | 0.50 | -142 | [-173, -111] | 74 | +51 |
| 24 | sonde-couleur-conforme-projetée vs sonde-couleur | 3.0 | 0.61 | -83 | [-121, -52] | 47 | +50 |
| 24 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.0 | 0.00 | 7 | [4, 11] | 150 | -0 |
| 24 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.03 | 0.02 | 7 | [4, 11] | 150 | -0 |
| 24 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.1 | 0.07 | 6 | [3, 10] | 150 | -0 |
| 24 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.3 | 0.20 | 8 | [3, 13] | 150 | -0 |
| 24 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 1.0 | 0.50 | -44 | [-76, -18] | 100 | -1 |
| 24 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 3.0 | 0.61 | -219 | [-294, -159] | 71 | -51 |
| 24 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.0 | 0.00 | -21 | [-31, -12] | 150 | -0 |
| 24 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.03 | 0.02 | -19 | [-29, -10] | 150 | -0 |
| 24 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.1 | 0.07 | -24 | [-35, -15] | 150 | -0 |
| 24 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.3 | 0.20 | -11 | [-21, -2] | 150 | -0 |
| 24 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 1.0 | 0.50 | -3 | [-19, 11] | 109 | -2 |
| 24 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 3.0 | 0.61 | -1 | [-16, 14] | 147 | +1 |
| 36 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.0 | 0.00 | -60 | [-62, -57] | 150 | -0 |
| 36 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.03 | 0.02 | -59 | [-61, -56] | 150 | -0 |
| 36 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.1 | 0.07 | -63 | [-66, -61] | 150 | -0 |
| 36 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 0.3 | 0.20 | -101 | [-106, -97] | 149 | +1 |
| 36 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 1.0 | 0.50 | -447 | [-536, -368] | 64 | -29 |
| 36 | sonde-couleur-conforme-3-sondes vs sonde-couleur | 3.0 | 0.61 | — | — | 4 | -93 |
| 36 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.0 | 0.00 | 1 | [-0, 2] | 150 | -0 |
| 36 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.03 | 0.02 | 1 | [-1, 2] | 150 | -0 |
| 36 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.1 | 0.07 | -4 | [-5, -2] | 150 | -0 |
| 36 | sonde-couleur-conforme-projetée vs sonde-couleur | 0.3 | 0.20 | -33 | [-36, -30] | 149 | +1 |
| 36 | sonde-couleur-conforme-projetée vs sonde-couleur | 1.0 | 0.50 | -462 | [-505, -421] | 119 | +19 |
| 36 | sonde-couleur-conforme-projetée vs sonde-couleur | 3.0 | 0.61 | -523 | [-583, -469] | 111 | -19 |
| 36 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.0 | 0.00 | 3 | [2, 5] | 150 | -0 |
| 36 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.03 | 0.02 | 3 | [1, 5] | 150 | -0 |
| 36 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.1 | 0.07 | 2 | [1, 3] | 150 | -0 |
| 36 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 0.3 | 0.20 | -1 | [-3, 2] | 150 | -0 |
| 36 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 1.0 | 0.50 | -90 | [-119, -64] | 99 | -4 |
| 36 | sonde-couleur-mémoire-conforme vs sonde-couleur-mémoire | 3.0 | 0.61 | -332 | [-415, -265] | 72 | -52 |
| 36 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.0 | 0.00 | -18 | [-23, -15] | 150 | -0 |
| 36 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.03 | 0.02 | -17 | [-21, -14] | 150 | -0 |
| 36 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.1 | 0.07 | -19 | [-22, -16] | 150 | -0 |
| 36 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 0.3 | 0.20 | -19 | [-23, -15] | 150 | -0 |
| 36 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 1.0 | 0.50 | -10 | [-18, -3] | 117 | +1 |
| 36 | sonde-couleur-mémoire-multiple vs sonde-couleur-mémoire | 3.0 | 0.61 | -5 | [-16, 4] | 150 | -0 |
