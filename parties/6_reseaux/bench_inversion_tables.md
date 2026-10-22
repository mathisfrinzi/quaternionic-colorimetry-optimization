# Partie 6.2 : inversion d'un réseau, tableaux

Succès [IC 95 %] / lectures médianes, 3 instances × 20 départs.


## σ = 0,0

| activation | n | défaut ℍ-lin. | sonde-couleur-mémoire-quaternionique | sonde-couleur-mémoire | Broyden-sans-couleur | Broyden-mêmes-mesures |
|---|---|---|---|---|---|---|
| tanh par composante | 2 | 0,39 | 45 % [33 ; 58] / 56 | 83 % [72 ; 91] / 25 | 90 % [80 ; 95] / 22 | 10 % [5 ; 20] / 90 |
| tanh par composante | 3 | 0,44 | 10 % [5 ; 20] / 220 | 28 % [19 ; 41] / 79 | 98 % [91 ; 100] / 31 | 10 % [5 ; 20] / 244 |
| tanh par composante | 4 | 0,55 | 0 % [0 ; 6] | 17 % [9 ; 28] / 84 | 82 % [70 ; 89] / 42 | 2 % [0 ; 9] / 479 |
| tanh par composante | 8 | 0,44 | 0 % [0 ; 6] | — | 65 % [52 ; 76] / 89 | 0 % [0 ; 6] |
| tanh radial | 2 | 0,35 | 42 % [30 ; 54] / 39 | 87 % [76 ; 93] / 28 | 95 % [86 ; 98] / 24 | 5 % [2 ; 14] / 135 |
| tanh radial | 3 | 0,41 | 3 % [1 ; 11] / 188 | 43 % [32 ; 56] / 74 | 88 % [78 ; 94] / 37 | 10 % [5 ; 20] / 230 |
| tanh radial | 4 | 0,44 | 3 % [1 ; 11] / 226 | 15 % [8 ; 26] / 90 | 72 % [59 ; 81] / 56 | 8 % [4 ; 18] / 439 |
| tanh radial | 8 | 0,38 | 0 % [0 ; 6] | — | 77 % [65 ; 86] / 114 | 0 % [0 ; 6] |
| saturation globale | 2 | 0,27 | 50 % [38 ; 62] / 26 | 95 % [86 ; 98] / 21 | 100 % [94 ; 100] / 18 | 5 % [2 ; 14] / 62 |
| saturation globale | 3 | 0,26 | 17 % [9 ; 28] / 60 | 82 % [70 ; 89] / 31 | 100 % [94 ; 100] / 24 | 12 % [6 ; 22] / 163 |
| saturation globale | 4 | 0,25 | 5 % [2 ; 14] / 101 | 48 % [36 ; 61] / 40 | 93 % [84 ; 97] / 29 | 7 % [3 ; 16] / 348 |
| saturation globale | 8 | 0,13 | 3 % [1 ; 11] / 62 | — | 100 % [94 ; 100] / 43 | 0 % [0 ; 6] |

## σ = 0,02

| activation | n | défaut ℍ-lin. | sonde-couleur-mémoire-quaternionique | sonde-couleur-mémoire | Broyden-sans-couleur | Broyden-mêmes-mesures |
|---|---|---|---|---|---|---|
| tanh par composante | 2 | 0,39 | 45 % [33 ; 58] / 63 | 77 % [65 ; 86] / 28 | 90 % [80 ; 95] / 22 | 22 % [13 ; 34] / 312 |
| tanh par composante | 3 | 0,44 | 2 % [0 ; 9] / 238 | 20 % [12 ; 32] / 80 | 98 % [91 ; 100] / 31 | 2 % [0 ; 9] / 600 |
| tanh par composante | 4 | 0,55 | 0 % [0 ; 6] | 5 % [2 ; 14] / 94 | 82 % [70 ; 89] / 42 | 0 % [0 ; 6] |
| tanh par composante | 8 | 0,44 | 0 % [0 ; 6] | — | 65 % [52 ; 76] / 89 | 0 % [0 ; 6] |
| tanh radial | 2 | 0,35 | 43 % [32 ; 56] / 66 | 80 % [68 ; 88] / 46 | 95 % [86 ; 98] / 24 | 18 % [11 ; 30] / 354 |
| tanh radial | 3 | 0,41 | 0 % [0 ; 6] | 28 % [19 ; 41] / 75 | 88 % [78 ; 94] / 37 | 5 % [2 ; 14] / 449 |
| tanh radial | 4 | 0,44 | 0 % [0 ; 6] | 5 % [2 ; 14] / 122 | 72 % [59 ; 81] / 56 | 0 % [0 ; 6] |
| tanh radial | 8 | 0,38 | 0 % [0 ; 6] | — | 77 % [65 ; 86] / 114 | 0 % [0 ; 6] |
| saturation globale | 2 | 0,27 | 50 % [38 ; 62] / 31 | 92 % [82 ; 96] / 25 | 100 % [94 ; 100] / 18 | 17 % [9 ; 28] / 402 |
| saturation globale | 3 | 0,26 | 2 % [0 ; 9] / 134 | 68 % [56 ; 79] / 36 | 100 % [94 ; 100] / 24 | 7 % [3 ; 16] / 538 |
| saturation globale | 4 | 0,25 | 5 % [2 ; 14] / 18 | 45 % [33 ; 58] / 58 | 93 % [84 ; 97] / 29 | 0 % [0 ; 6] |
| saturation globale | 8 | 0,13 | 2 % [0 ; 9] / 53 | — | 100 % [94 ; 100] / 43 | 0 % [0 ; 6] |

## Gain apparié de la sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur (σ = 0)

| activation | n = 2 | n = 3 | n = 4 | n = 8 |
|---|---|---|---|---|
| tanh par composante | -243 [-486 ; -114] (26) | -684 [-1004 ; -330] (6) | — | — |
| tanh radial | -86 [-182 ; -33] (25) | — | — | — |
| saturation globale | -50 [-124 ; +0] (30) | -153 [-268 ; -29] (10) | -288 [-996 ; -227] (3) | — |
