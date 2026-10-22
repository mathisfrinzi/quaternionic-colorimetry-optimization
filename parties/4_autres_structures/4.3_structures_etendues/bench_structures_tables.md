# Partie 4.3 : tableaux

Succès [IC 95 %] / lectures médianes, 3 instances × 20 départs, budget 600.


## σ = 0,0

| famille | n | défaut ℍ-lin. | modèle étendu | sonde-couleur-mémoire-quaternionique | sonde-couleur-mémoire | Broyden-sans-couleur | Broyden-mêmes-mesures |
|---|---|---|---|---|---|---|---|
| bilatère | 2 | 0,91 | 98 % [91 ; 100] / 32 | 0 % [0 ; 6] | 100 % [94 ; 100] / 24 | 100 % [94 ; 100] / 15 | 23 % [14 ; 35] / 124 |
| bilatère | 3 | 0,91 | 97 % [89 ; 99] / 64 | 0 % [0 ; 6] | 97 % [89 ; 99] / 42 | 100 % [94 ; 100] / 20 | 30 % [20 ; 43] / 219 |
| bilatère | 4 | 0,77 | 97 % [89 ; 99] / 104 | 0 % [0 ; 6] | 92 % [82 ; 96] / 80 | 100 % [94 ; 100] / 24 | 18 % [11 ; 30] / 271 |
| ℍ-linéaire + rang 2 | 2 | 0,44 | 85 % [74 ; 92] / 87 | 37 % [26 ; 49] / 70 | 98 % [91 ; 100] / 32 | 100 % [94 ; 100] / 21 | 12 % [6 ; 22] / 148 |
| ℍ-linéaire + rang 2 | 3 | 0,44 | 73 % [61 ; 83] / 80 | 55 % [42 ; 67] / 56 | 87 % [76 ; 93] / 44 | 100 % [94 ; 100] / 26 | 23 % [14 ; 35] / 262 |
| ℍ-linéaire + rang 2 | 4 | 0,46 | 97 % [89 ; 99] / 72 | 55 % [42 ; 67] / 114 | 82 % [70 ; 89] / 75 | 100 % [94 ; 100] / 31 | 18 % [11 ; 30] / 356 |
| inversion réseau (tanh) | 2 | 0,39 | 73 % [61 ; 83] / 50 | 45 % [33 ; 58] / 56 | 83 % [72 ; 91] / 25 | 90 % [80 ; 95] / 22 | 10 % [5 ; 20] / 90 |
| inversion réseau (tanh) | 3 | 0,44 | 43 % [32 ; 56] / 69 | 10 % [5 ; 20] / 220 | 28 % [19 ; 41] / 79 | 98 % [91 ; 100] / 31 | 10 % [5 ; 20] / 244 |
| inversion réseau (tanh) | 4 | 0,55 | 8 % [4 ; 18] / 78 | 0 % [0 ; 6] | 17 % [9 ; 28] / 84 | 82 % [70 ; 89] / 42 | 2 % [0 ; 9] / 479 |

## σ = 0,02

| famille | n | défaut ℍ-lin. | modèle étendu | sonde-couleur-mémoire-quaternionique | sonde-couleur-mémoire | Broyden-sans-couleur | Broyden-mêmes-mesures |
|---|---|---|---|---|---|---|---|
| bilatère | 2 | 0,91 | 100 % [94 ; 100] / 42 | 0 % [0 ; 6] | 100 % [94 ; 100] / 28 | 100 % [94 ; 100] / 15 | 48 % [36 ; 61] / 362 |
| bilatère | 3 | 0,91 | 97 % [89 ; 99] / 66 | 0 % [0 ; 6] | 90 % [80 ; 95] / 72 | 100 % [94 ; 100] / 20 | 0 % [0 ; 6] |
| bilatère | 4 | 0,77 | 88 % [78 ; 94] / 74 | 0 % [0 ; 6] | 97 % [89 ; 99] / 86 | 100 % [94 ; 100] / 24 | 0 % [0 ; 6] |
| ℍ-linéaire + rang 2 | 2 | 0,44 | 82 % [70 ; 89] / 86 | 52 % [39 ; 64] / 117 | 93 % [84 ; 97] / 52 | 100 % [94 ; 100] / 21 | 42 % [30 ; 54] / 394 |
| ℍ-linéaire + rang 2 | 3 | 0,44 | 62 % [49 ; 73] / 88 | 48 % [36 ; 61] / 52 | 77 % [65 ; 86] / 68 | 100 % [94 ; 100] / 26 | 2 % [0 ; 9] / 591 |
| ℍ-linéaire + rang 2 | 4 | 0,46 | 77 % [65 ; 86] / 78 | 35 % [24 ; 48] / 124 | 78 % [66 ; 87] / 87 | 100 % [94 ; 100] / 31 | 0 % [0 ; 6] |
| inversion réseau (tanh) | 2 | 0,39 | 55 % [42 ; 67] / 73 | 45 % [33 ; 58] / 63 | 77 % [65 ; 86] / 28 | 90 % [80 ; 95] / 22 | 22 % [13 ; 34] / 312 |
| inversion réseau (tanh) | 3 | 0,44 | 13 % [7 ; 24] / 67 | 2 % [0 ; 9] / 238 | 20 % [12 ; 32] / 80 | 98 % [91 ; 100] / 31 | 2 % [0 ; 9] / 600 |
| inversion réseau (tanh) | 4 | 0,55 | 2 % [0 ; 9] / 88 | 0 % [0 ; 6] | 5 % [2 ; 14] / 94 | 82 % [70 ; 89] / 42 | 0 % [0 ; 6] |

## Gains appariés du modèle étendu (σ = 0)

| famille | comparaison | n = 2 | n = 3 | n = 4 |
|---|---|---|---|---|
| bilatère | contre sonde-couleur-mémoire-quaternionique | — | — | — |
| bilatère | contre sonde-couleur-mémoire | -33 [-39 ; -25] (59) | -55 [-86 ; -17] (56) | -33 [-61 ; -9] (54) |
| bilatère | contre Broyden-sans-couleur | -113 [-120 ; -100] (59) | -222 [-262 ; -180] (58) | -331 [-375 ; -271] (58) |
| ℍ-linéaire + rang 2 | contre sonde-couleur-mémoire-quaternionique | -20 [-84 ; +13] (20) | -18 [-58 ; +23] (31) | +34 [-1 ; +55] (32) |
| ℍ-linéaire + rang 2 | contre sonde-couleur-mémoire | -169 [-253 ; -54] (50) | -90 [-188 ; +3] (41) | +5 [-51 ; +27] (47) |
| ℍ-linéaire + rang 2 | contre Broyden-sans-couleur | -314 [-419 ; -205] (51) | -222 [-356 ; -131] (44) | -132 [-169 ; -88] (58) |
| inversion réseau (tanh) | contre sonde-couleur-mémoire-quaternionique | +25 [-29 ; +68] (24) | — | — |
| inversion réseau (tanh) | contre sonde-couleur-mémoire | -102 [-173 ; -30] (38) | -29 [-355 ; +52] (5) | — |
| inversion réseau (tanh) | contre Broyden-sans-couleur | -116 [-202 ; -77] (40) | -126 [-260 ; -60] (26) | -129 [-265 ; -59] (5) |
