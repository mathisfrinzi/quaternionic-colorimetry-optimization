# Partie 3.1 — tableaux (succès % / lectures médianes)


## 24 bits par composante, bruit 0.0

| ε | méthode | n=1 | n=2 | n=3 | n=4 |
|---|---|---|---|---|---|
| 0.0 | sonde-couleur | 100 / 66 | 98 / 93 | 95 / 138 | 92 / 203 |
| 0.0 | sonde-couleur-quaternionique | 100 / 44 | 100 / 41 | 95 / 43 | 100 / 70 |
| 0.0 | sonde-couleur-mémoire | 100 / 16 | 100 / 24 | 97 / 36 | 90 / 48 |
| 0.0 | sonde-couleur-mémoire-quaternionique | 100 / 11 | 100 / 15 | 100 / 19 | 97 / 23 |
| 0.0 | Broyden-sans-couleur | 100 / 10 | 100 / 17 | 100 / 23 | 100 / 29 |
| 0.0 | Broyden-mêmes-mesures | 100 / 15 | 98 / 27 | 97 / 42 | 92 / 61 |
| 0.0 | Newton-sans-couleur | 100 / 21 | 100 / 37 | 100 / 53 | 100 / 69 |
| 0.0 | BFGS-sans-couleur | 100 / 46 | 100 / 163 | 100 / 339 | 60 / 545 |
| 0.3 | sonde-couleur | 100 / 68 | 100 / 94 | 98 / 168 | 92 / 235 |
| 0.3 | sonde-couleur-quaternionique | 100 / 44 | 92 / 51 | 100 / 72 | 100 / 91 |
| 0.3 | sonde-couleur-mémoire | 100 / 17 | 100 / 25 | 97 / 38 | 87 / 49 |
| 0.3 | sonde-couleur-mémoire-quaternionique | 100 / 11 | 98 / 19 | 93 / 28 | 98 / 34 |
| 0.3 | Broyden-sans-couleur | 100 / 12 | 100 / 19 | 100 / 25 | 100 / 30 |
| 0.3 | Broyden-mêmes-mesures | 98 / 15 | 97 / 28 | 93 / 48 | 95 / 61 |
| 0.3 | Newton-sans-couleur | 100 / 21 | 100 / 37 | 100 / 53 | 100 / 69 |
| 0.3 | BFGS-sans-couleur | 100 / 51 | 100 / 181 | 100 / 365 | 45 / 545 |
| 1.0 | sonde-couleur | 100 / 70 | 78 / 107 | 67 / 192 | 48 / 264 |
| 1.0 | sonde-couleur-quaternionique | 35 / 50 | 0 / — | 2 / 306 | 0 / — |
| 1.0 | sonde-couleur-mémoire | 100 / 18 | 92 / 50 | 72 / 63 | 65 / 86 |
| 1.0 | sonde-couleur-mémoire-quaternionique | 62 / 22 | 15 / 211 | 10 / 252 | 7 / 282 |
| 1.0 | Broyden-sans-couleur | 100 / 14 | 100 / 23 | 100 / 34 | 100 / 44 |
| 1.0 | Broyden-mêmes-mesures | 100 / 17 | 93 / 40 | 85 / 56 | 90 / 71 |
| 1.0 | Newton-sans-couleur | 100 / 21 | 98 / 46 | 98 / 79 | 100 / 86 |
| 1.0 | BFGS-sans-couleur | 100 / 61 | 100 / 244 | 98 / 469 | 2 / 596 |


### Gains appariés (24 bits, σ = 0.0) en % de lectures (+ = la 1ʳᵉ méthode est moins chère), IC 95 %

| ε | comparaison | n=1 | n=2 | n=3 | n=4 |
|---|---|---|---|---|---|
| 0.0 | sonde-couleur-quaternionique contre sonde-couleur | +38 [+34, +41] | +56 [+54, +58] | +70 [+68, +72] | +66 [+64, +68] |
| 0.0 | sonde-couleur-mémoire-quaternionique contre sonde-couleur-mémoire | +34 [+33, +36] | +40 [+35, +44] | +56 [+50, +61] | +57 [+53, +62] |
| 0.0 | sonde-couleur-mémoire-quaternionique contre Broyden-mêmes-mesures | +35 [+31, +38] | +48 [+44, +51] | +55 [+53, +58] | +59 [+55, +62] |
| 0.0 | sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur | -4 [-6, -3] | +6 [+5, +8] | +15 [+13, +17] | +17 [+13, +21] |
| 0.3 | sonde-couleur-quaternionique contre sonde-couleur | +39 [+35, +42] | +44 [+40, +48] | +58 [+56, +59] | +62 [+60, +63] |
| 0.3 | sonde-couleur-mémoire-quaternionique contre sonde-couleur-mémoire | +32 [+29, +34] | +30 [+24, +36] | +27 [+14, +37] | +29 [+16, +40] |
| 0.3 | sonde-couleur-mémoire-quaternionique contre Broyden-mêmes-mesures | +27 [+23, +30] | +37 [+32, +42] | +28 [+18, +37] | +32 [+22, +40] |
| 0.3 | sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur | +7 [+3, +11] | -2 [-7, +3] | -32 [-49, -19] | -36 [-56, -21] |
| 1.0 | sonde-couleur-quaternionique contre sonde-couleur | +31 [+22, +39] | — | — | — |
| 1.0 | sonde-couleur-mémoire-quaternionique contre sonde-couleur-mémoire | -120 [-204, -66] | -379 [-507, -259] | -191 [-458, -86] | — |
| 1.0 | sonde-couleur-mémoire-quaternionique contre Broyden-mêmes-mesures | -126 [-217, -68] | -447 [-591, -331] | -353 [-429, -289] | — |
| 1.0 | sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur | -181 [-279, -104] | -771 [-953, -620] | -734 [-868, -636] | — |


## 24 bits par composante, bruit 0.02

| ε | méthode | n=1 | n=2 | n=3 | n=4 |
|---|---|---|---|---|---|
| 0.0 | sonde-couleur | 98 / 69 | 98 / 101 | 83 / 178 | 60 / 254 |
| 0.0 | sonde-couleur-quaternionique | 100 / 44 | 98 / 44 | 98 / 57 | 85 / 76 |
| 0.0 | sonde-couleur-mémoire | 100 / 17 | 100 / 26 | 90 / 54 | 83 / 80 |
| 0.0 | sonde-couleur-mémoire-quaternionique | 100 / 11 | 100 / 16 | 98 / 22 | 97 / 29 |
| 0.0 | Broyden-sans-couleur | 100 / 10 | 100 / 17 | 100 / 23 | 100 / 29 |
| 0.0 | Broyden-mêmes-mesures | 100 / 17 | 100 / 44 | 100 / 72 | 100 / 114 |
| 0.0 | Newton-sans-couleur | 100 / 21 | 100 / 37 | 100 / 53 | 100 / 69 |
| 0.0 | BFGS-sans-couleur | 100 / 46 | 100 / 163 | 100 / 339 | 60 / 545 |
| 0.3 | sonde-couleur | 100 / 68 | 97 / 114 | 85 / 176 | 53 / 252 |
| 0.3 | sonde-couleur-quaternionique | 100 / 43 | 100 / 60 | 98 / 74 | 97 / 96 |
| 0.3 | sonde-couleur-mémoire | 100 / 17 | 100 / 30 | 92 / 64 | 90 / 76 |
| 0.3 | sonde-couleur-mémoire-quaternionique | 100 / 11 | 100 / 21 | 98 / 39 | 92 / 43 |
| 0.3 | Broyden-sans-couleur | 100 / 12 | 100 / 19 | 100 / 25 | 100 / 30 |
| 0.3 | Broyden-mêmes-mesures | 100 / 18 | 100 / 47 | 100 / 87 | 100 / 113 |
| 0.3 | Newton-sans-couleur | 100 / 21 | 100 / 37 | 100 / 53 | 100 / 69 |
| 0.3 | BFGS-sans-couleur | 100 / 51 | 100 / 181 | 100 / 365 | 45 / 545 |
| 1.0 | sonde-couleur | 97 / 69 | 50 / 125 | 32 / 220 | 15 / 289 |
| 1.0 | sonde-couleur-quaternionique | 67 / 88 | 0 / — | 0 / — | 0 / — |
| 1.0 | sonde-couleur-mémoire | 100 / 19 | 88 / 56 | 73 / 70 | 68 / 87 |
| 1.0 | sonde-couleur-mémoire-quaternionique | 58 / 22 | 5 / 141 | 0 / — | 0 / — |
| 1.0 | Broyden-sans-couleur | 100 / 14 | 100 / 23 | 100 / 34 | 100 / 44 |
| 1.0 | Broyden-mêmes-mesures | 100 / 22 | 100 / 56 | 100 / 92 | 100 / 131 |
| 1.0 | Newton-sans-couleur | 100 / 21 | 98 / 46 | 98 / 79 | 100 / 86 |
| 1.0 | BFGS-sans-couleur | 100 / 61 | 100 / 244 | 98 / 469 | 2 / 596 |


### Gains appariés (24 bits, σ = 0.02) en % de lectures (+ = la 1ʳᵉ méthode est moins chère), IC 95 %

| ε | comparaison | n=1 | n=2 | n=3 | n=4 |
|---|---|---|---|---|---|
| 0.0 | sonde-couleur-quaternionique contre sonde-couleur | +41 [+38, +45] | +58 [+56, +61] | +69 [+67, +70] | +70 [+68, +72] |
| 0.0 | sonde-couleur-mémoire-quaternionique contre sonde-couleur-mémoire | +36 [+34, +37] | +49 [+43, +55] | +56 [+51, +61] | +58 [+52, +63] |
| 0.0 | sonde-couleur-mémoire-quaternionique contre Broyden-mêmes-mesures | +45 [+41, +49] | +63 [+61, +66] | +70 [+68, +73] | +72 [+68, +75] |
| 0.0 | sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur | -6 [-9, -4] | -1 [-5, +3] | +0 [-6, +5] | -8 [-18, +1] |
| 0.3 | sonde-couleur-quaternionique contre sonde-couleur | +42 [+38, +45] | +44 [+39, +49] | +58 [+56, +60] | +61 [+58, +64] |
| 0.3 | sonde-couleur-mémoire-quaternionique contre sonde-couleur-mémoire | +30 [+25, +33] | +42 [+33, +50] | +25 [+13, +36] | +32 [+23, +39] |
| 0.3 | sonde-couleur-mémoire-quaternionique contre Broyden-mêmes-mesures | +38 [+32, +44] | +55 [+50, +59] | +50 [+42, +57] | +58 [+52, +63] |
| 0.3 | sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur | +3 [-4, +7] | -19 [-30, -9] | -79 [-105, -58] | -59 [-77, -44] |
| 1.0 | sonde-couleur-quaternionique contre sonde-couleur | -27 [-39, -18] | — | — | — |
| 1.0 | sonde-couleur-mémoire-quaternionique contre sonde-couleur-mémoire | -106 [-173, -58] | — | — | — |
| 1.0 | sonde-couleur-mémoire-quaternionique contre Broyden-mêmes-mesures | -77 [-133, -34] | — | — | — |
| 1.0 | sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur | -190 [-298, -121] | — | — | — |
