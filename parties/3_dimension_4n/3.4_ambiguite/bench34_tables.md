# Partie 3.4 : tableaux (générés par bench34_figs.py)

## Calibration (n = 32, graines 300–304 × 6 départs = 30 départs, jamais réutilisées)

| configuration | σ = 0 : succès / lectures méd. | σ = 0,02 : succès / lectures méd. |
|---|---|---|
| R1 k=2 tau=0.25 double=1 repli=1 maj=0 | 87 % / 140 | 70 % / 219 |
| R1 k=4 tau=0.25 double=1 repli=1 maj=0 | 87 % / 142 | 70 % / 207 |
| R1 k=2 tau=0.25 double=1 repli=0 maj=0 | 87 % / 141 | 67 % / 196 |
| R1 k=2 tau=0.4 double=1 repli=1 maj=0 | 93 % / 137 | 60 % / 207 |
| R1 k=4 tau=0.25 double=1 repli=0 maj=0 | 90 % / 137 | 63 % / 204 |
| R1 k=2 tau=0.4 double=1 repli=1 maj=1 | 93 % / 134 | 60 % / 208 |
| R1 k=2 tau=0.4 double=1 repli=0 maj=0 | 80 % / 135 | 70 % / 211 |
| R1 k=4 tau=0.4 double=1 repli=0 maj=0 | 87 % / 135 | 57 % / 205 |
| R1 k=4 tau=0.4 double=1 repli=1 maj=0 | 87 % / 138 | 53 % / 210 |
| R1 k=2 tau=0.25 double=0 repli=0 maj=1 | 40 % / 152 | – |
| R1 k=2 tau=0.25 double=0 repli=0 maj=0 | 37 % / 148 | – |
| R1 k=2 tau=0.25 double=0 repli=1 maj=0 | 33 % / 143 | – |
| R1 k=2 tau=0.25 double=0 repli=1 maj=1 | 20 % / 148 | – |

## Banc final (graines 100–109 × 10 départs = 100 départs par case)

### n = 16, σ = 0

| méthode | succès [IC 95 %] | lectures médianes [IC 95 %] | écart apparié à l'oracle [IC 95 %] (paires) | rejets ambiguïté / descente par lecture de pas | relectures par pas | pas acceptés avec une décade fausse muette |
|---|---|---|---|---|---|---|
| sonde-garde-ambiguïté | 82 % [73 ; 88] (100) | 100 [94 ; 112] | +45 [41 ; 54] (82) | 0.37 / 0.13 | 0.00 | 2.1 % |
| sonde-dépliage-vers-zéro | 98 % [93 ; 99] (100) | 63 [61 ; 65] | +10 [9.0 ; 11] (98) | 0.19 / 0.08 | 0.00 | 1.5 % |
| sonde-relecture-naïve | 64 % [54 ; 73] (100) | 87 [84 ; 92] | +34 [32 ; 35] (64) | 0.00 / 0.28 | 0.73 | 3.0 % |
| sonde-relecture | 99 % [95 ; 100] (100) | 87 [86 ; 89] | +34 [33 ; 34] (99) | 0.01 / 0.11 | 0.90 | 0.1 % |
| sonde-par-composante | 9 % [5 ; 16] (100) | 66 [64 ; 75] | +13 [8.0 ; 18] (9) | 0.00 / 0.85 | 0.00 | 28.4 % |
| sonde-relecture-par-composante | 96 % [90 ; 98] (100) | 87 [85 ; 90] | +33 [32 ; 34] (96) | 0.00 / 0.12 | 0.90 | 0.1 % |
| sonde-oracle-décade (triche) | 100 % [96 ; 100] (100) | 53 [52 ; 54] | réf. | 0.00 / 0.11 | 0.00 | 0.0 % |
| Broyden-sans-couleur | 100 % [96 ; 100] (100) | 84 [83 ; 85] | +31 [30 ; 32] (100) | – | – | – |

### n = 16, σ = 0.02

| méthode | succès [IC 95 %] | lectures médianes [IC 95 %] | écart apparié à l'oracle [IC 95 %] (paires) | rejets ambiguïté / descente par lecture de pas | relectures par pas | pas acceptés avec une décade fausse muette |
|---|---|---|---|---|---|---|
| sonde-garde-ambiguïté | 72 % [63 ; 80] (100) | 122 [112 ; 130] | +45 [34 ; 52] (72) | 0.30 / 0.47 | 0.00 | 1.9 % |
| sonde-dépliage-vers-zéro | 91 % [84 ; 95] (100) | 108 [101 ; 115] | +28 [22 ; 33] (91) | 0.20 / 0.39 | 0.00 | 1.2 % |
| sonde-relecture-naïve | 69 % [59 ; 77] (100) | 126 [118 ; 149] | +55 [43 ; 70] (69) | 0.00 / 0.70 | 0.48 | 6.2 % |
| sonde-relecture | 95 % [89 ; 98] (100) | 127 [123 ; 132] | +51 [44 ; 56] (95) | 0.01 / 0.44 | 0.68 | 0.3 % |
| sonde-par-composante | 0 % [0 ; 4] (100) | – [– ; –] | – | 0.00 / 0.96 | 0.00 | 56.5 % |
| sonde-relecture-par-composante | 97 % [92 ; 99] (100) | 126 [122 ; 133] | +47 [43 ; 55] (97) | 0.00 / 0.42 | 0.70 | 0.4 % |
| sonde-oracle-décade (triche) | 100 % [96 ; 100] (100) | 77 [74 ; 80] | réf. | 0.00 / 0.39 | 0.00 | 0.0 % |
| Broyden-sans-couleur | 100 % [96 ; 100] (100) | 84 [83 ; 85] | +11 [7.0 ; 13] (100) | – | – | – |

### n = 32, σ = 0

| méthode | succès [IC 95 %] | lectures médianes [IC 95 %] | écart apparié à l'oracle [IC 95 %] (paires) | rejets ambiguïté / descente par lecture de pas | relectures par pas | pas acceptés avec une décade fausse muette |
|---|---|---|---|---|---|---|
| sonde-garde-ambiguïté | 86 % [78 ; 91] (100) | 172 [167 ; 179] | +80 [73 ; 87] (81) | 0.45 / 0.17 | 0.00 | 1.6 % |
| sonde-dépliage-vers-zéro | 92 % [85 ; 96] (100) | 113 [109 ; 124] | +26 [20 ; 32] (88) | 0.29 / 0.19 | 0.00 | 2.3 % |
| sonde-relecture-naïve | 33 % [25 ; 43] (100) | 186 [148 ; 193] | +77 [58 ; 94] (32) | 0.06 / 0.48 | 0.66 | 8.9 % |
| sonde-relecture | 86 % [78 ; 91] (100) | 139 [136 ; 144] | +48 [46 ; 52] (84) | 0.01 / 0.26 | 0.82 | 0.5 % |
| sonde-par-composante | 0 % [0 ; 4] (100) | – [– ; –] | – | 0.00 / 0.94 | 0.00 | 61.1 % |
| sonde-relecture-par-composante | 85 % [77 ; 91] (100) | 137 [135 ; 140] | +47 [45 ; 50] (82) | 0.00 / 0.27 | 0.82 | 0.5 % |
| sonde-oracle-décade (triche) | 95 % [89 ; 98] (100) | 89 [88 ; 92] | réf. | 0.00 / 0.25 | 0.00 | 0.0 % |
| Broyden-sans-couleur | 100 % [96 ; 100] (100) | 152 [150 ; 153] | +65 [63 ; 69] (95) | – | – | – |

### n = 32, σ = 0.02

| méthode | succès [IC 95 %] | lectures médianes [IC 95 %] | écart apparié à l'oracle [IC 95 %] (paires) | rejets ambiguïté / descente par lecture de pas | relectures par pas | pas acceptés avec une décade fausse muette |
|---|---|---|---|---|---|---|
| sonde-garde-ambiguïté | 1 % [0 ; 5] (100) | 162 [– ; –] | – | 0.20 / 0.79 | 0.00 | 1.5 % |
| sonde-dépliage-vers-zéro | 27 % [19 ; 36] (100) | 172 [165 ; 183] | +46 [38 ; 57] (26) | 0.08 / 0.86 | 0.00 | 2.6 % |
| sonde-relecture-naïve | 18 % [12 ; 27] (100) | 222 [196 ; 248] | +95 [65 ; 122] (15) | 0.00 / 0.87 | 0.37 | 11.2 % |
| sonde-relecture | 59 % [49 ; 68] (100) | 219 [202 ; 232] | +88 [72 ; 102] (54) | 0.01 / 0.79 | 0.44 | 1.2 % |
| sonde-oracle-décade (triche) | 83 % [74 ; 89] (100) | 129 [126 ; 131] | réf. | 0.00 / 0.79 | 0.00 | 0.0 % |
| Broyden-sans-couleur | 100 % [96 ; 100] (100) | 152 [150 ; 154] | +27 [23 ; 31] (83) | – | – | – |
