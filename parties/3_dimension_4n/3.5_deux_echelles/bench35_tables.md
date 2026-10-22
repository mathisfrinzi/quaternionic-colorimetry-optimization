# Partie 3.5 : tableaux

Succès [IC 95 % de Wilson] / lectures médianes (succès seulement), 100 départs par case.

| méthode | n = 16, σ = 0,0 | n = 16, σ = 0,02 | n = 32, σ = 0,0 | n = 32, σ = 0,02 |
|---|---|---|---|---|
| sonde-garde-ambiguïté | 82 % [73 ; 88] / 100 | 72 % [63 ; 80] / 122 | 86 % [78 ; 91] / 172 | 1 % [0 ; 5] / 162 |
| sonde-dépliage-vers-zéro | 98 % [93 ; 99] / 63 | 91 % [84 ; 95] / 108 | 92 % [85 ; 96] / 113 | 27 % [19 ; 36] / 172 |
| sonde-relecture | 99 % [95 ; 100] / 87 | 95 % [89 ; 98] / 127 | 86 % [78 ; 91] / 139 | 59 % [49 ; 68] / 219 |
| sonde-relecture-économe | 99 % [95 ; 100] / 84 | 91 % [84 ; 95] / 119 | — | — |
| sonde-vote | 22 % [15 ; 31] / 82 | 5 % [2 ; 11] / 77 | — | — |
| sonde-deux-échelles | 97 % [92 ; 99] / 75 | 98 % [93 ; 99] / 89 | 96 % [90 ; 98] / 130 | 25 % [18 ; 34] / 152 |
| sonde-deux-échelles + test géométrique | 97 % [92 ; 99] / 71 | 100 % [96 ; 100] / 91 | 99 % [95 ; 100] / 122 | 44 % [35 ; 54] / 150 |
| sonde-deux-échelles + test pondéré | 100 % [96 ; 100] / 72 | 99 % [95 ; 100] / 84 | 100 % [96 ; 100] / 125 | — |
| sonde-oracle-décade | 100 % [96 ; 100] / 53 | 100 % [96 ; 100] / 77 | 95 % [89 ; 98] / 89 | 83 % [74 ; 89] / 129 |
| sonde-oracle-décade + test géométrique | 100 % [96 ; 100] / 52 | 100 % [96 ; 100] / 76 | 92 % [85 ; 96] / 86 | — |
| sonde-oracle-décade + test pondéré | 100 % [96 ; 100] / 52 | 100 % [96 ; 100] / 75 | — | — |
| Broyden-sans-couleur | 100 % [96 ; 100] / 84 | 100 % [96 ; 100] / 84 | 100 % [96 ; 100] / 152 | 100 % [96 ; 100] / 152 |

## Test de descente : refus et refus à tort

Part des lectures de pas refusées pour ambiguïté / par le test de descente ; entre parenthèses, part des refus du test de descente alors que la vraie moyenne des ‖ρᵢ‖² passait.

| méthode | n = 16, σ = 0,0 | n = 16, σ = 0,02 | n = 32, σ = 0,0 | n = 32, σ = 0,02 |
|---|---|---|---|---|
| sonde-relecture-économe | 0,01 / 0,12 (0,03) | 0,00 / 0,52 (0,02) | — | — |
| sonde-vote | 0,27 / 0,37 (0,00) | 0,09 / 0,82 (0,00) | — | — |
| sonde-deux-échelles | 0,00 / 0,44 (0,09) | 0,00 / 0,56 (0,03) | 0,00 / 0,57 (0,04) | 0,00 / 0,97 (0,00) |
| sonde-deux-échelles + test géométrique | 0,00 / 0,34 (0,09) | 0,00 / 0,46 (0,01) | 0,00 / 0,39 (0,02) | 0,00 / 0,92 (0,00) |
| sonde-deux-échelles + test pondéré | 0,00 / 0,22 (0,13) | 0,00 / 0,49 (0,08) | 0,00 / 0,38 (0,06) | — |
| sonde-oracle-décade | 0,00 / 0,11 (0,04) | 0,00 / 0,39 (0,03) | 0,00 / 0,25 (0,02) | 0,00 / 0,79 (0,01) |
| sonde-oracle-décade + test géométrique | 0,00 / 0,07 (0,00) | 0,00 / 0,30 (0,00) | 0,00 / 0,15 (0,00) | — |
| sonde-oracle-décade + test pondéré | 0,00 / 0,06 (0,15) | 0,00 / 0,30 (0,08) | — | — |

## Gains appariés (+ : la première est moins chère ; IC 95 % bootstrap)

| comparaison | n = 16, σ = 0,0 | n = 16, σ = 0,02 | n = 32, σ = 0,0 | n = 32, σ = 0,02 |
|---|---|---|---|---|
| sonde-deux-échelles contre sonde-relecture | +14 [+11 ; +16] (96) | +30 [+27 ; +35] (93) | +8 [+0 ; +13] (83) | +22 [+13 ; +36] (18) |
| sonde-deux-échelles contre sonde-dépliage-vers-zéro | -19 [-24 ; -15] (96) | +18 [+12 ; +23] (90) | -17 [-24 ; -4] (88) | +12 [-3 ; +19] (8) |
| sonde-deux-échelles + test géométrique contre sonde-deux-échelles | +5 [+2 ; +9] (94) | -2 [-7 ; +5] (98) | +4 [+0 ; +11] (95) | +1 [-8 ; +6] (17) |
| sonde-deux-échelles + test pondéré contre sonde-deux-échelles | +4 [+0 ; +7] (97) | +6 [+1 ; +9] (98) | +3 [-2 ; +9] (96) | — |
| sonde-deux-échelles + test géométrique contre Broyden-sans-couleur | +15 [+12 ; +18] (97) | -8 [-12 ; -1] (100) | +19 [+17 ; +22] (99) | +1 [-2 ; +5] (44) |
| sonde-relecture-économe contre sonde-relecture | +3 [+2 ; +6] (99) | +6 [-1 ; +11] (86) | — | — |
| sonde-oracle-décade + test géométrique contre sonde-oracle-décade | +1 [+0 ; +3] (100) | +1 [-3 ; +6] (100) | +4 [+2 ; +6] (88) | — |
