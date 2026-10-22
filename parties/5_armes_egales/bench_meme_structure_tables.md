# Partie 5.5 : même structure pour tous, tableaux

Succès [IC 95 %] / lectures médianes. Broyden-quaternionique-mêmes-mesures : meilleur des deux pas de différence finie (0,5 et 1) par case.

## Problème ℍ-linéaire des parties 3.1–3.5 (structure à gauche), 10 instances × 10 départs, départs à 1,5·√n, σ = 0,00

| structure / activation | taille | dimension | sonde de la partie d'origine (quaternionique / à droite / par sortie) | même sonde, oracle-décade | sonde-deux-échelles + test géométrique | Broyden-quaternionique-mêmes-mesures | Broyden-quaternionique-sans-couleur | Broyden-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| ℍ-linéaire (3.2) | n = 4 | 16 | 99 % [95 ; 100] / 25 | 100 % [96 ; 100] / 23 | 100 % [96 ; 100] / 26 | 100 % [96 ; 100] / 18 | 100 % [96 ; 100] / 14 | 100 % [96 ; 100] / 29 |
| ℍ-linéaire (3.2) | n = 8 | 32 | 89 % [81 ; 94] / 48 | 100 % [96 ; 100] / 34 | 99 % [95 ; 100] / 42 | 100 % [96 ; 100] / 29 | 100 % [96 ; 100] / 23 | 100 % [96 ; 100] / 49 |
| ℍ-linéaire (3.2) | n = 16 | 64 | 82 % [73 ; 88] / 100 | 100 % [96 ; 100] / 53 | 97 % [92 ; 99] / 71 | 100 % [96 ; 100] / 45 | 100 % [96 ; 100] / 34 | 100 % [96 ; 100] / 84 |
| ℍ-linéaire (3.2) | n = 32 | 128 | 86 % [78 ; 91] / 172 | 95 % [89 ; 98] / 89 | 99 % [95 ; 100] / 122 | 100 % [96 ; 100] / 70 | 100 % [96 ; 100] / 53 | 100 % [96 ; 100] / 152 |

## Problème ℍ-linéaire des parties 3.1–3.5 (structure à gauche), 10 instances × 10 départs, départs à 1,5·√n, σ = 0,02

| structure / activation | taille | dimension | sonde de la partie d'origine (quaternionique / à droite / par sortie) | même sonde, oracle-décade | sonde-deux-échelles + test géométrique | Broyden-quaternionique-mêmes-mesures | Broyden-quaternionique-sans-couleur | Broyden-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| ℍ-linéaire (3.2) | n = 4 | 16 | 99 % [95 ; 100] / 33 | 99 % [95 ; 100] / 26 | 100 % [96 ; 100] / 29 | 100 % [96 ; 100] / 20 | 100 % [96 ; 100] / 14 | 100 % [96 ; 100] / 29 |
| ℍ-linéaire (3.2) | n = 8 | 32 | 87 % [79 ; 92] / 72 | 100 % [96 ; 100] / 42 | 100 % [96 ; 100] / 50 | 100 % [96 ; 100] / 30 | 100 % [96 ; 100] / 23 | 100 % [96 ; 100] / 49 |
| ℍ-linéaire (3.2) | n = 16 | 64 | 72 % [63 ; 80] / 122 | 100 % [96 ; 100] / 77 | 100 % [96 ; 100] / 91 | 100 % [96 ; 100] / 47 | 100 % [96 ; 100] / 34 | 100 % [96 ; 100] / 84 |
| ℍ-linéaire (3.2) | n = 32 | 128 | 1 % [0 ; 5] / 162 | 83 % [74 ; 89] / 129 | 44 % [35 ; 54] / 150 | 100 % [96 ; 100] / 74 | 100 % [96 ; 100] / 53 | 100 % [96 ; 100] / 152 |

## Entraîner un neurone (partie 6.3, structure à droite), 3 instances × 20 départs, σ = 0,00

| structure / activation | taille | dimension | sonde de la partie d'origine (quaternionique / à droite / par sortie) | même sonde, oracle-décade | sonde-deux-échelles + test géométrique | Broyden-quaternionique-mêmes-mesures | Broyden-quaternionique-sans-couleur | Broyden-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| identité | n = 4 | 16 | 100 % [94 ; 100] / 12 | 100 % [94 ; 100] / 12 | — | 100 % [94 ; 100] / 8 | 100 % [94 ; 100] / 6 | 100 % [94 ; 100] / 18 |
| identité | n = 8 | 32 | 100 % [94 ; 100] / 17 | 100 % [94 ; 100] / 17 | — | 100 % [94 ; 100] / 13 | 100 % [94 ; 100] / 10 | 100 % [94 ; 100] / 34 |
| identité | n = 16 | 64 | 90 % [80 ; 95] / 33 | 100 % [94 ; 100] / 29 | — | 100 % [94 ; 100] / 23 | 100 % [94 ; 100] / 18 | 100 % [94 ; 100] / 66 |
| tanh radial, petit signal | n = 4 | 16 | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 11 | — | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 8 | 100 % [94 ; 100] / 20 |
| tanh radial, petit signal | n = 8 | 32 | 98 % [91 ; 100] / 16 | 100 % [94 ; 100] / 16 | — | 100 % [94 ; 100] / 13 | 100 % [94 ; 100] / 13 | 100 % [94 ; 100] / 36 |
| tanh radial, petit signal | n = 16 | 64 | 93 % [84 ; 97] / 34 | 100 % [94 ; 100] / 28 | — | 100 % [94 ; 100] / 23 | 100 % [94 ; 100] / 22 | 100 % [94 ; 100] / 68 |
| tanh radial | n = 4 | 16 | 60 % [47 ; 71] / 73 | 33 % [23 ; 46] / 25 | — | 100 % [94 ; 100] / 30 | 100 % [94 ; 100] / 30 | 100 % [94 ; 100] / 26 |
| tanh radial | n = 8 | 32 | 3 % [1 ; 11] / 36 | 7 % [3 ; 16] / 37 | — | 97 % [89 ; 99] / 70 | 100 % [94 ; 100] / 72 | 100 % [94 ; 100] / 46 |
| tanh radial | n = 16 | 64 | — | 0 % [0 ; 6] | — | 77 % [65 ; 86] / 132 | 93 % [84 ; 97] / 156 | — |

## Entraîner un neurone (partie 6.3, structure à droite), 3 instances × 20 départs, σ = 0,02

| structure / activation | taille | dimension | sonde de la partie d'origine (quaternionique / à droite / par sortie) | même sonde, oracle-décade | sonde-deux-échelles + test géométrique | Broyden-quaternionique-mêmes-mesures | Broyden-quaternionique-sans-couleur | Broyden-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| identité | n = 4 | 16 | 100 % [94 ; 100] / 13 | 100 % [94 ; 100] / 13 | — | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 6 | 100 % [94 ; 100] / 18 |
| identité | n = 8 | 32 | 93 % [84 ; 97] / 19 | 100 % [94 ; 100] / 18 | — | 100 % [94 ; 100] / 14 | 100 % [94 ; 100] / 10 | 100 % [94 ; 100] / 34 |
| identité | n = 16 | 64 | 65 % [52 ; 76] / 49 | 100 % [94 ; 100] / 34 | — | 100 % [94 ; 100] / 24 | 100 % [94 ; 100] / 18 | 100 % [94 ; 100] / 66 |
| tanh radial, petit signal | n = 4 | 16 | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 11 | — | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 8 | 100 % [94 ; 100] / 20 |
| tanh radial, petit signal | n = 8 | 32 | 97 % [89 ; 99] / 18 | 100 % [94 ; 100] / 17 | — | 100 % [94 ; 100] / 14 | 100 % [94 ; 100] / 13 | 100 % [94 ; 100] / 36 |
| tanh radial, petit signal | n = 16 | 64 | 72 % [59 ; 81] / 48 | 100 % [94 ; 100] / 33 | — | 100 % [94 ; 100] / 24 | 100 % [94 ; 100] / 22 | 100 % [94 ; 100] / 68 |
| tanh radial | n = 4 | 16 | 50 % [38 ; 62] / 96 | 47 % [35 ; 59] / 25 | — | 100 % [94 ; 100] / 32 | 100 % [94 ; 100] / 30 | 100 % [94 ; 100] / 26 |
| tanh radial | n = 8 | 32 | 2 % [0 ; 9] / 67 | 8 % [4 ; 18] / 32 | — | 75 % [63 ; 84] / 46 | 100 % [94 ; 100] / 72 | 100 % [94 ; 100] / 46 |
| tanh radial | n = 16 | 64 | — | 2 % [0 ; 9] / 98 | — | 42 % [30 ; 54] / 82 | 93 % [84 ; 97] / 156 | — |

## Entraîner une couche (partie 6.4, à droite + par sortie), 3 instances × 10 départs, σ = 0,00

| structure / activation | taille | dimension | sonde de la partie d'origine (quaternionique / à droite / par sortie) | même sonde, oracle-décade | sonde-deux-échelles + test géométrique | Broyden-quaternionique-mêmes-mesures | Broyden-quaternionique-sans-couleur | Broyden-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| identité | 4 × 4 | 64 | 100 % [89 ; 100] / 13 | 100 % [89 ; 100] / 13 | — | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 6 | 100 % [89 ; 100] / 66 |
| identité | 4 × 8 | 128 | 97 % [83 ; 99] / 18 | 100 % [89 ; 100] / 18 | — | 100 % [89 ; 100] / 13 | 100 % [89 ; 100] / 10 | 100 % [89 ; 100] / 130 |
| identité | 8 × 4 | 128 | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 14 | — | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 6 | 100 % [89 ; 100] / 130 |
| identité | 8 × 8 | 256 | 100 % [89 ; 100] / 18 | 100 % [89 ; 100] / 18 | — | 100 % [89 ; 100] / 13 | 100 % [89 ; 100] / 10 | 100 % [89 ; 100] / 258 |
| tanh radial, petit signal | 4 × 4 | 64 | 100 % [89 ; 100] / 12 | 100 % [89 ; 100] / 12 | — | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 68 |
| tanh radial, petit signal | 4 × 8 | 128 | 100 % [89 ; 100] / 17 | 100 % [89 ; 100] / 17 | — | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 13 | 100 % [89 ; 100] / 133 |
| tanh radial, petit signal | 8 × 4 | 128 | 100 % [89 ; 100] / 12 | 100 % [89 ; 100] / 12 | — | 100 % [89 ; 100] / 10 | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 133 |
| tanh radial, petit signal | 8 × 8 | 256 | 100 % [89 ; 100] / 18 | 100 % [89 ; 100] / 17 | — | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 261 |
| tanh radial | 4 × 4 | 64 | 0 % [0 ; 11] | — | — | 87 % [70 ; 95] / 48 | 93 % [79 ; 98] / 50 | 97 % [83 ; 99] / 80 |
| tanh radial | 4 × 8 | 128 | 0 % [0 ; 11] | — | — | 67 % [49 ; 81] / 88 | 83 % [66 ; 93] / 103 | 97 % [83 ; 99] / 164 |
| tanh radial | 8 × 4 | 128 | 0 % [0 ; 11] | — | — | 73 % [56 ; 86] / 68 | 90 % [74 ; 97] / 109 | 100 % [89 ; 100] / 152 |
| tanh radial | 8 × 8 | 256 | 0 % [0 ; 11] | — | — | 60 % [42 ; 75] / 143 | 97 % [83 ; 99] / 215 | 100 % [89 ; 100] / 286 |

## Entraîner une couche (partie 6.4, à droite + par sortie), 3 instances × 10 départs, σ = 0,02

| structure / activation | taille | dimension | sonde de la partie d'origine (quaternionique / à droite / par sortie) | même sonde, oracle-décade | sonde-deux-échelles + test géométrique | Broyden-quaternionique-mêmes-mesures | Broyden-quaternionique-sans-couleur | Broyden-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| identité | 4 × 4 | 64 | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 14 | — | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 6 | 100 % [89 ; 100] / 66 |
| identité | 4 × 8 | 128 | 93 % [79 ; 98] / 23 | 93 % [79 ; 98] / 20 | — | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 10 | 100 % [89 ; 100] / 130 |
| identité | 8 × 4 | 128 | 97 % [83 ; 99] / 15 | 100 % [89 ; 100] / 14 | — | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 6 | 100 % [89 ; 100] / 130 |
| identité | 8 × 8 | 256 | 87 % [70 ; 95] / 26 | 77 % [59 ; 88] / 20 | — | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 10 | 100 % [89 ; 100] / 258 |
| tanh radial, petit signal | 4 × 4 | 64 | 97 % [83 ; 99] / 12 | 100 % [89 ; 100] / 12 | — | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 68 |
| tanh radial, petit signal | 4 × 8 | 128 | 93 % [79 ; 98] / 24 | 100 % [89 ; 100] / 18 | — | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 13 | 100 % [89 ; 100] / 133 |
| tanh radial, petit signal | 8 × 4 | 128 | 97 % [83 ; 99] / 14 | 100 % [89 ; 100] / 13 | — | 100 % [89 ; 100] / 10 | 100 % [89 ; 100] / 9 | 100 % [89 ; 100] / 133 |
| tanh radial, petit signal | 8 × 8 | 256 | 97 % [83 ; 99] / 24 | 83 % [66 ; 93] / 19 | — | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 14 | 100 % [89 ; 100] / 261 |
| tanh radial | 4 × 4 | 64 | 0 % [0 ; 11] | — | — | 77 % [59 ; 88] / 40 | 93 % [79 ; 98] / 50 | 97 % [83 ; 99] / 80 |
| tanh radial | 4 × 8 | 128 | 0 % [0 ; 11] | — | — | 50 % [33 ; 67] / 80 | 83 % [66 ; 93] / 103 | 97 % [83 ; 99] / 164 |
| tanh radial | 8 × 4 | 128 | 0 % [0 ; 11] | — | — | 70 % [52 ; 83] / 89 | 90 % [74 ; 97] / 109 | 100 % [89 ; 100] / 152 |
| tanh radial | 8 × 8 | 256 | 0 % [0 ; 11] | — | — | 30 % [17 ; 48] / 91 | 97 % [83 ; 99] / 215 | 100 % [89 ; 100] / 286 |

## Gains appariés (+ : la première méthode est moins chère ; IC 95 % bootstrap ; départs réussis par les deux)

| famille | activation | taille | σ | comparaison | gain % |
|---|---|---|---|---|---|
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +40 [+38 ; +44] (99 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +22 [+18 ; +22] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +50 [+48 ; +52] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +52 [+50 ; +57] (89 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +15 [+12 ; +18] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +53 [+52 ; +55] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +66 [+64 ; +70] (82 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 16 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +15 [+13 ; +17] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +60 [+59 ; +60] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 32 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +69 [+68 ; +70] (86 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 32 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +22 [+20 ; +25] (95 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 32 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +65 [+64 ; +65] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +55 [+50 ; +60] (99 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +23 [+19 ; +23] (99 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +50 [+48 ; +52] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +68 [+66 ; +71] (87 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +29 [+26 ; +31] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +53 [+52 ; +55] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +72 [+70 ; +74] (72 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 16 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +39 [+36 ; +41] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +60 [+59 ; +60] (100 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 32 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 32 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +42 [+40 ; +43] (83 paires) |
| ℍ-linéaire | ℍ-linéaire (3.2) | n = 32 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +65 [+64 ; +65] (100 paires) |
| neurone | identité | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +52 [+50 ; +54] (60 paires) |
| neurone | identité | n = 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +33 [+33 ; +38] (60 paires) |
| neurone | identité | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +67 [+67 ; +67] (60 paires) |
| neurone | identité | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +41 [+41 ; +41] (60 paires) |
| neurone | identité | n = 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +24 [+24 ; +24] (60 paires) |
| neurone | identité | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +71 [+71 ; +71] (60 paires) |
| neurone | identité | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +45 [+44 ; +47] (54 paires) |
| neurone | identité | n = 16 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +21 [+18 ; +21] (60 paires) |
| neurone | identité | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +73 [+73 ; +73] (60 paires) |
| neurone | tanh radial, petit signal | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +27 [+27 ; +27] (60 paires) |
| neurone | tanh radial, petit signal | n = 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +18 [+18 ; +18] (60 paires) |
| neurone | tanh radial, petit signal | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +60 [+60 ; +60] (60 paires) |
| neurone | tanh radial, petit signal | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +19 [+19 ; +24] (59 paires) |
| neurone | tanh radial, petit signal | n = 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +19 [+19 ; +19] (60 paires) |
| neurone | tanh radial, petit signal | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +64 [+64 ; +64] (60 paires) |
| neurone | tanh radial, petit signal | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +34 [+33 ; +35] (56 paires) |
| neurone | tanh radial, petit signal | n = 16 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +18 [+18 ; +19] (60 paires) |
| neurone | tanh radial, petit signal | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +68 [+68 ; +68] (60 paires) |
| neurone | tanh radial | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +62 [+54 ; +69] (36 paires) |
| neurone | tanh radial | n = 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | -4 [-17 ; +10] (20 paires) |
| neurone | tanh radial | n = 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | -17 [-42 ; -8] (60 paires) |
| neurone | tanh radial | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| neurone | tanh radial | n = 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +12 [-31 ; +14] (4 paires) |
| neurone | tanh radial | n = 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | -57 [-113 ; -18] (60 paires) |
| neurone | tanh radial | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| neurone | tanh radial | n = 16 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| neurone | tanh radial | n = 16 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | — |
| neurone | identité | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +54 [+54 ; +54] (60 paires) |
| neurone | identité | n = 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +31 [+31 ; +31] (60 paires) |
| neurone | identité | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +67 [+67 ; +67] (60 paires) |
| neurone | identité | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +47 [+47 ; +50] (56 paires) |
| neurone | identité | n = 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +22 [+22 ; +26] (60 paires) |
| neurone | identité | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +71 [+71 ; +71] (60 paires) |
| neurone | identité | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +63 [+62 ; +65] (39 paires) |
| neurone | identité | n = 16 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +29 [+27 ; +31] (60 paires) |
| neurone | identité | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +73 [+73 ; +73] (60 paires) |
| neurone | tanh radial, petit signal | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +27 [+27 ; +33] (60 paires) |
| neurone | tanh radial, petit signal | n = 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +18 [+18 ; +25] (60 paires) |
| neurone | tanh radial, petit signal | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +60 [+60 ; +60] (60 paires) |
| neurone | tanh radial, petit signal | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +30 [+26 ; +32] (58 paires) |
| neurone | tanh radial, petit signal | n = 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +18 [+18 ; +22] (60 paires) |
| neurone | tanh radial, petit signal | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +64 [+64 ; +64] (60 paires) |
| neurone | tanh radial, petit signal | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +54 [+52 ; +57] (43 paires) |
| neurone | tanh radial, petit signal | n = 16 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +26 [+22 ; +27] (60 paires) |
| neurone | tanh radial, petit signal | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +68 [+68 ; +68] (60 paires) |
| neurone | tanh radial | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +71 [+62 ; +75] (30 paires) |
| neurone | tanh radial | n = 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | -8 [-22 ; +18] (28 paires) |
| neurone | tanh radial | n = 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | -17 [-42 ; -8] (60 paires) |
| neurone | tanh radial | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| neurone | tanh radial | n = 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | -12 [-34 ; +57] (4 paires) |
| neurone | tanh radial | n = 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | -57 [-113 ; -18] (60 paires) |
| neurone | tanh radial | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| neurone | tanh radial | n = 16 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| neurone | tanh radial | n = 16 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | — |
| couche | identité | 4 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +54 [+54 ; +54] (30 paires) |
| couche | identité | 4 × 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +31 [+31 ; +31] (30 paires) |
| couche | identité | 4 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +91 [+91 ; +91] (30 paires) |
| couche | identité | 4 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +44 [+44 ; +44] (29 paires) |
| couche | identité | 4 × 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +28 [+28 ; +28] (30 paires) |
| couche | identité | 4 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +92 [+92 ; +92] (30 paires) |
| couche | identité | 8 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +57 [+57 ; +57] (30 paires) |
| couche | identité | 8 × 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +36 [+33 ; +36] (30 paires) |
| couche | identité | 8 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +95 [+95 ; +95] (30 paires) |
| couche | identité | 8 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +44 [+44 ; +46] (30 paires) |
| couche | identité | 8 × 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +28 [+22 ; +28] (30 paires) |
| couche | identité | 8 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +96 [+96 ; +96] (30 paires) |
| couche | tanh radial, petit signal | 4 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +25 [+25 ; +25] (30 paires) |
| couche | tanh radial, petit signal | 4 × 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +25 [+25 ; +25] (30 paires) |
| couche | tanh radial, petit signal | 4 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +87 [+87 ; +87] (30 paires) |
| couche | tanh radial, petit signal | 4 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +24 [+18 ; +28] (30 paires) |
| couche | tanh radial, petit signal | 4 × 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +18 [+12 ; +18] (30 paires) |
| couche | tanh radial, petit signal | 4 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +90 [+89 ; +90] (30 paires) |
| couche | tanh radial, petit signal | 8 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +25 [+25 ; +31] (30 paires) |
| couche | tanh radial, petit signal | 8 × 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +17 [+17 ; +25] (30 paires) |
| couche | tanh radial, petit signal | 8 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +93 [+93 ; +93] (30 paires) |
| couche | tanh radial, petit signal | 8 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +22 [+18 ; +28] (30 paires) |
| couche | tanh radial, petit signal | 8 × 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +18 [+18 ; +18] (30 paires) |
| couche | tanh radial, petit signal | 8 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +95 [+95 ; +95] (30 paires) |
| couche | tanh radial | 4 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 4 × 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 4 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +40 [+10 ; +54] (27 paires) |
| couche | tanh radial | 4 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 4 × 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 4 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +35 [-40 ; +51] (25 paires) |
| couche | tanh radial | 8 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 8 × 4 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 8 × 4 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +28 [+14 ; +59] (27 paires) |
| couche | tanh radial | 8 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 8 × 8 | 0,00 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 8 × 8 | 0,00 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +25 [-39 ; +62] (29 paires) |
| couche | identité | 4 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +57 [+56 ; +57] (30 paires) |
| couche | identité | 4 × 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +36 [+31 ; +36] (30 paires) |
| couche | identité | 4 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +91 [+91 ; +91] (30 paires) |
| couche | identité | 4 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +57 [+53 ; +58] (28 paires) |
| couche | identité | 4 × 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +28 [+26 ; +30] (28 paires) |
| couche | identité | 4 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +92 [+92 ; +92] (30 paires) |
| couche | identité | 8 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +60 [+57 ; +62] (29 paires) |
| couche | identité | 8 × 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +36 [+29 ; +36] (30 paires) |
| couche | identité | 8 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +95 [+95 ; +95] (30 paires) |
| couche | identité | 8 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +62 [+58 ; +63] (26 paires) |
| couche | identité | 8 × 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +30 [+25 ; +30] (23 paires) |
| couche | identité | 8 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +96 [+96 ; +96] (30 paires) |
| couche | tanh radial, petit signal | 4 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +25 [+25 ; +31] (29 paires) |
| couche | tanh radial, petit signal | 4 × 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +25 [+17 ; +28] (30 paires) |
| couche | tanh radial, petit signal | 4 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +87 [+87 ; +87] (30 paires) |
| couche | tanh radial, petit signal | 4 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +46 [+42 ; +48] (28 paires) |
| couche | tanh radial, petit signal | 4 × 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +22 [+17 ; +26] (30 paires) |
| couche | tanh radial, petit signal | 4 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +90 [+89 ; +90] (30 paires) |
| couche | tanh radial, petit signal | 8 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +36 [+31 ; +40] (29 paires) |
| couche | tanh radial, petit signal | 8 × 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +23 [+19 ; +23] (30 paires) |
| couche | tanh radial, petit signal | 8 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +93 [+93 ; +93] (30 paires) |
| couche | tanh radial, petit signal | 8 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | +42 [+39 ; +48] (29 paires) |
| couche | tanh radial, petit signal | 8 × 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | +26 [+21 ; +26] (25 paires) |
| couche | tanh radial, petit signal | 8 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +95 [+95 ; +95] (30 paires) |
| couche | tanh radial | 4 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 4 × 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 4 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +40 [+10 ; +54] (27 paires) |
| couche | tanh radial | 4 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 4 × 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 4 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +35 [-40 ; +51] (25 paires) |
| couche | tanh radial | 8 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 8 × 4 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 8 × 4 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +28 [+14 ; +59] (27 paires) |
| couche | tanh radial | 8 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre sonde de la partie d'origine (quaternionique / à droite / par sortie) | — |
| couche | tanh radial | 8 × 8 | 0,02 | Broyden-quaternionique-mêmes-mesures contre même sonde, oracle-décade | — |
| couche | tanh radial | 8 × 8 | 0,02 | Broyden-quaternionique-sans-couleur contre Broyden-sans-couleur | +25 [-39 ; +62] (29 paires) |
