# Partie 6.5 : pourquoi la sonde s'effondre en régime saturé, tableaux

Succès [IC 95 %] / lectures médianes ; décade donnée (oracle), σ = 0, 3 instances × 20 départs.

## neurone tanh radial (défaut 0,27)

| variante | n = 4 | n = 8 |
|---|---|---|
| sonde-oracle-décade (réglage de la partie) | 33 % [23 ; 46] / 25 | 7 % [3 ; 16] / 37 |
| + test pondéré | 65 % [52 ; 76] / 30 | 22 % [13 ; 34] / 51 |
| + test géométrique | 42 % [30 ; 54] / 24 | 8 % [4 ; 18] / 37 |
| + pas de Newton entier (1,0) | 53 % [41 ; 65] / 26 | 25 % [16 ; 37] / 33 |
| + seuil d'acceptation 1,0 | 42 % [30 ; 54] / 26 | 12 % [6 ; 22] / 38 |
| + test pondéré + pas 1,0 | 70 % [57 ; 80] / 28 | 33 % [23 ; 46] / 32 |
| + test pondéré + pas 1,0 + seuil 1,0 | 63 % [51 ; 74] / 23 | 38 % [27 ; 51] / 35 |
| + coupure des valeurs singulières à 10⁻¹⁰ (écartée) | 55 % [42 ; 67] / 29 | 12 % [6 ; 22] / 42 |
| + taille des sondes fixe après échec (écartée) | 30 % [20 ; 43] / 24 | 7 % [3 ; 16] / 37 |
| + test pondéré + pas 1,0 + les deux précédentes | 80 % [68 ; 88] / 28 | 35 % [24 ; 48] / 33 |
| aucun test de descente (pas 1,0) | 82 % [70 ; 89] / 25 | 77 % [65 ; 86] / 36 |
| aucun test de descente (pas 0,65) | 72 % [59 ; 81] / 26 | 48 % [36 ; 61] / 39 |
| Broyden-quaternionique sans recherche linéaire (aucun test) | 100 % [94 ; 100] / 26 | 100 % [94 ; 100] / 43 |
| Broyden-quaternionique-sans-couleur | 100 % [94 ; 100] / 30 | 100 % [94 ; 100] / 72 |

## neurone tanh radial petit signal (défaut 0,05)

| variante | n = 4 | n = 8 |
|---|---|---|
| sonde-oracle-décade (réglage de la partie) | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 16 |
| + test pondéré | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 16 |
| + test géométrique | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 16 |
| + pas de Newton entier (1,0) | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 15 |
| + seuil d'acceptation 1,0 | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 16 |
| + test pondéré + pas 1,0 | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 15 |
| + test pondéré + pas 1,0 + seuil 1,0 | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 15 |
| + coupure des valeurs singulières à 10⁻¹⁰ (écartée) | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 16 |
| + taille des sondes fixe après échec (écartée) | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 16 |
| + test pondéré + pas 1,0 + les deux précédentes | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 15 |
| aucun test de descente (pas 1,0) | 100 % [94 ; 100] / 9 | 100 % [94 ; 100] / 15 |
| aucun test de descente (pas 0,65) | 100 % [94 ; 100] / 11 | 100 % [94 ; 100] / 16 |
| Broyden-quaternionique sans recherche linéaire (aucun test) | 100 % [94 ; 100] / 8 | 100 % [94 ; 100] / 13 |
| Broyden-quaternionique-sans-couleur | 100 % [94 ; 100] / 8 | 100 % [94 ; 100] / 13 |
