# Partie 1.4 — tableaux

## Succès % et lectures médianes (zéros isolés + Möbius ε ∈ {0, 1}, 3 instances × 30 départs)

| bits | σ | méthode | succès % | lectures méd. |
|---|---|---|---|---|
| 24 | 0.0 | sonde-couleur-mémoire | 97 | 17 |
| 24 | 0.0 | sonde-couleur-mémoire + dithering | 97 | 17 |
| 24 | 0.0 | sonde-couleur-mémoire + démarrage-conforme | 89 | 20 |
| 24 | 0.0 | sonde-couleur-mémoire + région-de-confiance | 90 | 20 |
| 24 | 0.0 | sonde-couleur-mémoire + les trois | 80 | 20 |
| 24 | 0.0 | Broyden-mêmes-mesures | 92 | 14 |
| 24 | 0.0 | Broyden-mêmes-mesures + dithering | 90 | 97 |
| 24 | 0.02 | sonde-couleur-mémoire | 94 | 18 |
| 24 | 0.02 | sonde-couleur-mémoire + dithering | 98 | 17 |
| 24 | 0.02 | sonde-couleur-mémoire + démarrage-conforme | 90 | 21 |
| 24 | 0.02 | sonde-couleur-mémoire + région-de-confiance | 93 | 20 |
| 24 | 0.02 | sonde-couleur-mémoire + les trois | 77 | 18 |
| 24 | 0.02 | Broyden-mêmes-mesures | 94 | 28 |
| 24 | 0.02 | Broyden-mêmes-mesures + dithering | 94 | 15 |
| 36 | 0.0 | sonde-couleur-mémoire | 97 | 15 |
| 36 | 0.0 | sonde-couleur-mémoire + dithering | 81 | 13 |
| 36 | 0.0 | sonde-couleur-mémoire + démarrage-conforme | 96 | 20 |
| 36 | 0.0 | sonde-couleur-mémoire + région-de-confiance | 94 | 19 |
| 36 | 0.0 | sonde-couleur-mémoire + les trois | 84 | 16 |
| 36 | 0.0 | Broyden-mêmes-mesures | 94 | 12 |
| 36 | 0.0 | Broyden-mêmes-mesures + dithering | 92 | 114 |
| 36 | 0.02 | sonde-couleur-mémoire | 95 | 20 |
| 36 | 0.02 | sonde-couleur-mémoire + dithering | 94 | 17 |
| 36 | 0.02 | sonde-couleur-mémoire + démarrage-conforme | 89 | 20 |
| 36 | 0.02 | sonde-couleur-mémoire + région-de-confiance | 89 | 16 |
| 36 | 0.02 | sonde-couleur-mémoire + les trois | 73 | 16 |
| 36 | 0.02 | Broyden-mêmes-mesures | 94 | 14 |
| 36 | 0.02 | Broyden-mêmes-mesures + dithering | 94 | 15 |


## Gains appariés en % de lectures (+ = la 1ʳᵉ méthode nommée est moins chère), IC 95 % bootstrap

| bits | σ | comparaison | gain % | n |
|---|---|---|---|---|
| 24 | 0.0 | sonde-couleur-mémoire + dithering contre sonde-couleur-mémoire | +10 [+6, +15] | 512 |
| 24 | 0.0 | sonde-couleur-mémoire + démarrage-conforme contre sonde-couleur-mémoire | -27 [-36, -19] | 470 |
| 24 | 0.0 | sonde-couleur-mémoire + région-de-confiance contre sonde-couleur-mémoire | -7 [-11, -3] | 483 |
| 24 | 0.0 | sonde-couleur-mémoire + les trois contre sonde-couleur-mémoire | -28 [-36, -21] | 427 |
| 24 | 0.0 | sonde-couleur-mémoire contre Broyden-mêmes-mesures | -36 [-43, -30] | 495 |
| 24 | 0.0 | sonde-couleur-mémoire + dithering contre Broyden-mêmes-mesures | -23 [-29, -19] | 490 |
| 24 | 0.02 | sonde-couleur-mémoire + dithering contre sonde-couleur-mémoire | -2 [-7, +3] | 505 |
| 24 | 0.02 | sonde-couleur-mémoire + démarrage-conforme contre sonde-couleur-mémoire | -45 [-55, -36] | 467 |
| 24 | 0.02 | sonde-couleur-mémoire + région-de-confiance contre sonde-couleur-mémoire | -18 [-24, -13] | 489 |
| 24 | 0.02 | sonde-couleur-mémoire + les trois contre sonde-couleur-mémoire | -28 [-37, -20] | 406 |
| 24 | 0.02 | sonde-couleur-mémoire contre Broyden-mêmes-mesures | +22 [+17, +26] | 488 |
| 24 | 0.02 | sonde-couleur-mémoire + dithering contre Broyden-mêmes-mesures | +21 [+16, +26] | 499 |
| 36 | 0.0 | sonde-couleur-mémoire + dithering contre sonde-couleur-mémoire | +12 [+10, +15] | 436 |
| 36 | 0.0 | sonde-couleur-mémoire + démarrage-conforme contre sonde-couleur-mémoire | -51 [-61, -41] | 513 |
| 36 | 0.0 | sonde-couleur-mémoire + région-de-confiance contre sonde-couleur-mémoire | -25 [-28, -22] | 506 |
| 36 | 0.0 | sonde-couleur-mémoire + les trois contre sonde-couleur-mémoire | -32 [-42, -22] | 448 |
| 36 | 0.0 | sonde-couleur-mémoire contre Broyden-mêmes-mesures | -22 [-26, -18] | 508 |
| 36 | 0.0 | sonde-couleur-mémoire + dithering contre Broyden-mêmes-mesures | -9 [-12, -7] | 434 |
| 36 | 0.02 | sonde-couleur-mémoire + dithering contre sonde-couleur-mémoire | +11 [+8, +14] | 505 |
| 36 | 0.02 | sonde-couleur-mémoire + démarrage-conforme contre sonde-couleur-mémoire | -18 [-24, -13] | 470 |
| 36 | 0.02 | sonde-couleur-mémoire + région-de-confiance contre sonde-couleur-mémoire | +10 [+6, +13] | 480 |
| 36 | 0.02 | sonde-couleur-mémoire + les trois contre sonde-couleur-mémoire | -3 [-8, +2] | 389 |
| 36 | 0.02 | sonde-couleur-mémoire contre Broyden-mêmes-mesures | -36 [-41, -31] | 502 |
| 36 | 0.02 | sonde-couleur-mémoire + dithering contre Broyden-mêmes-mesures | -21 [-25, -17] | 499 |
