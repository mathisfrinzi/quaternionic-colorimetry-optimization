# Partie 4.1 — tableaux

Défaut ℍ-linéaire de la jacobienne au zéro (moyenne des instances) :

| famille | n=2 | n=3 | n=4 |
|---|---|---|---|
| ℍ-linéaire (témoin) | 0.00 | 0.00 | 0.00 |
| bilatère a·d·b | 0.91 | 0.91 | 0.77 |
| polynôme non commutatif | 0.00 | 0.00 | 0.00 |
| inversion réseau quaternionique | 0.39 | 0.44 | 0.55 |
| point fixe réservoir quaternionique | 0.10 | 0.12 | 0.14 |
| générique (témoin) | 0.60 | 0.60 | 0.62 |


## 24 bits par composante, bruit 0.0 — succès % / lectures médianes

| famille | méthode | n=2 | n=3 | n=4 |
|---|---|---|---|---|
| ℍ-linéaire (témoin) | sonde-couleur-mémoire | 100 / 24 | 98 / 36 | 95 / 54 |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique | 100 / 16 | 98 / 20 | 98 / 23 |
| ℍ-linéaire (témoin) | sonde-projeter-puis-affiner | 98 / 18 | 100 / 24 | 87 / 35 |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | 100 / 16 | 97 / 20 | 100 / 24 |
| ℍ-linéaire (témoin) | Broyden-sans-couleur | 100 / 17 | 100 / 23 | 100 / 29 |
| ℍ-linéaire (témoin) | Broyden-mêmes-mesures | 98 / 31 | 93 / 48 | 88 / 61 |
| bilatère a·d·b | sonde-couleur-mémoire | 100 / 24 | 97 / 36 | 88 / 84 |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique | 0 / — | 0 / — | 0 / — |
| bilatère a·d·b | sonde-projeter-puis-affiner | 13 / 104 | 0 / — | 2 / 240 |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique-puis-libre | 10 / 127 | 3 / 108 | 2 / 164 |
| bilatère a·d·b | Broyden-sans-couleur | 100 / 15 | 100 / 20 | 100 / 24 |
| bilatère a·d·b | Broyden-mêmes-mesures | 97 / 31 | 97 / 52 | 92 / 73 |
| polynôme non commutatif | sonde-couleur-mémoire | 88 / 66 | 80 / 70 | 52 / 94 |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique | 92 / 17 | 95 / 20 | 88 / 32 |
| polynôme non commutatif | sonde-projeter-puis-affiner | 88 / 42 | 82 / 49 | 72 / 69 |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique-puis-libre | 88 / 17 | 97 / 20 | 82 / 31 |
| polynôme non commutatif | Broyden-sans-couleur | 87 / 35 | 88 / 43 | 75 / 69 |
| polynôme non commutatif | Broyden-mêmes-mesures | 70 / 44 | 72 / 59 | 57 / 86 |
| inversion réseau quaternionique | sonde-couleur-mémoire | 80 / 26 | 45 / 75 | 18 / 84 |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique | 42 / 48 | 3 / 137 | 0 / — |
| inversion réseau quaternionique | sonde-projeter-puis-affiner | 87 / 34 | 33 / 96 | 15 / 145 |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | 73 / 53 | 25 / 100 | 13 / 145 |
| inversion réseau quaternionique | Broyden-sans-couleur | 93 / 24 | 92 / 33 | 85 / 42 |
| inversion réseau quaternionique | Broyden-mêmes-mesures | 70 / 37 | 48 / 67 | 33 / 102 |
| point fixe réservoir quaternionique | sonde-couleur-mémoire | 100 / 17 | 100 / 21 | 100 / 25 |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique | 100 / 11 | 100 / 12 | 100 / 14 |
| point fixe réservoir quaternionique | sonde-projeter-puis-affiner | 100 / 11 | 100 / 12 | 100 / 14 |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | 100 / 11 | 100 / 12 | 100 / 14 |
| point fixe réservoir quaternionique | Broyden-sans-couleur | 100 / 13 | 100 / 17 | 100 / 21 |
| point fixe réservoir quaternionique | Broyden-mêmes-mesures | 95 / 17 | 98 / 23 | 93 / 30 |
| générique (témoin) | sonde-couleur-mémoire | 87 / 28 | 70 / 74 | 87 / 84 |
| générique (témoin) | sonde-couleur-mémoire-quaternionique | 25 / 167 | 12 / 201 | 13 / 266 |
| générique (témoin) | sonde-projeter-puis-affiner | 57 / 99 | 38 / 135 | 37 / 102 |
| générique (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | 58 / 127 | 37 / 162 | 18 / 138 |
| générique (témoin) | Broyden-sans-couleur | 100 / 15 | 100 / 20 | 100 / 24 |
| générique (témoin) | Broyden-mêmes-mesures | 90 / 42 | 88 / 56 | 87 / 87 |


### Gains appariés contre Broyden à lectures EXACTES (24 bits, σ = 0.0) ; + = la sonde est moins chère ; Δ succès en points

| famille | sonde | n=2 | n=3 | n=4 |
|---|---|---|---|---|
| ℍ-linéaire (témoin) | sonde-couleur-mémoire | -58 [-74, -47] (+0 pts) | -97 [-125, -75] (-2 pts) | -108 [-131, -89] (-5 pts) |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique | +5 [+3, +7] (+0 pts) | +13 [+9, +17] (-2 pts) | +18 [+16, +20] (-2 pts) |
| ℍ-linéaire (témoin) | sonde-projeter-puis-affiner | -31 [-47, -17] (-2 pts) | -22 [-35, -10] (+0 pts) | -26 [-41, -13] (-13 pts) |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | +5 [+4, +7] (+0 pts) | +15 [+12, +17] (-3 pts) | +16 [+12, +20] (+0 pts) |
| bilatère a·d·b | sonde-couleur-mémoire | -58 [-60, -56] (+0 pts) | -116 [-140, -96] (-3 pts) | -238 [-265, -215] (-12 pts) |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique | — | — | — |
| bilatère a·d·b | sonde-projeter-puis-affiner | -554 [-809, -369] (-87 pts) | — | — |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique-puis-libre | -722 [-1170, -441] (-90 pts) | — | — |
| polynôme non commutatif | sonde-couleur-mémoire | -69 [-93, -48] (+2 pts) | -41 [-66, -22] (-8 pts) | -43 [-69, -18] (-23 pts) |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique | +39 [+28, +46] (+5 pts) | +45 [+37, +52] (+7 pts) | +43 [+33, +52] (+13 pts) |
| polynôme non commutatif | sonde-projeter-puis-affiner | -19 [-41, -3] (+2 pts) | -8 [-28, +7] (-7 pts) | -9 [-28, +9] (-3 pts) |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique-puis-libre | +38 [+26, +48] (+2 pts) | +44 [+35, +51] (+8 pts) | +38 [+22, +50] (+7 pts) |
| inversion réseau quaternionique | sonde-couleur-mémoire | -50 [-83, -26] (-13 pts) | -111 [-152, -73] (-47 pts) | -94 [-139, -61] (-67 pts) |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique | -193 [-289, -126] (-52 pts) | — | — |
| inversion réseau quaternionique | sonde-projeter-puis-affiner | -84 [-138, -44] (-7 pts) | -178 [-280, -104] (-58 pts) | -287 [-445, -168] (-70 pts) |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | -166 [-219, -126] (-20 pts) | -211 [-308, -137] (-67 pts) | -326 [-465, -227] (-72 pts) |
| point fixe réservoir quaternionique | sonde-couleur-mémoire | -33 [-35, -31] (+0 pts) | -25 [-27, -23] (+0 pts) | -19 [-20, -19] (+0 pts) |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique | +11 [+10, +12] (+0 pts) | +28 [+27, +29] (+0 pts) | +34 [+32, +35] (+0 pts) |
| point fixe réservoir quaternionique | sonde-projeter-puis-affiner | +14 [+12, +15] (+0 pts) | +29 [+27, +30] (+0 pts) | +30 [+23, +34] (+0 pts) |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | +11 [+10, +12] (+0 pts) | +28 [+27, +29] (+0 pts) | +34 [+32, +35] (+0 pts) |
| générique (témoin) | sonde-couleur-mémoire | -113 [-138, -93] (-13 pts) | -230 [-272, -193] (-30 pts) | -265 [-290, -241] (-13 pts) |
| générique (témoin) | sonde-couleur-mémoire-quaternionique | -1031 [-1261, -817] (-75 pts) | -923 [-1171, -757] (-88 pts) | -1011 [-1319, -765] (-87 pts) |
| générique (témoin) | sonde-projeter-puis-affiner | -582 [-712, -479] (-43 pts) | -650 [-868, -489] (-62 pts) | -436 [-555, -346] (-63 pts) |
| générique (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | -706 [-877, -569] (-42 pts) | -677 [-851, -530] (-63 pts) | -447 [-563, -344] (-82 pts) |


## 24 bits par composante, bruit 0.02 — succès % / lectures médianes

| famille | méthode | n=2 | n=3 | n=4 |
|---|---|---|---|---|
| ℍ-linéaire (témoin) | sonde-couleur-mémoire | 100 / 28 | 97 / 60 | 87 / 84 |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique | 100 / 16 | 100 / 22 | 98 / 29 |
| ℍ-linéaire (témoin) | sonde-projeter-puis-affiner | 100 / 22 | 98 / 43 | 87 / 54 |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | 100 / 16 | 98 / 22 | 97 / 28 |
| ℍ-linéaire (témoin) | Broyden-sans-couleur | 100 / 17 | 100 / 23 | 100 / 29 |
| ℍ-linéaire (témoin) | Broyden-mêmes-mesures | 100 / 48 | 100 / 84 | 100 / 112 |
| bilatère a·d·b | sonde-couleur-mémoire | 100 / 27 | 100 / 70 | 93 / 82 |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique | 0 / — | 0 / — | 0 / — |
| bilatère a·d·b | sonde-projeter-puis-affiner | 8 / 80 | 2 / 111 | 0 / — |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique-puis-libre | 8 / 84 | 0 / — | 0 / — |
| bilatère a·d·b | Broyden-sans-couleur | 100 / 15 | 100 / 20 | 100 / 24 |
| bilatère a·d·b | Broyden-mêmes-mesures | 100 / 46 | 100 / 82 | 100 / 120 |
| polynôme non commutatif | sonde-couleur-mémoire | 77 / 60 | 72 / 81 | 42 / 104 |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique | 97 / 19 | 98 / 23 | 88 / 41 |
| polynôme non commutatif | sonde-projeter-puis-affiner | 85 / 41 | 80 / 54 | 60 / 66 |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique-puis-libre | 93 / 19 | 98 / 23 | 75 / 36 |
| polynôme non commutatif | Broyden-sans-couleur | 87 / 35 | 88 / 43 | 75 / 69 |
| polynôme non commutatif | Broyden-mêmes-mesures | 92 / 66 | 90 / 112 | 83 / 154 |
| inversion réseau quaternionique | sonde-couleur-mémoire | 80 / 48 | 23 / 92 | 2 / 111 |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique | 38 / 118 | 2 / 100 | 0 / — |
| inversion réseau quaternionique | sonde-projeter-puis-affiner | 80 / 50 | 13 / 136 | 3 / 102 |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | 70 / 60 | 15 / 194 | 2 / 116 |
| inversion réseau quaternionique | Broyden-sans-couleur | 93 / 24 | 92 / 33 | 85 / 42 |
| inversion réseau quaternionique | Broyden-mêmes-mesures | 90 / 70 | 78 / 114 | 48 / 162 |
| point fixe réservoir quaternionique | sonde-couleur-mémoire | 100 / 17 | 100 / 21 | 53 / 31 |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique | 100 / 11 | 100 / 13 | 100 / 15 |
| point fixe réservoir quaternionique | sonde-projeter-puis-affiner | 100 / 11 | 100 / 13 | 97 / 17 |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | 100 / 11 | 100 / 13 | 100 / 15 |
| point fixe réservoir quaternionique | Broyden-sans-couleur | 100 / 13 | 100 / 17 | 100 / 21 |
| point fixe réservoir quaternionique | Broyden-mêmes-mesures | 100 / 24 | 100 / 38 | 100 / 63 |
| générique (témoin) | sonde-couleur-mémoire | 87 / 51 | 78 / 78 | 70 / 91 |
| générique (témoin) | sonde-couleur-mémoire-quaternionique | 17 / 180 | 2 / 165 | 0 / — |
| générique (témoin) | sonde-projeter-puis-affiner | 38 / 103 | 22 / 110 | 18 / 154 |
| générique (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | 52 / 109 | 23 / 173 | 7 / 170 |
| générique (témoin) | Broyden-sans-couleur | 100 / 15 | 100 / 20 | 100 / 24 |
| générique (témoin) | Broyden-mêmes-mesures | 100 / 64 | 100 / 94 | 100 / 136 |


### Gains appariés contre Broyden à lectures EXACTES (24 bits, σ = 0.02) ; + = la sonde est moins chère ; Δ succès en points

| famille | sonde | n=2 | n=3 | n=4 |
|---|---|---|---|---|
| ℍ-linéaire (témoin) | sonde-couleur-mémoire | -108 [-134, -87] (+0 pts) | -163 [-198, -134] (-3 pts) | -168 [-195, -146] (-13 pts) |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique | -2 [-6, +1] (+0 pts) | -5 [-15, +2] (+0 pts) | -7 [-16, +1] (-2 pts) |
| ℍ-linéaire (témoin) | sonde-projeter-puis-affiner | -49 [-64, -36] (+0 pts) | -96 [-126, -70] (-2 pts) | -78 [-96, -63] (-13 pts) |
| ℍ-linéaire (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | -2 [-7, +1] (+0 pts) | -3 [-12, +4] (-2 pts) | -7 [-16, +1] (-3 pts) |
| bilatère a·d·b | sonde-couleur-mémoire | -141 [-176, -113] (+0 pts) | -247 [-275, -222] (+0 pts) | -242 [-257, -228] (-7 pts) |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique | — | — | — |
| bilatère a·d·b | sonde-projeter-puis-affiner | -446 [-651, -307] (-92 pts) | — | — |
| bilatère a·d·b | sonde-couleur-mémoire-quaternionique-puis-libre | -554 [-741, -433] (-92 pts) | — | — |
| polynôme non commutatif | sonde-couleur-mémoire | -67 [-94, -45] (-10 pts) | -74 [-99, -52] (-17 pts) | -56 [-83, -30] (-33 pts) |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique | +31 [+21, +39] (+10 pts) | +43 [+36, +49] (+10 pts) | +41 [+28, +51] (+13 pts) |
| polynôme non commutatif | sonde-projeter-puis-affiner | -34 [-60, -10] (-2 pts) | -9 [-27, +6] (-8 pts) | +14 [-3, +27] (-15 pts) |
| polynôme non commutatif | sonde-couleur-mémoire-quaternionique-puis-libre | +27 [+11, +38] (+7 pts) | +43 [+36, +49] (+10 pts) | +36 [+22, +47] (+0 pts) |
| inversion réseau quaternionique | sonde-couleur-mémoire | -85 [-122, -55] (-13 pts) | -132 [-191, -80] (-68 pts) | — |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique | -411 [-573, -281] (-55 pts) | — | — |
| inversion réseau quaternionique | sonde-projeter-puis-affiner | -100 [-157, -54] (-13 pts) | — | — |
| inversion réseau quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | -182 [-240, -139] (-23 pts) | -313 [-502, -174] (-77 pts) | — |
| point fixe réservoir quaternionique | sonde-couleur-mémoire | -35 [-37, -33] (+0 pts) | -33 [-42, -27] (+0 pts) | -45 [-49, -41] (-47 pts) |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique | +10 [+8, +11] (+0 pts) | +24 [+23, +25] (+0 pts) | +25 [+21, +28] (+0 pts) |
| point fixe réservoir quaternionique | sonde-projeter-puis-affiner | +8 [+4, +12] (+0 pts) | +22 [+20, +24] (+0 pts) | +11 [+1, +19] (-3 pts) |
| point fixe réservoir quaternionique | sonde-couleur-mémoire-quaternionique-puis-libre | +10 [+8, +11] (+0 pts) | +24 [+22, +25] (+0 pts) | +25 [+20, +28] (+0 pts) |
| générique (témoin) | sonde-couleur-mémoire | -222 [-264, -184] (-13 pts) | -287 [-317, -260] (-22 pts) | -278 [-304, -255] (-30 pts) |
| générique (témoin) | sonde-couleur-mémoire-quaternionique | -1108 [-1274, -954] (-83 pts) | — | — |
| générique (témoin) | sonde-projeter-puis-affiner | -578 [-736, -447] (-62 pts) | -481 [-661, -354] (-78 pts) | -458 [-569, -360] (-82 pts) |
| générique (témoin) | sonde-couleur-mémoire-quaternionique-puis-libre | -713 [-904, -559] (-48 pts) | -634 [-861, -451] (-77 pts) | — |
