# Partie 4.2 — tableaux

Cellules : succès % / lectures médianes (succès seulement). Défaut = défaut ℍ-linéaire de la jacobienne au point fixe.


## A. ρ = 0.8, β = 1

| n | activation | défaut | sonde-couleur-mémoire-quaternionique | sonde-couleur-mémoire-quaternionique-départ-identité | sonde-couleur-mémoire | Broyden-sans-couleur | Broyden-sans-couleur-départ-identité | Picard-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| 2 | tanh par composante | 0.18 | 37 % / 8 | 100 % / 11 | 83 % / 14 | 100 % / 12 | 100 % / 9 | 100 % / 18 |
| 2 | tanh radial | 0.15 | 97 % / 10 | 100 % / 10 | 90 % / 13 | 100 % / 12 | 100 % / 8 | 100 % / 11 |
| 2 | saturation globale | 0.12 | 73 % / 9 | 100 % / 10 | 90 % / 13 | 100 % / 12 | 100 % / 8 | 100 % / 12 |
| 8 | tanh par composante | 0.17 | 93 % / 19 | 97 % / 23 | 100 % / 41 | 100 % / 37 | 100 % / 14 | 100 % / 17 |
| 8 | tanh radial | 0.12 | 97 % / 17 | 100 % / 13 | 93 % / 39 | 100 % / 36 | 100 % / 10 | 100 % / 12 |
| 8 | saturation globale | 0.05 | 93 % / 16 | 100 % / 14 | 93 % / 40 | 100 % / 36 | 100 % / 11 | 100 % / 13 |
| 16 | tanh par composante | 0.12 | 100 % / 27 | 90 % / 28 | 93 % / 74 | 100 % / 69 | 100 % / 14 | 100 % / 20 |
| 16 | tanh radial | 0.10 | 100 % / 26 | 93 % / 18 | 97 % / 73 | 100 % / 68 | 100 % / 12 | 100 % / 15 |
| 16 | saturation globale | 0.03 | 100 % / 25 | 97 % / 19 | 93 % / 76 | 100 % / 68 | 100 % / 12 | 100 % / 14 |
| 32 | tanh par composante | 0.13 | 93 % / 58 | 100 % / 31 | 87 % / 144 | 100 % / 133 | 100 % / 14 | 100 % / 17 |
| 32 | tanh radial | 0.10 | 93 % / 47 | 97 % / 21 | 80 % / 144 | 100 % / 133 | 100 % / 12 | 100 % / 13 |
| 32 | saturation globale | 0.02 | 100 % / 44 | 100 % / 22 | 87 % / 163 | 100 % / 132 | 100 % / 12 | 100 % / 14 |


## A. ρ = 0.95, β = 1

| n | activation | défaut | sonde-couleur-mémoire-quaternionique | sonde-couleur-mémoire-quaternionique-départ-identité | sonde-couleur-mémoire | Broyden-sans-couleur | Broyden-sans-couleur-départ-identité | Picard-sans-couleur |
|---|---|---|---|---|---|---|---|---|
| 2 | tanh par composante | 0.22 | 70 % / 12 | 100 % / 12 | 77 % / 15 | 100 % / 12 | 100 % / 10 | 100 % / 37 |
| 2 | tanh radial | 0.17 | 87 % / 9 | 100 % / 10 | 93 % / 14 | 100 % / 12 | 100 % / 8 | 100 % / 14 |
| 2 | saturation globale | 0.14 | 57 % / 9 | 100 % / 10 | 90 % / 13 | 100 % / 12 | 100 % / 9 | 100 % / 14 |
| 8 | tanh par composante | 0.19 | 83 % / 22 | 77 % / 55 | 100 % / 41 | 100 % / 37 | 100 % / 18 | 100 % / 28 |
| 8 | tanh radial | 0.14 | 77 % / 20 | 100 % / 18 | 83 % / 42 | 100 % / 37 | 100 % / 12 | 100 % / 15 |
| 8 | saturation globale | 0.06 | 93 % / 16 | 100 % / 16 | 83 % / 40 | 100 % / 36 | 100 % / 12 | 100 % / 16 |
| 16 | tanh par composante | 0.14 | 87 % / 28 | 100 % / 34 | 90 % / 80 | 100 % / 69 | 100 % / 18 | 100 % / 34 |
| 16 | tanh radial | 0.11 | 100 % / 27 | 83 % / 24 | 97 % / 73 | 100 % / 69 | 100 % / 15 | 100 % / 20 |
| 16 | saturation globale | 0.03 | 97 % / 25 | 90 % / 24 | 100 % / 75 | 100 % / 68 | 100 % / 15 | 100 % / 19 |
| 32 | tanh par composante | 0.15 | 83 % / 55 | 97 % / 37 | 73 % / 147 | 100 % / 133 | 100 % / 18 | 100 % / 22 |
| 32 | tanh radial | 0.12 | 93 % / 54 | 97 % / 30 | 83 % / 150 | 100 % / 133 | 100 % / 14 | 100 % / 16 |
| 32 | saturation globale | 0.03 | 97 % / 55 | 100 % / 27 | 83 % / 145 | 100 % / 132 | 100 % / 15 | 100 % / 18 |


## B. défaut contrôlé (n = 8, ρ = 0,95)

| activation | β | part réelle t | défaut | sonde-couleur-mémoire-quaternionique | sonde-couleur-mémoire-quaternionique-départ-identité | sonde-couleur-mémoire | Broyden-sans-couleur | Broyden-sans-couleur-départ-identité | Picard-sans-couleur |
|---|---|---|---|---|---|---|---|---|---|
| saturation globale | 0.25 | 0.0 | 0.02 | 100 % / 18 | 100 % / 31 | 93 % / 44 | 100 % / 37 | 100 % / 26 | 100 % / 50 |
| saturation globale | 0.5 | 0.0 | 0.04 | 100 % / 17 | 93 % / 22 | 100 % / 42 | 100 % / 37 | 100 % / 18 | 100 % / 29 |
| saturation globale | 2.0 | 0.0 | 0.05 | 100 % / 15 | 100 % / 10 | 97 % / 39 | 100 % / 36 | 100 % / 8 | 100 % / 9 |
| saturation globale | 4.0 | 0.0 | 0.03 | 100 % / 14 | 97 % / 9 | 100 % / 41 | 100 % / 35 | 100 % / 6 | 100 % / 6 |
| tanh radial | 0.25 | 0.0 | 0.04 | 100 % / 18 | 100 % / 32 | 97 % / 41 | 100 % / 38 | 100 % / 28 | 100 % / 49 |
| tanh radial | 0.5 | 0.0 | 0.09 | 83 % / 20 | 90 % / 25 | 100 % / 42 | 100 % / 38 | 100 % / 20 | 100 % / 28 |
| tanh radial | 2.0 | 0.0 | 0.15 | 100 % / 17 | 100 % / 10 | 93 % / 40 | 100 % / 36 | 100 % / 8 | 100 % / 10 |
| tanh radial | 4.0 | 0.0 | 0.11 | 97 % / 16 | 100 % / 9 | 97 % / 39 | 100 % / 36 | 100 % / 6 | 100 % / 6 |
| tanh par composante | 0.25 | 0.0 | 0.04 | 77 % / 20 | 93 % / 33 | 90 % / 44 | 100 % / 37 | 100 % / 32 | 100 % / 64 |
| tanh par composante | 0.5 | 0.0 | 0.10 | 90 % / 18 | 77 % / 33 | 83 % / 44 | 100 % / 38 | 100 % / 24 | 100 % / 43 |
| tanh par composante | 2.0 | 0.0 | 0.27 | 93 % / 24 | 97 % / 23 | 100 % / 41 | 100 % / 37 | 100 % / 12 | 100 % / 15 |
| tanh par composante | 4.0 | 0.0 | 0.27 | 87 % / 19 | 90 % / 11 | 83 % / 39 | 100 % / 36 | 100 % / 8 | 100 % / 10 |
| tanh par composante | 1.0 | 0.1 | 0.22 | 83 % / 29 | 90 % / 44 | 90 % / 42 | 100 % / 37 | 100 % / 16 | 100 % / 22 |
| tanh par composante | 1.0 | 0.3 | 0.25 | 73 % / 26 | 30 % / 32 | 100 % / 41 | 100 % / 37 | 100 % / 16 | 100 % / 26 |
| tanh par composante | 1.0 | 0.6 | 0.32 | 70 % / 28 | 87 % / 52 | 100 % / 41 | 100 % / 37 | 100 % / 17 | 100 % / 36 |
