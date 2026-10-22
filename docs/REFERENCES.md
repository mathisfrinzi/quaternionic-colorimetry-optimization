# Références et positionnement

Revue non systématique : ce sont les travaux sur lesquels l'étude s'appuie ou qu'elle retrouve en route.

## Ce qui vient du dépôt complexe

- **Mathis Frinzi, [`complex-colorimetry-optimization`](https://github.com/mathisfrinzi/complex-colorimetry-optimization).** Le capteur (luminosité en dent de scie, teinte), la sonde à deux lectures, l'équation `W h = −1`, la sonde-Broyden, la sonde à cible pour les minima non nuls. Les fichiers repris tels quels sont dans `vendor/complex_colorimetry/` (licence MIT).

## Quaternions et fonctions de variable quaternionique

- **G. Gentili, D. C. Struppa (2007)**, *A new theory of regular functions of a quaternionic variable*, Advances in Mathematics 216. Les fonctions slice-régulières de la partie 1.1, et le fait que leurs zéros non réels forment des sphères.
- **D. Janovská, G. Opfer (2010)**, *A note on the computation of all zeros of simple quaternionic polynomials*, SIAM J. Numer. Anal. 48. Calcul des zéros de polynômes quaternioniques ; le Newton quaternionique est le Newton réel en dimension 4 (voir les [fausses pistes](../fausses-pistes/README.md)).
- **S. J. Sangwine (1996)**, *Fourier transforms of colour images using quaternion or hypercomplex numbers*, Electronics Letters 32. L'idée de coder une couleur RVB dans la partie imaginaire d'un quaternion, d'où l'axe en RVB `(u+1)/2` de la partie 0.

## Méthodes numériques

- **C. G. Broyden (1965)**, *A class of methods for solving nonlinear simultaneous equations*, Mathematics of Computation 19. La mise à jour de rang 1 de la sonde-couleur-mémoire (partie 1) et, dans sa version qui garde la structure, du Broyden quaternionique (parties 3 et 5.5).
- **L. K. Schubert (1970)**, [*Modification of a quasi-Newton method for nonlinear equations with a sparse Jacobian*](https://www.ams.org/mcom/1970-24-109/S0025-5718-1970-0258276-9/), Mathematics of Computation 24, 27–30. La mise à jour de Broyden restreinte à un motif creux : premier exemple de mise à jour sécante qui garde une structure imposée.
- **J. E. Dennis, H. F. Walker (1981)**, [*Convergence theorems for least-change secant update methods*](https://epubs.siam.org/doi/10.1137/0718067), SIAM J. Numer. Anal. 18, 949–987. Cadre général des mises à jour sécantes « de moindre changement » dans un sous-espace affine de matrices. Le Broyden quaternionique (blocs `L(a)`, `R(b)`, diagonale par sortie) est un cas de ce cadre, avec un sous-espace propre aux quaternions ; l'étude n'a pas trouvé de version quaternionique publiée, sans que la recherche ait été exhaustive.
- **J. E. Dennis, R. B. Schnabel (1983)**, *Numerical Methods for Unconstrained Optimization and Nonlinear Equations*, Prentice-Hall. Référence de manuel pour Broyden, les recherches linéaires et les mises à jour sécantes structurées.
- **R. B. Kearfott (1979)**, *An efficient degree-computation method for a generalized method of bisection*, Numerische Mathematik 32. Calcul du degré topologique par triangulation, base du certificat de la partie 1.5.
- **J. Bey (1995)**, *Tetrahedral grid refinement*, Computing 55. La subdivision des tétraèdres en 8 utilisée pour trianguler S³ (partie 1.5).
- **L. Schuchman (1964)**, *Dither signals and their effect on quantization noise*, IEEE Trans. Communication Technology 12. Le dithering soustractif de la partie 1.4.

## Suivi d'un zéro qui bouge (partie 5.8)

- **A. Simonetto, E. Dall'Anese, S. Paternain, G. Leus, G. B. Giannakis (2020)**, [*Time-varying convex optimization: time-structured algorithms and applications*](https://arxiv.org/abs/2006.08500), Proceedings of the IEEE 108. Revue des méthodes prédiction-correction pour suivre la solution d'un problème qui varie dans le temps.
- **R. Abbassi, H. Jerbi, M. Kchaou, T. E. Simos, S. D. Mourtas, V. N. Katsikis (2023)**, [*Towards higher-order zeroing neural networks for calculating quaternion matrix inverse with application to robotic motion tracking*](https://www.mdpi.com/2227-7390/11/12/2756), Mathematics 11(12), 2756. Les réseaux de neurones « à annulation » (*zeroing neural networks*) appliqués à un problème quaternionique variant dans le temps. La partie 5.8 n'utilise ni prédiction ni modèle de dérive : elle mesure seulement le retard de méthodes de Newton-Broyden sans test de descente.

## Réseaux de neurones quaternioniques

- **T. Parcollet, M. Morchid, G. Linarès (2020)**, *A survey of quaternion neural networks*, Artificial Intelligence Review 53. Couches à produit de Hamilton, activations « split » composante par composante : c'est le cas « presque ℍ-linéaire » des parties 4 et 6.
- **D. Sussillo, O. Barak (2013)**, *Opening the black box: low-dimensional dynamics in high-dimensional recurrent neural networks*, Neural Computation 25. Recherche des points fixes d'un RNN comme problème d'annulation : la partie 6.1.

## Ce qui est propre à cette étude

- La couleur d'un quaternion (luminosité, angle, axe) comme **seul** capteur, et la mesure de ce qu'elle coûte face aux méthodes classiques, à lectures et à bits égaux.
- Le degré topologique lu dans des couleurs quantifiées, comme certificat de zéro (partie 1.5) et de point fixe (partie 6.1).
- Un Broyden quaternionique à structure imposée (transposition directe de Broyden à l'arithmétique de ℍ, cas particulier du cadre de Dennis et Walker ; son gain est celui qu'on attend de 4n² paramètres au lieu de 16n²) (`L(a)` à gauche, `R(b)` à droite pour les poids d'un neurone, diagonale par sortie pour une couche), avec le choix du côté du quotient (le mauvais côté échoue à 0 %), dans une version qui lit la valeur et une version pilotée par les couleurs. La version couleur bat Broyden-sans-couleur quand la jacobienne est structurée, mais **pas** la version classique qui connaît la même structure (partie 5.5) : l'apport de la couleur est qu'elle suffit.
