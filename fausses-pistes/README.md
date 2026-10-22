# Fausses pistes (prévues par le plan, puis rencontrées ou écartées)

Cinq fausses pistes étaient annoncées dans le plan. Voici ce qu'il en est, plus celles trouvées en route. Les noms des méthodes sont ceux du rapport ([partie 0](../parties/0_capteurs_et_sondes/README.md)).

## Prévues

1. **« Le Newton quaternionique n'existe pas. »** Faux : il existe, et c'est Newton réel en 4D (Newton-sans-couleur, Newton-jacobienne-exacte). C'est lui, et Broyden, qui sont la référence, pas un Newton « complexe ». Résultat : il gagne à lectures égales sur les fonctions quelconques (parties 1.1–1.2).
2. **Non-commutativité de Hamilton.** Le quotient est ρ = F(q)⁻¹F(q+h), à gauche ; l'ordre compte. Garde-fou : `tests/test_quat.py` (associativité, L/R) et `tests/test_slice.py` (sur la tranche ℂ_i la sonde quaternionique redonne le complexe). La partie 1.1 sert de test de cohérence (critère d'arrêt passé).
3. **Visualiser les sphères de zéros en 3D.** Fait en partie 1.1 (`figures/03_spheres_de_zeros.png`), mais attention : un zéro de fonction slice-régulière est une *sphère* ; atteindre un point de la sphère suffit pour le succès.
4. **Ne jamais annoncer de victoire hors budget égal.** Toutes les tables comparent à lectures égales, avec départs identiques (comparaisons appariées) et intervalles de confiance quand l'effectif le permet.
5. **La sphère de zéros est une cible facile.** Chaque tableau étiquette le type de cible (point / sphère / courbe / surface / aucun zéro).

## Trouvées en dimension 4 (parties 0 à 2)

- **Degrés de liberté : 7, pas 4.** Une application conforme h ↦ αhβ a 7 paramètres effectifs (le « 4 » du plan vaut seulement pour q ↦ aq+b). Le gain annoncé « 2 lectures contre 4–5 » devient : 3 sondes + base = 4 lectures contre 5 pour une fonction quelconque.
- **La couleur-sans-axe non signée perd la conjugaison.** Avec θ ∈ [0, π], q et q̄ sont indistinguables : 59 % de succès sur les tranches contre 100 % en version signée (partie 1.1).
- **La couleur-sans-axe est aveugle au mauvais conditionnement.** Plateau à ~75 % quel que soit le nombre de bits, et 2–9 % sur la famille `ill` (partie 1.2).
- **L'axe en RVB est redondant** (|u| = 1) : à bits égaux, 2 angles battent RVB en précision (partie 0).
- **Le capteur-valeurs exige ≥ 48 bits** pour fonctionner (0 % à 24 bits).
- **Le bruit de réglage (≈ ±30 % sur le coût)** est du même ordre que beaucoup d'écarts entre méthodes proches : ne pas interpréter un écart de moins de 30 % sans comparaison appariée.
- **À 36 bits, la sonde-couleur-conforme-3-sondes est *moins* bonne que la sonde-couleur** (+60 % de lectures) ; le gain de 40 % à 24 bits n'est donc pas un effet de structure seul (partie 1.3).

## Trouvées en améliorant la sonde (partie 1.4)

- **Le dithering n'aide pas Broyden.** Avec des différences finies très fines, une erreur de quantification devenue aléatoire se transforme en bruit sur la jacobienne (97 à 114 lectures au lieu de 12–14). Il n'aide que la sonde, qui sonde à grande distance (η ≈ 0,3).
- **Région de confiance, première version** : le rapport « baisse lue / baisse prédite » en logarithme vaut −∞ pour un pas de Newton complet (la prédiction est |1 + W h| = 0), et le rayon s'effondrait (0 % de succès). Corrigé en baisses relatives bornées et en rayon minimal ; elle reste moins bonne que le simple retour en arrière.
- **Démarrage conforme sur des fonctions non conformes** : économiser une lecture au départ coûte une jacobienne moins bonne, et donc 18 à 51 % de lectures en plus.

## Trouvées en dimension 4n (partie 3)

- **Structure ℍ-linéaire hors de son domaine (partie 3.1)** : dès que la fonction s'écarte de la structure (ε = 1), la sonde structurée échoue presque toujours. Une sonde structurée ne doit pas être utilisée « par défaut ».
- **Amortir la mise à jour de la mémoire (partie 3.3).** Hypothèse : en grande dimension, le bruit de quantification de ρ, divisé par ‖h‖² dans la mise à jour de Broyden, abîme `W`. Diviser plutôt par ‖h‖² + a·thr² dégrade presque partout (n = 16, bruit 0,02 : 25 à 42 % de succès contre 75 %). La cause du décrochage est ailleurs : l'ambiguïté de la luminosité.
- **Ignorer le drapeau d'ambiguïté (sonde-sans-drapeau, partie 3.3).** Garder la décade la plus proche de la prédiction même quand elle est douteuse fait tomber le succès à 0–8 % : une décade fausse empoisonne la mémoire. Il faut vérifier la décade, pas la deviner.
- **Une première sonde « partielle » (partie 3.3)**, qui excluait les composantes ambiguës du test d'acceptation mais rebasait la mémoire avec leur valeur devinée : 0 % de succès à n = 8 (test rapide, 2 départs). Écartée ; l'idée de traiter les composantes une par une reste à reprendre proprement.

## Trouvées sur les autres structures (partie 4)

- **« Projeter sur les ℍ-linéaires puis affiner » par rang 1 (partie 4.1).** Partir d'une jacobienne ℍ-linéaire et la corriger par Broyden réel (sonde-projeter-puis-affiner), ou ne basculer qu'au premier pas refusé (sonde-couleur-mémoire-quaternionique-puis-libre), ne rattrape pas une structure fausse : 2 à 37 % de succès sur les familles générique et bilatère, contre 87–88 % sans structure. Il manque 12n² paramètres, une correction de rang 1 par pas ne suffit pas.
- **Inverser un réseau quaternionique à activation composante par composante (partie 4.1).** Le `tanh` appliqué à chaque composante réelle casse la structure ℍ-linéaire (défaut 0,4–0,55) : toutes les sondes-couleur échouent, alors que Broyden-sans-couleur réussit à 85–93 %.
- **Initialiser la sonde « comme Broyden-départ-identité », sans aucune lecture (partie 4.2).** La direction de chaque Fᵢ se lit dans la couleur, mais pas sa décade ; on a tenté de la fixer avec le premier pas lui-même, et d'autoriser des pas de Newton complets. Résultat : moins robuste que la sonde-couleur-mémoire-quaternionique-départ-identité (3/8 succès au lieu de 5/8 à n = 8 avec pas 0,8 ; 0/8 avec pas complets), sans gain de lectures. Raison : après un pas de Newton complet, |F| peut chuter de plusieurs décades, et la luminosité, en dent de scie, ne dit pas de combien. La sonde doit garder des pas amortis (|F| divisé par 5 à 10 au plus) pour que la baisse reste lisible : c'est ce qui lui coûte environ deux fois plus de lectures que Broyden-départ-identité, qui lit |F| directement.

## Trouvées en grande dimension et sur les réseaux (parties 3.4 à 6)

- **Mettre à jour la mémoire composante par composante (partie 3.4).** Accepter un pas sur les seules composantes sûres laisse entrer les décades fausses muettes : 0 à 9 % de succès.
- **La relecture sans double dépliage (partie 3.4).** Elle ne vérifie que les composantes qui s'écartent de la prédiction et laisse passer les erreurs de plus d'une demi-décade, qui se replient sans drapeau : 18 à 33 % à n = 32.
- **Le vote entre composantes (partie 3.5).** Les ambiguïtés sont corrélées (partie 3.3) ; on a corrigé les composantes douteuses par le décalage médian des composantes sûres. Résultat : 22 % de succès à n = 16 sans bruit, 5 % avec, contre 82 % et 72 % sans vote. Le décalage commun existe, mais il ne vaut pas une décade entière : redéplier autour de lui fait passer des décades fausses pour sûres, qui empoisonnent la mémoire.
- **Contraindre la pièce apprise d'un RNN à rester ℍ-linéaire (partie 6.1).** On espérait garder le défaut bas ; il ne baisse pas nettement (0,10–0,23 contre 0,15–0,17 pour une pièce libre ; plus haut avec `tanh` par composante et `tanh` radial à n = 4, plus bas seulement avec la saturation globale), parce qu'une pièce contrainte doit être plus grande pour apprendre la même tâche, et que le défaut vient de `Dσ·W`.
- **Changer l'activation pour inverser un réseau (partie 6.2).** La saturation globale baisse le défaut au zéro jusqu'à 0,13, mais la sonde échoue toujours (3 à 50 %) : loin du zéro l'activation sature, et le défaut au zéro ne dit rien du chemin.
- **La boîte grise pour l'inversion de réseau (partie 4.3).** Connaître la première couche et la forme diagonale de `Dσ` donne la jacobienne exacte en 2n sondes, mais la mise à jour ne suit pas la saturation : 73 % à n = 2, 8 % à n = 4.

## Trouvées en comparant la couleur à la valeur (parties 5.5 à 6.5)

- **« La couleur bat la valeur exacte » (parties 3, 4 et 6).** Les gains face à Broyden-sans-couleur venaient de la structure, que cette référence ignorait. Un Broyden quaternionique qui connaît la même structure et lit la valeur est moins cher que toutes les sondes (partie 5.5). La bonne formulation est « la couleur suffit à exploiter la structure ».
- **La dent de scie pour un système qui change d'échelle (partie 5.6).** La précision relative constante de la luminosité n'aide pas : à conditionnement 10³, toutes les sondes échouent, même avec la décade donnée, alors que le même capteur lu comme une valeur réussit 43 à 73 % des départs.
- **Les quotients contre un gain aléatoire (partie 5.7).** Ils n'éliminent qu'un gain commun aux deux lectures. Un gain indépendant à chaque lecture entre deux fois dans un quotient : la sonde décroche avant la valeur absolue bruitée.
- **Garder le test de descente en suivi (partie 5.8).** Il compare des valeurs lues à des instants différents et refuse presque tous les pas ; sans adaptation, aucune méthode ne suit un zéro mobile.
- **Trois hypothèses sur l'effondrement en régime saturé (partie 6.5).** Un défaut mal mesuré (recalculé normalisé : inchangé), la coupure des petites valeurs singulières, la taille des sondes après un échec : aucune n'explique l'écart. Corriger le test de descente en le pondérant par la vraie perte aide moins que le supprimer (22 % contre 77 % à n = 8).
