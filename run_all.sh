#!/bin/sh
# Tests puis bancs (compter plusieurs heures sur 2 cœurs). Chaque banc écrit son JSON, puis son script *_figs.py fait figures et tableaux.
# Ordre du rapport : 0 capteurs et sondes, 1 zéros, 2 minimums, 3 dimension 4n, 4 autres structures, 5 à armes égales, 6 réseaux.
set -e
cd "$(dirname "$0")"
P=parties
for t in tests/test_*.py; do python "$t"; done
# 0. Capteurs et sondes
python $P/0_capteurs_et_sondes/bench_capteurs.py && python $P/0_capteurs_et_sondes/bench_capteurs_figs.py
# 1. Zéros
python $P/1_zeros/1.1_tranche/bench_tranches.py && python $P/1_zeros/1.1_tranche/bench_tranches_figs.py
python $P/1_zeros/1.2_fonctions_quelconques/bench_generique.py && python $P/1_zeros/1.2_fonctions_quelconques/bench_bits.py && python $P/1_zeros/1.2_fonctions_quelconques/bench_generique_figs.py
python $P/1_zeros/1.3_structure_conforme/bench_mobius.py && python $P/1_zeros/1.3_structure_conforme/bench_mobius_figs.py
python $P/1_zeros/1.4_ameliorations/bench_ameliorations.py && python $P/1_zeros/1.4_ameliorations/bench_ameliorations_figs.py
python $P/1_zeros/1.5_certificat_degre/bench_degre.py && python $P/1_zeros/1.5_certificat_degre/bench_degre_figs.py
# 2. Minimums
python $P/2_minimums/bench_minimum.py && python $P/2_minimums/bench_suivi.py && python $P/2_minimums/bench_ablation.py && python $P/2_minimums/bench_ablation_figs.py
# 3. Dimension 4n
python $P/3_dimension_4n/3.1_jacobienne_H_lineaire/bench_hn.py && python $P/3_dimension_4n/3.1_jacobienne_H_lineaire/bench_hn_figs.py
python $P/3_dimension_4n/3.2_grande_dimension/bench_hn_grand.py
python $P/3_dimension_4n/3.3_decrochage/run_diag.py && python $P/3_dimension_4n/3.3_decrochage/run_diag.py --sigma 0.02
python $P/3_dimension_4n/3.4_ambiguite/bench34.py --eval && python $P/3_dimension_4n/3.4_ambiguite/bench34_figs.py   # plus de 2 h sur 2 cœurs, reprise possible
python $P/3_dimension_4n/3.5_deux_echelles/bench35.py && python src/figs_hn.py 35   # relit les cases communes dans le JSON de 3.4 ; plus d'une heure
# 4. Autres structures
python $P/4_autres_structures/4.1_hors_structure/bench_hors_structure.py && python $P/4_autres_structures/4.1_hors_structure/bench_hors_structure_figs.py
python $P/4_autres_structures/4.2_reservoirs/defaut_theorie.py && python $P/4_autres_structures/4.2_reservoirs/bench_reservoirs.py && python $P/4_autres_structures/4.2_reservoirs/bench_reservoirs_figs.py
python $P/4_autres_structures/4.3_structures_etendues/bench_structures.py && python src/figs_hn.py 43
# 5. À armes égales
python $P/5_armes_egales/bench_armes_egales.py && python $P/5_armes_egales/bench_armes_egales_figs.py
# 5.5 à 5.8 : même structure, boîte instable, amplitude peu fiable, zéro mobile (plusieurs heures ; bancs reprenables)
for b in bench_meme_structure bench_boite_instable bench_gain_amplitude bench_suivi; do python $P/5_armes_egales/$b.py; done
# 6. Réseaux quaternioniques
python $P/6_reseaux/bench_points_fixes.py && python $P/6_reseaux/bench_inversion.py && python $P/6_reseaux/bench_entrainement.py && python $P/6_reseaux/bench_couche.py
python src/figs_hn.py 61 62 63 64
python $P/6_reseaux/bench_diag_sature.py && python src/figs_hn.py 55 56 57 58 65     # 6.5 : diagnostic de la saturation
# Schémas et synthèses ; figures (i, j, k) avec / sans chaque élément de la sonde
python schemas/make_ijk.py
python schemas/make_schemas.py && python schemas/make_syntheses.py && python schemas/make_conclusion.py
