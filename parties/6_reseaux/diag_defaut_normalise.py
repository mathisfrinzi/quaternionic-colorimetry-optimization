"""Partie 6.5, hypothèse écartée : « le défaut qui compte pour la sonde est celui de la jacobienne NORMALISÉE par F, pas celui de la jacobienne brute ».
Mesure le défaut à droite de la jacobienne brute J et de la jacobienne normalisée R(F_i⁻¹)·J_i (ce que la sonde modélise), à 3 distances du zéro,
sur le neurone de la partie 6.3 (3 instances × 10 départs). Résultat : les deux défauts sont presque égaux (0,26–0,34 en tanh radial).
Usage : python diag_defaut_normalise.py   (écrit diag_defaut_normalise.md)"""
import sys,os,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,'..','..','src')); sys.path.insert(0,HERE)
from quat import Rmat, qinv
import banc_hn as B
from bench_entrainement import Neurone, defaut_cote
def defnorm(N,q):
    J=N.J(q); F=N.f(q).reshape(N.n,4); Jn=J.copy()
    for i in range(N.n):
        Ri=Rmat(qinv(F[i]))
        Jn[4*i:4*i+4,:]=Ri@J[4*i:4*i+4,:]
    return defaut_cote(J,"R"), defaut_cote(Jn,"R")
LIG=[]
for kind in ("radial-faible","radial"):
    for n in (4,8):
        raw=[];nor=[]
        for seed in range(3):
            N=Neurone(n,kind,seed)
            for i in range(10):
                q0=B.depart(N.qs,n,i,seed,0.5)
                for scale in (1.0,0.3,0.1):          # distance au zéro : 0.5√n, 0.15√n, 0.05√n
                    q=N.qs+(q0-N.qs)*scale
                    a,b=defnorm(N,q); raw.append((scale,a)); nor.append((scale,b))
        for scale in (1.0,0.3,0.1):
            ligne="| %s | %d | %.2f | %.2f | %.2f |"%(kind,n,scale,np.mean([a for s,a in raw if s==scale]),np.mean([b for s,b in nor if s==scale]))
            print(ligne); LIG.append(ligne)
open(os.path.join(HERE,"diag_defaut_normalise.md"),"w").write("# Défaut à droite, jacobienne brute et normalisée\n\n| activation | n | distance au zéro (× la distance de départ) | défaut brut | défaut normalisé |\n|---|---|---|---|---|\n"+"\n".join(LIG)+"\n")
