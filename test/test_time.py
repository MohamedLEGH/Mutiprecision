# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED
# MAIN 3


# Ajout du dossier src au path d'import
import sys
sys.path.insert(0, '/home/lucas/Documents/Projet_S6/Mutiprecision/src')

# Bibliothèques de calcul
import numpy as np
from flint import *
from math import sqrt

# Nos fonctions et la génération de matrice aléatoire
from random_dSS import random_dSS
from arb_resolution_naive import arb_lyap_naiv
from arb_mat_functions import arb_trace, convert_arb_mat, arb_transpose

# Gestion du temps
from time import clock

# Initialisation et choix du systême 

ctx.prec = 300

t = str( raw_input("Systême de matrices (random/last/\"nomfichier\"):"))

# Génération aléatoire
if( t == "random" or t == "r" or t=="") :
	n,p,q = 15,15,15
	A, B, C, D = random_dSS( n, p, q)
	
	print("Rayon spectral de A: " + str(np.amax( np.absolute( np.linalg.eigvals( A)))) + "\n")

# Dernier appel
elif( t == "last" or t == "l" ) :
	sys = np.load("log.npz")
	A = np.matrix( sys["arr_0"] )
	B = np.matrix( sys["arr_1"] )
	C = np.matrix( sys["arr_2"] )
	D = np.matrix( sys["arr_3"] )

# Depuis le fichier passé en argument
else :
	sys = np.load(t)
	A = np.matrix( sys["arr_0"] )
	B = np.matrix( sys["arr_1"] )
	C = np.matrix( sys["arr_2"] )
	D = np.matrix( sys["arr_3"] )

# version arb_mat des matrices
Ai = convert_arb_mat(A)
Bi = convert_arb_mat(B)
Ci = convert_arb_mat(C)
Di = convert_arb_mat(D)

# version arb_mat
Ait = arb_transpose(Ai)
Bit = arb_transpose(Bi)
Cit = arb_transpose(Ci)
Dit = arb_transpose(Di)


# Résolution de l'Equation de Lyapunov: W = A W At + B Bt
tps1 = clock()
Warb_naiv = arb_lyap_naiv(Ai, Bi*Bit)
tps2 = clock()
print( str(tps2-tps1)+"sec")

# Sauvegarde du systême testé
print ""
t = str( raw_input("Save in (default/\"nomfichier\"):"))
if (t=="default" or t=="d" or t==""):
	np.savez("log", A, B, C, D)
else :
	np.savez(t, A, B, C, D)
