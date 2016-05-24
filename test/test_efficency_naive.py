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

def isNaN(num):
    return num != num
    
# Initialisation et choix du systême
nbt = 1
ctx.prec = 700
n,p,q = 15, 15, 15
m_prec = arb(0)
m_tps = 0
m_succes = 0

# Résolution de l'Equation de Lyapunov: W = A W At + B Bt
for i in range(nbt):

	# Génération aléatoire
	A, B, C, D = random_dSS( n, p, q)

	# version arb_mat des matrices
	Ai = convert_arb_mat(A)
	Bi = convert_arb_mat(B)
	Ci = convert_arb_mat(C)
	Di = convert_arb_mat(D)

	# transposées
	Ait = arb_transpose(Ai)
	Bit = arb_transpose(Bi)
	Cit = arb_transpose(Ci)
	Dit = arb_transpose(Di)

	# Calculs
	try:
		tps1 = clock()
		Warb_naiv = arb_lyap_naiv(Ai, Bi*Bit)
		tps2 = clock()

		Harb_naiv = arb.sqrt( arb_trace( Ci*Warb_naiv*Cit + Di*Dit + Di*Dit))
		
		if(isNaN(Harb_naiv)):
			raise ValueError
			
		#print ("\tSuccès: ||H|| = " + str(Harb_naiv))
		m_prec = m_prec + Harb_naiv.rad()
		m_tps = m_tps + (tps2-tps1)
		m_succes = m_succes + 1
		print("\ttest"+str(i)+" succes"+str(m_succes))
	except:
		print("\ttest"+str(i)+" echec")

if(m_succes!=0):
	m_prec = m_prec/arb(m_succes)
	m_tps = m_tps/m_succes
	m_succes = float(m_succes)/float(nbt)
	print("Taux de succes:"+str(m_succes))
	print("Precision moyenne:"+str(m_prec.mid()))
	print("Temps moyen:"+str(m_tps))




