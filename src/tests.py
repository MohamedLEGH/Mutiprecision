# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED
# MAIN 3
# 09.02.16

# Question 1 
# calcule de la norme L2 d'un filtre modélisé par les matrices A, B, C, et D
# Le rayon spectral de A doit être inférieur à 1

#import
import numpy as np
from flint import *
from random_dSS import random_dSS
from scipy.linalg import solve_discrete_lyapunov
from pydare.dlyap import dlyap_schur, dlyap_slycot
from resolution_naive import lyap_naiv
from arb_resolution_naive import arb_lyap_naiv
from arb_mat_functions import arb_trace
from math import sqrt

def convert_arb_mat(M):
	n, m, = M.shape
	l = []
	for ligne in M.tolist():
		for coef in ligne:
			l.append(coef)
	return arb_mat( n, m, l)


# Initialisation
n,p,q = 10,5,5
A, B, C, D = random_dSS( n, p, q)
print("Rayon spectral de A: " + str(np.amax( np.absolute( np.linalg.eigvals( A)))) + "\n")

# version intervalle des matrices
Ai = convert_arb_mat(A)
Bi = convert_arb_mat(B)
Ci = convert_arb_mat(C)
Di = convert_arb_mat(D)

# Transposés des matrices
At = A.transpose()
Bt = B.transpose()
Ct = C.transpose()
Dt = D.transpose()

# version intervalle
Ait = convert_arb_mat(At)
Bit = convert_arb_mat(Bt)
Cit = convert_arb_mat(Ct)
Dit = convert_arb_mat(Dt)


#Calcul de W
# Equation de Lyapunov W = A W At + B Bt

# References
Wscipy = solve_discrete_lyapunov( A, B*Bt)
ref = A*Wscipy*At + B*Bt
print( "Précision résolution scipy: " + str(np.amax( Wscipy - ref)) )

Wschur = dlyap_schur( A, B*Bt)
ref = A*Wschur*At + B*Bt
print( "Précision résolution schur: " + str(np.amax( Wschur - ref)) )

Wslycot = dlyap_slycot( A, B*Bt)
ref = A*Wslycot*At + B*Bt
print( "Précision résolution slycot: " + str(np.amax( Wslycot - ref)) )

# Méthodes "maison"
Wnaiv = lyap_naiv(A, B*Bt)
ref = A*Wnaiv*At + B*Bt
print( "Précision résolution naive: " + str(np.amax( Wnaiv - ref)) )

Warb_naiv = arb_lyap_naiv(Ai, Bi*Bit)
ref = Ai*Warb_naiv*Ait + Bi*Bit
print( Warb_naiv)
#print( "Précision résolution arb_naive: " + str(np.amax( Warb_naiv - ref)) )

#Norme L2
print ""

Hscipy = sqrt( np.matrix.trace( C*Wscipy*Ct + D*Dt))
Hschur = sqrt( np.matrix.trace( C*Wschur*Ct + D*Dt + D*Dt))
Hslycot = sqrt( np.matrix.trace( C*Wslycot*Ct + D*Dt + D*Dt))

Hnaiv = sqrt( np.matrix.trace( C*Wnaiv*Ct + D*Dt + D*Dt))
Harb_naiv = arb.sqrt( arb_trace( Ci*Warb_naiv*Cit + Di*Dit + Di*Dit))

print ("Scipy: La norme L2 est " + str(Hscipy))
print ("Schur: La norme L2 est " + str(Hschur))
print ("slycot: La norme L2 est " + str(Hslycot))
print ("Naiv: La norme L2 est " + str(Hnaiv))
print ("Arb_naiv: La norme L2 est " + str(Harb_naiv))

