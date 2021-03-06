# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED
# MAIN 3


# Ajout du dossier src au path d'import
import sys
sys.path.insert(0, '/home/lucas/Documents/Projet_S6/Mutiprecision/src')

#imports
import numpy as np
from flint import *
from random_dSS import random_dSS
from math import sqrt
from scipy.linalg import solve_discrete_lyapunov
from pydare.dlyap import dlyap_schur, dlyap_slycot

from resolution_naive import lyap_naiv
from arb_resolution_naive import arb_lyap_naiv
from arb_mat_functions import arb_trace, convert_arb_mat
from lyap_bartels_stewart import lyap_bartel_stewart

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

# version numpy.mat+arb des matrices:
Ani = arb(1)*A
Bni = arb(1)*B
Cni = arb(1)*C
Dni = arb(1)*D

# Transposés des matrices
At = A.transpose()
Bt = B.transpose()
Ct = C.transpose()
Dt = D.transpose()

# version arb_mat
Ait = convert_arb_mat(At)
Bit = convert_arb_mat(Bt)
Cit = convert_arb_mat(Ct)
Dit = convert_arb_mat(Dt)


# Résolution de l'Equation de Lyapunov: W = A W At + B Bt

# References
#Wscipy = solve_discrete_lyapunov( A, B*Bt)
#ref = A*Wscipy*At + B*Bt
#print( "\"Précision\" résolution scipy: " + str(np.amax( Wscipy - ref)) )

Wschur = dlyap_schur( A, B*Bt)
#ref = A*Wschur*At + B*Bt
#print( "\"Précision\" résolution schur: " + str(np.amax( Wschur - ref)) )

Wslycot = dlyap_slycot( A, B*Bt)
ref = A*Wslycot*At + B*Bt
#print( "\"Précision\" résolution slycot: " + str(np.amax( Wslycot - ref)) )


# Méthodes "maison"
Wnaiv = lyap_naiv(A, B*Bt)
#ref = A*Wnaiv*At + B*Bt
#print( "\"Précision\" résolution naive: " + str(np.amax( Wnaiv - ref)) )

Wbartel = lyap_bartel_stewart( A, B*Bt)
#ref = A*Wbartel*At + B*Bt
#print( "\"Précision\" résolution bartel: " + str(np.amax( Wbartel - ref)) )

#Warb_naiv = arb_lyap_naiv(Ai, Bi*Bit)
#ref = Ai*Warb_naiv*Ait + Bi*Bit



# Calcul de la Norme L2
print ""

#Hscipy = sqrt( np.matrix.trace( C*Wscipy*Ct + D*Dt))
Hschur = sqrt( np.matrix.trace( C*Wschur*Ct + D*Dt + D*Dt))
Hslycot = sqrt( np.matrix.trace( C*Wslycot*Ct + D*Dt + D*Dt))

Hnaiv = sqrt( np.matrix.trace( C*Wnaiv*Ct + D*Dt + D*Dt))
Hbartel =  sqrt( np.matrix.trace( C*Wbartel*Ct + D*Dt + D*Dt))
#Harb_naiv = arb.sqrt( arb_trace( Ci*Warb_naiv*Cit + Di*Dit + Di*Dit))

#print ("Scipy: La norme L2 est " + str(Hscipy))
print ("Schur: La norme L2 est " + str(Hschur))
print ("Slycot: La norme L2 est " + str(Hslycot))
print ("Naiv: La norme L2 est " + str(Hnaiv))
print ("Bartel: La norme L2 est " + str(Hbartel))
#print ("Arb_naiv: La norme L2 est:" + str(Harb_naiv))


# Sauvegarde du systême testé
print ""
t = str( raw_input("Save in (default/\"nomfichier\"):"))
if (t=="default" or t=="d" or t==""):
	np.savez("log", A, B, C, D)
else :
	np.savez(t, A, B, C, D)
