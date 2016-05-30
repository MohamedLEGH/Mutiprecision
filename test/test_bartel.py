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

from lyap_bartels_stewart import lyap_bartel_stewart

# Initialisation et choix du systême 

ctx.prec = 300

t = str( raw_input("Systême de matrices (random/last/\"nomfichier\"):"))

# Génération aléatoire
if( t == "random" or t == "r" or t=="") :
	n,p,q = 5, 5, 5
	A, B, C, D = random_dSS( n, p, q)
	
	print("Rayon spectral de A: " + str(np.amax( np.absolute( np.linalg.eigvals( A)))) + "\n")
	
	# Sauvegarde du systême testé
	print ""
	t = str( raw_input("Save in (default/\"nomfichier\"):"))
	if (t=="default" or t=="d" or t==""):
		np.savez("log", A, B, C, D)
	else :
		np.savez(t, A, B, C, D)

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



# Transposés des matrices
At = A.transpose()
Bt = B.transpose()
Ct = C.transpose()
Dt = D.transpose()

# Résolution de l'Equation de Lyapunov: W = A W At + B Bt

Wbartel = lyap_bartel_stewart( A, B*Bt)
ref = A*Wbartel*At + B*Bt
print( "\"Précision\" résolution bartel: " + str(np.amax( Wbartel - ref)) )

# Calcul de la Norme L2
print ""
print ("W = " + str(Wbartel) )

Hbartel =  sqrt( np.matrix.trace( C*Wbartel*Ct + D*Dt + D*Dt))

print ("Bartel: La norme L2 est " + str(Hbartel))
