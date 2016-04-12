# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED
# MAIN 3
# 09.02.16

from math import sqrt
import numpy as np
	
def lyap_naiv(A, BBt):
	"""Méthode naive de résolution de l'équation discrete de lyapunov:
	
		vec(W) = ( I**2 - kronecker( A, A) )^(-1) * vec(B*Bt)
		
	This was inspired from "Lecture notes in numerical linear algebra" by Elias Jarlebring, 
	"""
	
	# Initialisation
	N = A.shape[1]
	Inn = np.eye(N**2) # Matrice identité de taille n^2
	Akro = np.kron(A,A) # Produit de kroneker entre A et lui même
	Bvec = BBt.flatten('F') # On transforme B*Bt en un vecteur, selon le style FORTRAN
	
	# Calcul de l'inverse de ( Id(n^2) - kronecker( A, A) )
	AIinf = np.linalg.inv(Inn-Akro) 
	Bvec = Bvec.transpose() # Transformation du vecteur ligne en colonne
	Wvec = AIinf*Bvec # Calcul de la vectorization de W
	return np.reshape(Wvec,(N,-1),order='F') # Dévectorization de W, toujours selon le style FORTRAN

