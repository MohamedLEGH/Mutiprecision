# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED
# MAIN 3


# Imports
import numpy as np
import mpmath
from flint import arb, arb_mat, acb, acb_mat

def convert_mpmath_mat( M):
	return mpmath.matrix([i for i in M.tolist()])

def schur_mp_to_arb( M):
	"""Computes the Schur decomposition of the matrix M
	
	The result is computed twice at different precisions using mpmath's schur function. Then the radius of the error interval is estimated to be half the difference beetween the two results.
	"""
	# Conversion de la matrice en mpmath.matrix
	A = convert_mpmath_mat( M)
	
	# Calcul en précision actuelle
	Q1, R1 = mpmath.schur( A)
	
	# Calcul en précision double
	mpmath.mp.prec *=2
	Q2, R2 = mpmath.schur( A)
	
	# Création de la matrice arb
	m, n = M.shape
	lq = []
	lr = []
	for i in range(m):
		for j in range(n):
			lq.append( acb( Q1[i,j].real, Q1[i,j].imag ) )
			lr.append( acb( R1[i,j].real, R1[i,j].imag ) )
	Q = acb_mat( m, n, lq)
	R = acb_mat( m, n, lr)
	
	mpmath.mp.prec /=2
	
	return( Q, R)







