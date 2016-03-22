# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED

import numpy as np
from math import sqrt
from arb_mat_functions import arb_kron, arb_vect, arb_devect
from flint import arb_mat
	
def arb_lyap_naiv(A, BBt):
	"""
	Computes the numerical solution of the Lyapunov equation:
		W = A W At + B Bt
	
	"""
	# Matrix dimension
	n, m = arb_mat.nrows(A), arb_mat.ncols(A)

	# I(n^2), vect(B*Bt) and kronecker( A, A)
	I = arb_mat(n**2,m**2)+1
	b = arb_vect(BBt)
	Ak = arb_kron(A,A)
	
	# ( I(n^2) - kronecker(A,A) )⁻¹ and A⁻¹ * vect(B*Bt)
#	print( "[arb_lyap_naiv] det(I-Ak) =" + str( arb_mat.det(I-Ak)) )
	#Ai = arb_mat.solve( I-Ak, I)
	Ai = arb_mat.inv( I - Ak)
	w = Ai*b
	
	return arb_devect(w, n, m)

