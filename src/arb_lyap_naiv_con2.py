# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED


import numpy as np
from math import sqrt
from arb_mat_functions import arb_kron, arb_vect, arb_devect
from flint import arb_mat

def arb_lyap_naiv_con2(A ,B , W):
#calcule vect(X)=(IxA+BxI)⁻¹vect(W)
	
	
# Matrix dimension
	n, m = arb_mat.nrows(A), arb_mat.ncols(A)
# I(n^2), vect(W) and kronecker( I, A) kronecker( B, I)
	I = arb_mat(n**2,m**2,1)
	w = arb_vect(W)
	Ak1 = arb_kron(I,A)
	Ak2 = arb_kron(B,I)
# ( kronecker(I,A) + kronecker(A,I) )⁻¹ 
	Ai = arb_mat.inv( Ak1 + Ak2)
	X = Ai*w
return arb_devect(X, n, m&)
