# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED


import numpy as np
from math import sqrt
from flint import arb_mat
	
def lyap_naiv_con(A, W):
#calcule vect(X)=(IxA+AxI)⁻¹vect(W)
	
	
# Matrix dimension
	N = A.shape[1]
# I(n^2), vect(W) and kronecker( I, A) kronecker( A, I)
	Inn=np.eye(N**2)
	w = W.flatten('F')
	Ak1 = np.kron(I,A)
	Ak2 = np.kron(A,I)
# ( kronecker(I,A) + kronecker(A,I) )⁻¹ 
	Ai = np.linalg.inv(Ak1 + Ak2)
	X = Ai*w
	return np.reshape(X,(N,-1),order='F') 
