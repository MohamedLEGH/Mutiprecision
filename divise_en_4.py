# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED


import numpy as np
from math import sqrt
from flint import arb_mat
	
def divise_en_4(T,A):
#divise en 4 A
# A1 A2
# A3 A4

# Matrix dimension
	N = A.shape[1]
	M = T.shape[1]
	
#test pour taille bloc
	i=M
	if (T[N-1,N-2] != 0) :
		i=M-1
	A1 = A[0:i-1,0:i-1]
	A2 = A[0:i-1,i-1:M]
	A3 = A[i-1:M,0:i-1]
	A4 = A[i-1:M,i-1:M]
	return [A1,A2,A3,A4]


