# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED

from flint import arb, arb_mat
from math import floor

def arb_transpose(M):
	"""	Computes the transpose of the matrix M """
	n, m = arb_mat.nrows(M), arb_mat.ncols(M)
	res = arb_mat( m, n)
	
	for i in range(n):
		for j in range(m):
			res[i,j] = M[j,i]
	
	return res


def arb_kron( A, B):
	""" Computes the Kronecker Product of the matrix A and B """
	
	# Matrix sizes
	n, m = arb_mat.nrows(A), arb_mat.ncols(A)
	p, q = arb_mat.nrows(B), arb_mat.ncols(B)
	
	# Initialization of the result
	res = arb_mat( n*p, m*q)
	
	# 
	for i in range(n*m) :
		for j in range(p*q) :
			res[i,j] = A[floor(i/n),floor(j/m)]*B[i%p,j%q]
					
	return res


def arb_vect(M):
	""" Computes the vectorization of the matrix M in FORTRAN style """
	
	# Initialization
	k = 0 # vector indice
	n, m = arb_mat.nrows(M), arb_mat.ncols(M) # dimensions of matrix M
	res = arb_mat( n*m, 1) # vector size n*m filled with zeros
	
	# Vectorization
	for j in range(m):
		for i in range(n):
			res[k,0] = M[i,j]
			k=k+1
			
	return res


def arb_devect(v, n, m):
	""" Computes the devectorization of the matrix M in FORTRAN style """
	
	# Initialization
	dim = arb_mat.nrows(v)
	res = arb_mat( n, m)
	
	# Computation
	for i in range(dim):
	 	res[ i%n , floor(i/m) ] = v[i,0]
	 
	return res


def arb_trace(M):
	""" Computes the trace of the matrix M """
	# Matrix size
	n, m = arb_mat.nrows(M), arb_mat.ncols(M)
	
	# Computation
	res = arb(0)
	for i in range(n):
		res = res + M[i,i]

	return res
	
