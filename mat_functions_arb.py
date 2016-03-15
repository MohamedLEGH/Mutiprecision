# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED
# MAIN 3
# 08.03.16

from flint import arb_mat
from math import floor


def arb_kron( A, B):
	"""
	
	"""
	# Taille des matrices
	n, m = arb_mat.nrows(A), arb_mat.ncols(A)
	p, q = arb_mat.nrows(B), arb_mat.ncols(B)
	
	# Initialistaion du résultat
	res = arb_mat( n*p, m*q)
	
	# traitement
	for i in range(n*m) :
		for j in range(p*q) :
			res[i,j] = A[floor(i/n),floor(j/m)]*B[i%p,j%q]
					
	return res

def arb_vect(A):
	"""
	Cette fonction permet de vectorizer une matrice : On transforme une matrice n*m sous la forme d'un vecteur n*m colonne 
	dans le style FORTRAN.
	"""
	n, m = arb_mat.nrows(A), arb_mat.ncols(A)
	
	k = 0 # l'indice pour notre vecteur
	
	a = arb_mat(n*m,1) # vecteur (matrice taille(A)*1 ) rempli avec des zéros
	
	for j in range(m):
		for i in range(n):
			a[k,0] = A[i,j]
			k=k+1
	return a
	
def arb_transpose(A):
	"""
	Cette fonction renvoie la transposé de la matrice A passée en argument
	"""
	
	
