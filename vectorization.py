# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA Mohamed
# MAIN 3
# 15.03.16

import python-flint 

def vec(A):
"""
Ce programme permet de vectorizer une matrice : On transforme une matrice n*m sous la forme d'un vecteur n*m colonne

"""
	k=0 # l'indice pour notre vecteur
	a=arb_mat(1,A.ncols*A.nrows) # vecteur (matrice 1xtaille(A) ) rempli par des zéros 
	for j in A.ncols:
		for i in A.nrows:
			a[k]=A[i,j]
			k=k+1
	return a

