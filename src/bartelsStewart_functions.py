# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED


import numpy as np

def part( M):
	"""
	Cette conction partitionne T une matrice triangulaire supérieure par blocs.
	
	Retourne r le nombre de blocs diagonaux et Mpart leurs tailles respectives.
	"""
	n, m = M.shape
	Mpart = []
	
	t = 1
	r = 0
	for i in range( n-1):
		if(M[i+1,i] != 0):
			t = t+1
		else:
			Mpart.append(t)
			r = r+1
			t = 1
	if(M[n-1,n-2]!=0):
		Mpart.append(2)
	else:
		Mpart.append(1)
	r = r+1
	
	return (r, Mpart)


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
	

def lyap_naiv_con(A, W):
	""" Calcule vect(X)=(IxA+AxI)⁻¹vect(W) """

	n = A.shape[1]
	
	w = W.flatten('F').transpose()
	
	In=np.eye(n)
	Ak1 = np.kron(In,A)
	Ak2 = np.kron(A,In)
	
	Ai = np.linalg.inv(Ak1 + Ak2)
	
	X = Ai*w
	
	return np.reshape(X,(n,-1),order='F') 
	

def lyap_naiv_con2(A ,B , W):
	"""	Calcule vect(X)=(IxA+BxI)⁻¹vect(W) """
	
	k = A.shape[0]
	m = B.shape[1]
		
	w = W.flatten('F').transpose()

	Im = np.eye(m)
	Ak = np.kron(Im,A)
	Ik = np.eye(k)
	Bk = np.kron(B,Ik)

	Ai = np.linalg.inv(Ak + Bk)

	X = Ai*w
	
	return np.reshape(X,(k,-1),order='F') 
	

