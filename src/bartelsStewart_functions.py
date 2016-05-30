# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED


import numpy as np

def part( M):
	"""
	Cette conction partitionne T une matrice triangulaire supérieure par blocs.
	
	Retourne r le nombre de blocs diagonaux et Tpart leurs tailles respectives.
	"""
	n, m = M.shape
	
	Mpart = []
	k = 0
	r = 1
	while (k<m-2) :
		r = r+1
		if( M[k+1,k] == 0 ):
			Mpart.append(1)
			k = k+1
		else :
			Mpart.append(2)
			k = k+2
	if( M[m-1,m-2] == 0 ):
		Mpart.append(1)
	else :
		Mpart.append(2)
	
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
#calcule vect(X)=(IxA+AxI)⁻¹vect(W)
	
# Matrix dimension
	N = A.shape[1]
# I(n^2), vect(W) and kronecker( I, A) kronecker( A, I)
	Inn=np.eye(N)
	print ("Inn=\n" + str(Inn) )
	print ""
	w = W.flatten('F')
	w = w.transpose()
	Ak1 = np.kron(Inn,A)
	Ak2 = np.kron(A,Inn)
# ( kronecker(I,A) + kronecker(A,I) )⁻¹ 
	Ai = np.linalg.inv(Ak1 + Ak2)
	print ("Ai=\n" + str(Ai) )
	print ""
	print ("w=\n" + str(w) )
	X = Ai*w
	return np.reshape(X,(N,-1),order='F') 
	

def lyap_naiv_con2(A ,B , W):
#calcule vect(X)=(IxA+BxI)⁻¹vect(W)
	
# Matrix dimension
	k = A.shape[0]
	m = B.shape[1]
# I(n^2), vect(W) and kronecker( I, A) kronecker( B, I)
	Ik = np.eye(k)
	Im = np.eye(m)
	w = W.flatten('F')
	print ("Ik=\n" + str(Ik) )
	print ""
	print ("Im=\n" + str(Im) )
	print ""
	print ("A=\n" + str(A) )
	print ""
	print ("B=\n" + str(B) )
	print ""
	
	Ak = np.kron(Im,A)
	Bk = np.kron(B,Ik)
# ( kronecker(I,A) + kronecker(A,I) )⁻¹ 
	Ai = np.linalg.inv(Ak + Bk)
	print ("Ai=\n" + str(Ai) )
	print ""
	print ("w=\n" + str(w) )
	X = Ai*w
	return np.reshape(X,(k,-1),order='F')
	
