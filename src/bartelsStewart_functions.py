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
#calcule vect(X)=(IxA+AxI)⁻¹vect(W)
	
# Matrix dimension
	n = A.shape[1]
# I(n^2), vect(W) and kronecker( I, A) kronecker( A, I)
	In=np.eye(n)
	w = W.flatten('F')
	w = w.transpose()
	Ak1 = np.kron(In,A)
	Ak2 = np.kron(A,In)
	print("\tAk1 = kron(In,A)")
	print("\tAk2 = kron(A,In)")
# ( kronecker(I,A) + kronecker(A,I) )⁻¹ 
	Ai = np.linalg.inv(Ak1 + Ak2)
	print ("\tAi=inv(Ak1+Ak2)\n" + str(Ai) )
	print ""
	print ("\tw=\n" + str(w) )
	X = Ai*w
	res = np.reshape(X,(n,-1),order='F') 
	print ""
	print ("\tres=\n"+str(res) )
	print ""
	return res
	

def lyap_naiv_con2(A ,B , W):
#calcule vect(X)=(IxA+BxI)⁻¹vect(W)
	
# Matrix dimension
	k = A.shape[0]
	m = B.shape[1]
# I(n^2), vect(W) and kronecker( I, A) kronecker( B, I)
	Ik = np.eye(k)
	Im = np.eye(m)
	w = W.flatten('F').transpose()
	print ("\tA=\n" + str(A) )
	print ""
	print ("\tB=\n" + str(B) )
	print ""
	
	print( "\tAk = kron(In,A)" )
	print( "\tBk= kron(B,Ik)" )
	Ak = np.kron(Im,A)
	Bk = np.kron(B,Ik)
# ( kronecker(I,A) + kronecker(A,I) )⁻¹ 
	Ai = np.linalg.inv(Ak + Bk)
	print ("\tAi=inv(Ak+Bk)\n" + str(Ai) )
	print ""
	print ("\tw=W.flatten\n" + str(w) )
	X = Ai*w
	res = np.reshape(X,(k,-1),order='F') 
	print ""
	print ("\tres=\n"+str(res) )
	print ""
	return res
	
