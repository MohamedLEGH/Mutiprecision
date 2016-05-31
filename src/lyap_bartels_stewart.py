# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED


import numpy as np
from scipy import linalg

from bartelsStewart_functions import lyap_naiv_con, lyap_naiv_con2, part, divise_en_4

def lyap_bartel_stewart(A,BBt):

	A = np.matrix(A)
	BBt = np.matrix(BBt)
	
	# Forme continue de l'équation
	n, m = A.shape
	I = np.eye(n)
	Acon = (A-I)*np.linalg.inv(A+I)
	Acon_t = Acon.transpose()
	D = 2*np.linalg.inv(A.transpose() + I)*BBt*np.linalg.inv(A+I)
	
	# Décomposition de schur de Acon et partition de T
	[T,Q] = linalg.schur(Acon)
	r, t = part(T) #r nombre de blocs diagonaux, nk leurs tailles resp.
	
	# Initialisation de Y et C
	Y = np.matrix(np.empty( (n, m),order='C'))
	C = Q.transpose()*D*Q

	# Traitement 
	for k in range(r-1,0,-1):
		
		#Partition de la matrice C
		[C11,C12,C21,C22] = divise_en_4( T[0:m,0:m], C[0:m,0:m])
		m = m - t[k]
		
		R22 = T[m:m+t[k],m:m+t[k]]
		R12 = T[0:m,m:m+t[k]]
		R11 = T[0:m,0:m]
	
		# Eq 5.9d :  calcul de Z22
		if( t[k] == 1) :
			Z22 = C22/(2*R22)
		else :
			Z22 = lyap_naiv_con( R22, C22)
		
		# 5.12a et 5.12b :  Calcul de Z12 et Z21
		C12p = C12 - R12*Z22
		
		mj = m
		
		# Initialisation de Z12 et Z21
		Z12 = np.mat( np.zeros((m,t[k])) )
		
		for j in range(k-1,-1,-1):
			# Calcul de Tjj
			mj = mj - t[j]
			Tjj = T[ mj:mj+t[j], mj:mj+t[j]]		
			
			# Selection du bon bloc
			C12j = C12p[mj:mj+t[j],:]
	
			# 5.15 : Calcul de Z12 et Z21
			Z12[mj:mj+t[j],:] = lyap_naiv_con2(Tjj,R22,C12j)

		# stockage des valeurs
		Z21 = Z12.transpose()
		
		Y[0:m,m:m+t[k]] = Z12
		Y[m:m+t[k],0:m] = Z21
		Y[m:m+t[k],m:m+t[k]] = Z22
		
		C = C11- R12*Z21 - Z12*R12.transpose()
	
	return Q*Y*Q.transpose()
		
		
		
