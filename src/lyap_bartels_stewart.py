# -*- coding: utf-8 -*-

# ALOUI Driss
# GAUDELET Lucas
# LEGHEBARA MOHAMED


import numpy as np
import scipy as sp

from bartelsStewart_functions import lyap_naiv_con, lyap_naiv_con2, part, divise_en_4

def lyap_bartel_stewart(A,BBt):
	
	# Forme continue de l'équation
	n, m = A.shape
	I = np.eye(n)
	Acon = (A-I)*np.linalg.inv(A+I)
	Acon_t = Acon.transpose()
	D = 2*np.linalg.inv(A.transpose() + I)*BBt*np.linalg.inv(A+I)
	
	# Décomposition de schur de Acon et partition de T
	[T,Q] = sp.linalg.schur(Acon)
	r, t = part(T) #r nombre de blocs diagonaux, nk leurs tailles resp.
	
	# Initialisation de Y et C
	Y = np.empty( (n, m),order='C')
	C = Q.transpose()*D*Q

	# Traitement 
	for k in range(r-1,0,-1):
		print ( "\t\tboucle for k:\t"+str(k)+"\n")
		print ("T=\n"+str(T[0:m,0:m]) )
		#Partition de la matrice C
		[C11,C12,C21,C22] = divise_en_4( T[0:m,0:m], C[0:m,0:m])
		m = m - t[k]
		
		
		# tests
		print ("C=\n" + str(C) )
		print ""
		print ("C11=\n" + str(C11) )
		print ""
		print ("C12=\n" + str(C12) )
		print ""
		print ("C21=\n" + str(C21) )
		print ""
		print ("C22=\n" + str(C22) )
		print ""
		
		R22 = T[m:m+t[k],m:m+t[k]]
		R12 = T[0:m,m:m+t[k]]
		R11 = T[0:m,0:m]
		
		# test
		print ("R22=\n" + str(R22) )
		print ""
		print ("R12=\n" + str(R12) )
		print ""
		print ("R11=\n" + str(R11) )
		print ""
		
		# Eq 5.9d :  calcul de Z22
		if( t[k] == 1) :
			Z22 = C22/(2*R22)
		else :
			Z22 = lyap_naiv_con( R22, C22)
		
		# 5.12a et 5.12b :  Calcul de Z12 et Z21
#		C12p = C12 - R12*Z22
		C12p = R12*Z22
		C12p = C12 - C12p
		C21p = C21.transpose() - R12*Z22.transpose()
		
		print ("C12p=\n" + str(C12p) )
		print ""
		print ("C21p=\n" + str(C21p) )
		print ""
		
		mj = m
		
		# Initialisation de Z12 et Z21
		Z12 = np.zeros((m,t[k])) 
		Z21t = np.zeros((t[k],m)) 
		
		for j in range(k-1,1,-1):
			print ( "\t\tboucle for j:\t"+str(j)+"\n" )
			# Calcul de Tjj
			mj = mj - t[j]
			Tjj = T[ mj:mj+t[j], mj:mj+t[j]]
			print ("Tjj=\n" + str(Tjj) )
			print ""
			
			# Selection du bon bloc
			C12j = C12p[0:mj+1,mj:mj+t[j]]
			C21j = C21p[mj:mj+t[j],0:mj+1]
			
			# test
			print ("C12j=\n" + str(C12j	) )
			print ""
			print ("C21jt=\n" + str(C21j) )
			print ""
			
			# 5.15 : Calcul de Z12 et Z21
			Z12[0:mj,mj:mj+t[j]] = lyap_naiv_con2(Tjj,R22,C12j)
			#Z12[] = lyap_naiv_con2(Tjj,R22,C12j)
			Z21t[mj:mj+t[j],0:mj] = lyap_naiv_con2(Tjj,R22,C21j)
			#Z21t[j] = lyap_naiv_con2(Tjj,R22,C21jt)
		
		# test
		print ("Z12=\n" + str(Z12) )
		print ""
		print ("Z21t=\n" + str(Z21t) )
		print ""
		
		# stockage des valeurs
		#Z21 = Z21t.transpose()
		Y[0:m,m:m+t[k]] = Z12
		Y[m:m+t[k],0:m] = Z21t
		#Y[1:m,m:m+t[k]-1] = Z21
		Y[m:m+t[k],m:m+t[k]] = Z22		
		
		Ctmp1 = R12*Z21t.transpose()
		Ctmp2 = Z12*R12.transpose()
		C = C11 - Ctmp1
		C = C - Ctmp2
		
		#C = C11 - R12*Z21t - Z12*R12.transpose()
		
	return Q*Y*Q.transpose()
		
		
		
