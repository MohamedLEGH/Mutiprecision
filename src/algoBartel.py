from math import sqrt
import numpy as np
from scipy import linalg
# Pas encore fini
def lyap_bartel(A,BBt):
	Ay = (A-I)*np.linalg.inv(A+I)
	D = 2*np.linalg.inv(At + I)*BBt*np.linalg.inv(A+I)
	Ayt = Ay.transpose()
	Y= np.empty(shape,order='C')
	[Q,T] = linalg.schur(Ay)
	C = Q.transpose()*D*Q
	m,n=A.shape
	r=0
	z=0
# r c'est le nombre de block de la matrice T
	for k in range(r,1,-1):
	#For k de r à 1 pas =-1
		m=m-1
		R22=T[k,k]
		R12=T[1:k-1,k]
		R11 = T - R22 - R12
		[C11,C12,C21,C22] = divise_en_4(C)
		nk,mk = R22.shape
		Z22=C22/(2*R22)
		C12p = C12 - R12*Z22
		C21p = C21t - R12*Z22.transpose()
		for j in range(k-1,1,-1):
			C12p(j) = C12(j) - np.sum(R11[j,1:k-1]*C12(1:k+1))
			C21p(j) = C21(j) - np.sum(R11[j,1:k-1]*C21(1:k+1))
			Z12(j) = naivContinu2(T[j,j],R22,C12p(j))
			Z21t(j) = naivContinu2(T[j,j],R22,C21p(j))
		Z21 = Z12t.transpose()
		Y(1:m,m+(1:nk)) = Z12
		Y(m+(1:nk),1:m) = Z21
		Y(m+(1:nk),m+(1:nk)) = Z22		
		C = C11 - R12*Z21 - Z12*R12.transpose()
	return Q*Y*Qt 	

