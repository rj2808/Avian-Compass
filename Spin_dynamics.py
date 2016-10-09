import qutip as qt
import math as mt
import numpy as np
import matplotlib.pyplot as plt
__author__ = 'Rakshit Jain'
import time
import Hamiltonian_Eigenvalues
# We will be using the singlet and triplet dynamics as calculated by Ritz et. al. 2000
# See Hamiltonian_Eigenvalues.py for more info

def Singlet_yield(theta, rate):
	I2_x = qt.jmat(1, 'x')
	I2_y = qt.jmat(1, 'y')
	I2_z = qt.jmat(1, 'z')
	# Defining Spin operators accoding to dimensions
	Sx_aq 	 = qt.tensor(qt.sigmax(), qt.identity(288))
	Sy_aq 	 = qt.tensor(qt.sigmay(), qt.identity(288))
	Sz_aq    = qt.tensor(qt.sigmaz(), qt.identity(288))
	Sx_bq 	 = qt.tensor(qt.sigmax(), qt.identity(192))
	Sy_bq 	 = qt.tensor(qt.sigmay(), qt.identity(192))
	Sz_bq    = qt.tensor(qt.sigmaz(), qt.identity(192))
	Sx_a  	 = Sx_aq.data
	Sy_a	 = Sy_aq.data
	Sz_a     = Sz_aq.data
	Sx_b 	 = Sx_bq.data
	Sy_b 	 = Sy_bq.data
	Sz_b     = Sz_bq.data
	S_a 	 = np.array([Sx_a, Sy_a, Sz_a])
	S_b      = np.array([Sx_b, Sy_b, Sz_b])




	(Ha_eigen_theta, Hb_eigen_theta) = Hamiltonian_Eigenvalues.Hamiltonian_eigenvalues(theta, 0)

	g_a1 = 0.0
	for p in range(0,3):
		for q in range(0,3):
			(rolpA,colpA) = S_a[p].nonzero()
			(rolqA,colqA) = S_a[q].nonzero()
			(rolpB,colpB) = S_b[p].nonzero()
			(rolqB,colqB) = S_b[q].nonzero()
			A = list(set(zip(rolpA, colpA)).intersection(zip(rolqA, colqA)))
			B = list(set(zip(rolpB, colpB)).intersection(zip(rolqB, colqB)))
			#print((10.0*p + q)/100)
			for n,m in A:
				#print (str(n/576.0*100.0))        		
				for r,s, in B:
					wa_mn = Ha_eigen_theta[m] - Ha_eigen_theta[n]						# Frequencies as deined in timmel et. al. 1998
					wb_rs = Hb_eigen_theta[r] - Hb_eigen_theta[s]
					g_a1 +=  S_a[p][n, m]*S_a[q][m, n]*S_b[p][r,s]*S_b[q][s, r]* (rate**2/(rate**2 + (wa_mn - wb_rs)**2))/4					

	singletyield = (g_a1/(576*384)) + .25
	return singletyield

A = np.zeros(6, dtype = np.complex128)
i1 = np.zeros(6)
for i in range(0,6):
	A[i] = Singlet_yield((i*90.0/5.0), 10**2)
	i1[i]	 =  i*90.0/5.0
	print (i1)
	print (A)
plt.plot(i1, np.absolute(A), 'ro')
print (np.absolute(A))
print (i1)
plt.show()
plt.pause(2**31-1)


#B = Singlet_yield(0, 10**4) 
#np.savetxt('yield0_4.txt', [B])
#C = Singlet_yield(60, 10**4) 
#np.savetxt('yield60_4.txt', [C])
#D = Singlet_yield(90, 10**4) 
#np.savetxt('yield90_4.txt', [D])












