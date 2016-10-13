import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import math as mt
__author__ = 'Rakshit Jain'

I2_x = qt.jmat(1, 'x')/2
I2_y = qt.jmat(1, 'y')/2
I2_z = qt.jmat(1, 'z')/2
I1_x = qt.sigmax()/2
I1_y = qt.sigmay()/2
I1_z = qt.sigmaz()/2
Sx 	 = qt.sigmax()/2
Sy 	 = qt.sigmay()/2
Sz   = qt.sigmaz()/2


def Hamiltonian_toymodel(theta, phi):

## Defining Constants
	omega = 1.4 * 10**6
	m_b = 5.7883818066*10**-5
	par=m_b*2.417990504024*10**11 # s-1/mT
	a1 = np.array([-.0989, -.0989, 1.7569]) * par # /mT
	a2 = np.array([0.0,0.0,1.0812]) * par

	H_a_zeeman_nt = omega * (Sx * mt.sin(theta) *  mt.cos(phi) + Sy *  mt.sin(theta) * mt. sin(phi) + Sz * mt.cos(theta))
	Ha_zeeman  = qt.tensor(H_a_zeeman_nt, qt.identity(3))

	H_a_hfi = a1[0] *  qt.tensor(Sx, I2_x) + a1[1] * qt.tensor(Sy, I2_y) + a1[2] * qt.tensor(Sz, I2_z)

	H_a = Ha_zeeman + H_a_hfi
	H_a_eigen, H_a_eigenstates = H_a.eigenstates()


	H_b_zeeman_nt = omega * (Sx * mt.sin(theta) *  mt.cos(phi) + Sy *  mt.sin(theta) * mt. sin(phi) + Sz * mt.cos(theta))
	Hb_zeeman  = qt.tensor(H_b_zeeman_nt, qt.identity(3))

	H_b_hfi = a2[0] *  qt.tensor(Sx, I2_x) + a2[1] * qt.tensor(Sy, I2_y) + a2[2] * qt.tensor(Sz, I2_z)

	H_b = Hb_zeeman + H_b_hfi
	H_b_eigen, H_b_eigenstates = H_b.eigenstates()

	return (H_a_eigen, H_a_eigenstates, H_b_eigen, H_b_eigenstates)

def Singlet_yield(theta, rate) :
	(Ha_eigen_theta,H_a_eigenstates, Hb_eigen_theta, H_b_eigenstates) = Hamiltonian_toymodel(theta,0)
	Sx_aq 	 = qt.tensor(Sx, qt.identity(3)).transform(H_a_eigenstates)
	Sy_aq 	 = qt.tensor(Sy, qt.identity(3)).transform(H_a_eigenstates)
	Sz_aq    = qt.tensor(Sz, qt.identity(3)).transform(H_a_eigenstates)
	Sx_bq 	 = qt.tensor(Sx, qt.identity(3)).transform(H_b_eigenstates)
	Sy_bq 	 = qt.tensor(Sy, qt.identity(3)).transform(H_b_eigenstates)
	Sz_bq    = qt.tensor(Sz, qt.identity(3)).transform(H_b_eigenstates)
	Sx_a  	 = Sx_aq.data
	Sy_a	 = Sy_aq.data
	Sz_a     = Sz_aq.data
	Sx_b 	 = Sx_bq.data
	Sy_b 	 = Sy_bq.data
	Sz_b     = Sz_bq.data
	S_a 	 = np.array([Sx_a, Sy_a, Sz_a])
	S_b      = np.array([Sx_b, Sy_b, Sz_b])
	g_a1 = 0.0
	for p in range(0,3):
		for q in range(0,3):
			(rolpA,colpA) = S_a[p].nonzero()
			(rolqA,colqA) = S_a[q].nonzero()
			(rolpB,colpB) = S_b[p].nonzero()
			(rolqB,colqB) = S_b[q].nonzero()
			A = list(set(zip(rolpA, colpA)).intersection(zip(rolqA, colqA))) # Techniques to iterate over sparse matrices
			B = list(set(zip(rolpB, colpB)).intersection(zip(rolqB, colqB))) # Taking Intersction over non-zero values
			#print((10.0*p + q)/100)
			for n,m in A:
				#print (str(n/576.0*100.0))        		
				for r,s in B:
					wa_mn = Ha_eigen_theta[m] - Ha_eigen_theta[n]						# Frequencies as deined in timmel et. al. 1998
					wb_rs = Hb_eigen_theta[r] - Hb_eigen_theta[s]
					g_a1 +=  S_a[p][n, m]*S_a[q][m, n]*S_b[p][r,s]*S_b[q][s, r]* (rate**2/(rate**2 + (wa_mn - wb_rs)**2))				

	singletyield = (g_a1/9.0) + .25
	return singletyield

A = np.zeros(101, dtype = np.complex128)
i1 = np.zeros(101)
for i in range(0,101):
	A[i] =  Singlet_yield((i*mt.pi/100.0), 10**4)
	i1[i]	 =  i*180.0/100.0
	#print (i1)
	#print (A)
plt.plot(i1, np.absolute(A), )
print (np.absolute(A))
print (i1)
plt.show()
plt.pause(2**31-1)














