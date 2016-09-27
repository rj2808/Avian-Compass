import qutip as qt
import math as mt
import numpy as np
import matplotlib.pyplot as plt
__author__ = 'Rakshit Jain'
# We will be using the singlet and triplet dynamics as calculated by Ritz et. al. 2000
# See Hamiltonian_Eigenvalues.py for more info
I2_x = qt.jmat(1, 'x')
I2_y = qt.jmat(1, 'y')
I2_z = qt.jmat(1, 'z')
# Defining Spin operators accoding to dimensions
Sx_a 	 = qt.tensor(qt.sigmax(), qt.identity(288))
Sy_a 	 = qt.tensor(qt.sigmay(), qt.identity(288))
Sz_a     = qt.tensor(qt.sigmaz(), qt.identity(288))
Sx_b 	 = qt.tensor(qt.sigmax(), qt.identity(192))
Sy_b 	 = qt.tensor(qt.sigmay(), qt.identity(192))
Sz_b     = qt.tensor(qt.sigmaz(), qt.identity(192))
g_a		 = np.zeros((576, 576), dtype = np.complex128)
g_b		 = np.zeros((384, 384), dtype = np.complex128)
# Defining various g_jk values for the calculation of singlet yields
for i in range(0, 576):
	for j in range(0, 576):
			g_a[i, j] += Sx_a[i, j]*Sx_a[j, i] + Sx_a[i, j]*Sy_a[j, i] + Sx_a[i, j]*Sz_a[j, i] 
			g_a[i, j] += Sy_a[i, j]*Sx_a[j, i] + Sy_a[i, j]*Sy_a[j, i] + Sy_a[i, j]*Sz_a[j, i]
			g_a[i, j] += Sz_a[i, j]*Sx_a[j, i] + Sz_a[i, j]*Sy_a[j, i] + Sz_a[i, j]*Sz_a[j, i]
	print(chr(27) + "[2J")
	print('g_a')
	print(i/576*100)
for i in range(0, 384):
	for j in range(0, 384):
			g_b[i, j] += Sx_b[i, j]*Sx_b[j, i] + Sx_b[i, j]*Sy_b[j, i] + Sx_b[i, j]*Sz_b[j, i] 
			g_b[i, j] += Sy_b[i, j]*Sx_b[j, i] + Sy_b[i, j]*Sy_b[j, i] + Sy_b[i, j]*Sz_b[j, i]
			g_b[i, j] += Sz_b[i, j]*Sx_b[j, i] + Sz_b[i, j]*Sy_b[j, i] + Sz_b[i, j]*Sz_b[j, i]
	print(chr(27) + "[2J")
	print('g_b')
	print(i/384*100)
def Singlet_yield(theta, rate):
	# This function calculates singlet yield according to Cintolesi et. al. 2003
	# Loading the file regarding eigenvalues
	opened_fileA =  str(theta*10 + 1)
	opened_fileB = 	str(theta*10 + 1)
	#Ha_eigen_theta = np.loadtxt('%s.txt' % opened_fileA)
	#Hb_eigen_theta = np.loadtxt('%s.txt' % opened_fileB)
	Ha_eigen_theta = np.full((576), 1, dtype = np.int)
	Hb_eigen_theta = np.full((384), 1, dtype = np.int)
	summation = 0
	# Now putting everything into the singlet relation
	for m in range(0, 576):
		for n in range(0, 576):
			for r in range(0, 384):
				for s in range(0, 384):
					# Summation is the summation term in the
					wa_mn = Ha_eigen_theta[m] - Ha_eigen_theta[n]						# Frequencies as deined in timmel et. al. 1998
					wb_rs = Hb_eigen_theta[r] - Hb_eigen_theta[s]						
					summation += (g_a[n,m] * g_b[r, s]) * (rate**2/(rate**2 + (wa_mn - wb_rs)**2))
		print(chr(27) + "[2J")
		print(m/576*100)


	singletyield = (summation/(576*384))
	return singletyield
A = Singlet_yield(30, 10**4) 
np.savetxt('yield30.txt', [A])
B = Singlet_yield(0, 10**4) 
np.savetxt('yield0.txt', [B])
C = Singlet_yield(60, 10**4) 
np.savetxt('yield0.txt', [C])
D = Singlet_yield(90, 10**4) 
np.savetxt('yield0.txt', [D])












