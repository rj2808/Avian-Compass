import qutip as qt
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import Hamiltonian_Eigenvalues

# We will be using the singlet and triplet dynamics as calculated by Ritz et. al. 2000
# See Hamiltonian_Eigenvalues.py for more info
def Singlet_yield():
	# This function calculates singlet yield according to Cintolesi et. al. 2003


	# Defining basic spin operators
	I2_x = qt.jmat(1, 'x')
	I2_y = qt.jmat(1, 'y')
	I2_z = qt.jmat(1, 'z')
	Sx 	 = qt.sigmax()
	Sy 	 = qt.sigmay()
	Sz   = qt.sigmaz()
	g 	 = np.zeros([3,3])
	# Defining various g_jk values for the calculation of singlet yields
		for i in range (0, 3):
			for j in range(0, 3):
				g[i, j]  = Sx[i, j]*Sx[j, i] + Sx[i, j]*Sy[j, i] + Sx[i, j]*Sz[j, i] 
				g[i, j] += Sy[i, j]*Sx[j, i] + Sy[i, j]*Sy[j, i] + Sy[i, j]*Sz[j, i]
				g[i, j] += Sz[i, j]*Sx[j, i] + Sz[i, j]*Sy[j, i] + Sz[i, j]*Sz[j, i]




