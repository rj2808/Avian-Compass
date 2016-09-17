# For Details see the supplements of Hancock et.al. 2016 Quantum Needle of Avaian Compass
__author__ = 'Rakshit Jain'
import qutip as qt
import math as mt
import numpy as np
import scipy as sc
# Basic 3*3 Tensor arays for Hyperfine interactions
# For FAD Radical
def Hamiltonian_eigenvalues(theta1, theta2, phi1, phi2): # Here o refers to the angle
	# Defining Spin Operators
	# I_2 denotes the spin operators for S = 1 (Nitrogen Atoms)
	# I_1 denotes the spin operators for S = 1/2.0 (Hydrogen Atoms)
	I2_x = qt.jmat(1, 'x')
	I2_y = qt.jmat(1, 'y')
	I2_z = qt.jmat(1, 'z')
	I1_x = qt.sigmax()
	I1_y = qt.sigmay()
	I1_z = qt.sigmaz()
	Sx 	 = qt.sigmax()
	Sy 	 = qt.sigmay()
	Sz   = qt.sigmaz()
	# Defining Constants here
	omega = 1.4 * 10**9				#Larmor Frequency in Hz for 50 microT
	g = 2  							# Electron Parameter
	# Isotropic Hyperfine Coupling Vector
	#For FAD Radical
	a_iso1 = np.array([.5233, .1887, -0.3872, 0.4399, 0.4070])
	T_N5n = np.array([[1.2336, 0, 0],[ 0, -0.6101, 0],[0, 0, -0.6234]])
	T_N5 = qt.Qobj(T_N5n) 
	T_N10n = np.array([[0.4159, 0, 0], [0, -0.2031, 0], [0, 0, 0.2128]])
	T_N10 = qt.Qobj(T_N10n)
	T_H6n = np.array([[0.1896, 0, 0], [0, -0.0464, 0], [0, 0, -0.1432]])
	T_H6 = qt.Qobj(T_H6n)
	T_H8n = np.zeros((3,3))
	T_H8 = qt.Qobj(T_H8n)
	T_Hbn = np.zeros((3,3))
	T_Hb = qt.Qobj(T_Hbn)
	# Similarly for TrpH Radical
	a_iso2 = np.array([.3215, -.5983, -.2780, -.488, -.0400, -.3636, 1.6046])
	T_N1n = np.array([[.7596, 0, 0], [0, -.3745, 0], [0, 0, -.3851]])
	T_N1 = qt.Qobj(T_N1n)
	T_H1n = np.array([[.5914, 0, 0], [0, -.1071, 0], [0, 0, -.4853]])
	T_H1 = qt.Qobj(T_H1n)
	T_H2n = np.array([[.2855, 0 , 0], [0, -.0919, 0], [0, 0, -.1936]])
	T_H2 = qt.Qobj(T_H2n)
	T_H4n = np.array([[.3001, 0, 0], [0, -.0480, 0], [0, 0, -.2520]])
	T_H4 = qt.Qobj(T_H4n)
	T_H5n = np.array([[-.0632, 0, 0], [0, -.0616, 0], [0, 0, .1248]])
	T_H5 = qt.Qobj(T_H5n)
	T_H7n = np.array([[.2540, 0, 0],[0, -.0594, 0], [0, 0, -.1945]])
	T_H7 =  qt.Qobj(T_H7n)
	T_Hb1n = np.array([[0.1521, 0, 0], [0, -.0456, 0], [0, 0, -.1065]])
	T_Hb1 = qt.Qobj(T_Hb1n)

	# Now since we have defined all the tensors required till now, Writing down the hamiltonian
	# H_a denotes the hamiltonian for FAD Radical
	# H_b denotes the Hamiltonian for TrpH radical

	Ha_zeeman_nt = omega * (Sx * mt.sin(theta1) *  mt.cos(phi1) + Sy *  mt.sin(theta1) * mt. sin(phi1) + Sz * mt.cos(theta1))			# Hamiltonian without any tensor product
	Ha_zeeman = qt.tensor(Ha_zeeman_nt, qt.identity(3**7))			# Compounding different states 
	# Using the anisotropic hyperfine tensors the hamiltonian comes down to this
	Ha_hfi_N5  = qt.tensor(a_iso1[0] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
	 				+ qt.tensor(Sx, T_N5 * I2_x) + qt.tensor(Sy, T_N5 * I2_y) + qt.tensor(Sz, T_N5 * I2_z), qt.identity(3**6))

	Ha_hfi_N10 = qt.tensor( qt.tensor(3**1), a_iso1[1] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_N10 * I2_x) + qt.tensor(Sy, T_N10 * I2_y) + qt.tensor(Sz, T_N10 * I2_z), qt.identity(3**5))

	Ha_hfi_H6  = qt.tensor( qt.tensor(3**2), (a_iso1[2] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z))) 
					+  qt.tensor(Sx, T_H6 * I2_x) + qt.tensor(Sy, T_H6 * I2_y) + qt.tensor(Sz, T_H6 * I2_z), qt.identity(3**4))

	Ha_hfi_H8  = qt.tensor( qt.tensor(3**3), a_iso1[3] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_H8 * I2_x) + qt.tensor(Sy, T_H8 * I2_y) + qt.tensor(Sz, T_H8 * I2_z), qt.identity(3**3))

	Ha_hfi_b1  = qt.tensor( qt.tensor(3**4), a_iso1[4] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_Hb * I2_x) + qt.tensor(Sy, T_Hb * I2_y) + qt.tensor(Sz, T_Hb * I2_z), qt.identity(3**2))

	Ha_hfi_b2  = qt.tensor( qt.tensor(3**5), a_iso1[4] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_Hb * I2_x) + qt.tensor(Sy, T_Hb * I2_y) + qt.tensor(Sz, T_Hb * I2_z), qt.identity(3**1))

	Ha_hfi_b3  = qt.tensor( qt.tensor(3**6), a_iso1[4] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_Hb * I2_x) + qt.tensor(Sy, T_Hb * I2_y) + qt.tensor(Sz, T_Hb * I2_z), qt.identity(3**0))



	H_a = Ha_zeeman + Ha_hfi_N5 + Ha_hfi_N10 + Ha_hfi_H6 + Ha_hfi_H8 + Ha_hfi_b1 + Ha_hfi_b2 + Ha_hfi_b3
	H_a_eigen = H_a.eigenenergies()

	# Now doing the same for the other TrpH molecule
	Hb_zeeman_nt = omega * (Sx * mt.sin(theta2) *  mt.cos(phi2) + Sy *  mt.sin(theta2) * mt. sin(phi2) + Sz * mt.cos(theta2))			# Hamiltonian without any tensor product
	Ha_zeeman = qt.tensor(Hb_zeeman_nt, qt.identity(3**7))			# Compounding different states 
	# Using the anisotropic hyperfine tensors the hamiltonian comes down to this
	Hb_hfi_N1  = qt.tensor(a_iso2[0] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_N1 * I2_x) + qt.tensor(Sy, T_N1 * I2_y) + qt.tensor(Sz, T_N1 * I2_z), qt.identity(3**6))

	Hb_hfi_H1  = qt.tensor( qt.tensor(3**1), a_iso2[1] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_H1 * I2_x) + qt.tensor(Sy, T_H1 * I2_y) + qt.tensor(Sz, T_H1 * I2_z), qt.identity(3**5))

	Hb_hfi_H2  = qt.tensor( qt.tensor(3**2), a_iso2[2] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_H2 * I2_x) + qt.tensor(Sy, T_H2 * I2_y) + qt.tensor(Sz, T_H2 * I2_z), qt.identity(3**4))

	Hb_hfi_H4  = qt.tensor( qt.tensor(3**3), a_iso2[3] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_H4 * I2_x) + qt.tensor(Sy, T_H4 * I2_y) + qt.tensor(Sz, T_H4 * I2_z), qt.identity(3**3))

	Hb_hfi_H5  = qt.tensor( qt.tensor(3**4), a_iso2[4] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_H5 * I2_x) + qt.tensor(Sy, T_H5 * I2_y) + qt.tensor(Sz, T_H5 * I2_z), qt.identity(3**2))

	Hb_hfi_H7  = qt.tensor( qt.tensor(3**5), a_iso2[5] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_H7 * I2_x) + qt.tensor(Sy, T_H7 * I2_y) + qt.tensor(Sz, T_H7 * I2_z), qt.identity(3**1))

	Hb_hfi_Hb1 = qt.tensor( qt.tensor(3**6), a_iso2[6] * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) 
					+  qt.tensor(Sx, T_Hb1 * I2_x) + qt.tensor(Sy, T_Hb1 * I2_y) + qt.tensor(Sz, T_Hb1 * I2_z), qt.identity(3**0))

	H_b = Hb_zeeman + Hb_hfi_N1 + Hb_hfi_H1 + Hb_hfi_H2 + Hb_hfi_H4 + Ha_hfi_H5 + Hb_hfi_H7 + Hb_hfi_Hb1
	H_b_eigen = H_b.eigenenergies()
	
	print H_b_eigen
	print H_a_eigen

Hamiltonian_eigenvalues(30, 30, 0, 0)




















