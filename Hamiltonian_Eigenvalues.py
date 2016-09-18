# For Details see the supplements of Hancock et.al. 2016 Quantum Needle of Avaian Compass
__author__ = 'Rakshit Jain'
import qutip as qt
import math as mt
import numpy as np
import scipy as sc
# Basic 3*3 Tensor arays for Hyperfine interactions
# For FAD Radical
def Hamiltonian_eigenvalues(theta,phi): # Here o refers to the angle
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
	# Aligning the axis of TrpH molecule's axis in order to have correct zeeman effects
	Sx_b = -.34554 * Sx + .6869 * Sy + -.6395 * Sz
	Sy_b = -.8269 * Sx +  .0995 * Sy +  .5535 * Sz
	Sz_b =  .4438 * Sx +  .7199 * Sy +  .5336 * Sz
	# Defining Constants here
	omega = 1.4 * 10**9				#Larmor Frequency in Hz for 50 microT
	g = 2  							# Electron Parameter
	# Isotropic Hyperfine Coupling Vector
	#For FAD Radical
	a_iso1 = np.array([.5233, .1887, -0.3872, 0.4399, 0.4070])
	T_N5n = np.array([1.2336,-0.6101,-0.6234])
	T_N5 = qt.Qobj(T_N5n) 
	T_N10n = np.array([0.4159, -0.2031, 0.2128])
	T_N10 = qt.Qobj(T_N10n)
	T_H6n = np.array([0.1896, -0.0464, 0-0.1432])
	T_H6 = qt.Qobj(T_H6n)
	T_H8n = np.zeros((3))
	T_H8 = qt.Qobj(T_H8n)
	T_Hbn = np.zeros((3))
	T_Hb = qt.Qobj(T_Hbn)
	# Similarly for TrpH Radical
	a_iso2 = np.array([.3215, -.5983, -.2780, -.488, -.0400, -.3636, 1.6046])
	T_N1n = np.array([.7596, -.3745, -.3851])
	T_N1 = qt.Qobj(T_N1n)
	T_H1n = np.array([.5914, -.1071,  -.4853])
	T_H1 = qt.Qobj(T_H1n)
	T_H2n = np.array([.2855, -.0919, -.1936])
	T_H2 = qt.Qobj(T_H2n)
	T_H4n = np.array([.3001, -.0480, -.2520])
	T_H4 = qt.Qobj(T_H4n)
	T_H5n = np.array([-.0632, -.0616, .1248])
	T_H5 = qt.Qobj(T_H5n)
	T_H7n = np.array([.2540, -.0594, -.1945])
	T_H7 =  qt.Qobj(T_H7n)
	T_Hb1n = np.array([0.1521, -.0456, -.1065])
	T_Hb1 = qt.Qobj(T_Hb1n)

	# Now since we have defined all the tensors required till now, Writing down the hamiltonian
	# H_a denotes the hamiltonian for FAD Radical
	# H_b denotes the Hamiltonian for TrpH radical

	Ha_zeeman_nt = omega * (Sx * mt.sin(theta) *  mt.cos(phi) + Sy *  mt.sin(theta) * mt. sin(phi) + Sz * mt.cos(theta))
	Ha_zeeman = qt.tensor(Ha_zeeman_nt, qt.identity(3**7))			# Compounding different states 
	# Using the anisotropic hyperfine tensors the hamiltonian comes down to this
	Ha_hfi_N5  = qt.tensor( qt.identity(3**0), float(a_iso1[0]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) + qt.tensor(Sx, T_N5n[0] * I2_x) + qt.tensor(Sy, T_N5n[1] * I2_y) + qt.tensor(Sz, T_N5n[2] * I2_z), qt.identity(3**6))

	Ha_hfi_N10 = qt.tensor( qt.identity(3**1),  float(a_iso1[1])*(qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_N10n[0] * I2_x) + qt.tensor(Sy, T_N10n[1] * I2_y) + qt.tensor(Sz, T_N10n[2] * I2_z), qt.identity(3**5))

	Ha_hfi_H6  = qt.tensor( qt.identity(3**2), float(a_iso1[2]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_H6n[0] * I2_x) + qt.tensor(Sy, T_H6n[1] * I2_y) + qt.tensor(Sz, T_H6n[2] * I2_z), qt.identity(3**4))

	Ha_hfi_H8  = qt.tensor( qt.identity(3**3), float(a_iso1[3]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_H8n[0] * I2_x) + qt.tensor(Sy, T_H8n[1] * I2_y) + qt.tensor(Sz, T_H8n[2] * I2_z), qt.identity(3**3))

	Ha_hfi_b1  = qt.tensor( qt.identity(3**4), float(a_iso1[4]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_Hbn[0] * I2_x) + qt.tensor(Sy, T_Hbn[1] * I2_y) + qt.tensor(Sz, T_Hbn[2] * I2_z), qt.identity(3**2))

	Ha_hfi_b2  = qt.tensor( qt.identity(3**5), float(a_iso1[4]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_Hbn[0] * I2_x) + qt.tensor(Sy, T_Hbn[1] * I2_y) + qt.tensor(Sz, T_Hbn[2] * I2_z), qt.identity(3**1))

	Ha_hfi_b3  = qt.tensor( qt.identity(3**6), float(a_iso1[4]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_Hbn[0] * I2_x) + qt.tensor(Sy, T_Hbn[1] * I2_y) + qt.tensor(Sz, T_Hbn[2] * I2_z), qt.identity(3**0))
	# Everything is pretty clear upto this point
	H_a1 = Ha_hfi_N5.data + Ha_hfi_N10.data + Ha_hfi_H6.data + Ha_hfi_H8.data + Ha_hfi_b1.data + Ha_hfi_b2.data + Ha_hfi_b3.data + Ha_zeeman.data
	H_a = qt.Qobj(H_a1)

	H_a_eigen = H_a.eigenenergies()
	

	# Now doing the same for the other TrpH molecule
	Hb_zeeman_nt = omega * (Sx_b * mt.sin(theta) *  mt.cos(phi) + Sy_b *  mt.sin(theta) * mt. sin(phi) + Sz_b * mt.cos(theta))			# Hamiltonian without any tensor product
	Hb_zeeman = qt.tensor(Hb_zeeman_nt, qt.identity(3**7))			# Compounding different states 
	# Using the anisotropic hyperfine tensors the hamiltonian comes down to this
	Hb_hfi_N1  = qt.tensor(float(a_iso2[0]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_N1n[0] * I2_x) + qt.tensor(Sy, T_N1n[1] * I2_y) + qt.tensor(Sz, T_N1n[2] * I2_z), qt.identity(3**6))

	Hb_hfi_H1  = qt.tensor( qt.identity(3**1), float(a_iso2[1]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_H1n[0] * I2_x) + qt.tensor(Sy, T_H1n[1] * I2_y) + qt.tensor(Sz, T_H1n[2] * I2_z), qt.identity(3**5))

	Hb_hfi_H2  = qt.tensor( qt.identity(3**2), float(a_iso2[2]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_H2n[0] * I2_x) + qt.tensor(Sy, T_H2n[1] * I2_y) + qt.tensor(Sz, T_H2n[2] * I2_z), qt.identity(3**4))

	Hb_hfi_H4  = qt.tensor( qt.identity(3**3), float(a_iso2[3]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_H4n[0] * I2_x) + qt.tensor(Sy, T_H4n[1] * I2_y) + qt.tensor(Sz, T_H4n[2] * I2_z), qt.identity(3**3))

	Hb_hfi_H5  = qt.tensor( qt.identity(3**4), float(a_iso2[4]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_H5n[0] * I2_x) + qt.tensor(Sy, T_H5n[1] * I2_y) + qt.tensor(Sz, T_H5n[2] * I2_z), qt.identity(3**2))

	Hb_hfi_H7  = qt.tensor( qt.identity(3**5), float(a_iso2[5]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_H7n[0] * I2_x) + qt.tensor(Sy, T_H7n[1] * I2_y) + qt.tensor(Sz, T_H7n[2] * I2_z), qt.identity(3**1))

	Hb_hfi_Hb1 = qt.tensor( qt.identity(3**6), float(a_iso2[6]) * (qt.tensor(Sx, I2_x) + qt.tensor(Sy, I2_y) + qt.tensor(Sz, I2_z)) +  qt.tensor(Sx, T_Hb1n[0] * I2_x) + qt.tensor(Sy, T_Hb1n[1] * I2_y) + qt.tensor(Sz, T_Hb1n[2] * I2_z), qt.identity(3**0))
				# Hamiltonian without any tensor product

	H_b1 = Hb_zeeman.data + Hb_hfi_N1.data + Hb_hfi_H1.data + Hb_hfi_H2.data + Hb_hfi_H4.data + Hb_hfi_H5.data + Hb_hfi_H7.data + Hb_hfi_Hb1.data
	H_b = qt.Qobj(H_b1)
	H_b_eigen = H_b.eigenenergies()

	print (H_b.isherm)
	print (H_a.isherm)


Hamiltonian_eigenvalues(30, 0)





















