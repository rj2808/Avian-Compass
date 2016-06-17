%matplotlib
import matplotlib.pyplot as plt
import numpy as num
from qutip import *

def qubit_integrate(epsilon, delta, g1, g2, solver):
	H = epsilon / 2.0 * sigmaz() + delta / 2.0 * sigmax()
	#collapse operators
 	c_ops = []
 	if g1 > 0.0 :
 		c_ops.append(num.sqrt(g1) * sigmam())
 	if g2 > 0.0 :
 		c_ops.append(num.sqrt(g2) * sigmax())
 	e_ops = [sigmax(), sigmay(), sigmaz()]

 	if solver == "me" :
 		output = mesolve(H, psi0, tlist, c_ops, e_ops)
 	elif solver == "es" :
 		output = essolve(H, psi0, tlist, c_ops, e_ops)
 	elif solver == "mc" :
 		ntraj = 250
 		output = mcsolve(H, psi0, tlist, ntraj, c_ops, [sigmax(), sigmay(), sigmaz()])
 	else:
 		 raise ValueError("Unknown Solver")
 	return output.expect[1], output.expect[2], output.expect[3]


