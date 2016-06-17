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





epsilon = 0.0 * 2 * num.pi
delta	 = 1.0 * 2 * num.pi
g2 = .15
g1 = 0.0

# intial state
psi0 = basis(2,0) # In State Up
tlist = num.linspace(0,5,200)

#analytics
sx_analytic = num.zeros(shape(tlist))
sy_analytic = -num.sin(2*np.pi*tlist) * num.exp(-tlist * g2)
sz_analytic = num.cos(2*np.pi*tlist) * num.exp(-tlist * g2)


fig, ax = plt.subplots(figsize=(12,6))
ax.plot(tlist, num.real(sx1), 'r')
ax.plot(tlist, num.real(sy1), 'b')
ax.plot(tlist, num.real(sz1), 'g')
ax.plot(tlist, sx_analytic, 'r*')
ax.plot(tlist, sy_analytic, 'g*')
ax.plot(tlist, sz_analytic, 'g*')
ax.legend(("sx", "sy", "sz"))
ax.set_xlabel('Time')
ax.set_ylabel('expectation value');





