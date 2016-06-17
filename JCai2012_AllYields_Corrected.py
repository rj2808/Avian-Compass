#
# Dynamics of radical pair subjected to Geomagnetic field
#
from qutip import *
from scipy import *
from pylab import *
import csv
import math

def run():

# problem parameters:
    pi=math.pi
    muB = 5.788*10**-5  # Bohr Magneton in eV/Tesla
    g=2 # g-factor for electron
    #gama=0.5*muB*g # in eV/Tesla
    gama=muB*g # in eV/Tesla
    gama=gama*0.1519756*10**16 # in sec-1/Tesla (Converted using natural units)
    B0=47*10**-6 # Geomagnetic field in Frankfurt (Tesla)
    #B0=46*10**-6 # Geomagnetic field in Frankfurt (Tesla)
    Brf=150*10**-9 # Disturbing RF field of strength 150 nanoTesla
    w = 2*pi*1.316*10**6      # Frequency of externally applied RF field
    sqrt2=math.sqrt(2)
##########################
# Initial State
    up=basis(2,0)
    down=basis(2,1)	
    singlet=(tensor(up,down)-tensor(down,up))/sqrt2 # Initial State of the radical pair
    trip0=(tensor(up,down)+tensor(down,up))/sqrt2 # Triplet state with zero spin
    tripu= tensor(up,up) # Tirplet state with spin +1
    tripd=tensor(down,down) # Triplet state with spin -1
    #s1=[[0 for x in xrange(1)] for x in xrange(10)]
    s1=[[0 for x in xrange(1)] for x in xrange(12)]
    #s2=[[0 for x in xrange(1)] for x in xrange(10)]
    s2=[[0 for x in xrange(1)] for x in xrange(12)]
    s1[1][0]=1/sqrt2;
    s1[2][0]=-1/sqrt2;
    s2[5][0]=1/sqrt2;
    s2[6][0]=-1/sqrt2;
    
    e_dm=singlet*singlet.dag()
    nuc_dm=0.5*qeye(2)
    rho=tensor(nuc_dm,e_dm)
    rho0=Qobj(pad(rho.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    
# Defining Hamiltonian operators

    S1x=tensor(qeye(2),sigmax(),qeye(2))
    S2x=tensor(qeye(2),qeye(2),sigmax())
    #S1x_10=Qobj(pad(S1x.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0))) # pad() function increases the size of array of S1x from
    S1x_12=Qobj(pad(S1x.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0))) # pad() function increases the size of array of S1x from  
    #8x8 to 10x10 and put 0s at the newly created places
    #S2x_10=Qobj(pad(S2x.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    S2x_12=Qobj(pad(S2x.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    #Sx=S1x_10+S2x_10
    Sx=S1x_12+S2x_12
    S1y=tensor(qeye(2),sigmay(),qeye(2))
    S2y=tensor(qeye(2),qeye(2),sigmay())
    #S1y_10=Qobj(pad(S1y.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    S1y_12=Qobj(pad(S1y.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    #S2y_10=Qobj(pad(S2y.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    S2y_12=Qobj(pad(S2y.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    #Sy=S1y_10+S2y_10
    Sy=S1y_12+S2y_12
    S1z=tensor(qeye(2),sigmaz(),qeye(2))
    S2z=tensor(qeye(2),qeye(2),sigmaz())
    #S1z_10=Qobj(pad(S1z.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    S1z_12=Qobj(pad(S1z.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    #S2z_10=Qobj(pad(S2z.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    S2z_12=Qobj(pad(S2z.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    #Sz=S1z_10+S2z_10
    Sz=S1z_12+S2z_12
    #ISx8=tensor(sigmax(),sigmax(),qeye(2))
    ISx8=tensor(sigmax(),qeye(2),sigmax())
    #ISy8=tensor(sigmay(),sigmay(),qeye(2))
    ISy8=tensor(sigmay(),qeye(2),sigmay())
    #ISz8=tensor(sigmaz(),sigmaz(),qeye(2))
    ISz8=tensor(sigmaz(),qeye(2),sigmaz())
    #ISx=Qobj(pad(ISx8.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    ISx=Qobj(pad(ISx8.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    #ISy=Qobj(pad(ISy8.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    ISy=Qobj(pad(ISy8.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))
    #ISz=Qobj(pad(ISz8.full(), ((0, 2), (0, 2)), 'constant', constant_values=(0)))
    ISz=Qobj(pad(ISz8.full(), ((0, 4), (0, 4)), 'constant', constant_values=(0)))

    # Decay rates [These three rates are simulated]
    k=5*10**5
    #k=10**3
    #print('k = ', k)    
    #lok=len(k)
    print('k = ', k)

    # Defining Projection operators    
   
    #P=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s

    #P1=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P1=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    P1[8][1]=1/sqrt2
    P1[8][2]=-1/sqrt2
    #P1=sqrt(k[ki])*Qobj(P1)
    P1=sqrt(k)*Qobj(P1)
    
    #P2=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P2=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    P2[8][5]=1/sqrt2
    P2[8][6]=-1/sqrt2
    #P2=sqrt(k[ki])*Qobj(P2)
    P2=sqrt(k)*Qobj(P2)
   
    #P3=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P3=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    P3[9][1]=1/sqrt2
    P3[9][2]=1/sqrt2
    #P3=sqrt(k[ki])*Qobj(P3)
    P3=sqrt(k)*Qobj(P3)

    #P4=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P4=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    P4[9][5]=1/sqrt2
    P4[9][6]=1/sqrt2
    #P4=sqrt(k[ki])*Qobj(P4)
    P4=sqrt(k)*Qobj(P4)

    #P5=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P5=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    #P5[9][0]=1
    P5[10][0]=1
    #P5=sqrt(k[ki])*Qobj(P5)
    P5=sqrt(k)*Qobj(P5)

    #P6=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P6=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    #P6[9][3]=1
    P6[10][4]=1
    #P6=sqrt(k[ki])*Qobj(P6)
    P6=sqrt(k)*Qobj(P6)

    #P7=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P7=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    #P7[9][4]=1
    P7[11][3]=1
    #P7=sqrt(k[ki])*Qobj(P7)
    P7=sqrt(k)*Qobj(P7)
	
    #P8=[[0 for x in xrange(10)] for x in xrange(10)] # 10x10 array containing only 0s
    P8=[[0 for x in xrange(12)] for x in xrange(12)] # 10x10 array containing only 0s
    #P8[9][7]=1
    P8[11][7]=1
    #P8=sqrt(k[ki])*Qobj(P8)
    P8=sqrt(k)*Qobj(P8)

    # Defining Collapse Operators
    c_op_list = [P1,P2,P3,P4,P5,P6,P7,P8]
    #c_ops=[]
    #print('c_ops', c_ops)

    # Time steps in Integration. It would take 6/k time to decay the population by almost 99%
    tlist=[]
    tlist = linspace(0.0, float(6)/k, 2000) # Defining time instants for simulation
    #tlist = linspace(0.0, float(6)/k, 10000) # Defining time instants for simulation
    noe=len(tlist)
    dt=tlist[1]-tlist[0]
    print("dt = ", dt)

    # Geomagnetic field angle (in radian)
    theta_list=linspace(0.0, pi/2, 40)
    #theta_list=linspace(pi/2, pi, 1)
    noa=len(theta_list)

# Hyperfine Coupling Constants (in eV)
    ax=0
    ay=0
    #a=logspace(log10(10**-2),log10(10**2),200)
    #a=logspace(log10(10**1),log10(10**2),50)
    #a=logspace(log10(float(1)/3),log10(3),1)
    a=logspace(log10(float(100.00)),log10(200),1)
    az=[x*gama*B0 for x in a]
    loa=len(az)
    print('No of hyperfines = ', loa)

    # Magnetic field made zero for calculating coherence
    #B0=0
    print('B0 = ', B0)
    print('a/B = ', a[0])

    for ki in range(0,loa):  # Loop over various hyperfine constant values
	singlet_yield=[] # Defining Singlet Yield variable
	singlet_yield_coh=[] # Defining Singlet Yield variable
	trip0_yield=[] # Defining Triplet Yield variable
	tripu_yield=[] # Defining Triplet Yield variable
	tripd_yield=[] # Defining Triplet Yield variable
	Cohr=0

    	for th in range(0,noa): # Loop over various angle values
	    Bx=B0*sin(theta_list[th])# Becomes most effective at 90 degrees
    	    By=0
  	    Bz=B0*cos(theta_list[th]) # Becomes most effective at 0 degrees
            Brfx = Brf*sin(3*pi/2+theta_list[th])
            Brfy = 0
            Brfz = Brf*cos(3*pi/2+theta_list[th])

	    H0 = gama*Bx*Sx/2 + gama*By*Sy/2 + gama*Bz*Sz/2 + ax*ISx/4 + ay*ISy/4 + az[ki]*ISz/4
	    #H0_coh = ax*ISx/4 + ay*ISy/4 + az[ki]*ISz/4 # Hyperfine Interaction Hmailtonian Only
	    #H0 = Bx*Sx/2 + By*Sy/2 + Bz*Sz/2 + ax*ISx/4 + ay*ISy/4 + az*ISz/4
            #H1 = gama*Brfx*Sx/2 + gama*Brfy*Sy/2 + gama*Brfz*Sz/2
	    #args={'w':w}
	    #H=[H0,[H1,'cos(w*t)']] # String method for writing time dependent Hamiltonian
    	    # Running the Master Equation solver
	    opts = Odeoptions(nsteps=5000)

	    result = mesolve(H0, rho0, tlist, c_op_list, [],options=opts)
	    #result_coh = mesolve(H0_coh, rho0_offd, tlist, c_op_list, [],options=opts)
	    #print('theta = ', theta_list[th],"Done")
	    #result = mesolve(H, rho0, tlist, c_op_list, [], args=args, options=None) # Hamiltonian with RF field also.
	    #rh_ss=steadystate(H, c_op_list,maxiter=500, tol=1e-06, method='mesolve')
	
    	    # Output states of the differential equation solver
	    rh=[]
	    rh_coh=[]
	    rh_static=[]
	    rh=result.states
	    singlet_yield.append(abs(rh[noe-1][8][0][8]))
	    trip0_yield.append(abs(rh[noe-1][9][0][9]))
	    tripu_yield.append(abs(rh[noe-1][10][0][10]))
	    tripd_yield.append(abs(rh[noe-1][11][0][11]))
	    '''
	    # Calculating L1-normed coherence from density matrix
	    for nn in range(0,noe):
		rh_array=rh[nn].full() # Converting quantum object into an array
		for ii in range(0,8):
		    for jj in range(0,8):
			if ii==jj:
			   Cohr=Cohr
			else:
			   Cohr=Cohr+abs(rh_array[ii][jj])
	    #rh_coh=result_coh.states
	    #print rh
	    #time=result.times
	    #noc=result.num_collapse
	    #print('noc',noc)
	    singlet_yield.append(abs(rh[noe-1][8][0][8]))
	    #singlet_yield_coh.append(abs(rh_coh[noe-1][8][0][8]))
 	    triplet_yield.append(abs(rh[noe-1][9][0][9]))		
	print('a/B = ', a[ki],"Done!")	
	Sens=max(singlet_yield)-min(singlet_yield)
        print('Sensitivity = ', Sens)
	sensitivity.append(Sens)
	#Cohr=max(singlet_yield_coh)	
	coherence.append(Cohr)
	'''
    #Sens=max(syield_theta11)-min(syield_theta11);

   # Saving data to a file
    ar_data=array([theta_list,singlet_yield,trip0_yield,tripu_yield,tripd_yield])
    fl = open('data/JCai2012_allyields_at_100.00B.csv', 'w')
    writer = csv.writer(fl)
    writer.writerow(['Theta', 'Singlet Yield','T0 Yield','TU Yield','TD Yield']) #if needed
    for values in ar_data:
        writer.writerow(values)
    fl.close()     

    #Plotting Sensitivity/Singlet Yield
    figure()
    plot(theta_list, singlet_yield, 'r',
 	 theta_list, trip0_yield, 'b',
	 theta_list, tripu_yield, 'g',
  	 theta_list, tripd_yield, 'k')

    ylim([0.00,1.00])
    xlabel("Angle "r'$\theta$',fontsize=18)
    ylabel('Yields', fontsize=18 )
    #legend(("No RF, k = $10^6 s^{-1}$", "With RF, k = $10^6 s^{-1}$","No RF, k = $10^5 s^{-1}$","With RF, k = $10^5 s^{-1}$"),loc='upper left')
    legend(("S", "T0","T+","T-"), loc='upper center',bbox_to_anchor=(0.50,1.10), ncol=4, fancybox=True, shadow=True)
    minorticks_on()
    xticks( theta_list, ('$0$', ' ',' ', ' ', ' ', r'$\frac{\pi}{8}$',' ', ' ', ' ', ' ', r'$\frac{\pi}{4}$',' ', ' ', ' ', ' ', r'$\frac{3\pi}{8}$', ' ', ' ', ' ',r'$\frac{\pi}{2}$'), fontsize=18 )
    savefig('JCai2012_allyields_at_100.00B.eps')
    savefig('JCai2012_allyields_at_100.00B.jpg')
    show()
if __name__ == '__main__':
    run()
