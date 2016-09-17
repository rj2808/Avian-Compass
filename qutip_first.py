from qutip import *
from scipy import *
from pylab import *
import math
import csv


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


    #Initial State
    up=basis(2,0)
    down=basis(2,1)	
    singlet=(tensor(up,down)-tensor(down,up))/sqrt2 # Initial State of the radical pair
    trip0=(tensor(up,down)+tensor(down,up))/sqrt2 # Triplet state with zero spin
    tripu= tensor(up,up) # Tirplet state with spin +1
    tripd=tensor(down,down) # Triplet state with spin -1
    


