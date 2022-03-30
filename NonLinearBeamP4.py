###
# Aidan Moretz
# MECH 408
# Geometric Nonlinear Beam 
# Problem 4
# Date Created: 3/6/2022
# Last Modified: 3/6/2022
###

import numpy as np
from numpy import pi
from scipy.special import ellipk

def findLoad(angle):
    k = np.sin(.5*angle)**2
    return EI*(ellipk(k)/L)**2

#LOWPRIO: input table for manual parameters
# parameters
s = 10E-3 #m
I = (s**4)/12 #m^4
E = 200E9 #Pa = N/m^2
EI = E*I #N*m^2
L = 1 #m
thetaL = -pi/3

criticalLoad = pi**2*EI/(4*L**2) #N
currentLoad = findLoad(thetaL)

changePercent = 100*np.abs(criticalLoad - currentLoad)/criticalLoad
print("To bend to {:.1f} deg, it takes {:.2f}% more than the buckling load.".format(np.degrees(thetaL), changePercent))