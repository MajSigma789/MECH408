#####
# Aidan Moretz
# MECH 408
# Geometric Non-linear Beam
# Problem 3
# Created: 3/3/2022
# Last Modified: 3/4/2022
#####

import numpy as np
from scipy.special import ellipk, ellipkinc, ellipe, ellipeinc

def findLoad(angle):
    """Calculates the load needed for an free end angle."""
    k = .5*(1 + np.sin(angle)) # -
    phi = np.arcsin((np.sqrt(2*k))**(-1)) #rad
    load = EI*((ellipk(k) - ellipkinc(phi,k))/L)**2
    return load

def nonlinearAngle(load,error=1E-6, maxIterations=50):
    """Uses the bisection method to find the free end angle for a given load."""
    angle1 = 0 + error
    angle2 = np.pi/2 - error
    count = 0
    residual1 = findLoad(angle1) - load
    residual2 = findLoad(angle2) - load
    #bisection to find 
    while np.abs(.5*(residual1 + residual2)) > error and count < maxIterations:
        angleNew = .5*(angle1 + angle2)
        residualNew = findLoad(angleNew) - load
        if residualNew == 0:
            break
        elif residualNew < 0:
            angle1 = angleNew
            residual1 = residualNew
        else:
            angle2 = angleNew
            residual2 = residualNew
        count += 1
    #error if takes too long
    if count == maxIterations:
        raise RuntimeError("Max Iterations exceeded!")
    return angleNew

def nonlinearDisplacement(load, angle):
    k = .5*(1 + np.sin(angle)) # -
    phi = np.arcsin((np.sqrt(2*k))**(-1)) #rad
    alpha = np.sqrt(load/E/I) #1/m
    return L - 2/alpha*(ellipe(k) - ellipeinc(phi, k)) #m

def linearDisplacement(load):
    return load*L**3/3/EI

#LOWPRIO: Take user defined inputs
s = 10E-3 #m
I = (s**4)/12 #m^4
E = 200E9 #Pa = N/m^2
EI = E*I #N*m^2
L = 1 #m
loads = (.5,5,50,100,250,500) #N

for P in loads:
    # find displacements
    angle = nonlinearAngle(load=P)
    displaceNL = nonlinearDisplacement(load=P,angle=angle)
    displaceL = linearDisplacement(load=P)
    error = (displaceL - displaceNL)/displaceNL*100
    # print displacements
    print(f"Force: {P}N")
    print("Nonlinear Displacement: {:5.6f}m".format(displaceNL))
    print("Linear Displacement:    {:5.6f}m".format(displaceL))
    print("Linear Error:           {:.2f}%".format(error))
