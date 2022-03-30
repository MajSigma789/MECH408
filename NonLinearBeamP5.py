###
# Aidan Moretz
# MECH 408
# Geometric Nonlinear Beam
# Problem 5
# Date Created: 3/6/2022
# Last Modified: 3/6/2022
###

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipeinc, ellipkinc

def curvePoints(angle, resolution=30):
    # form loading constants
    k = np.sin(angle/2)**2 #-
    beta = ellipk(k)/L #m^-1
    # find curve points
    phi = np.linspace(0,pi/2, resolution)
    x = (2*ellipeinc(phi,k) - ellipkinc(phi, k))/beta #m
    y = -2*np.sqrt(k)/beta*(1-np.cos(phi)) #m
    return x,y

#LOWPRIO: input table for manual parameters
# parameters
s = 10E-3 #m
I = (s**4)/12 #m^4
E = 200E9 #Pa = N/m^2
EI = E*I #N*m^2
L = 1 #m

angles = [-20, -40, -60, -80, -100, -120, -140, -160, -176]
thetaL = np.radians(np.array(angles)) #rad
# plot curves
plt.figure()
plt.title("Horizontal Beam Deflection for Various Free End Angles")
for angle, labelAngle in zip(thetaL,angles):
        x,y = curvePoints(angle=angle)
        plt.plot(x,y,'-',label="ThetaL: {:d}".format(labelAngle))
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.legend()
plt.show()
