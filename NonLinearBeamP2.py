#####
# Aidan Moretz
# MECH 408
# Geometric Non-linear Beam
# Problem 2
# Created: 2/27/2022
# Last Modified: 3/3/2022
#####

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipkinc, ellipe, ellipeinc

E = 200.0 #GPa = GN/m^2
I = 2.3974E-5 #m^4
L = 1.0 #m
increments = 20
# find loads for angle range
thetaL = np.array([-np.pi/4*k/increments for k in range(0,increments+1)]) #rad
k = .5*(1 - np.sin(thetaL)) # -
phi = np.arcsin((np.sqrt(2*k))**(-1)) #rad
load = E*I*((ellipk(k) - ellipkinc(phi,k))/L)**2 #GN
# find linear and non-linear displacements
alpha = np.sqrt(load/E/I) #1/m
yNLDisplace = L - 2/alpha*(ellipe(k) - ellipeinc(phi, k)) #m
yLDisplace = load/E/I*L**3 #m
load *= 1E3 #MN
# plot result
plt.figure()
plt.title("Load vs Displacement")
plt.plot(load, yNLDisplace, 'r.', label="Displacement (Nonlinear)")
plt.plot(load, yLDisplace, 'bs', label="Displacement (Small Deflection")
plt.xlabel("Load [MN = N*1E6]")
plt.ylabel("Y-Displacement [m]")
plt.legend()
plt.show()