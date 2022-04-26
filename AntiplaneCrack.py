###
# Aidan Moretz
# MECH 408
# Antiplane Crack
# Date Created: 4/11/2022
# Last Modified: 4/11/2022
###

import numpy as np
import matplotlib.pyplot as plt
from MechCustom import InputTable


parameters = {"mu":("GPa", 1), "alpha":("deg", 1), "R":("mm", 1), "load":("MPa", 1), "resolution":("-", 1), "terms":("-", 1)}
inTable = InputTable(rows=parameters, title="System Parmeters", button="Run")
[mu, alpha, R, load, resolution, terms] = [values[0] for values in inTable.getValues(entries=parameters.keys())]
inTable.destroy()
mu *= 1000
resolution = int(resolution)
terms = int(terms)
alpha = np.radians(alpha)

###mu = 80.E3 #MPa
###alpha = np.deg2rad(135.) #rad
###R = 250. #mm
###load = 30. #MPa
###posError = .001 #mm
###resolution = 200
###terms = 10

theta, radius = np.meshgrid(np.linspace(-alpha, alpha, resolution), np.linspace(1E-5, R, int(resolution/2)))
displacement = np.zeros_like(radius)
stressRZ = np.zeros_like(radius)
stressTZ = np.zeros_like(radius)


# find displacement and stresses until desired cutoff
for n in range(terms):
    lam = (2*n + 1)*np.pi/2/alpha
    a = 2*load/(mu*(2*n + 1)*np.pi*R**(lam - 1))
    sinTerm = np.sin(lam*theta)
    displaceCoeff = a*radius**lam

    displacement += displaceCoeff*sinTerm
    stressCoeff = displaceCoeff*mu*lam
    stressRZ += stressCoeff*sinTerm/radius
    stressTZ += stressCoeff*np.cos(lam*theta)

# normalize values
radiusNorm = radius/R
displacementNorm = displacement/R
stressRZNorm = stressRZ/load
stressTZNorm = stressTZ/load

maxDisplace = max((max(np.abs(rValue)) for rValue in displacement))

fig1, ax1 = plt.subplots(subplot_kw={'projection':'polar'})
#INPUT MUST BE IN RADIANS!!!
DisplacePlot = ax1.contourf(theta, radiusNorm, displacementNorm)
ax1.set_title("Normalized Displacement")
fig1.colorbar(DisplacePlot)

fig2, ax2 = plt.subplots(subplot_kw={'projection':'polar'})
#INPUT MUST BE IN RADIANS!!!
stressRZPlot = ax2.contourf(theta, radiusNorm, stressRZNorm)
ax2.set_title("Normalized R-Z Shear Stress")
fig2.colorbar(stressRZPlot)

fig3, ax3 = plt.subplots(subplot_kw={'projection':'polar'})
#INPUT MUST BE IN RADIANS!!!
stressTZPlot = ax3.contourf(theta, radiusNorm, stressTZNorm)
ax3.set_title("Normalized Theta-Z Shear Stress")
fig3.colorbar(stressTZPlot)

plt.show()

