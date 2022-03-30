###
# Aidan Moretz
# MECH 408
# Ch 5 Prob 15
# Compares different loadings for an infinite solid
# Date Created: 3/27/22
# Last Modified: 3/27/22
###

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi

def SigmaYY(x,y,load):
    return -2*load*y**3/(pi*(x**2 + y**2)**2)

def SigmaXY(x,y,load):
    return -2*load*x*y**2/(pi*(x**2 + y**2)**2)

P = 1 #kN
a = 1 #m
Resolution = 100
PlotWidth = 5
x = np.linspace(start=-PlotWidth, stop=PlotWidth, num=Resolution)

Heights = (1,10,100)
Fig1, AxialPlots = plt.subplots(len(Heights), sharex=True)
Fig1.suptitle("Axial Stress")
Fig2, ShearPlots = plt.subplots(len(Heights),sharex=True)
Fig2.suptitle("Shear Stress")
for c,pA,pS in zip(Heights,AxialPlots, ShearPlots):
    # find stresses
    # one force load
    StressAxialOneForce = SigmaYY(x=x,y=c,load=P)
    StressShearOneForce = SigmaXY(x=x, y=c, load=P)
    # two force load
    StressAxialTwoForce = SigmaYY(x=x-a/2, y=c, load=P/2) + SigmaYY(x=x+a/2, y=c, load=P/2)
    StressShearTwoForce = SigmaXY(x=x-a/2, y=c, load=P/2) + SigmaXY(x=x+a/2, y=c, load=P/2)
    #add to plot
    pA.set(title=f"y={c*a} m")
    pA.plot(x,StressAxialOneForce,label="Concentrated Load")
    pA.plot(x,StressAxialTwoForce,label="Two Loads")
    pA.legend()
    pA.xlim=(x[0], x[-1])

    pS.set(title=f"y={c*a} m")
    pS.plot(x,StressShearOneForce,label="Concentrated Load")
    pS.plot(x,StressShearTwoForce,label="Two Loads")
    pS.legend()
    pS.xlim=(x[0], x[-1])

pA.set(xlabel="horizontal position [m]", ylabel="Stress (per depth) [kN/m]")
pS.set(xlabel="horizontal position [m]", ylabel="Stress (per depth) [kN/m]")
plt.show()