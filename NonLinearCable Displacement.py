##
# Aidan Moretz
# MECH 408
# Nonlinear Cable Handout
# Created: 3/27/2022
# Modified: 4/7/2022
# Uses a corrected incremental approach to solve geometric nonlinearities
# in cable deformation with a diplacement controlled input
##

import numpy as np
import matplotlib.pyplot as plt
from MechCustom import InputTable, OutputTable

def tangentStiffness(position):
    """Forms the tangent stiffness matrix of the two cable system"""
    L1Offset = (L[0] + position[0])**(-1)
    L2SinOffset = L[1]*np.sin(alpha[1]) + position[1]
    L2CosOffset = L[1]*np.cos(alpha[1]) - position[0]
    w1 = np.sqrt(position[1]**2/(L1Offset)**(-2))
    w2 = (L2SinOffset)/(L2CosOffset)
    w1Term = (1 + w1**2)**(-1)
    w2Term = (1 + w2**2)**(-1)
    deltas = deltaCalculation(position=position)
    # linear material tangent stiffness matrix
    Kt = np.array([
        [k[0]*(w1Term + deltas[0]*w1*w1Term**(1.5)*position[1]*L1Offset**2) + k[1]*(w2Term + deltas[1]*w2*w2Term**1.5*(L2SinOffset)*(L2CosOffset)**(-2)), # dPx/du
        k[0]*(w1*w1Term - deltas[0]*w1*w1Term**(1.5)*L1Offset) - k[1]*(w2*w2Term - deltas[1]*w2*w2Term**1.5*(L[0]*np.cos(alpha[1]) - position[0])**(-1))], # dPx/dv
        [k[0]*(w1*w1Term + deltas[0]*w1Term**(1.5)*position[1]*L1Offset**2) - k[1]*(w2*w2Term - deltas[1]*w2Term**1.5*(L2SinOffset)*(L2CosOffset)**(-2)), # dPy/du
        k[0]*(w1**2*w1Term + deltas[0]*w1Term**(1.5)*L1Offset) + k[1]*(w2**2*w2Term + deltas[1]*w2Term**1.5*(L2CosOffset)**(-1))] # dPy/dv
        ])
    Kt = np.array([Kt[1][1] - Kt[1][0]*Kt[0][1]/Kt[0][0], -Kt[0][1]/Kt[0][0]])
    return Kt

def deltaCalculation(position):
    """Calculates the change in cable lengths for the current displacement"""
    delta1 = np.sqrt((L[0] + position[0])**2 + position[1]**2) - L[0]
    delta2 = np.sqrt((L[1]*np.cos(alpha[1]) - position[0])**2 + (L[1]*np.sin(alpha[1]) + position[1])**2) - L[1]
    return np.array((delta1, delta2))

#FIXME: something fucky comes this way in the NR iterations
def newtonRaphsonCorrection(position, load, stiffness=np.array([]), error=1E-6):
    """Uses the Newton-Raphson root finding method to enforce equillibrium"""
    # form stiffness if needed
    if np.array_equiv(stiffness,np.array([])):
        stiffness = tangentStiffness(position)
    
    count = 0
    global maxIterations
    while count < (maxIterations + 1)*5:
        # find residual forces
        updateGeometry(position)
        delta = deltaCalculation(position)
        reaction = k*delta
        rotation = np.array([[np.cos(beta[0]), -np.cos(beta[1])],[np.sin(beta[0]), np.sin(beta[1])]])
        residualForce = load - rotation@reaction
        # end if close to equillibrium
        if np.linalg.norm(residualForce) < error:
            maxIterations = max(count, maxIterations)
            return position
        # move towards equillibrium
        deltaCorrection = np.linalg.solve(stiffness, residualForce)
        position += deltaCorrection
        # prep for next iteration
        stiffness = tangentStiffness(position)
        count += 1

    raise RuntimeError("Failed to converge correction in time.")

def updateGeometry(position):
    """Repositions the cable angles for the current displacement."""
    global beta
    beta = [np.arcsin(position[1]/np.sqrt((L[0] + position[0])**2 + position[1]**2)), 
    np.arcsin((L[1]*np.sin(alpha[1]) + position[1])/np.sqrt((L[1]*np.cos(alpha[1]) - position[0])**2 + (L[1]*np.sin(alpha[1]) + position[1])**2))]

def incrementalLinearization(position, dPosition):
    """Advances the linearization by one step to find the change in position
    Parameters: dPosition [mm]
    Returns: deltaP [kN], deltaU [mm], Kt [kN/mm]"""
    Kt = tangentStiffness(position) #kN/mm
    deltaP, deltaU = Kt*dPosition
    return deltaP, deltaU, Kt

def getStresses(delta):
    """Calculates stress at a given deformation using the current element stiffness
    Parameters: delta [mm]
    Returns: stress [GPa = kN/mm^2]"""
    return k*delta/A #kN/mm^2 = GPa


# main method
data = {"L":("mm",2), "E":("GPa",2), "A":("mm^2",2), "alpha":("rad",2), "displacement":("mm",1), "increments":("-",1)}
# create input window
inTable = InputTable(rows=data, title="System Parameters", button="Run")

# get input parameters
L = np.array([float(entry.get()) for entry in inTable.values["L"]]) #mm
E = np.array([float(entry.get()) for entry in inTable.values["E"]]) #GPa = kN/mm^2
A = np.array([float(entry.get()) for entry in inTable.values["A"]]) #mm^2
alpha = np.array([float(entry.get()) for entry in inTable.values["alpha"]]) #rad
totalDisplacement = float(inTable.values["displacement"][0].get()) #mm
increments = int(inTable.values["increments"][0].get())
inTable.destroy()

# derived parameters
k = E*A/L #kN/mm
dPosition = totalDisplacement/increments #kN

# variables
beta = np.copy(alpha) #rad
positionIncrements = [dPosition]*increments #kN

# uncorrected solution
positionHistoryUncorrected = ([0.,], [0.,])
forceHistoryUncorrected = [0.,]
force = 0
for dV in positionIncrements:
    # step forward in loading amount
    position = np.array([history[-1] for history in positionHistoryUncorrected])
    deltaP, deltaU = incrementalLinearization(position=position, dPosition=dV)[0:2]
    position += np.array([deltaU,dV])
    force += deltaP
    # add component positions to history
    for history, displace in zip(positionHistoryUncorrected,position): history.append(displace)
    forceHistoryUncorrected.append(force)
finalDelta = deltaCalculation(position=position)
stressUC = 1000*getStresses(delta=finalDelta) #MPa

### corrected solution
##positionHistoryCorrected = ([0.,], [0.,])
##maxIterations = 0
##totalLoad = 0
##for dP in loadIncrements:
##    # step forward in loading amount
##    position = np.array([history[-1] for history in positionHistoryCorrected])
##    deltaPos, Kt = incrementalLinearization(position=position,load=dP)
##    totalLoad += dP
##    position = newtonRaphsonCorrection(position=position+deltaPos, load=totalLoad, stiffness=Kt)
##    for history, displace in zip(positionHistoryCorrected,position): history.append(displace)
##finalDelta = deltaCalculation(position=position)
##stressC = 1000*getStresses(delta=finalDelta) #MPa

# report uncorrected info
uHistoryUC, vHistoryUC = np.array(positionHistoryUncorrected[0]), np.array(positionHistoryUncorrected[1])
deltaHistoryUC = np.sqrt(uHistoryUC**2 + vHistoryUC**2)
rows = {"Final X":((f"{uHistoryUC[-1]:.3f}"), "mm"), "Final Y":((f"{vHistoryUC[-1]:.3f}"), "mm"), 
"Final Py":((f"{forceHistoryUncorrected[-1]:.3f}"), "kN"), 
"Stress 1":((f"{stressUC[0]:.2f}"), "MPa"), "Stress 2":((f"{stressUC[1]:.2f}"), "MPa")}
OutputTable(rows=rows, columnTitles=("Uncorrected"), title="Simulation Results")
###uHistoryC, vHistoryC = positionHistoryCorrected[0], positionHistoryCorrected[1]
###data = {"Final X":((f"{uHistoryUC[-1]:.3f}",f"{uHistoryC[-1]:.3f}"), "mm"), "Final Y":((f"{vHistoryUC[-1]:.3f}",f"{vHistoryC[-1]:.3f}"), "mm"),
###"Final Py":((f"{forceHistoryUncorrected[-1]:.3f}", f"{forceHistoryCorrected[-1]:.3f}"),"kN"), 
###"Stress 1":((f"{stressUC[0]:.2f}", f"{stressC[0]:.2f}"), "MPa"), "Stress 2":((f"{stressUC[1]:.2f}", f"{stressC[1]:.2f}"), "MPa"),
###"Max N-R":(("-", f"{maxIterations}"),"-")}
###OutputTable(rows=data, columnTitles=("Uncorrected","Corrected"), title="Simulation Results")
# plot various solutions
plt.figure()
plt.title("Iterated Displacements")
plt.plot(uHistoryUC, vHistoryUC, 'bs', label="Uncorrected Method")
##plt.plot(uHistoryC, vHistoryC, 'r.', label="Corrected Method")
plt.xlabel("x displacement [mm]")
plt.ylabel("y displacement [mm]")
plt.legend()
#force vs displacement
plt.figure()
plt.title("Force vs Displacement")
plt.plot(deltaHistoryUC, forceHistoryUncorrected, 'bs', label="Uncorrected Iteration")
plt.xlabel("Total Deformation [mm]")
plt.ylabel("Vertial Force [kN]")
plt.legend()
plt.show()
