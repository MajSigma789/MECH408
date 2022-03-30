##
# Aidan Moretz
# MECH 408
# Nonlinear Cable Handout
# Created: 1/26/2022
# Modified: 3/25/2022
# Uses a corrected incremental approach to solve geometric 
# and material nonlinearities in cable deformation
##

import numpy as np
import matplotlib.pyplot as plt
from tkinter.ttk import *

#TODO: reformat VariableTable to use given title list and dimensions
class VariableTable(Frame):
    def __init__(self, master= None):
        Frame.__init__(self,master)
        self.grid()
        self.create_inputs()

    def create_inputs(self):
        """Instantiates and lays out the input table fields and buttons"""
        titles = {"":None, "L":"[mm]", "E":"[GPa]", "A":"[mm^2]", "alpha":"[rad]","Force(x,y)":"[kN]"}
        tablewidth = 3
        # store entry values in dictionary for later
        self.values = {}
        row = -1
        for title, unit in titles.items():
            row += 1
            # top row with no entry points
            if title == "":
                Label(self, width=7, text="Cable 1").grid(row=row, column=1)
                Label(self, width=7, text="Cable 2").grid(row=row, column=2)
                continue
            # normal row
            Label(self, width = 15, text=" ".join([title,unit])).grid(row=row, column=0)
            self.values[title] = []
            for column in range(1,tablewidth):
                self.values[title].append(Entry(self, width=7))
                self.values[title][-1].grid(row=row, column=column)
        # last row with only one entry
        Label(self, width=15, text="increments").grid(row=len(titles), column=0)
        self.values["increments"] = Entry(self, width=7)
        self.values["increments"].grid(row=len(titles), column=1)
        # button to close
        Button(self, width=7, text="Run", command=lambda: self.quit()).grid(row=len(titles)+1, column=tablewidth-1)

#LOWPRIO: setup problem as a child of abstract Solvable, separate system from solver methods
def tangentStiffness(position):
    """Forms the tangent stiffness matrix of the two cable system"""
    Kt = np.zeros(shape=(2,2))
    L1Offset = (L[0] + position[0])**(-1)
    L2SinOffset = L[1]*np.sin(alpha[1]) + position[1]
    L2CosOffset = L[1]*np.cos(alpha[1]) - position[0]
    w1 = np.sqrt(position[1]**2/(L1Offset)**(-2))
    w2 = (L2SinOffset)/(L2CosOffset)
    w1Term = (1 + w1**2)**(-1)
    w2Term = (1 + w2**2)**(-1)
    deltas = deltaCalculation(position=position)
    # account for material non-linearities
    if materials is True:
        # find new stiffness derivative and stiffness constant
        stresses = getStresses(delta=deltas)
        updateElementStiffness(stress = stresses)
        stiffnessDerivatives = stiffnessDerivative(stress=stresses)
        # modify tangent stiffness for nonlinearity
        firstNonlinearTerm = stiffnessDerivatives[0]*deltas[0]*w1Term
        secondNonlinearTerm = stiffnessDerivatives[1]*deltas[1]*w2Term
        Kt += np.array([
            [firstNonlinearTerm + secondNonlinearTerm,
            firstNonlinearTerm*w1 - secondNonlinearTerm*w2],
            [firstNonlinearTerm*w1 - secondNonlinearTerm*w2,
            firstNonlinearTerm*w1*w1 + secondNonlinearTerm*w2*w2]
        ])
    # linear material tangent stiffness matrix
    Kt += np.array([
        [k[0]*(w1Term + deltas[0]*w1*w1Term**(1.5)*position[1]*L1Offset**2) + k[1]*(w2Term + deltas[1]*w2*w2Term**1.5*(L2SinOffset)*(L2CosOffset)**(-2)), # dPx/du
        k[0]*(w1*w1Term - deltas[0]*w1*w1Term**(1.5)*L1Offset) - k[1]*(w2*w2Term - deltas[1]*w2*w2Term**1.5*(L[0]*np.cos(alpha[1]) - position[0])**(-1))], # dPx/dv
        [k[0]*(w1*w1Term + deltas[0]*w1Term**(1.5)*position[1]*L1Offset**2) - k[1]*(w2*w2Term - deltas[1]*w2Term**1.5*(L2SinOffset)*(L2CosOffset)**(-2)), # dPy/du
        k[0]*(w1**2*w1Term + deltas[0]*w1Term**(1.5)*L1Offset) + k[1]*(w2**2*w2Term + deltas[1]*w2Term**1.5*(L2CosOffset)**(-1))] # dPy/dv
        ])
    return Kt

def deltaCalculation(position):
    """Calculates the change in cable lengths for the current displacement"""
    delta1 = np.sqrt((L[0] + position[0])**2 + position[1]**2) - L[0]
    delta2 = np.sqrt((L[1]*np.cos(alpha[1]) - position[0])**2 + (L[1]*np.sin(alpha[1]) + position[1])**2) - L[1]
    return np.array((delta1, delta2))

# TODO: try more iterations for correction with material non-linearity
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

def incrementalLinearization(position, load):
    """Advances the linearization by one step to find the change in position
    Parameters: position [mm], load [kN]
    Returns: deltaX [mm], Kt [kN/mm]"""
    Kt = tangentStiffness(position) #kN/mm
    deltaX = np.linalg.solve(Kt, load) #mm
    return deltaX, Kt

def getStresses(delta):
    """Calculates stress at a given deformation using the current element stiffness
    Parameters: delta [mm]
    Returns: stress [GPa = kN/mm^2]"""
    return k*delta/A #kN/mm^2 = GPa

def getStrainRBOG(stress):
    """Calculates the strain for a given stress using the Ramberg-Osgood relationship
    Parameters: stress [GPa]
    Returns: strain [-]"""
    return stress/E + .002*(stress/yieldStress)**(1/n) #-


def stiffnessDerivative(stress):
    """Calculates the instantaneous rate of change the element stiffness for a change in displacement
    Parameters: stress [GPa]
    Returns: dK/dDelta [kN/mm^2 = GPa]"""
    strain = getStrainRBOG(stress) # -
    if np.array_equiv(strain, np.zeros_like(strain)):
        return np.zeros_like(stress)
    strainDeriv = 1/E + .002/yieldStress/n*(stress/yieldStress)**(1/n - 1) #GPa^-1
    stiffDeriv = A/L/L/strain*(1/strainDeriv - stress/strain) #GPa = kN/mm^2
    return stiffDeriv

def updateElementStiffness(stress):
    """Updates the element stiffness for the current loading
    Parameters: stress [GPa]
    Returns: k [kN/mm] at reference"""
    global k
    # initial increment is elastic
    if np.array_equiv(stress, np.zeros_like(stress)):
        k = E*A/L
        return
    strain = getStrainRBOG(stress)
    k = A/L*stress/strain


# main method
#LOWPRIO: take material properties as input window
yieldStress = .11 #kN/mm^2
n = .1051
materials = True
# create input window
inTable = VariableTable()
inTable.master.title("System Parameters")
inTable.mainloop()

# get input parameters
L = np.array([float(entry.get()) for entry in inTable.values["L"]]) #mm
E = np.array([float(entry.get()) for entry in inTable.values["E"]]) #GPa = kN/mm^2
A = np.array([float(entry.get()) for entry in inTable.values["A"]]) #mm^2
alpha = np.array([float(entry.get()) for entry in inTable.values["alpha"]]) #rad
totalLoad = np.array([float(entry.get()) for entry in inTable.values["Force(x,y)"]]) #kN
increments = int(inTable.values["increments"].get())

# derived parameters
k = E*A/L #kN/mm
dLoad = totalLoad/increments #kN

# variables
beta = np.copy(alpha) #rad
loadIncrements = [dLoad]*increments #kN

# uncorrected solution
positionHistoryUncorrected = ([0.,], [0.,])
for dP in loadIncrements:
    # step forward in loading amount
    position = np.array([history[-1] for history in positionHistoryUncorrected])
    deltaPos = incrementalLinearization(position=position, load=dP)[0]
    position += deltaPos
    # add component positions to history
    for history, displace in zip(positionHistoryUncorrected,position): history.append(displace)
finalDelta = deltaCalculation(position=position)
stressUC = 1000*getStresses(delta=finalDelta) #MPa

# corrected solution
positionHistoryCorrected = ([0.,], [0.,])
maxIterations = 0
totalLoad = 0
for dP in loadIncrements:
    # step forward in loading amount
    position = np.array([history[-1] for history in positionHistoryCorrected])
    deltaPos, Kt = incrementalLinearization(position=position,load=dP)
    totalLoad += dP
    position = newtonRaphsonCorrection(position=position+deltaPos, load=totalLoad, stiffness=Kt)
    for history, displace in zip(positionHistoryCorrected,position): history.append(displace)
finalDelta = deltaCalculation(position=position)
stressC = 1000*getStresses(delta=finalDelta) #MPa

#LOWPRIO: TK formatted result table
# report uncorrected info
uHistoryUC, vHistoryUC = positionHistoryUncorrected[0], positionHistoryUncorrected[1]
uHistoryC, vHistoryC = positionHistoryCorrected[0], positionHistoryCorrected[1]
print("Final Uncorrected X: {:.3f}mm".format(uHistoryUC[-1]))
print("Final Uncorrected Y: {:.3f}mm".format(vHistoryUC[-1]))
print("Uncorrected Stress 1: {:.2f} MPa".format(stressUC[0]))
print("Uncorrected Stress 2: {:.2f} MPa".format(stressUC[1]))
# report corrected info
print("Final Corrected X: {:.3f}mm".format(uHistoryC[-1]))
print("Final Corrected Y: {:.3f}mm".format(vHistoryC[-1]))
print("Corrected Stress 1: {:.2f} MPa".format(stressC[0]))
print("Corrected Stress 2: {:.2f} MPa".format(stressC[1]))
print("Max N-R Iterations: {:d}".format(maxIterations))
# plot various solutions
plt.figure()
plt.title("Iterated Displacements")
plt.plot(uHistoryUC, vHistoryUC, 'bs', label="Uncorrected Method")
plt.plot(uHistoryC, vHistoryC, 'r.', label="Corrected Method")
plt.xlabel("x displacement [mm]")
plt.ylabel("y displacement [mm]")
plt.legend()
plt.show()
