###
# Aidan Moretz
# MECH 408
# Problem 7.13
# Date Created: 4/26/2022
# Last Modified: 4/26/2022
###

import numpy as np
import matplotlib.pyplot as plt
from MechCustom import InputTable
from typing import Optional


def main():
    # take inputs
    rows = {"nu":("-", 1), "y/l":("-", 1), "resolution":("-", 1)}
    inTable = InputTable(rows=rows, title="System Parameters", button="Run!")
    [nu, height, resolution] = [values[0] for values in inTable.getValues(entries=rows.keys())]
    resolution = int(resolution)

    # find vertical displacements
    x = np.linspace(-1,1, resolution)
    stress = planeStressVertical(position=x, nu=nu, height=height)
    strain = planeStrainVertical(position=x, nu=nu, height=height)

    # plot results
    plt.figure()
    plt.title(f"Normalized Vertical Displacement\ny* = {height}")
    plt.plot(x, stress, label="Plane Stress")
    plt.plot(x, strain, label="Plane Strain")
    plt.legend()
    plt.xlabel("Normalized Position [-]")
    plt.xlim((-1,1))
    plt.ylabel("Normalized Displacement [-]")
    plt.show()


def planeStressVertical(position, nu, height: Optional[float]=0):
    v = .5*(nu*height**2 + position**2 - 1)
    return v

def planeStrainVertical(position: float, nu: float, height: Optional[float]=0):
    v = .5*(1-nu**2)*(nu/(1-nu)*height**2 + position**2 - 1)
    return v

if __name__ == "__main__":
    main()