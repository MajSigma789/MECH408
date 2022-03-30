import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipe, ellipeinc, ellipk, ellipkinc

EI = 200E9*(10E-3)**4/12
L = 1

angle = np.linspace(-np.pi/100,np.pi/100, 100000)
sign = np.sign(angle)
#angle = np.abs(angle)
m = .5*(1 - np.cos(angle))
phi = np.arcsin((np.sqrt(2*m))**(-1)) #rad
mDerivTheta = -.5*np.cos(angle)
phiDerivM = -.5/np.sqrt(8*m**3 - 4*m**2)
kDerivM = ellipe(m)/(m*(1-m**2)) - ellipk(m)/m
kPhiDerivM = ellipeinc(phi, m)/(m*(1-m**2)) - ellipkinc(phi, m)/m
kPhiDerivPhi = (1 - m*(np.sin(phi)**2))**(-1)
derivative = 2*sign*EI/L/L*(ellipk(m) - ellipkinc(phi, m))*(kDerivM*mDerivTheta - kPhiDerivM*mDerivTheta - kPhiDerivPhi*phiDerivM*mDerivTheta)


#angle *= sign
plt.figure()
plt.plot(angle, derivative)
plt.show()

print(derivative)