import numpy as np
import matplotlib.pyplot as plt

# parameters
c = .1 #m
L = 1 #m
P = 1 #kN

x, y = np.meshgrid(np.linspace(0,L,30), np.linspace(-c,c,30))
stress = 3*P/4/c**3*(x*y + np.sqrt((x*y)**2 + (c**2 - y**2)**2))
levels = np.linspace(stress.min(), stress.max(), 30)

fig1, ax1 = plt.subplots(constrained_layout=True)
ax1.set_title("Max Principal Stress")
plt.style.use('_mpl-gallery-nogrid')
StressContour = ax1.contourf(x, y, stress, levels=levels)
ax1.set_xlabel("x position (m)")
ax1.set_ylabel("y position (m)")
CB = fig1.colorbar(StressContour)
CB.set_label("Stress (MPa)")
plt.show()