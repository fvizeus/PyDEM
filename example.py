import numpy as np
import matplotlib.pyplot as plt

# Cython implementation
# from cDEM import DEM
# Python implementation
from DEM import DEM

import Gassmann

# Matrix properties
Km = 77.0 # GPa
Gm = 32.0 # GPa
rhom = 2.71 # g/cm3

# Fluid properties
Kf = 3.0 # GPa
rhof = 1.0 # g/cm3

# Porosity
phimax = 0.4

# Inclusion properties
# In this example a mixture of three inclusion types are used:
# - 30% of 0.02 aspect ratio
# - 50% of 0.15 aspect ratio
# - 20% of 0.80 aspect ratio
alphas = np.array([0.01, 0.15, 0.8])
volumes = np.array([0.3, 0.5, 0.2])*phimax

# Dry inclusions
Kis = np.zeros(len(alphas), dtype=float)
Gis = np.zeros(len(alphas), dtype=float)

# The DEM function returns the bulk and shear moduli along with the porosity array to match them.
# The porosity array is not regularly spaced. If you need so, you should reinterpolate.
K, G, phi = DEM(Km, Gm, Kis, Gis, alphas, volumes)

rho = (1.0 - phi)*rhom + phi*rhof
Ks = Gassmann.Ks(K, Km, Kf, phi)

Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
Vs = np.sqrt(G/rho)

plt.plot(phi, Vp, 'b', label='Vp')
plt.plot(phi, Vs, 'g', label='Vs')

plt.show()
