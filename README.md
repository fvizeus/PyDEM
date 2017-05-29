# PyDEM

A Python implementation of DEM (Differential Effective Medium) for rock physics applications.

## What is DEM?

Differential Effective Medium theory allows the calculation of elastic properties (namely bulk and shear moduli) of composite media [\[1\]][1]. This theory "constructs" the model as a background solid medium with added ellipsoidal inclusions. The inclusions are added in an iteractive manner to avoid violating the Kuster and Tuksoz (1974) [\[2\]][2] premises.

DEM is widely used for estimating wave-velocity in porous media. For example Xu and Payne (2009) [\[3\]][3] proposed a rock physics model for calculating wave velocities in carbonates, extending the model proposed by Xu and White (1995) [\[4\]][4]. Both of this models use DEM theory.

[1]: http://asa.scitation.org/doi/pdf/10.1121/1.385172
[2]: http://geophysics.geoscienceworld.org/content/39/5/587
[3]: http://library.seg.org/doi/abs/10.1190/1.3064148?journalCode=leedff
[4]: http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2478.1995.tb00126.x/abstract

## Implementation

The goal of this package is to allow the modelling of bulk and shear moduli using only the matrix properties (bulk and shear moduli) and inclusion properties (bulk and shear moduli, aspect ratio, concentration). Additional tools for velocities modelling are also provided (i.e. Gassmann fluid substitution), since in practice the velocities are used instead of the elastic moduli.

## Example (example.py file)

``` Python
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
```

