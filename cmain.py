import matplotlib
matplotlib.use("TkAgg")
import numpy as np

# Try to use the Cython implementation, fallback to the Python implementation if it doesn't work
try:
    from cDEM import DEM
except:
    from DEM import DEM

def XuPayne():
    """Reproducing Figure 4 of Xu and Payne (2009)"""
    
    import matplotlib.pyplot as plt
    import Gassmann
    import time
    
    # Matrix properties
    # Source: Acoustic Properties of Carbonate Rocks p. 77
    Km_dol = 69.0 # GPa
    Gm_dol = 52.0 # GPa
    rhom_dol = 2.88 # g/cm3
    Km_lim = 77.0 # GPa
    Gm_lim = 32.0 # GPa
    rhom_lim = 2.71 # g/cm3
    
    # Fluid properties
    Kf = 3.0 # GPa
    rhof = 1.0 # g/cm3
    
    # Porosity
    phimax = 0.4
    
    # Aspect ratio of different pore types
    alpha_ref = 0.15
    alpha_crack = 0.02
    alpha_stiff = 0.8
    
    # Creating lists of inclusion configurations (alpha and volume for each inclusion type)
    crackandref_alphas = []
    crackandref_volumes = []
    stiffandref_alphas = []
    stiffandref_volumes = []
    fractions = [0.2, 0.4, 0.6, 0.8]
    for fraction in fractions:
        crackandref_alphas.append([alpha_crack, alpha_ref])
        crackandref_volumes.append([fraction*phimax, (1.0-fraction)*phimax])
        stiffandref_alphas.append([alpha_stiff, alpha_ref])
        stiffandref_volumes.append([fraction*phimax, (1.0-fraction)*phimax])
    
    alphas = [[alpha_crack]] + crackandref_alphas + [[alpha_ref]] + stiffandref_alphas + [[alpha_stiff]]
    volumes = [[phimax]] + crackandref_volumes + [[phimax]] + stiffandref_volumes + [[phimax]]
    
    # Adding one extra element for the dolomite line
    alphas += [[alpha_ref]]
    volumes += [[phimax]]
    
    # Creating lists of matrix proeprties
    Kms = len(alphas)*[Km_lim]
    Gms = len(alphas)*[Gm_lim]
    Kms[-1] = Km_dol
    Gms[-1] = Gm_dol
    
    # Creating style for each line
    crack_style = 'm:'
    ref_style = 'b-'
    stiff_style = 'r--'
    dol_style = 'g-'
    
    styles = [crack_style] + [crack_style]*len(crackandref_alphas) + [ref_style] + [stiff_style] + [stiff_style]*len(stiffandref_alphas)
    styles += [dol_style]
    
    # Plot configuration
    plt.figure(figsize=(7.0, 15.4/3.0), dpi=150)
    plt.subplots_adjust(left = 0.11, top = 0.88, bottom = 0.16, right = 0.96)
    plt.subplot(1, 1, 1)
    
    for inclusion_alphas, inclusion_volumes, style, Km, Gm in zip(alphas, volumes, styles, Kms, Gms):
        ni = len(inclusion_alphas)
        
        # Dry inclusions
        Kis = np.zeros(ni, dtype=float)
        Gis = np.zeros(ni, dtype=float)

        K, G, phi = DEM(Km, Gm, Kis, Gis, np.array(inclusion_alphas), np.array(inclusion_volumes))
        
        rho = (1.0 - phi)*rhom_lim + phi*rhof
        Ks = Gassmann.Ks(K, Km, Kf, phi)
        
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        
        plt.plot(phi, Vp, style)
    
    plt.ylabel("$V_p (km/s)$", fontsize=20)
    plt.xlabel("$\\phi$", fontsize=20)
    plt.ylim(2.0, 7.0)
    plt.xlim(0.0, phimax)
    
    plt.grid()
    plt.savefig("XuAndPayne2009_Figure5.png", dpi=150)

if __name__ == '__main__':
    XuPayne()