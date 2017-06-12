import matplotlib
matplotlib.use("TkAgg")
import numpy as np

from DEM import DEM

def XuPayne():
    # Reproducing Xu and Payne 2009
    import matplotlib.pyplot as plt
    import Gassmann
    import time
    # Source: Acoustic Properties of Carbonate Rocks p. 77
    Km_dol = 69.0 # GPa
    Gm_dol = 52.0 # GPa
    rhom_dol = 2.88 # g/cm3
    Km_lim = 77.0 # GPa
    Gm_lim = 32.0 # GPa
    rhom_lim = 2.71 # g/cm3
    
    Kf = 3.0 # GPa
    rhof = 1.0 # g/cm3
    
    phimax = 0.4
    
    alpha_ref = 0.15
    alpha_crack = 0.02
    alpha_stiff = 0.8
    
    cracks_alphas = []
    cracks_phis = []
    for frac in [0.2, 0.4, 0.6, 0.8]:
        cracks_alphas.append(np.array([alpha_crack, alpha_ref]))
        cracks_phis.append(np.array([frac*phimax, (1.0-frac)*phimax]))
    cracks_alphas.append(np.array([alpha_crack]))
    cracks_phis.append(np.array([phimax]))
    
    ref_alphas = [np.array([alpha_ref])]
    ref_phis = [np.array([phimax])]
    
    stiffs_alphas = []
    stiffs_phis = []
    for frac in [0.2, 0.4, 0.6, 0.8]:
        stiffs_alphas.append(np.array([alpha_stiff, alpha_ref]))
        stiffs_phis.append(np.array([frac*phimax, (1.0-frac)*phimax]))
    stiffs_alphas.append(np.array([alpha_stiff]))
    stiffs_phis.append(np.array([phimax]))
    
    plt.figure(figsize=(7.0, 14.0/3.0), dpi=150)
    plt.subplots_adjust(left = 0.12, top = 0.88, bottom = 0.15, right = 0.9)
    # plt.figure(figsize=(7.0, 15.4/3.0), dpi=150)
    # plt.subplots_adjust(left = 0.11, top = 0.88, bottom = 0.16, right = 0.96)
    plt.subplot(1, 1, 1)
    
    allalphas = cracks_alphas + ref_alphas + stiffs_alphas
    allphis = cracks_phis + ref_phis + stiffs_phis
    styles = ['m:']*len(cracks_alphas) + ['b-']*len(ref_alphas) + ['r--']*len(stiffs_alphas)
    
    for alphas, phis, style in zip(allalphas, allphis, styles):
        ni = len(alphas)
        K, G, phi = DEM(Km_lim, Gm_lim, np.zeros(ni, dtype=float), np.zeros(ni, dtype=float), alphas, phis)
        phistep = phi[1:] - phi[:-1]
        print phi[0], phistep[0]/phistep[-1], len(phistep)
        
        rho = (1.0 - phi)*rhom_lim + phi*rhof
        Ks = Gassmann.Ks(K, Km_lim, Kf, phi)
        
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/rho)
        
        plt.plot(phi, Vp, style)
    
    for alphas, phis in zip(ref_alphas, ref_phis):
        ni = len(alphas)
        K, G, phi = DEM(Km_dol, Gm_dol, np.zeros(ni, dtype=float), np.zeros(ni, dtype=float), alphas, phis)
        
        rho = (1.0 - phi)*rhom_lim + phi*rhof
        Ks = Gassmann.Ks(K, Km_lim, Kf, phi)
        
        Vp = np.sqrt((Ks + 4.0*G/3.0)/rho)
        Vs = np.sqrt(G/rho)
        
        plt.plot(phi, Vp, 'g-')
    
    # data = np.loadtxt("data.txt")
    # vplog, philog, litolog = data.T
    
    # plt.scatter(philog[litolog == 1.0], vplog[litolog == 1.0], c='b', s=5.0, zorder=-1000, linewidths=0.25)
    # plt.scatter(philog[litolog == 30.0], vplog[litolog == 30.0], c='g')
    
    plt.ylabel("$V_p (km/s)$", fontsize=20)
    plt.xlabel("$\\phi$", fontsize=20)
    # plt.ylim(2.0, 7.0)
    plt.ylim(1.0, 7.0)
    plt.xlim(0.0, phimax)
    
    plt.grid()
    plt.savefig("XuPayne1a.png", dpi=150)
    plt.show()

if __name__ == '__main__':
    XuPayne()