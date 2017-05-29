# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:35:15 2013

@author: fsantos
"""

import KusterToksoz as KT

_PFACTOR = 10


def XuWhite(Km, Gm, alpha, phi, curphi=0.0):
    nphi = int(phi*_PFACTOR/alpha + 0.5)
    if not nphi:
        nphi = 1
    dphi = phi/nphi

    K_ = Km
    G_ = Gm
    Kd = Km
    Gd = Gm

    for i in range(nphi):
        Kd = KT.Kd(K_, G_, alpha, dphi/(1.0 - curphi - i*dphi))
        Gd = KT.Gd(K_, G_, alpha, dphi/(1.0 - curphi - i*dphi))
        K_ = Kd
        G_ = Gd

    return Kd, Gd


if __name__ == '__maine__':
    import numpy as np

    Vpm = 5.85E3
    Vsm = 3.9E3
    rhom = 2.65E3

    Vpf = 1.6E3
    rhof = 1.1E3

    Km = rhom*(Vpm*Vpm - 4.0*Vsm*Vsm/3.0)
    Gm = rhom*Vsm*Vsm
    Kf = rhof*Vpf*Vpf

    phi = np.linspace(0.01, 0.99, 99)
    alpha = np.linspace(0.01, 0.99, 99)

    erroKd = np.zeros((99, 99))
    erroGd = np.zeros((99, 99))

    for i in range(99):
        for j in range(99):
            p = KT.T(Km, Gm, alpha[j])
            q = KT.F(Km, Gm, alpha[j])
            Kd_ = Km*(1.0 - phi[i])**p
            Gd_ = Gm*(1.0 - phi[i])**q
            Kd, Gd = XuWhite(Km, Gm, alpha[j], phi[i])
            erroKd[i][j] = 100.0*abs(Kd-Kd_)/Kd
            erroGd[i][j] = 100.0*abs(Gd-Gd_)/Gd
            print phi[i], alpha[j], erroKd[i][j]

    f = open('erroskd.txt', 'w')
    for i in range(99):
        for j in range(99):
            print >>f, erroKd[i][j],
        print >>f, ' '
    f.close()

    f = open('errosgd.txt', 'w')
    for i in range(99):
        for j in range(99):
            print >>f, erroGd[i][j],
        print >>f, ' '
    f.close()


if __name__ == '__main__':
# As seen on Xu & Keys (2002)
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams['figure.facecolor'] = 'white'

    Vpm = 5.85E3
    Vsm = 3.9E3
    rhom = 2.65E3

    Vpf = 1.6E3
    rhof = 1.1E3

    Km = rhom*(Vpm*Vpm - 4.0*Vsm*Vsm/3.0)
    Gm = rhom*Vsm*Vsm
    Kf = rhof*Vpf*Vpf

    phis = [0.01 + 0.01*i for i in range(40)]
    alphas = [0.01, 0.05, 0.10, 0.15]

    Vpss = []
    Vsss = []
    Vpss_ = []
    Vsss_ = []
    for alpha in alphas:
        Vps = []
        Vss = []
        Vps_ = []
        Vss_ = []
        for phi in phis:
            p = KT.T(Km, Gm, alpha)
            q = KT.F(Km, Gm, alpha)
            Kd_ = Km*(1.0 - phi)**p
            Gd_ = Gm*(1.0 - phi)**q
            Kd, Gd = XuWhite(Km, Gm, alpha, phi)

            rho = phi*rhof + (1.0 - phi)*rhom
            K_ = Kd_ + (1.0 - Kd_/Km)**2/(phi/Kf + (1.0 - phi)/Km - Kd_/Km/Km)
            K = Kd + (1.0 - Kd/Km)**2/(phi/Kf + (1.0 - phi)/Km - Kd/Km/Km)

            Vp = ((K + 4.0*Gd/3.0)/rho)**0.5
            Vs = (Gd/rho)**0.5
            Vp_ = ((K_ + 4.0*Gd/3.0)/rho)**0.5
            Vs_ = (Gd_/rho)**0.5

            Vps.append(Vp)
            Vss.append(Vs)
            Vps_.append(Vp_)
            Vss_.append(Vs_)

        Vpss.append(Vps)
        Vsss.append(Vss)
        Vpss_.append(Vps_)
        Vsss_.append(Vss_)

    plt.subplot(121)
    for V in Vpss:
        plt.plot(phis, V, 'k')

    for V in Vpss_:
        p2 = plt.plot(phis, V, 'ko', alpha=0.667)

    plt.xlabel('$\phi$', fontsize=24)
    plt.ylabel('$V_{p}\,[m/s]$', fontsize=24)
    bbox = {'fc': 'white', 'ec': 'white'}
    plt.annotate(r'$\alpha = 0.01$', (0.07, 3000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.annotate(r'$\alpha = 0.05$', (0.10, 4000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.annotate(r'$\alpha = 0.10$', (0.30, 3000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.annotate(r'$\alpha = 0.15$', (0.23, 4000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.legend([plt.gca().lines[1], plt.gca().lines[5]], ['Xu & White(1995)', 'Xu & Keys(2001)'], numpoints=1)
    plt.xlim(0, 0.4)
    plt.ylim(0, 6000)
    plt.grid()
    plt.xticks([0.1, 0.2, 0.3])
    plt.gca().set_aspect(0.4/6000)

    plt.subplot(122)
    for V in Vsss:
        plt.plot(phis, V, 'k')

    for V in Vsss_:
        plt.plot(phis, V, 'ko', alpha=0.667)

    plt.xlabel('$\phi$', fontsize=24)
    plt.ylabel('$V_{s}\,[m/s]$', fontsize=24)
    bbox = {'fc': 'white', 'ec': 'white'}
    plt.annotate(r'$\alpha = 0.01$', (0.06, 1000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.annotate(r'$\alpha = 0.05$', (0.14, 2000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.annotate(r'$\alpha = 0.10$', (0.25, 2000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.annotate(r'$\alpha = 0.15$', (0.33, 2000), horizontalalignment='center', fontsize=16, bbox=bbox)
    plt.legend([plt.gca().lines[1], plt.gca().lines[5]], ['Xu & White(1995)', 'Xu & Keys(2001)'], numpoints=1)
    plt.xlim(0, 0.4)
    plt.ylim(0, 6000)
    plt.grid()
    plt.xticks([0.1, 0.2, 0.3])
    plt.gca().set_aspect(0.4/6000)

    plt.show()
