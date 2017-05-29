"""
From Berryman 1980
"""

import numpy as np

def _theta_scalar(alpha):
    if alpha == 1.0:
        return 2.0/3.0
    else:
        return alpha*(np.arccos(alpha) - alpha*np.sqrt(1.0 - alpha*alpha))/(1.0 - alpha*alpha)**(3.0/2.0)

def _theta_array(alpha):
    not1 = alpha != 1.0
    thetaarray = np.empty_like(alpha)
    thetaarray[~not1] = 2.0/3.0
    alphanot1 = alpha[not1]
    thetaarray[not1] = alphanot1*(np.arccos(alphanot1) - alphanot1*np.sqrt(1.0 - alphanot1*alphanot1))/(1.0 - alphanot1*alphanot1)**(3.0/2.0)
    return thetaarray

def theta2(alpha):
    if isinstance(alpha, np.ndarray):
        return _theta_array(alpha)
    else:
        return _theta_scalar(alpha)

def _f_scalar(alpha, theta=None):
    if alpha == 1.0:
        return -2.0/5.0
    if theta is None:
        theta = _theta_scalar(alpha)
    return alpha*alpha*(3.0*theta - 2.0)/(1.0 - alpha*alpha)

def _f_array(alpha, theta=None):
    not1 = alpha != 1.0
    farray = np.empty_like(alpha)
    farray[~not1] = -2.0/5.0
    alphanot1 = alpha[not1]
    if theta is None:
        thetanot1 = _theta_array(alphanot1)
    else:
        thetanot1 = theta[not1]
    farray[not1] = alphanot1*alphanot1*(3.0*thetanot1 - 2.0)/(1.0 - alphanot1*alphanot1)
    return farray

def f2(alpha, theta=None):
    if isinstance(alpha, np.ndarray):
        return _f_array(alpha, theta)
    else:
        return _f_scalar(alpha, theta)

def theta(alpha):
    return alpha*(np.arccos(alpha) - alpha*np.sqrt(1.0 - alpha*alpha))/(1.0 - alpha*alpha)**(3.0/2.0)

def f(alpha, theta):
    return alpha*alpha*(3.0*theta - 2.0)/(1.0 - alpha*alpha)

def PQ(A, B, R, theta, f):
    F1 = 1.0 + A*(1.5*(f + theta) - R*(1.5*f + 2.5*theta - 4.0/3.0))
    F2 = 1.0 + A*(1.0 + 1.5*(f + theta) - R*(1.5*f + 2.5*theta)) + B*(3.0 - 4.0*R) + A*(A + 3.0*B)*(1.5 - 2.0*R)*(f + theta - R*(f - theta + 2.0*theta*theta))
    F3 = 1.0 + A*(1.0 - f - 1.5*theta + R*(f + theta))
    F4 = 1.0 + (A/4.0)*(f + 3.0*theta - R*(f - theta))
    F5 = A*(-f + R*(f + theta - 4.0/3.0)) + B*theta*(3.0 - 4.0*R)
    F6 = 1.0 + A*(1.0 + f - R*(f + theta)) + B*(1.0 - theta)*(3.0 - 4.0*R)
    F7 = 2.0 + (A/4.0)*(3.0*f + 9.0*theta - R*(3.0*f + 5.0*theta)) + B*theta*(3.0 - 4.0*R)
    F8 = A*(1.0 - 2.0*R + (f/2.0)*(R - 1.0) + (theta/2.0)*(5.0*R - 3.0)) + B*(1.0 - theta)*(3.0 - 4.0*R)
    F9 = A*((R - 1.0)*f - R*theta) + B*theta*(3.0 - 4.0*R)
    
    P = 3.0*F1/F2
    Q = 2.0/F3 + 1.0/F4 + (F4*F5 + F6*F7 - F8*F9)/(F2*F4)
    return P, Q

def KG(Km, Gm, Ki, Gi, ci, theta, f):
    A = Gi/Gm - 1.0
    B = (Ki/Km - Gi/Gm)/3.0
    R = Gm/(Km + (4.0/3.0)*Gm)
    Fm = (Gm/6.0)*(9.0*Km + 8.0*Gm)/(Km + 2.0*Gm)
    
    P, Q = PQ(A, B, R, theta, f)

    K = Km - (Km + (4.0/3.0)*Gm)*ci*(Km - Ki)*P/3.0/(Km + (4.0/3.0)*Gm + ci*(Km - Ki)*P/3.0)
    G = Gm - (Gm + Fm)*ci*(Gm - Gi)*Q/5.0/(Gm + Fm + ci*(Gm - Gi)*Q/5.0)
    
    return K, G

def DEM(Km, Gm, Ki, Gi, alphai, phii, curphi=0.0, minalphaphiratio=100000):
    ni = np.ceil(minalphaphiratio*phii/alphai).astype(np.int)
    dphii = phii/ni
    indexes = np.repeat(np.arange(len(ni)), ni)
    
    np.random.shuffle(indexes)
    
    thetai = theta2(alphai)
    fi = f2(alphai, thetai)
    
    K = np.empty(len(indexes))
    G = np.empty(len(indexes))
    phi = np.empty(len(indexes))

    K_ = Km
    G_ = Gm
    phi_ = curphi   
    
    for j in range(len(indexes)):
        i = indexes[j]
        ci = dphii[i]/(1.0 - phi_)
        phi_ += dphii[i]
        K_, G_ = KG(K_, G_, Ki[i], Gi[i], ci, thetai[i], fi[i])
        
        K[j] = K_
        G[j] = G_
        phi[j] = phi_
     
    return K, G, phi

def DEM2(Km, Gm, Ki, Gi, alpha, phi, phi0=0.0, r=100):
    n = np.log((1.0 - phi)/(1.0 - phi0))/np.log(1.0 - alpha/r)
    n_ = int(np.ceil(n))
    r_ = alpha/(1.0 - ((1.0 - phi)/(1.0 - phi0))**(1.0/n_))
    c = alpha/r_
    
    print alpha, r_, n_, phi
    
    K = np.empty(n_)
    G = np.empty(n_)
    phi = np.empty(n_)
    
    thetai = theta2(alpha)
    fi = f2(alpha, thetai)

    K_ = Km
    G_ = Gm
    phi_ = phi0 
    
    for j in range(n_):
        dphi = c*(1.0 - phi0)*(1.0 - c)**j
        phi_ += dphi
        K_, G_ = KG(K_, G_, Ki, Gi, c, thetai, fi)
        
        K[j] = K_
        G[j] = G_
        phi[j] = phi_
     
    return K, G, phi

def DEM3(Km, Gm, Ki, Gi, alphai, phii, phi0=0.0, r=100):
    from scipy.optimize import fsolve
    
    phi = np.sum(phii)
    fraci = phii/np.sum(phi)
    ci = fraci*alphai/r
    n = int(np.ceil((np.log(1.0-phi)-np.log(1.0-phi0))/np.sum(np.log(1.0-ci))))
    m = len(alphai)
    def func(r):
        return m*np.log(r) - np.sum(np.log(r - fraci*alphai)) + (np.log(1.0 - phi) - np.log(1.0 - phi0))/n
    
    def fprime(r):
        return m/r - np.sum(1.0/(r - fraci*alphai))
    
    r_ = fsolve(func, r, fprime=fprime)[0]
    
    ci = fraci*alphai/r_
    # WHY, GOD?! WHY!!!???
    # ci = fraci/r_
    
    thetai = theta2(alphai)
    fi = f2(alphai, thetai)
    
    K = np.empty(n+1)
    K[0] = Km
    G = np.empty(n+1)
    G[0] = Gm
    phi = np.empty(n+1)
    phi[0] = phi0

    K_ = Km
    G_ = Gm
    phi_ = phi0
    
    for i in range(1,n+1):
        dphi = ci[0]*(1.0 - phi_)
        K_, G_ = KG(K_, G_, Ki[0], Gi[0], ci[0]/alphai[0], thetai[0], fi[0])
        phi_ += dphi
        for j in range(1, m):
            dphi *= ci[j]*(1.0 - ci[j-1])/ci[j-1]
            K_, G_ = KG(K_, G_, Ki[j], Gi[j], ci[j]/alphai[j], thetai[j], fi[j])
            phi_ += dphi
        K[i] = K_
        G[i] = G_
        phi[i] = phi_
    
    print phi_
     
    return K, G, phi

from scipy.optimize import fsolve
def DEM4(Km, Gm, Ki, Gi, alphai, phii, phi0=0.0, r=1000, phitol=1.0E-10, gamma=0.01):
    phi = np.sum(phii)
    fraci = phii/np.sum(phi)
    ci = fraci*alphai/r
    n = int(np.ceil((np.log(1.0-phi)-np.log(1.0-phi0))/np.sum(np.log(1.0-ci))))
    m = len(alphai)
    
    def func(r):
        f = np.empty(m)
        f[0] = np.log(alphai[0]/r[0]) + np.log(1.0 - phi0/phi) - np.log(1 - ((1.0 - phi)/(1.0 - phi0))**(1.0/n))
        for j in range(1, m):
            f[j] = f[j-1] + np.log(alphai[j]/r[j]) + np.log(r[j-1]/alphai[j-1] - fraci[j-1])
        return f
    
    def fprime(r):
        jac = np.diag(-1.0/r)
        for j in range(0, m-1):
            jac[j+1:, j] = -1.0/r[j] + 1.0/(r[j] - fraci[j]*alphai[j])
        
        return jac
    
    r0 = r*np.ones(m)
    
    ri = fsolve(func, r0, fprime=fprime, factor=0.1)
    
    ci = fraci*alphai/ri
    
    # thetai = theta2(alphai)
    # fi = f2(alphai, thetai)
    thetai = theta(alphai)
    fi = f(alphai, thetai)
    
    K = np.empty(n)
    G = np.empty(n)
    phi = np.empty(n)

    K_ = Km
    G_ = Gm
    phi_ = phi0
    
    for i in range(n):
        dphi = ci[0]*(1.0 - phi_)
        K_, G_ = KG(K_, G_, Ki[0], Gi[0], ci[0], thetai[0], fi[0])
        phi_ += dphi
        for j in range(1, m):
            dphi *= ci[j]*(1.0 - ci[j-1])/ci[j-1]
            K_, G_ = KG(K_, G_, Ki[j], Gi[j], ci[j], thetai[j], fi[j])
            phi_ += dphi
        K[i] = K_
        G[i] = G_
        phi[i] = phi_
    
    return K, G, phi
    
def discrete_normal(zmin, zmax, n):
    from scipy.stats import norm
    z = np.linspace(zmin, zmax, n)
    x = np.empty(n + 1)
    x[0] = np.exp(-z[0]**2/2.0)/norm.cdf(z[0])
    x[-1] = np.exp(-z[-1]**2/2.0)/(1.0 - norm.cdf(z[-1]))
    x[1:-1] = (np.exp(-z[:-1]**2/2.0) - np.exp(-z[1:]**2/2.0))/(norm.cdf(z[1:]) - norm.cdf(z[:-1]))
    x /= np.sqrt(2.0*np.pi)

    p = np.empty(len(z) + 1)
    p[0] = norm.cdf(z[0])
    p[-1] = 1.0 - norm.cdf(z[-1])
    p[1:-1] = norm.cdf(z[1:]) - norm.cdf(z[:-1])
    
    return x, p
