from __future__ import division
import numpy as np
cimport numpy as np

DTYPE = np.double
ctypedef np.double_t DTYPE_t

cdef theta(np.ndarray[DTYPE_t] alpha):
    return alpha*(np.arccos(alpha) - alpha*np.sqrt(1.0 - alpha*alpha))/(1.0 - alpha*alpha)**(3.0/2.0)

cdef f(np.ndarray[DTYPE_t] alpha, np.ndarray[DTYPE_t] theta):
    return alpha*alpha*(3.0*theta - 2.0)/(1.0 - alpha*alpha)

cdef PQ(DTYPE_t A, DTYPE_t B, DTYPE_t R, DTYPE_t theta, DTYPE_t f):
    cdef DTYPE_t F1, F2, F3, F4, F5, F6, F7, F8, F9, P, Q
    
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

cdef KG(DTYPE_t Km, DTYPE_t Gm, DTYPE_t Ki, DTYPE_t Gi, DTYPE_t ci, DTYPE_t theta, DTYPE_t f):
    cdef DTYPE_t A, B, R, Fm, P, Q, K, G 
    
    A = Gi/Gm - 1.0
    B = (Ki/Km - Gi/Gm)/3.0
    R = Gm/(Km + (4.0/3.0)*Gm)
    Fm = (Gm/6.0)*(9.0*Km + 8.0*Gm)/(Km + 2.0*Gm)
    
    P, Q = PQ(A, B, R, theta, f)

    K = Km - (Km + (4.0/3.0)*Gm)*ci*(Km - Ki)*P/3.0/(Km + (4.0/3.0)*Gm + ci*(Km - Ki)*P/3.0)
    G = Gm - (Gm + Fm)*ci*(Gm - Gi)*Q/5.0/(Gm + Fm + ci*(Gm - Gi)*Q/5.0)
    
    return K, G

cpdef DEM(DTYPE_t Km, DTYPE_t Gm, np.ndarray[DTYPE_t] Ki, np.ndarray[DTYPE_t] Gi, np.ndarray[DTYPE_t] alphai, np.ndarray[DTYPE_t] phii, DTYPE_t phi0=0.0, DTYPE_t r=1000.0, DTYPE_t phitol=1.0E-10, DTYPE_t gamma=0.01):
    cdef np.ndarray[DTYPE_t] fval, fraci, ci, ri, fi, thetai, K, G, phi
    cdef np.ndarray[DTYPE_t, ndim=2] fjac
    cdef DTYPE_t phin, K_, G_, phi_, dphi, phic
    cdef Py_ssize_t n, m, i, j
    
    phin = np.sum(phii)
    fraci = phii/phin
    ci = fraci*alphai/r
    
    n = np.int(np.ceil((np.log(1.0-phin)-np.log(1.0-phi0))/np.sum(np.log(1.0-ci))))
    m = alphai.shape[0]
    
    phic = 1.0 - (1.0 - phi0)*np.exp(n*np.sum(np.log(1.0 - ci)))
    
    fval = np.zeros(m)
    fjac = np.zeros((m, m))
    ri = np.zeros(m)
    ri += r
    
    while np.abs(phic - phin) > phitol:
        fval[0] = np.log(alphai[0]/ri[0]) + np.log(1.0 - phi0/phin) - np.log(1 - ((1.0 - phin)/(1.0 - phi0))**(1.0/n))
        for j in range(1, m):
            fval[j] = fval[j-1] + np.log(alphai[j]/ri[j]) + np.log(ri[j-1]/alphai[j-1] - fraci[j-1])
        
        for j in range(m):
            fjac[j, j] = -1.0/ri[j]
        
        for j in range(m-1):
            fjac[j+1:, j] = -1.0/ri[j] + 1.0/(ri[j] - fraci[j]*alphai[j])
        
        step = gamma*np.dot(np.linalg.inv(fjac), fval)
        
        ri -= step
        
        ci = fraci*alphai/ri
        
        phic = 1.0 - (1.0 - phi0)*np.exp(n*np.sum(np.log(1.0 - ci)))
    
    thetai = theta(alphai)
    fi = f(alphai, thetai)
    
    K = np.zeros(n)
    G = np.zeros(n)
    phi = np.zeros(n)
    
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
