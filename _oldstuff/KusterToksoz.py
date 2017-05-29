# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:25:07 2013

@author: fsantos
"""

import numpy as np

__all__ = ['T', 'F', 'Kd', 'Gd']


def T2(K, G, a, Ki=0.0, Gi=0.0):
    A = Gi/G - 1.0
    B = (Ki/K - Gi/G)/3.0
    R = 1.0/(K/G + 4.0/3.0)
    C1 = 1.0 - a*a
    t = a*arccos(a)/C1/sqrt(C1) - a*a/(C1)
    f = a*a*(3.0*t - 2.0)/C1

    F1 = 1.0 + A*(1.5*(f + t) - R*(1.5*f + 2.5*t - 4.0/3.0))
    F2 = 1.0 + A*(1.0 + 1.5*(f + t) - R*(1.5*f + 2.5*t)) + B*(3.0 - 4.0*R) \
        + 0.5*A*(A + 3.0*B)*(3.0 - 4.0*R)*(f + t - R*(f - t + 2.0*t*t))

    return F1/F2


def F2(K, G, a, Ki=0.0, Gi=0.0):
    A = Gi/G - 1.0
    B = (Ki/K - Gi/G)/3.0
    R = 1.0/(K/G + 4.0/3.0)
    C1 = 1.0 - a*a
    t = a*arccos(a)/C1/sqrt(C1) - a*a/(C1)
    f = a*a*(3.0*t - 2.0)/C1

    F2 = 1.0 + A*(1.0 + 1.5*(f + t) - R*(1.5*f + 2.5*t)) + B*(3.0 - 4.0*R) \
        + 0.5*A*(A + 3.0*B)*(3.0 - 4.0*R)*(f + t - R*(f - t + 2.0*t*t))
    F3 = 1.0 + A*(1.0 - (f + 1.5*t) + R*(f + t))
    F4 = 1.0 + 0.25*A*(f + 3.0*t - R*(f - t))
    F5 = A*(-f + R*(f + t - 4.0/3.0)) + B*t*(3.0 - 4.0*R)
    F6 = 1.0 + A*(1.0 + f - R*(f + t)) + B*(1.0 - t)*(3.0 - 4.0*R)
    F7 = 2.0 + 0.25*A*(3.0*f + 9.0*t - R*(3.0*f + 5.0*t)) + B*t*(3.0 - 4.0*R)
    F8 = A*(1.0 - 2.0*R + 0.5*f*(R - 1.0) + 0.5*t*(5.0*R - 3.0)) \
        + B*(1.0 - t)*(3.0 - 4*R)
    F9 = A*((R - 1.0)*f - R*t) + B*t*(3.0 - 4.0*R)

    return (2.0/F3 + 1.0/F4 + (F4*F5 + F6*F7 - F8*F9)/F2/F4)/5.0


def T(K, G, a):  # Tiijj/3
    R = 1.0/(K/G + 4.0/3.0)
    C0 = 1.0 - a*a
    t = a*np.arccos(a)/C0/np.sqrt(C0) - a*a/(C0)
    f = a*a*(3.0*t - 2.0)/C0

    C1 = (R - 1.0)*(f + t)
    C2 = (R - 1.0)*(f - t)

    F1 = 1.0 - 4.0*R/3.0 + R*t + 1.5*C1
    F2 = R*(t*t*(4.0*R - 3.0) + 2.0*C2)

    return F1/F2


def F(K, G, a):  # (Tijij - Tiijj/3)/5
    R = 1.0/(K/G + 4.0/3.0)
    C0 = 1.0 - a*a
    t = a*np.arccos(a)/C0/np.sqrt(C0) - a*a/(C0)
    f = a*a*(3.0*t - 2.0)/C0

    C1 = (R - 1.0)*(f + t)
    C2 = (R - 1.0)*(f - t)

    F2 = R*(t*t*(4.0*R - 3.0) + 2.0*C2)
    F3 = 0.5*t - C1
    F4 = 1.0 - t + 0.25*C2
    F5 = 4.0*R/3.0 - t - C1
    F6 = t + C1
    F7 = 2.0 + t*(2.0*R - 3.0) + 0.75*C2
    F8 = (1.0 - t)*(2.0*R - 1.0) - 0.5*C1
    F9 = t - C2

    return (2.0/F3 + 1.0/F4 + (F4*F5 + F6*F7 - F8*F9)/F2/F4)/5.0


def Kd(Km, Gm, a, phi):
    p = phi*T2(Km, Gm, a)
    C1 = 3.0*Km + 4.0*Gm
    return Km*(C1 - 4.0*Gm*p)/(C1 + 3.0*Km*p)


def Gd(Km, Gm, a, phi):
    q = phi*F2(Km, Gm, a)
    C1 = 3.0*Km + 4.0*Gm
    return Gm*(1.0 - C1*q/(C1 + 1.2*(Km + 2.0*Gm)*q))
