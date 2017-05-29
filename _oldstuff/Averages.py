# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 17:46:15 2013

@author: fsantos
"""


def Voigt(*args):
    r = 0.0
    for i in range(0, len(args), 2):
        r += args[i]*args[i+1]
    return r


def Reuss(*args):
    r = 0.0
    for i in range(0, len(args), 2):
        r += args[i+1]/args[i]
    return 1.0/r


def VRH(*args):
    return 0.5*(Voigt(*args) + Reuss(*args))


Arithmetic = Voigt


Harmonic = Reuss


def Geometric(*args):
    r = 1.0
    for i in range(0, len(args), 2):
        r *= args[i]**args[i+1]
    return r
