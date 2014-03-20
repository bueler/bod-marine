#!/usr/bin/env python

from sys import exit
import numpy as np
import scipy.linalg as la

from exactsolns import *
# n, xg
# H, u, M = exactBod(x)
# T, B = exactBodBueler(x)
# H, u = exactVeen(x,M0)

def ddx_exactBod(x):
    hxx = - 2.0 * H0 / (L0 * L0)
    dHdx  = hxx * x
    dudx  = - hxx / k
    return dHdx, dudx

def ddx_exactVeen(x,M0):
    C = (rho * g * omega / (4.0 * Bg))**n
    tmp = ug * Hg + M0 * (x - xg)
    uspow = ug**(n+1) + (C / M0) * ( tmp**(n+1) - (ug * Hg)**(n+1) )
    if np.any(u <= 0):
        print "ERROR:  non-positive u detected in ddx_exactVeen()" 
        exit(1)
    H = tmp / u
    return H, u
    FIXME
    return dHdx, dudx

# linear ODE system for v = (u1(x),H1(x),T1(x)).T is of form
#    LL(v) v_x = RR(v) v
# where LL, RR are 3x3 matrices, equivalently
#    v_x = AA(v) v
# to understand stiffness we find eigs of AA at some x

for x in [xg-10.0e3, xg+10.0e3]:
    #LL = np.array([[ FIXME ]])
    #RR = np.array([[ FIXME ]])
    AA = linalg.solve(L,R)
    lam = linalg.eig(A)
    print lam

