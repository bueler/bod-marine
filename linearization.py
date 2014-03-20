#!/usr/bin/env python

from sys import exit
import numpy as np
import scipy.linalg as la

from exactsolns import H0, L0, rho, g, omega, Bg, ug, Hg, Mg, xg, n,
                       exactBod, exactBodBueler, exactVeen
# H, u, M = exactBod(x)
# T, B = exactBodBueler(x)
# H, u = exactVeen(x)

def ddx_exactBod(x):
    hxx = - 2.0 * H0 / (L0 * L0)
    dHdx  = hxx * x
    dudx  = - hxx / k
    return dHdx, dudx

def ddx_exactVeen(x):
    Cs = (rho * g * omega / (4.0 * Bg))**n
    Qg = ug * Hg
    tmp = Qg + Mg * (x - xg)
    u = ( ug**(n+1) + (Cs / Mg) * ( tmp**(n+1) - Qg**(n+1) ) )**(1.0/(n+1))
    if np.any(u <= 0):
        print "ERROR:  non-positive u detected in ddx_exactVeen()" 
        exit(1)
    dudx = C * tmp**n / u**n
    dHdx = (Mg * u - tmp * dudx) / u**2
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

