#!/usr/bin/env python

from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import exactsolns
from exactsolns import H0, L0, rho, g, omega, xa, Bg, ug, Hg, Mg, xg, n, k, \
                       exactBod, exactBodBueler, exactVeen

def ddx_exactBod(x):
    hxx = - 2.0 * H0 / (L0 * L0)
    dHdx  = hxx * x
    dudx  = - hxx / k
    return dHdx, dudx

def ddx_exactVeen(x):
    _, u = exactVeen(x)
    Qg = ug * Hg
    tmp = Qg + Mg * (x - xg)
    Cs = (rho * g * omega / (4.0 * Bg))**n
    dudx = Cs * tmp**n / u**n
    dHdx = (Mg * u - tmp * dudx) / u**2
    return dHdx, dudx

# linear ODE system for w = (tilu(x),tilH(x),tilT(x)).T is of form
#    LL(x) w_x = RR(x) w
# where LL, RR are 3x3 matrices, equivalently
#    v_x = AA(x) v
# where AA(x) = LL(x) \ RR(x)

q = (1.0/n) - 1.0
rg = rho * g

def getsides(x):
    if x <= xg:
       H, u, _ = exactBod(x)
       T, B = exactBodBueler(x)
       dH, du = ddx_exactBod(x)
       LL = np.array([[(2.0/n)*B*H*du**q, 0,     0],
                      [H,                 u,     0],
                      [0,                 -rg*H, 1]])
       RR = np.array([[0,       -2.0*B*du**(1.0/n), 1],
                      [-dH,     -du,                0],
                      [-k*rg*H, -k*rg*u+rg*dH,      0]])
    else:
       H, u = exactVeen(x)
       dH, du = ddx_exactVeen(x)
       LL = np.array([[(2.0/n)*Bg*H*du**q, 0,           0],
                      [H,                  u,           0],
                      [0,                  -omega*rg*H, 1]])
       RR = np.array([[0,   -2.0*Bg*du**(1.0/n), 1],
                      [-dH, -du,                 0],
                      [0,   omega*rg*dH,         0]])
    return LL, RR

# to understand stiffness we find eigs of AA(x) at each x
x = np.linspace(exactsolns.xa,exactsolns.xc,1001)
lam = np.zeros((len(x),3),'complex')
for j in range(len(x)):
    LL, RR = getsides(x[j])
    AA = la.solve(LL,RR)
    lam[j,:], X = la.eig(AA)
lam = 1000.0 * lam   # convert eigenvalues to km-1

# sort complex eigs so first is real (and neg) and last two are complex (w.
# small pos real part)
indx = np.argsort(np.real(lam),axis=1)
for j in range(len(x)):
    lam[j,:] = lam[j,indx[j]]

#print lam
#exit(0)

# show real parts of eigenvalues
fig = plt.figure(figsize=(6,4))
plt.subplot(2,1,1)
plt.semilogy((x-xa)/1000.0,-np.real(lam[:,0]),'k',lw=2.0)
plt.grid(True)
plt.ylabel(r"$-\lambda_1$")
plt.subplot(2,1,2)
plt.semilogy((x-xa)/1000.0,np.real(lam[:,1:2]),'k',lw=2.0)
plt.grid(True)
plt.ylabel(r"$\Re(\lambda_2)=\Re(\lambda_3)$")
plt.xlabel("x   (km)")

# absolute real parts
arl = np.zeros((len(x),3))
for j in range(len(x)):
    arl[j,:] = np.sort(abs(np.real(lam[j,:])))   # sort in increasing magnitude

# show "stiffness ratio" of largest to smallest real
fig = plt.figure(figsize=(6,4))
#plt.semilogy((x-xa)/1000.0,arl,'k--',lw=2.0)
plt.semilogy((x-xa)/1000.0,arl[:,-1]/arl[:,0],'k',lw=2.0)
plt.grid(True)
plt.xlabel("x   (km)")
nameroot = "exactmarine"
imagename = nameroot + '-stiffness-ratio.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()

