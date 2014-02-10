#!/usr/bin/env python

## full marine ice sheet by shooting

#from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.linalg import solve
from scipy.integrate import odeint

import exactN

afloat = -300.0 * exactN.a

def ground(H):  # decide flotation criterion
  return (exactN.rho * H >= exactN.rhow * exactN.bc)

def M(H):       # compute accumulation
  if ground(H):
      return exactN.a * (H - exactN.Hela)
  else:
      return afloat

def F(H):       # compute driving stress
  gover = exactN.rho * exactN.g * H  # grounded overburden pressure
  if ground(H):
      return - gover
  else:
      return - (1.0 - exactN.rho/exactN.rhow) * gover

def beta(H):    # compute basal resistance
  if ground(H):
      gover = exactN.rho * exactN.g * H  # grounded overburden pressure
      return exactN.k * gover
  else:
      return 0.0

def B(x):
  if x <= exactN.xc:
      _, _, _, _, B = exactN.exactN(x)
      return B
  else:
      _, _, _, _, Bc = exactN.exactN(exactN.xc)
      return Bc

def surfaces(H):
  hh = H
  flot = np.logical_not(ground(H))
  rr = exactN.rho / exactN.rhow
  hh[flot] = exactN.bc + (1.0 - rr) * H[flot]
  bb = np.zeros(np.shape(H))
  bb[flot] = exactN.bc - rr * H[flot]
  return hh, bb

# the ode system is
#
#   [ (2 B H)^n  0         0 ]  [ u' ]  =  [ T^n       ]
#   [ H          u         0 ]  [ H' ]  =  [ M(H)      ]
#   [ 0          F(H)      1 ]  [ T' ]  =  [ beta(H) u ]
#
# where
#   M(H)    = a (H - Hela)              GND
#             af                        FLOT
#   F(H)    = - rho g H                 GND
#             - (1 - rho/rhow) rho g H  FLOT
#   beta(H) = k rho g H                 GND
#             0                         FLOT
# and  GND = (rho H >= rhow b0)  and  FLOT = !GND

def f(v,x):
  """Defines the function for the ODE system, which is in form
  A(v) v' = g(v)  so  f(v,x) = A(v) \ g(v)."""
  # dependent var order:  v = [u, H, T]'
  A = np.array([[ (2.0 * B(x) * v[1])**exactN.n, 0.0,                0.0 ],
                [ v[1],                          v[0],               0.0 ],
                [ 0.0,                           F(v[1]),            1.0 ]  ]);
  g = np.array([[ v[2]**exactN.n    ],
                [ M(v[1])           ],
                [ beta(v[1]) * v[0] ]]);
  vp = solve(A,g)
  return vp.flatten()

xinit = 0.1 * exactN.L0
ux = (1.0/exactN.k) * (2.0 * exactN.H0 / exactN.L0**2)
uinit = ux * xinit
Hinit = exactN.H0 * (1.0 - (xinit / exactN.L0)**2.0)

Tinit = 2.0 * Hinit * B(xinit) * ux**(1.0/exactN.n)  # FIXME: a cheat to avoid shoot

print "initial conditions:"
print "  at xinit = %.3f km we have initial values" % xinit
print "  u = %7.2f m a-1, H = %.2f m, T %.3e UNITS" % (uinit, Hinit, Tinit)

x = np.linspace(xinit,1.2 * exactN.xc,1001)
dx = x[1] - x[0]

v0 = np.array([uinit, Hinit, Tinit])  # initial values at x[0] = xinit
v = odeint(f,v0,x)          # solve ODE system

u = v[:,0]
H = v[:,1]
T = v[:,2]

isgnd = (x <= exactN.xc)
xbod = x[isgnd]
Hbod = exactN.H0 * (1.0 - (xbod / exactN.L0)**2.0)
Herror = abs(H[isgnd] - Hbod).max()
print "results:"
print "  maximum error in H = %.2e m for grounded ice" % Herror

plt.figure(1)
hh, bb = surfaces(H)
plt.plot(x/1000.0,hh,'k',lw=2.0)
plt.plot(x/1000.0,bb,'k',lw=2.0)
#plt.plot(xbod/1000.0,Hbod,'go',lw=3.0)
plt.xlabel("x   [km]")
plt.ylabel("elevation   [m]")

plt.figure(2)
plt.subplot(3,1,1)
plt.plot(x/1000.0,u * exactN.secpera,'k',lw=2.0)
plt.grid(True)
plt.ylabel("velocity u(x)   [m]")
plt.subplot(3,1,2)
plt.plot(x/1000.0,T,'k',lw=2.0)
plt.grid(True)
plt.ylabel("T(x)")
BB = np.zeros(np.shape(x))
for j in range(len(x)):
   BB[j] = B(x[j])
plt.subplot(3,1,3)
plt.plot(x/1000.0,BB,'k',lw=2.0)
plt.grid(True)
plt.ylabel("B(x)")
plt.xlabel("x   [km]")

plt.show()

