#!/usr/bin/env python

## full marine ice sheet by shooting

import numpy as np
import matplotlib.pyplot as plt
from sys import exit
import argparse
from scipy.linalg import solve
from scipy.integrate import odeint

import exactsolns
from plottool import lw, plottwo

desc='''Solve ODE system for marine ice sheet, for which an exact solution
is known.  By default, does shooting to solve for the right value.
If options '--noshoot' are given then generates plots.'''

# process command-line arguments
parser = argparse.ArgumentParser(description=desc)
Tmult_default = 1.0
Mdrop_default = 0.1
parser.add_argument('--Tmult', default=Tmult_default,
                    help='when plotting, multiply exact value of T to use as initial value')
parser.add_argument('--Mdrop', default=Mdrop_default,
                    help='factor by which mass balance M(x) drops in magnitude across grounding line')
parser.add_argument('--saveroot', metavar='NAME', default='unnamed',
                    help='if given, save plots in NAME-geometry.pdf and NAME-other.pdf')
parser.add_argument('--noshoot', action='store_true',
                    help='do not do shooting with bisection and instead plot one case')
args = parser.parse_args()
Tmult = float(args.Tmult)
Mdrop = float(args.Mdrop)
doshoot = not(bool(args.noshoot))
if not doshoot:
    nameroot = str(args.saveroot)

# define the model

def ground(H):  # decide flotation criterion given thickness
  return (exactsolns.rho * H >= exactsolns.rhow * exactsolns.bedg)

def M(x):       # compute mass balance for location
  if x <= exactsolns.xg:
      _, _, _, Mexact = exactsolns.exactBod(x)
      return Mexact
  else:
      return Mdrop * exactsolns.Mg

def B(x):       # get contrived ice hardness for location
  if x <= exactsolns.xg:
      _, Bexact = exactsolns.exactBodBueler(x)
      return Bexact
  else:
      return exactsolns.Bg

def F(H):       # compute driving stress coefficient:    tau_d = - F(H) H_x
  gover = exactsolns.rho * exactsolns.g * H  # grounded overburden pressure
  if ground(H):
      return - gover
  else:
      return - (1.0 - exactsolns.rho/exactsolns.rhow) * gover

def beta(x,H):  # compute basal resistance coefficient:  tau_b = - beta(x,H) u
  if ground(H) & (x <= exactsolns.L0):
      Hexact, _, _, _ = exactsolns.exactBod(x)
      gover = exactsolns.rho * exactsolns.g * Hexact  # grounded overburden pressure
      return exactsolns.k * gover
  else:
      return 0.0

def surfaces(H):
  hh = H.copy()
  flot = np.logical_not(ground(H))
  rr = exactsolns.rho / exactsolns.rhow
  hh[flot] = exactsolns.bedg + (1.0 - rr) * H[flot]
  bb = np.zeros(np.shape(H))
  bb[flot] = exactsolns.bedg - rr * H[flot]
  return hh, bb

# the ode system is
#
#   [ (2 B(x) H)^n  0         0 ]  [ u' ]  =  [ T^n         ]
#   [ H             u         0 ]  [ H' ]  =  [ M(x)        ]
#   [ 0             F(H)      1 ]  [ T' ]  =  [ beta(x,H) u ]
#
# where
#   M(x)      = a (Hexact(x) - Hela)      x <= xg
#               Mg                        x >  xg
#   B(x)      = Bexact(x)                 x <= xg
#               Bg                        x >  xg
#   F(H)      = - rho g H                 GND
#               - omega rho g H           FLOT
#   beta(x,H) = k rho g Hexact(x)         GND
#               0                         FLOT
# and  GND = (rho H >= rhow b0)  and  FLOT = !GND

def f(v,x):
  """Defines the function for the ODE system, which is in form
  A(v) v' = g(v)  so  f(v,x) = A(v) \ g(v)."""
  # dependent var order:  v = [u, H, T]'
  A = np.array([[ (2.0 * B(x) * v[1])**exactsolns.n, 0.0,                0.0 ],
                [ v[1],                          v[0],               0.0 ],
                [ 0.0,                           F(v[1]),            1.0 ]  ]);
  g = np.array([[ v[2]**exactsolns.n    ],
                [ M(x)              ],
                [ beta(x,v[1]) * v[0] ]]);
  vp = solve(A,g)
  return vp.flatten()

xa = 0.1 * exactsolns.L0
xc = 1.4 * exactsolns.L0

rr = exactsolns.rho / exactsolns.rhow
Ha, _, ua, _ = exactsolns.exactBod(xa)

def objective(Tinit):
    v0 = np.array([ua, Ha, Tinit])  # initial values at x[0] = xa
    x = np.linspace(xa,xc,2)
    v = odeint(f,v0,x)              # solve ODE system to determine v[1,2] = T
    Tcalvc = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * v[1,1]**2
    return (v[1,2] - Tcalvc) / Tcalvc

uxcheat = (1.0/exactsolns.k) * (2.0 * exactsolns.H0 / exactsolns.L0**2)
Tcheat = 2.0 * Ha * B(xa) * uxcheat**(1.0/exactsolns.n)  # the cheat is to use B()

if doshoot:  # run bisection
    Tlow = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * 100.0**2
    Thigh = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * 600.0**2
    olow = objective(Tlow)
    ohigh = objective(Thigh)
    "  [Tlow,Thigh] = [%.7e,%.7e]  gives  [%2.3e,%2.3e]" % (Tlow,Thigh,olow,ohigh)
    for j in range(20):
        Tmid = (Tlow + Thigh) / 2.0
        omid = objective(Tmid)
        if omid * olow > 0.0:
            olow = omid
            Tlow = Tmid
        else:
            ohigh = omid
            Thigh = Tmid
        print "  [Tlow,Thigh] = [%.7e,%.7e]  gives  [%2.3e,%2.3e]    (omid = %10.3e)" \
            % (Tlow,Thigh,olow,ohigh,omid)
    print
    print "compare Tcheat = %.7e" % Tcheat
    exit(0)

# proceed to make figure if not doshoot

print "initial conditions:"
print "  at xa = %.3f km we have initial values" % (xa/1000.0)
print "  u = %7.2f m a-1, H = %.2f m" % (ua*exactsolns.secpera, Ha)

x = np.linspace(xa,xc,1001)
dx = x[1] - x[0]

Tinit = Tmult * Tcheat
print "  and  Tinit = %.6e Pa m = %f * Tcheat  where Tcheat = %.6e" % (Tinit, Tmult, Tcheat)

v0 = np.array([ua, Ha, Tinit])  # initial values at x[0] = xa
#v, info = odeint(f,v0,x,full_output=True)          # solve ODE system
v = odeint(f,v0,x)          # solve ODE system
#Tcalvc = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * v[-1,1]**2
#print "Tinit = %.6e,  (Tfinal - Tcalvc)/Tcalvc = %.6e" % ( Tinit, (v[-1,2]-Tcalvc)/Tcalvc )

u = v[:,0]
H = v[:,1]
T = v[:,2]

isgnd = (x <= exactsolns.xg)
xbod = x[isgnd]
Hbod, _, _, _ = exactsolns.exactBod(xbod)
Herror = abs(H[isgnd] - Hbod).max()
print "results:"
print "  maximum error in H = %.2e m for grounded ice" % Herror

print "calving front results:"
print "  at %.3f km we have values" % (x[-1]/1000.0)
print "  u = %7.2f m a-1, H = %.2f m, T = %.6e Pa m" % (u[-1]*exactsolns.secpera, H[-1], T[-1])
Tcalvc = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * H[-1]**2
print "COMPARE Tcalvc = %.6e Pa m" % Tcalvc

# geometry and velocity in figure 1 (two-axis)
hh, bb = surfaces(H)
fig = plt.figure(figsize=(6,4))
ax1fig1, ax2fig1 = plottwo(fig,x/1000.0,hh,u * exactsolns.secpera,
                           "z   (m, solid)",
                           "ice velocity   (m a-1, dashed)",
                           y1min=0.0)
ax1fig1.hold(True)
ax1fig1.plot(x/1000.0,bb,'k',lw=2.0)
ax1fig1.plot([xc/1000.0,xc/1000.0],[bb[-1],hh[-1]],'k',linewidth=lw)
# show water height as waves
xl = np.linspace(xc,1.14*xc,101)
ax1fig1.plot(xl/1000.0, exactsolns.bedg - 30.0 * abs(np.sin(xl / 4.0e3)),'k',linewidth=0.7*lw)
y1, y2 = ax1fig1.get_ylim()
ax1fig1.set_ylim(0.0,3000.0)
ax1fig1.set_yticks([0.0,500.0,1000.0,1500.0,2000.0,2500.0,3000.0])
ax1fig1.hold(False)
y1, y2 = ax2fig1.get_ylim()
ax2fig1.set_ylim(0.0,y2)
imagename = nameroot + '-geometry.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()

#exit(0)

# mass balance and ice hardness in figure 2 (two-axis)
MM = np.zeros(np.shape(x))
for j in range(len(x)):
   MM[j] = M(x[j])
BB = np.zeros(np.shape(x))
for j in range(len(x)):
   BB[j] = B(x[j])
fig = plt.figure(figsize=(6,4))
ax1fig2, ax2fig2 = plottwo(fig,x/1000.0,MM * exactsolns.secpera,BB/1.0e8,
                           "M(x)   (m a-1, solid)",
                           r"B(x)   ($10^8$ Pa s-1/3, dashed)")
ax1fig2.set_xlim(0.0,800.0)
ax2fig2.set_xlim(0.0,800.0)
ax1fig2.set_ylim(-5.0,5.0)
ax2fig2.set_ylim(0.0,5.0)
imagename = nameroot + '-M-B.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()

#exit(0)

fig = plt.figure(figsize=(6,4))
bbeta = np.zeros(np.shape(H))
for j in range(len(H)):
   if x[j] <= exactsolns.L0:
       Hexact, _, _, _ = exactsolns.exactBod(x[j])
       bbeta[j] = exactsolns.k * exactsolns.rho * exactsolns.g * Hexact
   else:
       bbeta[j] = 0.0
plt.plot(x/1000.0,bbeta/1.0e10,'k',lw=2.0)
plt.axis([0.0,800.0,-0.1,2.1])
plt.ylabel(r"$\beta_0$(x)   ($10^{10}$  Pa s m-1)")
plt.xlabel("x   (km)")
imagename = nameroot + '-beta.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()

