#!/usr/bin/env python

## full marine ice sheet by shooting

import numpy as np
import matplotlib.pyplot as plt
from sys import exit
import argparse
from scipy.linalg import solve
from scipy.integrate import odeint

import exactsolns
from plottwo import lw, plottwo

desc='''Solve ODE system for marine ice sheet, for which an exact solution
is known.  By default, does shooting to solve for the right value.
If options '--noshoot' are given then generates plots.'''

# process command-line arguments
parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
xocean_default = 1.05
parser.add_argument('--xocean', default=xocean_default,
                    help='end of ocean in plots, , as fraction of L0 = %.1f km' % (exactsolns.L0/1000.0))
parser.add_argument('--saveroot', metavar='NAME', default='unnamed',
                    help='if given, save plots in NAME-geometry.pdf and NAME-other.pdf')
parser.add_argument('--figures', action='store_true',
                    help='generate figures in addition to shooting')
args = parser.parse_args()
genfigs = bool(args.figures)
nameroot = str(args.saveroot)

# only relevant to plots
xoceanend = float(args.xocean) * exactsolns.L0

# define the model

def ground(H):  # decide flotation criterion given thickness
  return (exactsolns.rho * H >= exactsolns.rhow * exactsolns.bedg)

def M(x):       # compute mass balance for location
  if x <= exactsolns.xg:
      _, _, Mexact = exactsolns.exactBod(x)
      return Mexact
  else:
      return exactsolns.Mg

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
      Hexact, _, _ = exactsolns.exactBod(x)
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

rr = exactsolns.rho / exactsolns.rhow

def objective(Tinit):
    v0 = np.array([exactsolns.ua, exactsolns.Ha, Tinit])  # initial values at x[0] = xa
    x = np.linspace(exactsolns.xa,exactsolns.xc,2)
    v = odeint(f,v0,x)              # solve ODE system to determine v[1,2] = T
    Tcalvc = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * v[1,1]**2
    return (v[1,2] - Tcalvc) / Tcalvc

uxcheat = (1.0/exactsolns.k) * (2.0 * exactsolns.H0 / exactsolns.L0**2)
Tcheat = 2.0 * exactsolns.Ha * B(exactsolns.xa) * uxcheat**(1.0/exactsolns.n)  # the cheat is to use B()

# run bisection
Tlow = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * 100.0**2
Thigh = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * 600.0**2
olow = objective(Tlow)
ohigh = objective(Thigh)
"  [Tlow,Thigh] = [%.7e,%.7e]  gives  [%2.3e,%2.3e]" % (Tlow,Thigh,olow,ohigh)
for j in range(25):
    Tmid = (Tlow + Thigh) / 2.0
    omid = objective(Tmid)
    if omid * olow > 0.0:
        olow = omid
        Tlow = Tmid
    else:
        ohigh = omid
        Thigh = Tmid
    print "  T range [%.7e,%.7e]  gives  objective range [%2.3e,%2.3e]" \
        % (Tlow,Thigh,olow,ohigh)
Tmid = (Tlow + Thigh) / 2.0
omid = objective(Tmid)
print "final:   Tmid   = %.7e  gives  omid = %2.3e" % (Tmid,omid)
print "[compare Tcheat = %.7e]" % Tcheat

if not(genfigs):
    exit(0)

# proceed to make figures

print "initial conditions:"
print "  at xa = %.3f km we have initial values" % (exactsolns.xa/1000.0)
print "  u = %7.2f m a-1, H = %.2f m" % (exactsolns.ua*exactsolns.secpera, exactsolns.Ha)
print "  and  Tinit = %.6e Pa m (which is cheating, for plot purpose)" % Tcheat

# solve ODES: first with cheating initial, then with result of bisection
# get result on reasonably fine grid
x = np.linspace(exactsolns.xa,exactsolns.xc,1001)
v0 = np.array([exactsolns.ua, exactsolns.Ha, Tcheat])  # initial values at x[0] = xa
#v, info = odeint(f,v0,x,full_output=True)          # solve ODE system
v = odeint(f,v0,x)          # solve ODE system
v0mid = np.array([exactsolns.ua, exactsolns.Ha, Tmid])
vmid = odeint(f,v0mid,x)              # solve ODE system to determine v[1,2] = T

u = v[:,0]
H = v[:,1]
T = v[:,2]
umid = vmid[:,0]
Hmid = vmid[:,1]
Tmid = vmid[:,2]

# compute and report errors
xgnd = x[(x <= exactsolns.xg)]
Hbod, ubod, _ = exactsolns.exactBod(xgnd)
Hboderror = abs(H[(x <= exactsolns.xg)] - Hbod).max()
uboderror = abs(u[(x <= exactsolns.xg)] - ubod).max()
xveen = x[(x > exactsolns.xg)]
Hveen, uveen = exactsolns.exactVeen(xveen,exactsolns.Mg)
Hveenerror = abs(H[(x > exactsolns.xg)] - Hveen).max()
uveenerror = abs(u[(x > exactsolns.xg)] - uveen).max()
print "results:"
print "  max |H-Hexact| = %.2e m      for grounded ice" % Hboderror
print "  max |u-uexact| = %.2e m a-1  for grounded ice" % (uboderror * exactsolns.secpera)
print "  max |H-Hexact| = %.2e m      for floating ice" % Hveenerror
print "  max |u-uexact| = %.2e m a-1  for floating ice" % (uveenerror * exactsolns.secpera)

print "calving front results:"
print "  at %.3f km we have values" % (x[-1]/1000.0)
print "  u = %7.2f m a-1, H = %.2f m, T = %.6e Pa m" % (u[-1]*exactsolns.secpera, H[-1], T[-1])
Tcalvc = 0.5 * exactsolns.rho * exactsolns.g * (1.0 - rr) * H[-1]**2
print "COMPARE Tcalvc = %.6e Pa m" % Tcalvc

# geometry and velocity in figure 1 (two-axis)
hh, bb = surfaces(H)
fig = plt.figure(figsize=(6,4))
ax1fig1, ax2fig1 = plottwo(fig,(x-exactsolns.xa)/1000.0,hh,u * exactsolns.secpera,
                           "z   (m, solid)",
                           "ice velocity   (m a-1, dashed)",
                           y1min=0.0)
ax1fig1.set_xlim(0.0,(xoceanend-exactsolns.xa)/1000.0)
ax2fig1.set_xlim(0.0,(xoceanend-exactsolns.xa)/1000.0)
ax1fig1.hold(True)
ax1fig1.plot((x-exactsolns.xa)/1000.0,bb,'k',lw=2.0)
ax1fig1.plot([(exactsolns.xc-exactsolns.xa)/1000.0,(exactsolns.xc-exactsolns.xa)/1000.0],[bb[-1],hh[-1]],'k',linewidth=lw)
# show water height as waves
xl = np.linspace(exactsolns.xc,xoceanend,101)
ax1fig1.plot((xl-exactsolns.xa)/1000.0, exactsolns.bedg - 30.0 * abs(np.sin(xl / 2.0e3)),'k',linewidth=0.7*lw)
y1, y2 = ax1fig1.get_ylim()
ax1fig1.set_ylim(0.0,3000.0)
ax1fig1.set_yticks([0.0,500.0,1000.0,1500.0,2000.0,2500.0,3000.0])
ax1fig1.hold(False)
y1, y2 = ax2fig1.get_ylim()
ax2fig1.set_ylim(0.0,1.05*y2)
imagename = nameroot + '-geometry.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()
#exit(0)

# geometry and velocity DETAIL in figure 2 (two-axis)
fig = plt.figure(figsize=(6,4))
detail = (x > 0.95*exactsolns.xg)
ax1fig1, ax2fig1 = plottwo(fig,(x[detail]-exactsolns.xa)/1000.0,hh[detail],u[detail] * exactsolns.secpera,
                           "z   (m, solid)",
                           "ice velocity   (m a-1, dashed)",
                           y1min=0.0)
ax1fig1.hold(True)
ax1fig1.plot((x[detail]-exactsolns.xa)/1000.0,bb[detail],'k',lw=2.0)
ax1fig1.plot([(exactsolns.xc-exactsolns.xa)/1000.0,(exactsolns.xc-exactsolns.xa)/1000.0],[bb[-1],hh[-1]],'k',linewidth=lw)
# show water height as waves
xl = np.linspace(exactsolns.xc,1.01*exactsolns.L0,101)
ax1fig1.plot((xl-exactsolns.xa)/1000.0, exactsolns.bedg - 15.0 * abs(np.sin(xl / 0.8e3)),'k',linewidth=0.7*lw)
ax1fig1.hold(False)
plt.axis('tight')
imagename = nameroot + '-geometry-detail.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()
#exit(0)

# mass balance and ice hardness in figure 3 (two-axis)
MM = np.zeros(np.shape(x))
for j in range(len(x)):
   MM[j] = M(x[j])
BB = np.zeros(np.shape(x))
for j in range(len(x)):
   BB[j] = B(x[j])
fig = plt.figure(figsize=(6,4))
ax1fig2, ax2fig2 = plottwo(fig,(x-exactsolns.xa)/1000.0,MM * exactsolns.secpera,BB/1.0e8,
                           "M(x)   (m a-1, solid)",
                           r"B(x)   ($10^8$ Pa s-1/3, dashed)")
ax1fig2.set_xlim(0.0,(xoceanend-exactsolns.xa)/1000.0)
ax2fig2.set_xlim(0.0,(xoceanend-exactsolns.xa)/1000.0)
ax1fig2.set_ylim(-5.0,5.0)
ax2fig2.set_ylim(0.0,5.0)
imagename = nameroot + '-M-B.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()

# sliding resistance beta in figure 4
fig = plt.figure(figsize=(6,4))
bbeta = np.zeros(np.shape(H))
for j in range(len(H)):
   if x[j] <= exactsolns.L0:
       Hexact, _, _ = exactsolns.exactBod(x[j])
       bbeta[j] = exactsolns.k * exactsolns.rho * exactsolns.g * Hexact
   else:
       bbeta[j] = 0.0
ax1fig3, ax2fig3 = plottwo(fig,(x-exactsolns.xa)/1000.0,bbeta/1.0e10,T/1.0e8,
                           r"$\beta$(x)   ($10^{10}$  Pa s m-1, solid)",
                           r"$T$(x)   ($10^8$ Pa m, dashed)")
ax1fig3.set_ylim(-0.1,2.1)
ax2fig3.set_ylim(0.0,2.0)
bbeta[x > exactsolns.xg] = 0.0
ax1fig3.plot((x-exactsolns.xa)/1000.0,bbeta/1.0e10,'k:',lw=1.5)
imagename = nameroot + '-beta-T.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()

# thickness error and velocity error in figure 5
Hexact = np.zeros(np.shape(H))
uexact = np.zeros(np.shape(u))
for j in range(len(x)):
  if x[j] < exactsolns.xg:
    Hexact[j], uexact[j], _ = exactsolns.exactBod(x[j])
  else:
    Hexact[j], uexact[j] = exactsolns.exactVeen(x[j],exactsolns.Mg)
fig = plt.figure(figsize=(6,5))
plt.subplot(2,1,1)
plt.semilogy((x-exactsolns.xa)/1000.0,abs(H-Hexact),'k',lw=2.0)
plt.semilogy((x-exactsolns.xa)/1000.0,abs(Hmid-Hexact),'k--',lw=3.0)
ax1 = plt.gca()
ax1.set_xlim([0.0,(xoceanend-exactsolns.xa)/1000.0])
plt.grid(True)
plt.ylabel("H error   (m)", fontsize=14)
plt.subplot(2,1,2)
plt.semilogy((x-exactsolns.xa)/1000.0,abs(u-uexact) * exactsolns.secpera,'k',lw=2.0)
plt.semilogy((x-exactsolns.xa)/1000.0,abs(umid-uexact) * exactsolns.secpera,'k--',lw=3.0)
ax2 = plt.gca()
ax2.set_xlim([0.0,(xoceanend-exactsolns.xa)/1000.0])
plt.grid(True)
plt.ylabel("u error   (m a-1)", fontsize=14)
plt.xlabel("x   (km)")
imagename = nameroot + '-error.pdf'
plt.savefig(imagename)
print '  image file %s saved' % imagename
#plt.show()

