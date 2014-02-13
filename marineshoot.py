#!/usr/bin/env python

## full marine ice sheet by shooting

import numpy as np
import matplotlib.pyplot as plt
from sys import exit
import argparse
from scipy.linalg import solve
from scipy.integrate import odeint

import exactN
from plottool import lw, plottwo

desc='''Solve ODE system for marine ice sheet, for which an exact solution
is known.  By default, does shooting to solve for the right value.
If options '--noshoot --saveroot NAME' are given then generates plots.'''

parser = argparse.ArgumentParser(description=desc)

Tmultdefault = 1.0
parser.add_argument('--Tmult', default=Tmultdefault,
                    help='when plotting, multiply exact value of T to use as initial value')
parser.add_argument('--saveroot', metavar='NAME',
                    help='if given, save plots in NAME-geometry.pdf and NAME-other.pdf')
parser.add_argument('--noshoot', action='store_true',
                    help='do not do shooting with bisection and instead plot one case')

args = parser.parse_args()
Tmult = float(args.Tmult)
doshoot = not(bool(args.noshoot))
if not doshoot:
    nameroot = str(args.saveroot)

afloat = 0.0 * exactN.a   # FIXME: what do I want?

def ground(H):  # decide flotation criterion
  return (exactN.rho * H >= exactN.rhow * exactN.bg)

def M(H):       # compute accumulation
  if ground(H):
      return exactN.a * (H - exactN.Hela)
  else:
      return afloat

def F(H):       # compute driving stress coefficient:  tau_d = F(H) H_x
  gover = exactN.rho * exactN.g * H  # grounded overburden pressure
  if ground(H):
      return - gover
  else:
      return - (1.0 - exactN.rho/exactN.rhow) * gover

def beta(x,H):    # compute basal resistance coefficient:  tau_b = beta(x,H) u
  if ground(H) & (x <= exactN.L0):
      Hexact, _, _, _ = exactN.exactN(x)
      gover = exactN.rho * exactN.g * Hexact  # grounded overburden pressure
      return exactN.k * gover
  else:
      return 0.0

def B(x):       # get contrived ice hardness
  if x <= exactN.xg:
      _, B = exactN.exactNbueler(x)
      return B
  else:
      _, Bc = exactN.exactNbueler(exactN.xg)
      return Bc

def surfaces(H):
  hh = H.copy()
  flot = np.logical_not(ground(H))
  rr = exactN.rho / exactN.rhow
  hh[flot] = exactN.bg + (1.0 - rr) * H[flot]
  bb = np.zeros(np.shape(H))
  bb[flot] = exactN.bg - rr * H[flot]
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
                [ beta(x,v[1]) * v[0] ]]);
  vp = solve(A,g)
  return vp.flatten()

rr = exactN.rho / exactN.rhow
xinit = 0.1 * exactN.L0
ux = (1.0/exactN.k) * (2.0 * exactN.H0 / exactN.L0**2)
uinit = ux * xinit
Hinit = exactN.H0 * (1.0 - (xinit / exactN.L0)**2.0)

xc = 1.4 * exactN.L0

def objective(Tinit):
    v0 = np.array([uinit, Hinit, Tinit])  # initial values at x[0] = xinit
    x = np.linspace(xinit,xc,2)
    v = odeint(f,v0,x)                    # solve ODE system to determine v[1,2] = T
    Tcalvc = 0.5 * exactN.rho * exactN.g * (1.0 - rr) * v[1,1]**2
    #print "Tinit = %.6e,  (Tfinal - Tcalvc)/Tcalvc = %.6e" % ( Tinit, (v[-1,2]-Tcalvc)/Tcalvc )
    return (v[1,2] - Tcalvc) / Tcalvc

Tcheat = 2.0 * Hinit * B(xinit) * ux**(1.0/exactN.n)  # the cheat is to use B()

if doshoot:  # run bisection
    Tlow = 0.5 * exactN.rho * exactN.g * (1.0 - rr) * 100.0**2
    Thigh = 0.5 * exactN.rho * exactN.g * (1.0 - rr) * 600.0**2
    olow = objective(Tlow)
    ohigh = objective(Thigh)
    "  [Ta,Tb] = [%.7e,%.7e]  gives  [%2.3e,%2.3e]" % (Tlow,Thigh,olow,ohigh)
    for j in range(20):
        Tmid = (Tlow + Thigh) / 2.0
        omid = objective(Tmid)
        if omid * olow > 0.0:
            olow = omid
            Tlow = Tmid
        else:
            ohigh = omid
            Thigh = Tmid
        print "  [Ta,Tb] = [%.7e,%.7e]  gives  [%2.3e,%2.3e]    (omid = %10.3e)" \
            % (Tlow,Thigh,olow,ohigh,omid)
    print
    print "compare Tcheat = %.7e" % Tcheat
    exit(0)

# proceed to make figure if not doshoot

print "initial conditions:"
print "  at xinit = %.3f km we have initial values" % (xinit/1000.0)
print "  u = %7.2f m a-1, H = %.2f m" % (uinit*exactN.secpera, Hinit)

x = np.linspace(xinit,xc,1001)
dx = x[1] - x[0]

Tinit = Tmult * Tcheat
print "  and  Tinit = %.6e Pa m = %f * Tcheat  where Tcheat = %.6e" % (Tinit, Tmult, Tcheat)

v0 = np.array([uinit, Hinit, Tinit])  # initial values at x[0] = xinit
#v, info = odeint(f,v0,x,full_output=True)          # solve ODE system
v = odeint(f,v0,x)          # solve ODE system
#Tcalvc = 0.5 * exactN.rho * exactN.g * (1.0 - rr) * v[-1,1]**2
#print "Tinit = %.6e,  (Tfinal - Tcalvc)/Tcalvc = %.6e" % ( Tinit, (v[-1,2]-Tcalvc)/Tcalvc )

u = v[:,0]
H = v[:,1]
T = v[:,2]

isgnd = (x <= exactN.xg)
xbod = x[isgnd]
Hbod, _, _, _ = exactN.exactN(xbod)
Herror = abs(H[isgnd] - Hbod).max()
print "results:"
print "  maximum error in H = %.2e m for grounded ice" % Herror

print "calving front results:"
print "  at %.3f km we have values" % (x[-1]/1000.0)
print "  u = %7.2f m a-1, H = %.2f m, T = %.6e Pa m" % (u[-1]*exactN.secpera, H[-1], T[-1])
Tcalvc = 0.5 * exactN.rho * exactN.g * (1.0 - rr) * H[-1]**2
print "COMPARE Tcalvc = %.6e Pa m" % Tcalvc

# geometry and velocity in figure 1 (two-axis)
hh, bb = surfaces(H)
fig = plt.figure(figsize=(6,4))
ax1fig1, ax2fig1 = plottwo(fig,x/1000.0,hh,u * exactN.secpera,
                           "z   (m, solid)","ice velocity   (m a-1, dashed)",
                           y1min=0.0)
ax1fig1.hold(True)
ax1fig1.plot(x/1000.0,bb,'k',lw=2.0)
ax1fig1.plot([xc/1000.0,xc/1000.0],[bb[-1],hh[-1]],'k',linewidth=lw)
# show water height as waves
xl = np.linspace(xc,1.14*xc,101)
ax1fig1.plot(xl/1000.0, exactN.bg - 30.0 * abs(np.sin(xl / 4.0e3)),'k',linewidth=0.7*lw)
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

exit(0)

plt.figure(2)
plt.subplot(4,1,1)
plt.plot(x/1000.0,u * exactN.secpera,'k',lw=2.0)
plt.grid(True)
plt.ylabel("velocity u(x)   [m a-1]")
plt.subplot(4,1,2)
plt.plot(x/1000.0,T,'k',lw=2.0)
plt.grid(True)
plt.ylabel("T(x)   [Pa m]")
BB = np.zeros(np.shape(x))
for j in range(len(x)):
   BB[j] = B(x[j])
plt.subplot(4,1,3)
plt.plot(x/1000.0,BB,'k',lw=2.0)
plt.grid(True)
plt.ylabel("B(x)   [Pa s-(1/n)]")
bbeta = np.zeros(np.shape(H))
for j in range(len(H)):
   bbeta[j] = beta(x[j],H[j])
plt.subplot(4,1,4)
plt.plot(x/1000.0,bbeta,'k',lw=2.0)
plt.grid(True)
plt.ylabel(r"$\beta$(x)   [Pa s m-1]")
plt.xlabel("x   [km]")

plt.show()

