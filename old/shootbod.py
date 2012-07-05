#!/usr/bin/env python

## FIRST DRAFT will just solve the Bodvardsson (1955) problem by shooting; we
## know the exact solution, the parabola in equation (23)

from pylab import *
from sys import exit

from plhalfar import Params

p = Params()
L0 = p.R0  # = 750.0 m
H0 = p.H0  # = 3600.0 m
rg = p.rg

# parameters in Bodvardsson (1955); see equations (23) and (24)
He = H0 / 1.5   # m;  choose equilibrium-line altitude so H(0) = 1.5 He
a  = 0.001 / p.secpera # s-1;  = (1 mm/a) / m = 0.001 a-1
k  = 9.0 * He / (a * L0 * L0);   #  choose k so that eqn (24) gives our L0

def glen_viscosity(dudx):
  """Viscosity from Glen law;  nu in  tau = 2 nu u_x H."""
  eps0 = 0.0;  # unregularized for diagnosis purposes
  #eps0 = 1.0e-9 / p.secpera;  # = (.001 m/a) / (1000 km); regularizing strain rate
  D2 = dudx * dudx + eps0 * eps0;
  n = 3.0;
  A = 1.0e-16 / p.secpera; # 10-16 Pa-3 a-1
  q = ((1.0/n) - 1.0) / 2.0;
  return A**(-1.0/n) * D2**q;

def f(v,x):
  """Defines the function for the ODE system, which is in form
  M(v) v' = g(v)  so  f(v,x) = M(v) \ g(v)."""
  # dependent var order: v = [H, u, tau]'
  # 'mass' matrix:
  M = array([[ v[1],         v[0],              0.0 ],
             [ -rg*v[0],     0.0,               1.0 ],
             [ 0.0,          0.0,               1.0 ]  ]);
             #[ 0.0,          2.0 * nu0 * v[0],  0.0 ]  ]);
  g = array([[ a * (v[0] - He)      ],
             [ k * rg * v[0] * v[1] ],
             [ v[2]                 ]  ]);
  vp = solve(M,g)
  return vp.flatten()

x = linspace(0.0,L0,26) 
dx = x[1] - x[0]

from scipy.integrate import odeint

print "Bodvardsson theory:"
print "  equilibrium line altitude  H_ela = %7.2f m" % (He)

dudx0 = a / 3.0 # from (u H)' = a (H - H_ela) at zero
print "  initial strain rate dudx0 = %.5e s-1" % dudx0

nu0 = glen_viscosity(dudx0);  # Pa s; constant viscosity from Glen's law for
                              #   from exact strain rate from Bodvardsson soln
print "  initial viscosity (if Glen law used) = %.5e Pa s" % nu0

tau0 = 2.0 * nu0 * dudx0 * H0
print "  initial vertically-integrated longitudinal stress tau(0) = %.5e" % tau0

v0 = array([H0, 0.0, tau0])  # initial values at x[0] = 0.0
v = odeint(f,v0,x)          # solve ODE system

H = v[:,0]
u = v[:,1]
tau = v[:,2]

print "results:"
print "  velocity at margin is     u(L0)  = %7.2f m a-1" % (u[-1] * p.secpera)
print "  slope    at margin is     h'(L0) = %7.6f" % ((H[-1] - H[-2]) / dx)

H_bod = H0 * (1.0 - (x / L0)**2.0)
Herror = abs(H - H_bod).max()
print "  maximum error in H = %.2e m" % Herror

figure(1)
subplot(3,1,1)
plot(x/1000.0,H, x/1000.0,He * ones(shape(x)),'g--')
xlabel("x   [km]"), ylabel("elevation H(x)   [m]")
subplot(3,1,2)
plot(x/1000.0,a * (H - He) * p.secpera)
xlabel("x   [km]"), ylabel("accumulation M(x)   [m a-1]")
grid(True)
subplot(3,1,3)
plot(x/1000.0,u * p.secpera)
xlabel("x   [km]"), ylabel("velocity u(x)   [m a-1]")

dudx = (u[1:] - u[:-1]) / dx
print "  initial strain rate  dudx0 = %.5e s-1" % dudx[0]
nu = glen_viscosity(dudx)
print "  viscosity range, from Glen law: [%.5e, %.5e] Pa s" \
      % (nu.min(), nu.max())

show()

