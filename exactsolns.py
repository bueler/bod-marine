#!/usr/bin/env python

## the Bodvardsson (1955) problem rethought, as an exact SSA solution
## for a marine ice sheet

from sys import exit
import numpy as np

## the following code came from the PISM Test N form of the Bodvardsson solution
##   pism-dev/src/verif/tests/exactTestN.[h|c]
## see comments there

secpera = 31556926.0
g       = 9.81
rho     = 910.0
n       = 3.0

H0      = 3000.0
L0      = 500.0e3

a       = 0.003 / secpera

Hela    = H0 / 1.5
k       = 9.0 * Hela / (a * L0 * L0)


def exactBod(x):
  # Bodvardsson: get geometry and velocity from mass continuity and
  #   T_x = 0  and assumption  beta = k rho g H;  grounded only
  if np.any((x < 0.0) | (x > L0)):
    print "ERROR in exactBod(): x out of range [0,L0]"
    return 1
  hxx = - 2.0 * H0 / (L0 * L0)
  q   = (1.0 / n) - 1.0
  ux  = - hxx / k
  H   = H0 * (1.0 - (x / L0) * (x / L0))
  hx  = hxx * x
  u   = - hx / k
  M   = a * (H - Hela)
  return H, u, M

# now build marine ice sheet by deciding where the grounding line is
# and setting the ocean depth bg accordingly

rhow    = 1028.0
omega   = 1.0 - rho / rhow
xg      = 0.9 * L0

Hg, ug, Mg = exactBod(xg)
bedg    = rho * Hg / rhow

def exactBodBueler(x):
  # Bueler: additionally ...
  #   (1) set constant vertically-integrated longitudinal stress T
  #       to the value it would have at marine calving front at flotation
  #   (2) determine ice hardness B from T so that Bodvardsson's
  #       plug flow is actually grounded SSA
  hxx = - 2.0 * H0 / (L0 * L0)
  q   = (1.0 / n) - 1.0
  ux  = - hxx / k
  H, _, _ = exactBod(x)
  T = 0.5 * omega * rho * g * Hg**2
  B = T / ( 2.0 * H * (abs(ux)**q) * ux )
  return T, B

_, Bg   = exactBodBueler(xg)

def exactVeen(x):
  # van der Veen: get thickness and velocity for floating ice shelf, given
  #   Mg, Bg, xg, Hg, ug computed above from Bodvardsson
  Cs = (rho * g * omega / (4.0 * Bg))**n
  Qg = ug * Hg
  tmp = Qg + Mg * (x - xg)
  u = ( ug**(n+1) + (Cs / Mg) * ( tmp**(n+1) - Qg**(n+1) ) )**(1.0/(n+1))
  if np.any(u <= 0):
    print "ERROR:  non-positive u detected in exactVeen()" 
    exit(1)
  H = tmp / u
  return H, u


xa = 0.2 * L0
xc = 0.98 * L0
Ha, ua, _ = exactBod(xa)
Tg = 0.5 * rho * g * omega * Hg**2
Hc, uc = exactVeen(xc)
Tc = 0.5 * rho * g * omega * Hc**2

def printtable():
  print 'important constants and values for exact solution:'
  print '  H0    = %10.3f m' % H0
  print '  L0    = %10.3f km' % (L0 / 1000.0)
  print '  zo    = %10.3f m' % bedg
  print '  [xa   = %10.3f km]' % (xa / 1000.0)
  print '  H(0)  = %10.3f m' % Ha
  print '  u(0)  = %10.3f m/a' % (ua * secpera)
  print '  xg    = %10.3f km   (from start xa)' % ((xg - xa) / 1000.0)
  print '  H(xg) = %10.3f m' % Hg
  print '  u(xg) = %10.3f m/a' % (ug * secpera)
  print '  T(xg) = %10.3f 10^8 Pa m' % (1.0e-8 * Tg)
  print '  B(xg) = %10.3f 10^8 Pa s^(1/3)' % (1.0e-8 * Bg)
  print '  M(xg) = %10.3f m/a' % (Mg * secpera)
  print '  xc    = %10.3f km   (from start xa)' % ((xc - xa) / 1000.0)
  print '  H(xc) = %10.3f m' % Hc
  print '  u(xc) = %10.3f m/a' % (uc * secpera)
  print '  T(xc) = %10.3f 10^8 Pa m' % (1.0e-8 * Tc)

