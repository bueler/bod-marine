#!/usr/bin/env python

## the Bodvardsson (1955) problem rethought, as an exact SSA solution
## for a marine ice sheet

from sys import exit
import numpy as np

## the following code defining the Test N form of the Bodvardsson solution
## represents complete code duplication from
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

# Bodvardsson (1955) solution is purely grounded

def exactBod(x):
  # Bodvardsson: get geometry and velocity from mass continuity and
  #   T_x = 0  and assumption  beta = k rho g H
  if np.any((x < 0.0) | (x > L0)):
    print "ERROR in exactN(): x out of range [0,L0]"
    return 1
  hxx = - 2.0 * H0 / (L0 * L0)
  q   = (1.0 / n) - 1.0
  ux  = - hxx / k
  H   = H0 * (1.0 - (x / L0) * (x / L0))
  hx  = hxx * x
  u   = - hx / k
  M   = a * (H - Hela)
  return H, hx, u, M

# now build marine ice sheet by deciding where the grounding line is
# and setting the ocean depth bg accordingly

rhow    = 1028.0
omega   = 1.0 - rho / rhow
xg      = 0.9 * L0

Hg, _, ug, Mg = exactBod(xg)
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
  H, _, _, _ = exactBod(x)
  T = 0.5 * omega * rho * g * Hg**2
  B = T / ( 2.0 * H * (abs(ux)**q) * ux )
  return T, B

_, Bg   = exactBodBueler(xg)

def exactVeen(x,M0):
  # van der Veen: get thickness and velocity for floating ice shelf, given
  #   Mg, Bg, xg, Hg, ug computed above from Bodvardsson
  C = (rho * g * omega / (4.0 * Bg))**n
  tmp = ug * Hg + M0 * (x - xg)
  u = ( ug**(n+1) + (C / M0) * ( tmp**(n+1) - (ug * Hg)**(n+1) ) )**(1.0/(n+1))
  if np.any(u <= 0):
    print "ERROR:  non-positive u detected in exactVeen()" 
    exit(1)
  H = tmp / u
  return H, u

