#!/usr/bin/env python

## the Bodvardsson (1955) problem thought of as a
## verification test for a marine ice sheet

from sys import exit
from pylab import *

## the following code defining the Test N form of the Bodvardsson solution
## represents complete code duplication from
##   pism-dev/src/verif/tests/exactTestN.[h|c]
## see comments there

secpera = 31556926.0
g       = 9.81
rho     = 910.0
rhow    = 1028.0
n       = 3.0

H0      = 3000.0
L0      = 500.0e3

Hela    = H0 / 1.5
a       = 0.003 / secpera
k       = 9.0 * Hela / (a * L0 * L0)

xc      = 0.92 * L0
Hc      = H0 * (1.0 - (xc / L0) * (xc / L0))
bc      = rho * Hc / rhow

def exactN(x):
  if (x < 0.0) | (x > L0):
    print "ERROR in exactN()"
    return 1
  hxx = - 2.0 * H0 / (L0 * L0)
  q   = (1.0 / n) - 1.0
  ux  = - hxx / k
  T0  = 0.5 * (1.0 - rho / rhow) * rho * g * Hc * Hc
  H   = H0 * (1.0 - (x / L0) * (x / L0))
  hx  = hxx * x
  u   = - hx / k
  M   = a * (H - Hela)
  B   = T0 / ( 2.0 * H * (abs(ux)**q) * ux )
  return H, hx, u, M, B

