#!/usr/bin/env python

## generate figure showing the Bodvardsson (1955) problem thought of as a
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
n       = 3.0

H0      = 3000.0
L0      = 500.0e3

xc      = 0.9 * L0
Hela    = H0 / 1.5

hxx     = - 2.0 * H0 / (L0 * L0)
q       = (1.0 / n) - 1.0
a       = 0.003 / secpera
k       = 9.0 * Hela / (a * L0 * L0)
ux      = - hxx / k

rhow    = 1028.0

Hc      = H0 * (1.0 - (xc / L0) * (xc / L0))
T0      = 0.5 * (1.0 - rho / rhow) * rho * g * Hc * Hc

def exactN(x):
  if (x < 0.0) | (x > L0):
    print "ERROR in exactN()"
    return 1
  H = H0 * (1.0 - (x / L0) * (x / L0))
  hx = hxx * x
  u = - hx / k
  M = a * (H - Hela)
  B = T0 / ( 2.0 * H * (abs(ux)**q) * ux )
  return H, hx, u, M, B


lw = 2.0

# for two-axes technique, see
#   http://matplotlib.sourceforge.net/plot_directive/mpl_examples/api/two_scales.py
#   http://matplotlib.sourceforge.net/plot_directive/mpl_examples/api/fahrenheit_celcius_scales.py

def plottwo(x, y1, y2, label1, label2, y1min=-inf, y1max=inf, y2min=-inf, y2max=inf):
  import matplotlib.pyplot as plt

  fig = plt.figure(figsize=(6,3))
  ax1 = fig.add_subplot(111) # left axis
  ax2 = ax1.twinx()          # right axis

  ax1.plot(x,y1,'k',linewidth=lw)
  ax2.plot(x,y2,'k--',linewidth=lw)

  ax1.set_xlabel("x  (km)")

  #y1, y2 = ax1.get_ylim()
  
  if y1min == -inf:
    y1min = y1.min()
  if y1max == inf:
    y1max = y1.max()
  ax1.set_ylim(y1min, y1max)

  if y2min == -inf:
    y2min = y2.min()
  if y2max == inf:
    y2max = y2.max()
  ax2.set_ylim(y2min, y2max)

  # make the y-axis label and tick labels match the plot colors
  ax1.set_ylabel(label1, color='k')
  for tl in ax1.get_yticklabels():
      tl.set_color('k')
  ax2.set_ylabel(label2, color='k')
  for tl in ax2.get_yticklabels():
      tl.set_color('k')

  #ax1.set_xticks(array([0.0, 100.0e3, 200.0e3, 300.0e3, 400.0e3, xc]))
  #ax1.set_xticklabels(("0","100","200","300","400",r"$x_c$"))
  #ax2.set_xticks(array([0.0, 100.0e3, 200.0e3, 300.0e3, 400.0e3, xc]))
  #ax2.set_xticklabels(("0","100","200","300","400",r"$x_c$"))

  return ax1, ax2


N = 1001
xx = linspace(0.0,xc,N)
HH = xx.copy()
uu = xx.copy()
MM = xx.copy()
BB = xx.copy()
for j in range(N):
  HH[j], hx, uu[j], MM[j], BB[j] = exactN(xx[j])

hh = HH - (rho/rhow) * Hc  # surface elevation is so that ice reaches flotation
                           # exactly at xc

# figure 1: profile and velocity (two-axis)
ax1fig1, ax2fig1 = plottwo(xx/1000.0,HH,uu * secpera,
                           r"$H=$ ice thickness (m, solid)",r"$u=$ ice velocity (m a-1, dashed)",
                           y1min=0.0)
ax1fig1.hold(True)
xl = linspace(xc,L0,N)
ax1fig1.plot(xl/1000.0, (rho/rhow) * Hc - 30.0 * abs(sin(xl / 4.0e3)),'k',linewidth=0.7*lw)
y1, y2 = ax1fig1.get_ylim()
ax1fig1.set_ylim(0.0,y2)
ax1fig1.set_xlim(0.0,L0/1000.0)
ax1fig1.plot([xc/1000.0,xc/1000.0],[0.0,Hc],'k',linewidth=lw)
ax1fig1.text((xc-10.0e3)/1000.0,-300.0,r"$x_c$",fontsize=18)
ax1fig1.hold(False)
savefig('bodverifthickvel.pdf',bbox_inches='tight')

# figure 2: sliding \beta and hardness B (two-axis)
bbeta = k * rho * g * HH
ax1fig2, ax2fig2 = plottwo(xx/1000.0, bbeta * 1.0e-10, BB * 1.0e-8,
        r"$\beta$ ($10^{10}$ Pa s m-1, solid)",
        r"$B$ ($10^8$ Pa s-(1/3), dashed)")
ax1fig2.hold(True)
ax1fig2.plot(xl/1000.0, 0.0*xl,'k',alpha=0.0)
ax1fig2.text((xc-10.0e3)/1000.0,-0.25,r"$x_c$",fontsize=18)
ax1fig2.hold(False)
savefig('bodverifbetaB.pdf',bbox_inches='tight')

exit(0)

#figure 3: SMB   UNNEEDED
figure(3,figsize=(8.0,3.5))
plot(xx/1000.0, MM * secpera, 'k', linewidth=lw)
xlabel("x  (km)")
ylabel("surface mass balance (m a-1)")
savefig('bodverifsmb.pdf',bbox_inches='tight')

