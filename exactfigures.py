#!/usr/bin/env python

## generate figure showing the Bodvardsson solution

from sys import exit
from pylab import *

import exactN

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

  return ax1, ax2

N = 1001
xx = linspace(0.0,exactN.xc,N)
HH = xx.copy()
uu = xx.copy()
MM = xx.copy()
BB = xx.copy()
for j in range(N):
  HH[j], hx, uu[j], MM[j], BB[j] = exactN.exactN(xx[j])

# surface elevation is so that ice reaches flotation exactly at xc
hh = HH - (exactN.rho/exactN.rhow) * exactN.Hc

# figure 1: profile and velocity (two-axis)
ax1fig1, ax2fig1 = plottwo(xx/1000.0,HH,uu * exactN.secpera,
                           r"$H=$ ice thickness (m, solid)",r"$u=$ ice velocity (m a-1, dashed)",
                           y1min=0.0)
ax1fig1.hold(True)
# show water height as waves
xl = linspace(exactN.xc,exactN.L0,N)
ax1fig1.plot(xl/1000.0, (exactN.rho/exactN.rhow) * exactN.Hc - 30.0 * abs(sin(xl / 4.0e3)),'k',linewidth=0.7*lw)
y1, y2 = ax1fig1.get_ylim()
ax1fig1.set_ylim(0.0,y2)
ax1fig1.set_xlim(0.0,exactN.L0/1000.0)
ax1fig1.plot([exactN.xc/1000.0,exactN.xc/1000.0],[0.0,exactN.Hc],'k',linewidth=lw)
ax1fig1.text((exactN.xc-10.0e3)/1000.0,-200.0,r"$x_c$",fontsize=18)
ax1fig1.plot([0.0,8.0],[exactN.bc,exactN.bc],'k',linewidth=0.5*lw)
ax1fig1.text(+12.0,0.8*exactN.bc,r"$b_c$",fontsize=18)
ax1fig1.hold(False)
name = 'bodverifthickvel.pdf'
print '  generating figure %s ...' % name
savefig(name,bbox_inches='tight')

# figure 2: sliding \beta and hardness B (two-axis)
bbeta = exactN.k * exactN.rho * exactN.g * HH
ax1fig2, ax2fig2 = plottwo(xx/1000.0, bbeta * 1.0e-10, BB * 1.0e-8,
        r"$\beta$ ($10^{10}$ Pa s m-1, solid)",
        r"$B$ ($10^8$ Pa s-(1/3), dashed)")
y1, y2 = ax1fig2.get_ylim()
ax1fig2.set_ylim(0.0,y2)
y1, y2 = ax2fig2.get_ylim()
ax2fig2.set_ylim(0.0,y2)
ax1fig2.hold(True)
ax1fig2.plot(xl/1000.0, 0.0*xl,'k',alpha=0.0)
ax1fig2.text((exactN.xc-10.0e3)/1000.0,-0.25,r"$x_c$",fontsize=18)
ax1fig2.hold(False)
name = 'bodverifbetaB.pdf'
print '  generating figure %s ...' % name
savefig(name,bbox_inches='tight')

exit(0)

#figure 3: SMB   UNNEEDED
figure(3,figsize=(8.0,3.5))
plot(xx/1000.0, MM * exactN.secpera, 'k', linewidth=lw)
xlabel("x  (km)")
ylabel("surface mass balance (m a-1)")
savefig('bodverifsmb.pdf',bbox_inches='tight')

