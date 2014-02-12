#!/usr/bin/env python

## generate figure showing the Bodvardsson solution

from sys import exit
from pylab import *

import exactN

from plottool import lw, plottwo

N = 1001
xx = linspace(0.0,exactN.xc,N)

HH, hx, uu, MM = exactN.exactN(xx)
_, BB = exactN.exactNbueler(xx)

# surface elevation is so that ice reaches flotation exactly at xc
hh = HH - (exactN.rho/exactN.rhow) * exactN.Hc

# figure 1: profile and velocity (two-axis)
fig = plt.figure(figsize=(6,3))
ax1fig1, ax2fig1 = plottwo(fig, xx/1000.0,HH,uu * exactN.secpera,
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
fig = plt.figure(figsize=(6,3))
ax1fig2, ax2fig2 = plottwo(fig,xx/1000.0, bbeta * 1.0e-10, BB * 1.0e-8,
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

