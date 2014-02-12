#!/usr/bin/env python

from pylab import *
from sys import exit

H0 = 3000.0
L0 = 500.0e3

x = linspace(0.0,1.05*L0,1001)
xx = x / 1000.0

# get Weertman (1961) shape
xw = x.copy()
xw[x > L0] = L0
Hw = H0 * ( 1.0 - (xw / L0) )**0.5

# get Bodvardsson (1955), equation (23), shape
Hela = H0 / 1.5
Hb = 1.5 * Hela * ( 1.0 - (x / L0)**2.0 )
Hb[x > L0] = 0.0

figure(figsize=(6,3.5))
plot(xx,Hw,'k:',linewidth=1.5)
hold(True)
plot(xx,Hb,'k',linewidth=1.5)
hold(False)
axis([min(xx),max(xx),0.0,3.5])
setp(gca(),yticks=[0.0,500.0,1000.0,1500.0,2000.0,2500.0,3000.0])
ylabel('H   (m)')
xlabel('x   (km)')
savefig('twoparabolas.pdf',bbox_inches='tight')

