#!/usr/bin/env python

from pylab import *
from sys import exit

from plhalfar import Params, halfar_H, halfar_t0

p = Params()

nn = array([3, 10, 50, inf])  # values of n in Glen law to compare

r = linspace(0.0,1.15*p.R0,1001)
rr = r / 1000.0

# figure 1
figure(figsize=(6.0,2.0))
M = len(nn)
for j in range(0,M):
  if isinf(nn[j]):
    H = halfar_H(1.0,r,p,nn[j])
    plot(rr,H/1000.0,'k',linewidth=1.5)
  else:
    t0 = halfar_t0(p,nn[j])
    H = halfar_H(t0,r,p,nn[j])
    plot(rr,H/1000.0,'k',linewidth=1.5)
  hold(True)
hold(False)
axis([min(rr),max(rr),0.0,4.0])
setp(gca(),yticks=[0.0,1.0,2.0,3.0,4.0])
ylabel('H   (km)')
xlabel('r   (km)')
#show()
savefig('figure1.pdf',bbox_inches='tight')

