#!/usr/bin/env python

# Show Halfar solution (lambda=0) in time ranges for values of n.
# Formula and constants reference: Bueler et al 2005.
# Generates Figure 2 in paper, but gives tools for figures 1 and 3.

from pylab import *
from sys import exit

class Params:
    """Holds parameters for SIA profiles, but not including Glen's n."""
    def __init__(self):
      self.H0      = 3600.0         # m
      self.R0      = 750.0e3        # m;  = 750 km
      self.secpera = 31556926.0     # secs per year
      self.tau0    = 1.0e5          # 10^5 Pa = 1 bar
      self.edot0   = 0.1 / self.secpera  # 0.1 a-1; strain rate for A = 10^-16 Pa-3 a-1
      self.rg      = 910.0 * 9.81   # = rho * g; result units kg m-2 s-2 = Pa m-1

def halfar_t0(p, n):       # Characteristic time t_0 in the n < inf case.
  z = ( p.tau0 * p.R0 * (2.0*n+1.0) / (p.rg * p.H0**2.0 * (n+1.0)) )**n
  return z * p.R0 * (n+2.0) / (2.0 * p.H0 * p.edot0 * (5.0*n+3.0))

def halfar_H(t, r_in, p, n): # Profile H(t,r) for any 1 <= n <= inf
  r = r_in.copy()
  if isinf(n): # perfectly-plastic: ignor t in this case
    r[r > p.R0] = p.R0
    return p.H0 * ( 1.0 - (r / p.R0) )**0.5
  else:
    Tmbeta = (t / halfar_t0(p,n))**(- 1.0 / (5.0*n + 3.0))
    r[r > (p.R0/Tmbeta)] = p.R0/Tmbeta
    P = n / (2.0*n+1.0)   # outer power
    Q = (n+1.0) / n       # inner power
    return p.H0 * (Tmbeta**2.0) * ( 1.0 - (Tmbeta * r / p.R0)**Q )**P


if __name__ == "__main__":

  p = Params()
  
  #nn = array([3, 10, 50, inf])  # values of n in Glen law to compare
  nn = array([3, 10, 50])  # values of n in Glen law to compare in figure
  
  below = 0.001     # below * t0 < t < above * t0
  above = 1000.0
  
  r = linspace(0.0,1.5*p.R0,1001)
  rr = r / 1000.0
  
  # make figure 2 in paper
  figure(figsize=(6.0,9.0))
  M = len(nn)
  for j in range(0,M):
    subplot(M,1,j+1)
    if isinf(nn[j]):
      H = halfar_H(1.0,r,p,nn[j])
      plot(rr,H/1000.0,'k',linewidth=2.0)
      text((p.R0/10.0)/1000.0,p.H0/4000.0,'plastic (n = $\infty$)')
    else:
      t0 = halfar_t0(p,nn[j])
      H = halfar_H(t0,r,p,nn[j])
      Hbelow = halfar_H(below * t0,r,p,nn[j])
      Habove = halfar_H(above * t0,r,p,nn[j])
      plot(rr,Hbelow/1000.0,'g',linewidth=1.5)
      hold(True)
      fill(concatenate((rr,rr[::-1])),
           concatenate((Hbelow/1000.0,H[::-1]/1000.0)),'g',alpha=0.5)
      plot(rr,Habove/1000.0,'r',linewidth=1.5)
      fill(concatenate((rr,rr[::-1])),
           concatenate((H/1000.0,Habove[::-1]/1000.0)),'r',alpha=0.5)
      plot(rr,H/1000.0,'k',linewidth=2.0)
      hold(False)
      text((p.R0/10.0)/1000.0,p.H0/4000.0,'n = ' + str(nn[j]))
    axis([min(rr),max(rr),0.0,4.0])
    setp(gca(),yticks=[0.0,1.0,2.0,3.0,4.0])
    ylabel('H   (km)')
    if j == M-1:
      xlabel('r   (km)')
    if j < M-1:
      setp(gca(),xticklabels=[])

  #show()
  #savefig('plhalfar.png',dpi=300,bbox_inches='tight')
  savefig('plhalfar.pdf',bbox_inches='tight')

