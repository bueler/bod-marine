#!/usr/bin/env python

# run conv.sh in pism-dev/src/verif/tests/verifN/:
#   $ ./conv.sh >& conv.txt
# then link:
#   $ ln -s ~/pism-dev/src/verif/tests/verifN/conv.txt
# then run this script to get pyplot figs:
#   $ ./verifNfigure.py conv.txt verifN.pdf

from sys import stderr, exit, argv
from pylab import *

txtname = argv[1]
if len(argv) > 2:
  figname = argv[2]
else:
  figname = ""

print("reading text file %s ..." % txtname)

dx = []
errHinf = []
erruinf = []

f = file(txtname,'r')
count = 0
while True:
  try:
    line = f.next().split(' ')
  except:
    break
  if line[0].find('(dx,') >= 0:
    print line
    if len(line) < 4:
      print "ERROR: line needs 4"
      exit(1)
    dx.append(float(line[1]))
    errHinf.append(float(line[2]))
    erruinf.append(float(line[3]))
    count = count + 1
f.close()

dx = array(dx)
errHinf = array(errHinf)
erruinf = array(erruinf)

fig = figure(figsize=(7,4))
ax1 = fig.add_subplot(111)
loglog(dx,errHinf,'bo',markersize=9.0,label=r'$\max |H - H_{exact}|$  [m]')
hold(True)
loglog(dx,erruinf,'g*',markersize=12.0,label=r'$\max |u - u_{exact}|$  [m a-1]')
p = polyfit(log(dx),log(errHinf),1)
print "result:  thickness error decays at rate O(dx^%.2f)" % p[0]
loglog(dx,exp(polyval(p,log(dx))),'r:')
hold(False)
maxy = max(errHinf.max(),erruinf.max())
miny = min(errHinf.min(),erruinf.min())
dims = axis('tight')
axis((0.9*dx.min(),1.1*dx.max(),0.9*miny,1.1*maxy))
legend(loc='upper left')
xlabel(r'$\Delta x$  [m]')
xticks(array([300.0, 500.0, 1000.0, 2000.0, 3000.0, 5000.0, 10000.0, 30000.0]),
       ('300m','500m','1km','2km','3km','5km','10km','30km') )
grid(True)

if len(figname) == 0:
  show()
else:
  print "saving figure '%s'" % figname
  savefig(figname,bbox_inches='tight')

