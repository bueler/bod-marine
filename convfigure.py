#!/usr/bin/env python

# see convmarine-*.sh in petsc

# run as
#   $ ./convfigure.py petsc/convmarine-exactinit.sh petsc/convmarine-realistic.sh

from sys import stderr, exit, argv
from pylab import *

if len(argv) != 3:
  print('ERROR: requires two text files with convergence results to run ...')
  exit(1)

clen = zeros((2,1))
dx = zeros((2,12))   # at most 20 results in each convergence study
errHinf = dx.copy()
erruinf = dx.copy()

for fnum in range(2):
    print("reading text file %s ..." % argv[fnum+1]),
    f = file(argv[fnum+1],'r')
    count = 0
    while True:
      try:
        line = f.next().split(' ')
      except:
        break
      if line[0].find('#(dx,') >= 0:
        #print line
        if len(line) < 4:
          print "ERROR: line needs 4 objects separated by spaces"
          exit(1)
        dx[fnum][count] = float(line[1])
        errHinf[fnum][count] = float(line[2])
        erruinf[fnum][count] = float(line[3])
        count = count + 1
    f.close()
    clen[fnum] = count
    print("done")

fig = figure(figsize=(7,6))

sym = ('ko','k*')
symsize = (10.0,12.0)
symface = ('k','w')
dxmin = dx[0][:clen[0][0]].min()
dxmax = dx[1].max()

# upper subplot shows thickness errors
ax1 = fig.add_subplot(2,1,1)
for fnum in range(2):
    N = clen[fnum][0]
    loglog(dx[fnum][:N],errHinf[fnum][:N],sym[fnum],markersize=symsize[fnum],markerfacecolor=symface[fnum])
    hold(True)
    if fnum == 0:
      p = polyfit(log(dx[fnum][:N]),log(errHinf[fnum][:N]),1)
      print "result:  thickness error %d decays at rate O(dx^%.2f)" % (fnum,p[0])
      loglog(dx[fnum][:N],exp(polyval(p,log(dx[fnum][:N]))),'k:')
hold(False)
grid(True)
miny = errHinf[0][:clen[0][0]].min()
maxy = errHinf[1][:clen[1][0]].max()
axis((0.6*dxmin,1.7*dxmax,0.3*miny,3.0*maxy))
xlabel('')
xticks(array([10.0, 100.0, 1000.0, 10000.0]),
       ('', '', '', '') )
ylabel(r'$H$ error  (m)')

# lower subplot shows velocity errors  FIXME
ax2 = fig.add_subplot(2,1,2)
for fnum in range(2):
    N = clen[fnum][0]
    loglog(dx[fnum][:N],erruinf[fnum][:N],sym[fnum],markersize=symsize[fnum],markerfacecolor=symface[fnum])
    hold(True)
    if fnum == 0:
      p = polyfit(log(dx[fnum][:N]),log(erruinf[fnum][:N]),1)
      print "result:  velocity error %d decays at rate O(dx^%.2f)" % (fnum,p[0])
      loglog(dx[fnum][:N],exp(polyval(p,log(dx[fnum][:N]))),'k:')
hold(False)
grid(True)
miny = erruinf[0][:clen[0][0]].min()
maxy = erruinf[1][:clen[1][0]].max()
axis((0.6*dxmin,1.7*dxmax,0.3*miny,3.0*maxy))
xlabel('')
xticks(array([10.0, 100.0, 1000.0, 10000.0]),
       ('', '', '', '') )
ylabel(r'$u$ error  (m/a)')

# only lower x-axis gets labels
xlabel(r'$\Delta x$')
xticks(array([10.0, 100.0, 1000.0, 10000.0]),
       ('10m', '100m', '1km', '10km') )

figname = 'convmarine.pdf'
print "saving figure '%s'" % figname
savefig(figname,bbox_inches='tight')

