#!/bin/bash

# run as
#   $ ./convmarine-realistic.sh |'grep' errHinf
#(dx,errHinf,erruinf) 20000.000 1.0984e+02 1.3743e+02
#(dx,errHinf,erruinf) 9873.418 2.9679e+01 2.7486e+01
#(dx,errHinf,erruinf) 4968.153 5.7678e+00 4.7760e+00
#(dx,errHinf,erruinf) 1994.885 1.0727e+00 8.2240e-01
#(dx,errHinf,erruinf) 998.720 8.7900e+00 6.9828e+00

#MPIDO="mpiexec -n 4"
MPIDO="mpiexec -n 1"

for dx in 20000  10000  5000  2000  1000
do
  # use f.d. for jacobian, and realistic (slab) initial guess:
  $MPIDO ./marine -snes_fd_color -dx $dx -snes_converged_reason -snes_rtol 1.0e-8
done

