#!/bin/bash

# run as
#   $ ./convmarine-realistic.sh |'grep' errHinf
#(dx,errHinf,erruinf) 20000.000 1.0984e+02 1.3743e+02
#(dx,errHinf,erruinf) 9873.418 7.3900e+01 5.6537e+01
#(dx,errHinf,erruinf) 4968.153 4.2258e+01 3.2890e+01
#(dx,errHinf,erruinf) 1994.885 3.3542e+01 2.6632e+01
#(dx,errHinf,erruinf) 998.720 1.7569e+01 1.3896e+01

#MPIDO="mpiexec -n 4"
MPIDO="mpiexec -n 1"

for dx in 20000  10000  5000  2000  1000
do
  # use f.d. for jacobian, and realistic (slab) initial guess:
  $MPIDO ./marine -snes_fd -snes_max_funcs 100000 -dx $dx -snes_monitor -snes_rtol 1.0e-8
done

