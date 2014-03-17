#!/bin/bash

# run as
#   $ ./convmarine.sh |'grep' errHinf
#(dx,errHinf,erruinf) 9873.418 2.9679e+01 2.7486e+01
#(dx,errHinf,erruinf) 4968.153 5.7678e+00 4.7760e+00
#(dx,errHinf,erruinf) 1994.885 1.0727e+00 8.2240e-01
#(dx,errHinf,erruinf) 998.720 4.9633e-01 3.8600e-01
#(dx,errHinf,erruinf) 499.680 2.3843e-01 1.8681e-01
#(dx,errHinf,erruinf) 199.949 9.3042e-02 7.3229e-02

#MPIDO="mpiexec -n 4"
MPIDO="mpiexec -n 1"

#for dx in 10000  5000  2000  1000  500
for dx in 10000  5000  2000  1000  500 200 100 50 20 10
do
  # use exact soln as initial guess:
  #$MPIDO ./marine -snes_fd -snes_max_funcs 100000 -dx $dx -snes_monitor -snes_rtol 1.0e-8 -exactinit
  $MPIDO ./marine -snes_max_funcs 100000 -dx $dx -snes_monitor -snes_rtol 1.0e-8 -exactinit -noscale
done

