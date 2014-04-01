#!/bin/bash

# run as
#   $ ./convmarine-exactinit.sh |'grep' errHinf
#(dx,errHinf,erruinf) 4968.153 5.7678e+00 4.7760e+00
#(dx,errHinf,erruinf) 1994.885 1.0727e+00 8.2240e-01
#(dx,errHinf,erruinf) 998.720 4.9633e-01 3.8600e-01
#(dx,errHinf,erruinf) 499.680 2.3843e-01 1.8681e-01
#(dx,errHinf,erruinf) 199.949 9.3042e-02 7.3229e-02
#(dx,errHinf,erruinf) 99.987 4.6111e-02 3.6348e-02
#(dx,errHinf,erruinf) 49.997 2.2931e-02 1.8090e-02
#(dx,errHinf,erruinf) 19.999 9.1142e-03 7.1932e-03
#(dx,errHinf,erruinf) 10.000 4.5237e-03 3.5708e-03
#(dx,errHinf,erruinf) 5.000 2.2313e-03 1.7614e-03

#MPIDO="mpiexec -n 4"
MPIDO="mpiexec -n 1"

for dx in 5000  2000  1000  500  200  100  50  20  10  5
do
  # use exact, unscaled jacobian and exact soln as initial guess:
  $MPIDO ./marine -snes_max_funcs 100000 -dx $dx -snes_monitor -snes_rtol 1.0e-8 -exactinit -noscale
done

