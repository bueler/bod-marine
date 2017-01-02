#!/bin/bash

# run as
# $ ./convbod.sh |'grep' errHinf
#(dx,errHinf,erruinf) 10000.000 2.2529e+00 7.1132e-01
#(dx,errHinf,erruinf) 5000.000 5.8056e-01 1.7846e-01
#(dx,errHinf,erruinf) 2000.000 9.4591e-02 2.8663e-02
#(dx,errHinf,erruinf) 1000.000 2.3800e-02 7.1103e-03
#(dx,errHinf,erruinf) 500.000 5.9777e-03 1.7051e-03
#(dx,errHinf,erruinf) 200.000 9.6912e-04 1.9599e-04
#(dx,errHinf,erruinf) 100.000 2.7495e-04 2.7637e-05

MPIDO="mpiexec -n 1"

# dx =   10km  5km  2km  1km  500m  200m  100m
for N in   46   91  226  451   901  2251  4501
do
  # second-order upwinding (default) with exact soln as initial guess:
  $MPIDO ./bodvardsson -snes_fd_color -da_grid_x $N -snes_converged_reason -snes_rtol 1.0e-8 -ksp_rtol 1.0e-10 -bod_exact_init
done

