#!/bin/bash

# run as 
#   $ ./conv.sh >& conv.txt

MPIDO="mpiexec -n 1"

# dx =   10km  5km  2km  1km  500m  200m  100m
for N in   46   91  226  451   901  2251  4501
do
  # second-order upwinding (default)
  #./bodvardsson -snes_fd -da_grid_x $N -snes_monitor -snes_rtol 1.0e-8
  # first-order upwinding:
  #  ./bodvardsson -snes_fd -da_grid_x $N -snes_monitor -snes_rtol 1.0e-8 -bod_up_one

  # second-order upwinding (default) with exact soln as initial guess:
  $MPIDO ./bodvardsson -snes_fd -snes_max_funcs 100000 -da_grid_x $N -snes_monitor -snes_rtol 1.0e-8 -bod_exact_init
done

