/*
   Copyright (C) 2010-2014 Ed Bueler
*/

#ifndef __exactsolns_h
#define __exactsolns_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactBod() is a C implementation of the parabolic solution in 
! Bodvardsson (1955), treated here as a manufactured exact solution to
! a steady-state SSA flow problem, including the mass continuity equation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      The Bodvardsson (1955) can be thought of as solving mass continuity
      and SSA stress balance simultaneously:

         M(x) - (u H)_x = 0

         ( 2 H B(x) |u_x|^((1/n)-1) u_x )_x - beta0(x) u = rho g H h_x
         
      Here H = H(x) is ice thickness and u = u(x) is ice velocity.  Following
      Bodvardsson, the equilibrium line altitude, surface mass balance, and
      sliding coefficient are:

         Hela = H0 / 1.5
         M(x) = a (h(x) - Hela)
         beta0(x) = k rho g H(x)
      
      The boundary condition at x = xg implies that the calving front is
      exactly at the location where the ice sheet reaches flotation.

*/


int params_exactBod(double *H0, double *L0, double *xg,
                    double *a, double *Hela, double *k);
   /* outputs: H0   = dome thickness (m)
               L0   = full flow-line length from dome to margin where H->0 (m)
               xg   = in Bueler interpretation, the location of the calving front (m)
               a    = surface mass balance lapse rate, with elevation (s-1)
               Hela = elevation of equilibrium line (m)
               k    = coefficient for sliding                               */

int exactBod(double x, double *H, double *u, double *M);
   /* input    : x                   (m; 0.0 <= x <= L0)

      output   : H = H(x)            (m; ice thickness)
                 u = u(x)            (m s-1; ice velocity)
                 M = M(x)            (m s-1; surface mass balance)

      Assumes n = 3.

      return value =
         0 if successful
         1 if x < 0
         2 if x > L0                                                        */


int exactBodBueler(double x, double *T, double *B);
   /* input    : x                   (m; 0.0 <= x <= L0)

      output   : T = T(x)
                 B = B(x)            (Pa s^(1/3); ice hardness)
      
      In Bueler interpretation, T(x) and B(x) are constructed.

      return value =
         0 if successful
         1 if x < 0
         2 if x > L0                                                        */


int exactVeen(double x, double M0, double *H, double *u);
   /* input    : x                   (m; xg <= x)

      output   : H = H(x)            (m; thickness)
                 u = u(x)            (m s-1; ice velocity)

      van der Veen: get thickness and velocity for floating ice shelf, given
      Mg, Bg, xg, Hg, ug computed above from Bodvardsson

      return value =
         0 if successful
         1 if x < xg
         2 if u <= 0
         3 if H <= 0                                                        */

#ifdef __cplusplus
}
#endif

#endif  /* __exactsolns_h */

