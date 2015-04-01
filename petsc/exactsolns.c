/*
   Copyright (C) 2010-2014 Ed Bueler
*/

#include <math.h>
#include "exactsolns.h"

#define secpera  31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81
#define rho      910.0         /* ice density; kg m-3 */
#define rhow     1028.0        /* sea water density; kg m-3 */
#define n        3.0           /* Glen power */

int params_exactBod(double *H0, double *L0, double *xg,
                    double *a, double *Hela, double *k) {
  /* geometry */
  *H0   = 3000.0;
  *L0   = 500.0e3;
  *xg   = 0.9 * (*L0);
  /* mass balance */
  *a    = 0.003 / secpera;   /* s-1; mass balance gradient with elevation */
  *Hela = (*H0) / 1.5;       /* m;  H0 = 1.5 Hela  exactly */
  /* sliding */
  *k    = 9.0 * (*Hela) / ((*a) * (*L0) * (*L0)); /* s m-1; choose k so that eqn (24) gives our L0 */
  return 0;
}


int exactBod(double x, double *H, double *u, double *M) {

  double H0, L0, xg, a, Hela, k;
  double hx, hxx;
  params_exactBod(&H0, &L0, &xg, &a, &Hela, &k);

  if (x < 0.0) { return 1; }
  if (x > L0) { return 2; }

  hxx = - 2.0 * H0 / (L0 * L0);    /* constant concavity of h(x) */
  hx  = hxx * x;

  *H = H0 * (1.0 - (x / L0) * (x / L0));  /* eqn (23) in Bodvardsson */
  *u = - (hx) / k;                        /* eqn (10) in Bodvardson, once SIA is dropped */
  *M = a * ((*H) - Hela);                 /* page 6 in Bodvardsson, just before eqn (23) */
  return 0;
}


int exactBodBueler(double x, double *T, double *B) {

  const double omega = 1.0 - rho / rhow;
  double H0, L0, xg, a, Hela, k;
  double Hg, ug, Mg, H, u, M;
  double q, hxx, ux;
  int ierr;

  params_exactBod(&H0, &L0, &xg, &a, &Hela, &k);
  exactBod(xg, &Hg, &ug, &Mg);

  ierr = exactBod(x, &H, &u, &M);
  if (ierr)  return ierr;

  q   = (1.0 / n) - 1.0;           /* a useful power */
  hxx = - 2.0 * H0 / (L0 * L0);    /* constant concavity of h(x) */
  ux  = - hxx / k;                 /* constant strain rate */
  *T  = 0.5 * omega * rho * g * Hg * Hg;
  *B  = *T / ( 2.0 * H * pow(fabs(ux),q) * ux );
  return 0;
}


int exactVeen(double x, double M0, double *H, double *u) {

  const double omega = 1.0 - rho / rhow;
  double xg,Hg,ug,Bg,C,
         tmp1,tmp2,tmp3,tmp4,tmp5;

  params_exactBod(&tmp1, &tmp2, &xg, &tmp3, &tmp4, &tmp5);  // get xg
  if (x < xg) { return 1; }

  exactBod(xg, &Hg, &ug, &tmp1);
  exactBodBueler(xg, &tmp1, &Bg);

  C = pow(rho * g * omega / (4.0 * Bg), n);
  tmp1 = ug * Hg + M0 * (x - xg);
  tmp2 = pow(ug,n+1) + (C / M0) * ( pow(tmp1,n+1) - pow(ug * Hg,n+1) );
  *u = pow(tmp2,1.0/(n+1));
  if (*u <= 0.0) { return 2; }
  *H = tmp1 / *u;
  if (*H <= 0.0) { return 3; }
  return 0;
}

