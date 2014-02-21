static const char help[] =
"Solve steady flowline marine ice sheet equations using SNES, a dof=2 Vec\n"
"holding both thickness H and velocity u, and a 2nd-order finite difference\n"
"scheme to approximate the coupled mass continuity and SSA stress balance PDEs.\n"
"\n"
"IT WORKS!: this one is with a 'fair' initial condition:\n"
"  ./marine -snes_fd -dx 1000 -snes_monitor -snes_monitor_solution -draw_pause 0.1 -snes_max_funcs 1000000\n"
"A refinement path:\n"
"  for DX in 8000 4000 2000 1000 500 250; do ./marine -snes_fd -dx $DX -snes_max_funcs 100000 -exact_init; done\n"
"\n\n";

#include <petscdmda.h>
#include <petscsnes.h>

#include "exactsolns.h"


/* we will use dof=2 DMDA, and at each grid point have a thickness H and a velocity u */
typedef struct {
  PetscReal H, u;
} Node;

/* put info on exact solution in its own struct */
typedef struct {
  PetscReal   L0;  /* free boundary length in Bodvardsson solution */
  Vec         Hu;  /* exact thickness (Hu[i][0]) and exact velocity (Hu[i][1]) on regular grid */
  PetscReal   xg;  /* exact grounding line location */
  PetscReal   Hg;  /* exact thickness at grounding line location */
  PetscReal   Mg;  /* exact mass balance at grounding line location */
  PetscReal   Bg;  /* exact ice hardness at grounding line location */
} ExactCtx;

/* User-defined application context.  Filled by Fill...() and used by
   FunctionLocal() and JacobianLocal().  */
typedef struct {
  DM          da;      /* 1d,dof=2 distributed array for soln and residual */
  DM          stagda;  /* 1d,dof=1 distributed array for suitable for parameters at staggered grid points */
  PetscMPIInt rank;
  PetscInt    Mx, N, xs, xm;
  PetscReal   dx, secpera, n, rho, rhow, omega, g;
  PetscReal   xa;      /* location at which Dirichlet conditions are applied */
  PetscReal   xc;      /* calving front location at which stress (Neumann) condition applied to SSA eqn */
  PetscReal   Ha;      /* thickness at x=xa, for Dirichlet condition on mass cont */
  PetscReal   ua;      /* velocity at x=xa, for Dirichlet condition on mass cont */
  PetscReal   zocean;  /* surface elevation of ocean; bedrock is at zero elevation */
  PetscReal   k;       /* sliding parameter */
  PetscReal   Mfloat;  /* mass balance rate in floating ice */
  PetscReal   epsilon; /* regularization of viscosity, a strain rate */
  Vec         Mstag;   /* surface mass balance on staggered grid */
  Vec         Bstag;   /* ice hardness on staggered grid*/
} AppCtx;


extern PetscErrorCode FillExactSoln(ExactCtx*,AppCtx*);
extern PetscErrorCode FillInitial(AppCtx*,Vec*);
extern PetscErrorCode FillDistributedParams(ExactCtx*,AppCtx*);
extern PetscErrorCode FunctionLocalGLREG(DMDALocalInfo*,Node*,Node*,AppCtx*);
extern PetscErrorCode scshell(DMDALocalInfo*,Node*,Node*,AppCtx*);
/* extern PetscErrorCode JacobianMatrixLocal(DMDALocalInfo*,Node*,Mat,AppCtx*); */


int main(int argc,char **argv)
{
  PetscErrorCode         ierr;

  SNES                   snes;                 /* nonlinear solver */
  Vec                    Hu,r;                 /* solution, residual vectors */
  AppCtx                 user;                 /* user-defined work context */
  ExactCtx               exact;
  PetscInt               its;                  /* snes reports iteration count */
  SNESConvergedReason    reason;               /* snes reports convergence */
  PetscReal              tmp1, tmp2, tmp3,
                         errnorms[2], scaleNode[2], descaleNode[2];
  PetscInt               i;
  PetscBool              eps_set = PETSC_FALSE,
                         dump = PETSC_FALSE,
                         exactinitial = PETSC_FALSE,
                         snes_mf_set, snes_fd_set, dx_set;

  PetscInitialize(&argc,&argv,(char *)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &user.rank); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "MARINE solves for thickness and velocity in 1D, steady marine ice sheet\n"
    "  [run with -help for info and options]\n");CHKERRQ(ierr);

  user.n       = 3.0;          /* Glen flow law exponent */
  user.secpera = 31556926.0;
  user.rho     = 910.0;        /* kg m^-3 */
  user.rhow    = 1028.0;       /* kg m^-3 */
  user.omega   = 1.0 - user.rho / user.rhow;
  user.g       = 9.81;         /* m s^-2 */

  /* get parameters of exact solution */
  ierr = params_exactBod(&tmp1, &(exact.L0), &(exact.xg), &tmp2, &tmp3, &(user.k)); CHKERRQ(ierr);
  ierr = exactBod(exact.xg,&(exact.Hg),&tmp2,&(exact.Mg)); CHKERRQ(ierr);
  ierr = exactBodBueler(exact.xg,&tmp1,&(exact.Bg)); CHKERRQ(ierr);
  user.zocean = user.rho * exact.Hg / user.rhow;

/* see ../marineshoot.py: */
#define xa_default    0.2
#define xc_default    0.98
#define Mdrop_default 1.0

  /* define interval [xa,xc] */
  user.xa = xa_default * exact.L0;
  user.xc = xc_default * exact.L0;

  /* get Dirichlet boundary conditions, and mass balance on shelf */
  ierr = exactBod(user.xa, &(user.Ha), &(user.ua), &tmp1); CHKERRQ(ierr);
  user.Mfloat = Mdrop_default * exact.Mg;

  /* regularize using strain rate of 1/(length) per year */
  user.epsilon = (1.0 / user.secpera) / (user.xc - user.xa);
  /*user.epsilon = 0.0;*/

  user.dx = 10000.0;  /* default to coarse 10 km grid */

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to marine (steady marine ice sheet solver)","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsBool("-snes_mf","","",PETSC_FALSE,&snes_mf_set,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-snes_fd","","",PETSC_FALSE,&snes_fd_set,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-dx","grid spacing (m)","",
                            user.dx,&user.dx,&dx_set);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-exact_init",
             "initialize using exact solution instead of default linear function","",
             PETSC_FALSE,&exactinitial,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-dump_solns",
             "dump out exact and approximate solution and residual, as asci","",
             dump,&dump,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-epsilon","regularizing strain rate for stress computation (a-1)","",
                            user.epsilon * user.secpera,&user.epsilon,&eps_set);CHKERRQ(ierr);
    if (eps_set)  user.epsilon *= 1.0 / user.secpera;
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  if (!snes_mf_set && !snes_fd_set) {
    PetscPrintf(PETSC_COMM_WORLD,
       "\n***ERROR: marine needs one or zero of '-snes_mf', '-snes_fd'***\n\n"
       "USAGE FOLLOWS ...\n\n%s",help);
    PetscEnd();
  }

  if (snes_fd_set) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
       "  using approximate Jacobian; finite-differencing using coloring\n");
       CHKERRQ(ierr);
  } else if (snes_mf_set) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
       "  matrix free; no preconditioner\n"); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
       "  true Jacobian\n"); CHKERRQ(ierr);
  }

  if (dx_set && (user.dx <= 0.0)) {
    PetscPrintf(PETSC_COMM_WORLD,
       "\n***ERROR: -dx value unreasonable ... USAGE FOLLOWS ...\n\n%s",help);
    PetscEnd();
  }
  user.N = (int)ceil( ((user.xc - user.xa) / user.dx) - 0.5 );
  user.Mx = user.N + 2;
  user.dx = (user.xc - user.xa) / ((PetscReal)(user.N) + 0.5);  /* recompute so  dx * (N+1/2) = xc - xa  */

  /* Create machinery for parallel grid management (DMDA), nonlinear solver (SNES),
     and Vecs for fields (solution, RHS).  Note default grid points.  Also degrees
     of freedom = 2 (thickness and velocity at each point).  */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,-user.Mx,2,1,PETSC_NULL,&user.da);
            CHKERRQ(ierr);
  ierr = DMSetApplicationContext(user.da,&user);CHKERRQ(ierr);

  ierr = DMDASetFieldName(user.da,0,"ice thickness [non-dimensional]"); CHKERRQ(ierr);
  ierr = DMDASetFieldName(user.da,1,"ice velocity [non-dimensional]"); CHKERRQ(ierr);

  ierr = DMSetFromOptions(user.da); CHKERRQ(ierr);

  ierr = DMDAGetInfo(user.da,PETSC_IGNORE,&user.Mx,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(user.da,&user.xs,PETSC_NULL,PETSC_NULL,&user.xm,PETSC_NULL,PETSC_NULL);
                   CHKERRQ(ierr);

  /* another DMDA for scalar staggered parameters; one fewer point */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,user.Mx-1,1,1,PETSC_NULL,&user.stagda);CHKERRQ(ierr);

  /* establish geometry on grid; note xa = x_0 and xc = x_{N+1/2} */
  ierr = DMDASetUniformCoordinates(user.da,    user.xa,            user.xc+0.5*user.dx,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(user.stagda,user.xa+0.5*user.dx,user.xc,            0.0,1.0,0.0,1.0);CHKERRQ(ierr);

  /* report on current grid */
  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "  grid:  Mx = N+2 = %D regular points,  dx = %.3f m,  xa = %.2f km,  xc = %.2f km\n",
      user.Mx, user.dx, user.xa/1000.0, user.xc/1000.0);CHKERRQ(ierr);

  /* Extract/allocate global vectors from DMDAs and duplicate for remaining same types */
  ierr = DMCreateGlobalVector(user.da,&Hu);CHKERRQ(ierr);
  ierr = VecSetBlockSize(Hu,2);CHKERRQ(ierr);
  ierr = VecDuplicate(Hu,&r);CHKERRQ(ierr); /* inherits block size */
  ierr = VecDuplicate(Hu,&exact.Hu);CHKERRQ(ierr); /* ditto */

  ierr = DMCreateGlobalVector(user.stagda,&user.Mstag);CHKERRQ(ierr);
  ierr = VecDuplicate(user.Mstag,&user.Bstag);CHKERRQ(ierr);

  /* set up snes */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,user.da);CHKERRQ(ierr);

  ierr = DMDASNESSetFunctionLocal(user.da,INSERT_VALUES,(DMDASNESFunction)scshell,&user);CHKERRQ(ierr);

  /*ierr = DMDASNESSetJacobianLocal(user.da,(DMDASNESJacobian)JacobianMatrixLocal,&user);CHKERRQ(ierr);*/

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* the exact thickness and exact ice velocity (user.uHexact) are known */
  ierr = FillExactSoln(&exact, &user); CHKERRQ(ierr);
  /* the exact solution allows setting M(x), B(x) */
  ierr = FillDistributedParams(&exact, &user);CHKERRQ(ierr);

  if (exactinitial) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  using exact solution as initial guess\n");
             CHKERRQ(ierr);
    /* the initial guess is the exact continuum solution */
    ierr = VecCopy(exact.Hu,Hu); CHKERRQ(ierr);
  } else {
    /* the initial guess is a linear solution */
    ierr = FillInitial(&user, &Hu); CHKERRQ(ierr);
  }

  /************ SOLVE NONLINEAR SYSTEM  ************/
  /* recall that RHS  r  is used internally by KSP, and is set by the SNES */
  scaleNode[0] = 1000.0;
  scaleNode[1] = (100.0 / user.secpera);
  for (i = 0; i < 2; i++)  descaleNode[i] = 1.0 / scaleNode[i];
  ierr = VecStrideScaleAll(Hu,descaleNode); CHKERRQ(ierr); /* de-dimensionalize initial guess */
  ierr = SNESSolve(snes,PETSC_NULL,Hu);CHKERRQ(ierr);
  ierr = VecStrideScaleAll(Hu,scaleNode); CHKERRQ(ierr); /* put back in "real" scale */


  /* minimal report on solve */
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  %s Number of Newton iterations = %D\n",
           SNESConvergedReasons[reason],its);CHKERRQ(ierr);

  if (dump) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing combined result Hu\n");CHKERRQ(ierr);
    ierr = VecView(Hu,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing combined exact result exact.Hu\n");CHKERRQ(ierr);
    ierr = VecView(exact.Hu,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing final combined residual (at Hu)\n");CHKERRQ(ierr);
    ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  }

  /* evaluate error relative to exact solution */
  ierr = VecAXPY(Hu,-1.0,exact.Hu); CHKERRQ(ierr);  /* Hu = - Huexact + Hu */
  ierr = VecStrideNormAll(Hu,NORM_INFINITY,errnorms); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "(dx,errHinf,erruinf) %.3f %.4e %.4e\n",
           user.dx,errnorms[0],errnorms[1]*user.secpera);CHKERRQ(ierr);

  ierr = VecDestroy(&Hu);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = VecDestroy(&(exact.Hu));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.Mstag));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.Bstag));CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  ierr = DMDestroy(&(user.da));CHKERRQ(ierr);
  ierr = DMDestroy(&(user.stagda));CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}


/*  Compute the exact thickness and velocity on the regular grid, and
the staggered-grid ice hardness, over the locally-owned part of the grid.  */
PetscErrorCode FillExactSoln(ExactCtx *exact, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      dum1, *x;
  Node           *Hu;
  DM             coord_da;
  Vec            coord_x;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDM(user->da, &coord_da); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->da, &coord_x); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(coord_da,coord_x,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da,exact->Hu,&Hu);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    if (x[i] < exact->xg) {
      /* grounded part: Bodvardsson formula */
      ierr = exactBod(x[i], &(Hu[i].H), &(Hu[i].u), &dum1); CHKERRQ(ierr);
    } else {
      /* floating part: van der Veen formula */
      ierr = exactVeen(x[i], user->Mfloat, &(Hu[i].H), &(Hu[i].u));
      if (ierr) PetscPrintf(PETSC_COMM_WORLD,"WARNING:  exactVeen() returns error %d at i=%d, x[i]=%.3f km\n",
                            ierr,i,x[i]/1000.0);
    }
  }
  ierr = DMDAVecRestoreArray(coord_da,coord_x,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da,exact->Hu,&Hu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*  Put a not-unreasonable initial guess in Hu. */
PetscErrorCode FillInitial(AppCtx *user, Vec *vHu)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      *x;
  Node           *Hu;
  const PetscReal
      Hcguess = 300.0,
      ucguess = (100.0/user->secpera),
      Hslope  = (Hcguess - user->Ha) / (user->xc - user->xa),
      uslope  = (ucguess - user->ua) / (user->xc - user->xa);

  DM             coord_da;
  Vec            coord_x;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDM(user->da, &coord_da); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->da, &coord_x); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(coord_da,coord_x,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da,*vHu,&Hu);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    /* use linear function which uses upstream Dirichlet b.c. but guesses
       at calving-front values */
    Hu[i].H = user->Ha + Hslope * (x[i] - user->xa);
    Hu[i].u = user->ua + uslope * (x[i] - user->xa);
  }
  ierr = DMDAVecRestoreArray(coord_da,coord_x,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da,*vHu,&Hu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*  Compute the values of the surface mass balance and ice hardness. */
PetscErrorCode FillDistributedParams(ExactCtx *exact, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      dum1, dum2, *xstag, *Mstag, *Bstag;
  DM             coord_da;
  Vec            coord_x;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDM(user->stagda, &coord_da); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->stagda, &coord_x); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(coord_da,coord_x,&xstag);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->stagda,user->Mstag,&Mstag);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->stagda,user->Bstag,&Bstag);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm - 1; i++) {  /* note "-1" at end */
    if (xstag[i] < exact->xg) {
      ierr = exactBod(xstag[i], &dum1, &dum2, &(Mstag[i])); CHKERRQ(ierr);
      ierr = exactBodBueler(xstag[i], &dum1, &(Bstag[i])); CHKERRQ(ierr);
    } else  {
      Mstag[i] = user->Mfloat;
      Bstag[i] = exact->Bg;
    }
  }
  ierr = DMDAVecRestoreArray(coord_da,coord_x,&xstag);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->stagda,user->Mstag,&Mstag);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->stagda,user->Bstag,&Bstag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* define a power of the strain rate:   F  \approx  |u_x|^q u_x
   note  F(ul,ur) = - F(ur,ul) */
static inline PetscReal GetFSR(PetscReal dx, PetscReal eps, PetscReal n,
                                 PetscReal ul, PetscReal ur) {
  PetscReal dudx = (ur - ul) / dx,
            q    = (1.0 / n) - 1.0;
  return PetscPowScalar(dudx * dudx + eps * eps, q / 2.0) * dudx;
}


static inline PetscReal GLREG(PetscReal H, PetscReal Hg, PetscReal eps) {
  PetscReal tmp = exp(-(H-Hg)/eps);
  return 1.0 / (1.0 + tmp);
}

static inline PetscReal dGLREG(PetscReal H, PetscReal Hg, PetscReal eps) {
  /* this version has *no* regularization! */
  if (H > Hg) return 1.0;
  else        return 0.0;
}


/* Evaluate residual f at current iterate Hu, with regularization at grounding
lint.  */
PetscErrorCode FunctionLocalGLREG(DMDALocalInfo *info, Node *Hu, Node *f, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscReal      rg = user->rho * user->g,
                 omega = user->omega,
                 dx = user->dx,
                 Hg = user->zocean * (user->rhow / user->rho);
  PetscReal      *Mstag, *Bstag,
                 duH, ul, Fl, Fr, Tl, Tr, tmp, Hl, Hr, hl, hr, dhdx, beta;
  PetscInt       i, Mx = info->mx;
  Vec            locBstag, locMstag;

  PetscFunctionBegin;

  /* we need stencil width on Bstag and Mstag */
  ierr = DMGetLocalVector(user->stagda,&locBstag);CHKERRQ(ierr);  /* do NOT destroy it */
  ierr = DMGetLocalVector(user->stagda,&locMstag);CHKERRQ(ierr);  /* do NOT destroy it */

  ierr = DMGlobalToLocalBegin(user->stagda,user->Bstag,INSERT_VALUES,locBstag); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->stagda,user->Bstag,INSERT_VALUES,locBstag); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->stagda,user->Mstag,INSERT_VALUES,locMstag); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->stagda,user->Mstag,INSERT_VALUES,locMstag); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->stagda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->stagda,locMstag,&Mstag);CHKERRQ(ierr);
  for (i = info->xs; i < info->xs + info->xm; i++) {

    /* MASS CONT */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].H = Hu[0].H - user->Ha;
    } else {
      /* centered */
      if (i == 1) {
        duH = Hu[i].u * Hu[i].H - user->ua * user->Ha;
      } else {
        duH = Hu[i].u * Hu[i].H - Hu[i-1].u * Hu[i-1].H;
      }
      /**** residual for mass continuity ****/
      f[i].H = duH - dx * Mstag[i-1];
    }

    /* SSA */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].u = Hu[0].u - user->ua;
    } else if (i == Mx-1) {
      Fr = GetFSR(dx,user->epsilon,user->n, Hu[i-1].u,Hu[i].u);
      f[i].u = 0.25 * omega * rg * (Hu[i-1].H + Hu[i].H) - 2.0 * Bstag[i-1] * Fr;
    } else {
      /* residual: SSA eqn */
      /**** T = vertically-integrated longitudinal stress ****/
      ul = (i == 1) ? user->ua : Hu[i-1].u;
      Fl = GetFSR(dx,user->epsilon,user->n, ul,Hu[i].u);
      Tl = Bstag[i-1] * (Hu[i-1].H + Hu[i].H) * Fl;
      Fr = GetFSR(dx,user->epsilon,user->n, Hu[i].u,Hu[i+1].u);
      Tr = Bstag[i] * (Hu[i].H + Hu[i+1].H) * Fr;
      /**** sliding coefficient ****/
      beta = user->k * rg * Hu[i].H * GLREG(Hu[i].H,Hg,0.0);
      /**** gl-regularized surface slope ****/
      Hl = (i == 1) ? user->Ha : Hu[i-1].H;
      tmp = (1.0 - omega) * Hl + 0.0 - user->zocean;
      hl = omega * Hl + user->zocean + GLREG(Hl,Hg,0.0) * tmp;
      Hr = Hu[i+1].H;
      tmp = (1.0 - omega) * Hr + 0.0 - user->zocean;
      hr = omega * Hr + user->zocean + GLREG(Hr,Hg,0.0) * tmp;
      dhdx = (hr - hl) / (2.0 * dx);
      /**** residual for SSA ****/
      f[i].u = (Tr - Tl) / dx - beta * Hu[i].u - rg * Hu[i].H * dhdx;
    }
  }
  ierr = DMDAVecRestoreArray(user->stagda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->stagda,locMstag,&Mstag);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(user->stagda,&locBstag);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(user->stagda,&locMstag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* Apply scalings to variables and equations to make it all possible! */
PetscErrorCode scshell(DMDALocalInfo *info, Node *Hu, Node *f, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt i,
           start = PetscMax(info->xs - 1,            0),
           end   = PetscMin(info->xs + info->xm + 1, user->Mx);
  /* variable scaling coeffs */
  PetscReal scaleH = 1000.0,
            scaleu = 100.0 / user->secpera;

  /* residual scaling coeffs */
  PetscReal rscH      = 1.0 / scaleH,
            rscu      = 1.0 / scaleu,
            rscuH     = 1.0 / (scaleH * scaleu),
            rscstress = 1.0 / (user->rho * user->g);
  /* dimensionalize unknowns (put in "real" scale), including the ghost values */
  for (i = start; i < end; i++) {
    Hu[i].H *= scaleH;
    Hu[i].u *= scaleu;
  }
  /* compute residual in dimensional units */
  ierr = FunctionLocalGLREG(info, Hu, f, user); CHKERRQ(ierr);
  /* scale the residual to be O(1) */
  for (i = info->xs; i < info->xs + info->xm; i++) {
    if (i == 0) {
      f[0].H *= rscH;
      f[0].u *= rscu;
    } else {
      f[i].H *= rscuH;
      f[i].u *= rscstress;
    }
  }
  /* de-dimensionalize unknowns */
  for (i = start; i < end; i++) {
    Hu[i].H /= scaleH;
    Hu[i].u /= scaleu;
  }
  return 0;
}



#if 0
static inline PetscReal dFSRdleft(PetscReal dx, PetscReal eps, PetscReal n,
                                    PetscReal ul, PetscReal ur) {
  PetscReal dudx = (ur - ul) / dx,
              q    = (1.0 / n) - 1.0,
              D2   = dudx * dudx + eps * eps;
  return - (1.0 / dx) * PetscPowScalar(D2, (q / 2.0) - 1) * ( q * dudx * dudx + D2 );
}


static inline PetscReal dFSRdright(PetscReal dx, PetscReal eps, PetscReal n,
                                    PetscReal ul, PetscReal ur) {
  return - dFSRdleft(dx,eps,n,ul,ur);
}


/* Evaluate analytical Jacobian matrix. */
/* FIXME:  this is not fully implemented */
PetscErrorCode JacobianMatrixLocal(DMDALocalInfo *info, Node *Hu, Mat jac, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      v[6];
  MatStencil     row,col[6];

  PetscFunctionBegin;
  for (i=info->xs; i<info->xs+info->xm; i++) {

    /* MASS CONT */
    row.i = i; row.c = 0;
    if (i == 0) {
      col[0].i = i; col[0].c = 0;   v[0] = scentry(user,row,col[0],1.0);
      ierr = MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
    } else {
      if (user->upwind1) {
        /* 1st-order upwind */
        col[0].i = i; col[0].c = 0;   v[0] = scentry(user,row,col[0], - Hu[i].u);
        col[1].i = i; col[1].c = 1;   v[1] = scentry(user,row,col[1], - Hu[i].H);
        if (i == 1) {
          ierr = MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);      
        } else {
          col[2].i = i-1; col[2].c = 0;   v[2] = scentry(user,row,col[2], + Hu[i-1].u);
          col[3].i = i-1; col[3].c = 1;   v[3] = scentry(user,row,col[3], + Hu[i-1].H);
          ierr = MatSetValuesStencil(jac,1,&row,4,col,v,INSERT_VALUES);CHKERRQ(ierr);      
        }
      } else {
        /* 2nd-order upwind */
        if (i == 1) {
          col[0].i = i; col[0].c = 0;   v[0] = scentry(user,row,col[0], - 2.0 * Hu[i].u);
          col[1].i = i; col[1].c = 1;   v[1] = scentry(user,row,col[1], - 2.0 * Hu[i].H);
          ierr = MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);      
        } else {
          col[0].i = i;   col[0].c = 0;   v[0] = scentry(user,row,col[0], - 1.5 * Hu[i].u);
          col[1].i = i;   col[1].c = 1;   v[1] = scentry(user,row,col[1], - 1.5 * Hu[i].H);
          col[2].i = i-1; col[2].c = 0;   v[2] = scentry(user,row,col[2], + 2.0 * Hu[i-1].u);
          col[3].i = i-1; col[3].c = 1;   v[3] = scentry(user,row,col[3], + 2.0 * Hu[i-1].H);
          if (i == 2) {
            ierr = MatSetValuesStencil(jac,1,&row,4,col,v,INSERT_VALUES);CHKERRQ(ierr);
          } else {
            col[4].i = i-2; col[4].c = 0;   v[4] = scentry(user,row,col[4], - 0.5 * Hu[i-2].u);
            col[5].i = i-2; col[5].c = 1;   v[5] = scentry(user,row,col[5], - 0.5 * Hu[i-2].H);
            ierr = MatSetValuesStencil(jac,1,&row,6,col,v,INSERT_VALUES);CHKERRQ(ierr);
          }
        }
      }
    }  /* done with MASS CONT */

    /* SSA */
    row.i = i; row.c = 1;
    if (i == 0) {
      col[0].i = i; col[0].c = 1;   v[0] = scentry(user,row,col[0],1.0);
      ierr = MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
    } else {
      /* FIXME implement */
    }
  }

  /* assemble matrix, using the 2-step process */
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  /* tell matrix we will never add a new nonzero location; if we do then gives error  */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif
