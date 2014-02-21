static const char help[] =
"Solve steady flowline marine ice sheet equations using SNES, a dof=2 Vec\n"
"holding both thickness H and velocity u, and a 2nd-order finite difference\n"
"scheme to approximate the coupled mass continuity and SSA stress balance PDEs.\n"
"\n"
"TEST RUN DOES NOT WORK YET:\n"
"  ./marine -snes_fd -snes_monitor -snes_monitor_solution -draw_pause 1 -exact_init -ksp_diagonal_scale -snes_linesearch_monitor\n"
"FIXME This run shows success with finite-difference-Jacobian:\n"
"  ./marine -snes_fd -da_grid_x 181 -snes_monitor\n"
"FIXME Visualization:\n"
"  ./marine -snes_fd -da_grid_x 201 -snes_monitor -snes_monitor_solution -draw_pause 1\n\n";

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
/* extern PetscErrorCode JacobianMatrixLocal(DMDALocalInfo*,Node*,Mat,AppCtx*); */


int main(int argc,char **argv)
{
  PetscErrorCode         ierr;

  SNES                   snes;                 /* nonlinear solver */
  Vec                    Hu,r;                 /* solution, residual vectors */
  AppCtx                 user;                 /* user-defined work context */
  ExactCtx               exact;
  PetscInt               its, tmpxs, tmpxm; /* iteration count, index, etc. */
  PetscReal              tmp1, tmp2, tmp3,
                         errnorms[2];
  PetscBool              eps_set = PETSC_FALSE,
                         dump = PETSC_FALSE,
                         exactinitial = PETSC_FALSE,
                         snes_mf_set, snes_fd_set;
  SNESConvergedReason    reason;               /* Check convergence */

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

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,
           "","options to marine (steady marine ice sheet solver)","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsBool("-snes_mf","","",PETSC_FALSE,&snes_mf_set,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-snes_fd","","",PETSC_FALSE,&snes_fd_set,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-exact_init",
             "initialize using exact solution instead of default linear function","",
             PETSC_FALSE,&exactinitial,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-dump_solns",
             "dump out exact and approximate solution and residual, as asci","",
             dump,&dump,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-bod_epsilon","regularization (a strain rate in units of 1/a)","",
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

  /* Create machinery for parallel grid management (DMDA), nonlinear solver (SNES),
     and Vecs for fields (solution, RHS).  Note default grid points.  Also degrees
     of freedom = 2 (thickness and velocity at each point).  */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,-100,2,1,PETSC_NULL,&user.da);
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

  /* establish geometry on grid */
  /*user.dx = (user.xc-user.xa) / (PetscReal)(user.Mx-1);*/
  user.N = user.Mx - 2;
  user.dx = (user.xc - user.xa) / (PetscReal)(user.N+0.5);
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

  ierr = DMDASNESSetFunctionLocal(user.da,INSERT_VALUES,(DMDASNESFunction)FunctionLocalGLREG,&user);CHKERRQ(ierr);
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
  ierr = SNESSolve(snes,PETSC_NULL,Hu);CHKERRQ(ierr);

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
  ierr = VecDestroy(&(user.M));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.Bstag));CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  ierr = DMDestroy(&(user.da));CHKERRQ(ierr);
  ierr = DMDestroy(&(user.scalarda));CHKERRQ(ierr);

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
    if (x[i] < user->xa)
      SETERRQ2(PETSC_COMM_SELF,1,"ERROR on rank %d: x[i] = %f less than xa!",
               user->rank,x[i]);
    if (x[i] > user->xc)
      SETERRQ2(PETSC_COMM_SELF,2,"ERROR on rank %d: x[i] = %f greater than xc!",
               user->rank,x[i]);
    if (x[i] < exact->xg) {
      /* grounded part: Bodvardsson formula */
      ierr = exactBod(x[i], &(Hu[i].H), &(Hu[i].u), &dum1); CHKERRQ(ierr);
    } else {
      /* floating part: van der Veen formula */
      ierr = exactVeen(x[i], user->Mfloat, &(Hu[i].H), &(Hu[i].u)); CHKERRQ(ierr);
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
  PetscReal      dum1, dum2, xstag, *x, *M, *Bstag;
  DM             coord_da;
  Vec            coord_x;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDM(user->da, &coord_da); CHKERRQ(ierr);
  ierr = DMGetCoordinates(user->da, &coord_x); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(coord_da,coord_x,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->scalarda,user->Bstag,&Bstag);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    if (x[i] < user->xa)
      SETERRQ2(PETSC_COMM_SELF,1,"ERROR on rank %d: x[i] = %f less than xa!",
               user->rank,x[i]);
    if (x[i] > user->xc)
      SETERRQ2(PETSC_COMM_SELF,2,"ERROR on rank %d: x[i] = %f greater than xc!",
               user->rank,x[i]);
    if (x[i] < exact->xg) {
      ierr = exactBod(x[i], &dum1, &dum2, &(M[i])); CHKERRQ(ierr);
    } else  {
      M[i] = user->Mfloat;
    }
    xstag = x[i] + (user->dx/2.0);
    if (xstag < exact->xg) {
      ierr = exactBodBueler(xstag, &dum1, &(Bstag[i])); CHKERRQ(ierr);
    } else  {
      Bstag[i] = exact->Bg;
    }
  }
  ierr = DMDAVecRestoreArray(coord_da,coord_x,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->scalarda,user->Bstag,&Bstag);CHKERRQ(ierr);
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
  PetscReal tmp = exp(-(H-Hg)/eps);
  return tmp / (eps * (1.0 + tmp) * (1.0 + tmp));
}


/* Evaluate residual f at current iterate Hu, WITH ADDITIONAL REGULARIZATION
AT GROUNDING LINE.
For mass-continuity equation, note centered difference IS UNSTABLE:
   duH = Hu[i+1].u * Hu[i+1].H - ( (i == 1) ? 0.0 : Hu[i-1].u * Hu[i-1].H );
   duH *= 0.5;
*/
PetscErrorCode FunctionLocalGLREG(DMDALocalInfo *info, Node *Hu, Node *f, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscReal      rg = user->rho * user->g,
                 omega = user->omega,
                 dx = user->dx,
                 Hg = user->zocean * (user->rhow / user->rho),
                 epsH = 0.1 * Hg;
  PetscReal      *M, *Bstag,
                 duH, ul, Fl, Fr, Tl, Tr, dHdx, tmp, dhdx, beta;
  PetscInt       i, Mx = info->mx;
  Vec            locBstag;

  PetscFunctionBegin;

  /* we need stencil width on Bstag (but not for M, beta) */
  ierr = DMGetLocalVector(user->scalarda,&locBstag);CHKERRQ(ierr);  /* do NOT destroy it */
  ierr = DMGlobalToLocalBegin(user->scalarda,user->Bstag,INSERT_VALUES,locBstag); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->scalarda,user->Bstag,INSERT_VALUES,locBstag); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->scalarda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  for (i = info->xs; i < info->xs + info->xm; i++) {

    /* MASS CONT */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].H = user->ua * (Hu[0].H - user->Ha);  /* scale with u */
    } else {
      if (user->upwind1) {
        /* 1st-order upwind; leftward difference because u > 0 (because dH/dx < 0) */
        if (i == 1) {
          duH = Hu[i].u * Hu[i].H - user->ua * user->Ha;
        } else {
          duH = Hu[i].u * Hu[i].H - Hu[i-1].u * Hu[i-1].H;
        }
      } else {
        /* 2nd-order upwind; see Beam-Warming discussion in R. LeVeque, "Finite Volume ..." */
        if (i == 1) { /* first-order in this case */
          duH = Hu[i].u * Hu[i].H - user->ua * user->Ha;
        } else {
          if (i == 2) {
            duH = 3.0 * Hu[i].u * Hu[i].H - 4.0 * Hu[i-1].u * Hu[i-1].H + user->ua * user->Ha;
            duH *= 0.5;
          } else {
            duH = 3.0 * Hu[i].u * Hu[i].H - 4.0 * Hu[i-1].u * Hu[i-1].H + Hu[i-2].u * Hu[i-2].H;
            duH *= 0.5;
          }
        }
      }
      /**** residual for mass continuity ****/
      f[i].H = duH - dx * M[i];
    }

    /* SSA */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].u = Hu[0].u - user->ua;
    } else {
      /* residual: SSA eqn */
      /**** vertically-integrated longitudinal stress ****/
      ul = (i == 1) ? user->ua : Hu[i-1].u;
      Fl = GetFSR(dx,user->epsilon,user->n, ul,Hu[i].u);
      Tl = Bstag[i-1] * (Hu[i-1].H + Hu[i].H) * Fl;
      if (i == Mx-1) {
        /* here "staggered" T comes from calving-front boundary condition
           and introduced-point argument */
        Tr = 0.5 * user->omega * rg * Hu[i].H * Hu[i].H;
      } else {
        Fr = GetFSR(dx,user->epsilon,user->n, Hu[i].u,Hu[i+1].u);
        Tr = Bstag[i] * (Hu[i].H + Hu[i+1].H) * Fr;
      }
      /**** surface slope ****/
      if (i == 1) {
        dHdx  = (Hu[i+1].H - user->Ha) / (2.0 * dx);
      } else if (i == Mx-1) {
        /* nearly 2nd-order global convergence seems to occur even with this:
             dhdx  = (Hu[i].H - Hu[i-1].H) / dx; */
        dHdx  = (3.0*Hu[i].H - 4.0*Hu[i-1].H + Hu[i-2].H) / (2.0 * dx);
      } else { /* generic case */
        dHdx  = (Hu[i+1].H - Hu[i-1].H) / (2.0 * dx);
      }
      tmp = (1.0 - omega) * Hu[i].H + 0.0 - user->zocean;
      dhdx = dHdx * (omega + (1.0 - omega) * GLREG(Hu[i].H,Hg,epsH) + dGLREG(Hu[i].H,Hg,epsH) * tmp);
      /**** sliding coefficient ****/
      beta = user->k * rg * Hu[i].H * GLREG(Hu[i].H,Hg,epsH);
      /**** residual for SSA ****/
      if (i == Mx-1) {
        f[i].u = (2.0*Tr - 2.0*Tl) / dx;
      } else {
        f[i].u = (Tr - Tl) / dx;
      }
      f[i].u -= beta * Hu[i].u + rg * Hu[i].H * dhdx;
    }
  }
  ierr = DMDAVecRestoreArray(user->scalarda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(user->scalarda,&locBstag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
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
