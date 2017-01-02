static const char help[] = 
"Solve Bodvardsson equations (Bueler interpretation) using SNES, a dof=2 Vec\n"
"holding both thickness H and velocity u, and a 2nd-order finite difference\n"
"scheme to approximate the coupled mass continuity and SSA stress balance PDEs.\n"
"\n"
"This run shows success with finite-difference-Jacobian:\n"
"  ./bodvardsson -snes_fd -da_grid_x 181 -snes_monitor\n"
"This runs show success with matrix-free, no preconditioner:\n"
"  ./bodvardsson -snes_mf -da_grid_x 181 -snes_monitor -bod_exact_init\n"
"Visualization:\n"
"  ./bodvardsson -snes_fd -da_grid_x 201 -snes_monitor -snes_monitor_solution -draw_pause 1\n"
"See convbod.sh for convergence graph.  Add option -bod_up_one to see first-order upwinding.\n"
"Compare marine.c.\n\n";

#include <petscdmda.h>
#include <petscsnes.h>

#include "exactsolns.h"


/* we will use dof=2 DMDA, and at each grid point have a thickness H and a velocity u */
typedef struct {
  PetscReal H, u;
} Node;


/* User-defined application context.  Used especially by BodFunctionLocal().  */
typedef struct {
  DM          da;      /* 1d,dof=2 distributed array for soln and residual */
  DM          scalarda;/* 1d,dof=1 distributed array for suitable for parameters depending on x */
  PetscMPIInt rank;
  PetscInt    Mx, xs, xm, solnghostwidth;
  PetscBool   upwind1; /* if true, used low-order upwinding */
  PetscReal   dx, secpera, n, rho, rhow, g;
  PetscReal   H0;      /* thickness at x=0, for Dirichlet condition on mass cont */
  PetscReal   xc;      /* location at which stress (Neumann) condition applied to SSA eqn */
  PetscReal   Txc;     /* vertically-integrated longitudinal stress at xc, for Neumann cond:
                            T = 2 H B |u_x|^{(1/n)-1} u_x  */
  PetscReal   epsilon; /* regularization of viscosity, a strain rate */
  PetscReal   scaleNode[2]; /* scaling seen "inside" SNES is not what our function sees */
  Vec         Huexact; /* exact thickness (Huexact[i][0]) and exact velocity
                          (Huexact[i][1]) on regular grid */
  Vec         M;       /* surface mass balance on regular grid */
  Vec         beta;    /* sliding coefficient on regular grid */
  Vec         Bstag;   /* ice hardness on staggered grid*/
} AppCtx;


extern PetscErrorCode FillExactSoln(AppCtx*);
extern PetscErrorCode FillInitial(AppCtx*,Vec*);
extern PetscErrorCode FillDistributedParams(AppCtx*);
extern PetscErrorCode BodFunctionLocal(DMDALocalInfo*,Node*,Node*,AppCtx*);
extern PetscErrorCode scshell(DMDALocalInfo*,Node*,Node*,AppCtx*);


int main(int argc,char **argv)
{
  PetscErrorCode         ierr;

  SNES                   snes;                 /* nonlinear solver */
  Vec                    Hu,r;                 /* solution, residual vectors */
  AppCtx                 user;                 /* user-defined work context */
  PetscInt               i, tmpxs, tmpxm; /* iteration count, index, etc. */
  PetscReal              tmp1, tmp2, tmp3, tmp4,
                         errnorms[2], descaleNode[2];
  PetscBool              eps_set = PETSC_FALSE,
                         dump = PETSC_FALSE,
                         exactinitial = PETSC_FALSE;

  PetscInitialize(&argc,&argv,(char *)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &user.rank); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "BODVARDSSON solves for thickness and velocity in 1D, steady ice stream\n"
    "  [run with -help for info and options]\n");CHKERRQ(ierr);

  user.n       = 3.0;          /* Glen flow law exponent */
  user.secpera = 31556926.0;
  user.rho     = 910.0;        /* kg m^-3 */
  user.rhow    = 1028.0;       /* kg m^-3 */
  user.g       = 9.81;         /* m s^-2 */

  /* ask Bodvardsson soln for its parameters, but only those we need to solve */
  ierr = params_exactBod(&(user.H0), &tmp1, &(user.xc), &tmp2, &tmp3, &tmp4); CHKERRQ(ierr);
  /* regularize using strain rate of 1/xc per year */
  user.epsilon = (1.0 / user.secpera) / user.xc;
  /* tools for non-dimensionalizing to improve equation scaling */
  user.scaleNode[0] = 1000.0;  user.scaleNode[1] = 100.0 / user.secpera;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"bod_",
      "bodvardsson program options","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsBool("-up_one","","",PETSC_FALSE,&user.upwind1,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-exact_init","","",PETSC_FALSE,&exactinitial,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-dump",
      "dump out exact and approximate solution and residual, as asci","",
      dump,&dump,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-epsilon","regularization (a strain rate in units of 1/a)","",
                            user.epsilon * user.secpera,&user.epsilon,&eps_set);CHKERRQ(ierr);
    if (eps_set)  user.epsilon *= 1.0 / user.secpera;
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* Create machinery for parallel grid management (DMDA), nonlinear solver (SNES),
     and Vecs for fields (solution, RHS).  Note default Mx=46 grid points means
     dx=10 km.  Also degrees of freedom = 2 (thickness and velocity
     at each point) and stencil radius = ghost width = 2 for 2nd-order upwinding.  */
  user.solnghostwidth = 2;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,46,2,user.solnghostwidth,PETSC_NULL,&user.da);
            CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da); CHKERRQ(ierr);
  ierr = DMSetUp(user.da); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(user.da,&user);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(user.da,0.0,user.xc,0.0,1.0,0.0,1.0);CHKERRQ(ierr);

  ierr = DMDASetFieldName(user.da,0,"ice thickness [non-dimensional]"); CHKERRQ(ierr);
  ierr = DMDASetFieldName(user.da,1,"ice velocity [non-dimensional]"); CHKERRQ(ierr);

  ierr = DMSetFromOptions(user.da); CHKERRQ(ierr);

  ierr = DMDAGetInfo(user.da,PETSC_IGNORE,&user.Mx,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                     PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DMDAGetCorners(user.da,&user.xs,PETSC_NULL,PETSC_NULL,&user.xm,PETSC_NULL,PETSC_NULL);
                   CHKERRQ(ierr);
  user.dx = user.xc / (PetscReal)(user.Mx-1);

  /* another DMDA for scalar parameters, with same length */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,user.Mx,1,1,PETSC_NULL,&user.scalarda);CHKERRQ(ierr);
  ierr = DMSetUp(user.scalarda); CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(user.scalarda,0.0,user.xc,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
  /* check that parallel layout of scalar DMDA is same as dof=2 DMDA */
  ierr = DMDAGetCorners(user.scalarda,&tmpxs,PETSC_NULL,PETSC_NULL,&tmpxm,PETSC_NULL,PETSC_NULL);
                   CHKERRQ(ierr);
  if ((tmpxs != user.xs) || (tmpxm != user.xm)) {
    PetscPrintf(PETSC_COMM_SELF,
       "\n***ERROR: rank %d gets different ownership range for the two DMDAs!  ENDING ...***\n\n",
       user.rank);
    PetscEnd();
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "  Mx = %D points, dx = %.3f m\n  H0 = %.2f m, xc = %.2f km, Txc = %.5e Pa m\n",
      user.Mx, user.dx, user.H0, user.xc/1000.0, user.Txc);CHKERRQ(ierr);

  /* Extract/allocate global vectors from DMDAs and duplicate for remaining same types */
  ierr = DMCreateGlobalVector(user.da,&Hu);CHKERRQ(ierr);
  ierr = VecSetBlockSize(Hu,2);CHKERRQ(ierr);
  ierr = VecDuplicate(Hu,&r);CHKERRQ(ierr); /* inherits block size */
  ierr = VecDuplicate(Hu,&user.Huexact);CHKERRQ(ierr); /* ditto */

  ierr = DMCreateGlobalVector(user.scalarda,&user.M);CHKERRQ(ierr);
  ierr = VecDuplicate(user.M,&user.Bstag);CHKERRQ(ierr);
  ierr = VecDuplicate(user.M,&user.beta);CHKERRQ(ierr);

  /* set up snes */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,user.da);CHKERRQ(ierr);

  ierr = DMDASNESSetFunctionLocal(user.da,INSERT_VALUES,(DMDASNESFunction)scshell,&user);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* the Bodvardsson (1955) exact solution allows setting M(x), B(x), beta(x), T(xc) */
  ierr = FillDistributedParams(&user);CHKERRQ(ierr);

  /* the exact thickness and exact ice velocity (user.uHexact) are known from Bodvardsson (1955) */
  ierr = FillExactSoln(&user); CHKERRQ(ierr);

  if (exactinitial) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  using exact solution as initial guess\n");
             CHKERRQ(ierr);
    /* the initial guess is the exact continuum solution */
    ierr = VecCopy(user.Huexact,Hu); CHKERRQ(ierr);
  } else {
    ierr = FillInitial(&user, &Hu); CHKERRQ(ierr);
  }

  /************ SOLVE NONLINEAR SYSTEM  ************/
  /* recall that RHS  r  is used internally by KSP, and is set by the SNES */
  for (i = 0; i < 2; i++)  descaleNode[i] = 1.0 / user.scaleNode[i];
  ierr = VecStrideScaleAll(Hu,descaleNode); CHKERRQ(ierr); /* de-dimensionalize initial guess */
  ierr = SNESSolve(snes,PETSC_NULL,Hu);CHKERRQ(ierr);
  ierr = VecStrideScaleAll(Hu,user.scaleNode); CHKERRQ(ierr); /* put back in "real" scale */

  if (dump) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing combined result Hu\n");CHKERRQ(ierr);
    ierr = VecView(Hu,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing combined exact result Huexact\n");CHKERRQ(ierr);
    ierr = VecView(user.Huexact,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing final combined residual at Hu\n");CHKERRQ(ierr);
    ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  }

  /* evaluate error relative to exact solution */
  ierr = VecAXPY(Hu,-1.0,user.Huexact); CHKERRQ(ierr);  /* Hu = - Huexact + Hu */
  ierr = VecStrideNormAll(Hu,NORM_INFINITY,errnorms); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "(dx,errHinf,erruinf) %.3f %.4e %.4e\n",
           user.dx,errnorms[0],errnorms[1]*user.secpera);CHKERRQ(ierr);

  ierr = VecDestroy(&Hu);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = VecDestroy(&(user.Huexact));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.M));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.Bstag));CHKERRQ(ierr);
  ierr = VecDestroy(&(user.beta));CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  ierr = DMDestroy(&(user.da));CHKERRQ(ierr);
  ierr = DMDestroy(&(user.scalarda));CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

/*  Compute the exact thickness and velocity on the regular grid. */
PetscErrorCode FillExactSoln(AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      x, dum;
  Node           *Hu;

  PetscFunctionBegin;
  /* Compute regular grid exact soln and staggered-grid thickness over the
     locally-owned part of the grid */
  ierr = DMDAVecGetArray(user->da,user->Huexact,&Hu);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    x = user->dx * (PetscReal)i;  /* = x_i = distance from dome */
    ierr = exactBod(x, &(Hu[i].H), &(Hu[i].u), &dum); CHKERRQ(ierr);
  }
  ierr = DMDAVecRestoreArray(user->da,user->Huexact,&Hu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*  Put a not-unreasonable initial guess in Hu. */
PetscErrorCode FillInitial(AppCtx *user, Vec *vHu)
{
  PetscErrorCode ierr;
  PetscInt       i;
  Node           *Hu;

  PetscFunctionBegin;
  ierr = DMDAVecGetArray(user->da,*vHu,&Hu);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    Hu[i].H = 1000.0;
    Hu[i].u = 100.0 / user->secpera;
  }
  ierr = DMDAVecRestoreArray(user->da,*vHu,&Hu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*  Compute the values of the surface mass balance, ice hardness, and sliding coeff. */
PetscErrorCode FillDistributedParams(AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i, Mx=user->Mx;
  PetscReal      x, H, k, dum1, dum2, dum3, dum4, dum5, *M, *Bstag, *beta;

  PetscFunctionBegin;
  ierr = params_exactBod(&dum1, &dum2, &dum3, &dum4, &dum5, &k); CHKERRQ(ierr);
  /* Compute regular grid exact soln and staggered-grid thickness over the
     locally-owned part of the grid */
  ierr = DMDAVecGetArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->scalarda,user->Bstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    x = user->dx * (PetscReal)i;  /* = x_i = distance from dome; regular grid point */
    ierr = exactBod(x, &H, &dum1, &(M[i])); CHKERRQ(ierr);
    beta[i] = k * user->rho * user->g * H; 
    x = x + (user->dx/2.0);       /* = x_{i+1/2}; staggered grid point */
    if (i < Mx-1) {
      ierr = exactBodBueler(x, &dum1, &(Bstag[i])); CHKERRQ(ierr);
    } else {
      Bstag[i] = -9999.9999;  /* never accessed, and we don't have the value anyway */
    }
  }
  ierr = DMDAVecRestoreArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->scalarda,user->Bstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);
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


PetscErrorCode BodFunctionLocal(DMDALocalInfo *info, Node *Hu, Node *f, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscReal      rg = user->rho * user->g, dx = user->dx;
  PetscReal      *M, *Bstag, *beta,
                 duH, ul, u, ur, dHdx, Fl, Fr, Tl, Tr;
  PetscInt       i, Mx = info->mx;
  Vec            locBstag;

  PetscFunctionBegin;

  /* we need stencil width on Bstag (but not for M, beta) */
  ierr = DMGetLocalVector(user->scalarda,&locBstag);CHKERRQ(ierr);  /* do NOT destroy it */
  ierr = DMGlobalToLocalBegin(user->scalarda,user->Bstag,INSERT_VALUES,locBstag); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->scalarda,user->Bstag,INSERT_VALUES,locBstag); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(user->scalarda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);
  for (i = info->xs; i < info->xs + info->xm; i++) {
  
    /* MASS CONT */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].H = Hu[0].H - user->H0;
    } else {
      /* centered difference IS UNSTABLE:
        duH = Hu[i+1].u * Hu[i+1].H - ( (i == 1) ? 0.0 : Hu[i-1].u * Hu[i-1].H );
        duH *= 0.5; */

      if (user->upwind1) {
        /* 1st-order upwind; leftward difference because u > 0 (because dH/dx < 0) */
        duH = Hu[i].u * Hu[i].H - ( (i == 1) ? 0.0 : Hu[i-1].u * Hu[i-1].H );
      } else {
        /* 2nd-order upwind; see Beam-Warming discussion in R. LeVeque, "Finite Volume ..." */
        if (i == 1) { /* use PDE  M - (uH)_x = 0 to get quadratic poly, then diff that */
          duH = - dx * M[0] + 2.0 * Hu[1].u * Hu[1].H;
        } else {
          duH = 3.0 * Hu[i].u * Hu[i].H - 4.0 * Hu[i-1].u * Hu[i-1].H;
          /* if i == 2 then u=0 so uH = 0 */
          if (i >= 3)  duH += Hu[i-2].u * Hu[i-2].H;
          duH *= 0.5;
        }
      }

      f[i].H = dx * M[i] - duH;
    }

    /* SSA */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].u = Hu[0].u - 0.0;
    } else {
      /* residual: SSA eqn */
      /* consecutive values of u */
      ul = (i == 1) ? 0.0 : Hu[i-1].u;
      u  = Hu[i].u;
      ur = (i == Mx-1) ? -1.1e30 : Hu[i+1].u;
      /* surface slope */
      if (i == 1) { 
        dHdx  = (Hu[i+1].H - user->H0) / (2.0 * dx);
      } else if (i == Mx-1) {
        /* nearly 2nd-order global convergence seems to occur even with this:
             dHdx  = (Hu[i].H - Hu[i-1].H) / dx; */
        dHdx  = (3.0*Hu[i].H - 4.0*Hu[i-1].H + Hu[i-2].H) / (2.0 * dx);
      } else { /* generic case */
        dHdx  = (Hu[i+1].H - Hu[i-1].H) / (2.0 * dx);
      }
      /* vertically-integrated longitudinal stress */
      Fl = GetFSR(dx,user->epsilon,user->n, ul,u);
      if (i == Mx-1) {
        Tl = 2.0 * (Hu[i-1].H + Hu[i].H) * Bstag[i-1] * Fl;
        Tr = (1.0 - user->rho / user->rhow) * user->rho * user->g * Hu[i].H * Hu[i].H;
        /* exact value: Tr = 2.0 * user->Txc; */
      } else {
        Fr = GetFSR(dx,user->epsilon,user->n, u,ur);
        Tl = (Hu[i-1].H + Hu[i].H) * Bstag[i-1] * Fl;
        Tr = (Hu[i].H + Hu[i+1].H) * Bstag[i] * Fr;        
      }

      f[i].u = (Tr - Tl) - dx * beta[i] * u - dx * rg * Hu[i].H * dHdx; /* SSA */

    }
  }
  ierr = DMDAVecRestoreArray(user->scalarda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(user->scalarda,&locBstag);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/* Apply scalings to variables and equations to improve. */
PetscErrorCode scshell(DMDALocalInfo *info, Node *Hu, Node *f, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt i,
           start = PetscMax(info->xs - user->solnghostwidth,            0),
           end   = PetscMin(info->xs + info->xm + user->solnghostwidth, user->Mx);

  /* residual scaling coeffs */
  PetscReal rscH      = 1.0 / user->H0,
            rscu      = user->secpera / 100.0,
            rscuH     = user->secpera / user->H0,
            rscstress = 1.0 / (user->rho * user->g * user->H0 * user->dx * 0.001);
  /* dimensionalize unknowns (put in "real" scale), including the ghost values */
  for (i = start; i < end; i++) {
    Hu[i].H *= (user->scaleNode)[0];  Hu[i].u *= (user->scaleNode)[1];
  }
  /* compute residual in dimensional units */
  ierr = BodFunctionLocal(info, Hu, f, user); CHKERRQ(ierr);
  /* scale the residual to be O(1) */
  for (i = info->xs; i < info->xs + info->xm; i++) {
    if (i == 0) {
      f[0].H *= rscH;   f[0].u *= rscu;
    } else {
      f[i].H *= rscuH;  f[i].u *= rscstress; }
  }
  /* de-dimensionalize unknowns */
  for (i = start; i < end; i++) {
    Hu[i].H /= (user->scaleNode)[0];  Hu[i].u /= (user->scaleNode)[1];
  }
  return 0;
}

