/************************************************************************
 *                                                                      *
 * File       : cvfkx.c                                                 *
 * Programmers: Scott D. Cohen and Alan C. Hindmarsh and                *
 *              Radu Serban @ LLNL                                      *
 * Version of : 22 July 2003                                            *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * An ODE system is generated from the following 2-species diurnal      *
 * kinetics advection-diffusion PDE system in 2 space dimensions:       *
 *                                                                      *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)    *
 *                 + Ri(c1,c2,t)      for i = 1,2,   where              *
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,       *
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,                    *
 *   Kv(z) = Kv0*exp(z/5) ,                                             *
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)        *
 * vary diurnally.   The problem is posed on the square                 *
 *   0 <= x <= 20,    30 <= z <= 50   (all in km),                      *
 * with homogeneous Neumann boundary conditions, and for time t in      *
 *   0 <= t <= 86400 sec (1 day).                                       *
 * The PDE system is treated by central differences on a uniform        *
 * 10 x 10 mesh, with simple polynomial initial profiles.               *
 * The problem is solved with CVODES, with the BDF/GMRES method (i.e.   *
 * using the CVSPGMR linear solver) and the block-diagonal part of the  *
 * Newton matrix as a left preconditioner. A copy of the block-diagonal *
 * part of the Jacobian is saved and conditionally reused within the    *
 * Precond routine.                                                     *
 *                                                                      *
 * Optionally, CVODES can compute sensitivities with respect to the     *
 * problem parameters q1 and q2.                                        *
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and       *
 * STAGGERED1) can be used and sensitivities may be included in the     *
 * error test or not (error control set on FULL or PARTIAL,             *
 * respectively).                                                       *
 *                                                                      *
 * Execution:                                                           *
 *                                                                      *
 * If no sensitivities are desired:                                     *
 *    % cvskx -nosensi                                                  *
 * If sensitivities are to be computed:                                 *
 *    % cvskx -sensi sensi_meth err_con                                 *
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of    *
 * {full, partial}.                                                     *
 *                                                                      *
 ************************************************************************/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sundialstypes.h"     /* definitions of realtype, integertype            */
#include "cvodes.h"            /* main CVODES header file                         */
#include "cvodec.h"
#include "iterativ.h"          /* contains the enum for types of preconditioning  */
#include "cvsspgmr.h"          /* use CVSPGMR linear solver each internal step    */
#include "smalldense.h"        /* use generic DENSE solver for preconditioning    */
#include "nvector_serial.h"    /* definitions of type N_Vector, macro NV_DATA_S   */
#include "sundialsmath.h"      /* contains SQR macro                              */

/* Problem Constants */

#define NUM_SPECIES  2             /* number of species */
#define C1_SCALE     1.0e6         /* coefficients in initial profiles */
#define C2_SCALE     1.0e12

#define T0           0.0           /* initial time */
#define NOUT         12            /* number of output times */
#define TWOHR        7200.0        /* number of seconds in two hours  */
#define HALFDAY      4.32e4        /* number of seconds in a half day */
#define PI       3.1415926535898   /* pi */ 

#define XMIN          0.0          /* grid boundaries in x  */
#define XMAX         20.0           
#define ZMIN         30.0          /* grid boundaries in z  */
#define ZMAX         50.0
#define XMID         10.0          /* grid midpoints in x,z */          
#define ZMID         40.0

#define MX           15             /* MX = number of x mesh points */
#define MZ           15             /* MZ = number of z mesh points */
#define NSMX         NUM_SPECIES*MX /* NSMX = NUM_SPECIES*MX */
#define MM           (MX*MZ)        /* MM = MX*MZ */

/* CVodeMalloc Constants */
#define RTOL    1.0e-5            /* scalar relative tolerance */
#define FLOOR   100.0             /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */
#define NEQ     (NUM_SPECIES*MM)  /* NEQ = number of equations */

/* Sensitivity Constants */
#define NP    8
#define NS    2

#define ZERO  RCONST(0.0)

/* User-defined vector and matrix accessor macros: IJKth, IJth */

/* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MZ-1. The vdata array is obtained via
   the macro call vdata = NV_DATA_S(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])


/* Type : UserData 
   contains preconditioner blocks, pivot arrays, 
   problem parameters, and problem constants     */

typedef struct {
  realtype *p;
  realtype **P[MX][MZ], **Jbd[MX][MZ];
  integertype *pivot[MX][MZ];
  realtype q4, om, dx, dz, hdco, haco, vdco;
} *UserData;


/* Private Helper Functions */

static void WrongArgs(char *argv[]);
static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector y, realtype dx, realtype dz);
static void PrintOutput(void *cvode_mem, realtype t, N_Vector y);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(void *cvode_mem, booleantype sensi);

/* Functions Called by the CVODES Solver */

extern void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static int Precond(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int PSolve(realtype tn, N_Vector y, N_Vector fy,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *P_data, N_Vector vtemp);


static UserData AllocUserDataIm(void);
static void InitUserDataIm(UserData data);
static void FreeUserDataIm(UserData data);

extern void f_C_Cplx(realtype t_re,realtype t_im,
                     N_Vector y_Cplx_Re_,N_Vector y_Cplx_Im_,
                     N_Vector ydot_Cplx_Re_,N_Vector ydot_Cplx_Im_,
                     UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  int ia;
  NV_Spec nvSpec;
  void *cvode_mem;
  UserData data;
  realtype abstol, reltol, t, tout;
  N_Vector y;
  int iout, flag;

  realtype *pbar;
  int is, *plist;
  N_Vector *uS=NULL;
  booleantype sensi=FALSE;
  int sensi_meth=-1, err_con=-1;

  booleantype CSderivs;
  UserData data_im;

  sensi = FALSE;
  CSderivs = FALSE;

  /* Process arguments */


  ia = 1;
  while(ia < argc) {

    if (strcmp(argv[ia],"-x") == 0) {

      CSderivs = TRUE;
      ia++;

    }

    else if (strcmp(argv[ia],"-s") == 0) {

      if (argc < ia+3) WrongArgs(argv);

      sensi = TRUE;
      ia++;

      if (strcmp(argv[ia],"sim") == 0)
        sensi_meth = SIMULTANEOUS;
      else if (strcmp(argv[ia],"stg") == 0)
        sensi_meth = STAGGERED;
      else if (strcmp(argv[ia],"stg1") == 0)
        sensi_meth = STAGGERED1;
      else 
        WrongArgs(argv);

      ia++;

      if (strcmp(argv[ia],"full") == 0)
        err_con = FULL;
      else if (strcmp(argv[ia],"partial") == 0)
        err_con = PARTIAL;
      else
        WrongArgs(argv);
      
      ia++;

    }

    else 

      WrongArgs(argv);

  }


  /* Initialize serial vector specification */
  nvSpec = NV_SpecInit_Serial(NEQ);

  /* PROBLEM PARAMETERS */
  data = AllocUserData();
  InitUserData(data);

  /* INITIAL STATES */
  y = N_VNew(nvSpec);
  SetInitialProfiles(y, data->dx, data->dz);
  
  /* TOLERANCES */
  abstol=ATOL; 
  reltol=RTOL;

  /* CVODE_CREATE */
  cvode_mem = CVodeCreate(BDF, NEWTON);
  if (cvode_mem == NULL) { printf("CVodeCreate failed.\n"); return(1); }

  flag = CVodeSetFdata(cvode_mem, data);
  if (flag != SUCCESS) { printf("CVodeSetFdata failed.\n"); return(1); }

  flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
  if (flag != SUCCESS) { printf("CVodeSetMaxNumSteps failed.\n"); return(1); }

  if (CSderivs) {
    data_im = AllocUserDataIm();
    InitUserDataIm(data_im);
    flag = CVodeSetCSDerivs(cvode_mem, f_C_Cplx, data_im);
  }

  /* CVODE_MALLOC */
  flag = CVodeMalloc(cvode_mem, f, T0, y, SS, &reltol, &abstol, nvSpec);
  if (flag != SUCCESS) { printf("CVodeMalloc failed.\n"); return(1); }

  /* CVSPGMR */
  flag = CVSpgmr(cvode_mem, LEFT, 0);
  if (flag != SUCCESS) { printf("CVSpgmr failed."); return(1); }

  if (CSderivs)
    flag = CVSpgmrSetCSJacTimesVec(cvode_mem);

  flag = CVSpgmrSetPrecSetupFn(cvode_mem, Precond);
  flag = CVSpgmrSetPrecSolveFn(cvode_mem, PSolve);
  flag = CVSpgmrSetPrecData(cvode_mem, data);

  /* SENSITIVTY */
  if(sensi) {
    pbar = (realtype *) malloc(NP*sizeof(realtype));
    for(is=0; is<NP; is++) pbar[is] = data->p[is];
    plist = (int *) malloc(NS * sizeof(int));
    for(is=0; is<NS; is++) plist[is] = is+1;

    uS = N_VNew_S(NS, nvSpec);
    for(is=0;is<NS;is++)
      N_VConst(ZERO,uS[is]);

    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    flag = CVodeSetSensRho(cvode_mem, ZERO);
    flag = CVodeSetSensPbar(cvode_mem, pbar);

    if (CSderivs)
      flag = CVodeSetSensCSRhs(cvode_mem, data_im->p);

    flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, data->p, plist, uS);
    if (flag != SUCCESS) {printf("CVodeSensMalloc failed, flag=%d\n",flag);return(1);}
    
  }

  /* In loop over output points, call CVode, print results, test for error */
  
  printf("\n2-species diurnal advection-diffusion problem\n\n");

  printf("========================================================================\n");
  printf("     T     Q       H      NST                    Bottom left  Top right \n");
  printf("========================================================================\n");

  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    if (flag != SUCCESS) { 
      printf("CVode failed, flag=%d.\n", flag); 
      break; 
    }
    PrintOutput(cvode_mem, t, y);
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, t, uS);
      if (flag != OKAY) { printf("CVodeGetSens failed, flag=%d.\n", flag); break; }
      PrintOutputS(uS);
    }
    
    printf("------------------------------------------------------------------------\n");
    
  }

  /* Print final statistics */
  printf("\n\n========================================================");
  printf("\nFinal Statistics");
  printf("\nDerivative Calculation: ");
  if(CSderivs) printf("complex step");
  else         printf("finite differences");
  printf("\nSensitivity: ");
  if(sensi) {
    printf("YES ");
    if(sensi_meth == SIMULTANEOUS)   
      printf("( SIMULTANEOUS +");
    else 
      if(sensi_meth == STAGGERED) printf("( STAGGERED +");
      else                        printf("( STAGGERED1 +");   
    if(err_con == FULL) printf(" FULL ERROR CONTROL )");
    else                printf(" PARTIAL ERROR CONTROL )");
  } else {
    printf("NO");
  }

  PrintFinalStats(cvode_mem, sensi);
  printf("========================================================\n");

  /* Free memory */
  N_VFree(y);
  if(sensi) N_VFree_S(NS, uS);
  FreeUserData(data);
  if (CSderivs)
    FreeUserDataIm(data_im);
  CVodeFree(cvode_mem);
  NV_SpecFree_Serial(nvSpec);

  return(0);
}


/*********************** Private Helper Functions ************************/

/* ======================================================================= */
/* Exit if arguments are incorrect */

static void WrongArgs(char *argv[])
{
  printf("\nUsage: %s [OPTIONS]\n",argv[0]);
  printf("            -x          = use xs4c\n");
  printf("            -s meth err = perform SA\n");
  printf("               meth : sim, stg, or stg1\n");
  printf("               err  : full or partial\n");
  
  exit(0);
}

/* ======================================================================= */
/* Allocate memory for data structure of type UserData */

static UserData AllocUserData(void)
{
  int jx, jz;
  UserData data;

  data = (UserData) malloc(sizeof *data);

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      (data->P)[jx][jz] = denalloc(NUM_SPECIES);
      (data->Jbd)[jx][jz] = denalloc(NUM_SPECIES);
      (data->pivot)[jx][jz] = denallocpiv(NUM_SPECIES);
    }
  }

  data->p = (realtype *) malloc(NP*sizeof(realtype));

  return(data);
}

static UserData AllocUserDataIm(void)
{
  UserData data;

  data = (UserData) malloc(sizeof *data);

  data->p = (realtype *) malloc(NP*sizeof(realtype));

  return(data);
}

/* ======================================================================= */
/* Free data memory */

static void FreeUserData(UserData data)
{
  int jx, jz;

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      denfree((data->P)[jx][jz]);
      denfree((data->Jbd)[jx][jz]);
      denfreepiv((data->pivot)[jx][jz]);
    }
  }

  free(data->p);

  free(data);
}

static void FreeUserDataIm(UserData data)
{

  free(data->p);

  free(data);
}

/* ======================================================================= */
/* Load problem constants in data */

static void InitUserData(UserData data)
{
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Set problem parameters */
  Q1 = 1.63e-16; /* Q1  coefficients q1, q2, c3             */
  Q2 = 4.66e-16; /* Q2                                      */
  C3 = 3.7e16;   /* C3                                      */
  A3 = 22.62;    /* A3  coefficient in expression for q3(t) */
  A4 = 7.601;    /* A4  coefficient in expression for q4(t) */
  KH = 4.0e-6;   /* KH  horizontal diffusivity Kh           */ 
  VEL = 0.001;   /* VEL advection velocity V                */
  KV0 = 1.0e-8;  /* KV0 coefficient in Kv(z)                */  

  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/(MX-1);
  data->dz = (ZMAX-ZMIN)/(MZ-1);
  data->hdco = KH/SQR(data->dx);
  data->haco = VEL/(2.0*data->dx);
  data->vdco = (1.0/SQR(data->dz))*KV0;

  data->p[0] = Q1;
  data->p[1] = Q2;
  data->p[2] = C3;
  data->p[3] = A3;
  data->p[4] = A4;
  data->p[5] = KH;
  data->p[6] = VEL;
  data->p[7] = KV0;
}

static void InitUserDataIm(UserData data)
{
  int i;

  data->om   = 0.0;
  data->dx   = 0.0;
  data->dz   = 0.0; 
  data->hdco = 0.0; 
  data->haco = 0.0; 
  data->vdco = 0.0; 

  for(i=0; i<8; i++)
    data->p[i] = 0.0; 
}

/* ======================================================================= */
/* Set initial conditions in y */

static void SetInitialProfiles(N_Vector y, realtype dx, realtype dz)
{
  int jx, jz;
  realtype x, z, cx, cz;
  realtype *ydata;

  /* Set pointer to data array in vector y. */

  ydata = NV_DATA_S(y);

  /* Load initial profiles of c1 and c2 into y vector */

  for (jz=0; jz < MZ; jz++) {
    z = ZMIN + jz*dz;
    cz = SQR(0.1*(z - ZMID));
    cz = 1.0 - cz + 0.5*SQR(cz);
    for (jx=0; jx < MX; jx++) {
      x = XMIN + jx*dx;
      cx = SQR(0.1*(x - XMID));
      cx = 1.0 - cx + 0.5*SQR(cx);
      IJKth(ydata,1,jx,jz) = C1_SCALE*cx*cz; 
      IJKth(ydata,2,jx,jz) = C2_SCALE*cx*cz;
    }
  }
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector y)
{  
  int nst, qu;
  realtype hu;
  realtype *ydata;

  ydata = NV_DATA_S(y);
  
  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetLastOrder(cvode_mem, &qu);
  CVodeGetLastStep(cvode_mem, &hu);

  printf("%8.3e %2d  %8.3e %5d\n", t,qu,hu,nst);
  printf("                                Solution       ");
  printf("%12.4e %12.4e \n", IJKth(ydata,1,0,0), IJKth(ydata,1,MX-1,MZ-1)); 
  printf("                                               ");
  printf("%12.4e %12.4e \n", IJKth(ydata,2,0,0), IJKth(ydata,2,MX-1,MZ-1));
  
}

/* ======================================================================= */
/* Print sampled sensitivities */

static void PrintOutputS(N_Vector *uS)
{
  realtype *sdata;

  sdata = NV_DATA_S(uS[0]);

  printf("                                ----------------------------------------\n"); 
  printf("                                Sensitivity 1  ");
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
  printf("                                               ");
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));

  sdata = NV_DATA_S(uS[1]);

  printf("                                ----------------------------------------\n"); 
  printf("                                Sensitivity 2  ");
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
  printf("                                               ");
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));

}

/* ======================================================================= */
/* Print final statistics contained in iopt */

static void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
  int nst;
  int nfe, nsetups, nni, ncfn, netf;
  int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int nli, ncfl, npe, nps;

  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetNumRhsEvals(cvode_mem, &nfe);
  CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  CVodeGetNumErrTestFails(cvode_mem, &netf);
  CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  if (sensi) {
    CVodeGetNumSensRhsEvals(cvode_mem, &nfSe);
    CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    CVodeGetNumSensLinSolvSetups(cvode_mem, &nsetupsS);
    CVodeGetNumSensErrTestFails(cvode_mem, &netfS);
    CVodeGetNumSensNonlinSolvIters(cvode_mem, &nniS);
    CVodeGetNumSensNonlinSolvConvFails(cvode_mem, &ncfnS);
  }

  CVSpgmrGetNumLinIters(cvode_mem, &nli);
  CVSpgmrGetNumConvFails(cvode_mem, &ncfl);
  CVSpgmrGetNumPrecEvals(cvode_mem, &npe);
  CVSpgmrGetNumPrecSolves(cvode_mem, &nps);

  printf("\n\n");
  printf("nst     = %5d\n\n", nst);
  printf("nfe     = %5d\n",   nfe);
  printf("netf    = %5d    nsetups  = %5d\n", netf, nsetups);
  printf("nni     = %5d    ncfn     = %5d\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5d    nfeS     = %5d\n", nfSe, nfeS);
    printf("netfs   = %5d    nsetupsS = %5d\n", netfS, nsetupsS);
    printf("nniS    = %5d    ncfnS    = %5d\n", nniS, ncfnS);
  }

  printf("\n");
  printf("nli     = %5d    ncfl     = %5d\n", nli, ncfl);
  printf("npe     = %5d    nps      = %5d\n", npe, nps);
  

}

/* ======================================================================= */
/* Preconditioner setup routine. Generate and preprocess P. */

static int Precond(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype c1, c2, czdn, czup, diag, zdn, zup, q4coef, delz, verdco, hordco;
  realtype **(*P)[MZ], **(*Jbd)[MZ];
  integertype *(*pivot)[MZ];
  int ier, jx, jz;
  realtype *ydata, **a, **j;
  UserData data;
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Make local copies of pointers in P_data, and of pointer to y's data */
  data = (UserData) P_data;
  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  ydata = NV_DATA_S(y);

  /* Load problem coefficients and parameters */
  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];
  A3 = data->p[3];
  A4 = data->p[4];
  KH = data->p[5];
  VEL = data->p[6];
  KV0 = data->p[7];

  if (jok) {

  /* jok = TRUE: Copy Jbd to P */

    for (jz=0; jz < MZ; jz++)
      for (jx=0; jx < MX; jx++)
        dencopy(Jbd[jx][jz], P[jx][jz], NUM_SPECIES);

  *jcurPtr = FALSE;

  }

  else {
  /* jok = FALSE: Generate Jbd from scratch and copy to P */

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  delz = data->dz;
  verdco = data->vdco;
  hordco  = data->hdco;

  /* Compute 2x2 diagonal Jacobian blocks (using q4 values 
     computed on the last f call).  Load into P. */

    for (jz=0; jz < MZ; jz++) {
      zdn = ZMIN + (jz - .5)*delz;
      zup = zdn + delz;
      czdn = verdco*exp(0.2*zdn);
      czup = verdco*exp(0.2*zup);
      diag = -(czdn + czup + 2.0*hordco);
      for (jx=0; jx < MX; jx++) {
        c1 = IJKth(ydata,1,jx,jz);
        c2 = IJKth(ydata,2,jx,jz);
        j = Jbd[jx][jz];
        a = P[jx][jz];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + q4coef;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4coef) + diag;
        dencopy(j, a, NUM_SPECIES);
      }
    }

  *jcurPtr = TRUE;

  }

  /* Scale by -gamma */

    for (jz=0; jz < MZ; jz++)
      for (jx=0; jx < MX; jx++)
        denscale(-gamma, P[jx][jz], NUM_SPECIES);

  /* Add identity matrix and do LU decompositions on blocks in place. */

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      denaddI(P[jx][jz], NUM_SPECIES);
      ier = gefa(P[jx][jz], NUM_SPECIES, pivot[jx][jz]);
      if (ier != 0) return(1);
    }
  }

  return(0);
}


/* ======================================================================= */
/* Preconditioner solve routine */

static int PSolve(realtype tn, N_Vector y, N_Vector fy,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *P_data, N_Vector vtemp)
{
  realtype **(*P)[MZ];
  integertype *(*pivot)[MZ];
  int jx, jz;
  realtype *zdata, *v;
  UserData data;

  /* Extract the P and pivot arrays from P_data. */

  data = (UserData) P_data;
  P = data->P;
  pivot = data->pivot;
  zdata = NV_DATA_S(z);

  N_VScale(1.0, r, z);

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      v = &(IJKth(zdata, 1, jx, jz));
      gesl(P[jx][jz], NUM_SPECIES, pivot[jx][jz], v);
    }
  }

  return(0);
}
 
