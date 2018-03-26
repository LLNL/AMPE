#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_lapack.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>




#define LMM CP_ADAMS
#define NLS CP_FUNCTIONAL
#define TF 10.0



#define Ith(v,i)    NV_Ith_S(v,i-1)

/* Named constants */

#define ZERO   RCONST(0.0)
#define PT5    RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO    RCONST(2.0)
#define FOUR   RCONST(4.0)

/* Problem Constants */

#define RTOL  RCONST(1.0e-6)
#define ATOL  RCONST(1.0e-8)

typedef struct {
  realtype eps;
  realtype c;
  realtype e0;
} *UserData;


/* Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data);
static void PrintFinalStats(void *cpode_mem);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  UserData d;
  void *cpode_mem;
  N_Vector yy, yp, cc;
  realtype reltol, abstol, t, tout, Tout;
  realtype x, y, xd, yd, g;
  realtype pi;
  int iout, Nperiods, flag;

  pi = FOUR * atan(ONE);
  
  d = (UserData) malloc(sizeof(*d));
  d->eps = RCONST(0.01);
  d->c   = RCONST(0.6);
  d->e0 = -0.5 - 0.5*d->eps/(1-d->c)/(1-d->c)/(1-d->c);

  yy = N_VNew_Serial(4);
  yp = N_VNew_Serial(4);

  cc  = N_VNew_Serial(1);

  /* Initialize y */
  Ith(yy,1) = ONE - d->c;
  Ith(yy,2) = ZERO;
  Ith(yy,3) = ZERO;
  Ith(yy,4) = RSqrt((ONE+d->c)/(ONE-d->c));


  /* Check constraint at init. time */
  flag = cfun(ZERO, yy, cc, d);

  N_VPrint_Serial(cc);


  /* Set tolerances */
  reltol = RTOL; 
  abstol = ATOL;

  /* Initialize solver */
  cpode_mem = CPodeCreate(LMM, NLS);  
  flag = CPodeSetUserData(cpode_mem, d);
  flag = CPodeInitExpl(cpode_mem, (void *)f, ZERO, yy);
  flag = CPodeSStolerances(cpode_mem, reltol, abstol);
  flag = CPodeSetMaxNumSteps(cpode_mem, 50000);

  if (NLS == CP_NEWTON)
    flag = CPLapackDense(cpode_mem, 4);

  /* INTERNAL PROJECTION FUNCTION */  
  Ith(cc,1) = 1.0e-8;
  flag = CPodeProjInit(cpode_mem, CP_PROJ_L2NORM, CP_CNSTR_NONLIN, cfun, cc);
  flag = CPodeSetProjTestCnstr(cpode_mem, TRUE);
  flag = CPLapackDenseProj(cpode_mem, 1, 4, CPDIRECT_QR);

  /* Disable projection? */
  //CPodeSetProjFrequency(cpode_mem, 0);

 
  /* INTEGRATE THROUGH A SEQUENCE OF TIMES */

  Tout = TF;
  t = ZERO;
  while (t < Tout) {

    flag = CPode(cpode_mem, Tout, &t, yy, yp, CP_ONE_STEP);
    if (flag < 0) break;

    x  = Ith(yy,1);
    y  = Ith(yy,2);
    xd = Ith(yy,3);
    yd = Ith(yy,4);

    flag = cfun(t, yy, cc, d);
    g = Ith(cc,1);
    /*
    printf("%g  %14.10g  %14.10g  %14.10g  %14.10g    %14.10g    %d\n",
           t, x,y,xd,yd,g, flag);
    */
  }

  PrintFinalStats(cpode_mem);

  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(cc);
  CPodeFree(&cpode_mem);

  free(d);

  return(0);
}


static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  UserData d;
  realtype q1, q2, v1, v2, r, rr, rr2, rr3;

  d = (UserData)f_data;

  q1  = Ith(yy,1);
  q2  = Ith(yy,2);
  v1 = Ith(yy,3);
  v2 = Ith(yy,4);

  r = RSqrt(q1*q1+q2*q2);
  rr  = ONE/r;
  rr2 = rr * rr;
  rr3 = rr * rr2;

  Ith(fy,1) = v1;
  Ith(fy,2) = v2;
  Ith(fy,3) = - q1 * rr3 * (ONE + ONEPT5 * d->eps * rr2);
  Ith(fy,4) = - q2 * rr3 * (ONE + ONEPT5 * d->eps * rr2);

  return(0);
}

static int cfun(realtype t, N_Vector yy, N_Vector cout, void *c_data)
{
  UserData d;
  realtype q1, q2, v1, v2;
  realtype r, rr, rr2, rr3;
  realtype eps, e0;

  d = (UserData)c_data;
  eps = d->eps;
  e0 = d->e0;

  q1  = Ith(yy,1);
  q2  = Ith(yy,2);
  v1 = Ith(yy,3);
  v2 = Ith(yy,4);

  r = RSqrt(q1*q1+q2*q2);
  rr  = ONE/r;
  rr2 = rr * rr;
  rr3 = rr * rr2;

  Ith(cout,1) = PT5 * (v1*v1+v2*v2) - rr - PT5 * eps * rr3 - e0;

  return(0);
}


static void PrintFinalStats(void *cpode_mem)
{
  realtype h0u;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  long int nproj, nce, nsetupsP, nprf;
  int flag;

  flag = CPodeGetActualInitStep(cpode_mem, &h0u);
  flag = CPodeGetNumSteps(cpode_mem, &nst);
  flag = CPodeGetNumFctEvals(cpode_mem, &nfe);
  flag = CPodeGetNumLinSolvSetups(cpode_mem, &nsetups);
  flag = CPodeGetNumErrTestFails(cpode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cpode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  if (NLS == CP_NEWTON) {
    flag = CPDlsGetNumJacEvals(cpode_mem, &nje);
    flag = CPDlsGetNumFctEvals(cpode_mem, &nfeLS);
  } else {
    nje = 0;
    nfeLS = 0;
  }

  flag = CPodeGetProjStats(cpode_mem, &nproj, &nce, &nsetupsP, &nprf);

  flag = CPodeGetNumGEvals(cpode_mem, &nge);

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld\n",
	 nst, nfe, nsetups);
  printf("nfeLS = %-6ld nje = %ld\n",
	 nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld \n",
	 nni, ncfn, netf);
  printf("nproj = %-6ld nce = %-6ld nsetupsP = %-6ld nprf = %-6ld\n",
         nproj, nce, nsetupsP, nprf);
  printf("nge = %ld\n", nge);
}

