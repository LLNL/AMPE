/*
 * -----------------------------------------------------------------
 * $Revision: 1.27 $
 * $Date: 2005/11/08 23:41:52 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *   Modification of the cvfdx.c example to better illustrate the
 *   use of  CVSensRhs1Fn.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

/* Problem Constants */

#define NEQ   3             /* number of equations  */
#define Y1    RCONST(1.0)   /* initial y components */
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1e-4)  /* scalar relative tolerance */
#define ATOL1 RCONST(1e-8)  /* vector absolute tolerance components */
#define ATOL2 RCONST(1e-14)
#define ATOL3 RCONST(1e-6)
#define T0    RCONST(0.0)   /* initial time */
#define T1    RCONST(0.4)   /* first output time */
#define TMULT RCONST(10.0)  /* output time factor */
#define NOUT  12            /* number of output times */

#define NP    3             /* number of problem parameters */
#define NS    2             /* number of sensitivities computed */

#define ZERO  RCONST(0.0)

/* Type : UserData */

typedef struct {
  realtype p1;
  realtype p2;
  realtype p3;
} *UserData;

/* Prototypes of functions by CVODES */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static int Jac(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
              int iS, N_Vector yS, N_Vector ySdot, 
              void *fS_data, N_Vector tmp1, N_Vector tmp2);

/* Prototypes of private functions */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
static void PrintOutputS(N_Vector *uS, int Ns);
static void PrintFinalStats(void *cvode_mem);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *cvode_mem;
  UserData data;
  realtype t, tout, rtol;
  N_Vector y, atol;
  int iout, flag;

  realtype pbar[NS];
  int is; 
  N_Vector *yS;
  booleantype err_con;
  int sensi_meth;

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  data->p1 = RCONST(0.04);
  data->p2 = RCONST(1.0e4);
  data->p3 = RCONST(3.0e7);

  /* Initial conditions */
  y = N_VNew_Serial(NEQ);
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;

  /* Relative tolerance */
  rtol = RTOL;

  /* Absolute tolerance vector */
  atol = N_VNew_Serial(NEQ);
  Ith(atol,1) = ATOL1;
  Ith(atol,2) = ATOL2;
  Ith(atol,3) = ATOL3;

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  /* Allocate space for CVODES */
  flag = CVodeMalloc(cvode_mem, f, T0, y, CV_SV, rtol, atol);

  /* Attach user data */
  flag = CVodeSetFdata(cvode_mem, data);

  /* Attach linear solver */
  flag = CVDense(cvode_mem, NEQ);
  flag = CVDenseSetJacFn(cvode_mem, Jac, data);

  printf("\n3-species chemical kinetics problem\n");

  /* Sensitivity-related settings */

  sensi_meth = CV_SIMULTANEOUS;
  err_con = TRUE;

  /* We will compute sensitivities with respect to p1 and p3 */
  /* pbar is needed to estimate appropriate absolute tolerances for 
     the sensitivity variables */
  pbar[0] = data->p1;
  pbar[1] = data->p3;

  /* Initial conditions for sensitivity variables */  
  yS = N_VCloneVectorArray_Serial(NS, y);
  for (is=0;is<NS;is++) N_VConst(ZERO, yS[is]);

  /* Initialize sensitivity computations */
  flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, yS);

  flag = CVodeSetSensRhs1Fn(cvode_mem, fS, data);
  flag = CVodeSetSensErrCon(cvode_mem, err_con);
  flag = CVodeSetSensParams(cvode_mem, NULL, pbar, NULL);
  
  /* In loop over output points, call CVode, print results, test for error */
  
  printf("\n\n");
  printf("===================================================");
  printf("============================\n");
  printf("     T     Q       H      NST                    y1");
  printf("           y2           y3    \n");
  printf("===================================================");
  printf("============================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    PrintOutput(cvode_mem, t, y);

    flag = CVodeGetSens(cvode_mem, t, yS);
    PrintOutputS(yS, NS);

    printf("-------------------------------------------------");
    printf("------------------------------\n");

  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem);

  /* Free memory */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(atol);
  N_VDestroyVectorArray_Serial(yS, NS);
  free(data);
  CVodeFree(&cvode_mem);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/* f routine. Compute f(t,y) */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y1, y2, y3, yd1, yd3;
  UserData data;
  realtype p1, p2, p3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) f_data;
  p1 = data->p1; p2 = data->p2; p3 = data->p3;

  yd1 = Ith(ydot,1) = -p1*y1 + p2*y2*y3;
  yd3 = Ith(ydot,3) = p3*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;

  return(0);
}


/* Jacobian routine. Compute J(t,y) */

static int Jac(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2, p3;
 
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) jac_data;
  p1 = data->p1; p2 = data->p2; p3 = data->p3;
 
  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
                      IJth(J,3,2) = 2*p3*y2;

  return(0);
}
 
/* fS routine. Compute sensitivity r.h.s. */

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
              int iS, N_Vector yS, N_Vector ySdot, 
              void *fS_data, N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype sd1, sd2, sd3;

  data = (UserData) fS_data;
  p1 = data->p1; p2 = data->p2; p3 = data->p3;

  y1 = Ith(y,1);  y2 = Ith(y,2);  y3 = Ith(y,3);
  s1 = Ith(yS,1); s2 = Ith(yS,2); s3 = Ith(yS,3);

  sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
  sd3 = 2*p3*y2*s2;
  sd2 = -sd1-sd3;

  switch (iS) {
  case 0:    /* sensitivity w.r.t. p1 */
    sd1 += -y1;
    sd2 +=  y1;
    break;
  case 1:    /* sensitivity w.r.t. p3 */
    sd2 += -y2*y2;
    sd3 +=  y2*y2;
    break;
  }
  
  Ith(ySdot,1) = sd1;
  Ith(ySdot,2) = sd2;
  Ith(ySdot,3) = sd3;

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* Print current t, step count, order, stepsize, and solution */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u)
{
  long int nst;
  int qu, flag;
  realtype hu, *udata;
  
  udata = NV_DATA_S(u);

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  flag = CVodeGetLastOrder(cvode_mem, &qu);
  flag = CVodeGetLastStep(cvode_mem, &hu);

  printf("%8.3le %2d  %8.3le %5ld\n", t, qu, hu, nst);

  printf("                          Solution       ");

  printf("%12.4le %12.4le %12.4le \n", udata[0], udata[1], udata[2]);

}

/* Print sensitivities */

static void PrintOutputS(N_Vector *uS, int Ns)
{
  realtype *sdata;
  int iS;

  for (iS=0; iS<Ns; iS++) {
    sdata = NV_DATA_S(uS[iS]);
    printf("                          Sensitivity %1d  ",iS+1);
    printf("%12.4le %12.4le %12.4le \n", sdata[0], sdata[1], sdata[2]);
  }  
}

/* Print some final statistics from the CVODES memory */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nje, nfeLS;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  flag = CVodeGetNumSensRhsEvals(cvode_mem, &nfSe);
  flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
  flag = CVodeGetNumSensLinSolvSetups(cvode_mem, &nsetupsS);
  flag = CVodeGetNumSensErrTestFails(cvode_mem, &netfS);
  flag = CVodeGetNumSensNonlinSolvIters(cvode_mem, &nniS);
  flag = CVodeGetNumSensNonlinSolvConvFails(cvode_mem, &ncfnS);

  flag = CVDenseGetNumJacEvals(cvode_mem, &nje);

  flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeLS);

  printf("\nFinal Statistics\n\n");
  printf("nst     = %5ld\n\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  printf("\n");
  printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
  printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
  printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);

  printf("\n");
  printf("nje    = %5ld    nfeLS     = %5ld\n", nje, nfeLS);

}

