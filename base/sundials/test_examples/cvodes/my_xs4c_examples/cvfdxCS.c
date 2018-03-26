#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sundialstypes.h"   /* definitions of types realtype (set to double) and */
                             /* integertype (set to int), and the constant FALSE  */
#include "cvodes.h"          /* prototypes for CVodeMalloc, CVode, and CVodeFree, */
                             /* constants OPT_SIZE, BDF, NEWTON, SV, SUCCESS,     */
                             /* NST, NFE, NSETUPS, NNI, NCFN, NETF                */
#include "cvodec.h"
#include "cvsdense.h"        /* prototype for CVDense, constant DENSE_NJE         */
#include "nvector_serial.h"  /* definitions of type N_Vector and macro NV_Ith_S,  */
                             /* prototypes for N_VNew, N_VFree                    */
#include "dense.h"           /* definitions of type DenseMat, macro DENSE_ELEM    */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */
#define NEQ   3            /* number of equations  */
#define Y1    1.0          /* initial y components */
#define Y2    0.0
#define Y3    0.0
#define RTOL  1e-4         /* scalar relative tolerance            */
#define ATOL1 1e-8         /* vector absolute tolerance components */
#define ATOL2 1e-14
#define ATOL3 1e-6
#define T0    0.0          /* initial time           */
#define T1    0.4          /* first output time      */
#define TMULT 10.0         /* output time factor     */
#define NOUT  12           /* number of output times */

#define NP    3
#define NS    3

#define ZERO  0.0

/* Type : UserData */
typedef struct {
  realtype p[3];
} *UserData;

/* Private Helper Function */

static void WrongArgs(char *argv[]);
static void PrintFinalStats(void *cvode_mem, booleantype sensi);
static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
static void PrintOutputS(N_Vector *uS);

/* Functions Called by the CVODES Solver */

extern void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
extern void f_C_Cplx(realtype t_re,realtype t_im,
                     N_Vector y_Cplx_Re_,N_Vector y_Cplx_Im_,
                     N_Vector ydot_Cplx_Re_,N_Vector ydot_Cplx_Im_,
                     UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);

static void Jac(integertype N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
               int iS, N_Vector yS, N_Vector ySdot, 
               void *fS_data, N_Vector tmp1, N_Vector tmp2);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{

  int ia;

  NV_Spec nvSpec;
  void *cvode_mem;
  UserData data;
  realtype reltol, t, tout;
  N_Vector y, abstol;
  int iout, flag;

  realtype pbar[NP];
  int is, *plist; 
  N_Vector *yS=NULL;
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

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data);
  data->p[0] = 0.04;
  data->p[1] = 1.0e4;
  data->p[2] = 3.0e7;

  /* INITIAL STATES */
  y = N_VNew(nvSpec);
  abstol = N_VNew(nvSpec);

  /* Initialize y */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;
  Ith(y,3) = Y3;

  /* TOLERANCES */
  /* Set the scalar relative tolerance */
  reltol = RTOL;               
  /* Set the vector absolute tolerance */
  Ith(abstol,1) = ATOL1;       
  Ith(abstol,2) = ATOL2;
  Ith(abstol,3) = ATOL3;

  /* CVODE_CREATE */
  cvode_mem = CVodeCreate(BDF, NEWTON);

  flag = CVodeSetFdata(cvode_mem, data);


  if (CSderivs) {
    data_im = (UserData) malloc(sizeof *data);
    data_im->p[0] = 0.0;
    data_im->p[1] = 0.0;
    data_im->p[2] = 0.0;
    
    flag = CVodeSetCSDerivs(cvode_mem, f_C_Cplx, data_im);
  }

  /* CVODE_MALLOC */
  flag = CVodeMalloc(cvode_mem, f, T0, y, SV, &reltol, abstol, nvSpec);

  /* CVDENSE */
  flag = CVDense(cvode_mem, NEQ);

  if (CSderivs) {

    flag = CVDenseSetCSJac(cvode_mem);

  } else {

    flag = CVDenseSetJacFn(cvode_mem, Jac);
    flag = CVDenseSetJacData(cvode_mem, data);

  }


  /* SENSITIVITY */
  if(sensi) {
    pbar[0] = data->p[0];
    pbar[1] = data->p[1];
    pbar[2] = data->p[2];
    plist = (int *) malloc(NS * sizeof(int));
    for(is=0;is<NS;is++) plist[is] = is+1;

    yS = N_VNew_S(NS, nvSpec);
    for(is=0;is<NS;is++)
      N_VConst(0.0, yS[is]);

    if (CSderivs) {
      flag = CVodeSetSensCSRhs(cvode_mem, data_im->p);
    } else {
      flag = CVodeSetSensRhs1Fn(cvode_mem, fS);
      flag = CVodeSetSensFdata(cvode_mem, data);
    }

    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    flag = CVodeSetSensPbar(cvode_mem, pbar);

    flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, data->p, plist, yS);
    if (flag != SUCCESS) { printf("CVodeSensMalloc failed, flag=%d\n",flag); return(1); }
  }
  
  /* In loop over output points, call CVode, print results, test for error */

  printf("\n3-species chemical kinetics problem\n\n");
  printf("===================================================");
  printf("==================================\n");
  printf("     T     Q       H      NST                    y1");
  printf("           y2           y3    \n");
  printf("===================================================");
  printf("==================================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    if (flag != SUCCESS) {
      printf("CVode failed, flag=%d.\n", flag); 
      break; 
    }
    PrintOutput(cvode_mem, t, y);
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, t, yS);
      if (flag != OKAY) { printf("CVodeGetSens failed, flag=%d.\n", flag); break; }
      PrintOutputS(yS);
    } 
    printf("-------------------------------------------------");
    printf("------------------------------------\n");
  }

  /* Print final statistics */
  printf("\n\n========================================================");
  printf("\nFinal Statistics");
  printf("\nDerivative Calculation: ");
  if(CSderivs) printf("complex step");
  else         printf("analytical");
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
  N_VFree(y);                  /* Free the y and abstol vectors       */
  N_VFree(abstol);   
  if(sensi) N_VFree_S(NS, yS); /* Free the yS vectors                 */
  free(data);                  /* Free user data                      */
  if (CSderivs)
    free(data_im);
  CVodeFree(cvode_mem);        /* Free the CVODES problem memory      */
  NV_SpecFree_Serial(nvSpec);  /* Free the vector specification       */

  return(0);
}


/************************ Private Helper Function ************************/

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
/* Print current t, step count, order, stepsize, and solution  */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u)
{
  int nst, qu;
  realtype hu, *udata;
  
  udata = NV_DATA_S(u);

  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetLastOrder(cvode_mem, &qu);
  CVodeGetLastStep(cvode_mem, &hu);

  printf("%8.3e %2d  %8.3e %5d\n", t, qu, hu, nst);
  printf("                                Solution       ");
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
  
}
/* ======================================================================= */
/* Print sensitivities */

static void PrintOutputS(N_Vector *uS)
{

  realtype *sdata;

  sdata = NV_DATA_S(uS[0]);
  printf("                                Sensitivity 1  ");
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
  
  sdata = NV_DATA_S(uS[1]);
  printf("                                Sensitivity 2  ");
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);

  sdata = NV_DATA_S(uS[2]);
  printf("                                Sensitivity 3  ");
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);

}

/* ======================================================================= */
/* Print some final statistics located in the iopt array */

static void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
  int nst;
  int nfe, nsetups, nni, ncfn, netf;
  int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int njeD, nfeD;

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

  CVDenseGetNumJacEvals(cvode_mem, &njeD);
  CVDenseGetNumRhsEvals(cvode_mem, &nfeD);


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
  printf("njeD    = %5d    nfeD     = %5d\n", njeD, nfeD);


}


/***************** Functions Called by the CVODES Solver ******************/

/* ======================================================================= */
/* Jacobian routine. Compute J(t,y). */

static void Jac(integertype N, DenseMat J, realtype t,
                N_Vector y, N_Vector fy, void *jac_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2, p3;
 
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) jac_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
 
  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
                      IJth(J,3,2) = 2*p3*y2;

}
 
/* ======================================================================= */
/* fS routine. Compute sensitivity r.h.s. */

static void fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
               int iS, N_Vector yS, N_Vector ySdot, 
               void *fS_data, N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype sd1, sd2, sd3;

  data = (UserData) fS_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  y1 = Ith(y,1);  y2 = Ith(y,2);  y3 = Ith(y,3);
  s1 = Ith(yS,1); s2 = Ith(yS,2); s3 = Ith(yS,3);

  sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
  sd3 = 2*p3*y2*s2;
  sd2 = -sd1-sd3;

  switch (iS) {
  case 0:
    sd1 += -y1;
    sd2 +=  y1;
    break;
  case 1:
    sd1 +=  y2*y3;
    sd2 += -y2*y3;
    break;
  case 2:
    sd2 += -y2*y2;
    sd3 +=  y2*y2;
    break;
  }
  
  Ith(ySdot,1) = sd1;
  Ith(ySdot,2) = sd2;
  Ith(ySdot,3) = sd3;

}
