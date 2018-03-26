/*
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Type : UserData */

typedef struct {
  realtype l, alpha, epsilon, lambda, z, r0;
} *UserData;

static UserData InitUserData(void);

static void PrintOutput(realtype r, realtype p, realtype q);

static void PrintFinalStats(void *cvode_mem);

static int check_flag(void *flagvalue, char *funcname, int opt);

static int f(realtype r, N_Vector y, N_Vector ydot, void *f_data);

static int Jac(long int N, DenseMat J, realtype r,
               N_Vector y, N_Vector fy, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  UserData data;

  FILE *fout;

  long int Neq;
  void *cvode_mem;

  N_Vector y, abstol;
  realtype reltol, r, rout;

  realtype r0, l;

  int flag;


  y = abstol = NULL;
  cvode_mem = NULL;

  Neq = 2;

  data = InitUserData();

  /* Create serial vector of length NEQ for I.C. and abstol */

  y = N_VNew_Serial(Neq);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  abstol = N_VNew_Serial(Neq); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

  /* Initialize y */

  r0 = data->r0;
  l = data->l;

  Ith(y,1) = RPowerR(r0,l+1);
  //Ith(y,2) = 1.0/RPowerR(r0,l);
  Ith(y,2) = RPowerR(r0,l);

  printf("initial y:\n");
  N_VPrint_Serial(y);

  /* Set the scalar relative tolerance */

  reltol = 1.0e-6;

  /* Set the vector absolute tolerance */

  Ith(abstol,1) = 1.0e-10;
  Ith(abstol,2) = 1.0e-10;

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  
  flag = CVodeSetFdata(cvode_mem, data);
  if (check_flag((void *)cvode_mem, "CVodeSetFdata", 0)) return(1);

  flag = CVodeSetInitStep(cvode_mem, 1.0e-10);

  flag = CVodeMalloc(cvode_mem, f, r0, y, CV_SV, reltol, abstol);
  if (check_flag(&flag, "CVodeMalloc", 1)) return(1);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, Neq);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDenseSetJacFn(cvode_mem, Jac, data);
  if (check_flag(&flag, "CVDenseSetJacFn", 1)) return(1);

  fout = fopen("dirac.out","w");

  r = r0;
  rout = 10.0;

  while(r < rout) {

    flag = CVode(cvode_mem, rout, y, &r, CV_ONE_STEP);
    PrintOutput(r, Ith(y,1), Ith(y,2));

    fprintf(fout,"%g  %g  %g\n",r, Ith(y,1), Ith(y,2));

    if (check_flag(&flag, "CVode", 1)) break;
  }

  fclose(fout);

  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  free(data);


  return(0);
}


static UserData InitUserData(void)
{
  int jx, jy;
  UserData data;

  data = (UserData) malloc(sizeof *data);

  data->l = 1.0;
  data->alpha = 1.0/137.0;
  data->epsilon = 10.0;
  data->lambda = 3.0;
  data->z = 40.0;
  data->r0 = 1.0e-5;

  return(data);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y). 
 */

static int f(realtype r, N_Vector y, N_Vector ydot, void *f_data)
{
  UserData data;
  realtype p, q, pdot, qdot;
  realtype l, k, alpha, alpha2, epsilon, lambda, z, r0, cntrf, V;

  data = (UserData) f_data;

  //  printf("--------------------------\n");

  //  printf("r = %g \n",r);

  l = data->l;
  alpha = data->alpha;
  epsilon = data->epsilon;
  lambda = data->lambda;
  z = data->z;
  r0 = data->r0;
  k = l*(l+1);
  alpha2 = alpha*alpha;

  cntrf = k/r;
  V = -z/r * exp(-lambda*r);

  //  printf("cntrf = %g  V = %g\n",cntrf,V);

  p = Ith(y,1); 
  q = Ith(y,2);

  //  printf("p = %g  q = %g\n",p,q);
  //  printf("-cntrf*p = %g  alpha2*(epsilon-V) = %g (2.0 + alpha2*(epsilon-V)) = %g    (2.0 + alpha2*(epsilon-V)) * q = %g\n",
  //         -cntrf*p, alpha2*(epsilon-V), 2.0 + alpha2*(epsilon-V), (2.0 + alpha2*(epsilon-V)) * q);
  //  printf("pdot = %g\n",-cntrf * p + (2.0 + alpha2*(epsilon-V)) * q);

  pdot = -cntrf * p + (2.0 + alpha2*(epsilon-V)) * q;

  //  printf("pdot = %g\n",pdot);

  //  printf("epsilon-V = %g   (epsilon-V)*p = %g   p*cntrf = %g\n",epsilon-V,(epsilon-V)*p,p*cntrf);

  qdot = -(epsilon-V) * p + cntrf * q;

  //  printf("pdot = %g  qdot = %g\n",pdot,qdot);

  Ith(ydot,1) = pdot;
  Ith(ydot,2) = qdot;

  //  printf("y vector:\n");
  //  //  N_VPrint_Serial(y);
  //  printf("yd vector:\n");
  //  N_VPrint_Serial(ydot);

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(long int N, DenseMat J, realtype r,
               N_Vector y, N_Vector fy, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data;
  realtype p, q;
  realtype l, k, alpha, alpha2, epsilon, lambda, z, r0, cntrf, V;

  data = (UserData) jac_data;

  //  printf("==================================\n");

  //  printf("r = %g \n",r);

  l = data->l;
  alpha = data->alpha;
  epsilon = data->epsilon;
  lambda = data->lambda;
  z = data->z;
  r0 = data->r0;

  k = l*(l+1);
  alpha2 = alpha*alpha;

  cntrf = k/r;
  V = -z/r * exp(-lambda*r);

  p = Ith(y,1);
  q = Ith(y,2);

  //  printf("p = %g  q = %g\n",p,q);

  IJth(J,1,1) = -cntrf;
  IJth(J,1,2) = (2 + alpha2*(epsilon-V));
  IJth(J,2,1) = -(epsilon-V);
  IJth(J,2,2) = cntrf;


  //  DensePrint(J);

  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(realtype r, realtype p, realtype q)
{
  printf("At r = %0.4le      p = %14.6le  q = %14.6le\n", r, p, q);

  return;
}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDenseGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDenseGetNumJacEvals", 1);
  flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDenseGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
