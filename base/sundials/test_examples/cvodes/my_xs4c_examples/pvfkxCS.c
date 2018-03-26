/************************************************************************
 *                                                                      *
 * File       : pvfkx.c                                                 *
 * Programmers: S. D. Cohen, A. C. Hindmarsh, Radu Serban, and          *
 *              M. R. Wittman @ LLNL                                    *
 * Version of : 14 July 2003                                            *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * An ODE system is generated from the following 2-species diurnal      *
 * kinetics advection-diffusion PDE system in 2 space dimensions:       *
 *                                                                      *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)    *
 *                 + Ri(c1,c2,t)      for i = 1,2,   where              *
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,       *
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,                    *
 *   Kv(y) = Kv0*exp(y/5) ,                                             *
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)        *
 * vary diurnally.   The problem is posed on the square                 *
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),                      *
 * with homogeneous Neumann boundary conditions, and for time t in      *
 *   0 <= t <= 86400 sec (1 day).                                       *
 * The PDE system is treated by central differences on a uniform        *
 * mesh, with simple polynomial initial profiles.                       *
 *                                                                      *
 * The problem is solved by CVODES on NPE processors, treated as a      *
 * rectangular process grid of size NPEX by NPEY, with NPE = NPEX*NPEY. *
 * Each processor contains a subgrid of size MXSUB by MYSUB of the      *
 * (x,y) mesh.  Thus the actual mesh sizes are MX = MXSUB*NPEX and      *
 * MY = MYSUB*NPEY, and the ODE system size is neq = 2*MX*MY.           *
 *                                                                      *
 * The solution with CVODES is done with the BDF/GMRES method (i.e.     *
 * using the CVSPGMR linear solver) and the block-diagonal part of the  *
 * Newton matrix as a left preconditioner. A copy of the block-diagonal *
 * part of the Jacobian is saved and conditionally reused within the    *
 * Precond routine.                                                     *
 *                                                                      *
 * Performance data and sampled solution values are printed at selected *
 * output times, and all performance counters are printed on completion.*
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
 * NOTE: This version uses MPI for user routines, and the CVODES        *
 *       solver. In what follows, N is the number of processors,        *
 *       N = NPEX*NPEY (see constants below) and it is assumed that     *
 *       the MPI script mpirun is used to run a paralles application.   *
 * If no sensitivities are desired:                                     *
 *    % mpirun -np N pvfkx -nosensi                                     *
 * If sensitivities are to be computed:                                 *
 *    % mpirun -np N pvfkx -sensi sensi_meth err_con                    *
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of    *
 * {full, partial}.                                                     *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sundialstypes.h"    /* definitions of realtype, integertype            */
#include "cvodes.h"           /* main CVODES header file                         */
#include "cvodec.h"
#include "iterativ.h"         /* contains the enum for types of preconditioning  */
#include "cvsspgmr.h"         /* use CVSPGMR linear solver each internal step    */
#include "smalldense.h"       /* use generic DENSE solver in preconditioning     */
#include "nvector_parallel.h" /* definitions of type N_Vector, macro N_VDATA     */
#include "sundialsmath.h"     /* contains SQR macro                              */
#include "mpi.h"


/* Problem Constants */

#define NVARS        2             /* number of species                    */
#define C1_SCALE     1.0e6         /* coefficients in initial profiles     */
#define C2_SCALE     1.0e12

#define T0           0.0           /* initial time                         */
#define NOUT         12            /* number of output times               */
#define TWOHR        7200.0        /* number of seconds in two hours       */
#define HALFDAY      4.32e4        /* number of seconds in a half day      */
#define PI       3.1415926535898   /* pi                                   */ 

#define XMIN          0.0          /* grid boundaries in x                 */
#define XMAX         20.0           
#define YMIN         30.0          /* grid boundaries in y                 */
#define YMAX         50.0

#define NPEX         2              /* no. PEs in x direction of PE array  */
#define NPEY         2              /* no. PEs in y direction of PE array  */
                                    /* Total no. PEs = NPEX*NPEY           */
#define MXSUB        5              /* no. x points per subgrid            */
#define MYSUB        5              /* no. y points per subgrid            */

#define MX           (NPEX*MXSUB)   /* MX = number of x mesh points        */
#define MY           (NPEY*MYSUB)   /* MY = number of y mesh points        */
                                    /* Spatial mesh is MX by MY            */
/* CVodeMalloc Constants */

#define RTOL         1.0e-5         /* scalar relative tolerance            */
#define FLOOR        100.0          /* value of C1 or C2 at which tols.     */
                                    /* change from relative to absolute     */
#define ATOL         (RTOL*FLOOR)   /* scalar absolute tolerance            */

/* Sensitivity constants */
#define NP           8              /* number of problem parameters         */
#define NS           2              /* number of sensitivities              */

#define ZERO         RCONST(0.0)

/* User-defined matrix accessor macro: IJth */

/* IJth is defined in order to write code which indexes into small dense
   matrices with a (row,column) pair, where 1 <= row,column <= NVARS.   

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NVARS. The small matrix routines in dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */

#define IJth(a,i,j)        (a[j-1][i-1])

/* Type : UserData 
   contains problem constants, preconditioner blocks, pivot arrays, 
   grid constants, and processor indices */

typedef struct {
  realtype *p;
  realtype q4, om, dx, dy, hdco, haco, vdco;
  realtype uext[NVARS*(MXSUB+2)*(MYSUB+2)];
  integertype my_pe, isubx, isuby, nvmxsub, nvmxsub2;
  MPI_Comm comm;
} *UserData;

typedef struct {
  void *f_data;
  realtype **P[MXSUB][MYSUB], **Jbd[MXSUB][MYSUB];
  integertype *pivot[MXSUB][MYSUB];
} *PreconData;

/* Private Helper Functions */

static void WrongArgs(int my_pe, char *argv[]);
static PreconData AllocPreconData(UserData data);
static void InitUserData(int my_pe, MPI_Comm comm, UserData data);
static void FreePreconData(PreconData pdata);
static void SetInitialProfiles(N_Vector u, UserData data);
static void PrintOutput(void *cvode_mem, int my_pe, MPI_Comm comm,
                        realtype t, N_Vector u);
static void PrintOutputS(int my_pe, MPI_Comm comm, N_Vector *uS);
static void PrintFinalStats(void *cvode_mem, booleantype sensi);

/* Functions Called by the CVODES Solver */

extern void f(realtype t, N_Vector u, N_Vector udot, void *f_data);

static int Precond(realtype tn, N_Vector u, N_Vector fu,
                   booleantype jok, booleantype *jcurPtr, 
                   realtype gamma, void *P_data, 
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int PSolve(realtype tn, N_Vector u, N_Vector fu, 
                  N_Vector r, N_Vector z, 
                  realtype gamma, realtype delta,
                  int lr, void *P_data, N_Vector vtemp);


static void InitUserDataIm(UserData data);
void f_C_Cplx(realtype t_Cplx_Re_,realtype t_Cplx_Im_,N_Vector u_Cplx_Re_,
              N_Vector u_Cplx_Im_,N_Vector udot_Cplx_Re_,N_Vector udot_Cplx_Im_,
              UserData f_data_Cplx_Re_,UserData f_data_Cplx_Im_);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  int ia;
  NV_Spec nvSpec;
  realtype abstol, reltol, t, tout;
  N_Vector u;
  UserData data;
  PreconData predata;
  void *cvode_mem;
  int iout, flag, my_pe, npes;
  integertype neq, local_N;
  MPI_Comm comm;

  realtype *pbar;
  int is, *plist;
  N_Vector *uS=NULL;
  booleantype sensi=FALSE;
  int sensi_meth=-1, err_con=-1;

  booleantype CSderivs;
  UserData data_im;

  sensi = FALSE;
  CSderivs = FALSE;

  /* Set problem size neq */
  neq = NVARS*MX*MY;

  /* Get processor number and total number of pe's */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  if (npes != NPEX*NPEY) {
    if (my_pe == 0)
      printf("\n npes=%d is not equal to NPEX*NPEY=%d\n", npes,NPEX*NPEY);
    return(1);
  }

  /* Set local length */
  local_N = NVARS*MXSUB*MYSUB;

  /* Process arguments */

  ia = 1;
  while(ia < argc) {

    if (strcmp(argv[ia],"-x") == 0) {

      CSderivs = TRUE;
      ia++;

    }

    else if (strcmp(argv[ia],"-s") == 0) {

      if (argc < ia+3) WrongArgs(my_pe, argv);

      sensi = TRUE;
      ia++;

      if (strcmp(argv[ia],"sim") == 0)
        sensi_meth = SIMULTANEOUS;
      else if (strcmp(argv[ia],"stg") == 0)
        sensi_meth = STAGGERED;
      else if (strcmp(argv[ia],"stg1") == 0)
        sensi_meth = STAGGERED1;
      else 
        WrongArgs(my_pe, argv);

      ia++;

      if (strcmp(argv[ia],"full") == 0)
        err_con = FULL;
      else if (strcmp(argv[ia],"partial") == 0)
        err_con = PARTIAL;
      else
        WrongArgs(my_pe, argv);
      
      ia++;

    }

    else 

      WrongArgs(my_pe, argv);

  }

  /* Allocate and load user data block; allocate preconditioner block */
  data = (UserData) malloc(sizeof *data);
  data->p = (realtype *) malloc(NP*sizeof(realtype));
  InitUserData(my_pe, comm, data);
  predata = AllocPreconData (data);

  nvSpec = NV_SpecInit_Parallel(comm, local_N, neq, &argc, &argv);

  /* Allocate u, and set initial values and tolerances */ 
  u = N_VNew(nvSpec);
  SetInitialProfiles(u, data);
  abstol = ATOL; reltol = RTOL;

  /* CVODE_CREATE & CVODE_MALLOC */

  cvode_mem = CVodeCreate(BDF, NEWTON);
  if (cvode_mem == NULL) { 
    if (my_pe == 0) printf("CVodeCreate failed.\n"); 
    return(1); 
  }

  flag = CVodeSetFdata(cvode_mem, data);
  if (flag != SUCCESS) { 
    if (my_pe == 0) printf("CVodeSetFdata failed.\n"); 
    return(1); 
  }

  flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
  if (flag != SUCCESS) { 
    if (my_pe == 0) printf("CVodeSetMaxNumSteps failed.\n"); 
    return(1); 
  }


  if (CSderivs) {
    data_im = (UserData) malloc(sizeof *data);
    data_im->p = (realtype *) malloc(NP*sizeof(realtype));
    InitUserDataIm(data_im);
    flag = CVodeSetCSDerivs(cvode_mem, f_C_Cplx, data_im);
  }  

  flag = CVodeMalloc(cvode_mem, f, T0, u, SS, &reltol, &abstol, nvSpec);
  if (flag != SUCCESS) { 
    if (my_pe == 0) printf("CVodeMalloc failed.\n"); 
    return(1); 
  }

  /* CVSPGMR */
  flag = CVSpgmr(cvode_mem, LEFT, 0);
  if (flag != SUCCESS) { 
    if (my_pe == 0) printf("CVSpgmr failed."); 
    return(1); 
  }

  if (CSderivs)
    flag = CVSpgmrSetCSJacTimesVec(cvode_mem);

  flag = CVSpgmrSetPrecSetupFn(cvode_mem, Precond);
  flag = CVSpgmrSetPrecSolveFn(cvode_mem, PSolve);
  flag = CVSpgmrSetPrecData(cvode_mem, predata);

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
    if (flag != SUCCESS) {
      if (my_pe == 0) printf("CVodeSensMalloc failed, flag=%d\n",flag);
      return(1);
    }
  }

  if (my_pe == 0) {
    printf("\n2-species diurnal advection-diffusion problem\n\n");
    printf("========================================================================\n");
    printf("     T     Q       H      NST                    Bottom left  Top right \n");
    printf("========================================================================\n");
  }

  /* In loop over output points, call CVode, print results, test for error */
  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    if (flag != SUCCESS) {
      if (my_pe == 0) printf("CVode failed, flag=%d.\n", flag);
      break;
    }
    PrintOutput(cvode_mem, my_pe, comm, t, u);
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, t, uS);
      if (flag != SUCCESS) { 
        printf("CVodeSensExtract failed, flag=%d.\n", flag); 
        break; 
      }
      PrintOutputS(my_pe, comm, uS);
    }
    if (my_pe == 0)
      printf("------------------------------------------------------------------------\n");
  }

  /* Print final statistics */  
  if (my_pe == 0) {
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
  }

  /* Free memory */
  N_VFree(u);
  if(sensi) 
    N_VFree_S(NS, uS);
  free(data->p);  
  free(data);
  if (CSderivs) {
    free(data_im->p);  
    free(data_im);
  }
  FreePreconData(predata);
  CVodeFree(cvode_mem);
  NV_SpecFree_Parallel(nvSpec);
  MPI_Finalize();

  return(0);
}


/*********************** Private Helper Functions ************************/

/* ======================================================================= */
/* Exit if arguments are incorrect */

static void WrongArgs(int my_pe, char *argv[])
{
  if(my_pe == 0) {
    printf("\nUsage: %s [OPTIONS]\n",argv[0]);
    printf("            -x          = use xs4c\n");
    printf("            -s meth err = perform SA\n");
    printf("               meth : sim, stg, or stg1\n");
    printf("               err  : full or partial\n");
  }
  exit(0);
}

/* ======================================================================= */
/* Allocate memory for data structure of type UserData */

static PreconData AllocPreconData(UserData fdata)
{
  int lx, ly;
  PreconData pdata;

  pdata = (PreconData) malloc(sizeof *pdata);
  pdata->f_data = fdata;

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      (pdata->P)[lx][ly] = denalloc(NVARS);
      (pdata->Jbd)[lx][ly] = denalloc(NVARS);
      (pdata->pivot)[lx][ly] = denallocpiv(NVARS);
    }
  }

  return(pdata);
}

/* ======================================================================= */
/* Load constants in data */

static void InitUserData(int my_pe, MPI_Comm comm, UserData data)
{
  integertype isubx, isuby;
  realtype KH, VEL, KV0;

  /* Set problem parameters */
  data->p[0]  = 1.63e-16;       /* Q1  coeffs. q1, q2, c3             */
  data->p[1]  = 4.66e-16;       /* Q2                                 */
  data->p[2]  = 3.7e16;         /* C3                                 */
  data->p[3]  = 22.62;          /* A3  coeff. in expression for q3(t) */
  data->p[4]  = 7.601;          /* A4  coeff. in expression for q4(t) */
  KH  = data->p[5]  = 4.0e-6;   /* KH  horizontal diffusivity Kh      */ 
  VEL = data->p[6]  = 0.001;    /* VEL advection velocity V           */
  KV0 = data->p[7]  = 1.0e-8;   /* KV0 coeff. in Kv(z)                */ 

  /* Set problem constants */
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/((realtype)(MX-1));
  data->dy = (YMAX-YMIN)/((realtype)(MY-1));
  data->hdco = KH/SQR(data->dx);
  data->haco = VEL/(2.0*data->dx);
  data->vdco = (1.0/SQR(data->dy))*KV0;

  /* Set machine-related constants */
  data->comm = comm;
  data->my_pe = my_pe;
  /* isubx and isuby are the PE grid indices corresponding to my_pe */
  isuby = my_pe/NPEX;
  isubx = my_pe - isuby*NPEX;
  data->isubx = isubx;
  data->isuby = isuby;
  /* Set the sizes of a boundary x-line in u and uext */
  data->nvmxsub = NVARS*MXSUB;
  data->nvmxsub2 = NVARS*(MXSUB+2);

}

static void InitUserDataIm(UserData data)
{
  /* Set problem parameters */
  data->p[0]  =  0.0; 
  data->p[1]  =  0.0; 
  data->p[2]  =  0.0; 
  data->p[3]  =  0.0; 
  data->p[4]  =  0.0; 
  data->p[5]  =  0.0; 
  data->p[6]  =  0.0; 
  data->p[7]  =  0.0; 

  /* Set problem constants */
  data->om =  0.0; 
  data->dx =  0.0; 
  data->dy =  0.0; 
  data->hdco =  0.0; 
  data->haco =  0.0; 
  data->vdco =  0.0; 
}

/* ======================================================================= */
/* Free data memory */

static void FreePreconData(PreconData pdata)
{
  int lx, ly;

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      denfree((pdata->P)[lx][ly]);
      denfree((pdata->Jbd)[lx][ly]);
      denfreepiv((pdata->pivot)[lx][ly]);
    }
  }

  free(pdata);
}

/* ======================================================================= */
/* Set initial conditions in u */

static void SetInitialProfiles(N_Vector u, UserData data)
{
  integertype isubx, isuby, lx, ly, jx, jy, offset;
  realtype dx, dy, x, y, cx, cy, xmid, ymid;
  realtype *udata;

  /* Set pointer to data array in vector u */
  udata = NV_DATA_P(u);

  /* Get mesh spacings, and subgrid indices for this PE */
  dx = data->dx;         dy = data->dy;
  isubx = data->isubx;   isuby = data->isuby;

  /* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. */

  offset = 0;
  xmid = .5*(XMIN + XMAX);
  ymid = .5*(YMIN + YMAX);
  for (ly = 0; ly < MYSUB; ly++) {
    jy = ly + isuby*MYSUB;
    y = YMIN + jy*dy;
    cy = SQR(0.1*(y - ymid));
    cy = 1.0 - cy + 0.5*SQR(cy);
    for (lx = 0; lx < MXSUB; lx++) {
      jx = lx + isubx*MXSUB;
      x = XMIN + jx*dx;
      cx = SQR(0.1*(x - xmid));
      cx = 1.0 - cx + 0.5*SQR(cx);
      udata[offset  ] = C1_SCALE*cx*cy; 
      udata[offset+1] = C2_SCALE*cx*cy;
      offset = offset + 2;
    }
  }
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, int my_pe, MPI_Comm comm,
                        realtype t, N_Vector u)
{
  int nst, qu;
  realtype hu, *udata, tempu[2];
  integertype npelast, i0, i1;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;
  udata = NV_DATA_P(u);

  /* Send c at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&udata[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      tempu[0] = udata[i0];
      tempu[1] = udata[i1];
    }
  }

  /* On PE 0, receive c at top right, then print performance data
     and sampled solution values */ 
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&tempu[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    CVodeGetNumSteps(cvode_mem, &nst);
    CVodeGetLastOrder(cvode_mem, &qu);
    CVodeGetLastStep(cvode_mem, &hu);
    printf("%8.3e %2d  %8.3e %5d\n", t,qu,hu,nst);
    printf("                                Solution       ");
    printf("%12.4e %12.4e \n", udata[0], tempu[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", udata[1], tempu[1]);
  }
}

/* ======================================================================= */
/* Print sampled sensitivity values */

static void PrintOutputS(int my_pe, MPI_Comm comm, N_Vector *uS)
{
  realtype *sdata, temps[2];
  integertype npelast, i0, i1;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;

  sdata = NV_DATA_P(uS[0]);
  /* Send s1 at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&sdata[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      temps[0] = sdata[i0];
      temps[1] = sdata[i1];
    }
  }
  /* On PE 0, receive s1 at top right, then print sampled sensitivity values */ 
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&temps[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    printf("                                ----------------------------------------\n");
    printf("                                Sensitivity 1  ");
    printf("%12.4e %12.4e \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", sdata[1], temps[1]);
  }

  sdata = NV_DATA_P(uS[1]);
  /* Send s2 at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&sdata[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      temps[0] = sdata[i0];
      temps[1] = sdata[i1];
    }
  }
  /* On PE 0, receive s2 at top right, then print sampled sensitivity values */ 
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&temps[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    printf("                                ----------------------------------------\n");
    printf("                                Sensitivity 2  ");
    printf("%12.4e %12.4e \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", sdata[1], temps[1]);
  }

}

/* ======================================================================= */
/* Print final statistics contained in iopt */

static void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
  int nst;
  int nfe, nsetups, nni, ncfn, netf;
  int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;

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
}

/* ======================================================================= */
/* Preconditioner setup routine. Generate and preprocess P. */

static int Precond(realtype tn, N_Vector u, N_Vector fu,
                   booleantype jok, booleantype *jcurPtr, 
                   realtype gamma, void *P_data, 
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype c1, c2, cydn, cyup, diag, ydn, yup, q4coef, dely, verdco, hordco;
  realtype **(*P)[MYSUB], **(*Jbd)[MYSUB];
  int ier;
  integertype nvmxsub, *(*pivot)[MYSUB], offset;
  int lx, ly, jx, jy, isubx, isuby;
  realtype *udata, **a, **j;
  PreconData predata;
  UserData data;
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Make local copies of pointers in P_data, pointer to u's data,
     and PE index pair */
  predata = (PreconData) P_data;
  data = (UserData) (predata->f_data);
  P = predata->P;
  Jbd = predata->Jbd;
  pivot = predata->pivot;
  udata = NV_DATA_P(u);
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;

  /* Load problem coefficients and parameters */
  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];
  A3 = data->p[3];
  A4 = data->p[4];
  KH = data->p[5];
  VEL = data->p[6];
  KV0 = data->p[7];

  if (jok) {  /* jok = TRUE: Copy Jbd to P */

    for (ly = 0; ly < MYSUB; ly++)
      for (lx = 0; lx < MXSUB; lx++)
        dencopy(Jbd[lx][ly], P[lx][ly], NVARS);
    *jcurPtr = FALSE;

  } else {    /* jok = FALSE: Generate Jbd from scratch and copy to P */

    /* Make local copies of problem variables, for efficiency */
    q4coef = data->q4;
    dely = data->dy;
    verdco = data->vdco;
    hordco  = data->hdco;
    
    /* Compute 2x2 diagonal Jacobian blocks (using q4 values 
       computed on the last f call).  Load into P. */
    for (ly = 0; ly < MYSUB; ly++) {
      jy = ly + isuby*MYSUB;
      ydn = YMIN + (jy - .5)*dely;
      yup = ydn + dely;
      cydn = verdco*exp(0.2*ydn);
      cyup = verdco*exp(0.2*yup);
      diag = -(cydn + cyup + 2.0*hordco);
      for (lx = 0; lx < MXSUB; lx++) {
        jx = lx + isubx*MXSUB;
        offset = lx*NVARS + ly*nvmxsub;
        c1 = udata[offset];
        c2 = udata[offset+1];
        j = Jbd[lx][ly];
        a = P[lx][ly];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + q4coef;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4coef) + diag;
        dencopy(j, a, NVARS);
      }
    }
    
    *jcurPtr = TRUE;

  }

  /* Scale by -gamma */
  for (ly = 0; ly < MYSUB; ly++)
    for (lx = 0; lx < MXSUB; lx++)
      denscale(-gamma, P[lx][ly], NVARS);
  
  /* Add identity matrix and do LU decompositions on blocks in place */
  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      denaddI(P[lx][ly], NVARS);
      ier = gefa(P[lx][ly], NVARS, pivot[lx][ly]);
      if (ier != 0) return(1);
    }
  }
  
  return(0);
}


/* ======================================================================= */
/* Preconditioner solve routine */

static int PSolve(realtype tn, N_Vector u, N_Vector fu, 
                  N_Vector r, N_Vector z, 
                  realtype gamma, realtype delta,
                  int lr, void *P_data, N_Vector vtemp)
{
  realtype **(*P)[MYSUB];
  integertype nvmxsub, *(*pivot)[MYSUB];
  int lx, ly;
  realtype *zdata, *v;
  PreconData predata;
  UserData data;

  /* Extract the P and pivot arrays from P_data */
  predata = (PreconData) P_data;
  data = (UserData) (predata->f_data);
  P = predata->P;
  pivot = predata->pivot;

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z.
     First copy vector r to z. */
  N_VScale(1.0, r, z);

  nvmxsub = data->nvmxsub;
  zdata = NV_DATA_P(z);

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      v = &(zdata[lx*NVARS + ly*nvmxsub]);
      gesl(P[lx][ly], NVARS, pivot[lx][ly], v);
    }
  }

  return(0);
}
