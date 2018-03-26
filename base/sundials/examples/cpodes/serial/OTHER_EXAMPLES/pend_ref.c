#include <stdio.h>

#include <cpodes/cpodes.h>
#include <cpodes/cpodes_dense.h>
#include <nvector/nvector_serial.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Problem Constants */

#define RTOL  RCONST(1.0e-9)
#define ATOL  RCONST(1.0e-9)

/* Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static void PrintFinalStats(void *cpode_mem);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  void *cpode_mem;
  N_Vector yy, yp;
  realtype reltol, abstol, t, tout;
  realtype th, thd, x, y, xd, yd;
  int iout, flag;
  

  yy = N_VNew_Serial(2);
  yp = N_VNew_Serial(2);

  /* Initialize y */
  Ith(yy,1) = 0.0;  /* theta */
  Ith(yy,2) = 0.0;  /* thetad */

  /* Set tolerances */
  reltol = RTOL;
  abstol = ATOL;

  /* Initialize solver */
  cpode_mem = CPodeCreate(CP_EXPL, CP_BDF, CP_NEWTON);  
  flag = CPodeSetMaxNumSteps(cpode_mem, 100000);
  flag = CPodeInit(cpode_mem, f, NULL, 0.0, yy, yp, CP_SS, reltol, &abstol);
  flag = CPDense(cpode_mem, 2);


  flag = CPodeSetStopTime(cpode_mem, 20.0);
  flag = CPode(cpode_mem, 20.0, &t, yy, yp, CP_NORMAL_TSTOP);
  th  = Ith(yy,1);
  thd = Ith(yy,2);
  x = cos(th);
  y = sin(th);
  xd = -thd*sin(th);
  yd =  thd*cos(th);
  printf("%14.10e  %14.10e  %14.10e  %14.10e\n",
         x-1.0,y,xd,yd);




  /*
  t = 0.0;
  for(iout=1; iout<=50; iout++) {
    tout = iout*1.0;
    flag = CPode(cpode_mem, tout, &t, yy, yp, CP_NORMAL);
    if (flag < 0) break;
    th  = Ith(yy,1);
    thd = Ith(yy,2);

    x = cos(th);
    y = sin(th);

    xd = -thd*sin(th);
    yd =  thd*cos(th);

    printf("%lf  %14.10lf  %14.10lf  %14.10lf  %14.10lf    %14.10lf\n",  t, x,y,xd,yd,0.0);
  }
  */


  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  CPodeFree(&cpode_mem);

  return(0);
}


static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype th, thd, g;

  g = 13.7503716373294544;
  //  g = 0.1;

  th  = Ith(yy,1);
  thd  = Ith(yy,2);

  Ith(fy,1) = thd;
  Ith(fy,2) = -g*cos(th);

  return(0);
}

static void PrintFinalStats(void *cpode_mem)
{
  realtype h0u;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CPodeGetActualInitStep(cpode_mem, &h0u);
  flag = CPodeGetNumSteps(cpode_mem, &nst);
  flag = CPodeGetNumFctEvals(cpode_mem, &nfe);
  flag = CPodeGetNumLinSolvSetups(cpode_mem, &nsetups);
  flag = CPodeGetNumErrTestFails(cpode_mem, &netf);
  flag = CPodeGetNumNonlinSolvIters(cpode_mem, &nni);
  flag = CPodeGetNumNonlinSolvConvFails(cpode_mem, &ncfn);

  flag = CPDlsGetNumJacEvals(cpode_mem, &nje);
  flag = CPDlsGetNumFctEvals(cpode_mem, &nfeLS);

  flag = CPodeGetNumGEvals(cpode_mem, &nge);

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

