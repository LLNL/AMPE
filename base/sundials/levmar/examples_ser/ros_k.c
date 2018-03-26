#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"
#include "levmar.h"
#include "lmspgmr.h"
#include "nvector_serial.h"
#include "sundialsmath.h"


static realtype cf_fun(N_Vector p, void *f_data);
static void gr_fun(N_Vector p, N_Vector g, void *f_data);
static int hv_fun(N_Vector v, N_Vector Hv, N_Vector pp, booleantype *new_pp, 
                  void *H_data);

/***************************** Main Program ******************************/

int main()
{
  N_Vector p, s;
  int flag, maxl, nni;
  void *mem;
  realtype *pdata, cfval, grd_inf;

  /* Create serial vectors of length  */
  p = N_VNew_Serial(2);
  s = N_VNew_Serial(2);

  pdata = NV_DATA_S(p);
  pdata[0] = -1.2;
  pdata[1] = 1.0;

  N_VConst(1.0, s);

  mem = LMCreate();
  LMSetInitMu(mem, 0.577);
  LMSetEps1(mem, 1.0e-6);
  LMSetEps2(mem, 1.0e-10);
  LMSetMaxIters(mem, 50);
  flag = LMMalloc(mem, cf_fun, gr_fun, p);
  printf("flag LMMalloc = %d\n",flag);

  maxl = 2;
  flag = LMSpgmr(mem, maxl);
  printf("flag LMSpgmr = %d\n",flag);
  LMSpgmrSetHesTimesVecFn(mem, hv_fun);

  flag = LMSolve(mem, p, s, s);
  printf("flag LMSolve = %d\n",flag);

  LMGetNumNonlinSolvIters(mem, &nni);
  LMGetCost(mem, &cfval);
  LMGetGradNorm(mem, &grd_inf);

  printf("p =  %e  %e\n",pdata[0],pdata[1]);
  printf("nni = %d   cf = %e  ||grd|| = %e\n",nni,cfval,grd_inf);

  N_VDestroy(p);
  N_VDestroy(s);
  LMFree(mem);

  return(0);

}

static realtype cf_fun(N_Vector p, void *f_data)
{
  realtype *pdata, f1, f2, cf;

  pdata = NV_DATA_S(p);

  f1 = 10.0 * (pdata[1] - pdata[0]*pdata[0]);
  f2 = 1.0 - pdata[0];

  cf = 0.5*(f1*f1+f2*f2);

  //  printf("p:  %e  %e\n",pdata[0],pdata[1]);
  //  printf("cf: %e  %e   %e\n",f1,f2,cf);

  return(cf);

}
static void gr_fun(N_Vector p, N_Vector g, void *f_data)
{
  realtype *pdata, *gdata, f1, f2;

  pdata = NV_DATA_S(p);
  gdata = NV_DATA_S(g);

  f1 = 10.0 * (pdata[1] - pdata[0]*pdata[0]);
  f2 = 1.0 - pdata[0];

  gdata[0] = -20.0*pdata[0]*f1 - f2;
  gdata[1] = 10.0*f1;

  //  printf("g: %e  %e\n",gdata[0],gdata[1]);

}

static int hv_fun(N_Vector v, N_Vector Hv, N_Vector p, booleantype *new_p, 
                  void *H_data)
{
  realtype *pdata, *Hvdata, *vdata;

  pdata = NV_DATA_S(p);
  Hvdata = NV_DATA_S(Hv);
  vdata = NV_DATA_S(v);

  Hvdata[0] = 400.0*pdata[0]*pdata[0]*vdata[0] - 200.0*pdata[0]*vdata[1];
  Hvdata[1] = -200.0*pdata[0]*vdata[0] + 100.0*vdata[1];

  return(0);
}
