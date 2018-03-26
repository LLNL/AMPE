/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2007/08/20 20:56:24 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Rootfinding test
 *   y1' = y2
 *   y2' = -y1 + sin(t)
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>


static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

int main(int argc, char *argv[])
{
  void *cvode_mem;
  N_Vector y;
  realtype t;
  int flag;
  int rootsinfo[1];
  
  FILE *fout;

  y = N_VNew_Serial(2);
  NV_Ith_S(y,0) = 0.0;
  NV_Ith_S(y,1) = 1.0;

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeInit(cvode_mem, f, 0.0, y);
  flag = CVodeSStolerances(cvode_mem, 1.e-3, 1.e-3);
  flag = CVodeRootInit(cvode_mem, 1, g);
  flag = CVDense(cvode_mem, 2);

  fout = fopen("cvtestroot.out", "w");

  t = 0.0;
  while(t<10.0) {
    flag = CVode(cvode_mem, 10.0, y, &t, CV_ONE_STEP);
    fprintf(fout, "%g  %g  %g\n", t, NV_Ith_S(y,0), NV_Ith_S(y,1));
    if (flag == CV_ROOT_RETURN) {
      flag = CVodeGetRootInfo(cvode_mem, rootsinfo);
      printf("t = %g   root : %d\n", t, rootsinfo[0]);
    } else {
      printf("t = %g\n", t);
    }
  }

  fclose(fout);

  N_VDestroy_Serial(y);
  CVodeFree(&cvode_mem);

  return(0);
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y1, y2;

  y1 = NV_Ith_S(y,0);
  y2 = NV_Ith_S(y,1);

  NV_Ith_S(ydot,0) = y2;
  NV_Ith_S(ydot,1) = -y1 + sin(t);

  return(0);
}

static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
  gout[0] = NV_Ith_S(y,1);

  return(0);
}
