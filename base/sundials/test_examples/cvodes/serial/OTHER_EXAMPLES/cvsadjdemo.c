/*
 * -----------------------------------------------------------------
 * $Revision:$
 * $Date:$
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Integration of multiple adjoint systems example problem.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

typedef struct {
  realtype p[3];
} *UserData;

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int fB1(realtype t, N_Vector y, 
               N_Vector yB, N_Vector yBdot, void *f_dataB);
static int fB2(realtype t, N_Vector y, 
               N_Vector yB, N_Vector yBdot, void *f_dataB);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  void *cvode_mem;

  int NeqF;
  realtype tF;
  realtype reltolF, abstolF;
  N_Vector yF;

  realtype Tout;

  int steps, ncheck;

  int indexB1, indexB2, NeqB1, NeqB2;
  realtype tB1, tB2;
  realtype reltolB, abstolB;
  N_Vector yB1, yB2;

  realtype time, time1, time2;

  int flag, iout;
  long int nst;

  data = NULL;
  cvode_mem = NULL;
  yF = yB1 = yB2 = NULL;

  /* User data structure */

  data = (UserData) malloc(sizeof *data);
  data->p[0] = RCONST(1.0);
  data->p[1] = RCONST(2.0);
  data->p[2] = RCONST(3.0);

  /* Initialize yF */

  NeqF = 3;

  yF = N_VNew_Serial(NeqF);
  Ith(yF,1) = ONE;
  Ith(yF,2) = ONE;
  Ith(yF,3) = ONE;

  tF = 0.0;

  reltolF = 1.0e-6;
  abstolF = 1.0e-4;

  /* Create and allocate CVODES memory for forward run */

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  flag = CVodeMalloc(cvode_mem, f, tF, yF, CV_SS, reltolF, &abstolF);

  flag = CVodeSetFdata(cvode_mem, data);

  flag = CVDense(cvode_mem, NeqF);

  /* Allocate CVODEA memory */

  steps = 20;
  flag = CVodeAdjMalloc(cvode_mem, steps, CV_HERMITE);

  /* Perform forward run */

  printf("Forward integration ... ");

  Tout = 2.0;
  flag = CVodeF(cvode_mem, Tout, yF, &time, CV_NORMAL, &ncheck);
  flag = CVodeGetNumSteps(cvode_mem, &nst);

  printf("( nst = %ld, ncheck = %ld )\n", nst, ncheck);

  printf("  yF(%lf): %12.4le %12.4le %12.4le\n", 
         time, Ith(yF,1), Ith(yF,2), Ith(yF,3));

  /* Initialize yB1 and yB2  */

  NeqB1 = 3;
  tB1 = 2.0;
  yB1 = N_VNew_Serial(NeqB1);
  Ith(yB1,1) = ONE;
  Ith(yB1,2) = ZERO;
  Ith(yB1,3) = ZERO;

  NeqB2 = 2;
  tB2 = 1.0;
  yB2 = N_VNew_Serial(NeqB2);
  Ith(yB2,1) = ONE;
  Ith(yB2,2) = ZERO;


  reltolB = 1.0e-6;
  abstolB = 1.0e-4;

  /* Create and allocate CVODES memory for backward run */

  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB1);
  flag = CVodeMallocB(cvode_mem, indexB1, fB1, tB1, yB1, CV_SS, reltolB, &abstolB);
  flag = CVodeSetFdataB(cvode_mem, indexB1, data);
  flag = CVDenseB(cvode_mem, indexB1, NeqB1);

  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB2);
  flag = CVodeMallocB(cvode_mem, indexB2, fB2, tB2, yB2, CV_SS, reltolB, &abstolB);
  flag = CVodeSetFdataB(cvode_mem, indexB2, data);
  flag = CVDenseB(cvode_mem, indexB2, NeqB2);

  /* Backward Integration */

  printf("Backward integration ...\n");


  for (iout = 0; iout <= 20; iout++) {

    Tout = 2.0 - iout*0.1;

    printf("--------- Tout = %lf  --------- \n",Tout);

    flag = CVodeB(cvode_mem, Tout, CV_NORMAL);

    flag = CVodeGetB(cvode_mem, indexB1, &time1, yB1);
    
    flag = CVodeGetB(cvode_mem, indexB2, &time2, yB2);

    printf("  yB1(%lf): %12.4le %12.4le %12.4le\n", 
           time1, Ith(yB1,1), Ith(yB1,2), Ith(yB1,3));

    printf("  yB2(%lf): %12.4le %12.4le\n", 
           time2, Ith(yB2,1), Ith(yB2,2));

  }

  /*
  time = 2.0;
  while (1) {

    printf("-------------\n");

    flag = CVodeB(cvode_mem, 0.5, CV_ONE_STEP);

    flag = CVodeGetB(cvode_mem, indexB1, &time1, yB1);
    
    printf("  yB1(%lf): %12.4le %12.4le %12.4le\n", 
           time1, Ith(yB1,1), Ith(yB1,2), Ith(yB1,3));

    flag = CVodeGetB(cvode_mem, indexB2, &time2, yB2);

    printf("  yB2(%lf): %12.4le %12.4le\n", 
           time2, Ith(yB2,1), Ith(yB2,2));


    if ( (time1 < 0.5) && (time2<0.5) ) break;

  }
  */

  /* Free memory */

  printf("Free memory\n\n");


  CVodeFree(&cvode_mem);

  N_VDestroy_Serial(yF); 
  N_VDestroy_Serial(yB1);
  N_VDestroy_Serial(yB2);

  free(data);

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */


static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype y1, y2, y3, yd1, yd3;
  UserData data;
  realtype p1, p2, p3;

  data = (UserData) f_data;
  p1 = data->p[0]; 
  p2 = data->p[1]; 
  p3 = data->p[2];

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  Ith(ydot,1) = -p1*y1*y1 - p2*y3;
  Ith(ydot,2) = -y2;
  Ith(ydot,3) = -p3*y2*y3;

  return(0);
}


static int fB1(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *f_dataB)
{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  realtype l1, l2, l3;
  
  data = (UserData) f_dataB;
  p1 = data->p[0]; 
  p2 = data->p[1]; 
  p3 = data->p[2];

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  Ith(yBdot,1) = -p1*l1;
  Ith(yBdot,2) = -p2*y3*l1 - l2 - 2.0*p3*y2*l3;
  Ith(yBdot,3) = -p2*y2*l1;

  return(0);
}


static int fB2(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *f_dataB)
{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  realtype l1, l2;
  
  data = (UserData) f_dataB;
  p1 = data->p[0]; 
  p2 = data->p[1]; 
  p3 = data->p[2];

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 

  Ith(yBdot,1) = -p1*l1;
  Ith(yBdot,2) = -p2*y3*l1 - l2 - 2.0*p3*y2*l2;

  return(0);
}

