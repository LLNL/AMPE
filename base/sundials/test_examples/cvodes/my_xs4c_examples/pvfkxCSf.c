#include <math.h>
#include "sundialstypes.h"    /* definitions of realtype, integertype            */
//#include "cvodes.h"           /* main CVODES header file                         */
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
} UserDataStruct;

typedef UserDataStruct *UserData;

static void BSend(MPI_Comm comm, int my_pe, integertype isubx, 
                  integertype isuby, integertype dsizex, 
                  integertype dsizey, realtype udata[]);
static void BRecvPost(MPI_Comm comm, MPI_Request request[], int my_pe,
                      integertype isubx, integertype isuby,
                      integertype dsizex, integertype dsizey,
                      realtype uext[], realtype buffer[]);
static void BRecvWait(MPI_Request request[], integertype isubx, integertype isuby,
                      integertype dsizex, realtype uext[], realtype buffer[]);
static void ucomm(realtype t, N_Vector u, UserData data);
static void fcalc(realtype t, realtype udata[], realtype dudata[], UserData data);

/* Functions Called by the CVODES Solver */

void f(realtype t, N_Vector u, N_Vector udot, void *f_data);

/* ======================================================================= */
/* Routine to send boundary data to neighboring PEs */

static void BSend(MPI_Comm comm, int my_pe, integertype isubx, 
                  integertype isuby, integertype dsizex, integertype dsizey, 
                  realtype udata[])
{
  int i, ly;
  integertype offsetu, offsetbuf;
  realtype bufleft[NVARS*MYSUB], bufright[NVARS*MYSUB];

  /* If isuby > 0, send data from bottom x-line of u */
  if (isuby != 0)
    MPI_Send(&udata[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If isuby < NPEY-1, send data from top x-line of u */
  if (isuby != NPEY-1) {
    offsetu = (MYSUB-1)*dsizex;
    MPI_Send(&udata[offsetu], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If isubx > 0, send data from left y-line of u (via bufleft) */
  if (isubx != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = ly*dsizex;
      for (i = 0; i < NVARS; i++)
        bufleft[offsetbuf+i] = udata[offsetu+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);   
  }

  /* If isubx < NPEX-1, send data from right y-line of u (via bufright) */
  if (isubx != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = offsetbuf*MXSUB + (MXSUB-1)*NVARS;
      for (i = 0; i < NVARS; i++)
        bufright[offsetbuf+i] = udata[offsetu+i];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe+1, 0, comm);   
  }

}
 
/* ======================================================================= */
/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvPost(MPI_Comm comm, MPI_Request request[], int my_pe,
                      integertype isubx, integertype isuby,
                      integertype dsizex, integertype dsizey,
                      realtype uext[], realtype buffer[])
{
  integertype offsetue;
  /* Have bufleft and bufright use the same buffer */
  realtype *bufleft, *bufright;

  bufleft = buffer;
  bufright = buffer+NVARS*MYSUB;

  /* If isuby > 0, receive data for bottom x-line of uext */
  if (isuby != 0)
    MPI_Irecv(&uext[NVARS], dsizex, PVEC_REAL_MPI_TYPE,
              my_pe-NPEX, 0, comm, &request[0]);

  /* If isuby < NPEY-1, receive data for top x-line of uext */
  if (isuby != NPEY-1) {
    offsetue = NVARS*(1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&uext[offsetue], dsizex, PVEC_REAL_MPI_TYPE,
              my_pe+NPEX, 0, comm, &request[1]);
  }
  
  /* If isubx > 0, receive data for left y-line of uext (via bufleft) */
  if (isubx != 0) {
    MPI_Irecv(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE,
              my_pe-1, 0, comm, &request[2]);
  }
  
  /* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Irecv(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE,
              my_pe+1, 0, comm, &request[3]);
  }
  
}

/* ======================================================================= */
/* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvWait(MPI_Request request[], integertype isubx, integertype isuby,
                      integertype dsizex, realtype uext[], realtype buffer[])
{
  int i, ly;
  integertype dsizex2, offsetue, offsetbuf;
  realtype *bufleft, *bufright;
  MPI_Status status;

  bufleft = buffer;
  bufright = buffer+NVARS*MYSUB;
  
  
  dsizex2 = dsizex + 2*NVARS;
  
  /* If isuby > 0, receive data for bottom x-line of uext */
  if (isuby != 0)
    MPI_Wait(&request[0],&status);
  
  /* If isuby < NPEY-1, receive data for top x-line of uext */
  if (isuby != NPEY-1)
    MPI_Wait(&request[1],&status);

  /* If isubx > 0, receive data for left y-line of uext (via bufleft) */
  if (isubx != 0) {
    MPI_Wait(&request[2],&status);

    /* Copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetue = (ly+1)*dsizex2;
      for (i = 0; i < NVARS; i++)
        uext[offsetue+i] = bufleft[offsetbuf+i];
    }
  }

  /* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Wait(&request[3],&status);

    /* Copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetue = (ly+2)*dsizex2 - NVARS;
      for (i = 0; i < NVARS; i++)
        uext[offsetue+i] = bufright[offsetbuf+i];
    }
  }

}

/* ======================================================================= */
/* ucomm routine.  This routine performs all communication 
   between processors of data needed to calculate f. */

static void ucomm(realtype t, N_Vector u, UserData data)
{
  realtype *udata, *uext, buffer[2*NVARS*MYSUB];
  MPI_Comm comm;
  int my_pe;
  integertype isubx, isuby, nvmxsub, nvmysub;
  MPI_Request request[4];

  udata = NV_DATA_P(u);

  /* Get comm, my_pe, subgrid indices, data sizes, extended array uext */
  comm = data->comm;  my_pe = data->my_pe;
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;
  nvmysub = NVARS*MYSUB;
  uext = data->uext;

  /* Start receiving boundary data from neighboring PEs */
  BRecvPost(comm, request, my_pe, isubx, isuby, nvmxsub, nvmysub, uext, buffer);

  /* Send data from boundary of local grid to neighboring PEs */
  BSend(comm, my_pe, isubx, isuby, nvmxsub, nvmysub, udata);

  /* Finish receiving boundary data from neighboring PEs */
  BRecvWait(request, isubx, isuby, nvmxsub, uext, buffer);
}

/* ======================================================================= */
/* fcalc routine. Compute f(t,y).  This routine assumes that communication 
   between processors of data needed to calculate f has already been done,
   and this data is in the work array uext. */

static void fcalc(realtype t, realtype udata[], realtype dudata[], UserData data)
{
  realtype *uext;
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  realtype q4coef, dely, verdco, hordco, horaco;
  int i, lx, ly, jx, jy;
  integertype isubx, isuby, nvmxsub, nvmxsub2, offsetu, offsetue;
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Get subgrid indices, data sizes, extended work array uext */
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub; nvmxsub2 = data->nvmxsub2;
  uext = data->uext;

  /* Load problem coefficients and parameters */
  Q1  = data->p[0];
  Q2  = data->p[1];
  C3  = data->p[2];
  A3  = data->p[3];
  A4  = data->p[4];
  KH  = data->p[5];
  VEL = data->p[6];
  KV0 = data->p[7];

  /* Copy local segment of u vector into the working extended array uext */
  offsetu = 0;
  offsetue = nvmxsub2 + NVARS;
  for (ly = 0; ly < MYSUB; ly++) {
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = udata[offsetu+i];
    offsetu = offsetu + nvmxsub;
    offsetue = offsetue + nvmxsub2;
  }

  /* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of u to uext */

  /* If isuby = 0, copy x-line 2 of u to uext */
  if (isuby == 0) {
    for (i = 0; i < nvmxsub; i++) uext[NVARS+i] = udata[nvmxsub+i];
  }

  /* If isuby = NPEY-1, copy x-line MYSUB-1 of u to uext */
  if (isuby == NPEY-1) {
    offsetu = (MYSUB-2)*nvmxsub;
    offsetue = (MYSUB+1)*nvmxsub2 + NVARS;
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = udata[offsetu+i];
  }

  /* If isubx = 0, copy y-line 2 of u to uext */
  if (isubx == 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*nvmxsub + NVARS;
      offsetue = (ly+1)*nvmxsub2;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = udata[offsetu+i];
    }
  }

  /* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext */
  if (isubx == NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = (ly+1)*nvmxsub - 2*NVARS;
      offsetue = (ly+2)*nvmxsub2 - NVARS;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = udata[offsetu+i];
    }
  }

  /* Make local copies of problem variables, for efficiency */
  dely   = data->dy;
  verdco = data->vdco;
  hordco = data->hdco;
  horaco = data->haco;

  /* Set diurnal rate coefficients as functions of t, and save q4 in 
  data block for use by preconditioner evaluation routine */
  s = sin((data->om)*t);
  if (s > 0.0) {
    q3 = exp(-A3/s);
    q4coef = exp(-A4/s);
  } else {
    q3 = 0.0;
    q4coef = 0.0;
  }
  data->q4 = q4coef;

  /* Loop over all grid points in local subgrid */
  for (ly = 0; ly < MYSUB; ly++) {
    jy = ly + isuby*MYSUB;
    /* Set vertical diffusion coefficients at jy +- 1/2 */
    ydn = YMIN + (jy - .5)*dely;
    yup = ydn + dely;
    cydn = verdco*exp(0.2*ydn);
    cyup = verdco*exp(0.2*yup);
    for (lx = 0; lx < MXSUB; lx++) {
      jx = lx + isubx*MXSUB;
      /* Extract c1 and c2, and set kinetic rate terms */
      offsetue = (lx+1)*NVARS + (ly+1)*nvmxsub2;
      c1 = uext[offsetue];
      c2 = uext[offsetue+1];
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;
      /* Set vertical diffusion terms */
      c1dn = uext[offsetue-nvmxsub2];
      c2dn = uext[offsetue-nvmxsub2+1];
      c1up = uext[offsetue+nvmxsub2];
      c2up = uext[offsetue+nvmxsub2+1];
      vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);
      /* Set horizontal diffusion and advection terms */
      c1lt = uext[offsetue-2];
      c2lt = uext[offsetue-1];
      c1rt = uext[offsetue+2];
      c2rt = uext[offsetue+3];
      hord1 = hordco*(c1rt - 2.0*c1 + c1lt);
      hord2 = hordco*(c2rt - 2.0*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);
      /* Load all terms into dudata */
      offsetu = lx*NVARS + ly*nvmxsub;
      dudata[offsetu]   = vertd1 + hord1 + horad1 + rkin1; 
      dudata[offsetu+1] = vertd2 + hord2 + horad2 + rkin2;
    }
  }

}

/***************** Functions Called by the CVODES Solver ******************/

/* ======================================================================= */
/* f routine.  Evaluate f(t,y).  First call ucomm to do communication of 
   subgrid boundary data into uext.  Then calculate f by a call to fcalc. */

void f(realtype t, N_Vector u, N_Vector udot, void *f_data)
{
  realtype *udata, *dudata;
  UserDataStruct *data;

  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);
  data = (UserData) f_data;

  /* Call ucomm to do inter-processor communicaiton */
  ucomm (t, u, data);

  /* Call fcalc to calculate all right-hand sides */
  fcalc (t, udata, dudata, data);
}

