#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "levmar_impl.h"

#define loop for(;;)

#define ZERO     RCONST(0.0)
#define POINT01  RCONST(0.01)
#define ONETHIRD RCONST(.3333333333333333)
#define HALF     RCONST(0.5)
#define ONE      RCONST(1.0)
#define TWO      RCONST(2.0)



/*====================================================*/


static booleantype LMAllocVectors(LMMem lm_mem, N_Vector tmpl);
static void LMFreeVectors(LMMem lm_mem);

static int LMLinSolvDrv(LMMem lm_mem, N_Vector h, booleantype same_p);


/*====================================================*/

void *LMCreate(void)
{
  LMMem lm_mem;
  realtype uround;

  lm_mem = (LMMem) malloc(sizeof(struct LMMemRec));

  lm_mem->lm_uround = uround = UNIT_ROUNDOFF;

  /* set default values */
  lm_mem->lm_f_data = NULL;

  lm_mem->lm_mu0   = LM_INIT_MU;
  lm_mem->lm_eps1  = LM_STOP_EPS1;
  lm_mem->lm_eps2  = LM_STOP_EPS2;
  lm_mem->lm_maxit = LM_MAX_ITERS;

  /* eps should be changed */
  //lm_mem->lm_eps = POINT01 * RPowerR(uround,ONETHIRD);
  lm_mem->lm_eps = POINT01 * RPowerR(uround,HALF);

  lm_mem->lm_sqrt_relerr = RSqrt(uround);

  lm_mem->lm_msbpre   = LM_MSBPRE;

  return((void *) lm_mem);

}

int LMSetFdata(void *lmmem, void *f_data)
{
  LMMem lm_mem;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }

  lm_mem = (LMMem) lmmem;
  lm_mem->lm_f_data = f_data;

  return(LM_SUCCESS);

}

#define f_data (lm_mem->lm_f_data)


/* scale factor for initial mu */
int LMSetInitMu(void *lmmem, realtype mu0)
{
  LMMem lm_mem;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }

  lm_mem = (LMMem) lmmem;
  lm_mem->lm_mu0 = mu0;

  return(LM_SUCCESS);

}

#define mu0 (lm_mem->lm_mu0)


/* stoping threshold for ||J^T e||_inf */
int LMSetEps1(void *lmmem, realtype eps1)
{
  LMMem lm_mem;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }

  lm_mem = (LMMem) lmmem;
  lm_mem->lm_eps1 = eps1;

  return(LM_SUCCESS);

}

#define eps1 (lm_mem->lm_eps1)

/* stoping threshold for ||Dp||_2 */
int LMSetEps2(void *lmmem, realtype eps2)
{
  LMMem lm_mem;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }

  lm_mem = (LMMem) lmmem;
  lm_mem->lm_eps2 = eps2;

  return(LM_SUCCESS);

}

#define eps2 (lm_mem->lm_eps2)

/* max number of iterations */
int LMSetMaxIters(void *lmmem, int maxit)
{
  LMMem lm_mem;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }

  lm_mem = (LMMem) lmmem;
  lm_mem->lm_maxit = maxit;

  return(LM_SUCCESS);

}

#define maxit (lm_mem->lm_maxit)


int LMSetMaxPrecCalls(void *lmmem, int msbpre)
{
  LMMem lm_mem;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }

  lm_mem = (LMMem) lmmem;

  if (msbpre < 0) {
    return(LMS_ILL_INPUT);
  }

  lm_mem->lm_msbpre = msbpre;

  return(LM_SUCCESS);
}

#define msbpre (lm_mem->lm_msbpre)


int LMMalloc(void *lmmem, CostFn func, GradFn grad, N_Vector tmpl)
{
  long int liw1, lrw1;
  LMMem lm_mem;
  booleantype allocOK;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }
  lm_mem = (LMMem) lmmem;

  if (func == NULL) {
    return(LM_ILL_INPUT);
  }

  if (grad == NULL) {
    return(LM_ILL_INPUT);
  }

  if (tmpl->ops->nvspace != NULL) {
    N_VSpace(tmpl, &lrw1, &liw1);
    lm_mem->lm_lrw1 = lrw1;
    lm_mem->lm_liw1 = liw1;
  }
  else {
    lm_mem->lm_lrw1 = 0;
    lm_mem->lm_liw1 = 0;
  }

  /* allocate necessary vectors */

  allocOK = LMAllocVectors(lm_mem, tmpl);
  if (!allocOK) {
    free(lm_mem);
    return(LM_MEM_FAIL);
  }

  /* copy the input parameter into LMSol state */

  lm_mem->lm_func = func;
  lm_mem->lm_grad = grad;

  /* set the linear solver addresses to NULL */

  lm_mem->lm_linit  = NULL;
  lm_mem->lm_lsetup = NULL;
  lm_mem->lm_lsolve = NULL;
  lm_mem->lm_lfree  = NULL;

  lm_mem->lm_lmem   = NULL;
  
  /* problem memory has been successfully allocated */

  lm_mem->lm_MallocDone = TRUE;

  return(LM_SUCCESS);

}

/* Readability replacements */

#define func       (lm_mem->lm_func)
#define grad       (lm_mem->lm_grad)

#define nni        (lm_mem->lm_nni)
#define nnilsetup  (lm_mem->lm_nnilsetup)
#define nfe        (lm_mem->lm_nfe)
#define nge        (lm_mem->lm_nge)

#define uround     (lm_mem->lm_uround)

#define liw1       (lm_mem->lm_liw1)
#define lrw1       (lm_mem->lm_lrw1)
#define liw        (lm_mem->lm_liw)
#define lrw        (lm_mem->lm_lrw)

#define pscale     (lm_mem->lm_pscale)
#define gscale     (lm_mem->lm_gscale)

#define mu         (lm_mem->lm_mu)
#define nu         (lm_mem->lm_nu)

#define F          (lm_mem->lm_F)
#define grd_inf    (lm_mem->lm_grd_inf)

#define pp         (lm_mem->lm_pp)
#define gval        (lm_mem->lm_gval)
#define h          (lm_mem->lm_h)
#define p_new      (lm_mem->lm_p_new)

#define vtemp1     (lm_mem->lm_vtemp1)
#define vtemp2     (lm_mem->lm_vtemp2)

#define setupNonNull (lm_mem->lm_setupNonNull)
#define setupCurrent (lm_mem->lm_setupCurrent)

#define linit   (lm_mem->lm_linit)
#define lsetup  (lm_mem->lm_lsetup)
#define lsolve  (lm_mem->lm_lsolve)
#define lfree   (lm_mem->lm_lfree)

int LMSolve(void *lmmem, N_Vector p, N_Vector p_scale, N_Vector g_scale)
{

  LMMem lm_mem;

  int ret;

  realtype Fn, dF, dL;
  realtype h_L2, p_L2;
  realtype tmp;

  booleantype same_p;

  if (lmmem == NULL) {
    return(LM_MEM_NULL);
  }
  lm_mem = (LMMem) lmmem;

  pp = p;
  same_p = FALSE;

  pscale = p_scale;
  gscale = g_scale;

  /* Initializations */
  nni = 1;
  nnilsetup = 1;
  mu  = mu0;
  nu  = 2;
  nfe = 0;
  nge = 0;

  h_L2 = ZERO;

  /* initialize the linear solver */
  ret = linit(lm_mem);
  if(ret!=0) return(LM_LINIT_FAIL);

  /* compute the cost function F = 0.5 * [ f(p)^T * f(p) ] */
  F = func(p, f_data);   nfe++;

  /* compute gradient G = J^T e and the inf norm ||J^T e||_inf */
  grad(p, gval, f_data);    nge++;
  N_VAbs(gval, vtemp1);
  grd_inf = N_VMin(vtemp1);

  loop {

    if (grd_inf <= eps1) return(LM_SMALL_GRD);

    /* solve augmented equations ( J^T J + mu I ) h = -G */
    ret = LMLinSolvDrv(lm_mem, h, same_p);
    if (ret != 0) break;

    /* compute norm of update ||h||_2 */
    h_L2 = N_VDotProd(h,h);
    h_L2 = RSqrt(h_L2);

    /* ||p||_2 */
    p_L2 = N_VDotProd(p,p);
    p_L2 = RSqrt(p_L2);

    tmp = eps2 + p_L2;

    if(h_L2 <= eps2*tmp)        /* relative change in p is small, stop */
      return(LM_SMALL_H);
    
    if(h_L2 >= tmp/uround)    /* almost singular, stop */
      return(LM_ALMOST_SING);
    
    /* Find new p */
    N_VLinearSum(ONE,p,ONE,h,p_new);

    N_VLinearSum(mu,h,-ONE,gval,vtemp1);
    dL = N_VDotProd(h,vtemp1);
    dL = HALF * dL;
        
    /* compute the cost function at the new p */
    Fn = func(p_new, f_data);   nfe++; 
        
    dF = F-Fn;

    if( dL > ZERO && dF > ZERO ){   /* reduction in error, increment is accepted */
      
      same_p = FALSE;

      /* update p's estimate */
      N_VScale(ONE,p_new,p);
      F = Fn;

      /* compute gradient G = J^T e and the inf norm ||J^T e||_inf */
      grad(p, gval, f_data);    nge++;
      N_VAbs(gval, vtemp1);
      grd_inf = N_VMin(vtemp1);

      tmp = (TWO*dF/dL-ONE);
      tmp = ONE - tmp*tmp*tmp;
      mu = mu * ( (tmp>=ONETHIRD)? tmp : ONETHIRD );
      nu=2;

    } else {                        /* same p, increase mu */
      
      same_p = TRUE;

      mu*=nu;
      nu*=2;
  
    }

    nni++;

    if ( nni > maxit) return(LM_REACHED_MAXIT);

  }
  
  return(0);

}


int LMGetIntWorkSpace(void *lmmem, long int *leniw)
{
  LMMem lm_mem;

  if (lmmem == NULL) return(LM_MEM_NULL);

  lm_mem = (LMMem) lmmem;

  *leniw = liw;

  return(LM_OKAY);
}

int LMGetRealWorkSpace(void *lmmem, long int *lenrw)
{
  LMMem lm_mem;

  if (lmmem == NULL) return(LM_MEM_NULL);

  lm_mem = (LMMem) lmmem;

  *lenrw = lrw;

  return(LM_OKAY);
}

int LMGetNumNonlinSolvIters(void *lmmem, int *nniters)
{
  LMMem lm_mem;

  if (lmmem == NULL) return(LM_MEM_NULL);

  lm_mem = (LMMem) lmmem;

  *nniters = nni;

  return(LM_OKAY);
}

int LMGetNumFuncEvals(void *lmmem, int *nfevals)
{
  LMMem lm_mem;

  if (lmmem == NULL) return(LM_MEM_NULL);

  lm_mem = (LMMem) lmmem;

  *nfevals = nfe;

  return(LM_OKAY);
}

int LMGetNumGradEvals(void *lmmem, int *ngevals)
{
  LMMem lm_mem;

  if (lmmem == NULL) return(LM_MEM_NULL);

  lm_mem = (LMMem) lmmem;

  *ngevals = nge;

  return(LM_OKAY);
}

int LMGetCost(void *lmmem, realtype *cfval)
{
  LMMem lm_mem;

  if (lmmem == NULL) return(LM_MEM_NULL);

  lm_mem = (LMMem) lmmem;

  *cfval = F;

  return(LM_OKAY);
}

int LMGetGradNorm(void *lmmem, realtype *grad_inf_norm)
{
  LMMem lm_mem;

  if (lmmem == NULL) return(LM_MEM_NULL);

  lm_mem = (LMMem) lmmem;

  *grad_inf_norm = grd_inf;

  return(LM_OKAY);
}



void LMFree(void *lmmem)
{
  LMMem lm_mem;

  if (lmmem == NULL) return;

  lm_mem = (LMMem) lmmem;
  LMFreeVectors(lm_mem);
  lfree(lm_mem);
  free(lm_mem);
}


static booleantype LMAllocVectors(LMMem lm_mem, N_Vector tmpl)
{

  gval = N_VClone(tmpl);
  if (gval == NULL) {
    return(FALSE);
  }

  h = N_VClone(tmpl);
  if (h == NULL) {
    N_VDestroy(gval);
    return(FALSE);
  }

  p_new = N_VClone(tmpl);
  if (p_new == NULL) {
    N_VDestroy(gval);
    N_VDestroy(h);
    return(FALSE);
  }

  vtemp1 = N_VClone(tmpl);
  if (vtemp1 == NULL) {
    N_VDestroy(gval);
    N_VDestroy(h);
    N_VDestroy(p_new);
    return(FALSE);
  }

  vtemp2 = N_VClone(tmpl);
  if (vtemp2 == NULL) {
    N_VDestroy(gval);
    N_VDestroy(h);
    N_VDestroy(p_new);
    N_VDestroy(vtemp1);
    return(FALSE);
  }

  liw = 5*liw1;
  lrw = 5*lrw1;

  return(TRUE);

}

static void LMFreeVectors(LMMem lm_mem)
{
  N_VDestroy(gval);
  N_VDestroy(h);
  N_VDestroy(p_new);
  N_VDestroy(vtemp1);
  N_VDestroy(vtemp2);
}


static int LMLinSolvDrv(LMMem lm_mem, N_Vector x, booleantype same_p)
{
  int ret, diff_nni;
  N_Vector b;
  booleantype forceSetup;

  diff_nni = nni-nnilsetup;

  loop{

    if( setupNonNull ) {
      ret = lsetup(lm_mem, same_p, diff_nni, forceSetup, &setupCurrent);
      if (setupCurrent == TRUE) {
        nnilsetup = nni;
        forceSetup = FALSE;
      }
      if (ret != 0) return(LM_LSETUP_FAIL);
    }
    
    /* load vtemp2 with the current value of -grad */
    
    b = vtemp2;
    N_VScale(-ONE, gval, b);
    
    /* call the generic 'lsolve' routine to solve the system ( H + mu I ) x = b */
    
    ret = lsolve(lm_mem, x, b);
    
    /* If successful or unrecoverable, return now */
    if (ret == 0) return(LM_SUCCESS);
    if (ret < 0)  return(LM_LSOLVE_FAIL);

    /* This point can be reached only when using a Krylov solver */
    
    /* If we cannot call again the preconditioner setup, we failed */
    if (!setupNonNull || setupCurrent) return(LM_KRYLOV_FAIL);
    
    forceSetup = TRUE;

  }

}
