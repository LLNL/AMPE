#include <stdio.h>
#include <stdlib.h>
#include "lmspgmr_impl.h"
#include "levmar_impl.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "sundialsmath.h"
#include "iterative.h"
#include "spgmr.h"

/*
 * -----------------------------------------------------------------
 * private constants
 * -----------------------------------------------------------------
 */

#define MSBPRE 10
#define MAXRST 0

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/*
 * -----------------------------------------------------------------
 * function prototypes
 * -----------------------------------------------------------------
 */

/* LMSPGMR linit, lsetup, lsolve, and lfree routines */

static int LMSpgmrInit(LMMem lm_mem);
static int LMSpgmrSetup(LMMem lm_mem, booleantype same_p, int diff_nni,
                        booleantype forceSetup, booleantype *setupCurrent);
static int LMSpgmrSolve(LMMem lm_mem, N_Vector xx, N_Vector bb);
static int LMSpgmrFree(LMMem lm_mem);

/* LMSpgmr Atimes and PSolve routines called by generic SPGMR solver */

static int LMSpgmrAtimes(void *levmar_mem, N_Vector v, N_Vector z);
static int LMSpgmrPSolve(void *levmar_mem, N_Vector r, N_Vector z, int lr);

/* difference quotient approximation for Hessian times vector */

static int LMSpgmrDQHtimes(N_Vector v, N_Vector Hv,
                            N_Vector p, booleantype *new_p, 
                            void *hes_data);


/*
 * -----------------------------------------------------------------
 * readability replacements
 * -----------------------------------------------------------------
 */

#define lrw1           (lm_mem->lm_lrw1)
#define liw1           (lm_mem->lm_liw1)
#define uround         (lm_mem->lm_uround)
#define nni            (lm_mem->lm_nni)

#define grad           (lm_mem->lm_grad)
#define f_data         (lm_mem->lm_f_data)

#define linit          (lm_mem->lm_linit)
#define lsetup         (lm_mem->lm_lsetup)
#define lsolve         (lm_mem->lm_lsolve)
#define lfree          (lm_mem->lm_lfree)
#define lmem           (lm_mem->lm_lmem)

#define mu             (lm_mem->lm_mu)

#define pp             (lm_mem->lm_pp)
#define gval            (lm_mem->lm_p_new)
#define pscale         (lm_mem->lm_pscale)
#define gscale         (lm_mem->lm_gscale)
#define sqrt_relerr    (lm_mem->lm_sqrt_relerr)
#define setupNonNull   (lm_mem->lm_setupNonNull)

#define eps            (lm_mem->lm_eps)

#define vtemp1         (lm_mem->lm_vtemp1)
#define vec_tmpl       (lm_mem->lm_vtemp1)
#define vtemp2         (lm_mem->lm_vtemp2)

#define pretype   (lmspgmr_mem->g_pretype)
#define gstype    (lmspgmr_mem->g_gstype)
#define nli       (lmspgmr_mem->g_nli)
#define npe       (lmspgmr_mem->g_npe)
#define nps       (lmspgmr_mem->g_nps)
#define ncfl      (lmspgmr_mem->g_ncfl)
#define nhtimes   (lmspgmr_mem->g_nhtimes)
#define ngeSG     (lmspgmr_mem->g_ngeSG)
#define new_pp    (lmspgmr_mem->g_new_pp)
#define spgmr_mem (lmspgmr_mem->g_spgmr_mem)

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmr
 * -----------------------------------------------------------------
 * This routine allocates and initializes the memory record and
 * sets function fields specific to the SPGMR linear solver module.
 * LMSpgmr sets the lm_linit, lm_lsetup, lm_lsolve, and
 * lm_lfree fields in *lmmem to be LMSpgmrInit, LMSpgmrSetup,
 * LMSpgmrSolve, and LMSpgmrFree, respectively. It allocates
 * memory for a structure of type LMSpgmrMemRec and sets the
 * lm_lmem field in *lmmem to the address of this structure. It
 * also calls SpgmrMalloc to allocate memory for the module
 * SPGMR. In summary, LMSpgmr sets the following fields in the
 * LMSpgmrMemRec structure: 
 *                         
 *  pretype   = NONE
 *  gstype    = MODIFIED_GS
 *  g_maxl    = LMSPGMR_MAXL  if maxl <= 0             
 *            = maxl           if maxl > 0   
 *  g_maxlrst = 0 (default)
 *  g_pset    = NULL
 *  g_psolve  = NULL    
 *  g_htimes  = NULL                                    
 *  g_P_data  = NULL                                        
 * -----------------------------------------------------------------
 */

int LMSpgmr(void *lmmem, int maxl)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;
  int maxl1;

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);  
  lm_mem = (LMMem) lmmem;

  /* check for required vector operations */

  /* Note: do NOT need to check for N_VClone, N_VDestroy, N_VLinearSum, N_VDotProd,
     or N_VScale because they are required by LEVMAR */

  if ((vec_tmpl->ops->nvconst == NULL) ||
      (vec_tmpl->ops->nvprod == NULL) ||
      (vec_tmpl->ops->nvwl2norm == NULL) ||
      (vec_tmpl->ops->nvl1norm == NULL))
    return(LMSPGMR_ILL_INPUT);

  /* set four main function fields in lm_mem */

  linit  = LMSpgmrInit; 
  lsetup = LMSpgmrSetup;
  lsolve = LMSpgmrSolve;
  lfree  = LMSpgmrFree;

  /* get memory for LMSpgmrMemRec */

  lmspgmr_mem = (LMSpgmrMem) malloc(sizeof(LMSpgmrMemRec));
  if (lmspgmr_mem == NULL) return(LMSPGMR_MEM_FAIL);  

  /* set SPGMR parameters that were passed in call sequence */

  maxl1 = (maxl <= 0) ? LMSPGMR_MAXL : maxl;
  lmspgmr_mem->g_maxl = maxl1;  

  /* set default values for the rest of the SPGMR parameters */

  lmspgmr_mem->g_pretype   = NONE;
  lmspgmr_mem->g_gstype    = MODIFIED_GS;
  lmspgmr_mem->g_maxlrst   = MAXRST;
  lmspgmr_mem->g_msbpsetup = MSBPRE;
  lmspgmr_mem->g_P_data    = NULL;
  lmspgmr_mem->g_pset      = NULL;
  lmspgmr_mem->g_psolve    = NULL;
  lmspgmr_mem->g_htimes    = NULL;

  /* call SpgmrMalloc to allocate workspace for SPGMR */

  /* vec_templ passed as template vector */

  spgmr_mem = SpgmrMalloc(maxl1, vec_tmpl);
  if (spgmr_mem == NULL) {
    lmem = NULL;  /* set lmem to NULL and free that memory as a flag to a
                     later inadvertent Levmar call that SpgmrMalloc failed */
    free(lmem);
    return(LMSPGMR_MEM_FAIL);
  }

  /* attach linear solver memory to LEVMAR memory */

  lmem = lmspgmr_mem;

  return(LMSPGMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmrSetMaxRestarts
 * -----------------------------------------------------------------
 */

int LMSpgmrSetMaxRestarts(void *lmmem, int maxrs)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  /* check for legal maxrs */

  if (maxrs < 0) return(LMSPGMR_ILL_INPUT);
  lmspgmr_mem->g_maxlrst = maxrs;

  return(LMSPGMR_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSetMaxPrecCalls
 * -----------------------------------------------------------------
 */

int LMSpgmrSetMaxPrecCalls(void *lmmem, int msbpsetup)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  lmspgmr_mem->g_msbpsetup = msbpsetup;

  return(LMSPGMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSetPrecSetupFn
 * -----------------------------------------------------------------
 */

int LMSpgmrSetPrecSetupFn(void *lmmem, LMSpgmrPrecSetupFn pset)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  lmspgmr_mem->g_pset = pset;

  return(LMSPGMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSetPrecSolveFn
 * -----------------------------------------------------------------
 */

int LMSpgmrSetPrecSolveFn(void *lmmem, LMSpgmrPrecSolveFn psolve)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  lmspgmr_mem->g_psolve = psolve;

  return(LMSPGMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSetPrecData
 * -----------------------------------------------------------------
 */

int LMSpgmrSetPrecData(void *lmmem, void *P_data)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  lmspgmr_mem->g_P_data = P_data;

  return(LMSPGMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int LMSpgmrSetHesTimesVecFn(void *lmmem, LMSpgmrHesTimesVecFn htimes)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  lmspgmr_mem->g_htimes = htimes;

  return(LMSPGMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSetJacData
 * -----------------------------------------------------------------
 */

int LMSpgmrSetHesData(void *lmmem, void *H_data)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  lmspgmr_mem->g_H_data = H_data;

  return(LMSPGMR_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetIntWorkSpace
 * -----------------------------------------------------------------
 */

int LMSpgmrGetIntWorkSpace(void *lmmem, long int *leniwSG)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;
  int maxl;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  maxl = lmspgmr_mem->g_maxl;
  *leniwSG = liw1 * (maxl + 3);

  return(LMSPGMR_OKAY);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetRealWorkSpace
 * -----------------------------------------------------------------
 */

int LMSpgmrGetRealWorkSpace(void *lmmem, long int *lenrwSG)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;
  int maxl;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  maxl = lmspgmr_mem->g_maxl;
  *lenrwSG = lrw1 * (maxl + 3) + (maxl * (maxl + 4)) + 1;

  return(LMSPGMR_OKAY);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int LMSpgmrGetNumPrecEvals(void *lmmem, int *npevals)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  *npevals = npe;

  return(LMSPGMR_OKAY);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int LMSpgmrGetNumPrecSolves(void *lmmem, int *npsolves)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  *npsolves = nps;

  return(LMSPGMR_OKAY);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetNumLinIters
 * -----------------------------------------------------------------
 */

int LMSpgmrGetNumLinIters(void *lmmem, int *nliters)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  *nliters = nli;

  return(LMSPGMR_OKAY);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetNumConvFails
 * -----------------------------------------------------------------
 */

int LMSpgmrGetNumConvFails(void *lmmem, int *nlcfails)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  *nlcfails = ncfl;

  return(LMSPGMR_OKAY);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int LMSpgmrGetNumHtimesEvals(void *lmmem, int *nhvevals)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;

  *nhvevals = nhtimes;

  return(LMSPGMR_OKAY);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int LMSpgmrGetNumGradEvals(void *lmmem, int *ngevalsSG)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* return immediately if lmmem is NULL */

  if (lmmem == NULL) return(LMSPGMR_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMSPGMR_LMEM_NULL);
  lmspgmr_mem = (LMSpgmrMem) lmem;
  *ngevalsSG = ngeSG;

  return(LMSPGMR_OKAY);
}


/*
 * -----------------------------------------------------------------
 * additional readability replacements
 * -----------------------------------------------------------------
 */

#define maxl      (lmspgmr_mem->g_maxl)
#define maxlrst   (lmspgmr_mem->g_maxlrst)
#define msbpsetup (lmspgmr_mem->g_msbpsetup)
#define pset      (lmspgmr_mem->g_pset)
#define psolve    (lmspgmr_mem->g_psolve)
#define P_data    (lmspgmr_mem->g_P_data)
#define htimes    (lmspgmr_mem->g_htimes)
#define H_data    (lmspgmr_mem->g_H_data)

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrInit
 * -----------------------------------------------------------------
 * This routine initializes variables associated with the GMRES
 * linear solver. Memory allocation was done previously in
 * LMSpgmr.
 * -----------------------------------------------------------------
 */

static int LMSpgmrInit(LMMem lm_mem)
{
  LMSpgmrMem lmspgmr_mem;

  lmspgmr_mem = (LMSpgmrMem) lmem;

  /* initialize counters */

  npe = nli = nps = ncfl = 0;
  nhtimes = ngeSG = 0;

  /* set preconditioner type */

  if (psolve != NULL) {
    pretype = RIGHT;
  } else {
    pretype = NONE;
  }
  
  /* set setupNonNull to TRUE iff there is preconditioning with setup */

  setupNonNull = (psolve != NULL) && (pset != NULL);

  /* if htimes is NULL at this time, set it to private DQ routine */

  if (htimes == NULL) {
    htimes = LMSpgmrDQHtimes;
    H_data = lm_mem;
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the SPGMR linear
 * solver, that is, it is an interface to the user-supplied
 * routine pset.
 * -----------------------------------------------------------------
 */

static int LMSpgmrSetup(LMMem lm_mem, booleantype same_p, int diff_nni,
                        booleantype forceSetup, booleantype *setupCurrent)
{
  LMSpgmrMem lmspgmr_mem;
  int ret;

  lmspgmr_mem = (LMSpgmrMem) lmem;

  *setupCurrent = FALSE;

  if (forceSetup || (diff_nni>=msbpsetup)) {

    ret = pset(pp, pscale, gval, gscale, same_p, mu, P_data, vtemp1, vtemp2); 
    if (ret != 0) return(1);

    npe++;
    *setupCurrent = TRUE;

  }
  
  /* return the same value ret that pset returned */

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic SPGMR solver
 * SpgmrSolve for the solution of the linear system Ax = b.
 *
 * Appropriate variables are passed to SpgmrSolve and the counters
 * nli, nps, and ncfl are incremented, and the return value is set
 * according to the success of SpgmrSolve. The success flag is
 * returned if SpgmrSolve converged, or if this is the first Newton
 * iteration and the residual norm was reduced below its initial
 * value. Of the other error conditions, only preconditioner solver
 * failure is specifically returned. Otherwise a generic flag is
 * returned to denote failure of this routine.
 * -----------------------------------------------------------------
 */

static int LMSpgmrSolve(LMMem lm_mem, N_Vector xx, N_Vector bb)
{
  LMSpgmrMem lmspgmr_mem;
  int ret, nli_inc, nps_inc;
  realtype *res_norm;

  lmspgmr_mem = (LMSpgmrMem) lmem;

  /* Set initial guess to xx = 0. bb is set, by the routine
     calling LMSpgmrSolve, to the RHS vector for the system
     to be solved. */ 
 
  N_VConst(ZERO, xx);

  new_pp = TRUE;  /* set flag required for user Jacobian routine */

  /* call SpgmrSolve */

  ret = SpgmrSolve(spgmr_mem, lm_mem, xx, bb, pretype, gstype, eps, 
                   maxlrst, lm_mem, gscale, gscale, LMSpgmrAtimes,
                   LMSpgmrPSolve, res_norm, &nli_inc, &nps_inc);

  /* increment counters nli, nps, and ncfl */

  nli = nli + nli_inc;
  nps = nps + nps_inc;

  if (ret != 0) ncfl++;

  /* set return value to appropriate value */

  if ((ret == SPGMR_SUCCESS) || (ret == SPGMR_RES_REDUCED)) return(0);

  if (ret == SPGMR_PSOLVE_FAIL_REC) return(1);

  return(-1);  
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the SPGMR linear solver.
 * -----------------------------------------------------------------
 */

static int LMSpgmrFree(LMMem lm_mem)
{
  LMSpgmrMem lmspgmr_mem;

  lmspgmr_mem = (LMSpgmrMem) lmem;

  SpgmrFree(spgmr_mem);
  free(lmem);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrAtimes
 * -----------------------------------------------------------------
 * This routine coordinates the generation of the matrix-vector
 * product z = M*v, where M = H + mu * I 
 * It either calls LMSpgmrDQHtimes, which uses a difference quotient 
 * approximation for H*v, or it calls the user-supplied function 
 * LMSpgmrHesTimesVecFn if it is non-null.
 * -----------------------------------------------------------------
 */

static int LMSpgmrAtimes(void *levmar_mem, N_Vector v, N_Vector z)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;
  int ret;

  lm_mem = (LMMem) levmar_mem;
  lmspgmr_mem = (LMSpgmrMem) lmem;

  ret = htimes(v, z, pp, &new_pp, H_data);
  nhtimes++;
  if (ret!=0) return(ret);

  N_VLinearSum(ONE, z, mu, v, z);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic SpgmrSolve routine
 * and the user's psolve routine. It passes to psolve all required
 * state information from levmar_mem. Its return value is the same
 * as that returned by psolve. Note that the generic SPGMR solver
 * guarantees that LMSpgmrPSolve will not be called in the case in
 * which preconditioning is not done. This is the only case in which
 * the user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

static int LMSpgmrPSolve(void *levmar_mem, N_Vector r, N_Vector z, int lrdummy)
{
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;
  int ret;

  lm_mem = (LMMem) levmar_mem;
  lmspgmr_mem = (LMSpgmrMem) lmem;

  /* copy the rhs into z before the psolve call */   
  /* Note: z returns with the solution */

  N_VScale(ONE, r, z);

  /* this call is counted in nps within the LMSpgmrSolve routine */

  ret = psolve(pp, pscale, gval, gscale, z, mu, P_data, vtemp1);

  return(ret);     
}

/*
 * -----------------------------------------------------------------
 * Function : LMSpgmrDQHtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = H*v using a
 * difference quotient approximation. The approximation is 
 * H*v = [grd(pp + sigma*v) - grd(pp)]/sigma. Here sigma is based
 * on the dot products (pscale*pp, pscale*v) and
 * (pscale*v, pscale*v), the L1Norm(pscale*v), and on sqrt_relerr
 * (the square root of the relative error in the function). Note
 * that v in the argument list has already been both preconditioned
 * and unscaled.
 * -----------------------------------------------------------------
 */

static int LMSpgmrDQHtimes(N_Vector v, N_Vector Hv,
                           N_Vector p, booleantype *new_p, 
                           void *hes_data)
{
  realtype sigma, sigma_inv, sptsv, sq1norm, sign, vtv;
  LMMem lm_mem;
  LMSpgmrMem lmspgmr_mem;

  /* hes_data is lm_mem */

  lm_mem = (LMMem) hes_data;
  lmspgmr_mem = (LMSpgmrMem) lmem;

  /* scale the vector v and put Dp*v into vtemp1 */

  N_VProd(v, pscale, vtemp1);

  /* scale p and put into Hv (used as a temporary storage) */

  N_VProd(p, pscale, Hv);

  /* compute dot product (Dp*p).(Dp*v) */

  sptsv = N_VDotProd(Hv, vtemp1);

  /* compute dot product (Dp*v).(Dp*v) */

  vtv = N_VDotProd(vtemp1, vtemp1);

  sq1norm = N_VL1Norm(vtemp1);

  sign = (sptsv >= ZERO) ? ONE : -ONE ;
 
  /*  this expression for sigma is from p. 469, Brown and Saad paper */

  sigma = sign*sqrt_relerr*MAX(ABS(sptsv),sq1norm)/vtv; 

  sigma_inv = ONE/sigma;

  /* compute the p-prime at which to evaluate the gradient */

  N_VLinearSum(ONE, p, sigma, v, vtemp1);
 
  /* call the gradient function to calculate grad(p+sigma*v) */

  grad(vtemp1, vtemp2, f_data);    
  ngeSG++;

  /* finish the computation of the difference quotient */

  N_VLinearSum(sigma_inv, vtemp2, -sigma_inv, gval, Hv);

  return(0);
}
