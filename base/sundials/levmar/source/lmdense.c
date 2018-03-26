#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lmdense_impl.h"
#include "levmar_impl.h"

#include "sundialsmath.h"


#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* LMDENSE linit, lsetup, lsolve, and lfree routines */
 
static int LMDenseInit(LMMem lm_mem);

static int LMDenseSetup(LMMem lm_mem, booleantype same_p, int diff_nni,
                        booleantype forceSetup, booleantype *setupCurrent);

static int LMDenseSolve(LMMem lm_mem, N_Vector xx, N_Vector bb);

static int LMDenseFree(LMMem lm_mem);

static void LMDenseDQHes(long int Np, DenseMat H, N_Vector p, 
                         void *hes_data, N_Vector tmp1, N_Vector tmp2);

/* Readability Replacements */

#define lrw1           (lm_mem->lm_lrw1)
#define liw1           (lm_mem->lm_liw1)
#define uround         (lm_mem->lm_uround)

#define grad           (lm_mem->lm_grad)
#define f_data         (lm_mem->lm_f_data)

#define linit          (lm_mem->lm_linit)
#define lsetup         (lm_mem->lm_lsetup)
#define lsolve         (lm_mem->lm_lsolve)
#define lfree          (lm_mem->lm_lfree)
#define lmem           (lm_mem->lm_lmem)

#define mu             (lm_mem->lm_mu)

#define pp             (lm_mem->lm_pp)
#define gval           (lm_mem->lm_gval)
#define pscale         (lm_mem->lm_pscale)
#define gscale         (lm_mem->lm_gscale)
#define setupNonNull   (lm_mem->lm_setupNonNull)

#define vtemp1         (lm_mem->lm_vtemp1)
#define vec_tmpl       (lm_mem->lm_vtemp1)
#define vtemp2         (lm_mem->lm_vtemp2)

#define n         (lmdense_mem->d_n)
#define hes       (lmdense_mem->d_hes)
#define M         (lmdense_mem->d_M)
#define pivots    (lmdense_mem->d_pivots)
#define savedH    (lmdense_mem->d_savedH)
#define nhe       (lmdense_mem->d_nhe)
#define ngeD      (lmdense_mem->d_ngeD)
#define H_data    (lmdense_mem->d_H_data)

                  
int LMDense(void *lmmem, long int Np)
{
  LMMem lm_mem;
  LMDenseMem lmdense_mem;

  /* Return immediately if lmmem is NULL */
  if (lmmem == NULL) return(LMDENSE_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL)
    return(LMDENSE_ILL_INPUT);

  if (lfree !=NULL) lfree(lm_mem);

  /* Set four main function fields in lm_mem */
  linit  = LMDenseInit;
  lsetup = LMDenseSetup;
  lsolve = LMDenseSolve;
  lfree  = LMDenseFree;

  /* Get memory for LMDenseMemRec */
  lmdense_mem = (LMDenseMem) malloc(sizeof(LMDenseMemRec));
  if (lmdense_mem == NULL) return(LMDENSE_MEM_FAIL);

  /* Set default H function and H data */
  hes = LMDenseDQHes;
  H_data = lmmem;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = Np;

  /* Allocate memory for M, savedH, and pivot array */
  
  M = DenseAllocMat(Np);
  if (M == NULL) return(LMDENSE_MEM_FAIL);
  savedH = DenseAllocMat(Np);
  if (savedH == NULL) {
    DenseFreeMat(M);
    return(LMDENSE_MEM_FAIL);
  }
  pivots = DenseAllocPiv(Np);
  if (pivots == NULL) {
    DenseFreeMat(M);
    DenseFreeMat(savedH);
    return(LMDENSE_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = lmdense_mem;

  return(LMDENSE_SUCCESS);
}

int LMDenseSetHesFn(void *lmmem, LMDenseHesFn dhes)
{
  LMMem lm_mem;
  LMDenseMem lmdense_mem;

  /* Return immediately if lmmem is NULL */
  if (lmmem == NULL) return(LMDENSE_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMDENSE_LMEM_NULL);
  lmdense_mem = (LMDenseMem) lmem;

  hes = dhes;

  return(LMDENSE_SUCCESS);
}

int LMDenseSetHesData(void *lmmem, void *hes_data)
{
  LMMem lm_mem;
  LMDenseMem lmdense_mem;

  /* Return immediately if lmmem is NULL */
  if (lmmem == NULL) return(LMDENSE_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMDENSE_LMEM_NULL);
  lmdense_mem = (LMDenseMem) lmem;

  H_data = hes_data;

  return(LMDENSE_SUCCESS);
}

int LMDenseGetIntWorkSpace(void *lmmem, long int *leniwD)
{
  LMMem lm_mem;
  LMDenseMem lmdense_mem;

  /* Return immediately if lmmem is NULL */
  if (lmmem == NULL) return(LMDENSE_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMDENSE_LMEM_NULL);
  lmdense_mem = (LMDenseMem) lmem;

  *leniwD = n;

  return(LMDENSE_OKAY);
}

int LMDenseGetRealWorkSpace(void *lmmem, long int *lenrwD)
{
  LMMem lm_mem;
  LMDenseMem lmdense_mem;

  /* Return immediately if lmmem is NULL */
  if (lmmem == NULL) return(LMDENSE_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMDENSE_LMEM_NULL);
  lmdense_mem = (LMDenseMem) lmem;

  *lenrwD = 2*n*n;

  return(LMDENSE_OKAY);
}

int LMDenseGetNumHesEvals(void *lmmem, int *nhevalsD)
{
  LMMem lm_mem;
  LMDenseMem lmdense_mem;

  /* Return immediately if lmmem is NULL */
  if (lmmem == NULL) return(LMDENSE_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMDENSE_LMEM_NULL);
  lmdense_mem = (LMDenseMem) lmem;

  *nhevalsD = nhe;

  return(LMDENSE_OKAY);
}

int LMDenseGetNumGradEvals(void *lmmem, int *ngevalsD)
{
  LMMem lm_mem;
  LMDenseMem lmdense_mem;

  /* Return immediately if lmmem is NULL */
  if (lmmem == NULL) return(LMDENSE_MEM_NULL);
  lm_mem = (LMMem) lmmem;

  if (lmem == NULL) return(LMDENSE_LMEM_NULL);
  lmdense_mem = (LMDenseMem) lmem;

  *ngevalsD = ngeD;

  return(LMDENSE_OKAY);
}

static int LMDenseInit(LMMem lm_mem)
{
  LMDenseMem lmdense_mem;

  lmdense_mem = (LMDenseMem) lmem;
  
  nhe   = 0;
  ngeD  = 0;
  
  if (hes == NULL) {
    hes = LMDenseDQHes;
    H_data = lm_mem;
  }

  return(0);
}

static int LMDenseSetup(LMMem lm_mem, booleantype same_p, int diff_nni,
                        booleantype forceSetup, booleantype *setupCurrent)
{
  long int ier;

  long int i,j;

  LMDenseMem lmdense_mem;
  
  lmdense_mem = (LMDenseMem) lmem;
 
  if (same_p) {
    /* use saved copy of H */
    DenseCopy(savedH, M);
  } else {
    /* call jac function for new H */
    nhe++;
    DenseZero(M); 
    hes(n, M, pp, H_data, vtemp1, vtemp2);
    DenseCopy(M, savedH);
  }
  
  /* Add mu*I */
  DenseScale(ONE/mu, M);
  DenseAddI(M);
  DenseScale(mu, M);

  /* Do LU factorization of M */
  ier = DenseFactor(M, pivots); 
  
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);

  *setupCurrent = TRUE;

  return(0);
}

static int LMDenseSolve(LMMem lm_mem, N_Vector xx, N_Vector bb)
{
  LMDenseMem lmdense_mem;
  realtype *bd;

  lmdense_mem = (LMDenseMem) lmem;
  
  bd = N_VGetArrayPointer(bb);

  DenseBacksolve(M, pivots, bd);

  N_VScale(ONE, bb, xx);

  return(0);
}

static int LMDenseFree(LMMem lm_mem)
{
  LMDenseMem  lmdense_mem;

  lmdense_mem = (LMDenseMem) lmem;
  
  DenseFreeMat(M);
  DenseFreeMat(savedH);
  DenseFreePiv(pivots);
  free(lmdense_mem);

  return(0);
}

static void LMDenseDQHes(long int Np, DenseMat H, N_Vector p, 
                         void *hes_data, N_Vector tmp1, N_Vector tmp2)
{
  realtype gnorm, minInc, inc, inc_inv, pjsaved, srur;
  realtype *tmp2_data, *p_data, *pscale_data;
  N_Vector gtemp, jthCol;
  long int j;

  LMMem lm_mem;
  LMDenseMem  lmdense_mem;

  /* hes_data points to lmmem */
  lm_mem = (LMMem) hes_data;
  lmdense_mem = (LMDenseMem) lmem;


  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  gtemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for p */
  pscale_data = N_VGetArrayPointer(pscale);
  p_data   = N_VGetArrayPointer(p);

  srur = RSqrt(uround);
  for (j = 0; j < Np; j++) {

    /* Generate the jth col of H(p) */

    N_VSetArrayPointer(DENSE_COL(H,j), jthCol);

    pjsaved = p_data[j];
    inc = MIN_INC_MULT*srur*ABS(pjsaved);

    p_data[j] += inc;
    grad(p, gtemp, f_data);
    p_data[j] = pjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, gtemp, -inc_inv, gval, jthCol);

    DENSE_COL(H,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  /* Increment counter ngeD */
  ngeD += Np;
}
