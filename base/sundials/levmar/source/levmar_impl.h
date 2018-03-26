#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _levmar_impl_h
#define _levmar_impl_h

#include "levmar.h"
#include "sundialstypes.h"
#include "nvector.h"


typedef struct LMMemRec {

  realtype lm_uround;

  /* problem specification data */

  CostFn lm_func;
  GradFn lm_grad;
  void *lm_f_data;

  realtype lm_sqrt_relerr;         /* relative error bound for grad(p)           */

  realtype  lm_eps;
  
  booleantype lm_setupNonNull;     /* flag indicating if lin. solv. setup
				      function is non-null                       */
  booleantype lm_setupCurrent;     /* flag indicating if the lin.solv. setup is
                                      current                                    */

  /* Stopping test values */
  realtype lm_mu0;                 /* initial mu value */
  realtype lm_eps1;
  realtype lm_eps2;
  int lm_maxit;

  realtype lm_mu;
  realtype lm_nu;

  realtype lm_F;
  realtype lm_grd_inf;

  /* Counters */

  int lm_nni;       /* number of nonlinear iterations */
  int lm_nfe;       /* number of calls made to func routine */
  int lm_nge;       /* number of calls made to grad routine */
  int lm_nnilsetup; /* value of nni at last lsetup call */
  int lm_msbpre;    /* allowed no. of steps without a call to prec */

  /* N_Vectors */

  N_Vector lm_pscale;
  N_Vector lm_gscale;

  N_Vector lm_pp;
  N_Vector lm_p_new;

  N_Vector lm_gval;

  N_Vector lm_h;

  N_Vector lm_vtemp1;
  N_Vector lm_vtemp2;

  /* space requirements for vector storage */ 

  long int lm_lrw1;
  long int lm_liw1;
  long int lm_lrw;
  long int lm_liw;

  /* linear solver data */
 
  /* function prototypes (pointers) */

  int (*lm_linit)(struct LMMemRec *lm_mem);

  int (*lm_lsetup)(struct LMMemRec *lm_mem, booleantype same_p, int diff_nni,
                   booleantype forceSetup, booleantype *setupCurrent);

  int (*lm_lsolve)(struct LMMemRec *lm_mem, N_Vector xx, N_Vector bb);

  int (*lm_lfree)(struct LMMemRec *lm_mem);

  void *lm_lmem;

  booleantype lm_MallocDone;

} *LMMem;


#endif

#ifdef __cplusplus
}
#endif
