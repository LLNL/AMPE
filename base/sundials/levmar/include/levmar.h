#ifdef __cplusplus
extern "C" {
#endif
  
#ifndef _levmar_h_
#define _levmar_h_
  
#include "sundialstypes.h"
#include "nvector.h"
  
#define LM_MAX_ITERS     100
#define LM_INIT_MU    	 1.0E0
#define LM_STOP_EPS1	 1.0E-04
#define LM_STOP_EPS2     1.0E-08

#define LM_MSBPRE        10  
  
  typedef realtype (*CostFn)(N_Vector p, void *f_data );
  
  typedef void (*GradFn)(N_Vector p, N_Vector grd, void *f_data );
  
  void *LMCreate(void);
  
  int LMSetFdata(void *lmmem, void *f_data);
  int LMSetInitMu(void *lmmem, realtype mu0);
  int LMSetEps1(void *lmmem, realtype eps1);
  int LMSetEps2(void *lmmem, realtype eps2);
  int LMSetMaxIters(void *lmmem, int maxit);
  
  enum { SUCCESS = 0, LMS_NO_MEM = -1, LMS_ILL_INPUT = -2 };
  
  int LMMalloc(void *lmmem, CostFn func, GradFn grad, N_Vector tmpl);
  
  enum { LMM_NO_MEM = -1, LMM_ILL_INPUT = -2, LMM_MEM_FAIL = -3 };

  int LMSolve(void *lmmem, N_Vector p, N_Vector p_scale, N_Vector g_scale);

  int LMGetIntWorkSpace(void *lmmem, long int *leniw);
  int LMGetRealWorkSpace(void *lmmem, long int *lenrw);
  int LMGetNumNonlinSolvIters(void *lmmem, int *nniters);
  int LMGetNumFuncEvals(void *lmmem, int *nfevals);
  int LMGetNumGradEvals(void *lmmem, int *ngevals);
  int LMGetCost(void *lmmem, realtype *cfval);
  int LMGetGradNorm(void *lmmem, realtype *grad_inf_norm);

  void LMFree(void *lmmem);

#define LM_SUCCESS        0
#define LM_OKAY           0
#define LM_SMALL_GRD      1
#define LM_SMALL_H        2 
#define LM_MEM_NULL      -1
#define LM_ILL_INPUT     -2
#define LM_MEM_FAIL      -3  /* memory allocation failed */
#define LM_ALMOST_SING   -4 
#define LM_REACHED_MAXIT -5
#define LM_LINIT_FAIL    -6  /* lin. solv. init failed */
#define LM_LSETUP_FAIL   -7  /* lin. solv. setup failed  */
#define LM_LSOLVE_FAIL   -8  /* lin. solv. solution failed */
#define LM_KRYLOV_FAIL   -9  /* Krylov failed */
  
#endif
  
#ifdef __cplusplus
}
#endif

