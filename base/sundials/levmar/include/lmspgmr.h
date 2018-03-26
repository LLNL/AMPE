#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _lmspgmr_h
#define _lmspgmr_h

#include "levmar.h"
#include "spgmr.h"
#include "sundialstypes.h"
#include "nvector.h"
  
  typedef int (*LMSpgmrPrecSetupFn)(N_Vector pp, N_Vector pscale,
                                    N_Vector gval, N_Vector gscale,
                                    booleantype same_p, realtype mu, 
                                    void *P_data, 
                                    N_Vector vtemp1, N_Vector vtemp2);
  
  typedef int (*LMSpgmrPrecSolveFn)(N_Vector pp, N_Vector pscale, 
                                    N_Vector gval, N_Vector gscale, 
                                    N_Vector vv, realtype  mu,
                                    void *P_data, N_Vector vtemp);
  
  typedef int (*LMSpgmrHesTimesVecFn)(N_Vector v, N_Vector Hv,
                                      N_Vector pp, booleantype *new_pp, 
                                      void *H_data);
  
  int LMSpgmr(void *lmmem, int maxl);
  
  
  int LMSpgmrSetMaxRestarts(void *lmmem, int maxrs);
  int LMSpgmrSetPrecSetupFn(void *lmmem, LMSpgmrPrecSetupFn pset);
  int LMSpgmrSetPrecSolveFn(void *lmmem, LMSpgmrPrecSolveFn psolve);
  int LMSpgmrSetPrecData(void *lmmem, void *P_data);
  int LMSpgmrSetHesTimesVecFn(void *lmmem, LMSpgmrHesTimesVecFn htimes);
  int LMSpgmrSetHesData(void *lmmem, void *H_data);
  
  int LMSpgmrGetIntWorkSpace(void *lmmem, long int *leniwSG);
  int LMSpgmrGetRealWorkSpace(void *lmmem, long int *lenrwSG);
  int LMSpgmrGetNumPrecEvals(void *lmmem, int *npevals);
  int LMSpgmrGetNumPrecSolves(void *lmmem, int *npsolves);
  int LMSpgmrGetNumLinIters(void *lmmem, int *nliters);
  int LMSpgmrGetNumConvFails(void *lmmem, int *nlcfails);
  int LMSpgmrGetNumHtimesEvals(void *lmmem, int *nhvevals);
  int LMSpgmrGetNumGradEvals(void *lmmem, int *ngevalsSG); 
  
#define LMSPGMR_SUCCESS     0
#define LMSPGMR_OKAY        0
#define LMSPGMR_MEM_NULL   -1 
#define LMSPGMR_LMEM_NULL  -2 
#define LMSPGMR_ILL_INPUT  -3
#define LMSPGMR_MEM_FAIL   -4


#endif

#ifdef __cplusplus
}
#endif
