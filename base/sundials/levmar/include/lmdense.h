#ifdef __cplusplus
extern "C" {
#endif
  
#ifndef _lmdense_h
#define _lmdense_h
  
#include "dense.h"
#include "nvector.h"
#include "sundialstypes.h"
  
  typedef void (*LMDenseHesFn)(long int N, DenseMat H, N_Vector pp, 
                               void *hes_data, N_Vector tmp1, N_Vector tmp2);
  
  int LMDense(void *lmmem, long int Np); 
  int LMDenseSetHesFn(void *lmmem, LMDenseHesFn dhes);
  int LMDenseSetHesData(void *lmmem, void *H_data);
  
  int LMDenseGetIntWorkSpace(void *lmmem, long int *leniwD);
  int LMDenseGetRealWorkSpace(void *lmmem, long int *lenrwD);
  int LMDenseGetNumHesEvals(void *lmmem, int *nhevalsD);
  int LMDenseGetNumGradEvals(void *lmmem, int *ngevalsD);
  
#define LMDENSE_SUCCESS     0
#define LMDENSE_OKAY        0
#define LMDENSE_MEM_NULL   -1 
#define LMDENSE_LMEM_NULL  -2 
#define LMDENSE_ILL_INPUT  -3
#define LMDENSE_MEM_FAIL   -4

#endif

#ifdef __cplusplus
}
#endif
