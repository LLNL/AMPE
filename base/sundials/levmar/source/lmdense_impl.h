#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _lmdense_impl_h
#define _lmdense_impl_h

#include <stdio.h>

#include "lmdense.h"

#include "dense.h"
#include "nvector.h"
#include "sundialstypes.h"

typedef struct {

  long int d_n;       /* problem dimension                      */

  LMDenseHesFn d_hes; /* hes = Hessian routine to be called     */

  DenseMat d_M;       /* M = H + mu I                           */
  
  long int *d_pivots; /* pivots = pivot array for PM = LU       */
  
  DenseMat d_savedH;  /* savedH                                 */
  
  long int d_nhe;     /* nhe = no. of calls to hes              */

  long int d_ngeD;    /* nfeD = no. of calls to grad due to
                         difference quotient approximation of H */
  
  void *d_H_data;     /* H_data is passed to hes                */
  
} LMDenseMemRec, *LMDenseMem;

#endif

#ifdef __cplusplus
}
#endif
