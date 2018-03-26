#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _lmspgmr_impl_h
#define _lmspgmr_impl_h

#include "levmar.h"
#include "spgmr.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "lmspgmr.h"

/*
 * -----------------------------------------------------------------
 * LMSpgmr solver constant
 * -----------------------------------------------------------------
 * LMSPGMR_MAXL : maximum dimension of Krylov subspace allowed by
 *                 default
 * -----------------------------------------------------------------
 */

#define LMSPGMR_MAXL 10

/*
 * -----------------------------------------------------------------
 * Types : struct LMSpgmrMemRec and struct *LMSpgmrMem
 * -----------------------------------------------------------------
 * A variable declaration of type struct *LMSpgmrMem denotes a
 * pointer to a data structure of type struct LMSpgmrMemRec. The
 * LMSpgmrMemRec structure contains fields that must be accessible
 * by LMSPGMR/SPGMR solver module routines.
 * -----------------------------------------------------------------
 */

typedef struct {

  /* problem specification data */

  int  g_maxl;          /* maximum allowable dimension of Krylov subspace      */     
  int  g_pretype;       /* preconditioning type: NONE, RIGHT, LEFT or BOTH
			   (used by SPGMR module and defined in
			   shared/include/iterative.h)                         */
  int  g_gstype;        /* Gram-Schmidt orthogonalization procedure:
			   CLASSICAL_GS or MODIFIED_GS (used by SPGMR module
			   and defined in shared/include/iterative.h)          */
  booleantype g_new_pp; /* flag indicating if the iterate has been updated -
			   Hessian must be updated/reevaluated (meant to be
			   used by user-supplied htimes function)              */
  int g_maxlrst;        /* maximum number of times the SPGMR linear solver
			   can be restarted                                    */
  int g_msbpsetup;      /* maximum number of iterations before forcing a call
                           to the preconditioner setup                         */

  /* counters */

  int g_nli;     /* number of linear iterations performed                 */
  int g_npe;     /* number of preconditioner evaluations                  */
  int g_nps;     /* number of times preconditioner was applied to linear
	            system                                                */
  int g_ncfl;    /* number of linear convergence failures                 */
  int g_ngeSG;   /* number of evaluations of the system function F(u) or
                    number of calls made to func routine                  */    
  int g_nhtimes; /* number of times the matrix-vector product H(u)*v
                    was computed or number of calls made to htimes
                    routine                                               */

  /* functions */

  LMSpgmrPrecSetupFn g_pset;     /* routine called to compute preconditioner
				     matrix                                    */
  LMSpgmrPrecSolveFn g_psolve;   /* subroutine called to solve a
				     preconditioned linear system              */ 
  LMSpgmrHesTimesVecFn g_htimes; /* function called to compute matrix-vector
				     product H(p)*v                            */

  /* memory references (pointers) */

  void *g_P_data; /* pointer to user-allocated memory block that is passed
		     to pset and psolve                                        */
  void *g_H_data; /* pointer to user-allocated memory block that is passed
		     to jtimes (only required if using a user-supplied
		     htimes routine)                                           */

  SpgmrMem g_spgmr_mem; /* pointer to SPGMR memory block (allocated by
			   SpgmrMalloc routine)                                */

} LMSpgmrMemRec, *LMSpgmrMem;

#endif

#ifdef __cplusplus
}
#endif
