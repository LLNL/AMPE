/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006/11/08 01:07:05 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban  @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation of cmputation of consistent initial conditions
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cpodes_private.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

static int cpicDoProjection(CPodeMem cp_mem);
static int cpicProjLinear(CPodeMem cp_mem);
static int cpicProjNonlinear(CPodeMem cp_mem);

static int cpicLinSolvDrv(CPodeMem cp_mem);
static int cpicLineSearch(CPodeMem cp_mem);

static int cpicImplComputeYp(CPodeMem cp_mem);

static int cpicFailFlag(CPodeMem, int flag);

/* 
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define ode_type (cp_mem->cp_ode_type)
#define tn             (cp_mem->cp_tn)
#define zn             (cp_mem->cp_zn) 
#define tempv          (cp_mem->cp_tempv)

#define linit         (cp_mem->cp_linit)
#define lsetup        (cp_mem->cp_lsetup)
#define lsolve        (cp_mem->cp_lsolve)
#define lsetup_exists (cp_mem->cp_lsetup_exists)

#define proj_enabled   (cp_mem->cp_proj_enabled)
#define proj_type      (cp_mem->cp_proj_type)
#define proj_norm      (cp_mem->cp_proj_norm)
#define cnstr_type     (cp_mem->cp_cnstr_type)
#define pfun           (cp_mem->cp_pfun)
#define p_data         (cp_mem->cp_p_data)
#define cfun           (cp_mem->cp_cfun)
#define c_data         (cp_mem->cp_c_data)

#define prjcoef  (cp_mem->cp_prjcoef)
#define acorP          (cp_mem->cp_acorP)
#define ctemp          (cp_mem->cp_ctemp)
#define tempvP1        (cp_mem->cp_tempvP1)
#define tempvP2        (cp_mem->cp_tempvP2)

#define linitP         (cp_mem->cp_linitP)
#define lsetupP        (cp_mem->cp_lsetupP)
#define lsolveP        (cp_mem->cp_lsolveP)
#define lsetupP_exists (cp_mem->cp_lsetupP_exists)

#define icprjcoef  (cp_mem->cp_icprjcoef)
#define callSetup  (cp_mem->cp_callSetup)
#define jacCurrent (cp_mem->cp_jacCurrent)
#define yy0        (cp_mem->cp_yy0)
#define yp0        (cp_mem->cp_yp0)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

/*
 * CPodeCalcIC 
 *
 * This calculates corrected initial conditions that are consistent 
 * with the invariant constraints and (for implicit-form ODEs) with 
 * the ODE system itself. It first projects the initial guess for 
 * the state vector (given by the user through CPodeInit or CPodeReInit) 
 * and then, if necessary, computes a state derivative vector as 
 * solution of F(t0, y0, y0') = 0.
 *
 * For nonlinear constraints, when using the CPODES internal projection
 * functions, it uses Newton iteration combined with a line search
 * algorithm.
 */
int CPodeCalcIC(void *cpode_mem)
{
  CPodeMem cp_mem;
  int flag;

  /* Check if cpode_mem exists */
  if (cpode_mem==NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeCalcIC", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* Check if cpode_mem was allocated */
  if (cp_mem->cp_MallocDone == FALSE) {
    cpProcessError(cp_mem, CP_NO_MALLOC, "CPODES", "CPodeCalcIC", MSGCP_NO_MALLOC);
    return(CP_NO_MALLOC);
  }

  /* If projection is enabled, check if appropriate projection functions are available */
  if (proj_enabled) {

    switch (proj_type) {

    case CP_PROJ_USER:

      /* Check if user provided the projection function */
      if ( pfun == NULL ) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_NO_PFUN);
        return(CP_ILL_INPUT);
      }

      break;

    case CP_PROJ_INTERNAL:

      /* Check if user provided the constraint function */
      if (cfun == NULL) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_NO_CFUN);
        return(CP_ILL_INPUT);
      }
      /* Check that lsolveP exists */ 
      if ( lsolveP == NULL) {
        cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_PLSOLVE_NULL);
        return(CP_ILL_INPUT);
      }
      /* Call linitP if it exists */
      if ( linitP != NULL ) {
        flag = linitP(cp_mem);
        if (flag != 0) {
          cpProcessError(cp_mem, CP_PLINIT_FAIL, "CPODES", "CPodeCalcIC", MSGCP_PLINIT_FAIL);
          return(CP_PLINIT_FAIL);
        }
      }

      break;

    } 

  }

  /* For explicit ODE, if projection is not enabled, there is nothing to do */
  if ( (ode_type == CP_EXPL) && (!proj_enabled) ) {
    cpProcessError(cp_mem, CP_WARNING, "CPODES", "CPodeCalcIC", MSGCP_IC_NO_WORK);
    return(CP_SUCCESS);
  }

  /* For implicit ODE, we will always need a linear solver */
  if (ode_type == CP_IMPL) {
    /* Check that lsolve exists */ 
    if (lsolve == NULL) {
      cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeCalcIC", MSGCP_LSOLVE_NULL);
      return(CP_ILL_INPUT);
    }
    /* Call linit if it exists */
    if ( linit != NULL ) {
      flag = linit(cp_mem);
      if (flag != 0) {
        cpProcessError(cp_mem, CP_LINIT_FAIL, "CPODES", "CPodeCalcIC", MSGCP_LINIT_FAIL);
        return(CP_LINIT_FAIL);
      }
    }
  }

  /* Allocate space and initialize temporary vectors */
  yy0 = N_VClone(tempv);
  N_VScale(ONE, zn[0], yy0);
  if (ode_type == CP_IMPL) {
    yp0 = N_VClone(tempv);
    N_VScale(ONE, zn[1], yy0);
  }

  /* Compute y consistent with constraints */
  if (proj_enabled)
    flag = cpicDoProjection(cp_mem);
  
  /* Compute yp consistent with DE */
  if ( (flag = CP_SUCCESS) && (ode_type == CP_IMPL) )
    flag = cpicImplComputeYp(cp_mem);
  
  /* If successful, load yy0 and yp0 into Nordsieck array */
  if (flag == CP_SUCCESS) {
    N_VScale(ONE, yy0, zn[0]);
    if (ode_type == CP_IMPL)  N_VScale(ONE, yp0, zn[1]);
  }

  /* Free temporary space */
  N_VDestroy(yy0);
  if (ode_type == CP_IMPL)  N_VDestroy(yp0);

  /* On any type of failure, print message and return proper flag */
  if (flag != CP_SUCCESS) {
    cpicFailFlag(cp_mem, flag);
    return(flag);
  }

  /* Otherwise return success flag */
  return(CP_SUCCESS);

}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/*
 * cpicFailFlag
 *
 * This function is called if the calculation of consistent IC failed.
 * Depending on the flag value, it calls the error handler and returns
 * the proper flag.
 */
static int cpicFailFlag(CPodeMem cp_mem, int flag)
{
  switch(flag) {
  case CP_FIRST_CNSTRFUNC_ERR:
    cpProcessError(cp_mem, CP_FIRST_CNSTRFUNC_ERR, "CPODES", "CPodeCalcIC", MSGCP_IC_CNSTRFUNC_FIRST);
    break;
  case CP_CNSTRFUNC_FAIL:
    cpProcessError(cp_mem, CP_CNSTRFUNC_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_CNSTRFUNC_FAILED);
    break;
  case CP_PROJFUNC_FAIL:
    cpProcessError(cp_mem, CP_PROJFUNC_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_PROJFUNC_FAILED);
    break;
  case CP_PLSETUP_FAIL:
    cpProcessError(cp_mem, CP_PLSETUP_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_PLSETUP_FAILED);
    break;
  case CP_PLSOLVE_FAIL:
    cpProcessError(cp_mem, CP_PLSOLVE_FAIL, "CPODES", "CPodeCalcIC", MSGCP_IC_PLSOLVE_FAILED);
    break;
  case CP_NO_RECOVERY:
    cpProcessError(cp_mem, CP_NO_RECOVERY, "CPODES", "CPodeCalcIC", MSGCP_IC_NO_RECOVERY);
    break;
  }
}

/* 
 * cpicDoProjection
 *
 * This is the top-level function handling the projection onto the
 * invariant manifold. It either calls the user-supplied projection
 * function or, depending on the type of constraints, one of the
 * internal projection functions (cpicProjLinear or cpicProjNonlinear).
 *
 * For user supplied projection function, use tempv to store the
 * corection due to projection, acorP (tempv is not touched until it 
 * is potentially used in cpCompleteStep).
 *
 * For the internal projection function, space for acorP was allocated 
 * in CPodeProjInit.
 */
static int cpicDoProjection(CPodeMem cp_mem)
{
  int flag = CP_SUCCESS, retval;

  switch (proj_type) {

  case CP_PROJ_INTERNAL:

    /* Evaluate constraints at initial time and with the provided yy0 */
    retval = cfun(tn, yy0, ctemp, c_data);
    if (retval != 0) return(CP_FIRST_CNSTRFUNC_ERR);

    /* Perform projection step 
     * On a successful return, yy0 was updated. */
    if (cnstr_type == CP_CNSTR_NONLIN) flag = cpicProjNonlinear(cp_mem);
    else                               flag = cpicProjLinear(cp_mem);

    break;

  case CP_PROJ_USER:

    acorP = tempv;
    
    /* Call the user projection function (with err=NULL) */
    retval = pfun(tn, yy0, acorP, prjcoef, NULL, p_data);
    if (retval != 0) return(CP_PROJFUNC_FAIL);

    /* Update yy0 */
    N_VLinearSum(ONE, yy0, ONE, acorP, yy0);

    break;

  }

  return(flag);

}

/*
 * cpicProjLinear
 *
 * This function performs the projection onto a linear invariant manifold
 */
static int cpicProjLinear(CPodeMem cp_mem)
{
  int retval;

  /* Call the lsetupP function to evaluate and factorize the 
   * Jacobian of constraints */
  retval = lsetupP(cp_mem, yy0, ctemp, tempvP1, tempvP2, tempv);
  if (retval != 0) return(CP_PLSETUP_FAIL);

  /* Call lsolveP (rhs is ctemp; solution in acorP) */
  retval = lsolveP(cp_mem, ctemp, acorP, yy0, ctemp, tempvP1, tempv);
  if (retval != 0) return(CP_PLSOLVE_FAIL);

  /* Update yy0 */
  N_VLinearSum(ONE, yy0, -ONE, acorP, yy0);

  return(CP_SUCCESS);
}

/*
 * cpicProjNonlinear
 *
 * This function performs the projection onto a nonlinear invariant manifold.
 * It implements Newton iteration with line search.
 */
static int cpicProjNonlinear(CPodeMem cp_mem)
{
  int flag;

  icprjcoef = PT01 * prjcoef;

  callSetup = TRUE;

  loop {

    /* calculate the Newton step */
    flag = cpicLinSolvDrv(cp_mem);
    if (flag != CP_SUCCESS) return(flag);

    /* perform line search to calculate an acceptable step */
    flag = cpicLineSearch(cp_mem);
    if (flag != CP_SUCCESS) return(flag);

    /* check if convergence was obtained */

  }


  return(CP_SUCCESS);
}


/*
 * cpicLinSolDrv
 *
 * This routine handles the process of solving for the approximate
 * solution of the Newton equations in the Newton iteration.
 * Subsequent routines handle the nonlinear aspects of its
 * application. 
 */
static int cpicLinSolvDrv(CPodeMem cp_mem)
{
  int retval;

  /* This loop is traversed at most two times.
   * The second pass happens only if the linear solve function
   * fails recoverably and we can attempt recovery by updating
   * the Jacobian.
   */
  loop {

    jacCurrent = FALSE;

    if (lsetupP_exists && callSetup) {
      retval = lsetupP(cp_mem, yy0, ctemp, tempvP1, tempvP2, tempv);
      if (retval != 0) return(CP_PLSETUP_FAIL);
      jacCurrent = TRUE;
    }

    /* solve for Newton step (rhs is ctemp; solution in acorP) */
    retval = lsolveP(cp_mem, ctemp, acorP, yy0, ctemp, tempvP1, tempv);
    /* If successful, return */
    if (retval == 0) return(CP_SUCCESS);
    /* If lsolveP failed unrecoverably, return */ 
    if (retval < 0)  return(CP_PLSOLVE_FAIL);
    /* If lsolveP failed recoverably, but we cannot revcover, return */
    if (!lsetupP_exists || jacCurrent) return(CP_NO_RECOVERY);

    /* Loop back and try again with updated Jacobian */
    callSetup = TRUE;

  }

  return(CP_SUCCESS);
}


static int cpicLineSearch(CPodeMem cp_mem)
{
  realtype lambda;

  lambda = ONE;

  /* Loop on linesearch variable lambda */
  loop {

    /* Get new */

  }

  return(CP_SUCCESS);
}

/*
 * cpicImplComputeYp
 *
 * For implicit-form ODEs, this function computes y' values consistent 
 * with the ODE, given y values consistent with the constraints.
 * In other words, it solves F(t0,y0,y0') = 0 for y0' using a Newton
 * iteration.
 */
static int cpicImplComputeYp(CPodeMem cp_mem)
{

  return(CP_SUCCESS);
}
