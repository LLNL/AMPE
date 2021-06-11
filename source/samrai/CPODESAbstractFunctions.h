/*************************************************************************
 *
 * This file is adapted from the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE at https://github.com/LLNL/SAMRAI.
 *
 * Copyright:     (c) 1997-2021 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Poisson equation
 *
 ************************************************************************/
#ifndef included_CPODESAbstractFunctions
#define included_CPODESAbstractFunctions

#include "SAMRAI/SAMRAI_config.h"
#include "SundialsAbstractVector.h"

using namespace SAMRAI;

/**
 * Class CPODESAbstractFunctions is an abstract base class that defines
 * an interface for the user-supplied RHSFunction and preconditioner
 * routines to be used with CPODES and CPSpgmr via the C++ wrapper
 * class CPODESSolver.  To use CPODES with the C++ wrapper one must
 * derive a subclass of this base class and pass it into the CPODESSolver
 * constructor.  The pure virtual member functions in this interface are
 * used by CPODES and CPSpgmr during the ODE integration process.  The
 * complete argument lists in the function signatures defined by CPODES
 * for the user-supplied routines have been preserved for the most part.
 * In a few cases, some arguments do not appear in the function signatures
 * below since they are superfluous via this interface.
 *
 * @see solv::CPODESSolver
 * @see SundialsAbstractVector
 */

class CPODESAbstractFunctions
{
 public:
   /**
    * The constructor and destructor for CPODESAbstractFunctions
    * is empty.
    */
   CPODESAbstractFunctions();
   virtual ~CPODESAbstractFunctions();

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    *


    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    *


    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int evaluateRHSFunction(double t, SundialsAbstractVector* y,
                                   SundialsAbstractVector* y_dot,
                                   int fd_flag) = 0;

   /**
    * User-supplied function to project the current solution and estimated
    * error onto the constraint manifold.
    *
    * The function arguments are:
    *


    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b corr     (OUTPUT){correction such that y+corr is on the constraint
    manifold}
    * - \b epsProj  (INPUT){WRMS norm tolerance for a nonlinear solver
    iteration}
    * - \b err      (INPUT/OUTPUT){projection of err constraint manifold}
    *


    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int applyProjection(double t, SundialsAbstractVector* y,
                               SundialsAbstractVector* corr, double epsProj,
                               SundialsAbstractVector* err) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CPSpgmrPrecondSet(double t, SundialsAbstractVector* y,
                                 SundialsAbstractVector* fy, int jok,
                                 int* jcurPtr, double gamma,
                                 SundialsAbstractVector* vtemp1,
                                 SundialsAbstractVector* vtemp2,
                                 SundialsAbstractVector* vtemp3) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CPSpgmrPrecondSolve(double t, SundialsAbstractVector* y,
                                   SundialsAbstractVector* fy,
                                   SundialsAbstractVector* r,
                                   SundialsAbstractVector* z, double gamma,
                                   double delta, int lr,
                                   SundialsAbstractVector* vtemp) = 0;
};

#endif
