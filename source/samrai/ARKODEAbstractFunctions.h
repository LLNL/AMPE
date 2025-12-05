/*************************************************************************
 * Inspired by SAMRAI CVODEAbstractFunctions at
 * https://github.com/LLNL/SAMRAI
 * Adapted from PFiSM at https://github.com/ORNL/PFiSM
 ************************************************************************/

#ifndef included_ARKODEAbstractFunctions
#define included_ARKODEAbstractFunctions

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"


/**
 * Class ARKODEAbstractFunctions is an abstract base class that defines
 * an interface for the user-supplied RHSFunction and preconditioner
 * routines to be used with ARKODE and ARKSpgmr via the C++ wrapper
 * class ARKODESolver.  To use ARKODE with the C++ wrapper one must
 * derive a subclass of this base class and pass it into the ARKODESolver
 * constructor.  The pure virtual member functions in this interface are
 * used by ARKODE and ARKSpgmr during the ODE integration process.  The
 * complete argument lists in the function signatures defined by ARKODE
 * for the user-supplied routines have been preserved for the most part.
 * In a few cases, some arguments do not appear in the function signatures
 * below since they are superfluous via this interface.
 *
 * @see ARKODESolver
 * @see SundialsAbstractVector
 */

class ARKODEAbstractFunctions
{
 public:
   /**
    * The constructor and destructor for ARKODEAbstractFunctions
    * is empty.
    */
   ARKODEAbstractFunctions(){};
   virtual ~ARKODEAbstractFunctions(){};

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    *
    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int evaluateRHSFunction(
       double t, SAMRAI::solv::SundialsAbstractVector* y,
       SAMRAI::solv::SundialsAbstractVector* y_dot) = 0;

   /*
    * Implicit
    */
   virtual int evaluateRHSFunctionImp(
       double t, SAMRAI::solv::SundialsAbstractVector* y,
       SAMRAI::solv::SundialsAbstractVector* y_dot) = 0;

   /*
    * Explicit
    */
   virtual int evaluateRHSFunctionExp(
       double t, SAMRAI::solv::SundialsAbstractVector* y,
       SAMRAI::solv::SundialsAbstractVector* y_dot) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int ARKSpgmrPrecondSet(double t,
                                  SAMRAI::solv::SundialsAbstractVector* y,
                                  SAMRAI::solv::SundialsAbstractVector* fy,
                                  int jok, int* jcurPtr, double gamma) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int ARKSpgmrPrecondSolve(double t,
                                    SAMRAI::solv::SundialsAbstractVector* y,
                                    SAMRAI::solv::SundialsAbstractVector* fy,
                                    SAMRAI::solv::SundialsAbstractVector* r,
                                    SAMRAI::solv::SundialsAbstractVector* z,
                                    double gamma, double delta, int lr) = 0;
};

#endif
