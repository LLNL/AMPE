/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   KINSOL solver for use within a SAMRAI-based application.
 *
 ************************************************************************/

#ifndef included_solv_KINSOL_SAMRAIContext
#define included_solv_KINSOL_SAMRAIContext

#include "SAMRAI/SAMRAI_config.h"

/*
 ************************************************************************
 *  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT KINSOL
 ************************************************************************
 */
#ifdef HAVE_SUNDIALS

#include "SAMRAI/solv/NonlinearSolverStrategy.h"
#include "SAMRAI/solv/KINSOLSolver.h"
#include "SAMRAI/solv/KINSOLAbstractFunctions.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"

namespace SAMRAI {
namespace solv {

/*!
 * @brief Wraps the KINSOLSolver C++ wrapper class so that
 * KINSOL may be used in applications that require a nonlinear solver.
 *
 * Class KINSOL_SAMRAIContext wraps the KINSOLSolver
 * C++ wrapper class so that KINSOL may be used in applications that
 * require a nonlinear solver.  The class KINSOLSolver does
 * not depend on SAMRAI.  Making it derived from the SAMRAI nonlinear solver
 * interface would require the KINSOL C++ wrapper to depend on SAMRAI.
 *
 * Important note:  This class can only create a KINSOL C++ wrapper instance,
 * initialize it in a rudimentary way, and invoke the solution process.
 * All other interaction with the nonlinear solver (i.e., setting parameters,
 * retrieving solver statistics, etc.) must be done directly with the
 * KINSOL wrapper accessible via the routine getKINSOLSolver().
 * Alternatively, solver parameters may be set up at initialization time
 * using the SAMRAI input database.
 *
 * If no parameters are read from input, KINSOL defaults are used.  See KINSOL
 * documentation for default information.  Optional input keys and types are:
 *
 *  - @b residual_stop_tolerance
 *      double value for stopping tolerance on norm of scaled residual.
 *
 *  - @b max_nonlinear_iterations
 *      integer value for maximum number of nonlinear iterations (MXITER).
 *
 *  - @b max_krylov_dimension
 *      integer value for maximum dimension of Krylov space.
 *
 *  - @b global_newton_strategy
 *      integer flag for globalization strategy.  Choices are
 *      "INEXACT_NEWTON" (default) and "LINESEARCH".
 *
 *  - @b max_newton_step
 *      double value for maximum allowable Newton step (MXNEWTONSTEP).
 *
 *  - @b nonlinear_step_tolerance
 *      double value for stopping tolerance on maximum entry in
 *      scaled Newton step.
 *
 *  - @b relative_function_error
 *      double value for relative error in function evaluation (RELFUNC).
 *
 *  - @b solution_update_constraint
 *      double value for constraint on relative change in solution (RELU).
 *
 *  - @b linear_convergence_test
 *      integer flag for linear solver convergence tolerance (ETACHOICE).
 *      Choices are "ETACONSTANT" (default), "ETACHOICE1", ETACHOICE2".
 *
 *  - @b max_sub_setup_calls
 *      number of nonlinear iterations between checks by the
 *      nonlinear residual monitoring algorithm (specifies lenght of
 *      subinterval) NOTE: should be a multiple of
 *      MaxStepsWithNoPrecondSetup
 *
 *  - @b residual_monitoring_params
 *      values of omega_min and omega_max scalars used by nonlinear
 *      residual monitoring algorithm.
 *      Default is [0.00001 and 0.9]
 *
 *  - @b residual_monitoring_constant
 *      constant value used by residual monitoring algorithm. If
 *      omega=0, then it is estimated using omega_min and
 *      omega_max.
 *      Default is 0.0.
 *
 *  - @b no_min_eps flag
 *      control whether or not the value * of eps is bounded below
 *      by 0.01*fnormtol. FALSE = "constrain value of eps by setting to
 *      the following: eps = MAX{0.01*fnormtol, eps}" TRUE = "do
 *      notconstrain value of eps".  Default is FALSE
 *
 *  - @b eisenstat_walker_params
 *      array of two double values Eisenstat-Walker choice 2; i.e.,
 *      the values are given as ETAALPHA, followed by ETAGAMMA.  Note: the
 *      values only apply when linear convergence test is set to ETACHOICE2.
 *
 *  - @b linear_solver_relative_tolerance
 *      double value for constant linear solver relative tolerance
 *      (ETACONST).  Note: value only apply when convergence test is
 *      set to ETACONSTANT.
 *
 *  - @b precond_setup_flag
 *      integer flag for preconditioner setup strategy (PRECOND_NO_INIT).
 *
 *  - @b max_solves_no_precond_setup
 *      integer value for number of nonlinear steps separating successive
 *      calls to preconditioner setup routine.
 *
 *  - @b max_linear_solve_restarts
 *      integer value for maximum number of linear solver restarts allowed.
 *
 *  - @b KINSOL_log_filename
 *      string value for name of KINSOL log file; default is "kinsol.log".
 *
 *  - @b KINSOL_print_flag
 *      integer flag for KINSOL log file print options (PRINTFL).
 *
 *  - @b uses_preconditioner
 *      boolean flag indicating whether a preconditioner is supplied.
 *      Default is false.
 *
 *  - @b uses_jac_times_vector
 *      boolean flag indicating whether an analytic Jacobian-vector
 *      product is supplied.  Default is false.
 *
 * Note that all input values may override values read in from restart.  If
 * no new input value is given, the restart value is used.
 *
 * A sample input file entry might look like:
 *
 * @code
 *   residual_stop_tolerance  =  10.e-6
 *   max_nonlinear_iterations =  200
 *   max_newton_step          =  0.1
 *   KINSOL_log_filename      =  "mylogfile"
 *   KINSOL_print_flag        =  3   // print all output KINSOL has to offer
 * @endcode
 *
 * @see solv::NonlinearSolverStrategy
 */
class KINSOL_SAMRAIContext:
   public NonlinearSolverStrategy,
   public tbox::Serializable
{
public:
   /**
    * Constructor for algs::KINSOL_SAMRAIContext allocates the KINSOL
    * C++ wrapper object and initializes rudimentary state associated
    * with user-supplied solver components.  Then, it reads solver parameter
    * from input and restart which may override default values.
    *
    * When assertion checking is active, an unrecoverable assertion
    * will result if the name string is empty or the pointer to the
    * user-defined KINSOL functions object is null.
    */
   KINSOL_SAMRAIContext(
      const std::string& object_name,
      const boost::shared_ptr<tbox::Database>& input_db,
      KINSOLAbstractFunctions* my_functions);

   /**
    * Destructor for algs::KINSOL_SAMRAIContext destroys the KINSOL
    * C++ wrapper object and the KINSOL solution vector wrapper.
    */
   ~KINSOL_SAMRAIContext();

   /**
    * Initialize the state of KINSOL based on vector argument representing
    * the solution of the nonlinear system.  In general, this routine must
    * be called before the solve() routine is invoked.
    */
   void
   initialize(
      const boost::shared_ptr<SAMRAIVectorReal<double> >& solution);

   /**
    * Solve the nonlinear problem and return and integer value defined by
    * KINSOL.  A return value of 1 indicates success (i.e., KINSOL_SUCCESS).
    * Consult the KINSOL documentation, KINSOL header file kinsol.h, or the
    * header file for the class KINSOLSolver for more information
    * about KINSOL return codes.  In general, the initialize() routine must
    * be called before this solve function to set up the solver.
    */
   int
   solve();

   /**
    * Return pointer to KINSOL solver C++ wrapper object.
    */
   KINSOLSolver *
   getKINSOLSolver()
   {
      return d_KINSOL_solver;
   }

   /**
    * Read input parameters from given database.
    *
    * When assertion checking is active, an unrecoverable assertion
    * will result if the database pointer is null.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& db);

   /**
    * Retrieve solver parameters from restart database matching object name.
    *
    * When assertion checking is active, an unrecoverable assertion
    * will result if a restart database matching the object name does not
    * exist, or if the class version number does not match that in restart.
    */
   void
   getFromRestart();

   /**
    * Retrieve solver parameters from restart database matching object name.
    *
    * When assertion checking is active, an unrecoverable assertion
    * will result if database pointer is null.
    */
   void
   putToDatabase(
      const boost::shared_ptr<tbox::Database>& db) const;

   /**
    * Print out all members of integrator instance to given output stream.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /**
    * Returns the object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int SOLV_KINSOL_SAMRAI_CONTEXT_VERSION;

   std::string d_object_name;

   /*
    * KINSOL nonlinear solver object and solution vector for nonlinear system.
    */

   KINSOLSolver* d_KINSOL_solver;
   SundialsAbstractVector* d_solution_vector;

   /*
    * KINSOL state data maintained here for input/restart capabilities.
    */

   double d_residual_stop_tolerance;
   int d_max_nonlinear_iterations;
   int d_max_krylov_dimension;
   int d_global_newton_strategy;
   double d_max_newton_step;
   double d_nonlinear_step_tolerance;
   double d_relative_function_error;
   double d_solution_update_constraint;
   int d_linear_convergence_test;
   int d_max_subsetup_calls;
   double d_residual_monitoring_params[2];
   double d_residual_monitoring_constant;
   double d_eisenstat_walker_params[2];
   double d_linear_solver_constant_tolerance;
   int d_precond_setup_flag;
   int d_max_solves_no_precond_setup;
   int d_max_linear_solve_restarts;
   std::string d_KINSOL_log_filename;
   int d_KINSOL_print_flag;
   bool d_no_min_eps;
   bool d_uses_preconditioner;
   bool d_uses_jac_times_vector;
};

}
}

#endif
#endif
