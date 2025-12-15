/*************************************************************************
 * Inspired by SAMRAI CVODEAbstractFunctions at
 * https://github.com/LLNL/SAMRAI
 * Adapted from PFiSM at https://github.com/ORNL/PFiSM
 ************************************************************************/

#ifndef included_ARKODESolver
#define included_ARKODESolver

#include "SAMRAI/SAMRAI_config.h"

#include "ARKODEAbstractFunctions.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/Utilities.h"

extern "C" {
#include "sunlinsol/sunlinsol_spgmr.h"
}

// ARKODE includes
extern "C" {
#include "arkode/arkode.h"
#include "arkode/arkode_arkstep.h"
}


#include <string>

using namespace SAMRAI;

#ifndef LACKS_SSTREAM
#define ARKODE_ERROR(ierr)                                                 \
   do {                                                                    \
      if (ierr != ARK_SUCCESS) {                                           \
         std::ostringstream tboxos;                                        \
         tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__); \
      }                                                                    \
   } while (0)
#else
#define ARKODE_ERROR(ierr)                                         \
   do {                                                            \
      if (ierr != ARK_SUCCESS) {                                   \
         std::ostrstream tboxos;                                   \
         tbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__); \
      }                                                            \
   } while (0)
#endif

#define SUNABSVEC_CAST(v) \
   (static_cast<solv::SundialsAbstractVector*>(v->content))

/*!
 * @brief Class ARKODESolver serves as a C++ wrapper for the ARKODE
 * ordinary differential equation solver package.
 *
 * It is intended to be
 * sufficiently generic to be used independently of the SAMRAI framework.
 * This class declares one private static member function to link the
 * user-defined routine for right-hand side function evaluation and
 * two private statice member functions to link the user-defined
 * preconditioner setup and solve routines.  The implementation of these
 * functions is defined by the user in a subclass of the abstract base
 * class ARKODEAbstractFunctions.  The vector objects used within the
 * solver are given in a subclass of the abstract class
 * SundialsAbstractVector. The SundialsAbstractVector
 * class defines the vector kernel operations required by the ARKODE
 * package so that they may be easily supplied by a user who opts not
 * to use the vector kernel supplied by the ARKODE package.  (It should be
 * noted that the vector kernel used by ARKODE is the same as the one
 * used by the other packages in the Sundials of solvers).
 *
 * Note that this class provides no input or restart capabilities and
 * relies on ARKODE for output reporting.
 *
 * ARKODESolver Usage:
 *
 *    -  In order to use the ARKODESolver, the user must provide a
 *           concrete subclass of ARKODEAbstractFunctions abstract
 *           base class which defines the evaluateRHSFunction(),
 *           ARKSpgmrPrecondSet(), and ARKSpgmrPrecondSolve() methods.
 *
 *    -  Solving a system of ODEs using this ARKODE C++ interface
 *           requires four main stages.  First, a ARKODESolver
 *           object is created with a user-specified name and
 *           ARKODEAbstractFunctions object.  Second, the
 *           user must specify the integration parameters that s/he
 *           wishes to use.  Next, the user must call the ARKODESolver
 *           method initialize(solution_vector) with the
 *           SundialsAbstractVector that s/he wants to put the solution
 *           in.  Finally, the solve() method is invoked to solve the
 *           system of ODEs to the specified value of the independent
 *           variable.
 *
 *    -  The following is a list of integration parameters that
 *           must be specified by the user before calling the solve()
 *           method:
 *
 *            - Either relative or absolute tolerance must
 *                  be set - setRelativeTolerance(relative_tolerance),
 *                  setAbsoluteTolerance(absolute_tolerance)
 *
 *            - Initial value of independent variable -
 *                  setInitialValueOfIndependentVariable(init_time)
 *            - Final value of independent variable -
 *                  setFinalValueOfIndependentVariable(final_time
 *                      cvode_needs_initialization)
 *            - Initial condition vector -
 *                  setInitialConditionVector(ic_vector)
 *
 *    -  The following is a list of default values for integration
 *           parameters:
 *
 *           - @b Relative Tolerance
 *                0.0
 *
 *           - @b Scalar Absolute Tolerance
 *                0.0
 *
 *           - @b Vector Absolute Tolerance
 *                NULL
 *
 *           - @b Maximum Number of Internal Steps
 *                500
 *
 *           - @b Maximum Number of NIL Step Warnings
 *                10
 *
 *           - @b Initial Step Size
 *                determined by ARKODE
 *
 *           - @b Maximum Absolute Value of Step Size
 *                infinity
 *
 *           - @b Minimum Absolute Value of Step Size
 *                0.0
 *
 *           - @b ARKSpgmr Preconditioning Type
 *                NONE
 *
 *           - @b ARKSpgmr Gram Schmidt Algorithm
 *                MODIFIED_GS
 *
 *           - @b ARKSpgmr Maximum Krylov Dimension
 *                MIN(num_equations, ARKSPGMR_MAXL=5)
 *
 *           - @b ARKSpgmr Tolerance Scale Factor
 *                ARKSPGMR_DELT = 0.05.
 *
 * @see ARKODEAbstractFunctions
 * @see SundialsAbstractVector
 */

class ARKODESolver
{
 public:
   /**
    * Constructor for ARKODESolver sets default ARKODE parameters
    * and initializes the solver package with user-supplied functions
    * ARKODESolver parameters may be changed later using member
    * functions described below.
    *
    * Notes:
    *
    *        The solution vector is not passed into the constructor.
    *        Before the solver can be used, the initialize() function must
    *        be called.
    *
    * @pre !object_name.empty()
    * @pre my_functions != 0
    */
   ARKODESolver(const std::string& object_name,
                ARKODEAbstractFunctions* my_functions,
                const bool uses_preconditioner, const int im_ex);

   /**
    * Virtual destructor for ARKODESolver closes the
    * ARKODE log file and frees the memory allocated for the
    * ARKODE memory record.
    */
   virtual ~ARKODESolver();

   /**
    * Initialize solver with solution vector.  The solution vector is
    * required to initialize the memory record used internally within
    * ARKODE.  This routine must be called before the solver can be used.
    *
    * @pre solution != 0
    * @pre d_solution_vector == 0
    */
   void initialize(solv::SundialsAbstractVector* solution)
   {
      TBOX_ASSERT(solution != 0);
      TBOX_ASSERT(d_solution_vector == 0);
      d_solution_vector = solution;
      d_ARKODE_needs_initialization = true;
      initializeARKODE();
   }

   /**
    * Integrate ODE system specified t_f.  The integer return value is
    * a termination code defined by ARKODE.  The following is a table
    * of termination codes and a brief description of their meanings.
    *
    * If ARKODE or ARKSpgmr requires re-initialization, it is
    * automatically done before the solve.  This may be required if any
    * of the ARKODE or ARKSpgmr data parameters have changed since the
    * last call to the solver.
    *
    * @pre d_user_t_f > d_t_0
    */
   int solve()
   {
      initializeARKODE();

      /*
       * Check to make sure that user specified final value for t
       * is greater than initial value for t.
       */
      TBOX_ASSERT(d_user_t_f > d_t_0);

      /*
       * See cvode.h header file for definition of return types.
       */
      int retval = ARKStepEvolve(d_arkode_mem, d_user_t_f,
                                 d_solution_vector->getNVector(), &d_actual_t_f,
                                 d_stepping_method);
      return retval;
   }

   /**
    * Accessor function for setting ARKODE output log file name and output
    * printing options.  Output file name and options may be changed
    * throughout run as desired.
    */
   void setLogFileData(const std::string& log_fname = std::string())
   {
      if (log_fname != d_arkode_log_file_name) {
         if (log_fname.empty()) {
            d_arkode_log_file_name = "arkode.log";
         } else {
            d_arkode_log_file_name = log_fname;
         }
         d_ARKODE_needs_initialization = true;
      }
   }

   /**
    * Set ARKODESolver to use my_functions as the concrete subclass
    * of the ARKODEAbstractFunctions class that defines the
    * right-hand side evaluation and preconditioner functions.  The
    * uses_preconditioner argument indicates whether or not the
    * the user has defined preconditioner routines in their concrete
    * subclass of the ARKODEAbstractFunctions class.
    *
    * @pre my_functions != 0
    */
   void setARKODEFunctions(ARKODEAbstractFunctions* my_functions,
                           const bool uses_preconditioner)
   {
      TBOX_ASSERT(my_functions != 0);
      d_arkode_functions = my_functions;
      d_uses_preconditioner = uses_preconditioner;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Return pointer to object that provides user-defined functions for
    * ARKODE and ARKSpgmr.
    */
   ARKODEAbstractFunctions* getARKODEFunctions() const
   {
      return d_arkode_functions;
   }

   /**
    * Set the relative tolerance level.
    *
    * Note that pure absolute tolerance can be used by
    * setting the relative tolerance to 0.  However,
    * it is an error to simultaneously set relative and
    * absolute tolerances to 0.
    *
    * @pre relative_tolerance >= 0.0
    */
   void setRelativeTolerance(double relative_tolerance)
   {
      TBOX_ASSERT(relative_tolerance >= 0.0);
      d_relative_tolerance = relative_tolerance;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set the scalar absolute tolerance level.
    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However,
    * it is an error to simultaneously set relative and
    * absolute tolerances to 0.
    *
    * @pre absolute_tolerance >= 0.0
    */
   void setAbsoluteTolerance(double absolute_tolerance)
   {
      TBOX_ASSERT(absolute_tolerance >= 0.0);
      d_absolute_tolerance_scalar = absolute_tolerance;
      d_use_scalar_absolute_tolerance = true;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set the vector absolute tolerance level.
    *
    * Note that pure relative tolerance can be used by
    * setting the absolute tolerance to 0.  However,
    * it is an error to simultaneously set relative and
    * absolute tolerances to 0.
    *
    * @pre absolute_tolerance != 0
    * @pre absolute_tolerance->vecMin() >= 0.0
    */
   void setAbsoluteTolerance(solv::SundialsAbstractVector* absolute_tolerance)
   {
      TBOX_ASSERT(absolute_tolerance != 0);
      TBOX_ASSERT(absolute_tolerance->vecMin() >= 0.0);
      d_absolute_tolerance_vector = absolute_tolerance;
      d_use_scalar_absolute_tolerance = false;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set stepping method to use for integration.  There are
    * stepping methods: NORMAL and ONE_STEP.  The NORMAL
    * method has the solver take internal steps until
    * it has reached or just passed the user specified t_f
    * parameter. The solver then interpolates in order to
    * return an approximate value of y(t_f). The ONE_STEP
    * option tells the solver to just take one internal step
    * and return the solution at the point reached by that
    * step.
    *
    * Note: the enumeration constants NORMAL and ONE_STEP are
    * defined in cvode.h.
    */
   void setSteppingMethod(int stepping_method)
   {
      TBOX_ASSERT((stepping_method == ARK_NORMAL) ||
                  (stepping_method == ARK_ONE_STEP));
      d_stepping_method = stepping_method;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set initial value for independent variable.
    */
   void setInitialValueOfIndependentVariable(double t_0)
   {
      d_t_0 = t_0;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set final value for independent variable (i.e. the value of
    * independent variable to integrate the system to).  The boolean
    * argument specifies whether ARKODE should be re-initialized (i.e.
    * on first step) or if we are taking subsequent steps in a
    * sequence, in which case it is not initialized.
    */
   void setFinalValueOfIndependentVariable(double t_f,
                                           bool arkode_needs_initialization)
   {
      d_user_t_f = t_f;
      d_ARKODE_needs_initialization = arkode_needs_initialization;
   }

   /**
    * Set initial condition vector.
    *
    * @pre ic_vector != 0
    */
   void setInitialConditionVector(solv::SundialsAbstractVector* ic_vector)
   {
      TBOX_ASSERT(ic_vector != 0);
      d_ic_vector = ic_vector;
      d_ARKODE_needs_initialization = true;
   }

   void setMaximumNumberOfNonlinIters(int nls_max_iter)
   {
      TBOX_ASSERT(nls_max_iter >= 0);
      d_nls_max_iter = nls_max_iter;
      d_ARKODE_needs_initialization = true;
   }

   void setMethodOrder(int arkode_order)
   {
      TBOX_ASSERT(arkode_order >= 0);
      d_arkode_order = arkode_order;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set maximum number of warning messages issued by the solver
    * that (t + h == t) on the next internal step.  By default,
    * this is set to 10.
    *
    * @pre max_num_warnings >= 0
    */
   void setMaximumNumberOfNilStepWarnings(int max_num_warnings)
   {
      TBOX_ASSERT(max_num_warnings >= 0);
      d_max_num_warnings = max_num_warnings;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set initial step size.
    *
    * @pre init_step_size >= 0.0
    */
   void setInitialStepSize(double init_step_size)
   {
      TBOX_ASSERT(init_step_size >= 0.0);
      d_init_step_size = init_step_size;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set maximum absolute value of step size allowed.
    * By default, there is no upper bound on the absolute value
    * of step size.
    *
    * @pre max_step_size >= 0.0
    */
   void setMaximumAbsoluteStepSize(double max_step_size)
   {
      TBOX_ASSERT(max_step_size >= 0.0);
      d_max_step_size = max_step_size;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set minimum absolute value of step size allowed.
    * By default, this is set to 0.0.
    *
    * @pre min_step_size >= 0.0
    */
   void setMinimumAbsoluteStepSize(double min_step_size)
   {
      TBOX_ASSERT(min_step_size >= 0.0);
      d_min_step_size = min_step_size;
      d_ARKODE_needs_initialization = true;
   }

   // Methods for setting ARKSpgmr parameters.

   /**
    * Set the preconditioning type to be used by ARKSpgmr.
    * This must be one of the four enumeration constants
    * NONE, LEFT, RIGHT, or BOTH defined in iterativ.h.
    * These correspond to no preconditioning, left preconditioning only,
    * right preconditioning only, and both left and right
    * preconditioning, respectively.
    *
    * @pre (precondition_type == PREC_NONE) ||
    *      (precondition_type == PREC_LEFT) ||
    *      (precondition_type == PREC_RIGHT) ||
    *      (precondition_type == PREC_BOTH)
    */
   void setPreconditioningType(int precondition_type)
   {
      TBOX_ASSERT((precondition_type == PREC_NONE) ||
                  (precondition_type == PREC_LEFT) ||
                  (precondition_type == PREC_RIGHT) ||
                  (precondition_type == PREC_BOTH));
      d_precondition_type = precondition_type;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set the Gram-Schmidt orthogonalization type to be used by ARKSpgmr.
    * This must be one of the two enumeration constants MODIFIED_GS
    * or CLASSICAL_GS defined in iterativ.h. These correspond to
    * using modified Gram-Schmidt and classical Gram-Schmidt, respectively.
    *
    * @pre (gs_type == CLASSICAL_GS) || (gs_type == MODIFIED_GS)
    */
   void setGramSchmidtType(int gs_type)
   {
      TBOX_ASSERT((gs_type == CLASSICAL_GS) || (gs_type == MODIFIED_GS));
      d_gram_schmidt_type = gs_type;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set the maximum Krylov dimension to be used by ARKSpgmr.
    * This is an optional input to the ARKSPGMR solver. Pass 0 to
    * use the default value MIN(num_equations, ARKSPGMR_MAXL=5).
    *
    * @pre max_krylov_dim >= 0
    */
   void setMaxKrylovDimension(int max_krylov_dim)
   {
      TBOX_ASSERT(max_krylov_dim >= 0);
      d_max_krylov_dim = max_krylov_dim;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Set the factor by which the tolerance on the nonlinear
    * iteration is multiplied to get a tolerance on the linear iteration.
    * This is an optional input to the ARKSPGMR solver. Pass 0 to
    * use the default value ARKSPGMR_DELT = 0.05.
    *
    * @pre tol_scale_factor >= 0
    */
   void setARKSpgmrToleranceScaleFactor(double tol_scale_factor)
   {
      TBOX_ASSERT(tol_scale_factor >= 0);
      d_tol_scale_factor = tol_scale_factor;
      d_ARKODE_needs_initialization = true;
   }

   /**
    * Get solution vector.
    */
   solv::SundialsAbstractVector* getSolutionVector() const
   {
      return d_solution_vector;
   }

   /**
    * Get k-th derivative vector at the specified value of the
    * independent variable, t.  The integer return value is
    * return code the ARKODE CVodeDky() function.  The following is a table
    * of termination codes and a brief description of their meanings.
    *
    * CVodeDky Return Codes:
    *
    *    - @b OKAY (=0)
    *        CVodeDky succeeded.
    *
    *    - @b BAD_K (=-1)
    *
    *    - @b BAD_T (=-2)
    *
    *    - @b BAD_DKY (=-3)
    *
    *    - @b DKY_NO_MEM (=-4)
    *
    * Important Notes:
    *    -
    *       t must lie in the interval [t_cur - h, t_cur]
    *       where t_cur is the current internal time reached
    *       and h is the last internal step size successfully
    *       used by the solver.
    *
    *    -
    *       k may take on value 0, 1, . . . q where q is the order
    *       of the current linear multistep method being used.
    *
    *    -
    *       the dky vector must be allocated by the user.
    *
    *    -
    *       it is only leagal to call this method after a
    *       successful return from the solve() method.
    *
    */
   int getDkyVector(double t, int k, solv::SundialsAbstractVector* dky) const
   {
      int return_code = ARKStepGetDky(d_arkode_mem, t, k, dky->getNVector());
      return return_code;
   }

   /**
    * Get actual value of the independent variable that ARKODE integrated
    * to (i.e. the value of t that actually corresponds to the solution
    * vector y).
    */
   double getActualFinalValueOfIndependentVariable() const
   {
      return d_actual_t_f;
   }

   /**
    * Print ARKODE and ARKSpgmr statistics.
    */
   void printStatistics(std::ostream& os) const
   {
      printARKODEStatistics(os);
      if (d_im_ex > 0) printARKSpgmrStatistics(os);
   }

   /**
    * Print ARKODE statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    *
    *    - @b lenrw
    *       size (in double words) of memory used for doubles
    *
    *    - @b leniw
    *       size (in integer words) of memory used for integers
    *
    *    - @b nst
    *       cumulative number of internal steps taken by solver
    *
    *    - @b nfe
    *       number of right-hand side function evaluations
    *
    *    - @b nni
    *       number of NEWTON iterations performed
    *
    *    - @b nsetups
    *       number of calls made to linear solver's setup routine
    *
    *    - @b netf
    *       number of local error test failures
    *
    *    - @b ncfn
    *       number of nonlinear convergence failures
    *
    *    - @b qu
    *       order used during the last internal step
    *
    *    - @b qcur
    *       order to be used on the next internal step
    *
    *    - @b hu
    *       step size for the last internal step
    *
    *    - @b hcur
    *       step size to be attempted on the next internal step
    *
    *    - @b tcur
    *       current internal value of t reached by the solver
    *
    *    - @b tolsf
    *       suggested tolerance scaling factor
    */
   void printARKODEStatistics(std::ostream& os) const;

   // ARKODE optional return values.

   /**
    * Return the cumulative number of stability-limited steps taken by
    * the solver.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfExpStepsTaken() const
   {
      long int r;
      int ierr = ARKStepGetNumExpSteps(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of calls to the right-hand side function.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfRHSFunctionExCalls() const
   {
      long int nfe_evals;
      long int nfi_evals;
      int ierr = ARKStepGetNumRhsEvals(d_arkode_mem, &nfe_evals, &nfi_evals);
      ARKODE_ERROR(ierr);
      return static_cast<int>(nfe_evals);
   }

   int getNumberOfRHSFunctionImpCalls() const
   {
      long int nfe_evals;
      long int nfi_evals;
      int ierr = ARKStepGetNumRhsEvals(d_arkode_mem, &nfe_evals, &nfi_evals);
      ARKODE_ERROR(ierr);
      return static_cast<int>(nfi_evals);
   }

   /**
    * Return the number of calls made to linear solver setup
    * routines.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfLinearSolverSetupCalls() const
   {
      long int r;
      int ierr = ARKStepGetNumLinSolvSetups(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of NEWTON iterations performed.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfNonlinIters() const
   {
      long int r;
      int ierr = ARKStepGetNumNonlinSolvIters(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of nonlinear convergence failures that have
    * occurred.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfNonlinConvFails() const
   {
      long int r;
      int ierr = ARKStepGetNumNonlinSolvConvFails(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of local error test failures.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfLocalErrorTestFailures() const
   {
      long int r;
      int ierr = ARKStepGetNumErrTestFails(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the order of the linear multistep method used during
    * the last internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfSteps() const
   {
      long int r;
      int ierr = ARKStepGetNumSteps(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Returns the cumulative number of steps attempted by the solver
    * the next internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   int getNumberOfStepAttempts() const
   {
      long int r;
      int ierr = ARKStepGetNumStepAttempts(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the step size for the last internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   double getStepSizeForLastInternalStep() const
   {
      realtype r;
      int ierr = ARKStepGetLastStep(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return r;
   }

   /**
    * Return the step size to be used in the next internal step.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   double getStepSizeForNextInternalStep() const
   {
      realtype r;
      int ierr = ARKStepGetCurrentStep(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return r;
   }

   /**
    * Return the current internal value of the independent
    * variable reached by the solver.
    *
    * Note: if the solver was not set to collect statistics,
    * the minimum double value (as defined in float.h) is
    * returned.
    */
   double getCurrentInternalValueOfIndependentVariable() const
   {
      realtype r;
      int ierr = ARKStepGetTolScaleFactor(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return r;
   }

   /**
    * Return the suggested tolerance scaling factor.
    *
    * Note: if the solver was not set to collect statistics,
    * a value of -1 is returned.
    */
   double getARKODESuggestedToleranceScalingFactor() const
   {
      realtype r;
      int ierr = ARKStepGetTolScaleFactor(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return r;
   }

   // ARKSpgmr optional return values.

   /**
    * Print ARKSpgmr statistics to the stream.
    *
    * The abbreviations printed out refer to the following
    * quantities:
    *
    *    - @b spgmr_lrw
    *      size (in double words) of memory used for doubles
    *
    *    - @b spgmr_liw
    *       size (in integer words) of memory used for integers
    *
    *    - @b nli
    *       number of linear iterations
    *
    *    - @b ncfl
    *       number of linear convergence failures
    *
    *    - @b npe
    *       number of preconditioner evaluations
    *
    *    - @b nps
    *       number of calls to ARKSpgmrPrecondSolve()
    */
   void printARKSpgmrStatistics(std::ostream& os) const;

   /**
    * Return the number of preconditioner evaluations.
    */
   int getNumberOfPreconditionerEvaluations() const
   {
      long int r;
      int ierr = ARKStepGetNumPrecEvals(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of linear iterations.
    */
   int getNumberOfLinearIterations() const
   {
      long int r;
      int ierr = ARKStepGetNumLinIters(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of ARKSpgmrPrecondSolve() calls.
    */
   int getNumberOfPrecondSolveCalls() const
   {
      long int r;
      int ierr = ARKStepGetNumPrecSolves(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Return the number of linear convergence failures.
    */
   int getNumberOfLinearConvergenceFailures() const
   {
      long int r;
      int ierr = ARKStepGetNumLinConvFails(d_arkode_mem, &r);
      ARKODE_ERROR(ierr);
      return static_cast<int>(r);
   }

   /**
    * Print out all data members for this object.
    */
   virtual void printClassData(std::ostream& os) const;

   /**
    * Returns the object name.
    */
   const std::string& getObjectName() const { return d_object_name; }

 private:
   /*
    * Static member function for linkage with ARKODE routines.
    */
   static int ARKODERHSFuncEval(realtype t, N_Vector y, N_Vector y_dot,
                                void* my_solver)
   {
      return ((ARKODESolver*)my_solver)
          ->getARKODEFunctions()
          ->evaluateRHSFunction(t, SUNABSVEC_CAST(y), SUNABSVEC_CAST(y_dot));
   }

   /*
    * Static member function for linkage with ARKODE routines.
    */
   static int ARKODERHSFuncEvalImp(realtype t, N_Vector y, N_Vector y_dot,
                                   void* my_solver)
   {
      return ((ARKODESolver*)my_solver)
          ->getARKODEFunctions()
          ->evaluateRHSFunctionImp(t, SUNABSVEC_CAST(y), SUNABSVEC_CAST(y_dot));
   }

   /*
    * Static member function for linkage with ARKODE routines.
    */
   static int ARKODERHSFuncEvalExp(realtype t, N_Vector y, N_Vector y_dot,
                                   void* my_solver)
   {
      return ((ARKODESolver*)my_solver)
          ->getARKODEFunctions()
          ->evaluateRHSFunctionExp(t, SUNABSVEC_CAST(y), SUNABSVEC_CAST(y_dot));
   }

   /*
    * Static member functions for linkage with ARKSpgmr routines.
    */
   static int ARKSpgmrPrecondSet(realtype t, N_Vector y, N_Vector fy,
                                 booleantype jok, booleantype* jcurPtr,
                                 realtype gamma, void* my_solver)
   {
      int success =
          ((ARKODESolver*)my_solver)
              ->getARKODEFunctions()
              ->ARKSpgmrPrecondSet(t, SUNABSVEC_CAST(y), SUNABSVEC_CAST(fy),
                                   jok, jcurPtr, gamma);
      return success;
   }

   static int ARKSpgmrPrecondSolve(realtype t, N_Vector y, N_Vector fy,
                                   N_Vector r, N_Vector z, realtype gamma,
                                   realtype delta, int lr, void* my_solver)
   {
      int success =
          ((ARKODESolver*)my_solver)
              ->getARKODEFunctions()
              ->ARKSpgmrPrecondSolve(t, SUNABSVEC_CAST(y), SUNABSVEC_CAST(fy),
                                     SUNABSVEC_CAST(r), SUNABSVEC_CAST(z),
                                     gamma, delta, lr);
      return success;
   }

   /*
    * Open ARKODE log file, allocate main memory for ARKODE and initialize
    * ARKODE memory record.  ARKODE is initialized based on current state
    * of solver parameter data members.  If any solver parameters have
    * changed since last initialization, this function will be automatically
    * invoked at next call to the solve() method.  Also, if NEWTON iteration
    * is specified, this method also initializes the ARKSpgmr linear solver.
    *
    * @pre d_solution_vector != 0
    */
   void initializeARKODE();

   std::string d_object_name;

   /*
    * The following data members are input or set to default values in
    * the ARKODESolver constructor.  Many of these can be altered at
    * any time through class member functions.  When this occurs,
    * ARKODE may need to be re-initialized (e.g., if the linear solver
    * changes, ARKODE must change its memory record).  In this case,
    * the initializeARKODE() member function is invoked in the next
    * call to solve().
    */

   /*
    * Solution vector.
    */
   solv::SundialsAbstractVector* d_solution_vector;

   /*
    * Pointer to object which provides user-supplied functions to ARKODE
    * and ARKSpgmr.
    */
   ARKODEAbstractFunctions* d_arkode_functions;

   /*
    * ARKODE memory record.
    */
   void* d_arkode_mem;  // ARKODE memory structure

   /*
    * Linear solver for preconditioning
    */
   SUNLinearSolver d_linear_solver;


   /*
    * ARKODE log file information.
    */
   FILE* d_arkode_log_file;             // ARKODE message log file
   std::string d_arkode_log_file_name;  // ARKODE log file name

   /*
    * ODE parameters.
    */
   double d_t_0;         // initial value for independent variable
   double d_user_t_f;    // user-specified final value for independent variable
   double d_actual_t_f;  // actual final value of indep. variable after a step
   solv::SundialsAbstractVector* d_ic_vector;

   /*
    * ODE integration parameters.
    */
   int d_linear_multistep_method;
   double d_relative_tolerance;
   bool d_use_scalar_absolute_tolerance;
   double d_absolute_tolerance_scalar;
   solv::SundialsAbstractVector* d_absolute_tolerance_vector;
   int d_stepping_method;

   /*
    * Optional ARKODE parameters.
    */
   int d_arkode_order;
   int d_im_ex;
   int d_nls_max_iter;
   int d_max_num_warnings;
   double d_init_step_size;
   double d_max_step_size;
   double d_min_step_size;
   /*
    * ARKSpgmr parameters
    */
   int d_precondition_type;
   int d_gram_schmidt_type;
   int d_max_krylov_dim;
   double d_tol_scale_factor;

   /*
    * Boolean flag indicating whether ARKODE needs initialization
    * when solver is called.
    */
   bool d_ARKODE_needs_initialization;

   /*
    * Boolean flag indicating whether user-supplied preconditioner
    * routines are provided in the concrete subclass of
    * ARKODEAbstractFunctions.
    */
   bool d_uses_preconditioner;
};

#endif
