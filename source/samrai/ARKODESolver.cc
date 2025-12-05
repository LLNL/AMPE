/*************************************************************************
 * Inspired by SAMRAI CVODEAbstractFunctions at
 * https://github.com/LLNL/SAMRAI
 * Adapted from PFiSM at https://github.com/ORNL/PFiSM
 ************************************************************************/
#include "ARKODESolver.h"

/*
 * ARKODESolver constructor and destructor.
 */
ARKODESolver::ARKODESolver(const std::string& object_name,
                           ARKODEAbstractFunctions* my_functions,
                           const bool uses_preconditioner, const int im_ex)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(my_functions != 0);

   d_object_name = object_name;
   d_arkode_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_im_ex = im_ex;

   d_solution_vector = 0;

   /*
    * Set default parameters to safe values or to ARKODE/ARKSpgmr defaults.
    */

   /*
    * ARKODE memory record and log file.
    */
   d_arkode_mem = 0;
   d_linear_solver = 0;
   d_arkode_log_file = 0;
   d_arkode_log_file_name = "arkode.log";

   /*
    * ODE parameters.
    */
   d_t_0 = 0.0;
   d_user_t_f = 0.0;
   d_actual_t_f = 0.0;
   d_ic_vector = 0;

   /*
    * ODE integration parameters.
    */
   setRelativeTolerance(0.0);
   setAbsoluteTolerance(0.0);
   d_absolute_tolerance_vector = 0;
   setSteppingMethod(ARK_NORMAL);


   d_nls_max_iter = -1;
   d_max_num_warnings = -1;
   d_init_step_size = -1;
   d_max_step_size = -1;
   d_min_step_size = -1;
   d_arkode_order = -1;

   /*
    * ARKSpgmr parameters.
    *
    * Note that when the maximum krylov dimension and ARKSpgmr
    * tolerance scale factor are set to 0, ARKSpgmr uses its
    * internal default values.  These are described in the header for
    * this class.
    */
   setPreconditioningType(PREC_NONE);
   setGramSchmidtType(MODIFIED_GS);
   setMaxKrylovDimension(0);
   setARKSpgmrToleranceScaleFactor(0);

   d_ARKODE_needs_initialization = true;
}

ARKODESolver::~ARKODESolver()
{
   if (d_arkode_log_file) {
      fclose(d_arkode_log_file);
   }
   if (d_arkode_mem) {
      ARKStepFree(&d_arkode_mem);
   }
   if (d_linear_solver) {
      SUNLinSolFree(d_linear_solver);
   }
}

/*
 * Functions to initialize linear solver and reset ARKODE structure.
 */

void ARKODESolver::initializeARKODE()
{
   TBOX_ASSERT(d_solution_vector != 0);

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning(disable : 1572)
#endif

   if (d_ARKODE_needs_initialization) {

      /*
       * Set ARKODE log file.
       */
      if (d_arkode_log_file) {
         fclose(d_arkode_log_file);
      }
      d_arkode_log_file = fopen(d_arkode_log_file_name.c_str(), "w");

      /*
       * Make sure that either the relative tolerance or the
       * absolute tolerance has been set to a nonzero value.
       */
      bool tolerance_error = false;
      if (d_use_scalar_absolute_tolerance) {
         if ((d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_scalar == 0.0)) {
            tolerance_error = true;
         }
      } else {
         if ((d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_vector->maxNorm() == 0.0)) {
            tolerance_error = true;
         }
      }

      if (tolerance_error && d_arkode_log_file) {
         fprintf(d_arkode_log_file,
                 "%s: Both relative and absolute tolerance have value 0.0",
                 d_object_name.c_str());
      }

      /*
       * ARKODE function pointer.
       */
      ARKRhsFn RHSFunc = ARKODESolver::ARKODERHSFuncEval;

      /*
       * New function pointer for Implicit
       */
      ARKRhsFn RHSFuncImp = ARKODESolver::ARKODERHSFuncEvalImp;

      /*
       * New function pointer for Explicit
       */
      ARKRhsFn RHSFuncExp = ARKODESolver::ARKODERHSFuncEvalExp;

      /*
       * Free previously allocated ARKode memory.  Note that the
       * CVReInit() function is not used since the d_neq variable
       * might have been changed from the previous initializeARKODE()
       * call.
       */
      if (d_arkode_mem) ARKStepFree(&d_arkode_mem);

      /*
       * Allocate main memory for ARKODE package.
       */

      /*
       * Create ARKode member
       * im_ex = 0 is explicit
       * im_ex = 1 is implicit
       * im_ex = 2 IMEX
       */
      if (d_im_ex == 0) {
         d_arkode_mem =
             ARKStepCreate(RHSFunc, NULL, d_t_0, d_ic_vector->getNVector());
      } else {
         if (d_im_ex == 1) {
            d_arkode_mem =
                ARKStepCreate(NULL, RHSFunc, d_t_0, d_ic_vector->getNVector());
         } else {
            d_arkode_mem = ARKStepCreate(RHSFuncExp, RHSFuncImp, d_t_0,
                                         d_ic_vector->getNVector());
         }
      }

      //  ARKStepCreate takes place of CVodeInit

      int ierr = ARKStepSetUserData(d_arkode_mem, this);
      ARKODE_ERROR(ierr);

      ierr = ARKStepSStolerances(d_arkode_mem, d_relative_tolerance,
                                 d_absolute_tolerance_scalar);

      ierr = ARKStepSetOrder(d_arkode_mem, d_arkode_order);
      ARKODE_ERROR(ierr);

      if (d_im_ex != 0) {
         int pretype = d_uses_preconditioner ? PREC_LEFT : PREC_NONE;
         d_linear_solver = SUNLinSol_SPGMR(d_solution_vector->getNVector(),
                                           pretype, d_max_krylov_dim);

         ierr = ARKStepSetLinearSolver(d_arkode_mem, d_linear_solver, NULL);
         ARKODE_ERROR(ierr);
      }

      /*
       * Setup ARKSpgmr function pointers.
       */
      if (d_uses_preconditioner) {
         ARKLsPrecSetupFn precond_set = ARKODESolver::ARKSpgmrPrecondSet;
         ARKLsPrecSolveFn precond_solve = ARKODESolver::ARKSpgmrPrecondSolve;
         ierr =
             ARKStepSetPreconditioner(d_arkode_mem, precond_set, precond_solve);
         ARKODE_ERROR(ierr);
      }

      if (d_im_ex != 0 && !(d_nls_max_iter < 0)) {
         ierr = ARKStepSetMaxNonlinIters(d_arkode_mem, d_nls_max_iter);
         ARKODE_ERROR(ierr);
      }

      if (!(d_max_num_warnings < 0)) {
         ierr = ARKStepSetMaxNumConstrFails(d_arkode_mem, d_max_num_warnings);
         ARKODE_ERROR(ierr);
      }

      if (!(d_init_step_size < 0)) {
         ierr = ARKStepSetInitStep(d_arkode_mem, d_init_step_size);
         ARKODE_ERROR(ierr);
      }

      if (!(d_max_step_size < 0)) {
         ierr = ARKStepSetMaxStep(d_arkode_mem, d_max_step_size);
         ARKODE_ERROR(ierr);
      }

      if (!(d_min_step_size < 0)) {
         ierr = ARKStepSetMinStep(d_arkode_mem, d_min_step_size);
         ARKODE_ERROR(ierr);
      }
   }

   d_ARKODE_needs_initialization = false;
}

/*
 * Access methods for ARKODE statistics.
 */

void ARKODESolver::printARKODEStatistics(std::ostream& os) const
{

   char buf[256];

   os << "\nARKODESolver: ARKODE statistics... " << std::endl;

   sprintf(buf, "ns              = %5d     nsa              = %5d\n",
           getNumberOfSteps(), getNumberOfStepAttempts());
   os << buf;
   sprintf(buf, "nfi             = %5d     nfe              = %5d\n",
           getNumberOfRHSFunctionImpCalls(), getNumberOfRHSFunctionExCalls());
   os << buf;
   os << "nnli            = " << std::to_string(getNumberOfNonlinIters())
      << std::endl;
   os << "nnlf            = " << std::to_string(getNumberOfNonlinConvFails())
      << std::endl;
   os << "nsetups         = " << getNumberOfLinearSolverSetupCalls()
      << std::endl;
   os << "netf            = " << getNumberOfLocalErrorTestFailures()
      << std::endl;
   sprintf(buf, "hu              = %e      hcur             = %e\n",
           getStepSizeForLastInternalStep(), getStepSizeForNextInternalStep());
   os << buf;
   os << "tcur            = " << getCurrentInternalValueOfIndependentVariable()
      << std::endl;
}

/*
 * Access methods for ARKSpgmr statistics.
 */

void ARKODESolver::printARKSpgmrStatistics(std::ostream& os) const
{
   os << "ARKODESolver: ARKSpgmr statistics... " << std::endl;

   os << "nli             = " << std::to_string(getNumberOfLinearIterations())
      << std::endl;
   os << "nlf             = "
      << std::to_string(getNumberOfLinearConvergenceFailures()) << std::endl;
   os << "npe             = "
      << std::to_string(getNumberOfPreconditionerEvaluations()) << std::endl;
   os << "nps             = " << std::to_string(getNumberOfPrecondSolveCalls())
      << std::endl;
}

/*
 * Print ARKODESolver object data to given output stream.
 */
void ARKODESolver::printClassData(std::ostream& os) const
{
   os << "\nARKODESolver object data members..." << std::endl;
   os << "Object name = " << d_object_name << std::endl;

   os << "this = " << (ARKODESolver*)this << std::endl;
   os << "d_solution_vector = "
      << (solv::SundialsAbstractVector*)d_solution_vector << std::endl;

   os << "d_ARKODE_functions = " << (ARKODEAbstractFunctions*)d_arkode_functions
      << std::endl;

   os << "&d_arkode_mem = " << d_arkode_mem << std::endl;
   os << "d_arkode_log_file = " << (FILE*)d_arkode_log_file << std::endl;
   os << "d_arkode_log_file_name = " << d_arkode_log_file_name << std::endl;

   os << std::endl;
   os << "ARKODE parameters..." << std::endl;
   os << "d_t_0 = " << d_t_0 << std::endl;
   os << "d_ic_vector = " << (solv::SundialsAbstractVector*)d_ic_vector
      << std::endl;

   os << "d_linear_multistep_method = " << d_linear_multistep_method
      << std::endl;
   os << "d_relative_tolerance = " << d_relative_tolerance << std::endl;
   os << "d_use_scalar_absolute_tolerance = ";
   if (d_use_scalar_absolute_tolerance) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   }
   os << "d_absolute_tolerance_scalar = " << d_absolute_tolerance_scalar
      << std::endl;
   os << "d_absolute_tolerance_vector= " << std::endl;
   d_absolute_tolerance_vector->printVector();

   os << "Optional ARKODE inputs (see ARKODE docs for details):" << std::endl;

   os << "maximum number of non linear iterations = " << d_nls_max_iter
      << std::endl;
   os << "maximum number of nil internal step warnings = " << d_max_num_warnings
      << std::endl;

   os << "initial step size = " << d_init_step_size << std::endl;
   os << "maximum absolute value of step size = " << d_max_step_size
      << std::endl;
   os << "minimum absolute value of step size = " << d_min_step_size
      << std::endl;
   os << "last step size = " << getStepSizeForLastInternalStep() << std::endl;
   os << "...end of ARKODE parameters\n" << std::endl;

   os << std::endl;
   os << "ARKSpgmr parameters..." << std::endl;
   os << "d_precondition_type = " << d_precondition_type << std::endl;
   os << "d_gram_schmidt_type = " << d_gram_schmidt_type << std::endl;
   os << "d_max_krylov_dim = " << d_max_krylov_dim << std::endl;
   os << "d_tol_scale_factor = " << d_tol_scale_factor << std::endl;
   os << "...end of ARKSpgmr parameters\n" << std::endl;

   os << "d_ARKODE_needs_initialization = ";
   if (d_ARKODE_needs_initialization) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   }

   os << "...end of ARKODESolver object data members\n" << std::endl;
}
