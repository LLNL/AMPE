/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   FAC algorithm for solving linear equations on a hierarchy
 *
 ************************************************************************/

#ifndef included_solv_FACPreconditioner
#define included_solv_FACPreconditioner

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/solv/FACOperatorStrategy.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Utilities.h"

#include <algorithm>
#include <cctype>

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace solv {

/*!
 * @brief Implements the FAC iterative solution procedure
 * for a linear system of equations, Au=f, defined
 * on some subset of levels in an AMR patch hierarchy.
 *
 * The solution is found by applying an FAC algorithm
 * to the composite grid represented in the hierarchy.
 * After each FAC cycle the norm of the residual will be computed over
 * all levels involved.  The FAC iteration will stop when either the
 * maximum number of iterations is reached, or the residual norm on
 * all levels is below the given tolerance.
 *
 * The user must perform the following steps to use the FAC solver:
 * -# Create a FACPreconditioner, providing a valid concrete
 *    FACOperatorStrategy object.
 * -# Set the stopping criteria using the setStoppingParameters() function.
 * -# Set the number of smooting sweeps using the setSmoothingSweeps()
 *    function.  This is optional; if not used, the default is one sweep
 *    in each case.
 * -# Enable logging to record the FAC residual norms during the FAC
 *    itertion using the setFACLogFlag() function.  This is optional;
 *    the default is to turn logging off.  When loggin is turned on, the
 *    default mode is to send this information to the application log
 *    file (i.e., plog).
 * -# Invoke the FAC iteration process by calling solveSystem(),
 *    providing the vectors u and f, defined on a patch hierarchy
 *    containing the desired range of levels.
 * -# After solving, get solver statistics by viewing the log information
 *    and calling getNumberOfIterations(), getResidualNorm() functions
 *    if desired.
 */

class FACPreconditioner
{
public:
   /*!
    * Constructor.
    *
    * @param name Object name
    * @param user_ops Reference to user-specified FAC operator
    * @param database Input database with initialization parameters
    */
   FACPreconditioner(
      const std::string& name,
      FACOperatorStrategy& user_ops,
      const boost::shared_ptr<tbox::Database>& database =
         boost::shared_ptr<tbox::Database>());

   /*!
    * Virtual destructor.
    */
   virtual ~FACPreconditioner();

   /*!
    * @brief Solve linear system Au=f using the FAC algorithm.
    *
    * The return value is true if the solver
    * converged and false otherwise.
    * The problem-specific portions of the FAC procedure,
    * including the definitions of A are provided by the
    * FACOperatorStrategy object passed to the constructor.  More
    * information about the iteration can be found by calling the functions
    * getNumberOfIterations() and getResidualNorm() and by looking at the
    * log information.
    *
    * Before calling this function, the form of the solution and
    * right-hand-side quantities should be set properly by the user
    * on all patch interiors on the range of levels covered by the
    * FAC iteration.  All data in these vectors (that will be used
    * by the FACOperatorStrategy implementation) should be allocated.
    * Thus, the user is responsible for managing the
    * storage for the solution and right-hand-side.
    *
    * Conditions on arguments:
    * - vectors solution and rhs must have same hierarchy
    * - vectors solution and rhs must have same variables (except that
    *   solution can--and should--have enough ghost cells for computation).
    *
    * Upon return from this function,
    * the solution vector will contain the result of the solve.
    *
    * @param solution solution vector u
    * @param rhs right hand side vector f
    *
    * See initializeSolverState() and deallocateSolverState()
    * for opportunities to save overhead
    * when using multiple consecutive solves.
    *
    * @return whether solver converged to specified level
    *
    * @see initializeSolverState
    */
   bool
   solveSystem(
      SAMRAIVectorReal<double>& solution,
      SAMRAIVectorReal<double>& rhs);

   /*!
    * @brief Compute hierarchy-dependent data required for solving
    *
    * By default, the solveSystem() method
    * computes some required hierarchy-dependent data before
    * solving and removes that data after the solve.
    * For multiple solves using the same hierarchy configuration,
    * it is more efficient to manually compute, using
    * initializeSolverState(), and remove, using deallocateSolverState(),
    * the hierarchy-dependent data so that it is not done inside
    * solveSystem().  If solveSystem() detects that the solver state
    * is already initialized, it will @em NOT change the state.
    *
    * The vector arguments for solveSystem() need not match
    * those for initializeSolverState().  However, there must
    * be a certain degree of similarity, including
    * - hierarchy configuration (hierarchy pointer and level range)
    * - number, type and alignment of vector component data
    * - ghost cell width of data in the solution vector
    *
    * When assertion checking is enabled, limited checking is done
    * by solveSystem() to help ensure that of the vectors passed
    * to solveSystem() is compatible with the existing state.
    *
    * It is important to remember to reinitialize the solver state
    * when your hierarchy configuration changes.
    *
    * It is safe to initializeSolverState() when the state is
    * already initialized (the state is deallocated and reinitialized).
    *
    * Conditions on arguments:
    * - solution and rhs must have same hierarchy
    * - solution and rhs must have same structure, depth, etc.
    *   (except that u can--and should--have enough ghost cells
    *   for computation).
    * - coarsest_ln through finest_ln must exist in u.
    *
    * To unset the data set in this function,
    * see deallocateSolverState().
    *
    * After setting data for the current object, this function
    * calls the operator's corresponding function,
    * FACOperatorStrategy::initializeOperatorState()
    * so that the operator object can take steps to remain
    * in sync.
    *
    * @param solution solution vector u
    * @param rhs right hand side vector f
    */
   void
   initializeSolverState(
      const SAMRAIVectorReal<double>& solution,
      const SAMRAIVectorReal<double>& rhs);

   /*!
    * @brief Remove all hierarchy-dependent data computed by
    * initializeSolverState()
    *
    * Remove all hierarchy-dependent data set by initializeSolverState().
    * It is safe to call deallocateSolverState() even state is already
    * deallocated.
    *
    * After deallocating data for the current object, this function
    * calls the operator's corresponding function,
    * FACOperatorStrategy::deallocateOperatorState()
    * so that the operator object can take steps to remain
    * in sync.
    *
    * @see initializeSolverState
    */
   void
   deallocateSolverState();

   /*!
    * @brief Check compatibility of vectors with existing solver state.
    *
    * Check whether the solution and residual vectors given are
    * compatible with the existing solver state (solver state
    * must be initialized).  Compatibility means that the vectors
    * are sufficiently similar with the vectors with which the
    * state are initialized.  Compatibility implies that the
    * vectors may be used in solveSystem().
    *
    * The checking is not perfect!
    * Due to the possibility of user-defined patch data,
    * data-dependent checks cannot be performed.
    * It is possible that a false compatibility is returned.
    *
    * @return true if vectors are compatible with existing state
    */
   bool
   checkVectorStateCompatibility(
      const SAMRAIVectorReal<double>& solution,
      const SAMRAIVectorReal<double>& rhs) const;

   //@{
   //! @name Functions to set solving parameters.

   /*!
    * @brief Set the number of pre-smoothing sweeps during
    * FAC iteration process.
    *
    * Presmoothing is applied during the fine-to-coarse phase of the
    * iteration.  The default is to use one sweep.
    *
    * @param num_pre_sweeps Number of presmoothing sweeps
    */
   void
   setPresmoothingSweeps(
      int num_pre_sweeps)
   {
      d_presmoothing_sweeps = num_pre_sweeps;
   }

   /*!
    * @brief Set the number of post-smoothing sweeps during
    * FAC iteration process.
    *
    * Postsmoothing is applied during the coarse-to-fine phase of the
    * iteration.  The default is to use one sweep.
    *
    * @param num_post_sweeps Number of postsmoothing sweeps
    */
   void
   setPostsmoothingSweeps(
      int num_post_sweeps)
   {
      d_postsmoothing_sweeps = num_post_sweeps;
   }

   /*!
    * @brief Set the max number of iterations (cycles) to use per solve.
    */
   void
   setMaxCycles(
      int max_cycles)
   {
      d_max_iterations = max_cycles;
   }

   /*!
    * @brief Set the residual tolerance for stopping.
    *
    * The solution is considered converged if ||b-Ax|| <= residual_tol
    * @b or ||b-Ax|| <= relative_residual_tol*||b||.
    *
    * If you want the prescribed maximum number of cycles to always be taken,
    * set both residual tolerances to negative numbers.
    */
   void
   setResidualTolerance(
      double residual_tol,
      double relative_residual_tol = -1.0)
   {
      d_residual_tolerance = residual_tol;
      d_relative_residual_tolerance = relative_residual_tol;
   }

   /*!
    * @brief Set the choice of FAC cycling algorithm to use.
    *
    * For developer experimentation use only.
    * All others should use the default choice.
    *
    * @internal This function allows us to switch the cycling
    * algorithm to compare them.  This is mainly a debugging
    * feature and will be removed at some time.  Current
    * choices are:
    * - "default": the default recursive algorithm interpreted
    *        and coded by BTNG.
    * - "mccormick-s4.3": algorithm coded by BTNG, following Steve McCormick's
    *        section 4.3
    * - "pernice": algorithm coded by BTNG, interpretting the
    *        code originally written by Michael Pernice.
    */
   void
   setAlgorithmChoice(
      const std::string& choice);

   //@}

   //@{ @name Logging functions

   /*!
    * @brief Enable or disable logging.
    *
    * Set streams to NULL to turn off output.
    *
    * @param enabled Logging state.  true=on, false=off.
    */
   void
   enableLogging(
      bool enabled = true)
   {
      d_do_log = enabled;
   }

   //@}

   //@{
   //! @name Functions to get data on last solve.

   /*!
    * @brief Return FAC iteration count from last (or current
    * if there is one) FAC iteration process.
    */
   int
   getNumberOfIterations() const
   {
      return d_number_iterations;
   }

   /*!
    * @brief Get convergance rates of
    * the last (or current if there is one) FAC solve.
    *
    * The convergence factor is the factor to which the residual
    * has been reduced.  The final factor is that from the last
    * FAC cycle.
    *
    * @param avg_factor average convergence factor over FAC cycles
    *        from last solve.
    * @param final_factor convergence factor of the last FAC cycle
    */
   void
   getConvergenceFactors(
      double& avg_factor,
      double& final_factor) const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_number_iterations <= 0) {
         TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                                  << "a solve is invalid.\n");
      }
#endif
      avg_factor = d_avg_convergence_factor;
      final_factor = d_convergence_factor[d_number_iterations - 1];
   }

   /*!
    * @brief Get the net convergance rate of
    * the last (or current if there is one) FAC solve.
    *
    * The net factor is the factor to which the residual
    * has been reduced by the FAC cycles.
    * It is (current residual)/( initial residual + epsilon),
    * so it may not be accurate if the initial residual is very small.
    */
   double
   getNetConvergenceFactor() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_number_iterations <= 0) {
         TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                                  << "a solve is invalid.\n");
      }
#endif
      return d_net_convergence_factor;
   }

   /*!
    * @brief Get the average convergance rates of
    * the last (or current if there is one) FAC solve.
    *
    * The average factor is the net factor to the power of
    * 1/(number of FAC cycles).
    * It may not be accurate if the initial residual is very small.
    */
   double
   getAvgConvergenceFactor() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_number_iterations <= 0) {
         TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                                  << "a solve is invalid.\n");
      }
#endif
      return d_avg_convergence_factor;
   }

   /*!
    * @brief Get the final convergance rate of
    * the last (or current if there is one) FAC solve.
    *
    * The final factor is the factor to which the residual
    * has been reduced by the last FAC cycle.
    */
   double
   getFinalConvergenceFactor() const
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_number_iterations <= 0) {
         TBOX_ERROR(d_object_name << ": Seeking convergence factors before\n"
                                  << "a solve is invalid.\n");
      }
#endif
      return d_convergence_factor[d_number_iterations - 1];
   }

   /*!
    * @brief Return residual norm from the just-completed FAC iteration.
    *
    * The norm return value is computed as the maximum norm over all
    * patch levels involved in the solve.  The value corresponds to the
    * norm applied in the user-defined residual computation.
    *
    * The latest computed norm is the one returned.
    */
   double
   getResidualNorm() const
   {
      return d_residual_norm;
   }

   //@}

   /*!
    * @brief Print data members for debugging.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /*!
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   //@{
   //! @name Functions not implemented:
   FACPreconditioner(
      const FACPreconditioner&);
   void
   operator = (
      const FACPreconditioner&);
   //@}

   /*!
    * @brief Set state using database
    *
    * See the class description for the parameters that can be set
    * from a database.
    *
    * @param database Input database.  If a NULL pointer is given,
    * nothing is done.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& database);

   /*!
    * @brief Compute composite residual on all levels and
    * returns residual norm.
    *
    * Uses the FACOperatorStrategy::computeResidualOnLevel() function
    * provided by the operator object to compute per level residual.
    * Perform coarsen residual operation to get fine-grid approximation
    * of residual on to coarse grid.
    *
    * The residual is r=f-Au.
    *
    * @param residual residual vector r
    * @param solution solution vector u
    * @param rhs right hand side vector f
    */
   double
   computeFullCompositeResidual(
      SAMRAIVectorReal<double>& residual,
      SAMRAIVectorReal<double>& solution,
      SAMRAIVectorReal<double>& rhs);

   /*!
    * @brief Perform recursive FAC cycle iteration.
    *
    * Do one FAC iteration of Ae=r system.  The FAC algorithm
    * modifies Ae=r on coarser levels so that each coarser
    * level solves for the change in error of the next finer level,
    * as is expected in the FAC algorithm.
    *
    * The level number range lmax to lmin
    * must exist in the vectors e and r.
    *
    * Assumes:
    * - The error vector is preset to 0 on levels lmin to ln.
    *
    * @param error error vector e
    * @param residual residual vector r
    * @param solution solution vector u
    * @param lmax finest level number
    * @param lmin coarsest level number
    * @param ln current level number
    */
   void
   facCycle_Recursive(
      SAMRAIVectorReal<double>& error,
      SAMRAIVectorReal<double>& residual,
      SAMRAIVectorReal<double>& solution,
      int lmax,
      int lmin,
      int ln);

   /*!
    * @brief Perform recursive FAC cycle iteration from McCormick.
    *
    * Do one FAC iteration of Ae=r system.  The FAC algorithm
    * modifies Ae=r on coarser levels so that each coarser
    * level solves for the change in error of the next finer level,
    * as is expected in the FAC algorithm.
    *
    * The level number range lmax to lmin
    * must exist in the vectors e and r.
    *
    * Assumes:
    * - The error vector is preset to 0 on levels lmin to ln.
    *
    * @param error error vector e
    * @param residual residual vector r
    * @param solution solution vector u
    * @param lmax finest level number
    * @param lmin coarsest level number
    * @param ln current level number
    */
   void
   facCycle_McCormick(
      SAMRAIVectorReal<double>& error,
      SAMRAIVectorReal<double>& residual,
      SAMRAIVectorReal<double>& solution,
      int lmax,
      int lmin,
      int ln);

   /*!
    * @brief Perform FAC cycle iteration.
    *
    * Do one FAC iteration of Ae=r system.  The FAC algorithm
    * modifies Ae=r on coarser levels so that each coarser
    * level solves for the change in error of the next finer level,
    * as is expected in the FAC algorithm.
    *
    * The level number range lmax to lmin
    * must exist in the vectors e and r.
    *
    * Assumes:
    * - The error vector is preset to 0.
    *
    * @internal McCormick warned that cell-centered finite-volume
    * methods requires a W cycle, even in multigrid.  So this
    * function should be rewritten for that sort of flexibility.
    * Probably, a recursive function is needed.
    *
    * @param error error vector e
    * @param residual residual vector r
    * @param solution solution vector u
    * @param lmax finest level number
    * @param lmin coarsest level number
    */
   void
   facCycle(
      SAMRAIVectorReal<double>& error,
      SAMRAIVectorReal<double>& residual,
      SAMRAIVectorReal<double>& solution,
      int lmax,
      int lmin);

   /*!
    * @brief Name of this FAC solver object.
    */
   std::string d_object_name;

   /*!
    * @brief Object providing problem-specific routines.
    *
    * Reference is initialized by constructor @em never changes.
    */
   FACOperatorStrategy& d_fac_operator;

   //@{
   /*!
    * @name Solution vector-dependent data.
    *
    * These variables are set by
    * initializeSolverState and deallocateSolverState
    * and used only during the solve process.
    */

   boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
   int d_coarsest_ln;
   int d_finest_ln;

   /*!
    * @brief Clone of solution vector to store residual.
    */
   boost::shared_ptr<SAMRAIVectorReal<double> > d_residual_vector;

   /*!
    * @brief Clone of solution vector to store temporary residual.
    */
   boost::shared_ptr<SAMRAIVectorReal<double> > d_tmp_residual;

   /*!
    * @brief Error vector.
    */
   boost::shared_ptr<SAMRAIVectorReal<double> > d_error_vector;

   /*!
    * @brief Error vector for homogeneous boundary condition problem..
    */
   boost::shared_ptr<SAMRAIVectorReal<double> > d_tmp_error;

   //@}

   //@{
   /*!
    * @name Parameters for FAC iteration.
    */
   int d_max_iterations;
   double d_residual_tolerance;
   double d_relative_residual_tolerance;
   int d_presmoothing_sweeps;
   int d_postsmoothing_sweeps;
   std::string d_algorithm_choice;
   //@}

   //@{
   /*!
    * @name Status quantitities for FAC iteration.
    */
   int d_number_iterations;
   double d_residual_norm;

   /*!
    * @brief Norm of RHS, for computing relative residual.
    */
   double d_rhs_norm;
   /*!
    * @brief Convergence factor stack.
    *
    * The convergence factor stack is reset for each solve
    * and contains the convergence factors for each FAC cycle.
    */
   tbox::Array<double> d_convergence_factor;
   /*!
    * The average convergence factor computed from the current
    * values in d_convergence_factor.
    */
   double d_avg_convergence_factor;
   /*!
    * The net convergence factor computed from the current
    * values in d_convergence_factor.
    */
   double d_net_convergence_factor;
   //@}

   /*!
    * @brief Flag stating whether to log.
    */
   bool d_do_log;

   /*!
    * @brief Objects facilitating operations over a specific range
    * of levels.
    */
   tbox::Array<boost::shared_ptr<math::HierarchyDataOpsReal<double> > >
   d_controlled_level_ops;

   /*!
    * Timers for performance measurement.
    */
   boost::shared_ptr<tbox::Timer> t_solve_system;
};

}
}

#endif
