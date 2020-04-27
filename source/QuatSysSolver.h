// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#ifndef included_QuatSysSolver_h
#define included_QuatSysSolver_h

/*
 * This class is a high-level wrapper for a quaternion system
 * solver.
 *
 * This file was adapted from solv::CellPoissonFACSolver.h
 * in the SAMRAI library.
 */

#include "SAMRAI/SAMRAI_config.h"
#include "QuatFACOps.h"

#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"

using namespace SAMRAI;

/*
 * This class implements a linear solver for the quaternion
 * equation described in the working notes "Differential
 * Algebraic Formulation of the Pusztai, Bortel and Granasy
 * Quaternion Phase Field Equation of Motion" by M. Dorr.
 * In particular, this class implements a Fast Adaptive
 * Composite (FAC) algorithm on a structured AMR grid of a
 * linear system whose matrix is given by equations (11).
 * This class wraps up lower-level components (FAC cycling,
 * linear operator operations and boundary conditions) in a
 * single high-level interface.
 *
 * This class is a wrapper, providing a single class that coordinates
 * three major components: the FAC solver, the QuatFAC
 * and a default Robin bc coefficient implelemtation.  It is
 * perfectly acceptable to use those classes outside of this
 * class.
 *
 * The underlying solver is an FAC solver using cell-centered
 * discretization.  The difference scheme is second-order
 * central-difference.  On coarse-fine boundaries within the
 * solution levels, the composite grid operator uses, by default,
 * the discretization method of Ewing, Lazarov and Vassilevski
 * ("Local Refinement Techniques for Elliptic Problems on
 * Cell-Centered Grids, I. Error Analysis", Mathematics of
 * Computation, Vol. 56, No. 194, April 1991, pp. 437-461).
 *
 * Typical use of this class is:
 * -# Construct a QuatSysSolver object, declaring
 *    an object name  (used maily for prepending various output
 *    messages) and an input database.
 * -# Call setBoundaries() to state the types boundary conditions,
 *    along with supplemental data for setting those boundary
 *    conditions.  Note: this call is not neeed if all of the
 *    physical boundaries are periodic.
 * -# Call initializeSolverState() to set up information
 *    internal to the solver.  The purpose of this member function
 *    is to precompute and store as much information as possible
 *    about the "state" of the linear solve, defined as those
 *    underlying structures that do not depend on specific
 *    numerical data used to set boundary values, matrix
 *    coefficients or the right-hand side.  In other words,
 *    this function commits the object to the current hierarchy
 *    state, as well as the specific data and boundary condition
 *    types, but does NOT commit it to the specific boundary values
 *    or operator coefficients.  A hierarchy change (through
 *    adaption or other means) invalidates the state, so you
 *    must reinitialize or deallocateSolverState() the state
 *    before another solve.
 * -# Call setOperatorCoefficients() to set the matrix coefficients.
 *    The coefficients depend upon five data objects: (1) a
 *    scalar factor related to the current time step, (2)
 *    diffusion coefficients and (3) quaternion component
 *    gradients at cell faces, (4) quaterion components
 *    at cell centers and (5) Lagrange multiplier at cell
 *    centers.
 * -# Solve the equation with FACSolve().  Ids are input
 *    for the variables containing the right-hand sides of the
 *    quaternion and constraint equations, as well as the
 *    variables that will return the quaternion and Lagrange
 *    multiplier solutions.  The quaterion variable must
 *    include at least one ghost cell, and those ghost cells
 *    must contain the values on the boundary wherever Dirichlet
 *    boundary conditions are imposed.
 * -# Call deallocateSolverState() to free up internal resources
 *
 * After the solve, information on the solve can be obtained
 * by calling one of these functions:
 * - getNumberOfIterations() gives the number of FAC cycles used.
 * - getConvergenceFactors() gives the average and final convergence
 *   factors for the solve.
 * - getResidualNorm() gives the final residual
 *
 * Finer solver controls can be set using accessor member functions,
 * as well as through the input database (if supplied).  The following
 * parameters can be set.  Each is shown with its default value
 *
 * max_cycles = 10                         // Maximum number of FAC cycles
 * allowed residual_tol = 1.e-6                    // Desired residual tolerance
 * coarse_fine_discretization = "Ewing"    // Name of coarse-fine discretization
 * prolongation_method = "CONSTANT_REFINE" // Name of prolongation method
 * levelsolver_tolerance = 1e-8            // Level solver tolerance
 * levelsolver_max_iterations = 20         // Level solver max iterations
 * coarse_levelsolver_tolerance = 1e-8     // Coarse level solver tolerance
 * coarse_levelsolver_max_iterations = 20  // Coarse level solver max iterations
 *
 */
class QuatSysSolver
{

 public:
   /*
    * Constructor
    *
    * If the database is not NULL, initial settings will be set
    * using the database.
    * The solver is uninitialized until initializeSolverState()
    * is called.
    *
    * object_name:   Name of object used in outputs
    * database:      tbox::Database for initialization (may be NULL)
    */
   QuatSysSolver(const int qlen, const std::string& object_name,
                 std::shared_ptr<tbox::Database> database =
                     std::shared_ptr<tbox::Database>());

   /*
    * Destructor
    */
   ~QuatSysSolver(void);

   /*
    * Prepare the solver's internal state for solving
    *
    * This member function precomputes and stores some hierarchy-dependent
    * objects, which provides greater efficiency when solving multiple
    * systems with the same structure (i.e., hierarchy configuration,
    * variable and boundary condition type), but different numerical data
    * (boundary values, operator coefficients and right-hand side vectors).
    * The state must be reinitialized if the hierarchy or a boundary
    * condition type changes.
    *
    * To unset the state data set in this function,
    * see deallocateSolverState().
    *
    * The variable ids q_soln_id, l_soln_id, q_rhs_id and l_rhs_id
    * in the argument list are used to determine the type of the data
    * to be used in the solve (i.e., centering, number of components,
    * number of ghost cells).  They need not be the same objects
    * used for the solve, but they must have the same structure.  All
    * four variables are cell-centered, and the variable identified
    * by q_soln_id must have at least one ghost cell, though this is
    * not checked in the initialize phase, because data is not
    * required yet.
    *
    * q_soln_id = solution array for the quaternion q
    * q_rhs_id  = right hand side array for the quaternion equations
    * l_soln_id = solution array for the Lagrange multiplier lambda
    * l_rhs_id  = right hand side array for the constraint equation
    * hierarchy = patch hierarchy to solve on
    */
   void initializeSolverState(
       const int q_soln_id, const int q_rhs_id, const int weight_id,
       std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void resetSolverState(
       const int q_soln_id, const int q_rhs_id, const int weight_id,
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   /*
    * Remove the solver's internal state data
    *
    * Remove all hierarchy-dependent data set by initializeSolverState.
    * It is safe to call deallocateSolverState() even if the state is already
    * deallocated, but nothing is done in that case.
    *
    * See initializeSolverState()
    */
   void deallocateSolverState();

   /*
    * Specify the boundary conditions that are to be used at the
    * physical domain boundary.
    *
    * This method is used to set up the default SimpleCellRobinBcCoefs
    * object for specifying boundary conditions.  Note that you may
    * alternatively provide your own implementation of the Robin
    * boundary condition coefficients using the setBcObject() method.
    *
    * The boundary conditions specified as the
    * string argument "boundary_type."  The boundary type argument can be
    * "Dirichlet", "Neumann", or "Mixed".
    *
    * If using Dirichlet boundary conditions, then before the solver is
    * called, the storage for the unknown u
    * must have a layer of ghost cells at least one cell wide that includes
    * the Dirichlet boundary values.
    *
    * If using Neumann boundary conditions, then before the solver is called,
    * the outerface boundary flux data must be set for the Neumann conditions.
    * The fluxes argument gives the patch data index of this flux
    * data.
    *
    * The mixed boundary type is for a mixture of Dirichlet and Neumann
    * boundary conditions are used at the physical domain boundary.
    * The fluxes argument gives the patch data index of the outerface data
    * that specifies the flux data for the Neumann conditions.  The flags
    * array is an outerface data array of integer flags that specifies whether
    * Dirichlet (flag == zero) or Neumann (flag == one) conditions are to be
    * used at a particular cell boundary face.  Note that the flag data must
    * be set before the matrix entries can be computed and the flux data
    * must be set before the solver is called.  The bdry_types argument can
    * be used if the boundary conditions are mixed but one or more of the
    * faces of the physical boundary are entirely either Dirichlet or
    * Neumann boundaries.  The bdry_types argument should be an array of
    * 2*DIM integers, specifying the boundary conditions on each side of
    * the physical domain.  It should be ordered {x_lo, x_hi, y_lo, y_hi,
    * z_lo, z_hi}, with the values for each face being 0 for Dirichlet
    * conditions, 1 for Neumann conditions, and 2 for mixed boundary
    * conditions.  The bdry_type argument is never required, but if used
    * it can sometimes make the PoissonHYPRESolver class more efficient.
    */

   void setBoundaries(const std::string& boundary_type, const int fluxes = -1,
                      const int flags = -1, int* bdry_types = NULL);

   /*
    * Override internal implementation to set boundary condition
    * coefficients with user-provided implementation.
    *
    * This function is used to override the default internal
    * object for setting Robin boundary condition coefficients.
    * You should override when you need to avoid the limitations
    * of the SimpleCellRobinBcCoefs class or you prefer to
    * use your own implementation.
    *
    * Note that an important limitation of the SimpleCellRobinBcCoefs
    * class is the inability to support linear interpolation in
    * the prolongation step.
    *
    * Once the boundary condition object is overwritten by this
    * method, you must no longer call the setBoundaries() method.
    *
    * bc_object: Pointer to boundary condition object
    */
   void setBcObject(const solv::RobinBcCoefStrategy* bc_object)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!bc_object) {
         TBOX_ERROR(d_object_name << ": NULL pointer for boundary condition\n"
                                  << "object.\n");
      }
#endif
      d_bc_object = bc_object;
      d_fac_ops->setPhysicalBcCoefObject(d_bc_object);
   }

   void setOperatorCoefficients(const double time_step, const double epsilon_q,
                                const double quat_grad_floor,
                                const std::string quat_smooth_floor_type,
                                const int mobility_id,
                                const int mobility_deriv_id,
                                const int diff_coef_id,
                                const int diff_coef_deriv_id,
                                const int grad_q_id, const int q_id);

   bool solveSystem(const int q_soln_id, const int q_rhs_id);

   void evaluateRHS(const double epsilon_q, const double quat_grad_floor,
                    const std::string quat_smooth_floor_type,
                    const int diff_coef_id, const int grad_q_id,
                    const int grad_q_copy_id, const int rotations_id,
                    const int mobility_id, const int solution_id, int rhs_id,
                    const bool use_gradq_for_flux = false);

   void multiplyDQuatDPhiBlock(const int q_id, const int operator_q_id);

   void applyProjection(const int q_id, const int corr_id, const int err_id);

   /*
    * Enable printing of solver convergence information to the
    * tbox::pout stream.
    *
    * verbose: Print to tbox::pout?
    */
   void setVerbose(bool verbose)
   {
      d_verbose = verbose;
      d_fac_ops->setVerbose(verbose);
   }

   /*
    * Set the relative residual tolerance to which each
    * level solve will be performed.  For the coarsest level,
    * this tolerance is overridden by any value set by
    * setCoarsestLevelSolverTolerance().
    *
    * tol: Desired relative residual tolerance
    */
   void setLevelSolverTolerance(double tol)
   {
      d_fac_ops->setLevelSolverTolerance(tol);
   }

   /*
    * Set the maximum number of iterations allowed for
    * level solves.  For the coarsest level, this limit
    * is overridden by any value set by
    * setCoarsestLevelSolverMaxIterations().
    *
    * max_iterations: Maximum allowed number of level solve iterations
    */
   void setLevelSolverMaxIterations(int max_iterations)
   {
      d_fac_ops->setLevelSolverMaxIterations(max_iterations);
   }

   /*
    * Set the relative residual tolerance to which the
    * coarsest level solve will be performed.  This value
    * overrides values set by setLevelSolverTolerance().
    *
    * tol: Desired relative residual tolerance on the coarsest level
    */
   void setCoarsestLevelSolverTolerance(double tol)
   {
      d_fac_ops->setCoarsestLevelSolverTolerance(tol);
   }

   /*
    * Set the maximum number of iterations allowed for
    * coarsest level solve.  For the coarsest level, this
    * limit overrides values set by setLevelSolverMaxIterations().
    *
    * max_iterations: Maximum allowed number of level solve iterations
    *                 on the coarsest level
    */
   void setCoarsestLevelSolverMaxIterations(int max_iterations)
   {
      d_fac_ops->setCoarsestLevelSolverMaxIterations(max_iterations);
   }

   /*
    * Set the coarse-fine boundary discretization method.
    *
    * Specify the op_name string which will be passed to
    * xfer::Geometry::lookupRefineOperator() to get the operator
    * for setting fine grid ghost cells from the coarse grid.
    * Note that chosing this operator implicitly choses the
    * discretization method at the coarse-fine boundary.
    *
    * There is one important instance where this string is
    * not passed to xfer::Geometry::lookupRefineOperator().
    * If this variable is set to "Ewing", a constant refinement
    * method is used along with Ewing's correction.
    * For a reference to the correction method, see
    * "Local Refinement Techniques for Elliptic Problems on Cell-Centered
    * Grids, I. Error Analysis", Mathematics of Computation, Vol. 56, No. 194,
    * April 1991, pp. 437-461.
    *
    * coarsefine_method: String selecting the coarse-fine discretization method.
    */
   void setCoarseFineDiscretization(const std::string& coarsefine_method)
   {
      d_fac_ops->setCoarseFineDiscretization(coarsefine_method);
   }

   /*
    * Set the name of the prolongation method.
    *
    * Specify the op_name string which will be passed to
    * xfer::Geometry::lookupRefineOperator() to get the operator
    * for prolonging the coarse-grid correction.
    *
    * By default, "CONSTANT_REFINE" is used.  "LINEAR_REFINE" seems t
    * to lead to faster convergence, but it does NOT satisfy the Galerkin
    * condition.
    *
    * Prolonging using linear refinement requires a Robin bc
    * coefficient implementation that is capable of delivering
    * coefficients for non-hierarchy data, because linear refinement
    * requires boundary conditions to be set on temporary levels.
    *
    * prolongation_method: String selecting the coarse-fine discretization
    * method.
    */
   void setProlongationMethod(const std::string& prolongation_method)
   {
      d_fac_ops->setProlongationMethod(prolongation_method);
   }

   /*
    * Set the maximum number of FAC iterations (cycles) to use per solve.
    *
    * max_cycles: Maximum allowed number of FAC cycles
    */
   void setMaxCycles(int max_cycles);

   /*
    * Set the desired relative residual tolerance.
    *
    * If you want the prescribed maximum number of cycles to always be taken,
    * set the residual tolerance to a negative number.
    *
    * residual_tol: Residual tolerance
    */
   void setResidualTolerance(double residual_tol)
   {
      d_fac_solver.setResidualTolerance(residual_tol);
   }

   /*
    * Return the FAC iteration count from last (or current if there is one)
    * FAC iteration process.
    */
   int getNumberOfIterations() const
   {
      return d_fac_solver.getNumberOfIterations();
   }

   /*
    * Get the average convergance rate and convergence rate of
    * the last (or current if there is one) FAC solve.
    *
    * avg_factor:   Average convergence factor over current FAC cycles
    * final_factor: Convergence factor of the last FAC cycle
    */
   void getFACConvergenceFactors(double& avg_factor, double& final_factor) const
   {
      d_fac_solver.getConvergenceFactors(avg_factor, final_factor);
   }

   /*
    * Return the relative residual norm from the just-completed FAC iteration.
    *
    * The norm return value is computed as the maximum norm over all
    * patch levels involved in the solve.  The value corresponds to the
    * norm applied in the user-defined residual computation.
    *
    * The latest computed norm is the one returned.
    */
   double getResidualNorm() const { return d_fac_solver.getResidualNorm(); }

   /*
    * Print solver data
    *
    * solver_ret: solver return code
    */
   void printFACConvergenceFactors(const int solver_ret);

   int getFaceDiffCoeffId() { return d_fac_ops->getFaceDiffCoeffId(); }
   int getFaceDiffCoeffScratchId()
   {
      return d_fac_ops->getFaceDiffCoeffScratchId();
   }

 private:
   /*
    * Read parameters from the database
    *
    * See the above class description for the parameters that
    * can be set from a database.
    *
    * database: Input database.  If a NULL pointer is given,
    * nothing is done.
    */
   void getFromInput(const std::shared_ptr<tbox::Database>& database);

   /*
    * Set d_uv and d_fv to vectors wrapping the data
    * specified by patch data indices.
    */
   void createVectorWrappers(
       int q_u, int q_f, std::shared_ptr<solv::SAMRAIVectorReal<double> >& uv,
       std::shared_ptr<solv::SAMRAIVectorReal<double> >& fv);

   /*
    * Destroy the vector wrappers referenced to by d_uv and d_fv.
    */
   void destroyVectorWrappers(
       std::shared_ptr<solv::SAMRAIVectorReal<double> >& uv,
       std::shared_ptr<solv::SAMRAIVectorReal<double> >& fv);

   /*
    * Object name.
    */
   std::string d_object_name;

   /*
    * Context for all internally maintained data.
    */
   std::shared_ptr<hier::VariableContext> d_context;

   /*
    * FAC operator implementation
    */
   std::shared_ptr<QuatFACOps> d_fac_ops;

   /*
    * The FAC solver used for the preconditioner
    */
   FACPreconditioner d_fac_solver;

   /*
    * Robin bc object in use.
    */
   const solv::RobinBcCoefStrategy* d_bc_object;

   /*
    * Default implementation of RobinBcCoefStrategy
    */
   solv::SimpleCellRobinBcCoefs d_simple_bc;

   /*
    * Hierarchy and min and max levels
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_ln_min;
   int d_ln_max;

   /*
    * Descriptor index for an internally allocated variable
    * containing volume weights used in composite grid norm
    * calculations.
    */
   int d_weight_id;

   /*
    * Vector wrapper temporary for the solution.
    * See createVectorWrappers(), destroyVectorWrappers()
    */
   std::shared_ptr<solv::SAMRAIVectorReal<double> > d_uv;

   /*
    * Vector wrapper temporary for the right-hand side.
    * See createVectorWrappers(), destroyVectorWrappers()
    */
   std::shared_ptr<solv::SAMRAIVectorReal<double> > d_fv;

   /*
    * Solver initialization state
    */
   bool d_solver_is_initialized;

   /*
    * Boolean flag controlling the output of iteration data
    * to the tbox::pout stream.
    */
   bool d_verbose;

   /*
    * Timers for performance measurement.
    */
   std::shared_ptr<tbox::Timer> t_set_op_coef;
   std::shared_ptr<tbox::Timer> t_solve_system;

   /*
    * Relative residual tolerance for the outer GMRES iteration
    */
   double d_residual_tol;

   /*
    * Relative residual tolerance for the FAC preconditioner iteration
    */
   double d_precond_tol;

   /*
    * Maximum number of iterations for the FAC preconditioner iteration
    */
   int d_precond_maxiters;
};

#endif  // included_QuatSysSolver_h
