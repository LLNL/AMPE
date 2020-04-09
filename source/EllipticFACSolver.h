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
#ifndef included_EllipticFACSolver
#define included_EllipticFACSolver

#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/SAMRAI_config.h"

#include "EllipticFACOps.h"
#include "FACPreconditioner.h"

#include <boost/make_shared.hpp>
using namespace SAMRAI;

/*!
 * Note: This class is a generalization of the SAMRAI
 * solv::CellPoissonFACSolver. It solves an elliptic equation of the form
 *
 *     A div (D grad u) + C u = f
 *
 * whereas SAMRAI solv::CellPoissonFACSolver does not include the A factor.
 * Most of the code and the following documentation is the same, however.
 *
 * Class for solving scalar self-adjoint, elliptic equation on SAMR grid,
 * wrapping up lower-level components (FAC cycling, elliptic equation
 * operations and boundary conditions) in a single high-level interface.
 *
 * Note: this class provides a backward-compatible interface to
 * the soon-to-be obsolete PoissonHierarchySolver class.
 * Although this class hides the lower-level components (FAC cycling,
 * Poisson equation operations and boundary conditions), it is
 * perfectly acceptable to use those lower-level components directly.
 *
 * We solve the equation
 *    A div(D grad(u)) + Cu = f
 * where D is a side-centered array and C and A are cell-centered arrays,
 * u and f are also cell-centered.
 * Boundary conditions supported are Dirichlet, Neumann and mixed
 * (Dirichlet on some faces and Neumann on others).
 *
 * This class is a wrapper, providing a single class that coordinates
 * three major components: the FAC solver, the cell-centered Elliptic
 * FAC operator and a default Robin bc coefficient implelemtation.
 * It is perfectly acceptable to use those classes outside of this
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
 * -# Construct a EllipticFACSolver object, providing it
 *    the hierarchy and range of levels participating in the solve.
 * -# Set the parameters C, D and M using the functions named @c setC...
 *    and @c setD...  By default, D=1 and C=0 everywhere.
 * -# Call setBoundaries() to state the types boundary conditions,
 *    along with supplemental data for setting those boundary
 *    conditions.
 * -# Call initializeSolverState() to set up information
 *    internal to the solver.  This is step is not required
 *    but will save setup costs if you are making multiple
 *    solves.  This commits the object to the current hierarchy state
 *    and the specific @em types of boundary conditions you selected,
 *    It does NOT commit to the specific @em values of the boundary
 *    condition.  A hierarchy change (through adaption or other means)
 *    invalidates the state, thus you must reinitialize or
 *    deallocateSolverState() the state before another solve.
 * -# Solve the equation with solveSystem().  You provide the
 *    patch data indices for the solution u and the right hand
 *    side f.  u must have at least one ghost cell and where
 *    a Dirichlet boundary condition applies, those cells
 *    must be set to the value on the boundary.  If only Neumann
 *    boundary conditions are used, the ghost cell values
 *    do not matter.
 * -# Call deallocateSolverState() to free up internal resources,
 *    if initializeSolverState() was called before the solve.
 *
 * After the solve, information on the solve can be obtained
 * by calling one of these functions:
 * - getNumberOfIterations() gives the number of FAC cycles used.
 * - getConvergenceFactors() gives the average and final convergence
 *   factors for the solve.
 * - getResidualNorm() gives the final residual
 *
 * Finer solver controls can be set using the functions in this class.
 *
 * Object of this class can be set using input databases.
 * The following parameters can be set.  Each is shown with its
 * default value in the case where hypre is used.
 * @verbatim
 * enable_logging = TRUE // Bool flag to switch logging on/off
 * max_cycles = 10       // Integer number of max FAC cycles to use
 * residual_tol = 1.e-6  // Residual tolerance to solve for
 * num_pre_sweeps = 1    // Number of presmoothing sweeps to use
 * num_post_sweeps = 1   // Number of postsmoothing sweeps to use
 * coarse_fine_discretization = "Ewing" // Name of coarse-fine discretization
 * prolongation_method = "CONSTANT_REFINE" // Name of prolongation method
 * coarse_solver_choice = "hypre"  // Name of coarse level solver
 * coarse_solver_tolerance = 1e-10 // Coarse level tolerance
 * coarse_solver_max_iterations = 20 // Coarse level max iterations
 * use_smg = "FALSE"     // Whether to use hypre's smg solver
 *                       // (alternative is the pfmg solver)
 * @endverbatim
 *
 */
class EllipticFACSolver
{

 public:
   /*!
    * @brief Construct a solver.
    *
    * If the database is not NULL, initial settings will be set
    * using the database.
    * The solver is uninitialized until initializeSolverState()
    * is called.
    *
    * @param object_name Name of object used in outputs
    * @param database tbox::Database for initialization (may be NULL)
    */
   EllipticFACSolver(const std::string& object_name,
                     const boost::shared_ptr<EllipticFACOps> fac_ops,
                     const boost::shared_ptr<tbox::Database>& database =
                         boost::shared_ptr<tbox::Database>());

   /*!
    * @brief Destructor.
    */
   virtual ~EllipticFACSolver(void);

   /*
    * This function must be called after any of the C,D or A
    * coefficients are set or changed.  This lets the solver know
    * that no further changes will be made to these coefficients,
    * so that data that depends upon them can now be computed.
    */
   void finalizeCoefficients();

   /*!
    * @brief Solve Poisson's equation, assuming an uninitialized
    * solver state.
    *
    * Here, u is the "solution" patch data index and f is the
    * right hand side patch data index.
    * The return value is true if the solver converged and false otherwise.
    *
    * This function is a wrapper.
    * It simply initializes the solver state, call the
    * solveSystem(const int,const int) for the initialized solver then
    * deallocates the solver state.
    *
    * Upon return from this function,
    * solution will contain the result of the solve.
    *
    * See initializeSolverState() for opportunities to save overhead
    * when using multiple consecutive solves.
    *
    * @see solveSystem(const int,const int)
    *
    * @param solution hier::Patch data index for solution u
    * @param rhs hier::Patch data index for right hand side f
    * @param hierarchy The patch hierarchy to solve on
    * @param coarse_ln The coarsest level in the solve.
    * @param fine_ln The finest level in the solve.
    *
    * @return whether solver converged to specified level
    *
    * @see initializeSolverState
    */
   bool solveSystem(const int u_id, const int f_id, const int ew_id,
                    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                    int coarse_ln = -1, int fine_ln = -1);

   /*!
    * @brief Solve Poisson's equation using the current solver state
    * set by initializeSolverState().
    *
    * When the solver state has been initialized, this function may
    * be called repeadedly with different values on the rhs.
    * There is some cost savings for multiple solves when this
    * is done.
    *
    * Before calling this function, the solution and
    * right-hand-side quantities should be set properly by the user
    * on all patch interiors on the range of levels covered by the
    * FAC iteration.  All data for these patch data index should be allocated.
    * Thus, the user is responsible for managing the
    * storage for the solution and right-hand-side.
    *
    * @return whether solver converged to specified level
    *
    * @see solveSystem( const int, const int, boost::shared_ptr<
    * hier::PatchHierarchy >, int, int);
    */
   bool solveSystem(const int u_id, const int f_id, const int ew_id);

   /*!
    * @brief Specify the boundary conditions that are to be used at the
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

   /*!
    * @brief Override internal implementation to set boundary condition
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

   //@{ @name Functions for setting solver mathematic algorithm controls

   /*!
    * @brief Set max iterations for coarse level solve.
    *
    * If the coarse level solver requires a max iteration limit
    * (currently, they all do), the specified value is used.
    */
   void setCoarsestLevelSolverMaxIterations(int max_iterations);

#ifdef HAVE_HYPRE
   /*!
    * @brief Set whether to use HYPRe's PFMG algorithm instead of the
    * SMG algorithm.
    *
    * The flag is used to select which of HYPRE's linear solver algorithms
    * to use if true, the semicoarsening multigrid algorithm is used, and if
    * false, the ``PF'' multigrid algorithm is used.
    * By default, the SMG algorithm is used.
    *
    * This setting has effect only when HYPRe is chosen for the coarsest
    * level solver.  See setCoarsestLevelSolverChoice().
    *
    * Changing the algorithm must be done before setting up the matrix
    * coefficients.
    */
   void setUseSMG(bool use_smg);
#endif

   /*!
    * @brief Set the coarse-fine boundary discretization method.
    *
    * Specify the @c op_name string which will be passed to
    * xfer::Geometry::lookupRefineOperator() to get the operator
    * for setting fine grid ghost cells from the coarse grid.
    * Note that chosing this operator implicitly choses the
    * discretization method at the coarse-fine boundary.
    *
    * There is one important instance where this string is
    * @em not passed to xfer::Geometry::lookupRefineOperator().
    * If this variable is set to "Ewing", a constant refinement
    * method is used along with Ewing's correction.
    * For a reference to the correction method, see
    * "Local Refinement Techniques for Elliptic Problems on Cell-Centered
    * Grids, I. Error Analysis", Mathematics of Computation, Vol. 56, No. 194,
    * April 1991, pp. 437-461.
    *
    * @param coarsefine_method String selecting the coarse-fine discretization
    * method.
    */
   void setCoarseFineDiscretization(const std::string& coarsefine_method);

   void setVerbose(const bool verbose);

   /*!
    * @brief Set the number of pre-smoothing sweeps during
    * FAC iteration process.
    *
    * Presmoothing is applied during the fine-to-coarse phase of the
    * iteration.  The default is to use one sweep.
    *
    * @param num_pre_sweeps Number of presmoothing sweeps
    */
   void setPresmoothingSweeps(int num_pre_sweeps);

   /*!
    * @brief Set the number of post-smoothing sweeps during
    * FAC iteration process.
    *
    * Postsmoothing is applied during the coarse-to-fine phase of the
    * iteration.  The default is to use one sweep.
    *
    * @param num_post_sweeps Number of postsmoothing sweeps
    */
   void setPostsmoothingSweeps(int num_post_sweeps);

   /*!
    * @brief Set the residual tolerance for stopping.
    *
    * If you want the prescribed maximum number of cycles to always be taken,
    * set the residual tolerance to a negative number.
    *
    * Relative tolerance is used only if it is smaller than absolute tolerance
    */
   void setResidualTolerance(double residual_tol)
   {
      // tbox::pout << "  EllipticFACSolver::setResidualTolerance for precond to
      // "<<residual_tol<<endl;
      d_fac_precond.setResidualTolerance(residual_tol);
   }

   //@}

   /*!
    * @brief Prepare the solver's internal state for solving
    *
    * In the interest of efficiency, this class may prepare and
    * cache some hierarchy-dependent objects.  Though it is not required,
    * initializing the solver state makes for greater efficiency
    * when you are doing multiple solves on the same system of
    * equation.  If you do not initialize the state, it is initialized
    * and deallocated each time you call solveSystem(const int, const int).
    * The state must be reinitialized if the hierarchy or a boundary
    * condition type changes.
    *
    * To unset the data set in this function,
    * see deallocateSolverState().
    *
    * The @c solution and @c rhs patch data indices in the argument
    * list are used to determine the @em form of the data you
    * plan to use in the solve.  They need not be the same data
    * you solve on, but they should be similar.  Both must represent
    * cell-centered double data.  The solution must have at least one
    * ghost cell width, though this is not checked in the initialize
    * phase, because data is not required yet.
    *
    * @param solution solution patch data index for u
    * @param rhs right hand side patch data index for f
    * @param hierarchy The patch hierarchy to solve on
    * @param coarse_level The coarsest level in the solve
    * @param fine_level The finest level in the solve
    */
   void initializeSolverState(
       const int solution, const int rhs,
       const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int coarse_level = -1, const int fine_level = -1);

   /*!
    * @brief Remove the solver's internal state data
    *
    * Remove all hierarchy-dependent data set by initializeSolverState.
    * It is safe to call deallocateSolverState() even state is already
    * deallocated, but nothing is done in that case.
    *
    * @see initializeSolverState()
    */
   void deallocateSolverState();

   /*
    * -gamma*M*div(D*grad)+(I-gamma*M*C)
    */

   void resetSolverState(
       const int soln_id, const int rhs_id,
       const boost::shared_ptr<hier::PatchHierarchy> hierarchy);

   //@{
   //! @name Functions to get data on last solve.

   /*!
    * @brief Return FAC iteration count from last (or current
    * if there is one) FAC iteration process.
    */
   int getNumberOfIterations() const
   {
      return d_fac_precond.getNumberOfIterations();
   }


   /*!
    * @brief Get average convergance rate and convergence rate of
    * the last (or current if there is one) FAC solve.
    *
    * @param avg_factor average convergence factor over current FAC cycles
    * @param final_factor convergence factor of the last FAC cycle
    */
   void getConvergenceFactors(double& avg_factor, double& final_factor) const
   {
      d_fac_precond.getConvergenceFactors(avg_factor, final_factor);
      return;
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
   double getResidualNorm() const { return d_fac_precond.getResidualNorm(); }


   void printFACConvergenceFactors(const int solver_ret);

   //@}

 protected:
   /*!
    * @brief Set state using database
    *
    * See the class description for the parameters that can be set
    * from a database.
    *
    * @param database Input database.  If a NULL pointer is given,
    * nothing is done.
    */
   void getFromInput(const boost::shared_ptr<tbox::Database>& database);

   /*
    * @brief Set @c d_uv and @c d_fv to vectors wrapping the data
    * specified by patch data indices u and f.
    */
   void createVectorWrappers(int u, int f);

   /*
    * @brief Destroy vector wrappers referenced to by @c d_uv and @c d_fv.
    */
   void destroyVectorWrappers()
   {
      d_uv.reset();
      d_fv.reset();
   }

   /*!
    * @brief Object name.
    */
   std::string d_object_name;

   /*!
    * @brief FAC operator implementation corresponding to cell-centered
    * Poisson discretization.
    */
   boost::shared_ptr<EllipticFACOps> d_fac_ops;

   /*!
    * @brief FAC preconditioner algorithm.
    */
   FACPreconditioner d_fac_precond;

   /*!
    * @brief Robin bc object in use.
    */
   const solv::RobinBcCoefStrategy* d_bc_object;

   /*
    * @brief Default implementation of RobinBcCoefStrategy
    */
   solv::SimpleCellRobinBcCoefs d_simple_bc;

   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_ln_min;
   int d_ln_max;

   /*!
    * @brief Context for all internally maintained data.
    */
   boost::shared_ptr<hier::VariableContext> d_context;
   /*
    * @brief Vector wrapper for solution.
    * @see createVectorWrappers(), destroyVectorWrappers()
    */
   boost::shared_ptr<solv::SAMRAIVectorReal<double> > d_uv;
   /*
    * @brief Vector wrapper for source.
    * @see createVectorWrappers(), destroyVectorWrappers()
    */
   boost::shared_ptr<solv::SAMRAIVectorReal<double> > d_fv;

   bool d_solver_is_initialized;
   bool d_enable_logging;

   int d_vol_id;
   bool d_verbose;

   // boost::shared_ptr<math::HierarchyCellDataOpsReal<double> > d_hopscell;
};

#endif  // included_EllipticFACSolver
