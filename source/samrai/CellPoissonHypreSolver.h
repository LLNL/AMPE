/*************************************************************************
 *
 * This file is adapted from the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE at https://github.com/LLNL/SAMRAI.
 *
 * Copyright:     (c) 1997-2021 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Poisson equation
 *
 ************************************************************************/
#ifndef included_CellPoissonHypreSolver
#define included_CellPoissonHypreSolver

#include "SAMRAI/SAMRAI_config.h"

extern "C" {
#include "HYPRE_struct_ls.h"
}

#include "PoissonSpecifications.h"

#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/solv/GhostCellRobinBcCoefs.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <iostream>

using namespace SAMRAI;

/*!
 * @brief Use the HYPRE preconditioner library to solve (the cell-centered)
 * Poisson's equation on a single level in a hierarchy.
 *
 * Class CellPoissonHypreSolver uses the HYPRE preconditioner library
 * to solve linear equations of the form
 * @f$ M \nabla ( D \nabla u ) + C u = f @f$, where
 * M is a cell-centered array,
 * C is a cell-centered array, D is a face-centered array,
 * and u and f are cell-centered arrays
 * (see PoissonSpecifications).
 * The discretization is the standard second order
 * finite difference stencil.
 *
 * Robin boundary conditions are used through the
 * interface class RobinBcCoefStrategy.
 * Periodic boundary conditions are supported for powers of 2 only.
 *
 * The user must perform the following steps to use
 * CellPoissonHypreSolver:
 * - Create a CellPoissonHypreSolver object.
 * - Initialize CellPoissonHypreSolver object with a patch hierarchy,
 *   using the function initializeSolverState().
 * - Use the functions setPhysicalBcCoefObject()
 *   to provide implementations of RobinBcCoefStrategy.
 *   (For most problems you can probably find a suitable
 *   implementation to use without implementing the
 *   strategy yourself.  See for example
 *   SimpleCellRobinBcCoefs and GhostCellRobinBcCoefs.)
 * - Set the matrix coefficients in the linear system,
 *   using the function setMatrixCoefficients().
 * - Specify the stopping criteria using setStoppingCriteria().
 * - Solve the linear system, passing in u and f as the patch
 *   indices of the solution and the right hand side, respectively.
 *
 * Sample parameters for initialization from database (and their
 * default values):
 * @verbatim
 *     print_solver_info = FALSE      // Whether to print some data for
 * debugging max_iterations = 10            // Max iterations used by Hypre
 *     num_pre_relax_steps = 1        // # of presmoothing steps used by Hypre
 *     num_post_relax_steps = 1       // # of postsmoothing steps used by Hypre
 *     use_smg = FALSE                // Whether to use hypre's smg solver
 *                                    // (alternative is the pfmg solver)
 * @endverbatim
 * The parameter relative_residual_tol (Residual tolerance used by Hypre,
 * default value=1.0e-8) should be set using the interface function
 * setStoppingCriteria(const double residual_tol)
 */

class CellPoissonHypreSolver
{
 public:
   /*!
    * @brief Constructor.
    *
    * @param object_name Name of object.
    * @param database tbox::InputDatabase for input.
    */
   CellPoissonHypreSolver(const std::string &object_name,
                          std::shared_ptr<tbox::Database> database =
                              std::shared_ptr<tbox::Database>());

   /*!
    * The Poisson destructor releases all internally managed data.
    */
   ~CellPoissonHypreSolver();

   /*!
    * @brief Initialize to a given hierarchy.
    *
    * Initializer Poisson solver for a patch level in a hierarchy.
    *
    * @param hierarchy Hierarchy
    * @param ln Level number
    */
   void initializeSolverState(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                              int ln = 0);

   /*!
    * @brief Reset to an uninitialized state.
    */
   void deallocateSolverState();

   /*!
    * @brief Set the matrix coefficients
    *
    * For information describing the Poisson equation parameters,
    * see the light-weight PoissonSpecifications class where
    * you set the values of C and D.
    *
    * This method must be called before solveSystem().
    */
   void setMatrixCoefficients(const PoissonSpecifications &spec);

   /*!
    * @brief Set default depth of the solution data involved in the solve.
    *
    * If the solution data has multiple depths,
    * the solver uses just one depth at a time.
    * The default depth is the first depth.
    * Use this function to change it.
    * This is not used to set the depth of the data (which is not
    * controled by this class) but the depth used in the solve.
    *
    * Changing the depth after setting up the matrix is permissible,
    * as the solution data does not affect the matrix.
    */
   void setSolnIdDepth(const int depth) { d_soln_depth = depth; }

   /*!
    * @brief Set default depth of the rhs data involved in the solve.
    *
    * If the rhs data has multiple depths,
    * the solver uses just one depth at a time.
    * The default depth is the first depth.
    * Use this function to change it.
    * This is not used to set the depth of the data (which is not
    * controled by this class) but the depth used in the solve.
    *
    * Changing the depth after setting up the matrix is permissible,
    * as the rhs data does not affect the matrix.
    */
   void setRhsIdDepth(const int depth) { d_rhs_depth = depth; }

   /*!
    * @brief Set the stopping criteria (residual
    * tolerance) for the linear solver.
    *
    * @param relative_residual_tol the maximum error tolerance
    */
   void setStoppingCriteria(const double relative_residual_tol)
   {
      TBOX_ASSERT(relative_residual_tol >= 0.0);
      d_relative_residual_tol = relative_residual_tol;
   }

   /*!
    * @brief Solve the linear system Au=f.
    *
    * The solution u and the right hand side f are
    * specified via patch indices on the patch hierarchy.
    *
    * Member functions getNumberOfIterations() return the iterations
    * from the solver.
    * Note that the matrix coefficients and boundary condition object
    * must have been set up before this function is called.
    * As long as the matrix coefficients do not change,
    * this routine may be called repeatedly to solve any number of linear
    * systems (with the right-hand side varying).
    * If the boundary conditions or matrix coefficients are changed
    * then function setMatrixCoefficients() must be called again.
    *
    * When computing the matrix coefficients in setMatrixCoefficients(),
    * the inhomogeneous portion of the boundary condition (constant
    * terms, independent of u and thus having no effect on the matrix)
    * are saved and added to the source term, f,
    * before performing the matrix solve.  In some situations, it may be
    * useful to not add the inhomogeneous portion to f.  The flag argument
    * @c homoegneous_bc is used for this.  (This is a sort of optimization,
    * to avoid having to re-call setMatrixCoefficients() to change the
    * inhomogeneous portion.)
    *
    * @param u Descriptor of cell-centered unknown variable.
    * @param f Descriptor of cell-centered source variable.
    * @param homogeneous_bc Whether homogeneous boundary conditions
    *        are assumed.
    *
    * @return whether solver converged to specified level
    */
   int solveSystem(const int u, const int f, bool homogeneous_bc = false,
                   bool initial_zero = false);

   /*!
    * @brief Return the number of iterations taken by the solver to converge.
    *
    * @return number of iterations taken by the solver to converge
    */
   int getNumberOfIterations() const { return d_number_iterations; }

   /*!
    * @brief Return the final residual norm returned by the Hypre solve.
    * @return final residual norm returned by the Hypre solve.
    */
   double getRelativeResidualNorm() const { return d_relative_residual_norm; }

   /*!
    * @brief Specify boundary condition directly, without using
    * a RobinBcCoefStrategy object.
    *
    * Use @em either setBoundaries() @em or setPhysicalBcCoefObject(),
    * but not both.
    *
    * A SimpleCelBcCoef object is used to interpret and implement
    * the specified boundary conditions.
    * See SimpleCellRobinBcCoefs::setBoundaries()
    * for an explanation of the arguments.
    */
   void setBoundaries(const std::string &boundary_type, const int fluxes = -1,
                      const int flags = -1, int *bdry_types = NULL)
   {
      d_physical_bc_simple_case.setBoundaries(boundary_type, fluxes, flags,
                                              bdry_types);
      d_physical_bc_coef_strategy = &d_physical_bc_simple_case;
      d_physical_bc_variable.reset();
   }

   /*!
    * @brief Specify boundary condition through the use of a
    * Robin boundary condition object.
    *
    * Use @em either setBoundaries() @em or setPhysicalBcCoefObject(),
    * but not both.
    *
    * The Robin boundary condition object is used when setting
    * the matrix coefficient and when solving the system.
    * If your boundary conditions are fixed values at ghost
    * cell centers, use the GhostCellRobinBcCoefs
    * implementation of the RobinBcCoefStrategy strategy.
    *
    * @param physical_bc_coef_strategy std::shared_ptr a concrete
    *        implementation of the Robin bc strategy.
    * @param variable hier::Variable pointer to be passed
    *        to RobinBcCoefStrategy::setBcCoefs(),
    *        but otherwise unused by this class.
    */
   void setPhysicalBcCoefObject(
       const solv::RobinBcCoefStrategy *physical_bc_coef_strategy,
       const std::shared_ptr<hier::Variable> variable =
           std::shared_ptr<hier::Variable>())
   {
      d_physical_bc_coef_strategy = physical_bc_coef_strategy;
      d_physical_bc_variable = variable;
   }


   void printConvergenceFactors(std::ostream &os);

   /*!
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string &getObjectName() const { return d_object_name; }

 private:
   /*!
    * @biief Set state using database
    *
    * See the class description for the parameters that can be set
    * from a database.
    *
    * @param database Input database.  If a NULL pointer is given,
    * nothing is done.
    */
   void getFromInput(std::shared_ptr<tbox::Database> database);

   void setupHypreSolver();
   void destroyHypreSolver();
   void allocateHypreData();
   void deallocateHypreData();

   void copyToHypre(HYPRE_StructVector vect, pdat::CellData<double> &src,
                    int depth, const hier::Box &box);
   void copyFromHypre(pdat::CellData<double> &dst, int depth,
                      HYPRE_StructVector vect, const hier::Box box);


   /*!
    * @brief Add g*A*k0(a) from boundaries to rhs.
    *
    * Move the constant portion of the boundary condition
    * contribution to the right hand side and add it to the existing rhs.
    * This operation is done for physical as well as cf boundaries,
    * so it is placed in a function.
    *
    * The boundary boxes given must be to either the physical
    * boundary or coarse-fine boundary for the patch.  The
    * bc coefficient implementation should correspond to the
    * boundary being worked on.
    */
   void add_gAk0_toRhs(const hier::Patch &patch,
                       const std::vector<hier::BoundaryBox> &bdry_boxes,
                       const solv::RobinBcCoefStrategy *robin_bc_coef,
                       pdat::CellData<double> &rhs);


   //@{

   /*!
    * @name Dimension-independent functions to organize Fortran interface.
    */

   //! @brief Compute diagonal entries of the matrix when C is variable.
   void computeDiagonalEntries(
       pdat::CellData<double> &diagonal, const pdat::CellData<double> &C_data,
       const pdat::SideData<double> &variable_off_diagonal,
       const hier::Box &patch_box);
   //! @brief Compute diagonal entries of the matrix when C is constant.
   void computeDiagonalEntries(
       pdat::CellData<double> &diagonal, const double C,
       const pdat::SideData<double> &variable_off_diagonal,
       const hier::Box &patch_box);
   //! @brief Compute diagonal entries of the matrix when C is zero.
   void computeDiagonalEntries(
       pdat::CellData<double> &diagonal,
       const pdat::SideData<double> &variable_off_diagonal,
       const hier::Box &patch_box);
   /*!
    * @brief Adjust boundary entries for variable off-diagonals.
    *
    * At the same time, save information that are needed to adjust
    * the rhs.
    */
   void adjustBoundaryEntries(pdat::CellData<double> &diagonal,
                              pdat::SideData<double> &variable_off_diagonal,
                              const hier::Box &patch_box,
                              const pdat::ArrayData<double> &acoef_data,
                              const pdat::ArrayData<double> &bcoef_data,
                              const hier::Box bccoef_box,
                              pdat::ArrayData<double> &Ak0_data,
                              const hier::BoundaryBox &trimmed_boundary_box,
                              const double h[SAMRAI::MAX_DIM_VAL]);

   //@}

   //! @brief Free static variables at shutdown time.
   static void freeVariables();


   /*!
    * @brief Object dimension.
    */
   const tbox::Dimension d_dim;
   unsigned short d_actual_dim;

   /*!
    * @brief Object name.
    */
   std::string d_object_name;

   /*!
    * @brief Associated hierarchy.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Associated level number.
    *
    * Currently, this must be level number 0.
    */
   int d_ln;

   /*!
    * @brief Scratch context for this object.
    */
   std::shared_ptr<hier::VariableContext> d_context;

   //@{ @name Boundary condition handling

   /*!
    * @brief The coarse-fine boundary description for level d_ln.
    *
    * The coarse-fine boundary is computed when the operator
    * state is initialized.  It is used to allow solves on
    * levels that are not the coarsest in the hierarchy.
    */
   std::shared_ptr<hier::CoarseFineBoundary> d_cf_boundary;

   /*!
    * @brief Robin boundary coefficient object for physical
    * boundaries.
    *
    * If d_physical_bc_coef_strategy is set, use it, otherwise,
    * use d_physical_bc_simple_case.
    */
   const solv::RobinBcCoefStrategy *d_physical_bc_coef_strategy;
   std::shared_ptr<hier::Variable> d_physical_bc_variable;

   /*!
    * @brief Implementation of Robin boundary conefficients
    * for the case of simple boundary conditions.
    */
   solv::SimpleCellRobinBcCoefs d_physical_bc_simple_case;

   /*!
    * @brief Robin boundary coefficient object for coarse-fine
    * boundaries.
    *
    * This is a GhostCellRobinBcCoefs object because we
    * expect the users to have the correct ghost cell values
    * in the coarse-fine boundaries before solving.
    */
   solv::GhostCellRobinBcCoefs d_cf_bc_coef;
   std::shared_ptr<hier::Variable> d_coarsefine_bc_variable;

   //@}

   /*!
    * @brief hier::Patch index of A*k0(a) quantity
    *
    * A*k0(a) is the quantity that is saved for
    * later adding to the rhs.
    *
    * The Robin bc is expressed by the coefficients a and g
    * on the boundary (see RobinBcCoefStrategy).
    * This class uses a central difference approximation of
    * the Robin bc, which results in the value at a ghost cell,
    * uo, being writen as uo = g*k0(a) + k1(a)*ui, where ui is
    * the first interior cell value, k0 and k1 depend on a as
    * indicated.
    *
    * In setting up the Au=f system, the contribution of k1(a)*ui
    * is incorporated into the product Au.  The contribution of
    * A*g*k0(a) should be moved to the right hand side and saved for
    * later adding to f.  However, the value of g is not provided
    * until solve time.  Therefore, we save just A*k0(a) at the
    * patch data index d_Ak0_id.
    */
   int d_Ak0_id;

   static std::shared_ptr<pdat::OutersideVariable<double> > s_Ak0_var[3];

   /*!
    * @brief Depth of the solution variable.
    */
   int d_soln_depth;

   /*!
    * @brief Depth of the rhs variable.
    */
   int d_rhs_depth;

   int d_max_iterations;
   double d_relative_residual_tol;


   int d_number_iterations;          // iterations in solver
   int d_num_pre_relax_steps;        // pre-relax steps in solver
   int d_num_post_relax_steps;       // post-relax steps in solver
   double d_relative_residual_norm;  // norm from solver

   /*@
    * @brief Flag to use SMG or PFMG (default)
    *
    * The flag is used to select which of HYPRE's linear solver algorithms
    * to use if true, the semicoarsening multigrid algorithm is used, and if
    * false, the "PF" multigrid algorithm is used.
    * By default, the SMG algorithm is used.
    *
    * Changing the algorithm must be done before setting up the matrix
    * coefficients.
    */
   bool d_use_smg;

   //@{
   //! @name Hypre object
   //! @brief HYPRE grid
   HYPRE_StructGrid d_grid;
   //! @brief HYPRE stencil
   HYPRE_StructStencil d_stencil;
   //! @brief HYPRE structured matrix
   HYPRE_StructMatrix d_matrix;
   //! @brief Hypre RHS vector for linear solves
   HYPRE_StructVector d_linear_rhs;
   //! @brief Hypre solution vector
   HYPRE_StructVector d_linear_sol;
   //! @brief Hypre SMG solver data
   HYPRE_StructSolver d_mg_data;
   //@}

   int d_msqrt_id;
   std::shared_ptr<pdat::CellVariable<double> > d_msqrt_var;
   bool d_msqrt_transform;

   //@{

   //! @name Variables for debugging and analysis.

   /*!
    * @brief Flag to print solver info
    *
    * See setPrintSolverInfo().
    */
   bool d_print_solver_info;

   bool d_converged;

   //@}

   /*!
    * @brief Timers for performance measurement.
    */
   std::shared_ptr<tbox::Timer> t_solve_system;
   std::shared_ptr<tbox::Timer> t_set_matrix_coefficients;
   std::shared_ptr<tbox::Timer> t_copy_vectors;
   std::shared_ptr<tbox::Timer> t_initsolver;
   std::shared_ptr<tbox::Timer> t_setupsolver;
};

#endif
