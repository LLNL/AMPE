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
#ifndef included_QuatLevelSolver_h
#define included_QuatLevelSolver_h

/*
 * This class implements a level solver for a quaternion system
 * FAC solver.
 *
 * This file was adapted from solv::CellPoissonHypreSolver.C in
 * the SAMRAI library.
 */

#ifndef included_SAMRAI_config
#include "SAMRAI/SAMRAI_config.h"
#endif

#ifndef included_HYPRE_sstruct_ls
#define included_HYPRE_sstruct_ls
extern "C" {
#include "HYPRE_sstruct_ls.h"
}
#endif

#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/solv/GhostCellRobinBcCoefs.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/tbox/Database.h"

using namespace SAMRAI;

/*
 * Class QuatLevelSolver uses the HYPRE sysPFMG solver
 * to solve the quaternion system described in the working notes
 * "Differential Algebraic Formulation of the Pusztai, Bortel and
 * Granasy Quaternion Phase Field Equation of Motion" by M. Dorr
 * on a single level in a hierarchy.  The template parameter DIM
 * is the dimensionality of the underlying spatial domain.
 *
 * The solution and right-hand side are cell-centered, and the
 * discretization utilizes standard second order finite difference
 * stencils.
 *
 * Robin boundary conditions are used through the
 * interface class RobinBcCoefStrategy.
 *
 * The user must perform the following steps to use
 * QuatLevelSolver:
 * - Create a QuatLevelSolver object.
 * - Initialize QuatLevelSolver object with a patch hierarchy,
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
 *     print_solver_info = FALSE      // Whether to print some data for
 * debugging solver_type = 1                // Hypre solver type max_iterations
 * = 10            // Maximum number of iterations allowed relative_residual_tol
 * = 1.0e-8 // Relative residual tolerance desired num_pre_relax_steps = 1 // #
 * of presmoothing steps used by PFMG num_post_relax_steps = 1       // # of
 * postsmoothing steps used by PFMG
 */

class QuatLevelSolver
{
 public:
   /*
    * Constructor
    *
    * object_name: Name of object.
    * database:    tbox::Database for input.
    */
   QuatLevelSolver(const int qlen, const std::string& object_name,
                   std::shared_ptr<tbox::Database> database =
                       std::shared_ptr<tbox::Database>());

   /*
    * The destructor releases all internally managed data.
    */
   ~QuatLevelSolver();

   /*
    * Initialize with a given hierarchy and patch level with
    * that hierarchy
    *
    * hierarchy: Hierarchy
    * ln:        Level number
    */
   void initializeSolverState(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                              int ln = 0);

   /*
    * Reset to an uninitialized state
    */
   void deallocateSolverState();

   /*
    * Set the matrix coefficients
    *
    * This method must be called before solveSystem().
    */
   void setMatrixCoefficients(const double time_step,
                              const int sqrt_mobility_id,
                              const int face_coef_id);

   /*
    * Set the stopping criteria (max iterations and residual
    * tolerance) for the linear solver.
    *
    * max_iterations = the maximum number of iterations
    * residual_tol = the maximum error tolerance
    */
   void setStoppingCriteria(const int max_iterations = 10,
                            const double relative_residual_tol = 1.0e-6);

   /*
    * Solve the linear system Au=f.
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
    * homoegneous_bc is used for this.  (This is a sort of optimization,
    * to avoid having to re-call setMatrixCoefficients() to change the
    * inhomogeneous portion.)
    *
    * q_rhs_id:       Descriptor of rhs for quaternion equations
    * l_rhs_id:       Descriptor of rhs for constraint equation
    * q_solution_id:  Descriptor of quaternion solution
    * l_solution_id:  Descriptor of Lagrange multiplier solution
    * homogeneous_bc: Homogeneous boundary conditions assumed?
    *
    * Returns whether solver converged to specified level
    */
   int solveSystem(const int q_rhs_id, const int q_solution_id,
                   bool homogeneous_bc = false);

   /*
    * Return the number of iterations taken by the solver to converge.
    */
   int getNumberOfIterations() const;

   /*
    * Set the number of pre-relax steps used by the Hypre solve.
    */
   void setNumPreRelaxSteps(const int steps);

   /*
    * Set the number of post-relax steps used by the Hypre solve.
    */
   void setNumPostRelaxSteps(const int steps);

   /*
    * Return the final residual norm returned by the Hypre solve.
    */
   double getRelativeResidualNorm() const;

   /*
    * Specify boundary condition directly, without using
    * a RobinBcCoefStrategy object.
    *
    * Use either setBoundaries() @em or setPhysicalBcCoefObject(),
    * but not both.
    *
    * A SimpleCelBcCoef object is used to interpret and implement
    * the specified boundary conditions.
    * See SimpleCellRobinBcCoefs::setBoundaries()
    * for an explanation of the arguments.
    */
   void setBoundaries(const std::string& boundary_type, const int fluxes = -1,
                      const int flags = -1, int* bdry_types = NULL);

   /*
    * Specify boundary condition through the use of a
    * Robin boundary condition object.
    *
    * Use either setBoundaries() or setPhysicalBcCoefObject(),
    * but not both.
    *
    * The Robin boundary condition object is used when setting
    * the matrix coefficient and when solving the system.
    * If your boundary conditions are fixed values at ghost
    * cell centers, use the GhostCellRobinBcCoefs
    * implementation of the RobinBcCoefStrategy strategy.
    *
    * physical_bc_coef_strategy ==std::shared_ptr to a concrete
    *        implementation of the Robin bc strategy.
    * variable == hier::Variable pointer to be passed
    *        to RobinBcCoefStrategy::setBcCoefs(),
    *        but otherwise unused by this class.
    */
   void setPhysicalBcCoefObject(
       const solv::RobinBcCoefStrategy* physical_bc_coef_strategy,
       const std::shared_ptr<hier::Variable> variable =
           std::shared_ptr<hier::Variable>());

   /*
    * Set the flag for printing solver information.
    *
    * This optional function is used primarily for debugging.
    *
    * If set true, it will print the HYPRE matrix information
    * to the following files:
    *
    * - mat_bA.out - before setting matrix coefficients in matrix assemble
    * - mat_aA.out - after setting matrix coefficients in matrix assemble
    * - sol0.out   - u before solve (i.e. for system Au = b)
    * - sol.out    - u after solve
    * - mat0.out   - A before solve
    * - mat.out    - A after solve
    * - rhs.out    - b before and after solve
    *
    * If this method is not called, or the flag is set false, no printing
    * will occur.
    */
   void setPrintSolverInfo(const bool print);

   /*
    * Set the flag for printing iteration information to the
    * tbox::pout stream.
    *
    * verbose: Print?
    */
   void setVerbose(const bool verbose);

 private:
   /*
    * Set state using database
    *
    * See the class description for the parameters that can be set
    * from a database.
    *
    * database: Input database.  If a NULL pointer is given,
    *           nothing is done.
    */
   void getFromInput(std::shared_ptr<tbox::Database> database);

   /*
    * Allocate the data structures need by Hypre, except for
    * the actual solver itself, which is allocated by setupHypreSolver().
    */
   void allocateHypreData();

   /*
    * Deallocate the Hypre data structures allocated by allocateHypreData().
    */

   void deallocateHypreData();

   /*
    * Allocate and initialize the Hypre solver
    */
   void createHypreSolver(int var);

   /*
    * Deallocate the Hypre solver allocated by createHypreSolver()
    */
   void destroyHypreSolver(int var);

   /*
    * Copy data from the "component" component of the "src_data" CellData to the
    * "var" component of the Hypre vector "dst_vector"
    */

   void copyToHypre(const int component, const pdat::CellData<double>& src_data,
                    HYPRE_SStructVector dst_vector);

   /*
    * Copy data from the "var" component of the Hypre vector "src_vector" to the
    * "component" component of the "dst_data" CellData
    */

   void copyFromHypre(const int var, const HYPRE_SStructVector src_vector,
                      const int component, pdat::CellData<double>& dst_data);

   /*
    * Add g*A*k0(a) from boundaries to rhs.
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
   void add_gAk0_toRhs(const hier::Patch& patch,
                       const std::vector<hier::BoundaryBox>& bdry_boxes,
                       const solv::RobinBcCoefStrategy* robin_bc_coef,
                       pdat::CellData<double>& rhs);

   /*
    * Dimension-independent functions to organize Fortran interface.
    */

   /*
    * Adjust boundary entries for variable off-diagonals.
    *
    * At the same time, save information that are needed to adjust
    * the rhs.
    */
   void adjustBoundaryEntries(
       const double time_step, pdat::CellData<double>& sqrt_mobility_data,
       pdat::CellData<double>& diagonal,
       const pdat::SideData<double>& variable_off_diagonal,
       const hier::Box& patch_box, const pdat::ArrayData<double>& acoef_data,
       const pdat::ArrayData<double>& bcoef_data, const hier::Box bccoef_box,
       pdat::OutersideData<double>& Ak0,
       const hier::BoundaryBox& trimmed_boundary_box, const double h[NDIM]);

   // Free static variables at shutdown time.
   static void freeVariables();

   /*!
    * @brief Object dimension.
    */
   const tbox::Dimension d_dim;

   int d_qlen;

   /*
    * Object name.
    */
   std::string d_object_name;

   /*
    * Associated hierarchy.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*
    * Level number.
    */
   int d_ln;

   /*
    * Scratch context for this object.
    */
   std::shared_ptr<hier::VariableContext> d_context;

   /*
    * The coarse-fine boundary description for level d_ln.
    *
    * The coarse-fine boundary is computed when the operator
    * state is initialized.  It is used to allow solves on
    * levels that are not the coarsest in the hierarchy.
    */
   std::shared_ptr<hier::CoarseFineBoundary> d_cf_boundary;

   /*
    * Robin boundary coefficient object for physical
    * boundaries.
    *
    * If d_physical_bc_coef_strategy is set, use it, otherwise,
    * use d_physical_bc_simple_case.
    */
   const solv::RobinBcCoefStrategy* d_physical_bc_coef_strategy;
   std::shared_ptr<hier::Variable> d_physical_bc_variable;

   /*
    * Implementation of Robin boundary conefficients
    * for the case of simple boundary conditions.
    */
   solv::SimpleCellRobinBcCoefs d_physical_bc_simple_case;

   /*
    * Robin boundary coefficient object for coarse-fine
    * boundaries.
    *
    * This is a GhostCellRobinBcCoefs object because we
    * expect the users to have the correct ghost cell values
    * in the coarse-fine boundaries before solving.
    */
   solv::GhostCellRobinBcCoefs d_cf_bc_coef;
   std::shared_ptr<hier::Variable> d_coarsefine_bc_variable;

   /*
    * hier::Patch index of A*k0(a) quantity
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

   static std::shared_ptr<pdat::OutersideVariable<double> > s_Ak0_var[NDIM];

   /*
    * Maximum number of iterations allowed
    */
   int d_max_iterations;

   /*
    * Desired relative residual tolerance
    */
   double d_relative_residual_tol;

   /*
    * Number of pre-relax steps to be performed by Hypre solver
    */
   int d_num_pre_relax_steps;

   /*
    * Number of post-relax steps to be performed by Hypre solver
    */
   int d_num_post_relax_steps;

   /*
    * Number of iterations performed by solver (last solve)
    */
   int d_number_iterations;

   /*
    * Relative residual achieved by solver (last solve)
    */
   double d_relative_residual_norm;

   /*
    * Flag to indicate whether iteration data is to be
    * output to tbox::pout
    */
   bool d_verbose;

   int d_solver_type;

   // Hypre data objects
   HYPRE_SStructGrid d_grid;
   HYPRE_SStructGraph d_graph;
   HYPRE_SStructStencil d_stencil;

   HYPRE_SStructMatrix* d_matrix;
   HYPRE_SStructVector* d_linear_rhs;
   HYPRE_SStructVector* d_linear_sol;

   HYPRE_Solver* d_parcsr_solver;
   HYPRE_Solver* d_parcsr_precond;
   HYPRE_ParCSRMatrix* d_parcsr_A;
   HYPRE_ParVector* d_parcsr_b;
   HYPRE_ParVector* d_parcsr_x;

   HYPRE_StructSolver* d_struct_solver;
   HYPRE_StructSolver* d_struct_precond;
   HYPRE_StructMatrix* d_struct_A;
   HYPRE_StructVector* d_struct_b;
   HYPRE_StructVector* d_struct_x;

   int d_hypre_object_type;

   bool d_hypre_data_allocated;

   // Variables for debugging and analysis.

   /*
    * Flag to print solver info
    *
    * See setPrintSolverInfo().
    */
   bool d_print_solver_info;

   /*
    * Timers for performance measurement.
    */
   std::shared_ptr<tbox::Timer> t_solve_system;
   std::shared_ptr<tbox::Timer> t_set_matrix_coefficients;
   std::shared_ptr<tbox::Timer> t_create_hypre_solver;
   std::shared_ptr<tbox::Timer> t_hypre_solve;
};

#include "QuatLevelSolver.I"

#endif  // included_QuatLevelSolver_h
