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
#ifndef included_QuatFACOps_h
#define included_QuatFACOps_h

/*
 * This class provides operator specific functions supporting
 * a quaternion system FAC solver.
 *
 * This file was adapted from solv::CellPoissonFACOps.h in the
 * SAMRAI library.
 */

#ifndef included_SAMRAI_config
#include "SAMRAI/SAMRAI_config.h"
#endif

#include "FACPreconditioner.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideVariable.h"

#include "QuatLevelSolver.h"
#include "CartesianRobinBcHelperWithDepth.h"

#include <vector>

using namespace SAMRAI;

/*!
 * This class implements operator functions for use in the FAC
 * solution of the  quaternion system described in the working notes
 * "Differential Algebraic Formulation of the Pusztai, Bortel and
 * Granasy Quaternion Phase Field Equation of Motion" by M. Dorr
 * on a structured AMR hierarchy.
 *
 * This class provides operators that are used by the FAC solver
 * SAMRAI::solv::FACPreconditioner.  It is used to solve
 * the quaternion system using a cell-centered, second-order,
 * finite-volume discretization.  It is designed to provide all
 * operations required by the FAC algorithm that are specific to
 * the quaternion system.
 *
 * To use the SAMRAI::solv::FACPreconditioner solver with
 * this class, the user must provide:
 * -# The solution vector SAMRAIVectorReal,
 *    with appropriate norm weighting for the cell-centered AMR mesh.
 *    This class provides the function computeVectorWeights()
 *    to help with computing the appropriate weights.
 * -# The source vector SAMRAIVectorReal for f.
 * -# The boundary condition specifications in terms of the coefficients
 *    \alpha, \beta and \gamma in the
 *    Robin formula \alpha u + \beta u_n = \gamma applied on the
 *    boundary faces.  See RobinBcCoefStrategy.
 *
 * This class allocates and deallocates only its own scratch data.
 * Other data that it manipuates are passed in as function arguments.
 * Hence, it owns none of the solution vectors, error vectors,
 * coefficient data, or any such things.
 *
 * The following parameters are settable in the input database.  The default
 * values are shown
 *
 * levelsolver_tolerance = 1e-8            // see
 * setCoarsestLevelSolverTolerance() levelsolver_max_iterations = 10         //
 * see setCoarsestLevelSolverMaxIterations() coarse_levelsolver_tolerance = 1e-8
 * // see setCoarsestLevelSolverTolerance() coarse_levelsolver_max_iterations =
 * 10  // see setCoarsestLevelSolverMaxIterations() cf_discretization = "Ewing"
 * // see setCoarseFineDiscretization() prolongation_method = "CONSTANT_REFINE"
 * // see setProlongationMethod() enable_logging = FALSE                  // see
 * enableLogging()
 */
class QuatFACOps : public SAMRAI::solv::FACOperatorStrategy
{

 public:
   /*!
    * @brief Constructor
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    * @param qlen
    * @param object_name:  Ojbect name
    * @param database:     Input database
    */
   QuatFACOps(const int qlen, const std::string &object_name = std::string(),
              const std::shared_ptr<tbox::Database> &database =
                  std::shared_ptr<tbox::Database>());

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~QuatFACOps(void);

   /*
    * Enable verbosity.
    *
    * This flag controls the printing of solver information
    * to tbox::pout
    */
   void setVerbose(bool verbose);

   /*
    * Functions for setting solver mathematic algorithm controls
    */

   void setLevelSolverTolerance(double tol);

   void setLevelSolverMaxIterations(int max_iterations);

   /*
    * Set tolerance for coarse level solve.
    *
    * If the coarse level solver requires a tolerance (currently, they all do),
    * the specified value is used.
    */
   void setCoarsestLevelSolverTolerance(double tol);

   /*
    * Set max iterations for coarse level solve.
    *
    * If the coarse level solver requires a max iteration limit
    * (currently, they all do), the specified value is used.
    */
   void setCoarsestLevelSolverMaxIterations(int max_iterations);

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
    * not passed to xfer::Geometry::lookupRefineOperator.
    * If this variable is set to "Ewing", Ewing's coarse-fine
    * discretization is used (a constant refinement is performed,
    * and the flux is later corrected to result in Ewing's scheme).
    * For a reference to Ewing's discretization method, see
    * "Local Refinement Techniques for Elliptic Problems on Cell-Centered
    * Grids, I. Error Analysis", Mathematics of Computation, Vol. 56, No. 194,
    * April 1991, pp. 437-461.
    *
    * method: String selecting the coarse-fine discretization method.
    */
   void setCoarseFineDiscretization(const std::string &coarsefine_method);

   /*
    * Set the name of the prolongation method.
    *
    * Specify the op_name string which will be passed to
    * xfer::Geometry::lookupRefineOperator() to get the operator
    * for prolonging the coarse-grid correction.
    *
    * By default, "CONSTANT_REFINE" is used.  "LINEAR_REFINE" seems to
    * to lead to faster convergence, but it does NOT satisfy the Galerkin
    * condition.
    *
    * Prolonging using linear refinement requires a Robin bc
    * coefficient implementation that is capable of delivering
    * coefficients for non-hierarchy data, because linear refinement
    * requires boundary conditions to be set on temporary levels.
    *
    * prolongation_method: String selecting the coarse-fine
    *                      discretization method.
    */
   void setProlongationMethod(const std::string &prolongation_method);


   /*
    * Functions for setting patch data indices and coefficients
    */

   /*
    * Set the scratch patch data index for the flux.
    *
    * The use of this function is optional.
    * The patch data index should be a pdat::SideData type of variable.
    * If the flux id is -1 (the default initial value), scratch space
    * for the flux is allocated as needed and immediately deallocated
    * afterward, level by level.  If you have space preallocated for
    * flux and you would like that to be used, set flux id to the
    * patch data index of that space.
    */
   void setFluxId(int flux_id);


   /*
    * Provide an implementation for getting the physical bc coefficients
    *
    * If your solution is fixed at the physical boundary
    * ghost cell centers AND those cells have the correct
    * values before entering solveSystem(), you may use a
    * GhostCellRobinBcCoefs object.
    *
    * If your solution is not fixed at the ghost cell centers,
    * the ghost cell values will change as the interior
    * cell values change.  In those cases, the flexible
    * Robin boundary conditions are applied.  You must
    * call this function to provide the implementation for
    * determining the boundary condition coefficients.
    *
    * physical_bc_coef: std::shared_ptr to an object that can
    *                   set the Robin bc coefficients.
    */
   void setPhysicalBcCoefObject(
       const solv::RobinBcCoefStrategy *physical_bc_coef);

   /*
    * Functions for checking validity and correctness of state.
    */

   /*
    * Check validity of the flux patch data index.
    */

   void checkFluxPatchDataIndex() const;

   /*
    * Set the FAC preconditioner that will be using this object.
    *
    * The FAC Solver is accessed to get convergence data during
    * the cycle postprocessing step.  It is optional.
    */
   void setSolver(const FACPreconditioner *fac_solver);

   /*
    * Set the operator coefficients.
    */
   void setOperatorCoefficients(
       const double time_step, const double epsilon_q, const int mobility_id,
       const int mobility_deriv_id, const int diff_coef_id,
       const int diff_coef_deriv_id, const int grad_q_id, const int q_id,
       const double gradient_floor, const std::string grad_floor_type);

   // FACOperatorStrategy virtuals

   virtual void restrictSolution(const solv::SAMRAIVectorReal<double> &s,
                                 solv::SAMRAIVectorReal<double> &d,
                                 int dest_ln);

   virtual void restrictResidual(const solv::SAMRAIVectorReal<double> &s,
                                 solv::SAMRAIVectorReal<double> &d,
                                 int dest_ln);

   virtual void prolongErrorAndCorrect(const solv::SAMRAIVectorReal<double> &s,
                                       solv::SAMRAIVectorReal<double> &d,
                                       int dest_ln);

   virtual void smoothError(solv::SAMRAIVectorReal<double> &error,
                            const solv::SAMRAIVectorReal<double> &residual,
                            int ln, int num_sweeps);

   virtual int solveCoarsestLevel(
       solv::SAMRAIVectorReal<double> &error,
       const solv::SAMRAIVectorReal<double> &residual, int ln);

   virtual void computeCompositeResidualOnLevel(
       solv::SAMRAIVectorReal<double> &residual,
       const solv::SAMRAIVectorReal<double> &solution,
       const solv::SAMRAIVectorReal<double> &rhs, int ln,
       bool error_equation_indicator);

   virtual double computeResidualNorm(
       const solv::SAMRAIVectorReal<double> &residual, int fine_ln,
       int coarse_ln);

   virtual void initializeOperatorState(
       const solv::SAMRAIVectorReal<double> &solution,
       const solv::SAMRAIVectorReal<double> &rhs);

   virtual void deallocateOperatorState();

   virtual void postprocessOneCycle(
       int iteration_num, const solv::SAMRAIVectorReal<double> &current_soln,
       const solv::SAMRAIVectorReal<double> &residual);

   /*
    * return ID of the side-centered face diffusion coefficient data.
    */
   int getFaceDiffCoeffId() { return d_face_coef_id; }
   int getFaceDiffCoeffScratchId() { return d_face_coef_scratch_id; }

   void evaluateRHS(const double epsilon_q, const int diff_coef_id,
                    const int grad_q_id, const int grad_q_copy_id,
                    const double gradient_floor,
                    const std::string gradient_floor_type,
                    const int mobility_id, const int rotations_id,
                    const int q_id, int rhs_id, const bool use_gradq_for_flux);

   void accumulateOperatorOnLevel(const int mobility_id, const int face_coef_id,
                                  const int q_id, const int grad_q_id,
                                  int operator_q_id, int ln, bool project,
                                  bool error_equation_indicator);

   void multiplyDQuatDPhiBlock(const int q_id, const int operator_q_id);

   void applyProjectionOnLevel(const int q_id, const int corr_id,
                               const int err_id, const int ln);

   void multiplyMobilitySqrt(const int id);

   void divideMobilitySqrt(const int id);

 private:
   /*
    * Function to compute flux, using general diffusion
    * coefficient data.
    *
    * Recall that this solver class discretizes the PDE
    * [ \nabla \cdot D \nabla u + C u = f ] on an AMR grid.  This member
    * function allows users of this solver class to compute gradient
    * terms, D \nabla w , in their code in a manner consistent with the
    * solver discretization.   In particular, when solving PDE systems, it may
    * be necessary to discretize the gradient operator appearing in equations
    * not treated by the solver class in the same way as those treated by this
    * class.  These funtions allow users to do this easily.  The divergence
    * operator used in this solver is the standard sum of centered differences
    * involving flux terms on the cell sides computed by these routines.
    *
    * Note that the patch must exist on a level in an AMR hierarchy so that
    * the discretization can be computed properly at the coarse-fine interface.
    * Poisson coefficients C and D must exist on the patch, if they are
    * variable. Also, calling this function does not affect the internal solver
    * state in any way.  However, the solver must be fully initialized before it
    * is called and care should be exercised to pass arguments so that the
    * solver solution quantity and other internal solver quantities are not
    * adversely affected.
    *
    * patch:       patch on which computation will take place
    * ratio_to_coarser_level: refinement ratio from coarser level to level
    *                               on which patch lives; if current patch level
    *                               is level zero, this is ignored
    * w_data:      cell-centered data
    * Dgradw_data: side-centered flux data (i.e., D (grad w))
    */
   void computeFluxOnPatch(const hier::Patch &patch,
                           const hier::IntVector &ratio_to_coarser_level,
                           const pdat::SideData<double> &face_coef_data,
                           const pdat::CellData<double> &q_data,
                           pdat::SideData<double> &flux_data) const;

   /*!
    * compute flux using gradient of q at sides
    */
   void computeFluxOnPatch(const hier::Patch &patch,
                           const hier::IntVector &ratio_to_coarser_level,
                           const pdat::SideData<double> &face_coef_data,
                           const pdat::SideData<double> &grad_q_data,
                           pdat::SideData<double> &flux_data) const;

   void computeSymmetricFluxOnPatch(
       const hier::Patch &patch, const hier::IntVector &ratio_to_coarser_level,
       const pdat::SideData<double> &face_coef_data,
       const pdat::CellData<double> &sqrt_m_data,
       const pdat::CellData<double> &q_data,
       pdat::SideData<double> &flux_data) const;

   void computeLambdaOnPatch(
       const hier::Patch &patch, const pdat::SideData<double> &flux_data,
       const pdat::CellData<double> &q_data,
       std::shared_ptr<pdat::SideData<int> > rotation_index,
       pdat::CellData<double> &lambda_data) const;

   void computeFaceCoefs(const double epsilon_q, const int diff_coef_id,
                         const int grad_q_id, const double gradient_floor,
                         const std::string grad_floor_type,
                         const int face_coef_id);

   void computeFaceCoefsOnPatch(const hier::Patch &patch,
                                const double epsilon_q,
                                pdat::SideData<double> &diff_coef_data,
                                pdat::SideData<double> &grad_q_data,
                                pdat::SideData<double> &face_coef_data,
                                const double gradient_floor,
                                const std::string grad_floor_type) const;

   void computeDQuatDPhiFaceCoefs(const int dprime_id, const int phi_id,
                                  const int face_coef_id);

   void computeDQuatDPhiFaceCoefsOnPatch(
       const hier::Patch &patch, pdat::SideData<double> &dprime_data,
       pdat::CellData<double> &phi_data,
       pdat::SideData<double> &face_coef_data) const;

   /*
    * Solves the residual equation Ae=r on a level.
    *
    * error:              Error vector
    * residual:           Residual vector
    * ln:                 Level number
    * num_sweeps:         Number of sweeps
    * residual_tolerance: The maximum residual considered to be
    *                     converged
    */
   void doLevelSolve(solv::SAMRAIVectorReal<double> &error,
                     const solv::SAMRAIVectorReal<double> &residual, int ln,
                     int num_sweeps, double residual_tolerance = -1.0);

   /*
    * Fix flux per Ewing's coarse-fine boundary treatment.
    *
    * Ewing's coarse-fine boundary treatment can be implemented
    * using a constant refinement into the fine-grid ghost boundary,
    * naively computing the flux using the constant-refined data then
    * fixing up the flux to correct the error.
    *
    * To use this function
    * -# you must use constant refinement to fill the fine level ghost cells
    * -# the flux must first be computed and stored
    *
    * patch:            patch
    * soln_data:        cell-centered solution data
    * flux_data:        side-centered flux data
    * facecoef_data:    side-centered diffusion coefficient data
    * cfb:              coarse-fine boundary object for the level
    *                   in which patch resides
    * ratio_to_coarser: Refinement ratio to the next coarser level.
    */
   void ewingFixFlux(const hier::Patch &patch,
                     const pdat::CellData<double> &soln_data,
                     const pdat::SideData<double> &face_coef_data,
                     pdat::SideData<double> &flux_data,
                     const hier::IntVector &ratio_to_coarser) const;

   /*
    * AMR-unaware function to compute residual on a single patch,
    * with variable scalar field.
    *
    * patch: patch
    * flux_data:     side-centered flux data
    * soln_data:     cell-centered solution data
    * rhs_data:      cell-centered rhs data
    * residual_data: cell-centered residual data
    */
   void computeResidualOnPatch(
       const hier::Patch &patch, const pdat::SideData<double> &flux_data,
       std::shared_ptr<pdat::SideData<int> > rotations,
       const pdat::CellData<double> &sqrt_m_data,
       const pdat::CellData<double> &q_soln_data,
       const pdat::CellData<double> &q_rhs_data,
       pdat::CellData<double> &q_residual_data) const;

   void accumulateOperatorOnPatch(
       const hier::Patch &patch, const pdat::SideData<double> &flux_data,
       const pdat::CellData<double> &mobility_data,
       const pdat::CellData<double> &q_rhs_data) const;

   void accumulateProjectedOperatorOnPatch(
       const hier::Patch &patch, const pdat::SideData<double> &flux_data,
       const pdat::CellData<double> &mobility_data,
       const pdat::CellData<double> &q_soln_data,
       const pdat::CellData<double> &lambda_soln_data,
       std::shared_ptr<pdat::SideData<int> > rotations,
       const pdat::CellData<double> &q_rhs_data) const;

   void takeSquareRootOnPatch(pdat::CellData<double> &data);

   // For executing, caching and resetting communication schedules.

   /*
    * Execute a refinement schedule for prolonging cell data.
    *
    * General notes regarding internal objects for communication:
    * We maintain objects to support caching schedules to improve
    * efficiency.  Communication is needed in 5 distinct tasks.
    *   -# Prolongation
    *   -# Restriction
    *   -# Flux coarsening.  Changing the coarse grid flux to the
    *      composite grid flux by coarsening the fine grid flux
    *      at the coarse-fine boundaries.
    *   -# Fill boundary data from other patches in the same level
    *      and physical boundary condition.
    *   -# Fill boundary data from same level, coarser levels
    *      and physical boundary condition.
    *
    * For each task, we maintain a refine or coarsen operator,
    * a cache (array) of communication schedules and a list of keys
    * assocating the schedules to the exact data being transfered.
    * A key is a single integer uniquely identifying the data ids
    * and destination level of the transfer, allowing a simple 1D
    * lookup of the cached schedules.
    *
    * The 5 member functions named xeqSchedule... return
    * communication schedules appropriate for five specific tasks.
    * They use a cached schedule if possible or create and cache
    * a new schedule if needed.  These functions and the data
    * they manipulate are as follows:
    *
    *   xeqScheduleProlongation():
    *        d_prolongation_refine_operator
    *        d_prolongation_refine_schedules
    *   xeqScheduleURestriction():
    *        d_restriction_coarsen_operator,
    *        d_urestriction_coarsen_schedules.
    *   xeqScheduleRRestriction():
    *        d_restriction_coarsen_operator,
    *        d_rrestriction_coarsen_schedules.
    *   xeqScheduleFluxCoarsen():
    *        d_flux_coarsen_operator,
    *        d_flux_coarsen_schedules.
    *   xeqScheduleGhostFill():
    *        d_ghostfill_refine_operator,
    *        d_ghostfill_refine_schedules.
    *   xeqScheduleGhostFillNoCoarse():
    *        d_ghostfill_nocoarse_refine_operator,
    *        d_ghostfill_nocoarse_refine_schedules.
    *
    * Returns refinement schedule for prolongation
    */
   void xeqScheduleProlongation(int dst_id, int src_id, int scr_id,
                                int dest_ln);


   /*
    * Execute schedule for restricting solution to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * Returns coarsening schedule for restriction
    */
   void xeqScheduleURestriction(int dst_id, int src_id, int dest_ln);


   /*
    * Execute schedule for restricting residual to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * Returns coarsening schedule for restriction
    */
   void xeqScheduleRRestriction(int dst_id, int src_id, int dest_ln);


   /*
    * Execute schedule for coarsening flux to the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * Returns coarsening schedule for setting composite grid flux at
    * coarse-fine boundaries.
    */
   void xeqScheduleFluxCoarsen(int dst_id, int src_id, int dest_ln);


   /*
    * Execute schedule for filling ghosts on the specified
    * level or reregister an existing one.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * Returns refine schedule for filling ghost data from coarser level
    * and physical bc.
    */
   void xeqScheduleGhostFill(int dst_id, int dest_ln);


   /*
    * Execute schedule for filling ghosts on the specified level
    * or reregister an existing one.  This version does not get
    * data from coarser levels.
    *
    * See general notes for xeqScheduleProlongation().
    *
    * This function is used for the bottom solve level, since it does
    * not access data from any coarser level.  (Ghost data obtained
    * from coarser level must have been placed there before solve begins!)
    *
    * Returns refine schedule for filling ghost data from same level
    * and physical bc.
    */
   void xeqScheduleGhostFillNoCoarse(int dst_id, int dest_ln);

   // Return the patch data index for cell scratch data.
   int registerCellScratch() const;

   // Return the patch data index for flux scratch data.
   int registerFluxScratch() const;

   // Return the patch data index for outerflux scratch data.
   int registerOfluxScratch() const;

   // Free static variables at shutdown time.
   static void freeVariables();

   const int d_qlen;

   /*
    * Object name.
    */
   std::string d_object_name;


   // Hierarchy-dependent objects.

   /*
    * Reference hierarchy
    *
    * This variable is non-null between the initializeOperatorState()
    * and deallocateOperatorState() calls.  It is not truly needed,
    * because the hierarchy is obtainable through variables in most
    * function argument lists.  We use it to enforce working on one
    * hierarchy at a time.
    */
   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*
    * Coarsest level for solve.
    */
   int d_ln_min;

   /*
    * Finest level for solve.
    */
   int d_ln_max;

   /*
    * Description of coarse-fine boundaries.
    *
    * There is one coarse-fine boundary object for each level.
    * d_coarse_fine_boundary[i] is the description of
    * the coarse-fine boundary between level i and level i-1.
    * The coarse-fine boundary does not exist at the coarsest level,
    * although the hier::CoarseFineBoundary object still exists (it
    * should not contain any boxes).
    *
    * This array is initialized in initializeOperatorState() and
    * deallocated in deallocateOperatorState().  When allocated,
    * it is allocated for the index range [0,d_ln_max], though
    * the range [0,d_ln_min-1] is not used.  This is okay because
    * hier::CoarseFineBoundary is a light object before
    * it is set for a level.
    */
   std::vector<std::shared_ptr<hier::CoarseFineBoundary> > d_cf_boundary;


   /*
    * Private state variables for solution process.
    */

   /*
    * Coarse-fine discretization method.
    * See setCoarseFineDiscretization().
    */
   std::string d_cf_discretization;

   /*
    * Coarse-fine discretization method.
    *
    * The name of the refinement operator used to prolong the
    * coarse grid correction.
    *
    * See setProlongationMethod()
    */
   std::string d_prolongation_method;

   double d_levelsolver_tolerance;
   int d_levelsolver_max_iterations;

   /*
    * Tolerance specified to coarse solver
    * See setCoarsestLevelSolverTolerance()
    */
   double d_coarse_levelsolver_tolerance;

   /*
    * Coarse level solver iteration limit.
    * See setCoarsestLevelSolverMaxIterations()
    */
   int d_coarse_levelsolver_max_iterations;

   /*
    * Id of the flux.
    *
    * If set to -1, create and delete storage space on the fly.
    * Else, user has provided space for flux.
    *
    * See setFluxId
    */
   int d_flux_id;

   /*
    * Quaternion system solver objects on each level.
    */
   std::vector<QuatLevelSolver *> d_quat_level_solver;

   /*
    * Level solver input database.
    */
   std::shared_ptr<tbox::Database> d_levelsolver_database;

   /*
    * Externally provided physical boundary condition object.
    *
    * see setPhysicalBcCoefObject()
    */
   const solv::RobinBcCoefStrategy *d_physical_bc_coef;

   static std::shared_ptr<pdat::CellVariable<double> > s_cell_scratch_var;

   static std::shared_ptr<pdat::SideVariable<double> > s_flux_scratch_var;

   static std::shared_ptr<pdat::OutersideVariable<double> >
       s_oflux_scratch_var;

   static std::shared_ptr<pdat::SideVariable<double> > s_face_coef_var;

   static std::shared_ptr<pdat::SideVariable<double> > s_face_coef_deriv_var;

   static std::shared_ptr<pdat::CellVariable<double> > s_q_local_var;

   static std::shared_ptr<pdat::CellVariable<double> > s_residual_var;

   static std::shared_ptr<pdat::CellVariable<double> > s_sqrt_m_var;

   static std::shared_ptr<pdat::CellVariable<double> > s_m_deriv_var;

   /*
    * ID of the solution-like scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::CellVariable<double> named
    * d_object_name+"::cell_scratch".
    * Scratch data is allocated and removed as needed
    * to reduce memory usage.
    */
   int d_cell_scratch_id;

   /*
    * ID of the side-centered scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::SideVariable<double> named
    * d_object_name+"::flux_scratch".
    *
    * This data is allocated only as needed and deallocated
    * immediately after use.
    */
   int d_flux_scratch_id;

   /*
    * ID of the outerside-centered scratch data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::OutersideVariable<double> named
    * d_object_name+"::oflux_scratch".
    */
   int d_oflux_scratch_id;

   /*
    * ID of the side-centered face diffusion coefficient data.
    *
    * Set in setOperatorCoefficients.
    * Corresponds to a pdat::SideVariable<double> named
    * d_object_name+"::face_coef".
    */
   int d_face_coef_id;
   int d_face_coef_deriv_id;
   int d_face_coef_scratch_id;

   /*
    * ID of the cell-centered local q data.
    *
    * Set in constructor and never changed.
    * Corresponds to a pdat::CellVariable<double> named
    * d_object_name+"::face_q_scratch".
    *
    * This data is allocated only as needed and deallocated
    * immediately after use.
    */
   int d_q_local_id;

   int d_residual_id;

   /*
    * Various refine and coarsen objects used internally.
    */

   // Error prolongation (refinement) operator.
   std::shared_ptr<hier::RefineOperator> d_prolongation_refine_operator;
   std::shared_ptr<xfer::RefineAlgorithm> d_prolongation_refine_algorithm;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
       d_prolongation_refine_schedules;

   // Solution restriction (coarsening) operator.
   std::shared_ptr<hier::CoarsenOperator> d_urestriction_coarsen_operator;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_urestriction_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_urestriction_coarsen_schedules;

   // Residual restriction (coarsening) operator.
   std::shared_ptr<hier::CoarsenOperator> d_rrestriction_coarsen_operator;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_rrestriction_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_rrestriction_coarsen_schedules;

   // Coarsen operator for outerflux-to-flux
   std::shared_ptr<hier::CoarsenOperator> d_flux_coarsen_operator;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_flux_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_flux_coarsen_schedules;

   // Refine operator for cell-like data from coarser level.
   std::shared_ptr<hier::RefineOperator> d_ghostfill_refine_operator;
   std::shared_ptr<xfer::RefineAlgorithm> d_ghostfill_refine_algorithm;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
       d_ghostfill_refine_schedules;

   // Refine operator for cell-like data from same level.
   std::shared_ptr<hier::RefineOperator> d_ghostfill_nocoarse_refine_operator;
   std::shared_ptr<xfer::RefineAlgorithm>
       d_ghostfill_nocoarse_refine_algorithm;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
       d_ghostfill_nocoarse_refine_schedules;


   /*
    * Utility object employed in setting ghost cells and providing
    * xfer::RefinePatchStrategy implementation.
    *
    * Since this class deals only in scalar variables having
    * Robin boundary conditions, we take advantage of the corresponding
    * implementation in CartesianRobinBcHelper.  Whenever
    * we need an implementation of xfer::RefinePatchStrategy,
    * this object is used.  Note that in the code, before we
    * use this object to set ghost cell values, directly or
    * indirectly by calling xfer::RefineSchedule::fillData(),
    * we must tell d_bc_helper the patch data index we want
    * to set and whether we are setting data with homogeneous
    * boundary condition.
    */
   CartesianRobinBcHelperWithDepth d_bc_helper;


   /*
    * Non-essential objects used in outputs and debugging.
    */

   /*
    * Logging flag.
    */
   bool d_enable_logging;

   /*
    * Verbosity flag.
    */
   bool d_verbose;

   /*
    * Preconditioner using this object.
    *
    * See setPreconditioner().
    */
   const FACPreconditioner *d_solver;

   /*
    * Hierarchy cell operator
    */
   std::shared_ptr<math::HierarchyCellDataOpsReal<double> > d_hopscell;

   /*
    * Timers for performance measurement.
    */
   std::shared_ptr<tbox::Timer> t_restrict_solution;
   std::shared_ptr<tbox::Timer> t_restrict_residual;
   std::shared_ptr<tbox::Timer> t_prolong;
   std::shared_ptr<tbox::Timer> t_smooth_error;
   std::shared_ptr<tbox::Timer> t_solve_coarsest;
   std::shared_ptr<tbox::Timer> t_compute_composite_residual;
   std::shared_ptr<tbox::Timer> t_compute_rhs;
   std::shared_ptr<tbox::Timer> t_compute_residual_norm;

   double d_gamma;

   int d_sqrt_m_id;
   int d_m_deriv_id;

   int GetNumCellFacesInBox(const int *lower, const int *upper,
                            const int dim) const;

   /*!
    * rotation indexes used to rotate gradient on sides if
    * (i) fluxes are computed using gradients
    * (ii) symmetry is ON
    * Otherwise id set to -1 when calling evaluateRHS() function
    */
   int d_rotation_index_id;
};

#include "QuatFACOps.I"

#endif  // included_QuatFACOps_h
