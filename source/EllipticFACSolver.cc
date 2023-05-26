/*************************************************************************
 *
 * This file is adapted from the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE at https://github.com/LLNL/SAMRAI.
 *
 * Copyright:     (c) 1997-2021 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Poisson equation
 *
 ************************************************************************/
#include "EllipticFACSolver.h"
#include "EllipticFACOps.h"

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>


/*
*************************************************************************
*                                                                       *
* Constructor sets uninitialized solver state.                          *
* Set default iteration and convergence parameters.                     *
*                                                                       *
* By default settings:                                                  *
*   - Poisson equation specified has D=1, C=0.                          *
*   - State is uninitialized                                            *
*   - Logging is disabled                                               *
*   - Context for internal data is set based on object name.            *
*                                                                       *
*************************************************************************
*/
EllipticFACSolver::EllipticFACSolver(
    const std::string& object_name, std::shared_ptr<EllipticFACOps> fac_ops,
    const std::shared_ptr<tbox::Database>& database)
    : d_object_name(object_name),
      d_fac_ops(fac_ops),
      d_fac_precond(object_name + "::fac_precond", d_fac_ops, database),
      d_bc_object(NULL),
      d_simple_bc(tbox::Dimension(NDIM), object_name + "::bc"),
      d_ln_min(-1),
      d_ln_max(-1),
      d_context(hier::VariableDatabase::getDatabase()->getContext(object_name +
                                                                  "::CONTEXT")),
      d_solver_is_initialized(false),
      d_enable_logging(false),
      d_verbose(false)
{
   /*
    * The FAC operator optionally uses the preconditioner
    * to get data for logging.
    */
   d_fac_ops->setPreconditioner((const FACPreconditioner*)(&d_fac_precond));

   if (database) {
      getFromInput(database);
   }
}

/*
*************************************************************************
*                                                                       *
* Destructor for EllipticFACSolver.                                *
* Deallocate internal data.                                             *
*                                                                       *
*************************************************************************
*/
EllipticFACSolver::~EllipticFACSolver() { deallocateSolverState(); }


/*
********************************************************************
* Set state from database                                          *
*                                                                  *
* Do not allow FAC preconditioner and Poisson FAC operators to be  *
* set from database, as that may cause them to be inconsistent     *
* with this object if user does not coordinate the inputs          *
* correctly.  This is also why we don't allow direct access to     *
* those objects.  The responsibility for maintaining consistency   *
* lies in the public functions to set parameters, so use them      *
* instead of setting the parameters directly in this function.     *
********************************************************************
*/
void EllipticFACSolver::getFromInput(
    const std::shared_ptr<tbox::Database>& database)
{
   if (database->isBool("enable_logging")) {
      d_enable_logging = database->getBool("enable_logging");
   }
}

/*
*************************************************************************
*                                                                       *
* Prepare internal data for solve.                                      *
* Allocate scratch data.  Create std::vectors for u and f                    *
* required by the FACPreconditioner interface.                    *
* Set up internal boundary condition object.                            *
* Share data to coordinate with FAC preconditioner and                  *
* Poisson FAC operator.                                                 *
*                                                                       *
*************************************************************************
*/
void EllipticFACSolver::initializeSolverState(
    const int solution, const int rhs,
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarse_level, const int fine_level)
{
   TBOX_ASSERT(hierarchy);

   if (d_bc_object == NULL) {
      TBOX_ERROR(d_object_name << ": No BC coefficient strategy object!\n"
                               << "Use either setBoundaries or "
                                  "setPhysicalBcCoefObject\n"
                               << "to specify the boundary conidition.\n");
   }

   if (!d_solver_is_initialized) {
      if (solution < 0 || rhs < 0) {
         TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
      }

      if (!hierarchy) {
         TBOX_ERROR(d_object_name << ": NULL hierarchy pointer not allowed\n"
                                  << "in inititialization.");
      }
      d_hierarchy = hierarchy;

      d_ln_min = coarse_level;
      d_ln_max = fine_level;
      if (d_ln_min == -1) {
         d_ln_min = 0;
      }
      if (d_ln_max == -1) {
         d_ln_max = d_hierarchy->getFinestLevelNumber();
      }

      if (d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max) {
         TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
                                  << "inititialization.\n");
      }

      createVectorWrappers(solution, rhs);

      d_fac_precond.initializeSolverState(*d_uv, *d_fv);

      d_solver_is_initialized = true;
   }
}


void EllipticFACSolver::finalizeCoefficients()
{

   if (d_bc_object == &d_simple_bc) {
      d_simple_bc.setHierarchy(d_hierarchy, d_ln_min, d_ln_max);
      if (d_fac_ops->dIsConstant()) {
         d_simple_bc.setDiffusionCoefConstant(d_fac_ops->getDConstant());
      } else {
         d_simple_bc.setDiffusionCoefId(d_fac_ops->getDPatchDataId());
      }
   }

   d_fac_ops->finalizeCoefficients();
}

void EllipticFACSolver::deallocateSolverState()
{
   if (d_hierarchy) {

      d_fac_precond.deallocateSolverState();

      d_hierarchy.reset();
      d_ln_min = -1;
      d_ln_max = -1;
      d_solver_is_initialized = false;

      destroyVectorWrappers();
   }
   return;
}


void EllipticFACSolver::resetSolverState(
    const int soln_id, const int rhs_id,
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   if (d_solver_is_initialized) {
      assert(hierarchy);
      assert(soln_id >= 0);
      assert(rhs_id >= 0);
      deallocateSolverState();
      initializeSolverState(soln_id, rhs_id, hierarchy);
   }
}


void EllipticFACSolver::setBoundaries(const std::string& boundary_type,
                                      const int fluxes, const int flags,
                                      int* bdry_types)
{
   if (d_bc_object != NULL && d_bc_object != &d_simple_bc) {
      TBOX_ERROR(d_object_name << ": Bad attempt to set boundary condition\n"
                               << "by using default bc object after it has "
                                  "been overriden.\n");
   }
   d_simple_bc.setBoundaries(boundary_type, fluxes, flags, bdry_types);
   d_bc_object = &d_simple_bc;
   d_fac_ops->setPhysicalBcCoefObject(d_bc_object);
}


/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an initialized solver state.                      *
* Before solving, set the final piece of the boundary condition,        *
* which is not known until now, and initialize some internal            *
* solver quantities.                                                    *
*                                                                       *
*************************************************************************
*/
bool EllipticFACSolver::solveSystem(const int u_id, const int f_id)
{
   // tbox::pout<<"EllipticFACSolver::solveSystem() for object
   // "<<d_object_name<<endl;
   assert(u_id != -1);
   assert(f_id != -1);
   if (!d_solver_is_initialized) {
      TBOX_ERROR(d_object_name << ".solveSystem(int,int): uninitialized\n"
                               << "solver state.  You must call "
                                  "initializeSolverState()\n"
                               << "before using this function.  Or you can "
                                  "use\n"
                               << "solveSystem(int,int,...) to initialize the "
                                  "solver,\n"
                               << "solve and deallocate the solver.\n");
   }
   if (u_id < 0 || f_id < 0) {
      TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
   }
   if (d_bc_object == &d_simple_bc) {
      /*
       * Knowing that we are using the SimpelCellRobinBcCoefsX
       * implementation of RobinBcCoefStrategy, we must save
       * the ghost data in u before solving.
       * The solver overwrites it, but SimpleCellRobinBcCoefs
       * needs to get to access it repeatedly.
       */
      d_simple_bc.cacheDirichletData(u_id);
   }

   createVectorWrappers(u_id, f_id);

   bool solver_rval = d_fac_precond.solveSystem(*d_uv, *d_fv);

   if (d_bc_object == &d_simple_bc) {
      /*
       * Restore the Dirichlet cell data that were overwritten by the
       * solve process.  We do this to be backward compatible with the
       * user code.
       */
      d_simple_bc.restoreDirichletData(u_id);
   }

   if (d_verbose) printFACConvergenceFactors(solver_rval);

   return solver_rval;
}


/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an uninitialized solver state.                    *
* 1. Initialize the (currently uninitialized) solver state.             *
* 2. Solve.                                                             *
* 3. Deallocate the solver state.                                       *
*                                                                       *
*************************************************************************
*/
bool EllipticFACSolver::solveSystem(
    const int u_id, const int f_id,
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int coarse_ln,
    int fine_ln)
{
   TBOX_ASSERT(hierarchy);

   if (d_enable_logging) {
      tbox::plog << "EllipticFACSolver::solveSystem (" << d_object_name
                 << ")\n";
   }
   if (d_solver_is_initialized) {
      TBOX_ERROR(d_object_name << ".solveSystem(int,int,...): initialized\n"
                               << "solver state.  This function can only used "
                                  "when the\n"
                               << "solver state is uninitialized.  You should "
                                  "deallocate\n"
                               << "the solver state or use "
                                  "solveSystem(int,int).\n");
   }
   if (!hierarchy) {
      TBOX_ERROR(d_object_name << ".solveSystem(): Null hierarchy\n"
                               << "specified.\n");
   }
   initializeSolverState(u_id, f_id, hierarchy, coarse_ln, fine_ln);

   bool solver_rval = solveSystem(u_id, f_id);

   deallocateSolverState();

   return solver_rval;
}


void EllipticFACSolver::createVectorWrappers(int u, int f)
{
   // vol_id=-1 sets SAMRAIVectorReal without control volume array
   int vol_id = -1;

   hier::VariableDatabase& vdb(*hier::VariableDatabase::getDatabase());
   std::shared_ptr<hier::Variable> variable;

   if (!d_uv || d_uv->getComponentDescriptorIndex(0) != u) {
      d_uv.reset(new solv::SAMRAIVectorReal<double>(d_object_name + "::uv",
                                                    d_hierarchy, d_ln_min,
                                                    d_ln_max));
      vdb.mapIndexToVariable(u, variable);
      if (!variable) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index " << u
                                  << "\n");
      }
      std::shared_ptr<pdat::CellVariable<double> > cell_variable(
          SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(
              variable));
      if (!cell_variable) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << u
                                  << " is not a cell-double variable.\n");
      }
      d_uv->addComponent(variable, u, vol_id);
   }

   if (!d_fv || d_fv->getComponentDescriptorIndex(0) != f) {
      d_fv.reset(new solv::SAMRAIVectorReal<double>(d_object_name + "::fv",
                                                    d_hierarchy, d_ln_min,
                                                    d_ln_max));
      vdb.mapIndexToVariable(f, variable);
      if (!variable) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index " << f
                                  << "\n");
      }
      std::shared_ptr<pdat::CellVariable<double> > cell_variable(
          SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(
              variable));
      if (!cell_variable) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << f
                                  << " is not a cell-double variable.\n");
      }
      d_fv->addComponent(variable, f, vol_id);
   }
}


void EllipticFACSolver::printFACConvergenceFactors(const int solver_ret)
{
   d_fac_precond.printClassData(tbox::pout);
   double avg_factor, final_factor;
   getConvergenceFactors(avg_factor, final_factor);
   tbox::pout << "  EllipticFACSolver " << d_object_name << ", iteration ";
   tbox::pout << (solver_ret ? "" : "NOT ") << "converged "
              << "\n"
              << "     iterations: " << getNumberOfIterations() << "\n"
              << "     residual: " << getResidualNorm() << "\n"
              << "     average convergence: " << avg_factor << "\n"
              << "     final convergence: " << final_factor << "\n"
              << std::flush;
}
