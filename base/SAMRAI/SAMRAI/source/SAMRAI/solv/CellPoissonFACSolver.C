/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   High-level solver (wrapper) for scalar poisson equation.
 *
 ************************************************************************/
#ifndef included_solv_CellPoissonFACSolver_C
#define included_solv_CellPoissonFACSolver_C

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/solv/CellPoissonFACSolver.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include IOMANIP_HEADER_FILE

namespace SAMRAI {
namespace solv {

/*
 *************************************************************************
 *
 * Initialize the static data members.
 *
 *************************************************************************
 */

bool CellPoissonFACSolver::s_initialized = 0;
int CellPoissonFACSolver::s_weight_id[tbox::Dimension::
                                      MAXIMUM_DIMENSION_VALUE];
int CellPoissonFACSolver::s_instance_counter[tbox::Dimension::
                                             MAXIMUM_DIMENSION_VALUE];

/*
 *************************************************************************
 *
 * Constructor sets uninitialized solver state.
 * Set default iteration and convergence parameters.
 *
 * By default settings:
 *   - Poisson equation specified has D=1, C=0.
 *   - State is uninitialized
 *   - Logging is disabled
 *   - Context for internal data is set based on object name.
 *
 *************************************************************************
 */

CellPoissonFACSolver::CellPoissonFACSolver(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const boost::shared_ptr<tbox::Database>& database):
   d_dim(dim),
   d_object_name(object_name),
   d_poisson_spec(object_name + "::poisson_spec"),
   d_fac_ops(d_dim, object_name + "::fac_ops"),
   d_fac_precond(object_name + "::fac_precond", d_fac_ops),
   d_bc_object(NULL),
   d_simple_bc(d_dim, object_name + "::bc"),
   d_ln_min(-1),
   d_ln_max(-1),
   d_context(hier::VariableDatabase::getDatabase()->getContext(
                object_name + "::CONTEXT")),
   d_solver_is_initialized(false),
   d_enable_logging(false)
{

   if (!s_initialized) {
      initializeStatics();
   }

   setMaxCycles(10);
   setResidualTolerance(1e-6);
   setPresmoothingSweeps(1);
   setPostsmoothingSweeps(1);
   setCoarseFineDiscretization("Ewing");
#ifdef HAVE_HYPRE
   setCoarsestLevelSolverChoice("hypre");
   setCoarsestLevelSolverTolerance(1e-10);
   setCoarsestLevelSolverMaxIterations(20);
   setUseSMG(true);
#else
   setCoarsestLevelSolverChoice("redblack");
   setCoarsestLevelSolverTolerance(1e-8);
   setCoarsestLevelSolverMaxIterations(500);
#endif

   /*
    * Construct integer tag variables and add to variable database.  Note that
    * variables and patch data indices are shared among all instances.
    * The VariableDatabase holds the variables, once contructed and
    * registered via the VariableDatabase::registerInternalSAMRAIVariable()
    * function call.  Note that variables are registered and patch data indices
    * are made only for the first time through the constructor.
    */
   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   static std::string weight_variable_name("CellPoissonFACSolver_weight");

   boost::shared_ptr<pdat::CellVariable<double> > weight(
      var_db->getVariable(weight_variable_name),
      boost::detail::dynamic_cast_tag());
   if (!weight) {
      weight.reset(
         new pdat::CellVariable<double>(d_dim, weight_variable_name, 1));
   }

   if (s_weight_id[d_dim.getValue() - 1] < 0) {
      s_weight_id[d_dim.getValue() - 1] = var_db->registerInternalSAMRAIVariable(
            weight,
            hier::IntVector::getZero(d_dim));
   }

   /*
    * The default RobinBcCoefStrategy used,
    * SimpleCellRobinBcCoefs only works with constant refine
    * for prolongation.  So we use constant refinement
    * for prolongation by default.
    */
   setProlongationMethod("CONSTANT_REFINE");

   /*
    * The FAC operator optionally uses the preconditioner
    * to get data for logging.
    */
   d_fac_ops.setPreconditioner((const FACPreconditioner *)(&d_fac_precond));

   getFromInput(database);

   s_instance_counter[d_dim.getValue() - 1]++;
}

/*
 *************************************************************************
 *
 * Destructor for CellPoissonFACSolver.
 * Deallocate internal data.
 *
 *************************************************************************
 */

CellPoissonFACSolver::~CellPoissonFACSolver()
{
   s_instance_counter[d_dim.getValue() - 1]--;

   deallocateSolverState();

   if (s_instance_counter[d_dim.getValue() - 1] == 0) {
      hier::VariableDatabase::getDatabase()->
      removeInternalSAMRAIVariablePatchDataIndex(s_weight_id[d_dim.getValue() - 1]);
      s_weight_id[d_dim.getValue() - 1] = -1;
   }
}

/*
 ********************************************************************
 * Set state from database
 *
 * Do not allow FAC preconditioner and Poisson FAC operators to be
 * set from database, as that may cause them to be inconsistent
 * with this object if user does not coordinate the inputs
 * correctly.  This is also why we don't allow direct access to
 * those objects.  The responsibility for maintaining consistency
 * lies in the public functions to set parameters, so use them
 * instead of setting the parameters directly in this function.
 ********************************************************************
 */

void
CellPoissonFACSolver::getFromInput(
   const boost::shared_ptr<tbox::Database>& database)
{
   if (database) {
      if (database->isBool("enable_logging")) {
         bool logging = database->getBool("enable_logging");
         enableLogging(logging);
      }
      if (database->isInteger("max_cycles")) {
         int max_cycles = database->getInteger("max_cycles");
         setMaxCycles(max_cycles);
      }
      if (database->isDouble("residual_tol")) {
         double residual_tol = database->getDouble("residual_tol");
         setResidualTolerance(residual_tol);
      }
      if (database->isInteger("num_pre_sweeps")) {
         int num_pre_sweeps = database->getInteger("num_pre_sweeps");
         setPresmoothingSweeps(num_pre_sweeps);
      }
      if (database->isInteger("num_post_sweeps")) {
         int num_post_sweeps = database->getInteger("num_post_sweeps");
         setPostsmoothingSweeps(num_post_sweeps);
      }
      if (database->isString("coarse_fine_discretization")) {
         std::string s = database->getString("coarse_fine_discretization");
         setCoarseFineDiscretization(s);
      }
      if (database->isString("prolongation_method")) {
         std::string s = database->getString("prolongation_method");
         setProlongationMethod(s);
      }
      if (database->isString("coarse_solver_choice")) {
         std::string s = database->getString("coarse_solver_choice");
         setCoarsestLevelSolverChoice(s);
      }
      if (database->isDouble("coarse_solver_tolerance")) {
         double tol = database->getDouble("coarse_solver_tolerance");
         setCoarsestLevelSolverTolerance(tol);
      }
      if (database->isInteger("coarse_solver_max_iterations")) {
         int itr = database->getInteger("coarse_solver_max_iterations");
         setCoarsestLevelSolverMaxIterations(itr);
      }
#ifdef HAVE_HYPRE
      if (database->isBool("use_smg")) {
         bool smg = database->getBool("use_smg");
         setUseSMG(smg);
      }
#endif
   }
}

/*
 *************************************************************************
 *
 * Prepare internal data for solve.
 * Allocate scratch data.  Create vectors for u and f
 * required by the FACPreconditioner interface.
 * Set up internal boundary condition object.
 * Share data to coordinate with FAC preconditioner and
 * Poisson FAC operator.
 *
 *************************************************************************
 */

void
CellPoissonFACSolver::initializeSolverState(
   const int solution,
   const int rhs,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarse_level,
   const int fine_level)
{
   TBOX_ASSERT(hierarchy);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);

   if (d_bc_object == NULL) {
      TBOX_ERROR(
         d_object_name << ": No BC coefficient strategy object!\n"
                       << "Use either setBoundaries or setPhysicalBcCoefObject\n"
                       << "to specify the boundary conidition.\n");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (solution < 0 || rhs < 0) {
      TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
   }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!hierarchy) {
      TBOX_ERROR(d_object_name << ": NULL hierarchy pointer not allowed\n"
                               << "in inititialization.");
   }
#endif
   d_hierarchy = hierarchy;

   d_ln_min = coarse_level;
   d_ln_max = fine_level;
   if (d_ln_min == -1) {
      d_ln_min = 0;
   }
   if (d_ln_max == -1) {
      d_ln_max = d_hierarchy->getFinestLevelNumber();
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max) {
      TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
                               << "inititialization.\n");
   }
#endif

   int ln;
   for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
      d_hierarchy->getPatchLevel(ln)->allocatePatchData(
         s_weight_id[d_dim.getValue() - 1]);
   }

   d_fac_ops.computeVectorWeights(d_hierarchy,
      s_weight_id[d_dim.getValue() - 1],
      d_ln_min,
      d_ln_max);

   if (d_bc_object == &d_simple_bc) {
      d_simple_bc.setHierarchy(d_hierarchy,
         d_ln_min,
         d_ln_max);
      if (d_poisson_spec.dIsConstant()) {
         d_simple_bc.setDiffusionCoefConstant(d_poisson_spec.getDConstant());
      } else {
         d_simple_bc.setDiffusionCoefId(d_poisson_spec.getDPatchDataId());
      }
   }

   d_fac_ops.setPoissonSpecifications(d_poisson_spec);

   createVectorWrappers(solution, rhs);

   d_fac_precond.initializeSolverState(*d_uv, *d_fv);

   d_solver_is_initialized = true;
}

void
CellPoissonFACSolver::deallocateSolverState()
{
   if (d_hierarchy) {

      d_fac_precond.deallocateSolverState();

      /*
       * Delete internally managed data.
       */
      int ln;
      for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
         d_hierarchy->getPatchLevel(ln)->deallocatePatchData(
            s_weight_id[d_dim.getValue() - 1]);
      }

      d_hierarchy.reset();
      d_ln_min = -1;
      d_ln_max = -1;
      d_solver_is_initialized = false;

      destroyVectorWrappers();

   }
}

void
CellPoissonFACSolver::setBoundaries(
   const std::string& boundary_type,
   const int fluxes,
   const int flags,
   int* bdry_types)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_bc_object != NULL && d_bc_object != &d_simple_bc) {
      TBOX_ERROR(
         d_object_name << ": Bad attempt to set boundary condition\n"
                       << "by using default bc object after it has been overriden.\n");
   }
#endif
   d_simple_bc.setBoundaries(boundary_type,
      fluxes,
      flags,
      bdry_types);
   d_bc_object = &d_simple_bc;
   d_fac_ops.setPhysicalBcCoefObject(d_bc_object);
}

/*
 *************************************************************************
 *
 * Solve the linear system and report whether iteration converged.
 *
 * This version is for an initialized solver state.
 * Before solving, set the final piece of the boundary condition,
 * which is not known until now, and initialize some internal
 * solver quantities.
 *
 *************************************************************************
 */

bool
CellPoissonFACSolver::solveSystem(
   const int u,
   const int f)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!d_solver_is_initialized) {
      TBOX_ERROR(
         d_object_name << ".solveSystem(int,int): uninitialized\n"
                       << "solver state.  You must call initializeSolverState()\n"
                       << "before using this function.  Or you can use\n"
                       << "solveSystem(int,int,...) to initialize the solver,\n"
                       << "solve and deallocate the solver.\n");
   }
   if (u < 0 || f < 0) {
      TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
   }
#endif
   if (d_bc_object == &d_simple_bc) {
      /*
       * Knowing that we are using the SimpelCellRobinBcCoefsX
       * implementation of RobinBcCoefStrategy, we must save
       * the ghost data in u before solving.
       * The solver overwrites it, but SimpleCellRobinBcCoefs
       * needs to get to access it repeatedly.
       */
      d_simple_bc.cacheDirichletData(u);
   }

   createVectorWrappers(u, f);
   bool solver_rval;
   solver_rval = d_fac_precond.solveSystem(*d_uv, *d_fv);

   if (d_bc_object == &d_simple_bc) {
      /*
       * Restore the Dirichlet cell data that were overwritten by the
       * solve process.  We do this to be backward compatible with the
       * user code.
       */
      d_simple_bc.restoreDirichletData(u);
   }

   return solver_rval;
}

/*
 *************************************************************************
 *
 * Solve the linear system and report whether iteration converged.
 *
 * This version is for an uninitialized solver state.
 * 1. Initialize the (currently uninitialized) solver state.
 * 2. Solve.
 * 3. Deallocate the solver state.
 *
 *************************************************************************
 */

bool
CellPoissonFACSolver::solveSystem(
   const int u,
   const int f,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarse_ln,
   int fine_ln)
{
   TBOX_ASSERT(hierarchy);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);

   if (d_enable_logging) {
      tbox::plog << "CellPoissonFACSolver::solveSystem (" << d_object_name
                 << ")\n";
      d_poisson_spec.printClassData(tbox::plog);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_solver_is_initialized) {
      TBOX_ERROR(
         d_object_name << ".solveSystem(int,int,...): initialized\n"
                       << "solver state.  This function can only used when the\n"
                       << "solver state is uninitialized.  You should deallocate\n"
                       << "the solver state or use solveSystem(int,int).\n");
   }
   if (!hierarchy) {
      TBOX_ERROR(d_object_name << ".solveSystem(): Null hierarchy\n"
                               << "specified.\n");
   }
#endif
   initializeSolverState(u, f, hierarchy, coarse_ln, fine_ln);

   bool solver_rval;
   solver_rval = solveSystem(u, f);

   deallocateSolverState();

   return solver_rval;
}

void
CellPoissonFACSolver::createVectorWrappers(
   int u,
   int f)
{
   hier::VariableDatabase& vdb(*hier::VariableDatabase::getDatabase());
   boost::shared_ptr<hier::Variable> variable;

   if (!d_uv || d_uv->getComponentDescriptorIndex(0) != u) {
     d_uv.reset(new SAMRAIVectorReal<double>(d_object_name + "::uv",
                                             d_hierarchy,
                                             d_ln_min,
                                             d_ln_max));
      vdb.mapIndexToVariable(u, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!variable) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                                  << u << "\n");
      }
      boost::shared_ptr<pdat::CellVariable<double> > cell_variable(
         variable,
         boost::detail::dynamic_cast_tag());
      if (!cell_variable) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << u
                                  << " is not a cell-double variable.\n");
      }
#endif
      d_uv->addComponent(variable, u, s_weight_id[d_dim.getValue() - 1]);
   }

   if (!d_fv || d_fv->getComponentDescriptorIndex(0) != f) {
      d_fv.reset(new SAMRAIVectorReal<double>(d_object_name + "::fv",
                                              d_hierarchy,
                                              d_ln_min,
                                              d_ln_max));
      vdb.mapIndexToVariable(f, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!variable) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                                  << f << "\n");
      }
      boost::shared_ptr<pdat::CellVariable<double> > cell_variable(
         variable,
         boost::detail::dynamic_cast_tag());
      if (!cell_variable) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << f
                                  << " is not a cell-double variable.\n");
      }
#endif
      d_fv->addComponent(variable, f, s_weight_id[d_dim.getValue() - 1]);
   }
}

void
CellPoissonFACSolver::initializeStatics()
{
   for (int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      s_weight_id[d] = -1;
      s_instance_counter[d] = -1;
   }

   s_initialized = 1;
}

}
}
#endif
