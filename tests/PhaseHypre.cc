/*************************************************************************
 * Adapted from SAMRAI test suite
 ************************************************************************/
#include "PhaseHypre.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/SimpleCellRobinBcCoefs.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "FuncFort.h"

#include <string>


extern "C" {
void SAMRAI_F77_FUNC(phasesetexactandrhs2d, PHASESETEXACTANDRHS2D)(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1, double* exact, double* rhs, const double* dx,
    const double* xlower, const double& mobility, const double& delta,
    const double& omega, const double& gamma);

void SAMRAI_F77_FUNC(phasesetexactandrhs3d, PHASESETEXACTANDRHS3D)(
    const int& ifirst0, const int& ilast0, const int& ifirst1,
    const int& ilast1, const int& ifirst2, const int& ilast2, double* exact,
    double* rhs, const double* dx, const double* xlower, const double& mobility,
    const double& delta, const double& omega, const double& gamma);
}

/*
 *************************************************************************
 * Constructor creates a unique context for the object and register
 * all its internal variables with the variable database.
 *************************************************************************
 */
PhaseHypre::PhaseHypre(
    const std::string& object_name, const tbox::Dimension& dim,
    const std::shared_ptr<solv::LocationIndexRobinBcCoefs>& bc_coefs,
    std::shared_ptr<tbox::Database> database, const double epsilon,
    const double omega, const double delta, const double mobility_val,
    const double gamma)
    : d_object_name(object_name),
      d_dim(dim),
      d_poisson_solver(object_name + "::hypre_solver",
                       database && database->isDatabase("hypre_solver")
                           ? database->getDatabase("hypre_solver")
                           : std::shared_ptr<tbox::Database>()),
      d_epsilon(epsilon),
      d_omega(omega),
      d_delta(delta),
      d_mobility(mobility_val),
      d_gamma(gamma),
      d_poisson_spec(object_name + "::Poisson specs"),
      d_C_is_set(false),
      d_D_is_set(false),
      d_M_is_set(false)
{
   d_physical_bc_coef = bc_coefs;

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   /*
    * Get a unique context for variables owned by this object.
    */
   d_context = vdb->getContext(d_object_name + ":Context");

   /*
    * Register variables with hier::VariableDatabase
    * and get the descriptor indices for those variables.
    */
   std::shared_ptr<pdat::CellVariable<double> > comp_soln(
       new pdat::CellVariable<double>(dim, object_name + ":computed solution",
                                      1));
   d_comp_soln_id = vdb->registerVariableAndContext(
       comp_soln, d_context,
       hier::IntVector(dim, 1) /* ghost cell width is 1 for stencil widths */);

   std::shared_ptr<pdat::CellVariable<double> > exact_solution(
       new pdat::CellVariable<double>(dim, object_name + ":exact solution"));
   d_exact_id = vdb->registerVariableAndContext(
       exact_solution, d_context,
       hier::IntVector(dim, 1) /* ghost cell width is 1 in case needed */);

   std::shared_ptr<pdat::CellVariable<double> > rhs_variable(
       new pdat::CellVariable<double>(dim, object_name + ":linear system right "
                                                         "hand side"));
   d_rhs_id = vdb->registerVariableAndContext(
       rhs_variable, d_context,
       hier::IntVector(dim, 0) /* ghost cell width is 0 */);

   std::shared_ptr<pdat::CellVariable<double> > mobility(
       new pdat::CellVariable<double>(dim, object_name + ":mobility"));
   d_m_id = vdb->registerVariableAndContext(
       mobility, d_context,
       hier::IntVector(dim, 1) /* ghost cell width is 1 */);

   std::shared_ptr<pdat::CellVariable<double> > cdata(
       new pdat::CellVariable<double>(dim, object_name + ":cdata"));

   d_c_id = vdb->registerVariableAndContext(
       cdata, d_context, hier::IntVector(dim, 0) /* ghost cell width is 0 */);
}

/*
 *************************************************************************
 * Destructor does nothing interesting
 *************************************************************************
 */
PhaseHypre::~PhaseHypre() { d_poisson_solver.deallocateSolverState(); }

/*
 *************************************************************************
 * Initialize data on a level.
 *
 * Allocate the solution, exact solution and rhs memory.
 * Fill the rhs and exact solution.
 *************************************************************************
 */
void PhaseHypre::initializeLevelData(
    const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number, const double init_data_time,
    const bool can_be_refined, const bool initial_time,
    const std::shared_ptr<hier::PatchLevel>& old_level,
    const bool allocate_data)
{
   NULL_USE(init_data_time);
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(old_level);

   assert(d_m_id >= 0);
   assert(d_mobility > 0.);

   std::shared_ptr<hier::PatchHierarchy> hierarchy = patch_hierarchy;

   std::shared_ptr<hier::PatchLevel> level(
       hierarchy->getPatchLevel(level_number));

   if (allocate_data) {
      level->allocatePatchData(d_comp_soln_id);
      level->allocatePatchData(d_rhs_id);
      level->allocatePatchData(d_exact_id);
      level->allocatePatchData(d_m_id);
      level->allocatePatchData(d_c_id);
   }

   /*
    * Initialize data in all patches in the level.
    */
   for (hier::PatchLevel::iterator pi(level->begin()); pi != level->end();
        ++pi) {

      const std::shared_ptr<hier::Patch>& patch = *pi;
      if (!patch) {
         TBOX_ERROR(d_object_name << ": Cannot find patch.  Null patch "
                                     "pointer.");
      }
      hier::Box pbox = patch->getBox();
      std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                 hier::PatchGeometry>(
              patch->getPatchGeometry()));
      TBOX_ASSERT(patch_geom);

      std::shared_ptr<pdat::CellData<double> > mobility(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_m_id)));
      TBOX_ASSERT(mobility);
      mobility->fillAll(d_mobility);
#ifdef DEBUG_CHECK_ASSERTIONS
      math::PatchCellDataNormOpsReal<double> ops;
      const double norm_Mg = ops.maxNorm(mobility, mobility->getGhostBox());
      assert(norm_Mg == norm_Mg);
      assert(norm_Mg > 0.);
#endif

      std::shared_ptr<pdat::CellData<double> > exact_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_exact_id)));
      std::shared_ptr<pdat::CellData<double> > rhs_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_rhs_id)));
      TBOX_ASSERT(exact_data);
      TBOX_ASSERT(rhs_data);

      /*
       * Set source function and exact solution.
       */
#if (NDIM == 2)
      SAMRAI_F77_FUNC(phasesetexactandrhs2d, PHASESETEXACTANDRHS2D)
      (pbox.lower()[0], pbox.upper()[0], pbox.lower()[1], pbox.upper()[1],
       exact_data->getPointer(), rhs_data->getPointer(), patch_geom->getDx(),
       patch_geom->getXLower(), d_mobility, d_delta, d_omega, d_gamma);
#endif
#if (NDIM == 3)
      SAMRAI_F77_FUNC(phasesetexactandrhs3d, PHASESETEXACTANDRHS3D)
      (pbox.lower()[0], pbox.upper()[0], pbox.lower()[1], pbox.upper()[1],
       pbox.lower()[2], pbox.upper()[2], exact_data->getPointer(),
       rhs_data->getPointer(), patch_geom->getDx(), patch_geom->getXLower(),
       d_mobility, d_delta, d_omega, d_gamma);
#endif

   }  // End patch loop.

#ifdef DEBUG_CHECK_ASSERTIONS
   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
   double value = cellops.maxNorm(d_rhs_id);
   assert(value == value);
#endif
}

/*
 *************************************************************************
 * Reset the hierarchy-dependent internal information.
 *************************************************************************
 */
void PhaseHypre::resetHierarchyConfiguration(
    const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
    int coarsest_level, int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);

   d_hierarchy = new_hierarchy;
}

void PhaseHypre::setCtoZero()
{
   math::HierarchyCellDataOpsReal<double> mathops(d_hierarchy);

   mathops.setToScalar(d_c_id, 0., false);

   d_poisson_spec.setCPatchDataId(d_c_id);
   d_C_is_set = true;
}

// C = 1 + gamma * phi_mobility * (
//       phi_well_scale * phi_well_func'' +
//       eta_well_scale * eta_well_func * phi_interp_func''
void PhaseHypre::setC(const int phi_id, const double gamma,
                      const EnergyInterpolationType phi_interp_func_type,
                      const double phi_well_scale,
                      const std::string phi_well_func_type)
{
   assert(phi_id >= 0);
   assert(d_m_id >= 0);
   assert(d_c_id >= 0);
   assert(d_M_is_set);

   // tbox::pout<<"PhaseHypre::setC()..."<<endl;
   const char interpf = energyInterpChar(phi_interp_func_type);

   int ln = 0;
   {
      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::iterator pi(level->begin()); pi != level->end();
           ++pi) {

         std::shared_ptr<hier::Patch> patch = *pi;
         const hier::Box& patch_box = patch->getBox();

         std::shared_ptr<pdat::CellData<double> > phi_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phi_id)));

         std::shared_ptr<pdat::CellData<double> > mdata(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_m_id)));

         std::shared_ptr<pdat::CellData<double> > cdata(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_c_id)));

         setCOnPatchPrivate(phi_data, mdata, cdata, gamma, &interpf,
                            phi_well_scale, phi_well_func_type.c_str(),
                            patch_box);
      }
   }

   d_poisson_spec.setCPatchDataId(d_c_id);
   d_C_is_set = true;
}
void PhaseHypre::setCOnPatchPrivate(
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_m,
    std::shared_ptr<pdat::CellData<double> > cd_c, const double gamma,
    const char* phi_interp_func_type, const double phi_well_scale,
    const char* phi_well_func_type, const hier::Box& pbox)
{
   double* ptr_phi = cd_phi->getPointer();
   double* ptr_m = cd_m->getPointer();
   double* ptr_c = cd_c->getPointer();

   const hier::Box& c_gbox = cd_c->getGhostBox();
   int imin_c = c_gbox.lower(0);
   int jmin_c = c_gbox.lower(1);
   int jp_c = c_gbox.numberCells(0);
   int kmin_c = 0;
   int kp_c = 0;
#if (NDIM == 3)
   kmin_c = c_gbox.lower(2);
   kp_c = jp_c * c_gbox.numberCells(1);
#endif

   const hier::Box& m_gbox = cd_m->getGhostBox();
   int imin_m = m_gbox.lower(0);
   int jmin_m = m_gbox.lower(1);
   int jp_m = m_gbox.numberCells(0);
   int kmin_m = 0;
   int kp_m = 0;
#if (NDIM == 3)
   kmin_m = m_gbox.lower(2);
   kp_m = jp_m * m_gbox.numberCells(1);
#endif

   const hier::Box& pf_gbox = cd_phi->getGhostBox();
   int imin_pf = pf_gbox.lower(0);
   int jmin_pf = pf_gbox.lower(1);
   int jp_pf = pf_gbox.numberCells(0);
   int kmin_pf = 0;
   int kp_pf = 0;
#if (NDIM == 3)
   kmin_pf = pf_gbox.lower(2);
   kp_pf = jp_pf * pf_gbox.numberCells(1);
#endif

   int imin = pbox.lower(0);
   int imax = pbox.upper(0);
   int jmin = pbox.lower(1);
   int jmax = pbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = pbox.lower(2);
   kmax = pbox.upper(2);
#endif

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_c =
                (ii - imin_c) + (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;

            const int idx_m =
                (ii - imin_m) + (jj - jmin_m) * jp_m + (kk - kmin_m) * kp_m;

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;

            const double m = ptr_m[idx_m];
            const double phi = ptr_phi[idx_pf];

            const double g_phi_dbl_prime =
                SECOND_DERIV_WELL_FUNC(phi, phi_well_func_type);

            const double gamma_m = gamma * m;

            // C = 1 + gamma * phi_mobility * (
            //       phi_well_scale * phi_well_func'' +
            //       eta_well_scale * eta_well_func * phi_interp_func''

            ptr_c[idx_c] = 1.0 + gamma_m * phi_well_scale * g_phi_dbl_prime;
         }
      }
   }
}

/*
 *************************************************************************
 * Set up the initial guess and problem parameters
 * and solve the Poisson problem.  We explicitly initialize and
 * deallocate the solver state in this example.
 *************************************************************************
 */
int PhaseHypre::solve(EnergyInterpolationType phase_interp_func_type,
                      double phase_well_scale, std::string phase_well_func_type)
{
   assert(d_m_id >= 0);
   assert(d_comp_soln_id >= 0);

   if (!d_hierarchy) {
      TBOX_ERROR(d_object_name << "Cannot solve using an uninitialized "
                                  "object.\n");
   }

   /*
    * Fill in the initial guess.
    */
   for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;
         std::shared_ptr<pdat::CellData<double> > data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_comp_soln_id)));
         TBOX_ASSERT(data);
         data->fill(0.0);
      }
   }

   d_poisson_solver.initializeSolverState(d_hierarchy, 0);

   /*
    * Set the parameters for the equation
    *    -M*div(D grad(u)) + C*u = f
    */
   d_poisson_spec.setMPatchDataId(d_m_id);
   d_M_is_set = true;

#if 1
   setC(d_exact_id, d_gamma, phase_interp_func_type, d_omega,
        phase_well_func_type);
#else
   setCtoZero();
#endif

   d_poisson_spec.setDConstant(-d_gamma * d_epsilon * d_epsilon);

#ifdef DEBUG_CHECK_ASSERTIONS
   math::PatchCellDataNormOpsReal<double> ops;
   for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;
         std::shared_ptr<pdat::CellData<double> > data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_m_id)));
         TBOX_ASSERT(data);
         const double norm_Mg = ops.maxNorm(data, data->getGhostBox());
         assert(norm_Mg == norm_Mg);
      }
   }

   math::HierarchyCellDataOpsReal<double> cellops(d_hierarchy);
   double value = cellops.L1Norm(d_rhs_id);
   assert(value == value);
#endif


   d_poisson_solver.setPhysicalBcCoefObject(d_physical_bc_coef.get());
   d_poisson_solver.setMatrixCoefficients(d_poisson_spec);

   tbox::plog << "solving..." << std::endl;
   d_poisson_solver.setStoppingCriteria(1.e-8);
   int solver_ret = d_poisson_solver.solveSystem(d_comp_soln_id, d_rhs_id);
   /*
    * Present data on the solve.
    */
   tbox::plog << d_object_name << " Hypre solve " << (solver_ret ? "" : "NOT ")
              << "converged\n"
              << "\titerations: " << d_poisson_solver.getNumberOfIterations()
              << "\n"
              << "\tresidual: " << d_poisson_solver.getRelativeResidualNorm()
              << "\n";

   return !solver_ret;
}

#ifdef HAVE_HDF5
/*
 *************************************************************************
 * Set up external plotter to plot internal data from this class.
 * Register variables appropriate for plotting.
 *************************************************************************
 */
int PhaseHypre::setupPlotter(appu::VisItDataWriter& plotter) const
{
   if (!d_hierarchy) {
      TBOX_ERROR(d_object_name << ": No hierarchy in\n"
                               << " PhaseHypre::setupPlotter\n"
                               << "The hierarchy must be set before calling\n"
                               << "this function.\n");
   }
   plotter.registerPlotQuantity("Computed solution", "SCALAR", d_comp_soln_id);
   plotter.registerDerivedPlotQuantity("Error", "SCALAR",
                                       (appu::VisDerivedDataStrategy*)this);
   plotter.registerPlotQuantity("Exact solution", "SCALAR", d_exact_id);
   plotter.registerPlotQuantity("Poisson source", "SCALAR", d_rhs_id);

   return 0;
}
#endif

/*
 *************************************************************************
 * Write derived data to the given stream.
 *************************************************************************
 */
bool PhaseHypre::packDerivedDataIntoDoubleBuffer(
    double* buffer, const hier::Patch& patch, const hier::Box& region,
    const std::string& variable_name, int depth_id,
    double simulation_time) const
{
   NULL_USE(depth_id);
   NULL_USE(simulation_time);

   pdat::CellData<double>::iterator icell(pdat::CellGeometry::begin(region));
   pdat::CellData<double>::iterator icellend(pdat::CellGeometry::end(region));

   if (variable_name == "Error") {
      std::shared_ptr<pdat::CellData<double> > current_solution_(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_comp_soln_id)));
      std::shared_ptr<pdat::CellData<double> > exact_solution_(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_exact_id)));
      TBOX_ASSERT(current_solution_);
      TBOX_ASSERT(exact_solution_);
      pdat::CellData<double>& current_solution = *current_solution_;
      pdat::CellData<double>& exact_solution = *exact_solution_;
      for (; icell != icellend; ++icell) {
         double diff = (current_solution(*icell) - exact_solution(*icell));
         *buffer = diff;
         buffer = buffer + 1;
      }
   } else {
      // Did not register this name.
      TBOX_ERROR("Unregistered variable name '" << variable_name << "' in\n"
                                                << "FACPoissonX::"
                                                   "packDerivedDataIntoDoubleBu"
                                                   "ffer");
   }
   // Return true if this patch has derived data on it.
   // False otherwise.
   return true;
}

double PhaseHypre::compareSolutionWithExact()
{
   math::HierarchyCellDataOpsReal<double> mathops(d_hierarchy);

   mathops.subtract(d_exact_id, d_exact_id, d_comp_soln_id);
   return mathops.maxNorm(d_exact_id);
}
