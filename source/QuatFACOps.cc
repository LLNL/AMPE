// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
/*
 * This class provides operator specific functions supporting
 * a quaternion system FAC solver.  See the header file
 * QuatFACOps.h for additional documentation of the
 * class member functions and data.
 *
 * This file was adapted from solv::CellPoissonFACOps.C in the
 * SAMRAI library.
 */

#include "QuatFACOps.h"
#include "QuatFort.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cassert>

// Static class member definitions

std::shared_ptr<pdat::CellVariable<double> > QuatFACOps::s_cell_scratch_var;

std::shared_ptr<pdat::SideVariable<double> > QuatFACOps::s_flux_scratch_var;

std::shared_ptr<pdat::OutersideVariable<double> >
    QuatFACOps::s_oflux_scratch_var;

std::shared_ptr<pdat::SideVariable<double> > QuatFACOps::s_face_coef_var;

std::shared_ptr<pdat::SideVariable<double> > QuatFACOps::s_face_coef_deriv_var;

std::shared_ptr<pdat::CellVariable<double> > QuatFACOps::s_q_local_var;

std::shared_ptr<pdat::CellVariable<double> > QuatFACOps::s_residual_var;

std::shared_ptr<pdat::CellVariable<double> > QuatFACOps::s_sqrt_m_var;

std::shared_ptr<pdat::CellVariable<double> > QuatFACOps::s_m_deriv_var;


/*
********************************************************************
* Constructor.                                                     *
********************************************************************
*/

QuatFACOps::QuatFACOps(const int ql,
                       std::shared_ptr<QuatFaceCoeff> quat_face_coeff_strategy,
                       const std::string& object_name,
                       const std::shared_ptr<tbox::Database>& database)
    : d_qlen(ql),
      d_quat_face_coeff_strategy(quat_face_coeff_strategy),
      d_object_name(object_name),
      d_hierarchy(),
      d_ln_min(-1),
      d_ln_max(-1),
      d_cf_boundary(),
      d_cf_discretization("Ewing"),
      d_prolongation_method("CONSTANT_REFINE"),
      d_levelsolver_tolerance(1.e-8),
      d_levelsolver_max_iterations(10),
      d_coarse_levelsolver_tolerance(1.e-8),
      d_coarse_levelsolver_max_iterations(10),
      d_flux_id(-1),
      d_levelsolver_database(database && database->isDatabase("hypre_solver")
                                 ? database->getDatabase("hypre_solver")
                                 : std::shared_ptr<tbox::Database>()),
      d_physical_bc_coef(NULL),
      d_flux_scratch_id(-1),
      d_oflux_scratch_id(-1),
      d_bc_helper(tbox::Dimension(NDIM), d_object_name + "::bc helper"),
      d_enable_logging(false),
      d_verbose(false),
      d_solver(NULL),
      d_hopscell(),
      d_gamma(tbox::MathUtilities<double>::getSignalingNaN()),
      d_rotation_index_id(-1)
{
   // Since one can't initialize arrays in the initializer list, we do it here
   d_cell_scratch_id = -1;

   t_restrict_solution = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::restrictSolution()");
   t_restrict_residual = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::restrictResidual()");
   t_prolong = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::prolongErrorAndCorrect()");
   t_smooth_error = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::smoothError()");
   t_solve_coarsest = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::solveCoarsestLevel()");
   t_compute_composite_residual = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::computeCompositeResidualOnLevel()");
   t_compute_rhs = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::computeRHS()");
   t_compute_residual_norm = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatFACOps::computeResidualNorm()");

   if (NDIM == 1) {
      TBOX_ERROR(d_object_name << ": 1D is not implemented.\n");
   }

   if (!s_cell_scratch_var || !s_flux_scratch_var || !s_oflux_scratch_var ||
       !s_q_local_var || !s_residual_var || !s_sqrt_m_var || !s_m_deriv_var) {
      assert(!s_cell_scratch_var);
      assert(!s_flux_scratch_var);
      assert(!s_oflux_scratch_var);
      assert(!s_q_local_var);
      assert(!s_residual_var);
      assert(!s_sqrt_m_var);
      assert(!s_m_deriv_var);
      if (!s_cell_scratch_var) {
         s_cell_scratch_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            "QuatFACOps::private_cell_scratch",
                                            d_qlen));
      }

      if (!s_flux_scratch_var) {
         s_flux_scratch_var.reset(
             new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                            "QuatFACOps::private_flux_scratch",
                                            d_qlen));
      }

      if (!s_oflux_scratch_var) {
         s_oflux_scratch_var.reset(
             new pdat::OutersideVariable<double>(tbox::Dimension(NDIM),
                                                 "QuatFACOps::private_oflux_"
                                                 "scratch",
                                                 d_qlen));
      }

      if (!s_face_coef_var) {
         s_face_coef_var.reset(new pdat::SideVariable<double>(
             tbox::Dimension(NDIM), "QuatFACOps::private_face_coef", 1));
      }

      if (!s_face_coef_deriv_var) {
         s_face_coef_deriv_var.reset(
             new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                            "QuatFACOps::private_face_coef_"
                                            "deriv",
                                            2));
      }

      if (!s_q_local_var) {
         s_q_local_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), "QuatFACOps::private_q_local", d_qlen));
      }

      if (!s_residual_var) {
         s_residual_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), "QuatFACOps::private_residual", d_qlen));
      }

      if (!s_sqrt_m_var) {
         s_sqrt_m_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), "QuatFACOps::private_mobility_sqrt", 1));
      }

      if (!s_m_deriv_var) {
         s_m_deriv_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), "QuatFACOps::private_mobility_deriv", 1));
      }
   }

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
   std::shared_ptr<hier::VariableContext> private_context =
       vdb->getContext(object_name + "::PRIVATE_CONTEXT");
   std::shared_ptr<hier::VariableContext> scratch_context =
       vdb->getContext("SCRATCH");

   d_cell_scratch_id =
       vdb->registerVariableAndContext(s_cell_scratch_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       1));

   d_flux_scratch_id =
       vdb->registerVariableAndContext(s_flux_scratch_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));
   d_oflux_scratch_id =
       vdb->registerVariableAndContext(s_oflux_scratch_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));

   d_face_coef_scratch_id =
       vdb->registerVariableAndContext(s_face_coef_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));

   d_face_coef_id =
       vdb->registerVariableAndContext(s_face_coef_var, scratch_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));

   d_face_coef_deriv_id =
       vdb->registerVariableAndContext(s_face_coef_deriv_var, scratch_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));

   d_q_local_id =
       vdb->registerVariableAndContext(s_q_local_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       1));

   d_residual_id =
       vdb->registerVariableAndContext(s_residual_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));

   d_sqrt_m_id =
       vdb->registerVariableAndContext(s_sqrt_m_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       1));

   d_m_deriv_id =
       vdb->registerVariableAndContext(s_m_deriv_var, private_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));

   /*
    * Some variables initialized by default are overriden by input.
    */
   if (database) {

      d_levelsolver_tolerance =
          database->getDoubleWithDefault("levelsolver_tolerance",
                                         d_levelsolver_tolerance);
      d_levelsolver_max_iterations =
          database->getIntegerWithDefault("levelsolver_max_iterations",
                                          d_levelsolver_max_iterations);

      d_coarse_levelsolver_tolerance =
          database->getDoubleWithDefault("coarse_levelsolver_tolerance",
                                         d_coarse_levelsolver_tolerance);
      d_coarse_levelsolver_max_iterations =
          database->getIntegerWithDefault("coarse_levelsolver_max_iterations",
                                          d_coarse_levelsolver_max_iterations);

      d_cf_discretization = database->getStringWithDefault("cf_discretization",
                                                           d_cf_discretization);

      d_prolongation_method =
          database->getStringWithDefault("prolongation_method",
                                         d_prolongation_method);

      d_enable_logging =
          database->getBoolWithDefault("enable_logging", d_enable_logging);
   }

   /*
    * Check validity of the flux patch data id
    */
   checkFluxPatchDataIndex();
}


QuatFACOps::~QuatFACOps(void)
{
   for (int ln = 0; ln < (int)d_quat_level_solver.size(); ln++) {
      if (d_quat_level_solver[ln] != NULL) delete d_quat_level_solver[ln];
   }
}


/*
************************************************************************
* FACOperatorStrategy virtual initializeOperatorState function.   *
*                                                                      *
* Set internal variables to correspond to the solution passed in.      *
* Look up transfer operators.                                          *
************************************************************************
*/

void QuatFACOps::initializeOperatorState(
    const solv::SAMRAIVectorReal<double>& solution,
    const solv::SAMRAIVectorReal<double>& rhs)
{
   deallocateOperatorState();
   d_hierarchy = solution.getPatchHierarchy();
   d_ln_min = solution.getCoarsestLevelNumber();
   d_ln_max = solution.getFinestLevelNumber();
   d_hopscell.reset(new math::HierarchyCellDataOpsReal<double>(d_hierarchy,
                                                               d_ln_min,
                                                               d_ln_max));
   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();
#ifdef DEBUG_CHECK_ASSERTIONS

   if (d_physical_bc_coef == NULL) {
      /*
       * It's an error not to have bc object set.
       * Note that the bc object cannot be passed in through
       * the argument because the interface is inherited.
       */
      TBOX_ERROR(d_object_name << ": No physical bc object in\n"
                               << "QuatFACOps::initializeOperatorState\n"
                               << "You must use "
                               << "QuatFACOps::setPhysicalBcCoefObject\n"
                               << "to set one before calling "
                                  "initializeOperatorState\n");
   }

   if (solution.getNumberOfComponents() != 1) {
      TBOX_WARNING(d_object_name << ": Solution std::vector has invalid number "
                                    "of "
                                    "components.\n"
                                 << "Solver is for 1 components only.\n");
   }
   if (rhs.getNumberOfComponents() != 1) {
      TBOX_WARNING(d_object_name << ": RHS std::vector has invalid number of "
                                    "components.\n"
                                 << "Solver is for 1 components only.\n");
   }

   /*
    * Make sure that solution and rhs data
    *   are of correct type
    *   are allocated
    *   has sufficient ghost width
    */
   std::shared_ptr<hier::Variable> var;
   {
      vdb->mapIndexToVariable(rhs.getComponentDescriptorIndex(0), var);
      if (!var) {
         TBOX_ERROR(d_object_name << ": RHS component does not\n"
                                  << "correspond to a variable.\n");
      }
      std::shared_ptr<pdat::CellVariable<double> > cell_var(
          SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(
              var));
      if (!cell_var) {
         TBOX_ERROR(d_object_name << ": RHS component variable is not "
                                     "cell-centered double\n");
      }
   }
   {
      vdb->mapIndexToVariable(solution.getComponentDescriptorIndex(0), var);
      if (!var) {
         TBOX_ERROR(d_object_name << ": Solution component does not\n"
                                  << "correspond to a variable.\n");
      }
      std::shared_ptr<pdat::CellVariable<double> > cell_var(
          SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(
              var));
      if (!cell_var) {
         TBOX_ERROR(d_object_name << ": Solution component variable is not "
                                     "cell-centered double\n");
      }
   }
   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
           pi++) {
         hier::Patch& patch = **pi;
         {
            std::shared_ptr<hier::PatchData> fd =
                patch.getPatchData(rhs.getComponentDescriptorIndex(0));
            if (fd) {
               /*
               Some data checks can only be done if the data already exists.
             */
               std::shared_ptr<pdat::CellData<double> > cd(
                   SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                                          hier::PatchData>(fd));
               if (!cd) {
                  TBOX_ERROR(d_object_name << ": RHS data component is not "
                                              "cell-centered double\n");
               }
               if (cd->getDepth() > d_qlen) {
                  TBOX_WARNING(d_object_name << ": Depth of RHS data component "
                                                "is greater than "
                                             << d_qlen << "\n");
               }
            }
            std::shared_ptr<hier::PatchData> ud =
                patch.getPatchData(solution.getComponentDescriptorIndex(0));
            if (ud) {
               /*
               Some data checks can only be done if the data already exists.
             */
               std::shared_ptr<pdat::CellData<double> > cd(
                   SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                                          hier::PatchData>(ud));
               if (!cd) {
                  TBOX_ERROR(d_object_name << ": Solution data component is "
                                              "not cell-centered double\n");
               }
               if (cd->getDepth() > d_qlen) {
                  TBOX_WARNING(d_object_name << ": Depth of solution data "
                                                "component is greater than "
                                             << d_qlen << "\n");
               }
               if (cd->getGhostCellWidth() <
                   hier::IntVector(tbox::Dimension(NDIM), 1)) {
                  TBOX_ERROR(d_object_name << ": Solution component 0 data has "
                                              "insufficient ghost width\n");
               }
            }
         }
      }
   }

   /*
    * Solution and rhs must have some similar properties.
    */
   if (rhs.getPatchHierarchy() != d_hierarchy ||
       rhs.getCoarsestLevelNumber() != d_ln_min ||
       rhs.getFinestLevelNumber() != d_ln_max) {
      TBOX_ERROR(d_object_name << ": solution and rhs do not have\n"
                               << "the same set of patch levels.\n");
   }

#endif

   /*
    * Initialize the coarse-fine boundary description for the
    * hierarchy.
    */
   d_cf_boundary.resize(d_hierarchy->getNumberOfLevels());

   hier::IntVector max_gcw(tbox::Dimension(NDIM), 1);
   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      d_cf_boundary[ln].reset(
          new hier::CoarseFineBoundary(*d_hierarchy, ln, max_gcw));
   }

   d_quat_level_solver.resize(d_hierarchy->getNumberOfLevels());

   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {

      // Create new solver on this level
      d_quat_level_solver[ln - d_ln_min] =
          new QuatLevelSolver(d_qlen,
                              d_object_name + "::quaternion_level_solver-" +
                                  tbox::Utilities::intToString(ln, 2),
                              d_levelsolver_database);

      d_quat_level_solver[ln - d_ln_min]->setVerbose(d_verbose);

      d_quat_level_solver[ln - d_ln_min]->initializeSolverState(d_hierarchy,
                                                                ln);

      /*
       * Share the boundary condition object with the hypre solver
       * to make sure that boundary condition settings are consistent
       * between the two objects.
       */
      d_quat_level_solver[ln - d_ln_min]->setPhysicalBcCoefObject(
          d_physical_bc_coef);
   }

   /*
    * Get the transfer operators.
    * Flux coarsening is conservative.
    * Cell (solution, error, etc) coarsening is conservative.
    * Cell refinement from same level is constant refinement.
    * Cell refinement from coarser level is chosen by the
    *   choice of coarse-fine discretization, d_cf_discretization,
    *   which should be set to either "Ewing" or one of the
    *   acceptable std::strings for looking up the refine operator.
    */
   std::shared_ptr<geom::CartesianGridGeometry> geometry(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry,
                              hier::BaseGridGeometry>(
           d_hierarchy->getGridGeometry()));
   TBOX_ASSERT(geometry);

   std::shared_ptr<hier::Variable> variable;

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_prolongation_refine_operator =
       geometry->lookupRefineOperator(variable, d_prolongation_method);
   d_urestriction_coarsen_operator = d_rrestriction_coarsen_operator =
       geometry->lookupCoarsenOperator(variable, "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable(d_oflux_scratch_id, variable);
   d_flux_coarsen_operator =
       geometry->lookupCoarsenOperator(variable, "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_ghostfill_refine_operator =
       geometry->lookupRefineOperator(variable, d_cf_discretization == "Ewing"
                                                    ? "CONSTANT_REFINE"
                                                    : d_cf_discretization);

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_ghostfill_nocoarse_refine_operator =
       geometry->lookupRefineOperator(variable, "CONSTANT_REFINE");

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!d_prolongation_refine_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find prolongation refine operator");
   }
   if (!d_urestriction_coarsen_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find restriction coarsening "
                                  "operator");
   }
   if (!d_rrestriction_coarsen_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find restriction coarsening "
                                  "operator");
   }
   if (!d_flux_coarsen_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find flux coarsening operator");
   }
   if (!d_ghostfill_refine_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find ghost filling refinement "
                                  "operator");
   }
   if (!d_ghostfill_nocoarse_refine_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find ghost filling refinement "
                                  "operator");
   }
#endif

   for (int ln = d_ln_min + 1; ln <= d_ln_max; ++ln) {
      d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_oflux_scratch_id);
   }

   /*
     Make space for saving communication schedules.
     There is no need to delete the old schedules first
     because we have deallocated the solver state above.
   */
   d_prolongation_refine_schedules.resize(d_ln_max + 1);
   d_urestriction_coarsen_schedules.resize(d_ln_max + 1);
   d_rrestriction_coarsen_schedules.resize(d_ln_max + 1);
   d_ghostfill_refine_schedules.resize(d_ln_max + 1);
   d_ghostfill_nocoarse_refine_schedules.resize(d_ln_max + 1);
   d_flux_coarsen_schedules.resize(d_ln_max + 1);

   d_prolongation_refine_algorithm.reset(new xfer::RefineAlgorithm());
   d_urestriction_coarsen_algorithm.reset(
       new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
   d_rrestriction_coarsen_algorithm.reset(
       new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
   d_flux_coarsen_algorithm.reset(
       new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
   d_ghostfill_refine_algorithm.reset(new xfer::RefineAlgorithm());
   d_ghostfill_nocoarse_refine_algorithm.reset(new xfer::RefineAlgorithm());


   /*
     Register the transfer operations
   */
   assert(d_cell_scratch_id >= 0);
   d_prolongation_refine_algorithm->registerRefine(
       d_cell_scratch_id, solution.getComponentDescriptorIndex(0),
       d_cell_scratch_id, d_prolongation_refine_operator);
   d_urestriction_coarsen_algorithm->registerCoarsen(
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       d_urestriction_coarsen_operator);
   d_rrestriction_coarsen_algorithm->registerCoarsen(
       rhs.getComponentDescriptorIndex(0), rhs.getComponentDescriptorIndex(0),
       d_rrestriction_coarsen_operator);

   d_ghostfill_refine_algorithm->registerRefine(
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0), d_ghostfill_refine_operator);
   d_flux_coarsen_algorithm->registerCoarsen(
       ((d_flux_id != -1) ? d_flux_id : d_flux_scratch_id), d_oflux_scratch_id,
       d_flux_coarsen_operator);

   assert(d_ghostfill_nocoarse_refine_operator);
   d_ghostfill_nocoarse_refine_algorithm->registerRefine(
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       d_ghostfill_nocoarse_refine_operator);

   /*
     Create the refine schedules
   */
   for (int dest_ln = d_ln_min + 1; dest_ln <= d_ln_max; ++dest_ln) {
      xfer::RefinePatchStrategy* strategy = &d_bc_helper;
      d_prolongation_refine_schedules[dest_ln] =
          d_prolongation_refine_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln),
              std::shared_ptr<hier::PatchLevel>(), dest_ln - 1, d_hierarchy,
              strategy);
      if (!d_prolongation_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                                     "prolongation!\n");
      }
      d_ghostfill_refine_schedules[dest_ln] =
          d_ghostfill_refine_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln), dest_ln - 1, d_hierarchy,
              &d_bc_helper);
      if (!d_ghostfill_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                                     "ghost filling!\n");
      }
      d_ghostfill_nocoarse_refine_schedules[dest_ln] =
          d_ghostfill_nocoarse_refine_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln), &d_bc_helper);
      if (!d_ghostfill_nocoarse_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                                     "ghost filling on bottom level!\n");
      }
   }

   /*
     Create the coarsen schedules
   */

   for (int dest_ln = d_ln_min; dest_ln < d_ln_max; ++dest_ln) {
      d_urestriction_coarsen_schedules[dest_ln] =
          d_urestriction_coarsen_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln),
              d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_urestriction_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for U "
                                     "restriction!\n");
      }
      d_rrestriction_coarsen_schedules[dest_ln] =
          d_rrestriction_coarsen_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln),
              d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_rrestriction_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for R "
                                     "restriction!\n");
      }
      d_flux_coarsen_schedules[dest_ln] =
          d_flux_coarsen_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln),
              d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_flux_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for "
                                     "flux transfer!\n");
      }
   }
   d_ghostfill_nocoarse_refine_schedules[d_ln_min] =
       d_ghostfill_nocoarse_refine_algorithm->createSchedule(
           d_hierarchy->getPatchLevel(d_ln_min), &d_bc_helper);
   if (!d_ghostfill_nocoarse_refine_schedules[d_ln_min]) {
      TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for ghost "
                                  "filling on bottom level!\n");
   }

   assert(d_face_coef_id);
   assert(d_face_coef_deriv_id);
   assert(d_face_coef_scratch_id);
   assert(d_q_local_id);
   assert(d_residual_id);
   assert(d_sqrt_m_id);
   assert(d_m_deriv_id);
   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level_ptr =
          d_hierarchy->getPatchLevel(ln);

      // Allocate the local patch data
      if (!level_ptr->checkAllocated(d_face_coef_id))
         level_ptr->allocatePatchData(d_face_coef_id);
      if (!level_ptr->checkAllocated(d_face_coef_deriv_id))
         level_ptr->allocatePatchData(d_face_coef_deriv_id);
      if (!level_ptr->checkAllocated(d_face_coef_scratch_id))
         level_ptr->allocatePatchData(d_face_coef_scratch_id);
      if (!level_ptr->checkAllocated(d_q_local_id))
         level_ptr->allocatePatchData(d_q_local_id);
      if (!level_ptr->checkAllocated(d_residual_id))
         level_ptr->allocatePatchData(d_residual_id);
      if (!level_ptr->checkAllocated(d_sqrt_m_id))
         level_ptr->allocatePatchData(d_sqrt_m_id);
      if (!level_ptr->checkAllocated(d_m_deriv_id))
         level_ptr->allocatePatchData(d_m_deriv_id);
   }
}


void QuatFACOps::computeDQuatDPhiFaceCoefs(const int dprime_id,
                                           const int phi_id,
                                           const int face_coef_id)
{
   // Note: This function assumes that the ghost cells of phi have already been
   // filled

   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
           pi++) {
         std::shared_ptr<hier::Patch> patch = *pi;

         std::shared_ptr<pdat::SideData<double> > dprime_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(dprime_id)));
         std::shared_ptr<pdat::CellData<double> > phi_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phi_id)));
         std::shared_ptr<pdat::SideData<double> > face_coef_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(face_coef_id)));

         computeDQuatDPhiFaceCoefsOnPatch(*patch, *dprime_data, *phi_data,
                                          *face_coef_data);
      }
   }
}


void QuatFACOps::takeSquareRootOnPatch(pdat::CellData<double>& data)
{
   const hier::Box& gbox = data.getGhostBox();
   const hier::Index& glower = gbox.lower();
   const hier::Index& gupper = gbox.upper();

#if NDIM == 2
   TAKE_SQUARE_ROOT2D(glower[0], gupper[0], glower[1], gupper[1],
                      data.getPointer(), glower[0], gupper[0], glower[1],
                      gupper[1]);
#endif
#if NDIM == 3
   TAKE_SQUARE_ROOT3D(glower[0], gupper[0], glower[1], gupper[1], glower[2],
                      gupper[2], data.getPointer(), glower[0], gupper[0],
                      glower[1], gupper[1], glower[2], gupper[2]);
#endif
}


void QuatFACOps::setOperatorCoefficients(const double gamma,
                                         const int mobility_id,
                                         const int mobility_deriv_id,
                                         const int phase_id, const int temp_id,
                                         const int face_coef_deriv_id,
                                         const int grad_q_id, const int q_id)
{
   d_gamma = gamma;

   // Check for non-positive mobility
   double mobility_min = d_hopscell->min(mobility_id);

   if (mobility_min <= 0.) {
      TBOX_ERROR(d_object_name << ": Non-positive mobility passed to "
                                  "setOperatorCoefficients().");
   }

   // Compute the face coefficients

   d_quat_face_coeff_strategy->computeFaceCoefs(d_hierarchy, phase_id, temp_id,
                                                grad_q_id, d_face_coef_id);

   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
           pi++) {
         std::shared_ptr<hier::Patch> patch = *pi;

         // Copy q solution to "local" array member

         std::shared_ptr<pdat::CellData<double> > q_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(q_id)));
         std::shared_ptr<pdat::CellData<double> > q_local_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_q_local_id)));

         q_local_data->copy(*q_data);

         // Copy mobility into sqrt_m_data (including ghost values assumed to be
         // filled)
         std::shared_ptr<pdat::CellData<double> > mobility_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(mobility_id)));
         std::shared_ptr<pdat::CellData<double> > sqrt_m_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_sqrt_m_id)));

         sqrt_m_data->copy(*mobility_data);
         takeSquareRootOnPatch(*sqrt_m_data);

         if (mobility_deriv_id >= 0) {

            // Copy mobility derivatives to local array

            std::shared_ptr<pdat::CellData<double> > mobility_deriv_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(mobility_deriv_id)));
            std::shared_ptr<pdat::CellData<double> > m_deriv_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_m_deriv_id)));

            m_deriv_data->copy(*mobility_deriv_data);
         }

         if (face_coef_deriv_id >= 0) {

            // Copy face coef derivatives to local array

            std::shared_ptr<pdat::SideData<double> > face_coef_deriv_data(
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(face_coef_deriv_id)));
            std::shared_ptr<pdat::SideData<double> > local_face_coef_deriv_data(
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(d_face_coef_deriv_id)));

            local_face_coef_deriv_data->copy(*face_coef_deriv_data);
         }
      }

      // Set the matrix coefficients
      d_quat_level_solver[ln - d_ln_min]->setMatrixCoefficients(d_gamma,
                                                                d_sqrt_m_id,
                                                                d_face_coef_id);
   }
}


/*
********************************************************************
* FACOperatorStrategy virtual deallocateOperatorState         *
* function.  Deallocate internal hierarchy-dependent data.         *
* State is allocated iff hierarchy is set.                         *
********************************************************************
*/
void QuatFACOps::deallocateOperatorState()
{
   if (d_hierarchy) {
      int ln;
      for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
         std::shared_ptr<hier::PatchLevel> level_ptr =
             d_hierarchy->getPatchLevel(ln);
         if (ln > d_ln_min) {
            level_ptr->deallocatePatchData(d_oflux_scratch_id);
         }
         level_ptr->deallocatePatchData(d_m_deriv_id);
         level_ptr->deallocatePatchData(d_sqrt_m_id);
         level_ptr->deallocatePatchData(d_residual_id);
         level_ptr->deallocatePatchData(d_q_local_id);
         level_ptr->deallocatePatchData(d_face_coef_scratch_id);
         level_ptr->deallocatePatchData(d_face_coef_deriv_id);
         level_ptr->deallocatePatchData(d_face_coef_id);
      }
      d_cf_boundary.resize(0);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_quat_level_solver.size() > 0);
#endif
      for (int i = 0; i < (int)d_quat_level_solver.size(); i++) {
         d_quat_level_solver[i]->deallocateSolverState();
         delete d_quat_level_solver[i];
      }
      d_quat_level_solver.clear();
      d_hierarchy.reset();
      d_ln_min = -1;
      d_ln_max = -1;

      d_prolongation_refine_algorithm.reset();
      d_urestriction_coarsen_algorithm.reset();
      d_rrestriction_coarsen_algorithm.reset();

      d_flux_coarsen_algorithm.reset();
      d_ghostfill_refine_algorithm.reset();
      d_ghostfill_nocoarse_refine_algorithm.reset();
   }
}


/*
********************************************************************
* FACOperatorStrategy virtual postprocessOneCycle function.  *
********************************************************************
*/
void QuatFACOps::postprocessOneCycle(
    int fac_cycle_num, const solv::SAMRAIVectorReal<double>& current_soln,
    const solv::SAMRAIVectorReal<double>& residual)
{
   (void)current_soln;
   (void)residual;

   if (d_enable_logging) {
      if (d_solver) {
         /*
          * Output convergence progress.  This is probably only appropriate
          * if the solver is NOT being used as a preconditioner.
          */
         double avg_factor, final_factor;
         d_solver->getConvergenceFactors(avg_factor, final_factor);
         tbox::plog << "iter=" << std::setw(4) << fac_cycle_num
                    << " resid=" << d_solver->getResidualNorm()
                    << " net conv=" << d_solver->getNetConvergenceFactor()
                    << " final conv=" << d_solver->getNetConvergenceFactor()
                    << " avg conv=" << d_solver->getAvgConvergenceFactor()
                    << std::endl;
      }
   }
}

/*
********************************************************************
* FACOperatorStrategy virtual restrictSolution function.      *
* After restricting solution, update ghost cells of the affected   *
* level.                                                           *
********************************************************************
*/
void QuatFACOps::restrictSolution(const solv::SAMRAIVectorReal<double>& s,
                                  solv::SAMRAIVectorReal<double>& d,
                                  int dest_ln)
{
   t_restrict_solution->start();

   xeqScheduleURestriction(d.getComponentDescriptorIndex(0),
                           s.getComponentDescriptorIndex(0), dest_ln);

   // Fill ghost cells.  Only component 0 has them.
   int id = d.getComponentDescriptorIndex(0);

   d_bc_helper.setHomogeneousBc(false);
   d_bc_helper.setTargetDataId(id);

   if (dest_ln == d_ln_min) {
      xeqScheduleGhostFillNoCoarse(id, dest_ln);
   } else {
      xeqScheduleGhostFill(id, dest_ln);
   }

   t_restrict_solution->stop();
}

/*
********************************************************************
* FACOperatorStrategy virtual restrictresidual function.      *
********************************************************************
*/
void QuatFACOps::restrictResidual(const solv::SAMRAIVectorReal<double>& s,
                                  solv::SAMRAIVectorReal<double>& d,
                                  int dest_ln)
{
   t_restrict_residual->start();

   xeqScheduleRRestriction(d.getComponentDescriptorIndex(0),
                           s.getComponentDescriptorIndex(0), dest_ln);

   t_restrict_residual->stop();
}


/*
***********************************************************************
* FACOperatorStrategy virtual prolongErrorAndCorrect function.   *
* After the prolongation, we set the physical boundary condition      *
* for the correction, which is zero.  Other ghost cell values,        *
* which are preset to zero, need not be set.                          *
***********************************************************************
*/
void QuatFACOps::prolongErrorAndCorrect(const solv::SAMRAIVectorReal<double>& s,
                                        solv::SAMRAIVectorReal<double>& d,
                                        int dest_ln)
{
   t_prolong->start();

   if (s.getPatchHierarchy() != d_hierarchy ||
       d.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                                  "internal state hierarchy.");
   }

   std::shared_ptr<hier::PatchLevel> fine_level(
       d_hierarchy->getPatchLevel(dest_ln));

   /*
    * Data is prolonged into the scratch space corresponding
    * to index d_cell_scratch_id and allocated here.
    */
   math::HierarchyCellDataOpsReal<double> hierarchy_math_ops(d_hierarchy,
                                                             dest_ln, dest_ln);

   fine_level->allocatePatchData(d_cell_scratch_id);

   /*
    * Refine solution into scratch space to fill the fine level
    * interior in the scratch space, then use that refined data
    * to correct the fine level error.
    */
   d_bc_helper.setTargetDataId(d_cell_scratch_id);
   d_bc_helper.setHomogeneousBc(true);

   const int src_index = s.getComponentDescriptorIndex(0);
   xeqScheduleProlongation(d_cell_scratch_id, src_index, d_cell_scratch_id,
                           dest_ln);

   /*
    * Add the refined error in the scratch space
    * to the error currently residing in the destination level.
    */
   const int dst_index = d.getComponentDescriptorIndex(0);
   hierarchy_math_ops.add(dst_index, dst_index, d_cell_scratch_id);

   fine_level->deallocatePatchData(d_cell_scratch_id);

   t_prolong->stop();
}


void QuatFACOps::smoothError(solv::SAMRAIVectorReal<double>& error,
                             const solv::SAMRAIVectorReal<double>& residual,
                             int ln, int num_sweeps)
{
   t_smooth_error->start();

   /*
    * Since we are calling a Hypre multigrid solver on each level,
    * rather an some simple iterative method like Gauss-Seidel,
    * we're not really smoothing here, but instead doing a full
    * solve of the residual equation.  Therefore, the num_sweeps
    * parameter is somewhat redundant.  There is no point to sweeping
    * more than once, since that one "sweep" actually corresponds
    * to multigrid solve performed with tolerances and iteration
    * limits set by setLevelSolverTolerance() and
    * setLevelSolverMaxIterations(), whose settings can be changed
    * in order to modify how hard the solver works on each level.
    * The main use of the num_sweep sweep parameter is therefore to
    * turn off either the pre or post smooth steps by setting the
    * respective parameter to 0, for which we test here in order to
    * determine whether or not to perform the level solve.
    */
   if (num_sweeps > 0) {
      doLevelSolve(error, residual, ln, d_levelsolver_max_iterations,
                   d_levelsolver_tolerance);
   }

   t_smooth_error->stop();
}


void QuatFACOps::doLevelSolve(solv::SAMRAIVectorReal<double>& solution,
                              const solv::SAMRAIVectorReal<double>& residual,
                              int ln, int max_iterations,
                              double residual_tolerance)
{
   assert(d_quat_level_solver[ln] != NULL);

   int q_solution_id = solution.getComponentDescriptorIndex(0);
   int q_residual_id = residual.getComponentDescriptorIndex(0);

   d_bc_helper.setTargetDataId(q_solution_id);
   d_bc_helper.setHomogeneousBc(true);

   // Fill q ghost cells
   xeqScheduleGhostFillNoCoarse(q_solution_id, ln);

   if (ln > d_ln_min) {
      /*
       * Perform a one-time transfer of data from coarser level,
       * to fill ghost boundaries that will not change through
       * the smoothing loop.
       */
      xeqScheduleGhostFill(q_solution_id, ln);
   }

   d_quat_level_solver[ln]->setStoppingCriteria(max_iterations,
                                                residual_tolerance);

   // Solve
   d_quat_level_solver[ln]->solveSystem(q_residual_id, q_solution_id);

   // Update q ghost cells
   xeqScheduleGhostFillNoCoarse(q_solution_id, ln);

   /*
    * Present data on the solve.
    * The QuatLevelSolver returns 0 if converged.
    */
   if (d_enable_logging)
      tbox::plog << d_object_name << " Quaterion solve at level " << ln
                 << "\titerations: "
                 << d_quat_level_solver[ln]->getNumberOfIterations() << "\n"
                 << "\tresidual: "
                 << d_quat_level_solver[ln]->getRelativeResidualNorm() << "\n";
}


/*
********************************************************************
* Fix flux on coarse-fine boundaries computed from a               *
* constant-refine interpolation of coarse level data.              *
********************************************************************
*/
void QuatFACOps::ewingFixFlux(const hier::Patch& patch,
                              const pdat::CellData<double>& soln_data,
                              const pdat::SideData<double>& face_coef_data,
                              pdat::SideData<double>& flux_data,
                              const hier::IntVector& ratio_to_coarser) const
{
   (void)patch;
   (void)soln_data;
   (void)face_coef_data;
   (void)flux_data;
   (void)ratio_to_coarser;

   TBOX_ERROR("QuatFACOps::ewingFixFlux) not implemented!!!\n");
#if 0
   const int patch_ln = patch.getPatchLevelNumber();
   const int pn = patch.getPatchNumber();
   std::shared_ptr< geom::CartesianPatchGeometry > patch_geom ( 
      SAMRAI_SHARED_PTR_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(patch.getPatchGeometry()) );
   const double * dx = patch_geom->getDx();
   const hier::Box & patch_box( patch.getBox() );
   const hier::Index & plower = patch_box.lower();
   const hier::Index & pupper = patch_box.upper();

   const std::vector< hier::BoundaryBox > & bboxes = d_cf_boundary[patch_ln].getBoundaries(pn,1);
   int nboxes=bboxes.getSize();

   for (int bn=0; bn<nboxes; ++bn ) {
     const hier::BoundaryBox &boundary_box=bboxes[bn];
#ifdef DEBUG_CHECK_ASSERTIONS
     assert( boundary_box.getBoundaryType() == 1 );
#endif
     const hier::Box &bdry_box = boundary_box.getBox();
     const hier::Index &blower = bdry_box.lower();
     const hier::Index &bupper = bdry_box.upper();
     const int location_index = boundary_box.getLocationIndex();
     const int depth = d_qlen;
#if NDIM == 2
     fixflux2D(flux_data.getPointer(0), flux_data.getPointer(1),
                flux_data.getGhostCellWidth()[0],
                flux_data.getGhostCellWidth()[1],
                depth,
                face_coef_data.getPointer(0), face_coef_data.getPointer(1),
                face_coef_data.getGhostCellWidth()[0],
                face_coef_data.getGhostCellWidth()[1],
                depth,
                soln_data.getPointer(),
                soln_data.getGhostCellWidth()[0],
                soln_data.getGhostCellWidth()[1],
                depth,
                plower[0], pupper[0], plower[1], pupper[1], depth,
                location_index,
                ratio_to_coarser,
                blower, bupper,
                dx);
#endif
#if NDIM == 3
     fixflux3D(flux_data.getPointer(0), flux_data.getPointer(1), flux_data.getPointer(2),
                flux_data.getGhostCellWidth()[0],
                flux_data.getGhostCellWidth()[1],
                flux_data.getGhostCellWidth()[2],
                depth,
                face_coef_data.getPointer(0), face_coef_data.getPointer(1), face_coef_data.getPointer(2),
                face_coef_data.getGhostCellWidth()[0],
                face_coef_data.getGhostCellWidth()[1],
                face_coef_data.getGhostCellWidth()[2],
                depth,
                soln_data.getPointer() ,
                soln_data.getGhostCellWidth()[0],
                soln_data.getGhostCellWidth()[1],
                soln_data.getGhostCellWidth()[2],
                depth,
                plower[0], pupper[0], plower[1], pupper[1], plower[2], pupper[2], depth,
                location_index,
                ratio_to_coarser,
                blower, bupper,
                dx);
#endif
   }
#endif
}


/*
********************************************************************
* FACOperatorStrategy virtual solveCoarsestLevel             *
* function                                                         *
********************************************************************
*/
int QuatFACOps::solveCoarsestLevel(
    solv::SAMRAIVectorReal<double>& data,
    const solv::SAMRAIVectorReal<double>& residual, int coarsest_ln)

{
   t_solve_coarsest->start();

   doLevelSolve(data, residual, coarsest_ln,
                d_coarse_levelsolver_max_iterations,
                d_coarse_levelsolver_tolerance);

   t_solve_coarsest->stop();

   return 0;  // Always successful
}


/*
********************************************************************
* FACOperatorStrategy virtual computeResidualNorm             *
* function                                                         *
********************************************************************
*/
double QuatFACOps::computeResidualNorm(
    const solv::SAMRAIVectorReal<double>& residual, int fine_ln, int coarse_ln)
{
   if (coarse_ln != residual.getCoarsestLevelNumber() ||
       fine_ln != residual.getFinestLevelNumber()) {
      TBOX_ERROR("QuatFACOps::computeResidualNorm() is not\n"
                 << "set up to compute residual except on the range of\n"
                 << "levels defining the std::vector.\n");
   }
   t_compute_residual_norm->start();
   /*
    * The residual std::vector was cloned from std::vectors that has
    * the proper weights associated with them, so we do not
    * have to explicitly weight the residuals.
    */

   double norm = residual.RMSNorm();

   t_compute_residual_norm->stop();

   return norm;
}


void QuatFACOps::computeLambdaOnPatch(
    const hier::Patch& patch, const pdat::SideData<double>& flux_data,
    const pdat::CellData<double>& q_data,
    std::shared_ptr<pdat::SideData<int> > rotation_index,
    pdat::CellData<double>& lambda_data) const
{
   assert(patch.inHierarchy());
   assert(q_data.getDepth() == d_qlen);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& q_gbox = q_data.getGhostBox();
   const hier::Index& qlower = q_gbox.lower();
   const hier::Index& qupper = q_gbox.upper();

   const hier::Box& f_gbox = flux_data.getGhostBox();
   const hier::Index& flower = f_gbox.lower();
   const hier::Index& fupper = f_gbox.upper();

   const hier::Box& l_gbox = lambda_data.getGhostBox();
   const hier::Index& llower = l_gbox.lower();
   const hier::Index& lupper = l_gbox.upper();

#if NDIM == 2
   if (rotation_index)
      COMPUTE_LAMBDA_FLUX2D_SYMM(lower[0], upper[0], lower[1], upper[1], d_qlen,
                                 flux_data.getPointer(0), flower[0],
                                 fupper[0] + 1, flower[1], fupper[1],
                                 flux_data.getPointer(1), flower[0], fupper[0],
                                 flower[1], fupper[1] + 1, q_data.getPointer(),
                                 qlower[0], qupper[0], qlower[1], qupper[1], dx,
                                 lambda_data.getPointer(0), llower[0],
                                 lupper[0], llower[1], lupper[1],
                                 rotation_index->getPointer(0),
                                 rotation_index->getPointer(1),
                                 rotation_index->getGhostCellWidth()[0]);
   else
      COMPUTE_LAMBDA_FLUX2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                            flux_data.getPointer(0), flower[0], fupper[0] + 1,
                            flower[1], fupper[1], flux_data.getPointer(1),
                            flower[0], fupper[0], flower[1], fupper[1] + 1,
                            q_data.getPointer(), qlower[0], qupper[0],
                            qlower[1], qupper[1], dx, lambda_data.getPointer(0),
                            llower[0], lupper[0], llower[1], lupper[1]);
#endif
#if NDIM == 3
   assert(!rotation_index);
   COMPUTE_LAMBDA_FLUX3D(lower[0], upper[0], lower[1], upper[1], lower[2],
                         upper[2], d_qlen, flux_data.getPointer(0), flower[0],
                         fupper[0] + 1, flower[1], fupper[1], flower[2],
                         fupper[2], flux_data.getPointer(1), flower[0],
                         fupper[0], flower[1], fupper[1] + 1, flower[2],
                         fupper[2], flux_data.getPointer(2), flower[0],
                         fupper[0], flower[1], fupper[1], flower[2],
                         fupper[2] + 1, q_data.getPointer(), qlower[0],
                         qupper[0], qlower[1], qupper[1], qlower[2], qupper[2],
                         dx, lambda_data.getPointer(0), llower[0], lupper[0],
                         llower[1], lupper[1], llower[2], lupper[2]);
#endif
}


/*
********************************************************************
* Check the validity of the flux variable id.                      *
********************************************************************
*/
void QuatFACOps::checkFluxPatchDataIndex() const
{
   if (d_flux_id != -1) {
      hier::VariableDatabase& vdb(*hier::VariableDatabase::getDatabase());
      std::shared_ptr<hier::Variable> var;
      vdb.mapIndexToVariable(d_flux_id, var);
      std::shared_ptr<pdat::SideVariable<double> > flux_var(
          SAMRAI_SHARED_PTR_CAST<pdat::SideVariable<double>, hier::Variable>(
              var));
      assert(flux_var);
   }
}


/*
*******************************************************************
*                                                                 *
* AMR-unaware patch-centered computational kernels.               *
*                                                                 *
*******************************************************************
*/

void QuatFACOps::computeDQuatDPhiFaceCoefsOnPatch(
    const hier::Patch& patch, pdat::SideData<double>& dprime_data,
    pdat::CellData<double>& phi_data,
    pdat::SideData<double>& face_coef_data) const
{
   assert(patch.inHierarchy());
   assert(dprime_data.getDepth() == 2);
   assert(face_coef_data.getDepth() == 1);

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& dp_gbox = dprime_data.getGhostBox();
   const hier::Index& dplower = dp_gbox.lower();
   const hier::Index& dpupper = dp_gbox.upper();

   const hier::Box& p_gbox = phi_data.getGhostBox();
   const hier::Index& plower = p_gbox.lower();
   const hier::Index& pupper = p_gbox.upper();

   const hier::Box& fc_gbox = face_coef_data.getGhostBox();
   const hier::Index& fclower = fc_gbox.lower();
   const hier::Index& fcupper = fc_gbox.upper();

#if NDIM == 2
   COMPUTE_DQUATDPHI_FACE_COEF2D(
       lower[0], upper[0], lower[1], upper[1], d_qlen,
       dprime_data.getPointer(0), dplower[0], dpupper[0] + 1, dplower[1],
       dpupper[1], dprime_data.getPointer(1), dplower[0], dpupper[0],
       dplower[1], dpupper[1] + 1, phi_data.getPointer(0), plower[0], pupper[0],
       plower[1], pupper[1], face_coef_data.getPointer(0), fclower[0],
       fcupper[0] + 1, fclower[1], fcupper[1], face_coef_data.getPointer(1),
       fclower[0], fcupper[0], fclower[1], fcupper[1] + 1);
#endif
#if NDIM == 3
   COMPUTE_DQUATDPHI_FACE_COEF3D(
       lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], d_qlen,
       dprime_data.getPointer(0), dplower[0], dpupper[0] + 1, dplower[1],
       dpupper[1], dplower[2], dpupper[2], dprime_data.getPointer(1),
       dplower[0], dpupper[0], dplower[1], dpupper[1] + 1, dplower[2],
       dpupper[2], dprime_data.getPointer(2), dplower[0], dpupper[0],
       dplower[1], dpupper[1], dplower[2], dpupper[2] + 1,
       phi_data.getPointer(0), plower[0], pupper[0], plower[1], pupper[1],
       plower[2], pupper[2], face_coef_data.getPointer(0), fclower[0],
       fcupper[0] + 1, fclower[1], fcupper[1], fclower[2], fcupper[2],
       face_coef_data.getPointer(1), fclower[0], fcupper[0], fclower[1],
       fcupper[1] + 1, fclower[2], fcupper[2], face_coef_data.getPointer(2),
       fclower[0], fcupper[0], fclower[1], fcupper[1], fclower[2],
       fcupper[2] + 1);
#endif
}


void QuatFACOps::computeFluxOnPatch(
    const hier::Patch& patch, const hier::IntVector& ratio_to_coarser_level,
    const pdat::SideData<double>& face_coef_data,
    const pdat::CellData<double>& q_data,
    pdat::SideData<double>& flux_data) const
{
   assert(patch.inHierarchy());
   assert(q_data.getDepth() == d_qlen);
   assert(flux_data.getDepth() == d_qlen);
   assert(face_coef_data.getDepth() == 1);

   // tbox::pout<<"QuatFACOps::computeFluxOnPatch() NOT using grad_q
   // data..."<<endl;

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& q_gbox = q_data.getGhostBox();
   const hier::Index& qlower = q_gbox.lower();
   const hier::Index& qupper = q_gbox.upper();

   const hier::Box& d_gbox = face_coef_data.getGhostBox();
   const hier::Index& dlower = d_gbox.lower();
   const hier::Index& dupper = d_gbox.upper();

   const hier::Box& f_gbox = flux_data.getGhostBox();
   const hier::Index& flower = f_gbox.lower();
   const hier::Index& fupper = f_gbox.upper();

#if NDIM == 2
   COMPUTE_FLUX2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                  face_coef_data.getPointer(0), dlower[0], dupper[0] + 1,
                  dlower[1], dupper[1], face_coef_data.getPointer(1), dlower[0],
                  dupper[0], dlower[1], dupper[1] + 1, q_data.getPointer(),
                  qlower[0], qupper[0], qlower[1], qupper[1], dx,
                  flux_data.getPointer(0), flower[0], fupper[0] + 1, flower[1],
                  fupper[1], flux_data.getPointer(1), flower[0], fupper[0],
                  flower[1], fupper[1] + 1);
#endif
#if NDIM == 3
   COMPUTE_FLUX3D(lower[0], upper[0], lower[1], upper[1], lower[2], upper[2],
                  d_qlen, face_coef_data.getPointer(0), dlower[0],
                  dupper[0] + 1, dlower[1], dupper[1], dlower[2], dupper[2],
                  face_coef_data.getPointer(1), dlower[0], dupper[0], dlower[1],
                  dupper[1] + 1, dlower[2], dupper[2],
                  face_coef_data.getPointer(2), dlower[0], dupper[0], dlower[1],
                  dupper[1], dlower[2], dupper[2] + 1, q_data.getPointer(),
                  qlower[0], qupper[0], qlower[1], qupper[1], qlower[2],
                  qupper[2], dx, flux_data.getPointer(0), flower[0],
                  fupper[0] + 1, flower[1], fupper[1], flower[2], fupper[2],
                  flux_data.getPointer(1), flower[0], fupper[0], flower[1],
                  fupper[1] + 1, flower[2], fupper[2], flux_data.getPointer(2),
                  flower[0], fupper[0], flower[1], fupper[1], flower[2],
                  fupper[2] + 1);
#endif

   const int patch_ln = patch.getPatchLevelNumber();

   if (d_cf_discretization == "Ewing" && patch_ln > d_ln_min) {
      ewingFixFlux(patch, q_data, face_coef_data, flux_data,
                   ratio_to_coarser_level);
   }
}

void QuatFACOps::computeFluxOnPatch(
    const hier::Patch& patch, const hier::IntVector& ratio_to_coarser_level,
    const pdat::SideData<double>& face_coef_data,
    const pdat::SideData<double>& gradq_data,
    pdat::SideData<double>& flux_data) const
{
   (void)ratio_to_coarser_level;

   // tbox::pout<<"check array sizes in
   // QuatFACOps::computeFluxOnPatch()..."<<endl;
   assert(patch.inHierarchy());
   assert(gradq_data.getDepth() == d_qlen * NDIM);
   assert(flux_data.getDepth() == d_qlen);
   assert(face_coef_data.getDepth() == 1);

   // tbox::pout<<"QuatFACOps::computeFluxOnPatch() using grad_q data..."<<endl;

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& d_gbox = face_coef_data.getGhostBox();
   const hier::Index& dlower = d_gbox.lower();
   const hier::Index& dupper = d_gbox.upper();

   const hier::Box& f_gbox = flux_data.getGhostBox();
   const hier::Index& flower = f_gbox.lower();
   const hier::Index& fupper = f_gbox.upper();

#if NDIM == 2
   COMPUTE_FLUX2D_FROM_GRADQ(
       lower[0], upper[0], lower[1], upper[1], d_qlen,
       face_coef_data.getPointer(0), dlower[0], dupper[0] + 1, dlower[1],
       dupper[1], face_coef_data.getPointer(1), dlower[0], dupper[0], dlower[1],
       dupper[1] + 1,
       gradq_data.getPointer(0, 0 * d_qlen),  // side 0, depth 0 (x component)
       gradq_data.getPointer(1, 1 * d_qlen),  // side 1, depth 1 (y component)
       flux_data.getPointer(0), flower[0], fupper[0] + 1, flower[1], fupper[1],
       flux_data.getPointer(1), flower[0], fupper[0], flower[1], fupper[1] + 1);
#endif
#if NDIM == 3
   COMPUTE_FLUX3D_FROM_GRADQ(
       lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], d_qlen,
       face_coef_data.getPointer(0), dlower[0], dupper[0] + 1, dlower[1],
       dupper[1], dlower[2], dupper[2], face_coef_data.getPointer(1), dlower[0],
       dupper[0], dlower[1], dupper[1] + 1, dlower[2], dupper[2],
       face_coef_data.getPointer(2), dlower[0], dupper[0], dlower[1], dupper[1],
       dlower[2], dupper[2] + 1,
       gradq_data.getPointer(0, 0 * d_qlen),  // side 0, depth 0 (x component)
       gradq_data.getPointer(1, 1 * d_qlen),  // side 1, depth 1 (y component)
       gradq_data.getPointer(2, 2 * d_qlen),  // side 2, depth 2 (z component)
       flux_data.getPointer(0), flower[0], fupper[0] + 1, flower[1], fupper[1],
       flower[2], fupper[2], flux_data.getPointer(1), flower[0], fupper[0],
       flower[1], fupper[1] + 1, flower[2], fupper[2], flux_data.getPointer(2),
       flower[0], fupper[0], flower[1], fupper[1], flower[2], fupper[2] + 1);
#endif
}


void QuatFACOps::computeSymmetricFluxOnPatch(
    const hier::Patch& patch, const hier::IntVector& ratio_to_coarser_level,
    const pdat::SideData<double>& face_coef_data,
    const pdat::CellData<double>& sqrt_m_data,
    const pdat::CellData<double>& q_data,
    pdat::SideData<double>& flux_data) const
{
   assert(patch.inHierarchy());
   assert(q_data.getDepth() == d_qlen);
   assert(flux_data.getDepth() == d_qlen);
   assert(face_coef_data.getDepth() == 1);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& q_gbox = q_data.getGhostBox();
   const hier::Index& qlower = q_gbox.lower();
   const hier::Index& qupper = q_gbox.upper();

   const hier::Box& m_gbox = sqrt_m_data.getGhostBox();
   const hier::Index& mlower = m_gbox.lower();
   const hier::Index& mupper = m_gbox.upper();

   const hier::Box& d_gbox = face_coef_data.getGhostBox();
   const hier::Index& dlower = d_gbox.lower();
   const hier::Index& dupper = d_gbox.upper();

   const hier::Box& f_gbox = flux_data.getGhostBox();
   const hier::Index& flower = f_gbox.lower();
   const hier::Index& fupper = f_gbox.upper();

#if NDIM == 2
   COMPUTE_SYM_FLUX2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                      face_coef_data.getPointer(0), dlower[0], dupper[0] + 1,
                      dlower[1], dupper[1], face_coef_data.getPointer(1),
                      dlower[0], dupper[0], dlower[1], dupper[1] + 1,
                      sqrt_m_data.getPointer(), mlower[0], mupper[0], mlower[1],
                      mupper[1], q_data.getPointer(), qlower[0], qupper[0],
                      qlower[1], qupper[1], dx, flux_data.getPointer(0),
                      flower[0], fupper[0] + 1, flower[1], fupper[1],
                      flux_data.getPointer(1), flower[0], fupper[0], flower[1],
                      fupper[1] + 1);
#endif
#if NDIM == 3
   COMPUTE_SYM_FLUX3D(
       lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], d_qlen,
       face_coef_data.getPointer(0), dlower[0], dupper[0] + 1, dlower[1],
       dupper[1], dlower[2], dupper[2], face_coef_data.getPointer(1), dlower[0],
       dupper[0], dlower[1], dupper[1] + 1, dlower[2], dupper[2],
       face_coef_data.getPointer(2), dlower[0], dupper[0], dlower[1], dupper[1],
       dlower[2], dupper[2] + 1, sqrt_m_data.getPointer(), mlower[0], mupper[0],
       mlower[1], mupper[1], mlower[2], mupper[2], q_data.getPointer(),
       qlower[0], qupper[0], qlower[1], qupper[1], qlower[2], qupper[2], dx,
       flux_data.getPointer(0), flower[0], fupper[0] + 1, flower[1], fupper[1],
       flower[2], fupper[2], flux_data.getPointer(1), flower[0], fupper[0],
       flower[1], fupper[1] + 1, flower[2], fupper[2], flux_data.getPointer(2),
       flower[0], fupper[0], flower[1], fupper[1], flower[2], fupper[2] + 1);
#endif

   const int patch_ln = patch.getPatchLevelNumber();

   if (d_cf_discretization == "Ewing" && patch_ln > d_ln_min) {
      ewingFixFlux(patch, q_data, face_coef_data, flux_data,
                   ratio_to_coarser_level);
   }
}


/*
********************************************************************
* FACOperatorStrategy virtual                                *
* computeCompositeResidualOnLevel function                         *
********************************************************************
*/

void QuatFACOps::computeCompositeResidualOnLevel(
    solv::SAMRAIVectorReal<double>& residual,
    const solv::SAMRAIVectorReal<double>& solution,
    const solv::SAMRAIVectorReal<double>& rhs, int ln,
    bool error_equation_indicator)
{
   t_compute_composite_residual->start();

   checkFluxPatchDataIndex();
#ifdef DEBUG_CHECK_ASSERTIONS
   if (residual.getPatchHierarchy() != d_hierarchy ||
       solution.getPatchHierarchy() != d_hierarchy ||
       rhs.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                                  "internal hierarchy.");
   }
#endif
   std::shared_ptr<hier::PatchLevel> level(d_hierarchy->getPatchLevel(ln));

   /*
    * Set up the bc helper so that when we use a refine schedule
    * to fill ghosts, the correct data is operated on.
    */
   const int q_id = solution.getComponentDescriptorIndex(0);
   d_bc_helper.setTargetDataId(q_id);
   d_bc_helper.setHomogeneousBc(error_equation_indicator);

   const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

   /*
    * Assumptions:
    * 1. Data does not yet exist in ghost boundaries.
    * 2. Residual data on next finer grid (if any)
    *    has been computed already.
    * 3. Flux data from next finer grid (if any) has
    *    been computed but has not been coarsened to
    *    this level.
    *
    * Steps:
    * S1. Fill solution ghost data by refinement
    *     or setting physical boundary conditions.
    *     This also brings in information from coarser
    *     to form the composite grid flux.
    * S2. Compute flux on ln.
    * S3. If next finer is available,
    *     Coarsen flux data on next finer level,
    *     overwriting flux computed from coarse data.
    * S4. Compute residual data from flux.
    */

   /* S1. Fill ghost data for quaternions (component 0). */
   {
      if (ln > d_ln_min) {
         /* Fill from current, next coarser level and physical boundary */
         xeqScheduleGhostFill(q_id, ln);
      } else {
         /* Fill from current and physical boundary */
         xeqScheduleGhostFillNoCoarse(q_id, ln);
      }
   }

   /*
    * For the whole level, make sure the internal
    * side-centered data is allocated and note
    * whether that data should be deallocated when done.
    * We do this for the whole level because the data
    * undergoes transfer operations which require the
    * whole level data.
    */
   bool deallocate_flux_data_when_done = false;
   if (flux_id == d_flux_scratch_id) {
      if (!level->checkAllocated(flux_id)) {
         level->allocatePatchData(flux_id);
         deallocate_flux_data_when_done = true;
      }
   }

   /*
    * S2. Compute flux on patches in level.
    */
   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        pi++) {
      std::shared_ptr<hier::Patch> patch = *pi;

      std::shared_ptr<pdat::CellData<double> > q_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              solution.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));
      std::shared_ptr<pdat::SideData<double> > face_coef_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(d_face_coef_id)));
      std::shared_ptr<pdat::CellData<double> > sqrt_m_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_sqrt_m_id)));

      computeSymmetricFluxOnPatch(*patch, level->getRatioToCoarserLevel(),
                                  *face_coef_data, *sqrt_m_data, *q_data,
                                  *flux_data);
   }

   /*
    * S3. Coarsen oflux data from next finer level so that
    * the computed flux becomes the composite grid flux.
    */
   if (ln < d_ln_max) {
      xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
   }

   /*
    * S4. Compute residual on patches in level.
    */
   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        pi++) {
      std::shared_ptr<hier::Patch> patch = *pi;
      std::shared_ptr<pdat::CellData<double> > q_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              solution.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::CellData<double> > q_rhs_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              rhs.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::CellData<double> > q_residual_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              residual.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::CellData<double> > sqrt_m_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_sqrt_m_id)));
      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));
      std::shared_ptr<pdat::SideData<int> > rotation_index;
      if (d_rotation_index_id >= 0) {
         rotation_index =
             std::dynamic_pointer_cast<pdat::SideData<int>, hier::PatchData>(
                 patch->getPatchData(d_rotation_index_id));
         assert(rotation_index);
      }

      computeResidualOnPatch(*patch, *flux_data, rotation_index, *sqrt_m_data,
                             *q_data, *q_rhs_data, *q_residual_data);

      if (ln > d_ln_min) {
         /*
          * Save outerflux data so that next coarser level
          *  can compute its coarse-fine composite flux.
          *  This is not strictly needed in this "compute residual"
          *  loop through the patches, but we put it here to
          *  avoid writing another loop for it.
          */
         std::shared_ptr<pdat::OutersideData<double> > oflux_data(
             SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>,
                                    hier::PatchData>(
                 patch->getPatchData(d_oflux_scratch_id)));
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(oflux_data);
#endif
         oflux_data->copy(*flux_data);
      }
   }

   if (deallocate_flux_data_when_done) {
      level->deallocatePatchData(flux_id);
   }

   t_compute_composite_residual->stop();
}


void QuatFACOps::computeResidualOnPatch(
    const hier::Patch& patch, const pdat::SideData<double>& flux_data,
    std::shared_ptr<pdat::SideData<int> > rotation_index,
    const pdat::CellData<double>& sqrt_m_data,
    const pdat::CellData<double>& q_soln_data,
    const pdat::CellData<double>& q_rhs_data,
    pdat::CellData<double>& q_residual_data) const
{
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));

   const double* dx = patch_geom->getDx();

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& m_gbox = sqrt_m_data.getGhostBox();
   const hier::Index& mlower = m_gbox.lower();
   const hier::Index& mupper = m_gbox.upper();

   const hier::Box& f_gbox = flux_data.getGhostBox();
   const hier::Index& flower = f_gbox.lower();
   const hier::Index& fupper = f_gbox.upper();

   const hier::Box& q_gbox = q_soln_data.getGhostBox();
   const hier::Index& qlower = q_gbox.lower();
   const hier::Index& qupper = q_gbox.upper();

   const hier::Box& qrh_gbox = q_rhs_data.getGhostBox();
   const hier::Index& qrhlower = qrh_gbox.lower();
   const hier::Index& qrhupper = qrh_gbox.upper();

   const hier::Box& qr_gbox = q_residual_data.getGhostBox();
   const hier::Index& qrlower = qr_gbox.lower();
   const hier::Index& qrupper = qr_gbox.upper();

   if (rotation_index) {
#if NDIM == 2
      COMPUTE_Q_RESIDUAL2D_SYMM(
          lower[0], upper[0], lower[1], upper[1], d_qlen,
          sqrt_m_data.getPointer(), mlower[0], mupper[0], mlower[1], mupper[1],
          flux_data.getPointer(0), flower[0], fupper[0] + 1, flower[1],
          fupper[1], flux_data.getPointer(1), flower[0], fupper[0], flower[1],
          fupper[1] + 1, q_soln_data.getPointer(), qlower[0], qupper[0],
          qlower[1], qupper[1], dx, d_gamma, q_rhs_data.getPointer(),
          qrhlower[0], qrhupper[0], qrhlower[1], qrhupper[1],
          q_residual_data.getPointer(), qrlower[0], qrupper[0], qrlower[1],
          qrupper[1], rotation_index->getPointer(0),
          rotation_index->getPointer(1),
          rotation_index->getGhostCellWidth()[0]);
#endif
#if NDIM == 3
      COMPUTE_Q_RESIDUAL3D_SYMM(
          lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], d_qlen,
          sqrt_m_data.getPointer(), mlower[0], mupper[0], mlower[1], mupper[1],
          mlower[2], mupper[2], flux_data.getPointer(0), flower[0],
          fupper[0] + 1, flower[1], fupper[1], flower[2], fupper[2],
          flux_data.getPointer(1), flower[0], fupper[0], flower[1],
          fupper[1] + 1, flower[2], fupper[2], flux_data.getPointer(2),
          flower[0], fupper[0], flower[1], fupper[1], flower[2], fupper[2] + 1,
          q_soln_data.getPointer(), qlower[0], qupper[0], qlower[1], qupper[1],
          qlower[2], qupper[2], dx, d_gamma, q_rhs_data.getPointer(),
          qrhlower[0], qrhupper[0], qrhlower[1], qrhupper[1], qrhlower[2],
          qrhupper[2], q_residual_data.getPointer(), qrlower[0], qrupper[0],
          qrlower[1], qrupper[1], qrlower[2], qrupper[2],
          rotation_index->getPointer(0), rotation_index->getPointer(1),
          rotation_index->getPointer(2),
          rotation_index->getGhostCellWidth()[0]);
#endif
   } else {
#if NDIM == 2
      COMPUTE_Q_RESIDUAL2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                           sqrt_m_data.getPointer(), mlower[0], mupper[0],
                           mlower[1], mupper[1], flux_data.getPointer(0),
                           flower[0], fupper[0] + 1, flower[1], fupper[1],
                           flux_data.getPointer(1), flower[0], fupper[0],
                           flower[1], fupper[1] + 1, q_soln_data.getPointer(),
                           qlower[0], qupper[0], qlower[1], qupper[1], dx,
                           d_gamma, q_rhs_data.getPointer(), qrhlower[0],
                           qrhupper[0], qrhlower[1], qrhupper[1],
                           q_residual_data.getPointer(), qrlower[0], qrupper[0],
                           qrlower[1], qrupper[1]);
#endif
#if NDIM == 3
      COMPUTE_Q_RESIDUAL3D(
          lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], d_qlen,
          sqrt_m_data.getPointer(), mlower[0], mupper[0], mlower[1], mupper[1],
          mlower[2], mupper[2], flux_data.getPointer(0), flower[0],
          fupper[0] + 1, flower[1], fupper[1], flower[2], fupper[2],
          flux_data.getPointer(1), flower[0], fupper[0], flower[1],
          fupper[1] + 1, flower[2], fupper[2], flux_data.getPointer(2),
          flower[0], fupper[0], flower[1], fupper[1], flower[2], fupper[2] + 1,
          q_soln_data.getPointer(), qlower[0], qupper[0], qlower[1], qupper[1],
          qlower[2], qupper[2], dx, d_gamma, q_rhs_data.getPointer(),
          qrhlower[0], qrhupper[0], qrhlower[1], qrhupper[1], qrhlower[2],
          qrhupper[2], q_residual_data.getPointer(), qrlower[0], qrupper[0],
          qrlower[1], qrupper[1], qrlower[2], qrupper[2]);
#endif
   }
}


void QuatFACOps::evaluateRHS(
    const int phase_id, const int temp_id, const int grad_q_id,
    const int grad_q_copy_id,  // for computation of diffusion coefficient
    const int mobility_id, const int rotation_index_id, const int q_id,
    int rhs_id, const bool use_gradq_for_flux)
{
   t_compute_rhs->start();

   if (use_gradq_for_flux) assert(grad_q_id >= 0);
   assert(grad_q_copy_id >= 0);

   d_rotation_index_id = rotation_index_id;

   const int gq_id = use_gradq_for_flux ? grad_q_id : -1;

   // Initialize the output array
   d_hopscell->setToScalar(rhs_id, 0., false);

   d_quat_face_coeff_strategy->computeFaceCoefs(d_hierarchy, phase_id, temp_id,
                                                grad_q_copy_id,
                                                d_face_coef_scratch_id);

   for (int ln = d_ln_max; ln >= d_ln_min; ln--) {
      accumulateOperatorOnLevel(mobility_id, d_face_coef_scratch_id, q_id,
                                gq_id, rhs_id, ln, true, false);
   }

   t_compute_rhs->stop();
}


void QuatFACOps::multiplyDQuatDPhiBlock(const int phase_id, const int out_id)
{
   // Initialize the output array
   d_hopscell->setToScalar(out_id, 0., false);

   // Compute the product of the dquatdphi preconditioner block containing the
   // mobility derivative times the input phi
   for (int ln = d_ln_max; ln >= d_ln_min; ln--) {

      accumulateOperatorOnLevel(d_m_deriv_id, d_face_coef_id, d_q_local_id, -1,
                                out_id, ln, false, false);

      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::Iterator pi(level->begin());
      for (; pi != level->end(); pi++) {
         hier::Patch& patch = **pi;

         std::shared_ptr<pdat::CellData<double> > out_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(out_id)));
         std::shared_ptr<pdat::CellData<double> > phi_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(phase_id)));

         const hier::Box& box = patch.getBox();
         const hier::Index& lower = box.lower();
         const hier::Index& upper = box.upper();

         const hier::Box& o_gbox = out_data->getGhostBox();
         const hier::Index& olower = o_gbox.lower();
         const hier::Index& oupper = o_gbox.upper();

         const hier::Box& p_gbox = phi_data->getGhostBox();
         const hier::Index& plower = p_gbox.lower();
         const hier::Index& pupper = p_gbox.upper();

#if NDIM == 2
         MULTICOMPONENT_MULTIPLY2D(lower[0], upper[0], lower[1], upper[1],
                                   phi_data->getPointer(), plower[0], pupper[0],
                                   plower[1], pupper[1], out_data->getPointer(),
                                   olower[0], oupper[0], olower[1], oupper[1],
                                   d_qlen);
#endif
#if NDIM == 3
         MULTICOMPONENT_MULTIPLY3D(lower[0], upper[0], lower[1], upper[1],
                                   lower[2], upper[2], phi_data->getPointer(),
                                   plower[0], pupper[0], plower[1], pupper[1],
                                   plower[2], pupper[2], out_data->getPointer(),
                                   olower[0], oupper[0], olower[1], oupper[1],
                                   olower[2], oupper[2], d_qlen);
#endif
      }
   }

   // Add the product of the dquatdphi preconditioner block with the diffusion
   // coefficient derivative times the input phi

   computeDQuatDPhiFaceCoefs(d_face_coef_deriv_id, phase_id,
                             d_face_coef_scratch_id);

   for (int ln = d_ln_max; ln >= d_ln_min; ln--) {
      accumulateOperatorOnLevel(d_sqrt_m_id, d_face_coef_scratch_id,
                                d_q_local_id, -1, out_id, ln, false, false);
   }
}


void QuatFACOps::accumulateOperatorOnLevel(const int mobility_id,
                                           const int face_coef_id,
                                           const int q_id, const int grad_q_id,
                                           int rhs_q_id, int ln, bool project,
                                           bool error_equation_indicator)
{
   checkFluxPatchDataIndex();
   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

   /*
    * Set up the bc helper so that when we use a refine schedule
    * to fill ghosts, the correct data is operated on.
    */

   d_bc_helper.setTargetDataId(q_id);
   d_bc_helper.setHomogeneousBc(error_equation_indicator);

   const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

   /*
    * Assumptions:
    * 1. Data does not yet exist in ghost boundaries.
    * 2. RHS data on next finer grid (if any)
    *    has been computed already.
    * 3. Flux data from next finer grid (if any) has
    *    been computed but has not been coarsened to
    *    this level.
    *
    * Steps:
    * S1. Fill solution ghost data by refinement
    *     or setting physical boundary conditions.
    *     This also brings in information from coarser
    *     to form the composite grid flux.
    * S2. Compute flux on ln.
    * S3. If next finer is available,
    *     Coarsen flux data on next finer level,
    *     overwriting flux computed from coarse data.
    * S4. Compute RHS data from flux.
    */

   /* S1. Fill ghost data for quaternions (component 0). */
   if (grad_q_id == -1) {
      if (ln > d_ln_min) {
         /* Fill from current, next coarser level and physical boundary */
         xeqScheduleGhostFill(q_id, ln);
      } else {
         /* Fill from current and physical boundary */
         xeqScheduleGhostFillNoCoarse(q_id, ln);
      }
   }

   /*
    * For the whole level, make sure the internal
    * side-centered data is allocated and note
    * whether that data should be deallocated when done.
    * We do this for the whole level because the data
    * undergoes transfer operations which require the
    * whole level data.
    */
   bool deallocate_flux_data_when_done = false;
   if (flux_id == d_flux_scratch_id) {
      if (!level->checkAllocated(flux_id)) {
         level->allocatePatchData(flux_id);
         deallocate_flux_data_when_done = true;
      }
   }

   /*
    * S2. Compute flux on patches in level.
    */
   hier::PatchLevel::Iterator pi(level->begin());
   for (; pi != level->end(); pi++) {
      std::shared_ptr<hier::Patch> patch = *pi;

      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));
      std::shared_ptr<pdat::SideData<double> > face_coef_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(face_coef_id)));

      if (grad_q_id == -1) {
         std::shared_ptr<pdat::CellData<double> > q_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(q_id)));

         computeFluxOnPatch(*patch, level->getRatioToCoarserLevel(),
                            *face_coef_data, *q_data, *flux_data);
      } else {
         // tbox::pout<<"call computeFluxOnPatch()..."<<endl;
         std::shared_ptr<pdat::SideData<double> > grad_q_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(grad_q_id)));
         computeFluxOnPatch(*patch, level->getRatioToCoarserLevel(),
                            *face_coef_data, *grad_q_data, *flux_data);
      }
   }

   /*
    * S3. Coarsen oflux data from next finer level so that
    * the computed flux becomes the composite grid flux.
    */
   if (ln < d_ln_max) {
      xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
   }

   /*
    * S4. Accumulate operator on patches in level.
    */
   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        pi++) {
      std::shared_ptr<hier::Patch> patch = *pi;

      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));
      std::shared_ptr<pdat::CellData<double> > mobility_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(mobility_id)));
      std::shared_ptr<pdat::CellData<double> > q_rhs_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(rhs_q_id)));

      if (project && d_qlen != 1) {

         std::shared_ptr<pdat::CellData<double> > q_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(q_id)));

         // Temporary to store lambda
         std::shared_ptr<pdat::CellData<double> > lambda_data(
             new pdat::CellData<double>(pi->getBox(), 1,
                                        hier::IntVector(tbox::Dimension(NDIM),
                                                        0)));
         assert(lambda_data);

         std::shared_ptr<pdat::SideData<int> > rotation_index;
         if (d_rotation_index_id >= 0) {
            rotation_index =
                std::dynamic_pointer_cast<pdat::SideData<int>, hier::PatchData>(
                    patch->getPatchData(d_rotation_index_id));
            assert(rotation_index);
         }

         computeLambdaOnPatch(*patch, *flux_data, *q_data, rotation_index,
                              *lambda_data);

         accumulateProjectedOperatorOnPatch(*patch, *flux_data, *mobility_data,
                                            *q_data, *lambda_data,
                                            rotation_index, *q_rhs_data);
      } else {
         accumulateOperatorOnPatch(*patch, *flux_data, *mobility_data,
                                   *q_rhs_data);
      }

      if (ln > d_ln_min) {
         /*
          * Save outerflux data so that next coarser level
          *  can compute its coarse-fine composite flux.
          *  This is not strictly needed in this "compute RHS"
          *  loop through the patches, but we put it here to
          *  avoid writing another loop for it.
          */
         std::shared_ptr<pdat::OutersideData<double> > oflux_data(
             SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>,
                                    hier::PatchData>(
                 patch->getPatchData(d_oflux_scratch_id)));
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(oflux_data);
#endif
         oflux_data->copy(*flux_data);
      }
   }

   if (deallocate_flux_data_when_done) {
      level->deallocatePatchData(flux_id);
   }
}


void QuatFACOps::accumulateOperatorOnPatch(
    const hier::Patch& patch, const pdat::SideData<double>& flux_data,
    const pdat::CellData<double>& mobility_data,
    const pdat::CellData<double>& q_rhs_data) const
{
   assert(d_rotation_index_id == -1);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));

   const double* dx = patch_geom->getDx();

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& m_gbox = mobility_data.getGhostBox();
   const hier::Index& mlower = m_gbox.lower();
   const hier::Index& mupper = m_gbox.upper();

   const hier::Box& f_gbox = flux_data.getGhostBox();
   const hier::Index& flower = f_gbox.lower();
   const hier::Index& fupper = f_gbox.upper();

   const hier::Box& qrh_gbox = q_rhs_data.getGhostBox();
   const hier::Index& qrhlower = qrh_gbox.lower();
   const hier::Index& qrhupper = qrh_gbox.upper();

#if NDIM == 2
   ADD_QUAT_OP2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                 mobility_data.getPointer(), mlower[0], mupper[0], mlower[1],
                 mupper[1], flux_data.getPointer(0), flower[0], fupper[0] + 1,
                 flower[1], fupper[1], flux_data.getPointer(1), flower[0],
                 fupper[0], flower[1], fupper[1] + 1, dx,
                 q_rhs_data.getPointer(), qrhlower[0], qrhupper[0], qrhlower[1],
                 qrhupper[1]);
#endif
#if NDIM == 3
   ADD_QUAT_OP3D(lower[0], upper[0], lower[1], upper[1], lower[2], upper[2],
                 d_qlen, mobility_data.getPointer(), mlower[0], mupper[0],
                 mlower[1], mupper[1], mlower[2], mupper[2],
                 flux_data.getPointer(0), flower[0], fupper[0] + 1, flower[1],
                 fupper[1], flower[2], fupper[2], flux_data.getPointer(1),
                 flower[0], fupper[0], flower[1], fupper[1] + 1, flower[2],
                 fupper[2], flux_data.getPointer(2), flower[0], fupper[0],
                 flower[1], fupper[1], flower[2], fupper[2] + 1, dx,
                 q_rhs_data.getPointer(), qrhlower[0], qrhupper[0], qrhlower[1],
                 qrhupper[1], qrhlower[2], qrhupper[2]);
#endif
}


void QuatFACOps::accumulateProjectedOperatorOnPatch(
    const hier::Patch& patch, const pdat::SideData<double>& flux_data,
    const pdat::CellData<double>& mobility_data,
    const pdat::CellData<double>& q_soln_data,
    const pdat::CellData<double>& lambda_soln_data,
    std::shared_ptr<pdat::SideData<int> > rotation_index,
    const pdat::CellData<double>& q_rhs_data) const
{
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));

   const double* dx = patch_geom->getDx();

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& m_gbox = mobility_data.getGhostBox();
   const hier::Index& mlower = m_gbox.lower();
   const hier::Index& mupper = m_gbox.upper();

   const hier::Box& f_gbox = flux_data.getGhostBox();
   const hier::Index& flower = f_gbox.lower();
   const hier::Index& fupper = f_gbox.upper();

   const hier::Box& q_gbox = q_soln_data.getGhostBox();
   const hier::Index& qlower = q_gbox.lower();
   const hier::Index& qupper = q_gbox.upper();

   const hier::Box& l_gbox = lambda_soln_data.getGhostBox();
   const hier::Index& llower = l_gbox.lower();
   const hier::Index& lupper = l_gbox.upper();

   const hier::Box& qrh_gbox = q_rhs_data.getGhostBox();
   const hier::Index& qrhlower = qrh_gbox.lower();
   const hier::Index& qrhupper = qrh_gbox.upper();

#if NDIM == 2
   if (rotation_index)
      ADD_QUAT_PROJ_OP2D_SYMM(
          lower[0], upper[0], lower[1], upper[1], d_qlen,
          mobility_data.getPointer(), mlower[0], mupper[0], mlower[1],
          mupper[1], flux_data.getPointer(0), flower[0], fupper[0] + 1,
          flower[1], fupper[1], flux_data.getPointer(1), flower[0], fupper[0],
          flower[1], fupper[1] + 1, q_soln_data.getPointer(), qlower[0],
          qupper[0], qlower[1], qupper[1], lambda_soln_data.getPointer(),
          llower[0], lupper[0], llower[1], lupper[1], dx,
          q_rhs_data.getPointer(), qrhlower[0], qrhupper[0], qrhlower[1],
          qrhupper[1], rotation_index->getPointer(0),
          rotation_index->getPointer(1),
          rotation_index->getGhostCellWidth()[0]);
   else
      ADD_QUAT_PROJ_OP2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                         mobility_data.getPointer(), mlower[0], mupper[0],
                         mlower[1], mupper[1], flux_data.getPointer(0),
                         flower[0], fupper[0] + 1, flower[1], fupper[1],
                         flux_data.getPointer(1), flower[0], fupper[0],
                         flower[1], fupper[1] + 1, q_soln_data.getPointer(),
                         qlower[0], qupper[0], qlower[1], qupper[1],
                         lambda_soln_data.getPointer(), llower[0], lupper[0],
                         llower[1], lupper[1], dx, q_rhs_data.getPointer(),
                         qrhlower[0], qrhupper[0], qrhlower[1], qrhupper[1]);
#endif
#if NDIM == 3
   assert(!rotation_index);
   ADD_QUAT_PROJ_OP3D(lower[0], upper[0], lower[1], upper[1], lower[2],
                      upper[2], d_qlen, mobility_data.getPointer(), mlower[0],
                      mupper[0], mlower[1], mupper[1], mlower[2], mupper[2],
                      flux_data.getPointer(0), flower[0], fupper[0] + 1,
                      flower[1], fupper[1], flower[2], fupper[2],
                      flux_data.getPointer(1), flower[0], fupper[0], flower[1],
                      fupper[1] + 1, flower[2], fupper[2],
                      flux_data.getPointer(2), flower[0], fupper[0], flower[1],
                      fupper[1], flower[2], fupper[2] + 1,
                      q_soln_data.getPointer(), qlower[0], qupper[0], qlower[1],
                      qupper[1], qlower[2], qupper[2],
                      lambda_soln_data.getPointer(), llower[0], lupper[0],
                      llower[1], lupper[1], llower[2], lupper[2], dx,
                      q_rhs_data.getPointer(), qrhlower[0], qrhupper[0],
                      qrhlower[1], qrhupper[1], qrhlower[2], qrhupper[2]);
#endif
}


void QuatFACOps::multiplyMobilitySqrt(const int id)
{
   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::Iterator pi(level->begin());
      for (; pi != level->end(); pi++) {
         hier::Patch& patch = **pi;

         std::shared_ptr<pdat::CellData<double> > data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(id)));
         std::shared_ptr<pdat::CellData<double> > sqrt_m_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(d_sqrt_m_id)));
         int ncomp = data->getDepth();

         const hier::Box& gbox = data->getGhostBox();
         const hier::Index& lower = gbox.lower();
         const hier::Index& upper = gbox.upper();

         const hier::Box& m_gbox = sqrt_m_data->getGhostBox();
         const hier::Index& mlower = m_gbox.lower();
         const hier::Index& mupper = m_gbox.upper();

#if NDIM == 2
         MULTICOMPONENT_MULTIPLY2D(lower[0], upper[0], lower[1], upper[1],
                                   sqrt_m_data->getPointer(), mlower[0],
                                   mupper[0], mlower[1], mupper[1],
                                   data->getPointer(), lower[0], upper[0],
                                   lower[1], upper[1], ncomp);
#endif
#if NDIM == 3
         MULTICOMPONENT_MULTIPLY3D(lower[0], upper[0], lower[1], upper[1],
                                   lower[2], upper[2],
                                   sqrt_m_data->getPointer(), mlower[0],
                                   mupper[0], mlower[1], mupper[1], mlower[2],
                                   mupper[2], data->getPointer(), lower[0],
                                   upper[0], lower[1], upper[1], lower[2],
                                   upper[2], ncomp);
#endif
      }
   }
}


void QuatFACOps::divideMobilitySqrt(const int id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   math::ArrayDataNormOpsReal<double> ops;
#endif

   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
      hier::PatchLevel::Iterator pi(level->begin());
      for (; pi != level->end(); pi++) {
         hier::Patch& patch = **pi;

         std::shared_ptr<pdat::CellData<double> > data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(id)));
         std::shared_ptr<pdat::CellData<double> > sqrt_m_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(d_sqrt_m_id)));
         int ncomp = data->getDepth();

#ifdef DEBUG_CHECK_ASSERTIONS
         double nb = ops.maxNorm(data->getArrayData(), patch.getBox());
         assert(nb == nb);
#endif

         const hier::Box& box = patch.getBox();
         const hier::Index& lower = box.lower();
         const hier::Index& upper = box.upper();

         const hier::Box& d_gbox = data->getGhostBox();
         const hier::Index& dlower = d_gbox.lower();
         const hier::Index& dupper = d_gbox.upper();

         const hier::Box& m_gbox = sqrt_m_data->getGhostBox();
         const hier::Index& mlower = m_gbox.lower();
         const hier::Index& mupper = m_gbox.upper();

#if NDIM == 2
         MULTICOMPONENT_DIVIDE2D(lower[0], upper[0], lower[1], upper[1],
                                 sqrt_m_data->getPointer(), mlower[0],
                                 mupper[0], mlower[1], mupper[1],
                                 data->getPointer(), dlower[0], dupper[0],
                                 dlower[1], dupper[1], ncomp);
#endif
#if NDIM == 3
         MULTICOMPONENT_DIVIDE3D(lower[0], upper[0], lower[1], upper[1],
                                 lower[2], upper[2], sqrt_m_data->getPointer(),
                                 mlower[0], mupper[0], mlower[1], mupper[1],
                                 mlower[2], mupper[2], data->getPointer(),
                                 dlower[0], dupper[0], dlower[1], dupper[1],
                                 dlower[2], dupper[2], ncomp);
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
         nb = ops.maxNorm(data->getArrayData(), patch.getBox());
         assert(nb == nb);
#endif
      }
   }
}


void QuatFACOps::applyProjectionOnLevel(const int q_id, const int corr_id,
                                        const int err_id, const int ln)
{
   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
   hier::PatchLevel::Iterator pi(level->begin());
   for (; pi != level->end(); pi++) {
      hier::Patch& patch = **pi;

      std::shared_ptr<pdat::CellData<double> > q_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(q_id)));
      std::shared_ptr<pdat::CellData<double> > corr_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(corr_id)));
      std::shared_ptr<pdat::CellData<double> > err_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(err_id)));

      const hier::Box& box = patch.getBox();
      const hier::Index& lower = box.lower();
      const hier::Index& upper = box.upper();

      const hier::Box& q_gbox = q_data->getGhostBox();
      const hier::Index& qlower = q_gbox.lower();
      const hier::Index& qupper = q_gbox.upper();

      const hier::Box& c_gbox = corr_data->getGhostBox();
      const hier::Index& clower = c_gbox.lower();
      const hier::Index& cupper = c_gbox.upper();

      const hier::Box& e_gbox = err_data->getGhostBox();
      const hier::Index& elower = e_gbox.lower();
      const hier::Index& eupper = e_gbox.upper();

#if NDIM == 2
      PROJECT2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                q_data->getPointer(), qlower[0], qupper[0], qlower[1],
                qupper[1], corr_data->getPointer(), clower[0], cupper[0],
                clower[1], cupper[1], err_data->getPointer(), elower[0],
                eupper[0], elower[1], eupper[1]);
#endif
#if NDIM == 3
      PROJECT3D(lower[0], upper[0], lower[1], upper[1], lower[2], upper[2],
                d_qlen, q_data->getPointer(), qlower[0], qupper[0], qlower[1],
                qupper[1], qlower[2], qupper[2], corr_data->getPointer(),
                clower[0], cupper[0], clower[1], cupper[1], clower[2],
                cupper[2], err_data->getPointer(), elower[0], eupper[0],
                elower[1], eupper[1], elower[2], eupper[2]);
#endif
   }
}

/*
*******************************************************************
*                                                                 *
* Communication members.                                          *
*                                                                 *
*******************************************************************
*/

void QuatFACOps::xeqScheduleProlongation(int dst_id, int src_id, int scr_id,
                                         int dest_ln)
{
   // tbox::plog<<"xeqScheduleProlongation for component "<<component<<",
   // dest_ln="<<dest_ln<<endl;

   if (!d_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm refiner;

   refiner.registerRefine(dst_id, src_id, scr_id,
                          d_prolongation_refine_operator);
   refiner.resetSchedule(d_prolongation_refine_schedules[dest_ln]);
   d_prolongation_refine_schedules[dest_ln]->fillData(0.0);
   d_prolongation_refine_algorithm->resetSchedule(
       d_prolongation_refine_schedules[dest_ln]);
}


void QuatFACOps::xeqScheduleURestriction(int dst_id, int src_id, int dest_ln)
{
   if (!d_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm coarsener(tbox::Dimension(NDIM));
   coarsener.registerCoarsen(dst_id, src_id, d_urestriction_coarsen_operator);
   coarsener.resetSchedule(d_urestriction_coarsen_schedules[dest_ln]);
   d_urestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_urestriction_coarsen_algorithm->resetSchedule(
       d_urestriction_coarsen_schedules[dest_ln]);
}


void QuatFACOps::xeqScheduleRRestriction(int dst_id, int src_id, int dest_ln)
{
   if (!d_rrestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm coarsener(tbox::Dimension(NDIM));
   coarsener.registerCoarsen(dst_id, src_id, d_rrestriction_coarsen_operator);
   coarsener.resetSchedule(d_rrestriction_coarsen_schedules[dest_ln]);
   d_rrestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_rrestriction_coarsen_algorithm->resetSchedule(
       d_rrestriction_coarsen_schedules[dest_ln]);
}


void QuatFACOps::xeqScheduleFluxCoarsen(int dst_id, int src_id, int dest_ln)
{
   if (!d_flux_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm coarsener(tbox::Dimension(NDIM));
   coarsener.registerCoarsen(dst_id, src_id, d_flux_coarsen_operator);
   coarsener.resetSchedule(d_flux_coarsen_schedules[dest_ln]);
   d_flux_coarsen_schedules[dest_ln]->coarsenData();
   d_flux_coarsen_algorithm->resetSchedule(d_flux_coarsen_schedules[dest_ln]);
}


void QuatFACOps::xeqScheduleGhostFill(int dst_id, int dest_ln)
{
   // tbox::plog<<"QuatFACOps::xeqScheduleGhostFill() for
   // dest_ln="<<dest_ln<<endl;
   if (!d_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm refiner;

   refiner.registerRefine(dst_id, dst_id, dst_id, d_ghostfill_refine_operator);
   refiner.resetSchedule(d_ghostfill_refine_schedules[dest_ln]);
   d_ghostfill_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_refine_algorithm->resetSchedule(
       d_ghostfill_refine_schedules[dest_ln]);
}


/*
 * Function called even with single level
 */
void QuatFACOps::xeqScheduleGhostFillNoCoarse(int dst_id, int dest_ln)
{
   // tbox::plog<<"QuatFACOps::xeqScheduleGhostFillNoCoarse()..."<<endl;
   if (!d_ghostfill_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   assert(d_ghostfill_nocoarse_refine_operator);

   xfer::RefineAlgorithm refiner;
   refiner.registerRefine(dst_id, dst_id, dst_id,
                          d_ghostfill_nocoarse_refine_operator);
   refiner.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dest_ln]);
   d_ghostfill_nocoarse_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_nocoarse_refine_algorithm->resetSchedule(
       d_ghostfill_nocoarse_refine_schedules[dest_ln]);
}


/*
*******************************************************************
*                                                                 *
* Miscellaneous helper members.                                   *
*                                                                 *
*******************************************************************
*/

void QuatFACOps::freeVariables()
{
   s_cell_scratch_var.reset();
   s_flux_scratch_var.reset();
   s_oflux_scratch_var.reset();
   s_face_coef_var.reset();
   s_q_local_var.reset();
   s_residual_var.reset();
}


int QuatFACOps::GetNumCellFacesInBox(const int* lower, const int* upper,
                                     const int dim) const
{
   // Return the number of cell faces in the box (lower, upper)

   int num_cell_faces = 1;
   for (int i = 0; i < NDIM; i++) {
      int extra = (i == dim) ? 1 : 0;
      num_cell_faces *= (upper[i] - lower[i] + 1 + extra);
   }

   return num_cell_faces;
}
