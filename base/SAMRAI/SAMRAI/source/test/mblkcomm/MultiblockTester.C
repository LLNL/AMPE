/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Manager class for patch data communication tests.
 *
 ************************************************************************/

#include "MultiblockTester.h"

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "PatchMultiblockTestStrategy.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/MultiblockGriddingTagger.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

using namespace SAMRAI;

/*
 *************************************************************************
 *
 * The constructor initializes object state.  The destructor is empty.
 *
 *************************************************************************
 */

MultiblockTester::MultiblockTester(
   const string& object_name,
   const tbox::Dimension& dim,
   boost::shared_ptr<tbox::Database>& main_input_db,
   boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   PatchMultiblockTestStrategy* data_test,
   bool do_refine,
   bool do_coarsen,
   const string& refine_option):
   xfer::CoarsenPatchStrategy(dim),
   xfer::RefinePatchStrategy(dim),
   d_object_name(object_name),
   d_dim(dim),
   d_data_test_strategy(data_test),
   d_do_refine(do_refine),
   d_do_coarsen(false),
   d_refine_option(refine_option),
   d_patch_hierarchy(hierarchy),
   d_fake_time(0.0),
   d_source(
     hier::VariableDatabase::getDatabase()->getContext("SOURCE")),
   d_destination(
      hier::VariableDatabase::getDatabase()->getContext("DESTINATION")),
   d_refine_scratch(
      hier::VariableDatabase::getDatabase()->getContext("REFINE_SCRATCH")),
   d_reset_source(
      hier::VariableDatabase::getDatabase()->getContext("SOURCE")),
   d_reset_destination(
      hier::VariableDatabase::getDatabase()->getContext("DESTINATION")),
   d_reset_refine_scratch(
      hier::VariableDatabase::getDatabase()->getContext("REFINE_SCRATCH")),
   d_reset_refine_algorithm(dim),
   d_reset_coarsen_algorithm(dim),
   d_is_reset(false)
{
   NULL_USE(main_input_db);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(data_test != (PatchMultiblockTestStrategy *)NULL);
#endif

   if (!do_refine) {
      d_do_coarsen = do_coarsen;
   }

   if (!((d_refine_option == "INTERIOR_FROM_SAME_LEVEL")
         || (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL"))) {
      TBOX_ERROR(object_name << " input error: illegal refine_option = "
                             << d_refine_option << endl);
   }

   d_patch_data_components.clrAllFlags();

   d_data_test_strategy->registerVariables(this);
}

MultiblockTester::~MultiblockTester()
{

}

/*
 *************************************************************************
 *
 * Add variable with associated attributes to set of test variables.
 *
 *************************************************************************
 */

void MultiblockTester::registerVariable(
   const boost::shared_ptr<hier::Variable> src_variable,
   const boost::shared_ptr<hier::Variable> dst_variable,
   const hier::IntVector& src_ghosts,
   const hier::IntVector& dst_ghosts,
   const boost::shared_ptr<hier::BaseGridGeometry> xfer_geom,
   const string& operator_name)
{
   TBOX_ASSERT(src_variable);
   TBOX_ASSERT(dst_variable);
   TBOX_ASSERT(xfer_geom);
   TBOX_ASSERT(!operator_name.empty());

   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   int src_id = variable_db->registerVariableAndContext(src_variable,
         d_source,
         src_ghosts);

   int dst_id = variable_db->registerVariableAndContext(dst_variable,
         d_destination,
         dst_ghosts);

   TBOX_ASSERT(src_id != -1);
   TBOX_ASSERT(dst_id != -1);

   d_patch_data_components.setFlag(src_id);
   d_patch_data_components.setFlag(dst_id);

   boost::shared_ptr<hier::RefineOperator> refine_operator;
   boost::shared_ptr<hier::CoarsenOperator> coarsen_operator;

   if (d_do_refine) {
      refine_operator = xfer_geom->lookupRefineOperator(src_variable,
            operator_name);

      d_mblk_refine_alg.reset(new xfer::RefineAlgorithm(d_dim));

      hier::IntVector scratch_ghosts =
         hier::IntVector::max(src_ghosts, dst_ghosts);
      scratch_ghosts.max(hier::IntVector(d_dim, 1));
      if (refine_operator) {
         scratch_ghosts.max(refine_operator->getStencilWidth());
      }
      int scratch_id =
         variable_db->registerVariableAndContext(src_variable,
            d_refine_scratch,
            scratch_ghosts);
      TBOX_ASSERT(scratch_id != -1);

      d_patch_data_components.setFlag(scratch_id);

      d_mblk_refine_alg->registerRefine(dst_id,
         src_id,
         scratch_id,
         refine_operator);

   } else if (d_do_coarsen) {
      coarsen_operator = xfer_geom->lookupCoarsenOperator(src_variable,
            operator_name);
      d_coarsen_algorithm->registerCoarsen(dst_id,
         src_id,
         coarsen_operator);

   }

   registerVariableForReset(src_variable, dst_variable,
      src_ghosts, dst_ghosts, xfer_geom,
      operator_name);

}

void MultiblockTester::registerVariableForReset(
   const boost::shared_ptr<hier::Variable> src_variable,
   const boost::shared_ptr<hier::Variable> dst_variable,
   const hier::IntVector& src_ghosts,
   const hier::IntVector& dst_ghosts,
   const boost::shared_ptr<hier::BaseGridGeometry> xfer_geom,
   const string& operator_name)
{
   TBOX_ASSERT(src_variable);
   TBOX_ASSERT(dst_variable);
   TBOX_ASSERT(xfer_geom);
   TBOX_ASSERT(!operator_name.empty());

   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   int src_id = variable_db->registerVariableAndContext(src_variable,
         d_reset_source,
         src_ghosts);

   int dst_id = variable_db->registerVariableAndContext(dst_variable,
         d_reset_destination,
         dst_ghosts);

   d_patch_data_components.setFlag(src_id);
   d_patch_data_components.setFlag(dst_id);

   boost::shared_ptr<hier::RefineOperator> refine_operator;
   boost::shared_ptr<hier::CoarsenOperator> coarsen_operator;

   if (d_do_refine) {
      refine_operator = xfer_geom->lookupRefineOperator(src_variable,
            operator_name);

      hier::IntVector scratch_ghosts =
         hier::IntVector::max(src_ghosts, dst_ghosts);

      scratch_ghosts.max(hier::IntVector(d_dim, 1));
      if (refine_operator) {
         scratch_ghosts.max(refine_operator->getStencilWidth());
      }
      int scratch_id =
         variable_db->registerVariableAndContext(src_variable,
            d_reset_refine_scratch,
            scratch_ghosts);

      d_patch_data_components.setFlag(scratch_id);

      d_reset_refine_algorithm.registerRefine(dst_id,
         src_id,
         scratch_id,
         refine_operator);

   } else if (d_do_coarsen) {
      coarsen_operator = xfer_geom->lookupCoarsenOperator(src_variable,
            operator_name);
      d_reset_coarsen_algorithm.registerCoarsen(dst_id,
         src_id,
         coarsen_operator);
   }

}

/*
 *************************************************************************
 *
 * Create refine and coarsen communication schedules for hierarchy.
 *
 *************************************************************************
 */

void MultiblockTester::createRefineSchedule(
   const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif

   boost::shared_ptr<hier::PatchLevel> level(
      d_patch_hierarchy->getPatchLevel(level_number));

   if (d_do_refine) {

      d_refine_schedule.resizeArray(
         d_patch_hierarchy->getFinestLevelNumber() + 1);
      d_refine_schedule[level_number].reset();

      if (level_number == 0) {
         d_refine_schedule[level_number] =
            d_mblk_refine_alg->createSchedule(level,
               this);
      } else if (d_refine_option == "INTERIOR_FROM_SAME_LEVEL") {
         d_refine_schedule[level_number] =
            d_mblk_refine_alg->createSchedule(level,
               level_number - 1,
               d_patch_hierarchy,
               this);
      } else if (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") {
         d_refine_schedule[level_number] =
            d_mblk_refine_alg->createSchedule(level,
               boost::shared_ptr<hier::PatchLevel>(),
               level_number - 1,
               d_patch_hierarchy,
               this);
      }

   }

}

void MultiblockTester::resetRefineSchedule(
   const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif

   if (d_do_refine) {

      d_reset_refine_algorithm.resetSchedule(d_refine_schedule[level_number]);

   }

   d_is_reset = true;
}

void MultiblockTester::createCoarsenSchedule(
   const int level_number)
{
   NULL_USE(level_number);
/*
 * if (d_do_coarsen && (level_number > 0)) {
 *
 *    d_coarsen_schedule.resizeArray(
 *       d_patch_hierarchy->getFinestLevelNumber()+1);
 *    d_coarsen_schedule[level_number].reset();
 *
 *
 *
 *    boost::shared_ptr<hier::PatchLevel > level =
 *       d_patch_hierarchy->getPatchLevel(level_number);
 *    boost::shared_ptr<hier::PatchLevel > coarser_level =
 *       d_patch_hierarchy->getPatchLevel(level_number-1);
 *
 *    tbox::Array< boost::shared_ptr< hier::Connector > > fine_to_coarse;
 *    tbox::Array< boost::shared_ptr< hier::Connector > > coarse_to_fine;
 *
 *    const hier::Connector *fine_to_coarse =
 *       &lh.getConnector(level_number, level_number-1);
 *    const hier::Connector *coarse_to_fine =
 *       &lh.getConnector(level_number-1, level_number);
 *    const hier::Connector::TransposePair
 *       coarse_fine_pair( coarse_to_fine, fine_to_coarse );
 *
 *    if ( dlbg_schedule ) {
 *       d_coarsen_schedule[level_number] =
 *          d_coarsen_algorithm.createSchedule(coarser_level,
 *                                             level,
 *                                             coarse_fine_pair,
 *                                             this);
 *    } else {
 *       TBOX_ERROR("The following must be replaced with the DLBG version.");
 * #if 0
 *    d_coarsen_schedule[level_number] =
 *       d_coarsen_algorithm->createSchedule(coarser_level, level, this);
 * #endif
 *    }
 *
 * }
 */
}

void MultiblockTester::resetCoarsenSchedule(
   const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_patch_hierarchy->getFinestLevelNumber()));
#endif

   if (d_do_coarsen && (level_number > 0)) {

      d_reset_coarsen_algorithm.resetSchedule(
         d_coarsen_schedule[level_number]);

   }

   d_is_reset = true;
}

/*
 *************************************************************************
 *
 * Perform data refine and coarsen operations.
 *
 *************************************************************************
 */

void MultiblockTester::performRefineOperations(
   const int level_number)
{
   if (d_do_refine) {
      if (d_is_reset) {
         d_data_test_strategy->setDataContext(d_reset_refine_scratch);
      } else {
         d_data_test_strategy->setDataContext(d_destination);
      }
      if (d_refine_schedule[level_number]) {
         d_refine_schedule[level_number]->fillData(d_fake_time);
      }
      d_data_test_strategy->clearDataContext();
   }
}

void MultiblockTester::performCoarsenOperations(
   const int level_number)
{
   if (d_do_coarsen) {
      if (d_is_reset) {
         d_data_test_strategy->setDataContext(d_reset_source);
      } else {
         d_data_test_strategy->setDataContext(d_source);
      }
      if (d_coarsen_schedule[level_number]) {
         d_coarsen_schedule[level_number]->coarsenData();
      }
      d_data_test_strategy->clearDataContext();
   }
}

/*
 *************************************************************************
 *
 * Verify results of communication operations.
 *
 *************************************************************************
 */

bool MultiblockTester::verifyCommunicationResults() const
{
   bool success = true;
   if (d_is_reset) {
      d_data_test_strategy->setDataContext(d_reset_destination);
   } else {
      d_data_test_strategy->setDataContext(d_destination);
   }
   for (int ln = 0;
        ln <= d_patch_hierarchy->getFinestLevelNumber(); ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_patch_hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator mi(level->begin());
           mi != level->end(); ++mi) {

         success = d_data_test_strategy->verifyResults(
               **mi, d_patch_hierarchy, ln,
               mi->getBox().getBlockId());
      }

   }
   d_data_test_strategy->clearDataContext();

   return success;
}

/*
 *************************************************************************
 *
 * Cell tagging and patch level data initialization routines declared
 * in the GradientDetectorStrategy interface.  They are used to
 * construct the hierarchy initially.
 *
 *************************************************************************
 */
void MultiblockTester::initializeLevelData(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const boost::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(old_level);
   NULL_USE(allocate_data);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(level_number >= 0);
#endif

   boost::shared_ptr<hier::PatchHierarchy> mblk_hierarchy(hierarchy);

   boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   level->allocatePatchData(d_patch_data_components, time);

   for (hier::PatchLevel::iterator p(level->begin());
        p != level->end(); ++p) {
      const boost::shared_ptr<hier::Patch>& patch = *p;
      const hier::BlockId& block_id = patch->getBox().getBlockId();

      int level_num = level->getLevelNumber();

      d_data_test_strategy->setDataContext(d_source);
      d_data_test_strategy->initializeDataOnPatch(*patch,
         mblk_hierarchy,
         level_num, block_id,
         's');
      d_data_test_strategy->clearDataContext();

      d_data_test_strategy->setDataContext(d_reset_source);
      d_data_test_strategy->initializeDataOnPatch(*patch,
         mblk_hierarchy,
         level_num, block_id,
         's');
      d_data_test_strategy->clearDataContext();

      if (d_do_coarsen) {

         d_data_test_strategy->setDataContext(d_destination);
         d_data_test_strategy->initializeDataOnPatch(*patch,
            mblk_hierarchy,
            level_num, block_id,
            'd');
         d_data_test_strategy->clearDataContext();

         d_data_test_strategy->setDataContext(d_reset_destination);
         d_data_test_strategy->initializeDataOnPatch(*patch,
            mblk_hierarchy,
            level_num, block_id,
            'd');
         d_data_test_strategy->clearDataContext();

      }
   }

}

void MultiblockTester::resetHierarchyConfiguration(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
}

void MultiblockTester::applyGradientDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double dt_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(dt_time);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif

   boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   d_data_test_strategy->setDataContext(d_source);

   for (hier::PatchLevel::iterator p(level->begin());
        p != level->end(); ++p) {
      const boost::shared_ptr<hier::Patch>& patch = *p;

      d_data_test_strategy->tagCellsToRefine(*patch,
         hierarchy,
         level_number,
         tag_index);
   }

   d_data_test_strategy->clearDataContext();

}

/*
 *************************************************************************
 *
 * Physical boundary condition and user-defined coarsen and refine
 * operations declared in RefinePatchStrategy and CoarsenPatchStrategy.
 * They are passed off to patch data test object.
 *
 *************************************************************************
 */

void MultiblockTester::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw)
{
   NULL_USE(time);
   boost::shared_ptr<hier::VariableContext> save_context(
      d_data_test_strategy->getDataContext());

   d_data_test_strategy->setDataContext(d_refine_scratch);

   d_data_test_strategy->setPhysicalBoundaryConditions(patch,
      d_fake_time,
      gcw);

   d_data_test_strategy->setDataContext(save_context);

}

void MultiblockTester::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   const hier::Connector& dst_to_encon,
   const double time,
   const hier::Box& fill_box,
   const hier::BoundaryBox& boundary_box,
   const boost::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
   NULL_USE(grid_geometry);
   NULL_USE(time);

   boost::shared_ptr<hier::VariableContext> save_context(
      d_data_test_strategy->getDataContext());

//   if (d_filling_coarse_scratch) {
   d_data_test_strategy->setDataContext(d_refine_scratch);
//   } else {
//      d_data_test_strategy->setDataContext(d_destination);
//   }

   d_data_test_strategy->fillSingularityBoundaryConditions(
      patch,
      encon_level,
      dst_to_encon,
      fill_box,
      boundary_box,
      grid_geometry);

   d_data_test_strategy->setDataContext(save_context);
}

hier::IntVector MultiblockTester::getRefineOpStencilWidth() const
{
   return hier::IntVector(d_dim, 0);
}

void MultiblockTester::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->preprocessRefine(fine, coarse, d_refine_scratch,
      fine_box, ratio);
}

void MultiblockTester::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->postprocessRefine(fine, coarse, d_refine_scratch,
      fine_box, ratio);
}

hier::IntVector MultiblockTester::getCoarsenOpStencilWidth() const
{
   return hier::IntVector(d_dim, 0);
}

void MultiblockTester::preprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->preprocessCoarsen(coarse, fine,
      boost::shared_ptr<hier::VariableContext>(),
      coarse_box, ratio);
}

void MultiblockTester::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   d_data_test_strategy->postprocessCoarsen(coarse, fine,
      boost::shared_ptr<hier::VariableContext>(),
      coarse_box, ratio);
}

/*
 *************************************************************************
 *
 * Create and configure gridding objects used to build the hierarchy.
 * Then, create hierarchy and initialize data.  Note this routine
 * must be called after variables are registered with this tester object.
 *
 *************************************************************************
 */

void MultiblockTester::setupHierarchy(
   boost::shared_ptr<tbox::Database> main_input_db,
   boost::shared_ptr<mesh::StandardTagAndInitialize> cell_tagger)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(main_input_db);
#endif

   boost::shared_ptr<mesh::BergerRigoutsos> box_generator(
      new mesh::BergerRigoutsos(d_dim));

   boost::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
      new mesh::TreeLoadBalancer(d_dim,
         "TreeLoadBalancer",
         main_input_db->getDatabase("TreeLoadBalancer")));
   load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

   boost::shared_ptr<mesh::GriddingAlgorithm> gridding_alg(
      new mesh::GriddingAlgorithm(
         d_patch_hierarchy,
         "GriddingAlgorithm",
         main_input_db->getDatabase("GriddingAlgorithm"),
         cell_tagger,
         box_generator,
         load_balancer,
         load_balancer,
         true));

   int fake_tag_buffer = 0;

   gridding_alg->makeCoarsestLevel(d_fake_time);

   bool initial_time = true;
   for (int ln = 0; d_patch_hierarchy->levelCanBeRefined(ln); ln++) {
      gridding_alg->makeFinerLevel(d_fake_time,
         initial_time, fake_tag_buffer,
         d_fake_time);
   }

   tbox::plog << "\n\nHierarchy:\n";
   d_patch_hierarchy->recursivePrint(tbox::plog, "", 2);

   for (int ln = 1; ln < d_patch_hierarchy->getNumberOfLevels(); ln++) {
      hier::CoarseFineBoundary cf_bndry(*d_patch_hierarchy, ln,
                                        hier::IntVector::getOne(d_dim));
   }
}
