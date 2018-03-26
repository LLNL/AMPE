/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines.
 *
 ************************************************************************/

#ifndef included_mesh_GriddingAlgorithm_C
#define included_mesh_GriddingAlgorithm_C

#include "SAMRAI/mesh/GriddingAlgorithm.h"

#include "SAMRAI/tbox/IEEE.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/math/PatchCellDataOpsInteger.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <iomanip>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int GriddingAlgorithm::ALGS_GRIDDING_ALGORITHM_VERSION = 3;

/*
 *************************************************************************
 *
 * Initialize the static data members.
 *
 *************************************************************************
 */

tbox::Array<int> * GriddingAlgorithm::s_tag_indx = 0;
tbox::Array<int> * GriddingAlgorithm::s_buf_tag_indx = 0;

tbox::StartupShutdownManager::Handler
GriddingAlgorithm::s_startup_shutdown_handler(
   0,
   GriddingAlgorithm::startupCallback,
   GriddingAlgorithm::shutdownCallback,
   0,
   tbox::StartupShutdownManager::priorityListElements);

tbox::StartupShutdownManager::Handler
GriddingAlgorithm::s_initialize_handler(
   GriddingAlgorithm::initializeCallback,
   0,
   0,
   GriddingAlgorithm::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *
 * Constructor and destructor for GriddingAlgorithm.
 *
 *************************************************************************
 */
GriddingAlgorithm::GriddingAlgorithm(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const std::string& object_name,
   const boost::shared_ptr<tbox::Database>& input_db,
   const boost::shared_ptr<TagAndInitializeStrategy>& tag_init_strategy,
   const boost::shared_ptr<BoxGeneratorStrategy>& generator,
   const boost::shared_ptr<LoadBalanceStrategy>& balancer,
   const boost::shared_ptr<LoadBalanceStrategy>& balancer0,
   bool register_for_restart):
   GriddingAlgorithmStrategy(),
   d_hierarchy(hierarchy),
   d_connector_width_requestor(),
   d_dim(hierarchy->getDim()),
   d_buf_tag_ghosts(d_dim, 0),
   d_object_name(object_name),
   d_registered_for_restart(register_for_restart),
   d_tag_init_strategy(tag_init_strategy),
   d_box_generator(generator),
   d_load_balancer(balancer),
   d_load_balancer0(balancer0 ? balancer0 : balancer),
   d_mb_tagger_strategy(NULL),
   d_true_tag(1),
   d_false_tag(0),
   d_base_ln(-1),
   d_check_nonrefined_tags('w'),
   d_check_overlapping_patches('w'),
   d_check_nonnesting_user_boxes('e'),
   d_check_boundary_proximity_violation('e'),
   d_sequentialize_patch_indices(true),
   d_log_metadata_statistics(false),
   d_enforce_proper_nesting(true),
   d_extend_to_domain_boundary(true),
   d_load_balance(true),
   d_barrier_and_time(false),
   d_check_overflow_nesting(false),
   d_check_proper_nesting(false),
   d_check_connectors(false),
   d_print_steps(false)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);
   TBOX_ASSERT(tag_init_strategy);
   TBOX_ASSERT(generator);
   TBOX_ASSERT(balancer);

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(d_object_name, this);
   }

   d_hierarchy->registerConnectorWidthRequestor(
      d_connector_width_requestor);

   /*
    * Construct integer tag variables and add to variable database.
    * Note that variables and patch data indices are shared among
    * all GriddingAlgorithm instances.  The VariableDatabase
    * holds the variables, once contructed and registered via the
    * VariableDatabase::makeInternalSAMRAIWorkVariablePatchDataIndex()
    * function call.  Note that variables are registered and patch data
    * indices are made only for the first time through the constructor.
    */
   hier::VariableDatabase* var_db =
      hier::VariableDatabase::getDatabase();

   std::string tag_interior_variable_name("GriddingAlgorithm__tag-interior");
   std::string tag_buffer_variable_name("GriddingAlgorithm__tag-buffer");

   std::ostringstream dim_extension;
   dim_extension << "_" << d_dim.getValue();

   tag_interior_variable_name += dim_extension.str();
   tag_buffer_variable_name += dim_extension.str();

   d_tag = boost::dynamic_pointer_cast<pdat::CellVariable<int>, hier::Variable>(
       var_db->getVariable(tag_interior_variable_name));
   if (!d_tag) {
      d_tag.reset(
         new pdat::CellVariable<int>(d_dim, tag_interior_variable_name, 1));
   }

   d_buf_tag = boost::dynamic_pointer_cast<pdat::CellVariable<int>, hier::Variable>(
      var_db->getVariable(tag_buffer_variable_name));
   if (!d_buf_tag) {
      d_buf_tag.reset(new pdat::CellVariable<int>(d_dim,
                                                  tag_buffer_variable_name,
                                                  1));
   }

   if ((*s_tag_indx)[d_dim.getValue() - 1] < 0) {
      (*s_tag_indx)[d_dim.getValue() - 1] =
         var_db->registerInternalSAMRAIVariable(d_tag,
            hier::IntVector::getZero(d_dim));
   }
   if ((*s_buf_tag_indx)[d_dim.getValue() - 1] < 0) {
      (*s_buf_tag_indx)[d_dim.getValue() - 1] =
         var_db->registerInternalSAMRAIVariable(d_buf_tag,
            hier::IntVector::getOne(d_dim));
      d_buf_tag_ghosts = hier::IntVector::getOne(d_dim); 
   }

   d_tag_indx = (*s_tag_indx)[d_dim.getValue() - 1];
   d_buf_tag_indx = (*s_buf_tag_indx)[d_dim.getValue() - 1];

   if (d_hierarchy->getGridGeometry()->getNumberBlocks() > 1) {
      d_mb_tagger_strategy = new MultiblockGriddingTagger(d_dim);
      d_mb_tagger_strategy->setScratchTagPatchDataIndex(d_buf_tag_indx);
   }

   /*
    * Initialize communication algorithm for exchanging buffered tag data.
    */
   d_bdry_fill_tags.reset(new xfer::RefineAlgorithm(d_dim));

   d_bdry_fill_tags->
   registerRefine(d_buf_tag_indx,
      d_buf_tag_indx,
      d_buf_tag_indx,
      boost::shared_ptr<hier::RefineOperator>());

   d_efficiency_tolerance.resizeArray(d_hierarchy->getMaxNumberOfLevels());
   d_combine_efficiency.resizeArray(d_hierarchy->getMaxNumberOfLevels());

   for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ln++) {
      d_efficiency_tolerance[ln] = 0.8e0;
      d_combine_efficiency[ln] = 0.8e0;
   }

   d_proper_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels(),
      hier::BoxLevel(d_dim));
   d_to_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels());
   d_from_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels());

   /*
    * Initialize object with data read from input and restart databases.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart && d_registered_for_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);

   if (d_hierarchy->getMaxNumberOfLevels() > 1) {
      tbox::Array<hier::IntVector> ratio_to_coarser(d_hierarchy->getMaxNumberOfLevels(),
                                                    hier::IntVector::getOne(d_dim));
      for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
         ratio_to_coarser[ln] = d_hierarchy->getRatioToCoarserLevel(ln);
      }
      d_tag_init_strategy->checkCoarsenRatios(ratio_to_coarser);
      // Note: checkCoarsenRatios() must precede getErrorCoarsenRatio().
   }
   if (d_tag_init_strategy->getErrorCoarsenRatio() > 1) {
      boost::shared_ptr<StandardTagAndInitialize> std_tag_init(
         d_tag_init_strategy,
         boost::detail::dynamic_cast_tag());
      if (std_tag_init) {
         d_hierarchy->registerConnectorWidthRequestor(
            std_tag_init->getConnectorWidthRequestor());
      }
   }

   d_bdry_sched_tags.resizeArray(d_hierarchy->getMaxNumberOfLevels());

   warnIfDomainTooSmallInPeriodicDir();

   allocateTimers();

#ifdef GA_RECORD_STATS

   d_boxes_stat.resizeArray(d_hierarchy->getMaxNumberOfLevels());
   d_cells_stat.resizeArray(d_hierarchy->getMaxNumberOfLevels());
   d_timestamp_stat.resizeArray(d_hierarchy->getMaxNumberOfLevels());

   for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
      std::string ln_text = tbox::Utilities::intToString(ln, 2);
      const std::string num_boxes_str = std::string("GA_BoxesL") + ln_text;
      const std::string num_cells_str = std::string("GA_CellsL") + ln_text;
      const std::string timestamp_str = std::string("GA_TimeL") + ln_text;
      d_boxes_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(num_boxes_str, "PROC_STAT");
      d_cells_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(num_cells_str, "PROC_STAT");
      d_timestamp_stat[ln] = tbox::Statistician::getStatistician()->
         getStatistic(timestamp_str, "PROC_STAT");
   }

#endif

}

/*
 *************************************************************************
 *
 * Destructor tells the tbox::RestartManager to remove this object from
 * the list of restart items.
 *
 *************************************************************************
 */
GriddingAlgorithm::~GriddingAlgorithm()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
   delete d_mb_tagger_strategy;
   d_mb_tagger_strategy = NULL;
}

boost::shared_ptr<TagAndInitializeStrategy>
GriddingAlgorithm::getTagAndInitializeStrategy() const
{
   return d_tag_init_strategy;
}

boost::shared_ptr<LoadBalanceStrategy>
GriddingAlgorithm::getLoadBalanceStrategy() const
{
   return d_load_balancer;
}

boost::shared_ptr<LoadBalanceStrategy>
GriddingAlgorithm::getLoadBalanceStrategyZero() const
{
   return d_load_balancer0;
}

/*
 *************************************************************************
 *
 * Construct coarsest level in the patch hierarchy (i.e., level 0).
 * This routine operates in two modes:
 *
 * (1) If level 0 does not exist in the hierarchy, then a new level
 *  will be created based on the domain description provided by
 *  the grid geometry associated with the hierarchy.
 * (2) If level 0 exists in the hierarchy, then a new level is made
 *  by re-applying the load balancing routines to the existing level.
 *  The pre-existing level will be removed from the hierarchy.
 *
 * In either case, the new level is placed in the hierarchy as level 0.
 * If level 0 can be refined (i.e., it will be subject to tagging cells
 * for refinement), the domain boxes are checked against the constraints
 * of regridding process.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::makeCoarsestLevel(
   const double level_time)
{

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *d_hierarchy);

   if (d_barrier_and_time) {
      t_make_coarsest->barrierAndStart();
   }

   const hier::IntVector& zero_vec(hier::IntVector::getZero(d_dim));

   const hier::BoxLevelConnectorUtils dlbg_edge_utils;

   TBOX_ASSERT(d_hierarchy->getMaxNumberOfLevels() > 0);

   t_make_domain->start();

   const int ln = 0;

   bool level_zero_exists = false;
   if ((d_hierarchy->getNumberOfLevels() > 0)) {
      if (d_hierarchy->getPatchLevel(ln)) {
         level_zero_exists = true;
      }
   }

   hier::IntVector smallest_patch(d_dim);
   hier::IntVector largest_patch(d_dim);
   hier::IntVector extend_ghosts(d_dim);
   {
      hier::IntVector smallest_box_to_refine(d_dim);
      // "false" argument: for_building_finer level = false
      getGriddingParameters(
         smallest_patch,
         smallest_box_to_refine,
         largest_patch,
         extend_ghosts,
         ln,
         false);
   }

   /*
    * If making coarsest level for the first time, check domain boxes
    * for violations of user constraints.
    */
   if (!level_zero_exists) {
      for (int b = 0; b < d_hierarchy->getGridGeometry()->getNumberBlocks();
           b++) {
         hier::BoxContainer domain_boxes(
            d_hierarchy->getGridGeometry()->getPhysicalDomain(),
            hier::BlockId(b));
         checkDomainBoxes(domain_boxes);
      }
   }

   const hier::BoxLevel& domain_mapped_box_level(d_hierarchy->getDomainBoxLevel());

   t_make_domain->stop();

   /*
    * domain_to_domain is a temporary connector for bridging the
    * balanced BoxLevel to itself.  Since domain is small, owned
    * by just proc 0 and already globalized, computing
    * domain_to_domain is fast and does not require communication.
    */

   hier::Connector domain_to_domain(
      domain_mapped_box_level,
      domain_mapped_box_level,
      hier::IntVector::max(
         hier::IntVector::getOne(d_dim),
         d_hierarchy->getRequiredConnectorWidth(0, 0)));

   const hier::OverlapConnectorAlgorithm oca;
   oca.findOverlaps(domain_to_domain);

   if (d_barrier_and_time) {
      t_load_balance0->barrierAndStart();
   }

   t_load_balance_setup->start();

   // hier::IntVector patch_cut_factor(d_tag_init_strategy-> getErrorCoarsenRatio());
   const hier::IntVector patch_cut_factor(d_dim, 1);

   /*
    * FIXME: The code for generating the coarsest level's boxes is not
    * scalable because we use the domain BoxLevel in bridges and
    * Connector modifications, forcing proc 0 (which owns all domain
    * nodes) to do all of the searches.  This problem will show up in
    * performance data if we re-partition the coarsest level
    * regularly.
    */

   hier::BoxLevel new_mapped_box_level(domain_mapped_box_level);
   hier::Connector domain_to_new(domain_to_domain);
   domain_to_new.setHead(new_mapped_box_level, true);
   hier::Connector new_to_domain(domain_to_domain);
   new_to_domain.setBase(new_mapped_box_level, true);

   t_load_balance_setup->stop();

   d_load_balancer0->loadBalanceBoxLevel(
      new_mapped_box_level,
      new_to_domain,
      domain_to_new,
      d_hierarchy,
      ln,
      hier::Connector(),
      hier::Connector(),
      smallest_patch,
      largest_patch,
      domain_mapped_box_level,
      extend_ghosts,
      patch_cut_factor);

   if (d_barrier_and_time) {
      t_load_balance0->stop();
   }

   if (d_sequentialize_patch_indices) {
      renumberBoxes(new_mapped_box_level,
         domain_to_new,
         new_to_domain,
         false /* sort_by_corners */,
         true /* sequentialize_global_indices */);
   }

   dlbg_edge_utils.addPeriodicImagesAndRelationships(
      new_mapped_box_level,
      new_to_domain,
      domain_to_new,
      d_hierarchy->getGridGeometry()->getDomainSearchTree(),
      domain_to_domain);

   if (d_check_connectors) {
      TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_domain) == 0);
      TBOX_ASSERT(oca.checkOverlapCorrectness(domain_to_new) == 0);
   }

   hier::Connector* new_to_new = new hier::Connector;
   if (domain_mapped_box_level.getLocalNumberOfBoxes(0) ==
       (size_t)domain_mapped_box_level.getGlobalNumberOfBoxes()) {
      /*
       * If proc 0 owns all new boxes, it is faster find new<==>new by
       * globalizing the new boxes.
       *
       * The standard approach of bridging basically does the same,
       * but forces proc 0 to compute all the overlaps and send that
       * info to each processor one at a time.
       */

      if (d_barrier_and_time) {
         t_find_new_to_new->barrierAndStart();
      }

      new_to_new->setBase(new_mapped_box_level);
      new_to_new->setHead(new_mapped_box_level);
      new_to_new->setWidth(d_hierarchy->getRequiredConnectorWidth(0, 0), true);
      oca.findOverlaps(*new_to_new);

      if (d_barrier_and_time) {
         t_find_new_to_new->stop();
      }

   } else {

      if (d_barrier_and_time) {
         t_bridge_new_to_new->barrierAndStart();
      }

      oca.bridgeWithNesting(
         *new_to_new,
         *new_to_new,
         new_to_domain,
         domain_to_new,
         new_to_domain,
         domain_to_new,
         zero_vec,
         zero_vec,
         d_hierarchy->getRequiredConnectorWidth(0, 0));

      if (d_barrier_and_time) {
         t_bridge_new_to_new->stop();
      }

      TBOX_ASSERT(new_to_new->getConnectorWidth() ==
         d_hierarchy->getRequiredConnectorWidth(0, 0));
      TBOX_ASSERT(&new_to_new->getBase() == &new_mapped_box_level);
      TBOX_ASSERT(&new_to_new->getHead() == &new_mapped_box_level);

   }

   if (d_check_overlapping_patches != 'i') {
      checkOverlappingPatches(*new_to_new);
   }

   t_make_new->start();
   if (!level_zero_exists) {

      d_hierarchy->makeNewPatchLevel(ln,
         new_mapped_box_level);
      /*
       * Add computed Connectors to new level's collection of
       * persistent overlap Connectors.
       */
      boost::shared_ptr<hier::PatchLevel> new_level(
         d_hierarchy->getPatchLevel(ln));
      new_level->getBoxLevel()->getPersistentOverlapConnectors().
      cacheConnector(
         *new_level->getBoxLevel(),
         new_to_new);

      d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
         *d_hierarchy->getPatchLevel(0));

      // "true" argument: const bool initial_time = true;
      d_tag_init_strategy->initializeLevelData(d_hierarchy,
         ln,
         level_time,
         d_hierarchy->levelCanBeRefined(ln),
         true);

   } else {

      /*
       * Save old data before they are overwritten by the new BoxLevel.
       */
      hier::BoxLevel old_mapped_box_level = *d_hierarchy->getBoxLevel(0);

      boost::shared_ptr<hier::PatchLevel> old_level(
         d_hierarchy->getPatchLevel(ln));

      d_hierarchy->removePatchLevel(ln);

      d_hierarchy->makeNewPatchLevel(ln,
         new_mapped_box_level);

      /*
       * Compute old<==>new.  Doing it this way is not scalable, but
       * we only do this for the coarsest level.  The old approach of
       * bridging across the domain BoxLevel is probably not
       * scalable anyway, because the domain isusually owned by just
       * one processor.
       */
      old_mapped_box_level.getPersistentOverlapConnectors().createConnector(
         new_mapped_box_level,
         d_hierarchy->getRequiredConnectorWidth(0, 0));
      new_mapped_box_level.getPersistentOverlapConnectors().createConnector(
         old_mapped_box_level,
         d_hierarchy->getRequiredConnectorWidth(0, 0));

      d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
         *d_hierarchy->getPatchLevel(0));

      // "false" argument: const bool initial_time = false;
      d_tag_init_strategy->initializeLevelData(d_hierarchy,
         ln,
         level_time,
         d_hierarchy->levelCanBeRefined(ln),
         false,
         old_level);

      old_level.reset();

   }
   t_make_new->stop();

   if (d_barrier_and_time) {
      t_reset_hier->barrierAndStart();
   }
   d_tag_init_strategy->resetHierarchyConfiguration(d_hierarchy, ln, ln);
   if (d_barrier_and_time) {
      t_reset_hier->stop();
   }

   if (d_log_metadata_statistics) {
      logMetadataStatistics("makeCoarsestLevel", 0, false, false);
   }

#ifdef GA_RECORD_STATS
   recordStatistics(level_time);
#endif

   if (d_barrier_and_time) {
      t_make_coarsest->stop();
   }
}

/*
 *************************************************************************
 *
 * Perform error estimation process on the finest hierarchy level to
 * determine if and where a new finest level is needed.  If it is
 * determined  that a new finest level should exist, it is created and
 * its patch data is allocated and initialized.  The algorithm is
 * summarized as follows:
 *
 * (1) Compute boxes whose union constitutes the region within the level
 *  in which the next finer level may reside (i.e., proper nesting).
 *
 * (2) Set tags on the level to indicate which cells should be refined.
 *
 * (3) Buffer the tags.  This prevents disturbances from moving off
 *  refined regions before the next remesh occurs.
 *
 * (4) Determine boxes for new patches that will cover the tagged cells.
 *  This step includes load balancing of the patches.
 *
 * (5) If there exist some regions to refine, construct the next finer
 *  level and insert it in the hierarchy.  Then, initialize the data
 *  on that level.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::makeFinerLevel(
   const double level_time,
   const bool initial_time,
   const int tag_buffer,
   const double regrid_start_time)
{
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::makeFinerLevel entered with finest ln = "
      << d_hierarchy->getFinestLevelNumber() << "\n";
   }

   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT(d_hierarchy->getPatchLevel(
                    d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(tag_buffer >= 0);

   if (d_barrier_and_time) {
      t_make_finer->barrierAndStart();
   }

   const hier::OverlapConnectorAlgorithm oca;

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   const hier::BoxLevelConnectorUtils dlbg_edge_utils;

   const int tag_ln = d_hierarchy->getFinestLevelNumber();
   const int new_ln = tag_ln + 1;

   if (d_hierarchy->levelCanBeRefined(tag_ln)) {

      t_make_finer_setup->start();

      /*
       * d_base_ln is used by private methods during regrid.
       */
      d_base_ln = tag_ln;

      /*
       * Compute nesting data at d_base_ln for use in constructing
       * level d_base_ln+1;
       */
      computeProperNestingData(d_base_ln);

      const boost::shared_ptr<hier::PatchLevel> tag_level(
         d_hierarchy->getPatchLevel(tag_ln));

      hier::BoxLevel new_mapped_box_level(d_dim);
      hier::Connector* tag_to_new = new hier::Connector;
      hier::Connector* new_to_tag = new hier::Connector;
      hier::Connector* new_to_new = new hier::Connector;

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *   1) only user supplied refine boxes are used
       *   2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       */
      bool do_tagging = true;
      if (d_tag_init_strategy->refineUserBoxInputOnly()) do_tagging = false;
      t_make_finer_setup->stop();

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         t_make_finer_tagging->start();

         /*
          * Initialize integer tag arrays on level to false.
          */
         tag_level->allocatePatchData(d_tag_indx);
         fillTags(d_false_tag, tag_level, d_tag_indx);

         /*
          * Perform pre-processing of error estimation data.
          */
         d_tag_init_strategy->
         preprocessErrorEstimation(d_hierarchy,
            tag_ln,
            level_time,
            regrid_start_time,
            initial_time);

         /*
          * Determine cells needing refinement on level and set tags to true.
          * Because we are constructing a new level, not regridding the level,
          * the coarsest_sync_level argument is always false.
          */
         bool coarsest_sync_level = false;
         d_tag_init_strategy->
         tagCellsForRefinement(d_hierarchy,
            tag_ln,
            level_time,
            d_tag_indx,
            initial_time,
            coarsest_sync_level,
            d_hierarchy->levelCanBeRefined(tag_ln),
            regrid_start_time);

         /*
          * Check for user-tagged cells that violate proper nesting.
          * except if user specified that the violating tags be ignored.
          */
         if (d_check_nonrefined_tags != 'i') {
            checkNonrefinedTags(*tag_level, tag_ln);
         }

         /*
          * Buffer true tagged cells by specified amount which should be
          * sufficient to keep disturbance on refined region until next regrid
          * of the level occurs.
          */
         hier::IntVector max_descriptor_ghosts(
            d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_dim));

         /*
          * If the tag buffer value passed into this method is greater than the
          * current ghost width of the data that handles tag buffering, then
          * the call to resetTagBufferingData resets that data to have a ghost
          * width equal to the tag buffer.
          */

         if (tag_buffer > d_buf_tag_ghosts.max()) {
            resetTagBufferingData(tag_buffer);
         }

         /*
          * Create communication schedule for buffer tags on this level.
          */
         t_bdry_fill_tags_create->start();
         d_bdry_sched_tags[tag_ln] =
            d_bdry_fill_tags->createSchedule(tag_level, d_mb_tagger_strategy);
         t_bdry_fill_tags_create->stop();


         tag_level->allocatePatchData(d_buf_tag_indx);
         bufferTagsOnLevel(d_true_tag, tag_level, tag_buffer);
         tag_level->deallocatePatchData(d_buf_tag_indx);
         t_make_finer_tagging->stop();

         /*
          * We cannot leave this method with the tag buffering data having
          * ghosts greater than any other data managed by the patch descriptor,
          * so if that is the case, we reset it to the default value of 1.
          */

         if (tag_buffer > max_descriptor_ghosts.max()) {
            resetTagBufferingData(1);
         }

         /*
          * Determine Boxes for new fine level.
          */
         findRefinementBoxes(new_mapped_box_level,
            *tag_to_new,
            *new_to_tag,
            tag_ln);

         if (new_mapped_box_level.isInitialized()) {

            if (d_check_proper_nesting) {
               /*
                * Check that the new BoxLevel nests inside the tag level.
                *
                * SAMRAI convention (my best understanding of it) supported
                * (or should have been supported) by grid generation:
                * - L0 must be equivalent to the domain.
                * - L1 must nest in L0 by the max ghost width in L1 index space,
                *   except where L1 touches physical boundary.
                * - L(n) must nest in L(n-1) by getProperNestingBuffer(n-1),
                *   except where L(n) touches physical boundary.
                */

               hier::IntVector required_nesting(d_dim);
               if (tag_ln > 0) {
                  required_nesting = hier::IntVector(d_dim,
                        d_hierarchy->getProperNestingBuffer(tag_ln));
               } else {
                  required_nesting = max_descriptor_ghosts;
               }

               bool locally_nests = false;
               const bool new_nests_in_tag =
                  dlbg_edge_utils.baseNestsInHead(
                     &locally_nests,
                     new_mapped_box_level,
                     *d_hierarchy->getBoxLevel(tag_ln),
                     required_nesting,
                     zero_vector,
                     zero_vector,
                     &d_hierarchy->getGridGeometry()->getPeriodicDomainSearchTree());

               if (!new_nests_in_tag) {
                  hier::BoxLevel violator(d_dim);
                  hier::Connector new_to_violator;
                  hier::BoxLevelConnectorUtils edge_utils;
                  t_compute_external_parts->start();
                  edge_utils.computeExternalParts(
                     violator,
                     new_to_violator,
                     *new_to_tag,
                     hier::IntVector(d_dim, -d_hierarchy->getProperNestingBuffer(tag_ln)),
                     d_hierarchy->getGridGeometry()->getDomainSearchTree());
                  t_compute_external_parts->stop();

                  TBOX_ERROR(
                     "Internal library error: Failed to produce proper nesting."
                     << "GriddingAlgorithm::makeFinerLevel:\n"
                     << "tag_ln=" << tag_ln << ":\n"
                     << "new mapped_box_level does not properly nest\n"
                     << "in tag mapped_box_level by the required nesting buffer of "
                     << d_hierarchy->getProperNestingBuffer(tag_ln)
                     << ".\nLocal nestingness: " << locally_nests
                     << ".\nProper nesting violation with new_mapped_box_level of\n"
                     << new_mapped_box_level.format("N->", 2)
                     << "Proper nesting violation with tag mapped_box_level of\n"
                     << d_hierarchy->getBoxLevel(tag_ln)->format("F->", 2)
                     << "tag_to_new:\n" << tag_to_new->format("N->", 2)
                     << "new_to_tag:\n" << new_to_tag->format("N->", 2)
                     << "violator:\n" << violator.format("N->", 2)
                     << "new_to_violator:\n" << new_to_violator.format("N->", 2));
               }
            }

         }

         /*
          * Deallocate tag arrays and schedule -- no longer needed.
          */
         tag_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[tag_ln].reset();

      } else { /* do_tagging == false */

         /*
          * If tagging is not necessary (do_tagging = false) we simply
          * need to access the level boxes, either from a dumped file or
          * from user-supplied refine boxes, and load balance them before
          * constructing the finer level.
          */
         bool remove_old_fine_level = false;
         readLevelBoxes(new_mapped_box_level,
            *tag_to_new,
            *new_to_tag,
            tag_ln,
            level_time,
            remove_old_fine_level);

         /*
          * Check for user-specified boxes that violate nesting requirements.
          */

         if (d_check_nonnesting_user_boxes != 'i') {
            checkNonnestingUserBoxes(
               *new_to_tag,
               new_to_tag->getRatio() * d_hierarchy->getProperNestingBuffer(tag_ln));
         }

         if (d_check_boundary_proximity_violation != 'i') {
            checkBoundaryProximityViolation(tag_ln, new_mapped_box_level);
         }

      } /* do_tagging == false */

      if (new_mapped_box_level.isInitialized()) {

         /*
          * If we made a new_mapped_box_level, proceed to make a
          * PatchLevel from it.
          */

         // Bridge for new<==>new.
         t_bridge_links->start();
         t_bridge_new_to_new->start();
         oca.bridgeWithNesting(
            *new_to_new,
            *new_to_new,
            *new_to_tag,
            *tag_to_new,
            *new_to_tag,
            *tag_to_new,
            zero_vector,
            zero_vector,
            d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));
         t_bridge_new_to_new->stop();
         t_bridge_links->stop();

         if (d_check_overlapping_patches != 'i') {
            checkOverlappingPatches(*new_to_new);
         }

         t_make_finer_create->start();

         d_hierarchy->makeNewPatchLevel(new_ln, new_mapped_box_level);

         boost::shared_ptr<hier::PatchLevel> new_level(
            d_hierarchy->getPatchLevel(new_ln));

         new_level->getBoxLevel()->getPersistentOverlapConnectors().
         cacheConnector(
            *new_level->getBoxLevel(),
            new_to_new);

         new_level->getBoxLevel()->getPersistentOverlapConnectors().
         cacheConnector(
            *tag_level->getBoxLevel(),
            new_to_tag);

         tag_level->getBoxLevel()->getPersistentOverlapConnectors().
         cacheConnector(
            *new_level->getBoxLevel(),
            tag_to_new);

         d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
            *d_hierarchy->getPatchLevel(new_ln));

         d_tag_init_strategy->initializeLevelData(d_hierarchy,
            new_ln,
            level_time,
            d_hierarchy->levelCanBeRefined(new_ln),
            initial_time);

         t_reset_hier->barrierAndStart();
         d_tag_init_strategy->resetHierarchyConfiguration(d_hierarchy,
            new_ln,
            new_ln);
         t_reset_hier->stop();
         t_make_finer_create->stop();

         if (d_log_metadata_statistics) {
            logMetadataStatistics("makeFinerLevel", d_hierarchy->getFinestLevelNumber(), false, true);
         }
      }
      else {
         delete tag_to_new;
         delete new_to_tag;
         delete new_to_new;
      }

      d_base_ln = -1;

   }  // if level cannot be refined, the routine drops through...

   if (d_barrier_and_time) {
      t_make_finer->stop();
   }

}

/*
 *************************************************************************
 *
 * Regrid each level in the hierarchy which is finer than the
 * specified level.  If the regridding procedure employs time
 * integration, we perform any pre-processing necessary to regrid the
 * levels.  Then, each level finer than the specified level is
 * regridded from fine to coarse.  The recursive regridding procedure
 * is performed by the function regridFinerLevel().  Finally, after
 * the new hierarchy configuration is set, the application-specific
 * operations for resetting hierarchy-dependent infomation is called.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::regridAllFinerLevels(
   const int level_number,
   const double regrid_time,
   const tbox::Array<int>& tag_buffer,
   const tbox::Array<double> regrid_start_time,
   const bool level_is_coarsest_sync_level)
{
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(d_hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(tag_buffer.getSize() >= level_number + 1);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < tag_buffer.getSize(); i++) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif

   if (d_barrier_and_time) {
      t_regrid_all_finer->barrierAndStart();
   }
   t_misc5->start();

   if (d_hierarchy->levelCanBeRefined(level_number)) {

      if (d_print_steps) {
         tbox::plog
         << "GriddingAlgorithm::regridAllFinerLevels: regridding finer than "
         << level_number << std::endl;
      }

      /*
       * d_base_ln is used by private methods during regrid.
       */
      d_base_ln = level_number;

      t_process_error->start();
      /*
       * Perform pre-processing of error estimation data, if
       * appropriate.
       */
      if (d_tag_init_strategy->usesTimeIntegration()) {
         for (int ln = level_number;
              ln <= d_hierarchy->getFinestLevelNumber(); ln++) {
            if (d_hierarchy->levelCanBeRefined(ln)) {
               bool initial_time = false;
               double level_regrid_start_time = 0.;
               if (regrid_start_time.getSize() < ln + 1) {
                  tbox::IEEE::setNaN(level_regrid_start_time);
               } else {
                  level_regrid_start_time = regrid_start_time[ln];
               }

               d_tag_init_strategy->
               preprocessErrorEstimation(d_hierarchy,
                  ln,
                  regrid_time,
                  level_regrid_start_time,
                  initial_time);
            }
         }
      }
      t_process_error->stop();

      t_misc5->stop();
      /*
       * Recursively regrid each finer level.
       */
      const int finest_level_not_regridded = level_number;
      regridFinerLevel(
         level_number,
         regrid_time,
         finest_level_not_regridded,
         level_is_coarsest_sync_level,
         tag_buffer,
         regrid_start_time);

      t_misc5->start();

      /*
       * Invoke application-specific routines to reset information for those
       * levels which have been modified.
       */

      if (d_hierarchy->getFinestLevelNumber() >= (level_number + 1)) {
         if (d_barrier_and_time) {
            t_reset_hier->barrierAndStart();
         }
         d_tag_init_strategy->
         resetHierarchyConfiguration(d_hierarchy,
            level_number + 1,
            d_hierarchy->getFinestLevelNumber());
         if (d_barrier_and_time) {
            t_reset_hier->stop();
         }
      }

      d_base_ln = -1;

      if (d_print_steps) {
         tbox::plog
         << "GriddingAlgorithm::regridAllFinerLevels: regridded finer than "
         << level_number << std::endl;
      }

   } //  if level cannot be refined, the routine drops through...

   else {
      if (d_print_steps) {
         tbox::plog
         << "GriddingAlgorithm::regridAllFinerLevels: level "
         << level_number << " cannot be refined." << std::endl;
      }
   }

#ifdef GA_RECORD_STATS
   // Verified that this does not use much time.
   recordStatistics(regrid_time);
#endif

   t_misc5->stop();

   if (d_barrier_and_time) {
      t_regrid_all_finer->stop();
   }

}

/*
 *************************************************************************
 *
 * Recursively, regrid each AMR hierarchy level finer than the
 * specified level (indexed by tag_ln).  The process is as follows:
 *
 * (1) Initialize tags to false on the level.
 *
 * (2) If a finer level exists, set tag to true on level for each cell
 * that is refined.
 *
 * (3) Tag cells for refinement on level by applying application-
 * specific error estimation routines.
 *
 * (4) If a finer level exists, invoke process recursively (i.e.,
 * invoke step 1 on next finer level).
 *
 * (5) (Note we have popped out of recursion at this point).  Buffer
 * true tags on current level to keep disturbances on fine grids until
 * regridding occurs next.
 *
 * (6) If level tag_ln+2 exists, tag under its boxes to ensure level
 * tag_ln+1 properly nests it.
 *
 * (7) Determine box configuration for new finer level, by calling
 * findRefinementBoxes() function.
 *
 * (8) If a finer level should exist in the hierarchy, create its
 * patches from the box description and initialize its data.  If
 * necessary, discard old level.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::regridFinerLevel(
   const int tag_ln,
   const double regrid_time,
   const int finest_level_not_regridded,
   const bool level_is_coarsest_sync_level,
   const tbox::Array<int>& tag_buffer,
   const tbox::Array<double>& regrid_start_time)
{
   TBOX_ASSERT((tag_ln >= 0)
      && (tag_ln <= d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(d_hierarchy->getPatchLevel(tag_ln));
   TBOX_ASSERT(finest_level_not_regridded >= 0
      && finest_level_not_regridded <= tag_ln);
   TBOX_ASSERT(tag_buffer.getSize() >= tag_ln + 1);
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < tag_buffer.getSize(); i++) {
      TBOX_ASSERT(tag_buffer[i] >= 0);
   }
#endif

   const hier::OverlapConnectorAlgorithm oca;

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel: entered with tag_ln = "
      << tag_ln << "\n";
   }

   hier::BoxLevelConnectorUtils dlbg_edge_utils;

   if (d_hierarchy->levelCanBeRefined(tag_ln)) {

      if (d_barrier_and_time) {
         t_misc4->start();
      }

      int new_ln = tag_ln + 1;

      boost::shared_ptr<hier::PatchLevel> tag_level(
         d_hierarchy->getPatchLevel(tag_ln));

      /*
       * Compute nesting data at tag_ln for use in constructing
       * level tag_ln+1;
       */
      computeProperNestingData(tag_ln);

      /*
       * The boolean "do_tagging" specifies whether or not tagging will
       * be performed. This will be true except in two circumstances:
       *   1) only user supplied refine boxes are used
       *   2) the boxes are read from a previously dumped file.
       *
       * If either of these circumstances is true, tagging operations
       * are NOT necessary so do_tagging will be set to false.
       *
       * The old level is generally removed when regridding, but
       * some circumstances may warrant keeping the old level.  For
       * example, if the refine region has not changed, there is no
       * need to regenerate the finer level.  The boolean
       * "remove_old_fine_level" specifies if the old level should
       * be removed.
       */
      bool do_tagging = true;
      if (d_tag_init_strategy->refineUserBoxInputOnly()) do_tagging = false;
      bool remove_old_fine_level = true;

      hier::BoxLevel new_mapped_box_level(d_dim);
      hier::Connector* tag_to_new = new hier::Connector;
      hier::Connector* new_to_tag = new hier::Connector;

      if (d_barrier_and_time) {
         t_misc4->stop();
      }

      /*
       * tag_to_finer is [tag_ln]->[tag_ln+2].
       * finer_to_tag is [tag_ln+2]->[tag_ln].
       *
       * These are declared in this scope, computed and cached if
       * [tag_ln+2] exists.  They are used to tag the footprint of
       * [tag_ln+2] on level tag_ln.  Later on, they are used to
       * bridge for [new_ln] <-> [tag_ln+2].
       */
      hier::Connector tag_to_finer, finer_to_tag;

      /*
       * Tag cells, determine refine boxes from tagged regions, and
       * load balance in preparation for constructing new refined level.
       */
      if (do_tagging) {

         /*
          * Tagging stuff have been factored out to shorten this method.
          */
         regridFinerLevel_doTaggingBeforeRecursiveRegrid(
            tag_ln,
            level_is_coarsest_sync_level,
            regrid_start_time,
            regrid_time);

         /*
          * Perform regridding recursively on finer levels, if appropriate.
          */
         if (d_hierarchy->finerLevelExists(tag_ln)
             && d_hierarchy->levelCanBeRefined(new_ln)) {

            if (d_print_steps) {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel: recursing to tag_ln = "
               << new_ln << "\n";
            }

            regridFinerLevel(
               new_ln,
               regrid_time,
               finest_level_not_regridded,
               false,
               tag_buffer,
               regrid_start_time);

            if (d_print_steps) {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel: recursion returned to tag_ln = "
               << tag_ln << "\n";
            }

         }

         /*
          * Tagging stuff have been factored out to shorten this method.
          */
         regridFinerLevel_doTaggingAfterRecursiveRegrid(
            tag_to_finer,
            finer_to_tag,
            tag_ln,
            tag_buffer);

         /*
          * Determine boxes containing cells on level with a true tag
          * value.
          */
         findRefinementBoxes(
            new_mapped_box_level,
            *tag_to_new,
            *new_to_tag,
            tag_ln);

         if (d_barrier_and_time) {
            t_misc3->start();
         }

         if (d_print_steps) {
            if (new_mapped_box_level.isInitialized()) {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel got inititalized new_mapped_box_level\n";
            } else {
               tbox::plog
               << "GriddingAlgorithm::regridFinerLevel got un-inititalized new_mapped_box_level\n";
            }
         }

         /*
          * Deallocate tag arrays and schedule; no longer needed on current
          * level.
          */

         tag_level->deallocatePatchData(d_tag_indx);
         d_bdry_sched_tags[tag_ln].reset();

         if (d_barrier_and_time) {
            t_misc3->stop();
         }

      } else { /* do_tagging == false */

         if (d_hierarchy->finerLevelExists(tag_ln)
             && d_hierarchy->levelCanBeRefined(new_ln)) {
            regridFinerLevel(
               new_ln,
               regrid_time,
               finest_level_not_regridded,
               false,
               tag_buffer,
               regrid_start_time);
         }

         /*
          * If tagging is not necessary (do_tagging = false) we simply
          * need to access the user-supplied refine boxes, and load
          * balance them before constructing the finer level.
          */
         readLevelBoxes(new_mapped_box_level,
            *tag_to_new,
            *new_to_tag,
            tag_ln,
            regrid_time,
            remove_old_fine_level);

      } /* end do_tagging == false */

      /*
       * Make new finer level (new_ln) if necessary, or remove
       * next finer level if it is no longer needed.
       */

      if (new_mapped_box_level.isInitialized()) {

         /*
          * Create the new PatchLevel from the new_mapped_box_level.
          */
         regridFinerLevel_createAndInstallNewLevel(
            tag_ln,
            regrid_time,
            tag_to_new,
            new_to_tag,
            tag_to_finer,
            finer_to_tag);

         if (d_log_metadata_statistics) {
            // Don't log the coarse Connector, if the coarse level will be updated.
            logMetadataStatistics("regridFinerLevel", new_ln, new_ln<d_hierarchy->getFinestLevelNumber(), tag_ln==d_base_ln);
            tbox::plog << "GriddingAlgorithm::regridFinerLevel finished logging level stats." << std::endl;
         }

      } else {

         /*
          * The new level has no boxes, so we don't generate it.
          * Remove the pre-existing fine level if it exists.
          */

         if (d_hierarchy->finerLevelExists(tag_ln)
             && remove_old_fine_level) {
            d_hierarchy->removePatchLevel(new_ln);
         }

         delete tag_to_new;
         delete new_to_tag;

      } // if we are not re-regenerating level new_ln.

   } //  if level cannot be refined, the routine drops through...

}

/*
 *************************************************************************
 * Various tagging stuff done before recursively regridding a finer level.
 *************************************************************************
 */
void
GriddingAlgorithm::regridFinerLevel_doTaggingBeforeRecursiveRegrid(
   const int tag_ln,
   const bool level_is_coarsest_sync_level,
   const tbox::Array<double>& regrid_start_time,
   const double regrid_time)
{
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel_doTaggingBeforeRecursiveRegrid: entered with tag_ln = "
      << tag_ln << "\n";
   }

   const hier::IntVector& zero_vec(hier::IntVector::getZero(d_dim));
   const boost::shared_ptr<hier::PatchLevel>& tag_level(
      d_hierarchy->getPatchLevel(tag_ln));

   /*
    * Create communication schedule for buffer tags and set tags to
    * false.
    */

   if (d_barrier_and_time) {
      t_misc1->start();
   }

   tag_level->allocatePatchData(d_tag_indx);
   fillTags(d_false_tag, tag_level, d_tag_indx);

   /*
    * Set tags to true for cells that currently cover next finer level.
    * Note that this is not needed for all regridding strategies.  But
    * knowledge of currently refined cells is generally useful to avoid
    * repeated refining and coarsening of cells near boundaries of
    * refined regions where error estimates may hover around the error
    * tolerance. For example, the regridding scheme may require that
    * the error in a cell that is currently refined fall below a
    * certain tolerance  (generally different than the tolerance to
    * refine the cell in the first place) before the cell will be
    * de-refined.
    */

   if (d_hierarchy->finerLevelExists(tag_ln)) {
      fillTagsFromBoxLevel(
         d_true_tag,
         tag_level,
         d_tag_indx,
         d_hierarchy->getConnector(tag_ln, tag_ln + 1),
         true,
         zero_vec);
   }

   /*
    * Determine cells needing refinement according to a specific
    * error estimation procedure and set to true.
    *
    * The "level_is_coarsest_sync_level" is provided as an argument
    * to this method.  Provide the additional check of whether the
    * level is not the coarsest level and that it is not a new level
    * in the hierarchy.  If all three conditions are true, the
    * "coarsest_sync_level" argument passed into the tagCells method
    * will be true.  Otherwise, it will be false.
    */

   bool coarsest_sync_level =
      level_is_coarsest_sync_level &&
      tag_ln > 0 &&
      tag_ln <= d_base_ln;

   bool initial_time = false;
   double level_regrid_start_time = 0.;
   if (regrid_start_time.getSize() < tag_ln + 1) {
      tbox::IEEE::setNaN(level_regrid_start_time);
   } else {
      level_regrid_start_time = regrid_start_time[tag_ln];
   }

   if (d_barrier_and_time) {
      t_misc1->stop();
   }

   if (d_barrier_and_time) {
      t_tag_cells_for_refinement->barrierAndStart();
   }
   d_tag_init_strategy->
   tagCellsForRefinement(d_hierarchy,
      tag_ln,
      regrid_time,
      d_tag_indx,
      initial_time,
      coarsest_sync_level,
      d_hierarchy->levelCanBeRefined(tag_ln),
      level_regrid_start_time);
   if (d_barrier_and_time) {
      t_tag_cells_for_refinement->stop();
   }

   /*
    * Check for user-tagged cells that violate proper nesting,
    * except if user specified that the violating tags be ignored.
    */
   if (d_check_nonrefined_tags != 'i') {
      checkNonrefinedTags(*tag_level, tag_ln);
   }
}

/*
 *************************************************************************
 * Various tagging stuff done after recursively regridding a finer level.
 *************************************************************************
 */
void
GriddingAlgorithm::regridFinerLevel_doTaggingAfterRecursiveRegrid(
   hier::Connector& tag_to_finer,
   hier::Connector& finer_to_tag,
   const int tag_ln,
   const tbox::Array<int>& tag_buffer)
{
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel_doTaggingAfterRecursiveRegrid: entered with tag_ln = "
      << tag_ln << "\n";
   }

   const int new_ln = tag_ln + 1;
   const boost::shared_ptr<hier::PatchLevel>& tag_level(
      d_hierarchy->getPatchLevel(tag_ln));
   const hier::OverlapConnectorAlgorithm oca;
   const hier::BoxLevelConnectorUtils mblc_utils;

   if (d_barrier_and_time) {
      t_misc2->barrierAndStart();
   }

   /*
    * Buffer true tagged cells by specified amount which should be
    * sufficient to keep disturbance on refined region until next
    * regrid of the level occurs.
    */

   hier::IntVector max_descriptor_ghosts(
      d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_dim));

   /*
    * If the tag buffer value passed into this method is greater than the
    * current ghost width of the data that handles tag buffering, then the
    * call to resetTagBufferingData resets that data to have a ghost width
    * equal to the tag buffer.
    */ 
   if (tag_buffer[tag_ln] > d_buf_tag_ghosts.max()) {
      resetTagBufferingData(tag_buffer[tag_ln]);
   }

   t_bdry_fill_tags_create->start();
   d_bdry_sched_tags[tag_ln] =
      d_bdry_fill_tags->createSchedule(tag_level, d_mb_tagger_strategy);
   t_bdry_fill_tags_create->stop();

   tag_level->allocatePatchData(d_buf_tag_indx);
   bufferTagsOnLevel(d_true_tag, tag_level, tag_buffer[tag_ln]);

   /*
    * We cannot leave this method with the tag buffering data having ghosts
    * greater than any other data managed by the patch descriptor, so if that
    * is the case, we reset it to the default value of 1.
    */
   if (tag_buffer[tag_ln] > max_descriptor_ghosts.max()) {
      resetTagBufferingData(1);
   }

   if (d_hierarchy->finerLevelExists(new_ln)) {

      /*
       * Add tags to ensure new_ln properly nests level
       * new_ln+1.  This means tagging under level new_ln+1
       * plus an appropriate nesting buffer equal to how
       * much new_ln+1 should nest in new_ln.
       *
       * To determine where to tag, we compute Connector
       * [tag_ln] ---> [new_ln+1], aka tag_to_finer.  Then use
       * its neighbor data to determine where to tag.
       */

      if (d_barrier_and_time) {
         t_second_finer_tagging->start();
      }

      const hier::Connector& tag_to_old = d_hierarchy->getConnector(tag_ln, new_ln);
      const hier::Connector& old_to_finer = d_hierarchy->getConnector(new_ln, new_ln + 1);
      const hier::Connector& finer_to_old = d_hierarchy->getConnector(new_ln + 1, new_ln);
      const hier::Connector& old_to_tag = d_hierarchy->getConnector(new_ln, tag_ln);
      oca.bridge(
         tag_to_finer,
         finer_to_tag,
         tag_to_old,
         old_to_finer,
         finer_to_old,
         old_to_tag);

      // Nesting buffer in resolution of level new_ln+1.
      const hier::IntVector nesting_buffer =
         d_hierarchy->getRatioToCoarserLevel(new_ln + 1)
         * d_hierarchy->getProperNestingBuffer(tag_ln + 1);

#ifdef DEBUG_CHECK_ASSERTIONS
      oca.assertOverlapCorrectness(tag_to_finer, false, true, true);
      oca.assertOverlapCorrectness(finer_to_tag, false, true, true);
      TBOX_ASSERT(
         tag_to_finer.getConnectorWidth()
         * d_hierarchy->getRatioToCoarserLevel(tag_ln + 1)
         * d_hierarchy->getRatioToCoarserLevel(new_ln + 1) >= nesting_buffer);
#endif

      // Add periodic relationships in tag_to_finer.
      hier::BoxLevel dummy_finer_mapped_box_level = finer_to_old.getBase();
      mblc_utils.addPeriodicImagesAndRelationships(
         dummy_finer_mapped_box_level,
         finer_to_tag,
         tag_to_finer,
         d_hierarchy->getGridGeometry()->getDomainSearchTree(),
         d_hierarchy->getConnector(tag_ln, tag_ln));

#ifdef DEBUG_CHECK_ASSERTIONS
      oca.assertOverlapCorrectness(tag_to_finer, false, true, false);
      oca.assertOverlapCorrectness(finer_to_tag, false, true, false);
#endif

      fillTagsFromBoxLevel(
         d_true_tag,
         tag_level,
         d_tag_indx,
         tag_to_finer,
         true,
         nesting_buffer);

      if (d_barrier_and_time) {
         t_second_finer_tagging->stop();
      }

   } // End tagging under level new_ln+1.

   tag_level->deallocatePatchData(d_buf_tag_indx);

   if (d_barrier_and_time) {
      t_misc2->stop();
   }
}

/*
 *************************************************************************
 * After the metadata describing the new level is computed,
 * this method creates and installs new PatchLevel in the hierarchy.
 *************************************************************************
 */
void
GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel(
   const int tag_ln,
   const double regrid_time,
   hier::Connector* tag_to_new,
   hier::Connector* new_to_tag,
   const hier::Connector& tag_to_finer,
   const hier::Connector& finer_to_tag)
{
   /*
    * Compute self-Connector for the new level.
    */
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel: bridge links\n";
   }

   const int new_ln = tag_ln + 1;
   const boost::shared_ptr<hier::PatchLevel>& tag_level(
      d_hierarchy->getPatchLevel(tag_ln));
   const hier::BoxLevel& new_mapped_box_level(tag_to_new->getHead());

   const hier::OverlapConnectorAlgorithm oca;
   const hier::BoxLevelConnectorUtils mblc_utils;
   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   hier::Connector* new_to_new = new hier::Connector;

   t_bridge_links->start();
   t_bridge_new_to_new->start();

   oca.bridgeWithNesting(
      *new_to_new,
      *new_to_new,
      *new_to_tag,
      *tag_to_new,
      *new_to_tag,
      *tag_to_new,
      zero_vector,
      zero_vector,
      d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));
   new_to_new->setConnectorType(hier::Connector::COMPLETE_OVERLAP);

   t_bridge_new_to_new->stop();
   t_bridge_links->stop();

   TBOX_ASSERT(new_to_new->getConnectorWidth() ==
      d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));

   if (d_check_overlapping_patches != 'i') {
      checkOverlappingPatches(*new_to_new);
   }

   if (d_barrier_and_time) {
      t_regrid_finer_create->barrierAndStart();
   }

   /*
    * Either remove pre-existing fine level from hierarchy and make
    * a new level, or just make a new fine level for hierarchy.
    */

   /*
    * Save references to old objects before hierarchy removes them.
    * We need this while installing new objects.
    */

   boost::shared_ptr<hier::PatchLevel> old_fine_level;
   boost::shared_ptr<const hier::BoxLevel> old_mapped_box_level;
   const hier::Connector* old_to_tag = NULL;
   const hier::Connector* tag_to_old = NULL;

   hier::IntVector ratio(tag_level->getRatioToLevelZero()
                         * d_hierarchy->getRatioToCoarserLevel(new_ln));

   if (d_hierarchy->finerLevelExists(tag_ln)) {

      old_mapped_box_level = d_hierarchy->getBoxLevel(new_ln);
      old_to_tag = &d_hierarchy->getConnector(new_ln, tag_ln);
      tag_to_old = &d_hierarchy->getConnector(tag_ln, new_ln);

      old_fine_level = d_hierarchy->getPatchLevel(new_ln);
      d_hierarchy->removePatchLevel(new_ln);
      // ratio = old_fine_level->getRatioToLevelZero();
      TBOX_ASSERT(ratio == old_fine_level->getRatioToLevelZero());

   }

   if (d_hierarchy->levelExists(new_ln + 1)) {

      if (d_check_proper_nesting) {

         /*
          * Check that the new_mapped_box_level nests the next finer
          * level (new_ln+1).
          */

         hier::IntVector required_nesting(d_dim, d_hierarchy->getProperNestingBuffer(new_ln));
         required_nesting *= d_hierarchy->getRatioToCoarserLevel(new_ln + 1);

         bool locally_nests = false;
         const bool finer_nests_in_new =
            mblc_utils.baseNestsInHead(
               &locally_nests,
               *d_hierarchy->getBoxLevel(new_ln + 1),
               new_mapped_box_level,
               required_nesting,
               zero_vector,
               zero_vector,
               &d_hierarchy->getGridGeometry()->getPeriodicDomainSearchTree());

         if (!finer_nests_in_new) {

            tbox::perr << "GriddingAlgorithm: new mapped_box_level\n"
                       << "at ln=" << new_ln
                       << " does not properly nest\n"
                       << "existing finer mapped_box_level at ln="
                       << new_ln + 1
                       << " by the required nesting buffer of "
                       << required_nesting << " in fine resolution.\n"
                       << "Local nestingness: " << locally_nests
                       << ".\nWriting BoxLevels out to log file."
                       << std::endl;
            tbox::plog
            << "Proper nesting violation with new_mapped_box_level of\n"
            << new_mapped_box_level.format("N->", 2)
            << "Proper nesting violation with finer mapped_box_level of\n"
            << d_hierarchy->getBoxLevel(new_ln + 1)->format("F->", 2);

            hier::BoxLevel external(d_dim);
            hier::Connector finer_to_external;
            hier::Connector finer_to_new(
               *d_hierarchy->getBoxLevel(new_ln + 1),
               new_mapped_box_level,
               required_nesting);
            const hier::OverlapConnectorAlgorithm oca;
            oca.findOverlaps(finer_to_new);
            tbox::plog << "Finer to new:\n" << finer_to_new.format("FN->", 3);
            mblc_utils.computeExternalParts(
               external,
               finer_to_external,
               finer_to_new,
               -required_nesting,
               d_hierarchy->getGridGeometry()->getDomainSearchTree());
            tbox::plog << "External parts:\n" << finer_to_external.format("FE->", 3);

            TBOX_ERROR(
               "Internal library error: Failed to produce proper nesting.");

         } /* !finer_nests_in_new */

      } /* d_check_proper_nesting */

   } /* d_hierarchy->levelExists(new_ln + 1) */

   d_hierarchy->makeNewPatchLevel(new_ln, new_mapped_box_level);

   /*
    * Cache Connectors for new level.
    */
   boost::shared_ptr<hier::PatchLevel> new_level(
      d_hierarchy->getPatchLevel(new_ln));
   new_level->getBoxLevel()->getPersistentOverlapConnectors().
   cacheConnector(
      *new_level->getBoxLevel(),
      new_to_new);
   new_level->getBoxLevel()->getPersistentOverlapConnectors().
   cacheConnector(
      *tag_level->getBoxLevel(),
      new_to_tag);
   tag_level->getBoxLevel()->getPersistentOverlapConnectors().
   cacheConnector(
      *new_level->getBoxLevel(),
      tag_to_new);

   if (d_hierarchy->levelExists(new_ln + 1)) {
      /*
       * There is a level finer than new_ln.  Connect the new level to
       * the finer level.
       */
      hier::Connector* new_to_finer = new hier::Connector;
      hier::Connector* finer_to_new = new hier::Connector;

      t_bridge_new_to_finer->start();
      oca.bridgeWithNesting(
         *new_to_finer,
         *finer_to_new,
         *new_to_tag,
         tag_to_finer,
         finer_to_tag,
         *tag_to_new,
         zero_vector,
         -hier::IntVector::getOne(d_dim),
         d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln + 1));
      new_to_finer->setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      finer_to_new->setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      t_bridge_new_to_finer->stop();

      TBOX_ASSERT(
         new_to_finer->getConnectorWidth() ==
         d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln + 1));
      TBOX_ASSERT(
         finer_to_new->getConnectorWidth() ==
         d_hierarchy->getRequiredConnectorWidth(new_ln + 1, new_ln));

#ifdef DEBUG_CHECK_ASSERTIONS
      oca.assertOverlapCorrectness(*new_to_finer, false, true, true);
      oca.assertOverlapCorrectness(*finer_to_new, false, true, true);
#endif

      boost::shared_ptr<hier::PatchLevel> finer_level(
         d_hierarchy->getPatchLevel(new_ln + 1));

      new_level->getBoxLevel()->getPersistentOverlapConnectors().
      cacheConnector(
         *finer_level->getBoxLevel(),
         new_to_finer);

      finer_level->getBoxLevel()->getPersistentOverlapConnectors().
      cacheConnector(
         *new_level->getBoxLevel(),
         finer_to_new);
   }

   if (old_mapped_box_level) {

      /*
       * Connect old to new by bridging.
       *
       * Cache these Connectors for use when creating schedules to
       * transfer data from old to new.
       */

      hier::Connector* old_to_new = new hier::Connector;
      hier::Connector* new_to_old = new hier::Connector;
      t_bridge_new_to_old->start();
      oca.bridgeWithNesting(
         *old_to_new,
         *new_to_old,
         *old_to_tag,
         d_hierarchy->getConnector(tag_ln, tag_ln + 1),
         d_hierarchy->getConnector(tag_ln + 1, tag_ln),
         *tag_to_old,
         zero_vector,
         zero_vector,
         d_hierarchy->getRequiredConnectorWidth(new_ln, new_ln));
      old_to_new->setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      new_to_old->setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      t_bridge_new_to_old->stop();

      new_level->getBoxLevel()->getPersistentOverlapConnectors().
      cacheConnector(
         *old_fine_level->getBoxLevel(),
         new_to_old);
      old_fine_level->getBoxLevel()->getPersistentOverlapConnectors().
      cacheConnector(
         *new_level->getBoxLevel(),
         old_to_new);

   }

   d_hierarchy->getGridGeometry()->adjustMultiblockPatchLevelBoundaries(
      *d_hierarchy->getPatchLevel(new_ln));

   // "false" argument": const bool initial_time = false;
   d_tag_init_strategy->initializeLevelData(d_hierarchy,
      new_ln,
      regrid_time,
      d_hierarchy->levelCanBeRefined(new_ln),
      false,
      old_fine_level);

   /*
    * Destroy old patch level, if such a level existed prior to regrid.
    */
   old_fine_level.reset();
   if (d_barrier_and_time) {
      t_regrid_finer_create->stop();
   }
}

/*
 *************************************************************************
 * Check for boundary proximity violations.  Boxes may not be within
 * extend_ghosts of a physical boundary without touching the physical
 * boundary.
 *************************************************************************
 */
size_t
GriddingAlgorithm::checkBoundaryProximityViolation(
   const hier::BoxLevel& mapped_box_level,
   const hier::IntVector& extend_ghosts) const
{
   /*
    * 1. Compute the boundary regions of the boxes.
    *    a. Grow temporary boxes by the max ghost width.
    *    b. Remove orig boxes from the grown boxes to get ghost regions.
    *    c. Remove domain from the ghost region.
    * 2. Physical boundary regions which should not
    *    be less wide than the max ghost width.  If they are,
    *    it means they are partially inside the domain.
    */

   const hier::BaseGridGeometry& grid_geometry(
      *d_hierarchy->getGridGeometry()); 

   const hier::BoxContainer& periodic_domain_search_tree(
      grid_geometry.getPeriodicDomainSearchTree());

   hier::BoxContainer refined_periodic_domain_search_tree(
      periodic_domain_search_tree);
   refined_periodic_domain_search_tree.refine(mapped_box_level.getRefinementRatio());
   refined_periodic_domain_search_tree.makeTree(&grid_geometry);

   size_t nerr(0);

   for (hier::RealBoxConstIterator bi(mapped_box_level.getBoxes().realBegin());
        bi != mapped_box_level.getBoxes().realEnd(); ++bi) {

      hier::BoxContainer external_parts(*bi);
      external_parts.grow(extend_ghosts);
      external_parts.removeIntersections(
         mapped_box_level.getRefinementRatio(),
         refined_periodic_domain_search_tree);

      for (hier::BoxContainer::iterator bli(external_parts);
           bli != external_parts.end();
           ++bli) {
         hier::IntVector leftover_size((*bli).numberCells());
         for (int d = 0; d < d_dim.getValue(); ++d) {
            if (leftover_size(d) != 0 && leftover_size(d) < extend_ghosts(d)) {
               ++nerr;
               TBOX_WARNING("GriddingAlgorithm::makeFinerLevel:\n"
                  << "User-specified box (refined) " << *bi
                  << " violates boundary proximity.\n"
                  << "In dimension " << d << ", it is "
                  << extend_ghosts(d) - leftover_size(d)
                  << " cells from a physical domain boundary.\n"
                  << "All boxes must be at least " << extend_ghosts
                  << " from physical boundaries or touching the physical boundary.");
            }
         }
      }

   }

   return nerr;
}

/*
 *******************************************************************
 * Check domain boxes for violations of user constraints.
 *******************************************************************
 */
void
GriddingAlgorithm::checkDomainBoxes(const hier::BoxContainer& domain_boxes) const {

   hier::IntVector smallest_patch(d_dim);
   hier::IntVector largest_patch(d_dim);
   hier::IntVector extend_ghosts(d_dim);
   {
      hier::IntVector smallest_box_to_refine(d_dim);
      // "false" argument: for_building_finer level = false
      getGriddingParameters(
         smallest_patch,
         smallest_box_to_refine,
         largest_patch,
         extend_ghosts,
         0,
         false);
   }

   /*
    * Check minimum size violations.
    */
   int i = 0;
   for (hier::BoxContainer::const_iterator itr(domain_boxes);
        itr != domain_boxes.end(); ++itr, ++i) {

      hier::Box test_box = *itr;
      for (int dir = 0; dir < d_dim.getValue(); dir++) {

         if (test_box.numberCells(dir) < smallest_patch(dir)) {

            int error_coarsen_ratio =
               d_tag_init_strategy->getErrorCoarsenRatio();
            if (error_coarsen_ratio > 1) {
               TBOX_ERROR(
                  d_object_name << ": " << "\ndomain Box " << i << ", " << test_box
                                << ", violates the minimum patch size constraints."
                                << "\nVerify that boxes are larger than"
                                << "the maximum ghost width and/or"
                                << "\nthe specified minimum patch size."
                                << "\nNOTE: to assure the constraints are"
                                << "properly enforced during coarsening for"
                                << "\nerror computation, the minimum patch"
                                << "size is the smallest patch size multiplied"
                                << "\nby the error coarsen ratio, which is "
                                << error_coarsen_ratio
                                << " in this case."
                                << std::endl);
            } else {
               TBOX_ERROR(
                  d_object_name << ": "
                                << "\ndomain Box " << i << ", " << test_box
                                << ", violates the minimum patch size constraints."
                                << "\nVerify that boxes are larger than"
                                << "the maximum ghost width and/or"
                                << "\nthe specified minimum patch size."
                                << std::endl);
            }
         }
      }
   }

   /*
    * Check for overlapping boxes.
    * TODO: This check only works for single-block.
    */
   if (domain_boxes.boxesIntersect()) {
      TBOX_ERROR(d_object_name << ":  "
                               << "Boxes specified for coarsest level "
                               << "contain intersections with each other!");
   }

   /*
    * Check for violations of implementation of TagAndInitStrategy.
    */
   if ((d_hierarchy->getMaxNumberOfLevels() > 1)
       && (!d_tag_init_strategy->coarsestLevelBoxesOK(domain_boxes))) {
      TBOX_ERROR(d_object_name << ":  "
                               << "level gridding strategy encountered"
                               << " a problem with the domain boxes!");
   }
}

/*
 *******************************************************************
 * Check for non-nesting user-specified boxes.
 *******************************************************************
 */
void
GriddingAlgorithm::checkNonnestingUserBoxes(
   const hier::Connector& new_to_tag,
   const hier::IntVector& nesting_buffer) const
{

   const hier::BoxLevel& new_mapped_box_level(new_to_tag.getBase());

   hier::BoxLevel violating_parts(d_dim);
   hier::Connector new_to_violating_parts;

   hier::BoxLevelConnectorUtils mblc_utils;
   mblc_utils.computeExternalParts(
      violating_parts,
      new_to_violating_parts,
      new_to_tag,
      -nesting_buffer,
      d_hierarchy->getGridGeometry()->getDomainSearchTree());

   if (violating_parts.getGlobalNumberOfBoxes() > 0) {

      tbox::perr << "GriddingAlgorihtm: user-specified refinement boxes\n"
                 << "violates nesting requirement.  Diagnostics will be\n"
                 << "writen to log files." << std::endl;
      const std::string left_margin("ERR: ");
      tbox::plog
      << left_margin << "Tag BoxLevel:\n" << new_to_tag.getHead().format(left_margin, 2)
      << left_margin << "User-specified boxes:\n" << new_mapped_box_level.format(left_margin, 2)
      << left_margin << "Violating parts:\n" << violating_parts.format(left_margin, 2)
      << left_margin << "User-specified boxes and their violating parts:\n"
      << new_to_violating_parts.format(left_margin, 2);

      if (d_check_nonnesting_user_boxes == 'e') {
         TBOX_ERROR("Exiting due to above error");
      }
      if (d_check_nonnesting_user_boxes == 'w') {
         TBOX_WARNING("Proceeding with nesting violation as requested.\n"
            << "SAMRAI is not guaranteed to work with nesting"
            << "violations!");
      }

   }
}

/*
 *******************************************************************
 * Check for non-nesting user-specified boxes.
 *******************************************************************
 */
void
GriddingAlgorithm::checkBoundaryProximityViolation(
   const int tag_ln,
   const hier::BoxLevel& new_mapped_box_level) const
{
   hier::IntVector extend_ghosts(d_dim);
   hier::IntVector smallest_patch(d_dim);
   hier::IntVector smallest_box_to_refine(d_dim);
   hier::IntVector largest_patch(d_dim);
   getGriddingParameters(smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts,
      tag_ln,
      false);

   const size_t nerr(
      checkBoundaryProximityViolation(
         new_mapped_box_level,
         extend_ghosts));
   if (nerr > 0 && d_check_boundary_proximity_violation == 'e') {
      TBOX_ERROR("GriddingAlgorithm::makeFinerLevel: User error:\n"
         << "Making level " << tag_ln + 1 << ".\n"
         << "New boxes violate boundary proximity.\n"
         << "All boxes must be at least " << extend_ghosts
         << " from physical boundaries or touching the physical boundary.");
   }
}

/*
 *************************************************************************
 *                                                                       *
 *                                                                       *
 *************************************************************************
 */
void
GriddingAlgorithm::recordStatistics(
   double current_time)
{
#ifdef GA_RECORD_STATS
// GA_RECORD_STATS is defined in GriddingAlgorithm.h
/*
 * For statistics, record number of cells and patches on new level.
 */
   for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
      int level_gridcells = 0;
      int level_local_patches = 0;
      if (ln < d_hierarchy->getNumberOfLevels()) {
         const boost::shared_ptr<hier::PatchLevel>& patch_level =
            d_hierarchy->getPatchLevel(ln);
         level_gridcells = patch_level->getLocalNumberOfCells();
         level_local_patches = patch_level->getLocalNumberOfPatches();
      }
      d_boxes_stat[ln]->recordProcStat(double(level_local_patches));
      d_cells_stat[ln]->recordProcStat(double(level_gridcells));
      d_timestamp_stat[ln]->recordProcStat(double(current_time));
   }
#endif
}

/*
 *************************************************************************
 *                                                                       *
 *                                                                       *
 *************************************************************************
 */
void
GriddingAlgorithm::printStatistics(
   std::ostream& s) const
{
#ifdef GA_RECORD_STATS
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());
   /*
    * Output statistics.
    */
   // Collect statistic on mesh size.
   tbox::Statistician* statn = tbox::Statistician::getStatistician();

   statn->finalize(false);
   // statn->printLocalStatData(s);
   if (mpi.getRank() == 0) {
      // statn->printAllGlobalStatData(s);
      for (int ln = 0; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
         tbox::Statistic& cstat = *d_cells_stat[ln];
         tbox::Statistic& bstat = *d_boxes_stat[ln];
         tbox::Statistic& tstat = *d_timestamp_stat[ln];
         s << "statistic " << cstat.getName() << ":" << std::endl;
         if (0) {
            s << "Global: \n";
            statn->printGlobalProcStatDataFormatted(cstat.getInstanceId(), s);
         }
         s
         <<
         " Seq#  SimTime       C-Sum      C-Avg      C-Min   ->    C-Max  C-NormDiff  B-Sum B-Avg B-Min -> B-Max B-NormDiff C/B-Avg\n";
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
         for (int sn = 0; sn < cstat.getStatSequenceLength(); ++sn) {
            double csum = statn->getGlobalProcStatSum(cstat.getInstanceId(), sn);
            double cmax = statn->getGlobalProcStatMax(cstat.getInstanceId(), sn);
            double cmin = statn->getGlobalProcStatMin(cstat.getInstanceId(), sn);
            double cdiffnorm = cmax != 0 ? 1.0 - cmin / cmax : 0;
            double bsum = statn->getGlobalProcStatSum(bstat.getInstanceId(), sn);
            double bmax = statn->getGlobalProcStatMax(bstat.getInstanceId(), sn);
            double bmin = statn->getGlobalProcStatMin(bstat.getInstanceId(), sn);
            double bdiffnorm = bmax != 0 ? 1.0 - bmin / bmax : 0;
            double stime = statn->getGlobalProcStatMin(
                  tstat.getInstanceId(), sn);
            s << std::setw(3) << sn << "  "
              << std::scientific << std::setprecision(6) << std::setw(12)
              << stime
              << " "
              << std::fixed << std::setprecision(0)
              << std::setw(10) << csum << " "
              << std::setw(10) << csum / mpi.getSize() << " "
              << std::setw(10) << cmin << " -> "
              << std::setw(10) << cmax
              << "  " << std::setw(4) << std::setprecision(4) << cdiffnorm
              << "  "
              << std::fixed << std::setprecision(0)
              << std::setw(6) << bsum << " "
              << std::setw(5) << bsum / mpi.getSize() << " "
              << std::setw(5) << bmin << "  ->"
              << std::setw(5) << bmax
              << "   " << std::setw(4) << std::setprecision(4) << bdiffnorm
              << std::setw(10) << std::setprecision(0)
              << (bsum != 0 ? csum / bsum : 0)
              << std::endl;
         }
      }
   }
#endif
}

/*
 *************************************************************************
 * All tags reside in the tag level.  But due to nesting restrictions,
 * not all cells in the level may be refined.
 *
 * We look for the portions of the level that would violate nesting if
 * refined.  Any tags there are nonnesting tags.
 *
 * TODO: remove level from this interface.  tag_ln is sufficient.
 *************************************************************************
 */
void
GriddingAlgorithm::checkNonrefinedTags(
   const hier::PatchLevel& level,
   int tag_ln) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, level);

   const hier::BoxLevel& tag_mapped_box_level = *d_hierarchy->getBoxLevel(tag_ln);
   hier::BoxLevel violator(d_dim);
   hier::Connector tag_to_violator;
   const hier::Connector& tag_mapped_box_level_to_self = d_hierarchy->getConnector(tag_ln,
         tag_ln);
   computeNestingViolator(
      violator,
      tag_to_violator,
      tag_mapped_box_level,
      tag_mapped_box_level_to_self,
      tag_mapped_box_level_to_self,
      tag_ln);

   /*
    * Check for user-tagged cells in the violating parts of the tag level.
    */
   math::PatchCellDataBasicOps<int> dataop;
   math::PatchCellDataOpsInteger dataopi;
   int maxval = 0;
   for (hier::Connector::ConstNeighborhoodIterator ei = tag_to_violator.begin();
        ei != tag_to_violator.end(); ++ei) {
      const hier::BoxId& mapped_box_id = *ei;
      boost::shared_ptr<pdat::CellData<int> > tag_data(
         level.getPatch(mapped_box_id)->getPatchData(d_tag_indx),
         boost::detail::dynamic_cast_tag());
      for (hier::Connector::ConstNeighborIterator na = tag_to_violator.begin(ei);
           na != tag_to_violator.end(ei); ++na) {
         const hier::Box& vio_mapped_box = *na;
         maxval = dataop.max(tag_data, vio_mapped_box);
         if (maxval > 0) {
            break;
         }
      }
      if (maxval > 0) {
         break;
      }
   }
   const tbox::SAMRAI_MPI mpi(tag_mapped_box_level.getMPI());
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&maxval, 1, MPI_MAX);
   }

   if (maxval > 0) {
      if (d_check_nonrefined_tags == 'w') {
         TBOX_WARNING("User code has tagged cells in\n"
            << "violation of nesting requirements.\n"
            << "Violating tags will be discarded.\n"
            << "See GriddingAlgorithm::checkNonrefinedTags()\n");
      } else if (d_check_nonrefined_tags == 'e') {
         TBOX_ERROR("User code has tagged cells in\n"
            << "violation of nesting requirements.\n"
            << "See GriddingAlgorithm::checkNonrefinedTags()\n");
      }
   }
}

/*
 *******************************************************************
 * Reset tag buffering data to be able to handle given buffer size.
 *******************************************************************
 */

void GriddingAlgorithm::resetTagBufferingData(const int tag_buffer)
{
   d_buf_tag_ghosts = hier::IntVector(d_dim, tag_buffer);

   d_bdry_fill_tags.reset();

   hier::VariableDatabase* var_db =
      hier::VariableDatabase::getDatabase();

   /*
    * Remove d_buf_tag from the VariableDatabase and re-register it with
    * the new ghost width.
    */
   var_db->removeInternalSAMRAIVariablePatchDataIndex(d_buf_tag_indx);

   (*s_buf_tag_indx)[d_dim.getValue() - 1] =
      var_db->registerInternalSAMRAIVariable(d_buf_tag,
         d_buf_tag_ghosts);

   d_buf_tag_indx = (*s_buf_tag_indx)[d_dim.getValue() - 1];

   if (d_hierarchy->getGridGeometry()->getNumberBlocks() > 1) {
      TBOX_ASSERT(d_mb_tagger_strategy); 
      d_mb_tagger_strategy->setScratchTagPatchDataIndex(d_buf_tag_indx);
   }

   d_bdry_fill_tags.reset(new xfer::RefineAlgorithm(d_dim));

   d_bdry_fill_tags->
   registerRefine(d_buf_tag_indx,
      d_buf_tag_indx,
      d_buf_tag_indx,
      boost::shared_ptr<hier::RefineOperator>());
}

/*
 *************************************************************************
 *************************************************************************
 */
void
GriddingAlgorithm::checkOverlappingPatches(
   const hier::Connector& mapped_box_level_to_self) const
{
   bool has_overlap = false;
   const hier::BoxLevel& mapped_box_level = mapped_box_level_to_self.getBase();
   const hier::BaseGridGeometry& grid_geom =
      *mapped_box_level.getGridGeometry();
   const hier::IntVector& ratio = mapped_box_level.getRefinementRatio();

   for (hier::Connector::ConstNeighborhoodIterator ei = mapped_box_level_to_self.begin();
        ei != mapped_box_level_to_self.end() && !has_overlap; ++ei) {

      const hier::Box& mapped_box = *mapped_box_level.getBoxStrict(*ei);

      for (hier::Connector::ConstNeighborIterator na = mapped_box_level_to_self.begin(ei);
           na != mapped_box_level_to_self.end(ei) && !has_overlap;
           ++na) {
         const hier::Box& nabr = *na;

         if (!nabr.isIdEqual(mapped_box)) {
            if (nabr.getBlockId() == mapped_box.getBlockId()) {
               has_overlap = nabr.intersects(mapped_box);
            } else {
               hier::Box nabr_box(nabr);
               grid_geom.transformBox(nabr_box,
                  ratio,
                  mapped_box.getBlockId(),
                  nabr.getBlockId());
               has_overlap = nabr_box.intersects(mapped_box);
            }
         }
      }
   }

   if (has_overlap) {
      if (d_check_overlapping_patches == 'w') {
         TBOX_WARNING(
            "PatchLevel has patches which overlap in index space\n"
            << "See GriddingAlgorithm::checkOverlappingPatches().\n"
            <<
            "Note that setting allow_patches_smaller_than_minimum_size_to_prevent_overlaps = FALSE\n"
            <<
            "in the PatchHierarchy can allow some patches to violate min size in order to prevent overlap.\n");
      } else if (d_check_overlapping_patches == 'e') {
         TBOX_ERROR(
            "PatchLevel has patches which overlap in index space\n"
            << "See GriddingAlgorithm::checkOverlappingPatches().\n"
            <<
            "Note that setting allow_patches_smaller_than_minimum_size_to_prevent_overlaps = FALSE\n"
            <<
            "in the PatchHierarchy can allow some patches to violate min size in order to prevent overlap.\n");
      }
   }
}

/*
 *************************************************************************
 *
 * For cases where tagging is not performed read the new level boxes
 * either from user input or from stored level boxes.
 *
 *************************************************************************
 */
void
GriddingAlgorithm::readLevelBoxes(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& coarser_to_new,
   hier::Connector& new_to_coarser,
   const int tag_ln,
   const double regrid_time,
   bool& remove_old_fine_level)
{
   TBOX_ASSERT((tag_ln >= 0)
      && (tag_ln <= d_hierarchy->getFinestLevelNumber()));

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, *d_hierarchy, new_mapped_box_level);

   const hier::BoxLevel& coarser_mapped_box_level = *d_hierarchy->getBoxLevel(
         tag_ln);

   int fine_level_number = tag_ln + 1;
   hier::BoxContainer boxes_to_refine;

   /*
    * Access the user supplied refine boxes.  The
    * "new_level_has_new_boxes" boolean specifies whether the
    * level boxes have changed from the last time
    * getUserSuppliedRefineBoxes() was called.  If they have changed,
    * it returns true.  If they are unchanged, it returns false.
    */
   bool new_level_has_new_boxes = true;
   if (d_tag_init_strategy->refineUserBoxInputOnly()) {

      new_level_has_new_boxes = d_tag_init_strategy->
         getUserSuppliedRefineBoxes(boxes_to_refine,
            tag_ln,
            regrid_time);

   }

   /*
    * If "new_level_has_new_boxes" is false we wish to keep the
    * existing fine level intact.  Avoid further work by setting
    * the parameter "compute_load_balanced_level_boxes" to false
    * and indicate that we want to avoid removing the old fine level
    * by setting "remove_old_fine_level" to false.
    */
   bool compute_load_balanced_level_boxes = true;
   if (!new_level_has_new_boxes) {
      compute_load_balanced_level_boxes = false;
      remove_old_fine_level = false;
   }

   /*
    * If we are using the nonuniform load balance option, we
    * still need to redo the load balance and construct a new level,
    * even if the level boxes have not changed.
    */

   if (d_load_balancer->getLoadBalanceDependsOnPatchData(fine_level_number)
       && boxes_to_refine.size() > 0) {
      compute_load_balanced_level_boxes = true;
      remove_old_fine_level = true;
   }

   /*
    * If the boxes_to_refine are empty, this implies that no
    * refinement is desired so a new finer level will NOT be
    * constructed.  In this case, avoid load balance steps and
    * specify that we want to remove the old fine level.
    */
   if (boxes_to_refine.size() == 0) {
      compute_load_balanced_level_boxes = false;
      remove_old_fine_level = true;
   }

   if (compute_load_balanced_level_boxes) {

      hier::BoxLevel unbalanced_mapped_box_level(d_dim);
      unbalanced_mapped_box_level.initialize(
         coarser_mapped_box_level.getRefinementRatio(),
         coarser_mapped_box_level.getGridGeometry(),
         d_hierarchy->getMPI(),
         hier::BoxLevel::GLOBALIZED);
      hier::LocalId i(0);
      for (hier::BoxContainer::iterator itr(boxes_to_refine);
           itr != boxes_to_refine.end(); ++itr, ++i) {
         hier::Box unbalanced_mapped_box(*itr, i, 0);
         unbalanced_mapped_box.setBlockId(hier::BlockId(0));
         unbalanced_mapped_box_level.addBox(unbalanced_mapped_box);
      }

      const hier::IntVector& ratio =
         d_hierarchy->getRatioToCoarserLevel(fine_level_number);

      new_mapped_box_level = unbalanced_mapped_box_level;
      coarser_to_new.clearNeighborhoods();
      coarser_to_new.setBase(coarser_mapped_box_level);
      coarser_to_new.setHead(new_mapped_box_level);
      coarser_to_new.setWidth(
         d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln + 1),
         true);
      new_to_coarser.clearNeighborhoods();
      new_to_coarser.setBase(new_mapped_box_level);
      new_to_coarser.setHead(coarser_mapped_box_level);
      new_to_coarser.setWidth(hier::IntVector::ceilingDivide(
            d_hierarchy->getRequiredConnectorWidth(tag_ln + 1, tag_ln), ratio),
            true);
      const hier::OverlapConnectorAlgorithm oca;
      oca.findOverlaps(coarser_to_new);
      oca.findOverlaps(new_to_coarser);

      hier::IntVector smallest_patch(d_dim);
      hier::IntVector largest_patch(d_dim);
      hier::IntVector extend_ghosts(d_dim);
      {
         hier::IntVector smallest_box_to_refine(d_dim);
         // "false" argument: for_building_finer level = false
         getGriddingParameters(smallest_patch,
            smallest_box_to_refine,
            largest_patch,
            extend_ghosts,
            fine_level_number,
            false);
      }

      hier::IntVector patch_cut_factor(d_dim, 1);

      t_load_balance0->start();
      d_load_balancer0->loadBalanceBoxLevel(
         new_mapped_box_level,
         new_to_coarser,
         coarser_to_new,
         d_hierarchy,
         tag_ln,
         hier::Connector(),
         hier::Connector(),
         smallest_patch,
         largest_patch,
         d_hierarchy->getDomainBoxLevel(),
         extend_ghosts,
         patch_cut_factor);
      t_load_balance0->stop();

      refineNewBoxLevel(new_mapped_box_level,
         coarser_to_new,
         new_to_coarser,
         ratio);
      if (d_sequentialize_patch_indices) {
         renumberBoxes(new_mapped_box_level,
            coarser_to_new,
            new_to_coarser,
            false,
            true);
      }

      oca.findOverlaps(coarser_to_new);
      oca.findOverlaps(new_to_coarser);

      /*
       * Periodic relationships exist in new_to_coarser, but are not
       * complete because new doesn't have any periodic images yet.
       * Remove these relationships to make new<==>coarser proper
       * transposes.
       */
      new_to_coarser.removePeriodicRelationships();

      const hier::Connector& coarser_to_coarser =
         d_hierarchy->getConnector(tag_ln, tag_ln);
      const hier::BoxLevelConnectorUtils dlbg_edge_utils;
      dlbg_edge_utils.addPeriodicImagesAndRelationships(
         new_mapped_box_level,
         new_to_coarser,
         coarser_to_new,
         d_hierarchy->getGridGeometry()->getDomainSearchTree(),
         coarser_to_coarser);
      new_mapped_box_level.finalize();
   }
}

/*
 *************************************************************************
 * Set all tags on a level to tag_value.
 *************************************************************************
 */

void
GriddingAlgorithm::fillTags(
   const int tag_value,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const int tag_index) const
{
   TBOX_ASSERT((tag_value == d_true_tag) || (tag_value == d_false_tag));
   TBOX_ASSERT(tag_level);
   TBOX_ASSERT(tag_index == d_tag_indx || tag_index == d_buf_tag_indx);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *tag_level);

   t_fill_tags->start();

   for (hier::PatchLevel::iterator ip(tag_level->begin());
        ip != tag_level->end(); ++ip) {

      const boost::shared_ptr<hier::Patch>& patch = *ip;
      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_index),
         boost::detail::dynamic_cast_tag());
      TBOX_ASSERT(tag_data);

      tag_data->fill(tag_value);

   }
   t_fill_tags->stop();
}

/*
 *************************************************************************
 *
 * Set each integer value in specified tag array to tag_value where
 * patch level intersects given box array.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::fillTagsFromBoxLevel(
   const int tag_value,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const int tag_index,
   const hier::Connector& tag_level_to_fill_mapped_box_level,
   const bool interior_only,
   const hier::IntVector& fill_box_growth) const
{
   TBOX_ASSERT((tag_value == d_true_tag) || (tag_value == d_false_tag));
   TBOX_ASSERT(tag_level);
   TBOX_ASSERT(tag_index == d_tag_indx || tag_index == d_buf_tag_indx);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, *tag_level, fill_box_growth);

   /*
    * This method assumes fill is finer than tag, but that is easy to
    * change, if needed.
    */
   TBOX_ASSERT(tag_level_to_fill_mapped_box_level.getHeadCoarserFlag() == false);

   t_fill_tags->start();

   const hier::OverlapConnectorAlgorithm oca;

   const boost::shared_ptr<const hier::BaseGridGeometry>& grid_geom(
      d_hierarchy->getGridGeometry());

   const hier::IntVector& ratio = tag_level_to_fill_mapped_box_level.getRatio();

   const hier::IntVector growth_in_tag_resolution =
      hier::IntVector::ceilingDivide(fill_box_growth,
         tag_level_to_fill_mapped_box_level.getRatio());

   for (hier::PatchLevel::iterator ip(tag_level->begin());
        ip != tag_level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& patch = *ip;

      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_index),
         boost::detail::dynamic_cast_tag());

      TBOX_ASSERT(tag_data);

      const hier::BoxId& mapped_box_id(patch->getBox().getId());

      NeighborSet neighbors;

      oca.extractNeighbors(
         neighbors,
         tag_level_to_fill_mapped_box_level,
         mapped_box_id,
         growth_in_tag_resolution);

      for (NeighborSet::const_iterator
           ni = neighbors.begin(); ni != neighbors.end(); ++ni) {
         const hier::Box& neighbor(*ni);
         hier::Box box = neighbor;
         box.grow(fill_box_growth);
         box.coarsen(ratio);
         if (neighbor.getBlockId() != patch->getBox().getBlockId()) {
            grid_geom->transformBox(box,
               tag_level->getRatioToLevelZero(),
               patch->getBox().getBlockId(),
               neighbor.getBlockId());
         }
         if (interior_only) {
            box = box * tag_data->getBox();
         }
         tag_data->fill(tag_value, box);
      }

   }
   t_fill_tags->stop();
}

/*
 *************************************************************************
 *
 * Buffer each integer tag with given value on the patch level by the
 * specified buffer size.  Note that the patch data indexed by
 * d_buf_tag_indx is used temporarily to buffer the tag data. The
 * communication of ghost cell (i.e., buffer) information forces all
 * tags on all patch interiors to represent a consistent buffering of
 * the original configuration of tagged cells.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::bufferTagsOnLevel(
   const int tag_value,
   const boost::shared_ptr<hier::PatchLevel>& level,
   const int buffer_size) const
{
   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::bufferTagsOnLevel: entered with tag_ln = "
      << level->getLevelNumber() << "\n";
   }

   TBOX_ASSERT((tag_value == d_true_tag) || (tag_value == d_false_tag));
   TBOX_ASSERT(level);
   TBOX_ASSERT(buffer_size >= 0);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *level);
   /*
    * Start timer for this method.
    */
   t_buffer_tags->start();

   /*
    * Set temporary buffered tags based on buffer width and
    * distance from actual tags.
    */
   const int not_tag = ((tag_value == d_true_tag) ? d_false_tag : d_true_tag);
   for (hier::PatchLevel::iterator ip1(level->begin());
        ip1 != level->end(); ++ip1) {
      const boost::shared_ptr<hier::Patch>& patch = *ip1;

      boost::shared_ptr<pdat::CellData<int> > buf_tag_data(
         patch->getPatchData(d_buf_tag_indx),
         boost::detail::dynamic_cast_tag());
      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(d_tag_indx),
         boost::detail::dynamic_cast_tag());

      buf_tag_data->fillAll(not_tag);

      const hier::Box& interior(patch->getBox());

      pdat::CellIterator icend(interior, false);
      for (pdat::CellIterator ic(interior, true); ic != icend; ++ic) {
         if ((*tag_data)(*ic) == tag_value) {
            (*buf_tag_data)(*ic) = d_true_tag;
         }
      }
   }

   /*
    * Communicate boundary data for buffered tag array so that tags
    * near patch boundaries will become buffered properly.
    */
   const double dummy_time = 0.0;

   t_bdry_fill_tags_comm->start();
   d_bdry_sched_tags[level->getLevelNumber()]->fillData(dummy_time, false);
   t_bdry_fill_tags_comm->stop();

   /*
    * Buffer tags on patch interior according to buffered tag data.
    */
   for (hier::PatchLevel::iterator ip2(level->begin());
        ip2 != level->end(); ++ip2) {
      const boost::shared_ptr<hier::Patch>& patch = *ip2;

      boost::shared_ptr<pdat::CellData<int> > buf_tag_data(
         patch->getPatchData(d_buf_tag_indx),
         boost::detail::dynamic_cast_tag());
      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(d_tag_indx),
         boost::detail::dynamic_cast_tag());

      const hier::Box& tag_box(tag_data->getBox());
      const hier::BlockId& tag_box_block_id = tag_box.getBlockId();
      hier::Box buf_tag_box(tag_box);
      buf_tag_box.grow(hier::IntVector(d_dim, buffer_size));

      tag_data->fillAll(not_tag);

      pdat::CellIterator icend(buf_tag_box, false);
      for (pdat::CellIterator ic(buf_tag_box, true); ic != icend; ++ic) {
         if ((*buf_tag_data)(*ic) == d_true_tag) {
            hier::Box buf_box(*ic - buffer_size,
               *ic + buffer_size,
               tag_box_block_id);
            tag_data->fill(tag_value, buf_box);
         }
      }

   }

   t_buffer_tags->stop();
}

/*
 *************************************************************************
 *
 * Given a patch level, determine an appropriate array of boxes from
 * which a new finer level may be constructed.  That is, find an array
 * of boxes that covers all tags having the specified tag value.  Note
 * that it is assumed that the integer tag arrays have been set
 * properly; i.e., cells have been tagged through error estimation and
 * the tags have been buffered to ensure disturbances remain on fine
 * level until next regrid occurs.  Note that load balancing is
 * performed once an appropriate list of boxes containing the tags is
 * found.  This procedure massages the list of boxes further and then
 * assigns each to a single processor (i.e., the mapping).
 *
 *************************************************************************
 */

void
GriddingAlgorithm::findRefinementBoxes(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const int tag_ln) const
{
   TBOX_ASSERT((tag_ln >= 0)
      && (tag_ln <= d_hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(d_hierarchy->getPatchLevel(tag_ln));
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, new_mapped_box_level, *d_hierarchy);

   if (d_print_steps) {
      tbox::plog
      << "GriddingAlgorithm::findRefinementBoxes entered with tag_ln = "
      << tag_ln << "\n";
   }

   const hier::OverlapConnectorAlgorithm oca;
   const hier::BoxLevelConnectorUtils dlbg_edge_utils;

   /*
    * Start timer for this method.
    */
   if (d_barrier_and_time) {
      t_find_refinement->barrierAndStart();
   }

   const hier::BoxLevel& tag_mapped_box_level = *d_hierarchy->getBoxLevel(tag_ln);

   const int new_ln = tag_ln + 1;

   /*
    * Construct list of boxes covering the true tags on the level.
    * Note that box list will be contained in the bounding box
    * but will not be contained in the list of proper nesting boxes,
    * in general.  So we intersect the box list against the list of
    * nesting boxes.  Note that this may produce boxes which are too
    * small.  Thus, boxes are regrown later.
    */

   hier::IntVector smallest_patch(d_dim);
   hier::IntVector smallest_box_to_refine(d_dim);
   hier::IntVector largest_patch(d_dim);
   hier::IntVector extend_ghosts(d_dim);
   // "true" argument: for_building_finer level = true
   getGriddingParameters(smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts,
      new_ln,
      true);

   const hier::IntVector smallest_patch_in_tag_space =
      hier::IntVector::ceilingDivide(smallest_patch,
         d_hierarchy->getRatioToCoarserLevel(new_ln));
   const hier::IntVector largest_patch_in_tag_space =
      largest_patch / d_hierarchy->getRatioToCoarserLevel(new_ln);
   const hier::IntVector extend_ghosts_in_tag_space =
      hier::IntVector::ceilingDivide(extend_ghosts,
         d_hierarchy->getRatioToCoarserLevel(new_ln));

   boost::shared_ptr<hier::PatchLevel> level(
      d_hierarchy->getPatchLevel(tag_ln));

   if (d_print_steps) {
      tbox::plog << "GriddingAlgorithm::findRefinementBoxes: clustering\n";
   }

   t_find_boxes_containing_tags->barrierAndStart();
   hier::IntVector ratio = d_hierarchy->getRatioToCoarserLevel(new_ln);

   /*
    * Compute the width for tag<==>cluster.  This width be wide enough to
    * guarantee completeness for tag<==>new when we massage the
    * cluster into the new level.  In the massage step, we grow the new boxes.
    * The width of tag<==>cluster must be big enough to see new overlaps
    * generated by the growths.
    */
   hier::IntVector tag_to_cluster_width =
      d_hierarchy->getRequiredConnectorWidth(tag_ln, tag_ln + 1);
   if (d_extend_to_domain_boundary) {
      // For extending boxes to domain boundary by amount of extend_ghosts_in_tag_space.
      tag_to_cluster_width += extend_ghosts_in_tag_space;
   }
   // For growing within domain by amount of smallest_box_to_refine.
   tag_to_cluster_width += smallest_box_to_refine;

   const int nblocks = d_hierarchy->getGridGeometry()->getNumberBlocks();
   if (nblocks == 1) {
      hier::Box bounding_box(d_dim);
      bounding_box =
         d_hierarchy->getBoxLevel(tag_ln)->getGlobalBoundingBox(0);
      d_box_generator->findBoxesContainingTags(
         new_mapped_box_level,
         tag_to_new,
         new_to_tag,
         level, d_tag_indx, d_true_tag, bounding_box,
         smallest_box_to_refine,
         getEfficiencyTolerance(tag_ln),
         getCombineEfficiency(tag_ln),
         tag_to_cluster_width,
         hier::BlockId(0),
         hier::LocalId(0));
   } else {
      /*
       * TODO: The following loop is an inefficient work-around to
       * handle multiple blocks through the single-block
       * findBoxesContainingTags interface.  The interfaces should be
       * changed to support multiblock.
       */
      hier::BoxContainer accumulated_mapped_boxes;
      hier::LocalId first_local_id(0);
      for (int bn = 0; bn < nblocks; ++bn) {
         /*
          * Determine single smallest bounding box for all nesting boxes.
          */
         hier::Box bounding_box(d_dim);
         bounding_box =
            d_hierarchy->getBoxLevel(tag_ln)->getGlobalBoundingBox(bn);

         if (!bounding_box.isEmpty()) {
            d_box_generator->findBoxesContainingTags(
               new_mapped_box_level,
               tag_to_new,
               new_to_tag,
               level, d_tag_indx, d_true_tag, bounding_box,
               smallest_box_to_refine,
               getEfficiencyTolerance(tag_ln),
               getCombineEfficiency(tag_ln),
               tag_to_cluster_width,
               hier::BlockId(bn),
               first_local_id);
            accumulated_mapped_boxes.insert(
               new_mapped_box_level.getBoxes().begin(),
               new_mapped_box_level.getBoxes().end());
            if (accumulated_mapped_boxes.size() > 0) {
               first_local_id =
                  accumulated_mapped_boxes.back().getId().getLocalId() + 1;
            }
         }

      }
#ifdef DEBUG_CHECK_ASSERTIONS
      std::set<int> local_ids;
      for (hier::BoxContainer::iterator
           ac_itr = accumulated_mapped_boxes.begin();
           ac_itr != accumulated_mapped_boxes.end(); ++ac_itr) {
         local_ids.insert(ac_itr->getId().getLocalId().getValue());
      }
      TBOX_ASSERT(static_cast<int>(local_ids.size()) == accumulated_mapped_boxes.size());
#endif
      const hier::BoxLevel& tag_mapped_box_level(tag_to_new.getBase());
      new_mapped_box_level.swapInitialize(
         accumulated_mapped_boxes,
         new_mapped_box_level.getRefinementRatio(),
         new_mapped_box_level.getGridGeometry(),
         new_mapped_box_level.getMPI());

      /*
       * Set up tag<==>new Connectors.  We cannot use the neighborhood sets
       * from findBoxesContainingTags because those do not include cross-block
       * neighbors.
       */
      tag_to_new.clearNeighborhoods();
      tag_to_new.setBase(tag_mapped_box_level);
      tag_to_new.setHead(new_mapped_box_level);
      tag_to_new.setWidth(tag_to_cluster_width, true);
      new_to_tag.clearNeighborhoods();
      new_to_tag.setBase(new_mapped_box_level);
      new_to_tag.setHead(tag_mapped_box_level);
      new_to_tag.setWidth(tag_to_cluster_width, true);
      oca.findOverlaps(tag_to_new);
      oca.findOverlaps(new_to_tag);

   }
   t_find_boxes_containing_tags->stop();

   if (new_mapped_box_level.getGlobalNumberOfBoxes() > 0) {

      if (d_check_connectors) {
         /*
          * At this stage, there are no edges to periodic images yet, so
          * don't check for them.
          */
         if (d_print_steps) {
            tbox::plog
            <<
            "GriddingAlgorithm::findRefinementBoxes: checking new-->tag from findBoxesContainingTags\n";
         }
         TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag,
               false,
               true,
               true) == 0);
         if (d_print_steps) {
            tbox::plog
            <<
            "GriddingAlgorithm::findRefinementBoxes: checking tag-->new from findBoxesContainingTags\n";
         }
         TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new,
               false,
               true,
               true) == 0);
      }

      t_box_massage->start();

      {
         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes: enforcing overflow nesting\n";
         }

         if (d_barrier_and_time) {
            t_limit_overflow->barrierAndStart();
         }
         /*
          * Do not allow the new mapped_box_level to overflow the tag mapped_box_level.
          * If we want to allow the overflow, we have to add the
          * overflow ammount to width of tag->new.  Such additions
          * may make the ABR algorithm slower, because more
          * non-contributing processors would have to be included
          * in the contributing group (unless ABR keep track of
          * the non-contributing processors and don't seek tag
          * histogram from them).
          */
         hier::BoxLevel nested_mapped_box_level(d_dim);
         hier::Connector unnested_to_nested;
         makeOverflowNestingMap(
            nested_mapped_box_level,
            unnested_to_nested,
            new_mapped_box_level,
            new_to_tag);
         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes applying overflow nesting map."
            << std::endl;
            tbox::plog << "Overflow nesting map:\n" << unnested_to_nested.format("", 3);
         }
         t_use_overflow_map->start();
         t_modify_connector->start();
         hier::MappingConnectorAlgorithm mca;
         mca.modify(tag_to_new,
            new_to_tag,
            unnested_to_nested,
            &new_mapped_box_level);
         t_modify_connector->stop();
         t_use_overflow_map->stop();
         if (d_barrier_and_time) {
            t_limit_overflow->stop();
         }
         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes finished applying overflow nesting map."
            << std::endl;
         }

         if (d_check_overflow_nesting) {
            if (d_print_steps) {
               tbox::plog << "GriddingAlgorithm::findRefinementBoxes checking overflow."
                          << std::endl;
            }
            bool locally_nested = false;
            bool nested = dlbg_edge_utils.baseNestsInHead(
                  &locally_nested,
                  new_mapped_box_level,
                  tag_mapped_box_level,
                  hier::IntVector::getZero(d_dim),
                  hier::IntVector::getZero(d_dim),
                  hier::IntVector::getZero(d_dim),
                  &d_hierarchy->getGridGeometry()->getDomainSearchTree());
            if (!nested) {
               TBOX_ERROR(
                  "Failed overflow nesting: new mapped_box_level does not nest in tagged mapped_box_level.\n"
                  << "Local nestedness = " << locally_nested << std::endl
                  << "tag_mapped_box_level:\n" << tag_mapped_box_level.format("", 2)
                  << "new_mapped_box_level:\n" << new_mapped_box_level.format("", 2)
                  << "tag_to_new:\n" << tag_to_new.format("", 2)
                  << "new_to_tag:\n" << new_to_tag.format("", 2));
            }
         }
      }

      if (d_enforce_proper_nesting) {

         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes: enforcing proper nesting\n";
         }
         if (d_barrier_and_time) {
            t_enforce_nesting->barrierAndStart();
         }

         hier::Connector unnested_to_nested;
         hier::BoxLevel nested_mapped_box_level(d_dim);

         makeProperNestingMap(
            nested_mapped_box_level,
            unnested_to_nested,
            new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            new_ln);
         if (d_print_steps) {
            tbox::plog << "GriddingAlgorithm::findRefinementBoxes applying proper nesting map.\n"
                       << "Proper nesting map:\n" << unnested_to_nested.format("", 3);
         }
         t_use_nesting_map->start();
         t_modify_connector->start();
         const hier::MappingConnectorAlgorithm mca;
         mca.modify(tag_to_new,
            new_to_tag,
            unnested_to_nested,
            &new_mapped_box_level);
         t_modify_connector->stop();
         t_use_nesting_map->stop();

         if (d_barrier_and_time) {
            t_enforce_nesting->stop();
         }

         if (tag_ln == d_base_ln && d_check_proper_nesting) {
            /*
             * Tag level will be regridded when we exit the current
             * recursion if tag_ln is not d_base_ln, so do not check
             * proper nesting in that case.
             *
             * Check that the new mapped_box_level nest in the tag
             * level (tag_ln).
             */
            hier::IntVector required_nesting(d_dim);
            if (tag_ln > 0) {
               required_nesting =
                  hier::IntVector(d_dim, d_hierarchy->getProperNestingBuffer(tag_ln));
            } else {
               required_nesting =
                  d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_dim);
            }
            bool locally_nests = false;
            const bool new_nests_in_tag =
               dlbg_edge_utils.baseNestsInHead(
                  &locally_nests,
                  new_mapped_box_level,
                  tag_mapped_box_level,
                  required_nesting,
                  hier::IntVector::getZero(d_dim),
                  hier::IntVector::getZero(d_dim),
                  &d_hierarchy->getGridGeometry()->getPeriodicDomainSearchTree());
            if (!new_nests_in_tag) {
               tbox::perr << "GriddingAlgorithm: new BoxLevel\n"
                          << "at ln=" << new_ln
                          << " does not properly nest in\n"
                          << "tag level at tag_ln=" << tag_ln
                          << " by the required nesting buffer of "
                          << required_nesting
                          << ".\nLocal nestingness: " << locally_nests
                          << ".\nWriting BoxLevels out to log file."
                          << std::endl;
               tbox::plog
               << "Proper nesting violation with new BoxLevel of\n"
               << new_mapped_box_level.format("N->", 2)
               << "Proper nesting violation with tag BoxLevel of\n"
               << tag_mapped_box_level.format("T->", 2);
               hier::BoxLevel external(d_dim);
               hier::Connector tmp_new_to_tag(
                  new_mapped_box_level,
                  tag_mapped_box_level,
                  required_nesting);
               oca.findOverlaps(tmp_new_to_tag);
               tbox::plog << "tmp_new_to_tag:\n" << tmp_new_to_tag.format("NT->", 3);
               hier::Connector new_to_external;
               dlbg_edge_utils.computeExternalParts(
                  external,
                  new_to_external,
                  tmp_new_to_tag,
                  -required_nesting,
                  d_hierarchy->getGridGeometry()->getDomainSearchTree());
               tbox::plog << "External parts:\n" << new_to_external.format("NE->", 3);
               TBOX_ERROR(
                  "Internal library error: Failed to produce proper nesting.");
            }
         }
      }

      if (d_extend_to_domain_boundary) {

         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes: extending nodes\n";
         }

         t_extend_to_domain_boundary->barrierAndStart();
         extendBoxesToDomainBoundary(
            new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            level->getPhysicalDomainArray(),
            extend_ghosts_in_tag_space);
         t_extend_to_domain_boundary->stop();
      }

      bool allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
         d_hierarchy->allowPatchesSmallerThanMinimumSize();

      // BTNG: these if-else blocks can be significantly simplified by factoring.
      if (!allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes: growing boxes\n";
         }
         t_extend_within_domain->start();
         growBoxesWithinNestingDomain(
            new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            smallest_box_to_refine,
            tag_ln);
         t_extend_within_domain->stop();
      } else {
         const hier::IntVector periodic_dirs(
            d_hierarchy->getGridGeometry()->getPeriodicShift(
               hier::IntVector::getOne(d_dim)));

         bool need_to_grow = false;
         hier::IntVector min_size(hier::IntVector::getOne(d_dim));
         for (int i = 0; i < d_dim.getValue(); i++) {
            if (periodic_dirs(i)) {
               need_to_grow = true;
               min_size(i) = smallest_box_to_refine(i);
            }
         }

         if (need_to_grow) {
            t_extend_within_domain->start();
            growBoxesWithinNestingDomain(
               new_mapped_box_level,
               tag_to_new,
               new_to_tag,
               min_size,
               tag_ln);
            t_extend_within_domain->stop();
         } else {
            /*
             * Had we need to grow, the growth would have shrunken the widths.
             * We must manually shrink the widths as if we used the growing map
             * (matching the result of an empty map).
             */
            tag_to_new.shrinkWidth(
               tag_to_new.getConnectorWidth() - smallest_box_to_refine);
            new_to_tag.shrinkWidth(
               new_to_tag.getConnectorWidth() - smallest_box_to_refine);
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      std::set<int> new_local_ids;
      const hier::BoxContainer& new_boxes = new_mapped_box_level.getBoxes();
      for (hier::BoxContainer::const_iterator new_itr = new_boxes.begin();
           new_itr != new_boxes.end(); ++new_itr) {
         new_local_ids.insert(new_itr->getId().getLocalId().getValue());
      }
      TBOX_ASSERT(static_cast<int>(new_local_ids.size()) == new_boxes.size());
#endif

      t_box_massage->stop();

      if (d_load_balance) {
         if (d_print_steps) {
            tbox::plog
            << "GriddingAlgorithm::findRefinementBoxes: load balancing\n";
         }

         t_load_balance->barrierAndStart();
         t_load_balance_setup->start();

         hier::IntVector patch_cut_factor(d_dim, 1);

         t_load_balance_setup->stop();

         d_load_balancer->loadBalanceBoxLevel(
            new_mapped_box_level,
            new_to_tag,
            tag_to_new,
            d_hierarchy,
            new_ln,
            new_to_tag, // FIXME: try using finer as the attractor.
            tag_to_new, // FIXME: try using finer as the attractor.
            smallest_patch_in_tag_space,
            largest_patch_in_tag_space,
            d_hierarchy->getDomainBoxLevel(),
            extend_ghosts_in_tag_space,
            patch_cut_factor);

         t_load_balance->stop();

         if (d_check_connectors) {
            tbox::plog << "GriddingAlgorithm checking new-tag" << std::endl;
            int errs = 0;
            if (oca.checkOverlapCorrectness(new_to_tag, false, true, true)) {
               ++errs;
               tbox::perr << "Error found in new_to_tag!\n";
            }
            if (oca.checkOverlapCorrectness(tag_to_new, false, true, true)) {
               ++errs;
               tbox::perr << "Error found in tag_to_new!\n";
            }
            if (new_to_tag.checkTransposeCorrectness(tag_to_new)) {
               ++errs;
               tbox::perr << "Error found in new-tag transpose!\n";
            }
            if (errs != 0) {
               TBOX_ERROR(
                  "Errors found after using load balance map."
                  << "new_mapped_box_level:\n" << new_mapped_box_level.format("", 2)
                  << "tag_mapped_box_level:\n" << tag_mapped_box_level.format("", 2)
                  << "new_to_tag:\n" << new_to_tag.format("", 2)
                  << "tag_to_new:\n" << tag_to_new.format("", 2));
            }
         }
      }

      if (d_sequentialize_patch_indices) {
         if (d_print_steps) {
            tbox::plog << "GriddingAlgorithm begin sorting nodes." << std::endl;
         }
         renumberBoxes(new_mapped_box_level,
            tag_to_new,
            new_to_tag,
            false,
            true);
         if (d_print_steps) {
            tbox::plog << "GriddingAlgorithm end sorting nodes." << std::endl;
         }
      }

      /*
       * Add periodic image Boxes to new_mapped_box_level and add edges
       * incident on those nodes.
       */
      const hier::Connector& tag_to_tag =
         d_hierarchy->getConnector(tag_ln, tag_ln);
      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm begin adding periodic images."
                    << std::endl;
      }
      dlbg_edge_utils.addPeriodicImagesAndRelationships(
         new_mapped_box_level,
         new_to_tag,
         tag_to_new,
         d_hierarchy->getGridGeometry()->getDomainSearchTree(),
         tag_to_tag);
      if (d_print_steps) {
         tbox::plog << "GriddingAlgorithm begin adding periodic images."
                    << std::endl;
      }

      if (d_check_connectors) {
         TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag) == 0);
         TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new) == 0);
      }

      /*
       * We have been working with new_mapped_box_level in the
       * tag_mapped_box_level's index space.  Now, refine it so we can
       * build the new level.
       */
      refineNewBoxLevel(new_mapped_box_level,
         tag_to_new,
         new_to_tag,
         ratio);

      if (d_check_connectors) {
         TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag) == 0);
         TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new) == 0);
      }

   } else {

      /*
       * On return, new_mapped_box_level should be initialized if we
       * generated boxes, deallocated if we didnt.
       */
      new_mapped_box_level.clear();

   }

   d_hierarchy->getMPI().Barrier();

   if (d_barrier_and_time) {
      t_find_refinement->stop();
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
GriddingAlgorithm::renumberBoxes(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   bool sort_by_corners,
   bool sequentialize_global_indices) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, new_mapped_box_level);

   t_sort_nodes->start();

   hier::MappingConnectorAlgorithm mca;

   hier::Connector sorting_map;
   hier::BoxLevel seq_mapped_box_level(d_dim);
   hier::BoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.makeSortingMap(
      seq_mapped_box_level,
      sorting_map,
      new_mapped_box_level,
      sort_by_corners,
      sequentialize_global_indices);

   t_modify_connector->start();
   mca.modify(tag_to_new,
      new_to_tag,
      sorting_map,
      &new_mapped_box_level);
   t_modify_connector->stop();

   t_sort_nodes->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
GriddingAlgorithm::refineNewBoxLevel(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const hier::IntVector& ratio) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, new_mapped_box_level, ratio);

   new_mapped_box_level.refineBoxes(new_mapped_box_level,
      ratio,
      new_mapped_box_level.getRefinementRatio()*ratio);
   new_mapped_box_level.finalize();

   const hier::IntVector& new_to_tag_width =
      ratio * new_to_tag.getConnectorWidth();
   new_to_tag.setBase(new_mapped_box_level);
   new_to_tag.setWidth(new_to_tag_width, true);

   tag_to_new.setHead(new_mapped_box_level, true);
   tag_to_new.refineLocalNeighbors(ratio);
#if 0
   const hier::OverlapConnectorAlgorithm oca;
   TBOX_ASSERT(oca.checkOverlapCorrectness(tag_to_new) == 0);
   TBOX_ASSERT(oca.checkOverlapCorrectness(new_to_tag) == 0);
#endif
}

/*
 *************************************************************************
 * Extend Boxes to domain boundary if they are too close.
 *************************************************************************
 */

void
GriddingAlgorithm::extendBoxesToDomainBoundary(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const tbox::Array<hier::BoxContainer>& physical_domain_array,
   const hier::IntVector& extend_ghosts) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, new_mapped_box_level);

   tbox::SAMRAI_MPI mpi(new_mapped_box_level.getMPI());

   /*
    * Extend boxes to domain boundary if they are too close to one.  I
    * think there is no need to modify connectivities when a box is
    * extended to the domain boundary.  There is no increased overlap
    * to a finer level, because the finer level already nests in the
    * nodes being extended.  There is increased overlap with the
    * coarser level, but it should not result in additional relationships,
    * because there would not be enough room for an unseen coarse Box
    * to live in the small gap across which the Box is being extended.
    */
   const hier::BoxContainer& before_nodes =
      new_mapped_box_level.getBoxes();

   hier::BoxLevel after_mapped_box_level(d_dim);
   after_mapped_box_level.initialize(
      new_mapped_box_level.getRefinementRatio(),
      new_mapped_box_level.getGridGeometry(),
      new_mapped_box_level.getMPI());

   hier::Connector before_to_after(
      new_mapped_box_level,
      after_mapped_box_level,
      extend_ghosts);
   before_to_after.setConnectorType(hier::Connector::MAPPING);

   for (hier::BoxContainer::const_iterator nn = before_nodes.begin();
        nn != before_nodes.end(); ++nn) {
      const hier::Box& before_mapped_box = *nn;
      hier::Box after_mapped_box = before_mapped_box;
      hier::BoxUtilities::extendBoxToDomainBoundary(
         after_mapped_box,
         physical_domain_array[before_mapped_box.getBlockId().getBlockValue()],
         extend_ghosts);
      after_mapped_box_level.addBox(after_mapped_box);
      before_to_after.insertLocalNeighbor(
         after_mapped_box,
         before_mapped_box.getId());
   }

   const hier::MappingConnectorAlgorithm mca;
   mca.modify(tag_to_new,
      new_to_tag,
      before_to_after,
      &new_mapped_box_level);

#if 0
   TBOX_WARNING("Performing extensive error checking due to using new code!");
   TBOX_ASSERT(new_to_tag.checkOverlapCorrectness() == 0);
   TBOX_ASSERT(tag_to_new.checkOverlapCorrectness() == 0);
#endif
}

/*
 *************************************************************************
 * Make a map that can be used to enforce overflow nesting.
 *************************************************************************
 */

void
GriddingAlgorithm::makeOverflowNestingMap(
   hier::BoxLevel& nested_mapped_box_level,
   hier::Connector& unnested_to_nested,
   const hier::BoxLevel& unnested_mapped_box_level,
   const hier::Connector& unnested_to_reference) const
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(unnested_mapped_box_level);
#endif

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim,
      unnested_mapped_box_level,
      nested_mapped_box_level);

   t_make_overflow_map->start();

   hier::BoxLevel violator_mapped_box_level(d_dim);
   hier::Connector unnested_to_violator;
   if (d_print_steps) {
      tbox::plog << " GriddingAlgorithm::makeOverflowNestingMap computing external parts."
                 << std::endl;
   }
   t_make_overflow_map_compute->start();
   hier::BoxLevelConnectorUtils edge_utils;
   t_compute_external_parts->start();
   edge_utils.computeExternalParts(
      violator_mapped_box_level,
      unnested_to_violator,
      unnested_to_reference,
      hier::IntVector::getZero(d_dim),
      d_hierarchy->getGridGeometry()->getDomainSearchTree());
   t_compute_external_parts->stop();
   t_make_overflow_map_compute->stop();

   TBOX_ASSERT(unnested_to_violator.isLocal());

   if (d_print_steps) {
      tbox::plog << " GriddingAlgorithm::makeOverflowNestingMap making remainer map." << std::endl;
   }
   t_make_overflow_map_convert->start();
   edge_utils.makeRemainderMap(
      nested_mapped_box_level,
      unnested_to_nested,
      unnested_to_violator);
   t_make_overflow_map_convert->stop();
   t_make_overflow_map->stop();

   if (d_print_steps) {
      tbox::plog << " GriddingAlgorithm::makeOverflowNestingMap finished." << std::endl;
   }
}

/*
 *************************************************************************
 * Make a map that can be used to enforce proper nesting.
 * @param unnested_ln Level number for refinenement ratio of unnested_mapped_box_level.
 *************************************************************************
 */

void
GriddingAlgorithm::makeProperNestingMap(
   hier::BoxLevel& nested_mapped_box_level,
   hier::Connector& unnested_to_nested,
   const hier::BoxLevel& unnested_mapped_box_level,
   const hier::Connector& hierarchy_to_unnested,
   const hier::Connector& unnested_to_hierarchy,
   const int unnested_ln) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim,
      unnested_mapped_box_level,
      nested_mapped_box_level,
      nested_mapped_box_level);

   t_make_nesting_map->start();

   hier::Connector unnested_to_violator;
   hier::BoxLevel violator(d_dim);
   t_make_nesting_map_compute->start();
   computeNestingViolator(
      violator,
      unnested_to_violator,
      unnested_mapped_box_level,
      unnested_to_hierarchy,
      hierarchy_to_unnested,
      unnested_ln - 1);
   t_make_nesting_map_compute->stop();

   /*
    * unnested_to_violator is the Connector from the nodes
    * that violate nesting to their violating parts.
    * Convert it to the mapping from unnested to nested.
    */
   const hier::BoxLevelConnectorUtils edge_utils;
   t_make_nesting_map_convert->start();
   edge_utils.makeRemainderMap(
      nested_mapped_box_level,
      unnested_to_nested,
      unnested_to_violator);
   t_make_nesting_map_convert->stop();

   t_make_nesting_map->stop();
}

/*
 *************************************************************************
 * Make a map from a BoxLevel to parts of that BoxLevel
 * that violate proper nesting.
 *
 * The violating Boxes are found by comparing candidate
 * Boxes to d_to_nesting_complement's head BoxLevel.
 * Boxes inside the nesting complement violate nesting.
 *************************************************************************
 */
void
GriddingAlgorithm::computeNestingViolator(
   hier::BoxLevel& violator,
   hier::Connector& candidate_to_violator,
   const hier::BoxLevel& candidate,
   const hier::Connector& candidate_to_hierarchy,
   const hier::Connector& hierarchy_to_candidate,
   const int tag_ln) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, candidate, violator);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   // Check requirements on arguments.
   TBOX_ASSERT(candidate_to_hierarchy.getRatio() ==
      hier::IntVector::getOne(d_dim));
   TBOX_ASSERT(hierarchy_to_candidate.getRatio() ==
      hier::IntVector::getOne(d_dim));

   const hier::BaseGridGeometry& grid_geometry(
      *d_hierarchy->getGridGeometry());

   t_compute_nesting_violator->start();

   const hier::BoxLevelConnectorUtils edge_utils;
   const hier::OverlapConnectorAlgorithm oca;

   const hier::BoxContainer& candidate_mapped_boxes = candidate.getBoxes();
   /*
    * Bridge candidate to d_nesting_complement.  Any part of the
    * candidate BoxLevel that overlaps d_nesting_complement
    * violates nesting.
    */
   hier::Connector candidate_to_complement;
   oca.bridge(candidate_to_complement,
      candidate_to_hierarchy,
      d_to_nesting_complement[tag_ln],
      d_from_nesting_complement[tag_ln],
      hierarchy_to_candidate);

   edge_utils.computeInternalParts(
      violator,
      candidate_to_violator,
      candidate_to_complement,
      zero_vector,
      grid_geometry.getDomainSearchTree());
   /*
    * Above step ignored the domain complement components of nesting
    * definition (by necessity).  Where the candidate falls outside
    * the domain, it violates nesting.
    */

   hier::BoxContainer refined_domain_search_tree(
      d_hierarchy->getGridGeometry()->getDomainSearchTree());
   refined_domain_search_tree.refine(candidate.getRefinementRatio());
   refined_domain_search_tree.makeTree(&grid_geometry);

   for (hier::BoxContainer::const_iterator ni = candidate_mapped_boxes.begin();
        ni != candidate_mapped_boxes.end(); ++ni) {
      const hier::Box& cmb = *ni;
      hier::BoxContainer addl_violators(cmb);
      addl_violators.removeIntersections(
         candidate.getRefinementRatio(),
         refined_domain_search_tree);
      if (!addl_violators.isEmpty()) {
         /*
          * Non-periodic BoxId needed for NeighborhoodSet::find()
          */
         hier::BoxId cmb_non_per_id(cmb.getGlobalId(),
                                    hier::PeriodicId::zero());
         if (candidate_to_violator.hasNeighborSet(cmb_non_per_id)) {
            /*
             * Remove parts that we already know, through
             * candidate_to_violator, are non-nesting.  Leftovers are
             * non-nesting parts not found using
             * candidate_to_complement.
             */
            hier::Connector::NeighborhoodIterator base_box_itr =
               candidate_to_violator.makeEmptyLocalNeighborhood(cmb_non_per_id);
            hier::Connector::ConstNeighborhoodIterator current_violators =
               candidate_to_violator.find(cmb_non_per_id);
            for (hier::Connector::ConstNeighborIterator na = candidate_to_violator.begin(current_violators);
                 na != candidate_to_violator.end(current_violators) && !addl_violators.isEmpty();
                 ++na) {
               addl_violators.removeIntersections(*na);
            }
            if (!addl_violators.isEmpty()) {
               for (hier::BoxContainer::iterator bi(addl_violators);
                    bi != addl_violators.end(); ++bi) {
                  hier::BoxContainer::const_iterator new_violator = violator.addBox(
                        *bi, cmb.getBlockId());
                  candidate_to_violator.insertLocalNeighbor(*new_violator,
                     base_box_itr);
               }
            }
         }
      }
   }

   t_compute_nesting_violator->stop();
}

/*
 *************************************************************************
 * Precompute data used to define proper nesting.  Data is associated
 * with level number ln, to be used for constructing level number ln+1.
 *
 * Data computed: d_proper_nesting_complement[ln],
 * d_from_nesting_complement[ln], d_to_proper_nesting_complement[ln].
 *
 * If ln > d_base_ln, assume data at ln-1 is already set.
 *************************************************************************
 */

void
GriddingAlgorithm::computeProperNestingData(
   const int ln)
{
   TBOX_ASSERT(d_base_ln >= 0 && ln >= d_base_ln);

   const hier::BoxLevelConnectorUtils edge_utils;
   const hier::OverlapConnectorAlgorithm oca;

   if (ln == d_base_ln) {
      /*
       * At the base level, nesting domain is level d_base_ln,
       * shrunken by d_proper_nesting_buffer[d_base_ln].
       */
      hier::BoxLevel& proper_nesting_complement =
         d_proper_nesting_complement[ln];
      const hier::Connector& self_connector =
         d_hierarchy->getConnector(ln, ln);

      // This assert shoud pass due to GriddingAlgorithmConnectorWidthRequestor.
      TBOX_ASSERT( self_connector.getConnectorWidth() >=
                   hier::IntVector(d_dim, -d_hierarchy->getProperNestingBuffer(ln)) );

      edge_utils.computeExternalParts(
         proper_nesting_complement,
         d_to_nesting_complement[ln],
         self_connector,
         hier::IntVector(d_dim, -d_hierarchy->getProperNestingBuffer(ln)),
         d_hierarchy->getGridGeometry()->getDomainSearchTree());

      d_from_nesting_complement[ln].initializeToLocalTranspose(
         d_to_nesting_complement[ln]);

   } else {

      TBOX_ASSERT(d_to_nesting_complement[ln - 1].isFinalized());

      /*
       * How to build d_proper_nesting_complement[ln] and connect it to level ln:
       *
       * In the left column are the BoxLevels in the hierarchy.
       * In the right are their proper nesting complements for that level number.
       *
       *                   (new
       *                 Connector)
       *                     |
       *             Mapped  |   Proper
       *               box   |   nesting
       *             levels  | complements
       *             ======  | ===========
       *                     v
       *               ln <-----> ln
       *                ^        ^
       *                |       /
       *                |      /
       *                |     /
       * (existing      |    / <--(temporary
       *  Connector)--> |   /      Connector)
       *                |  /
       *                | /
       *                |/
       *                v
       *             ln-1 <-----> ln-1
       *                     ^
       *                     |
       *                 (existing
       *                  Connector)
       *
       * We have existing Connectors between ln and ln-1 and also between
       * level ln-1 and the nesting complement at ln-1.
       *
       * 1. Build the complement at ln (d_proper_nesting_complement[ln])
       *    from the complement at ln-1 (d_proper_nesting_complement[ln-1])
       *    by refining and growing the complement boxes.
       *
       * 2. Build the temporary Connector from level ln-1 to
       *    complements at ln by using the fact that the complement at ln
       *    is similar to the one at ln-1.
       *
       * 3. Bridge for the new Connector from level ln to the
       *    complement at ln, using the temporary Connector.
       */

      /*
       * 1. Build the complement at ln (d_proper_nesting_complement[ln])
       *    from the complement at ln-1 (d_proper_nesting_complement[ln-1]).
       */
      d_proper_nesting_complement[ln].initialize(
         d_hierarchy->getBoxLevel(ln)->getRefinementRatio(),
         d_hierarchy->getGridGeometry(),
         d_to_nesting_complement[ln - 1].getMPI());
      const hier::BoxContainer& lnm1_complement_mapped_boxes =
         d_proper_nesting_complement[ln - 1].getBoxes();
      for (hier::BoxContainer::const_iterator ni =
           lnm1_complement_mapped_boxes.begin();
           ni != lnm1_complement_mapped_boxes.end(); ++ni) {
         hier::Box tmp_mapped_box = *ni;
         TBOX_ASSERT(!tmp_mapped_box.isPeriodicImage());
         tmp_mapped_box.refine(d_hierarchy->getRatioToCoarserLevel(ln));
         tmp_mapped_box.grow(
            hier::IntVector(d_dim, d_hierarchy->getProperNestingBuffer(ln)));
         d_proper_nesting_complement[ln].addBox(tmp_mapped_box);
      }

      /*
       * 2. Temporarily connect level ln-1 and d_proper_nesting_complement[ln].
       */
      hier::Connector lnm1_to_ln_complement(
         *d_hierarchy->getBoxLevel(ln - 1),
         d_proper_nesting_complement[ln],
         d_to_nesting_complement[ln - 1].getConnectorWidth());
      for (hier::Connector::ConstNeighborhoodIterator ei =
              d_to_nesting_complement[ln - 1].begin();
           ei != d_to_nesting_complement[ln - 1].end(); ++ei) {
         for (hier::Connector::ConstNeighborIterator na =
                 d_to_nesting_complement[ln - 1].begin(ei);
              na != d_to_nesting_complement[ln - 1].end(ei); ++na) {
            hier::Box tmp_mapped_box = *na;
            tmp_mapped_box.refine(d_hierarchy->getRatioToCoarserLevel(ln));
            tmp_mapped_box.grow(
               hier::IntVector(d_dim, d_hierarchy->getProperNestingBuffer(ln)));
            lnm1_to_ln_complement.insertLocalNeighbor(tmp_mapped_box, *ei);
         }
      }
      hier::Connector ln_complement_to_lnm1(d_from_nesting_complement[ln - 1]);
      ln_complement_to_lnm1.setBase(d_proper_nesting_complement[ln]);
      ln_complement_to_lnm1.setHead(*d_hierarchy->getBoxLevel(ln - 1));
      ln_complement_to_lnm1.setWidth(
         d_from_nesting_complement[ln - 1].getConnectorWidth() *
            d_hierarchy->getRatioToCoarserLevel(ln),
         true);

      /*
       * 3. Bridge for Connector between level ln and d_proper_nesting_complement[ln].
       */
      oca.bridge(
         d_to_nesting_complement[ln],
         d_from_nesting_complement[ln],
         d_hierarchy->getConnector(ln, ln - 1),
         lnm1_to_ln_complement,
         ln_complement_to_lnm1,
         d_hierarchy->getConnector(ln - 1, ln),
         d_hierarchy->getRequiredConnectorWidth(ln - 1, ln));
   }
}

/*
 *************************************************************************
 * Make a mapping Connector that can be used to grow boxes within
 * nesting domain by the minimum amount needed to make all boxes in a
 * BoxLevel satisfy the min_size requirement.
 *
 * Apply the map.
 *************************************************************************
 */

void
GriddingAlgorithm::growBoxesWithinNestingDomain(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const hier::IntVector& min_size,
   const int tag_ln) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim,
      new_mapped_box_level,
      min_size);

   const hier::OverlapConnectorAlgorithm oca;

   const hier::BaseGridGeometry& grid_geometry(
      *d_hierarchy->getGridGeometry());
   const int nblocks = grid_geometry.getNumberBlocks();
   hier::IntVector current_min_size(d_dim, tbox::MathUtilities<int>::getMax());
   for (int bn = 0; bn < nblocks; ++bn) {
      current_min_size.min(new_mapped_box_level.getGlobalMinBoxSize(bn));
   }

   if (current_min_size >= min_size) {
      /*
       * No box growing is needed.  Just shrink the Connector widths
       * to mimic expected the side-effect of applying a map with a
       * width of min_size.  The code should give the same result
       * without this special bypass.
       */
      tag_to_new.shrinkWidth(tag_to_new.getConnectorWidth() - min_size);
      new_to_tag.shrinkWidth(new_to_tag.getConnectorWidth() - min_size);
      return;
   }

   const hier::BoxContainer& new_mapped_boxes =
      new_mapped_box_level.getBoxes();

   const hier::Connector& tag_to_nesting_complement =
      d_to_nesting_complement[tag_ln];

   const hier::Connector& nesting_complement_to_tag =
      d_from_nesting_complement[tag_ln];

   /*
    * Connect new_mapped_box_level to the nesting complement so we
    * know where it cannot exist.  Use a Connector width of zero
    * because we don't need it any bigger and we don't want
    * inter-block neighbors.  (The new Boxes are already confined to
    * their own blocks before entering this method, so a zero
    * Connector width eliminates inter-block neighbors.)
    */

   hier::Connector new_to_nesting_complement;
   oca.bridge(
      new_to_nesting_complement,
      new_to_tag,
      tag_to_nesting_complement,
      nesting_complement_to_tag,
      tag_to_new,
      hier::IntVector::getZero(d_dim));

   /*
    * Set up the empty grown_mapped_box_level to be populated as we
    * determine whether each box needs to be grown.
    */
   hier::BoxLevel grown_mapped_box_level(
      new_mapped_box_level.getRefinementRatio(),
      new_mapped_box_level.getGridGeometry(),
      new_mapped_box_level.getMPI());

   // Create the mapping Connector from new to grown.
   hier::Connector new_to_grown(
      new_mapped_box_level,
      grown_mapped_box_level,
      min_size);
   new_to_grown.setConnectorType(hier::Connector::MAPPING);

   hier::BoxContainer refined_domain_search_tree(
      grid_geometry.getDomainSearchTree());
   refined_domain_search_tree.refine(new_mapped_box_level.getRefinementRatio());
   refined_domain_search_tree.makeTree(&grid_geometry);

   std::vector<hier::Box> tmp_mapped_box_vector;
   tmp_mapped_box_vector.reserve(10);

   /*
    * Loop through the new Boxes and grow if needed.  For each
    * Box, determine a sufficient view of the domain where it
    * can grow into.  This view includes domain boxes overlapping the
    * Box minus parts removed to satisfy nesting requirements.
    */

   for (hier::BoxContainer::const_iterator ni = new_mapped_boxes.begin();
        ni != new_mapped_boxes.end(); ++ni) {
      const hier::Box& omb = *ni;
      TBOX_ASSERT(!omb.isPeriodicImage());

      if (omb.numberCells() <= min_size) {
         // This box does not need growing.
         grown_mapped_box_level.addBox(omb);
         continue;
      }

      hier::BoxContainer nesting_domain;

      refined_domain_search_tree.findOverlapBoxes(
         nesting_domain,
         omb,
         new_mapped_box_level.getRefinementRatio());

      if (new_to_nesting_complement.hasNeighborSet(omb.getId())) {
         hier::Connector::ConstNeighborhoodIterator neighbors =
            new_to_nesting_complement.find(omb.getId());
         for (hier::Connector::ConstNeighborIterator na = new_to_nesting_complement.begin(neighbors);
              na != new_to_nesting_complement.end(neighbors); ++na) {
            nesting_domain.removeIntersections(*na);
         }
      }

      hier::Box grown_mapped_box = omb;
      hier::BoxUtilities::growBoxWithinDomain(
         grown_mapped_box,
         nesting_domain,
         min_size);

      /*
       * If the box is grown, generate the mapping for it.  If not,
       * keep the old box and don't generate a mapping.
       */
      if (!omb.isSpatiallyEqual(grown_mapped_box)) {
         grown_mapped_box_level.addBox(grown_mapped_box);
         new_to_grown.insertLocalNeighbor(grown_mapped_box, omb.getId());
      } else {
         grown_mapped_box_level.addBox(omb);
      }

   }

   /*
    * Use the mapping Connector.
    */

   t_modify_connector->start();
   const hier::MappingConnectorAlgorithm mca;
   mca.modify(tag_to_new,
      new_to_tag,
      new_to_grown,
      &new_mapped_box_level);
   t_modify_connector->stop();
}

void
GriddingAlgorithm::getGriddingParameters(
   hier::IntVector& smallest_patch,
   hier::IntVector& smallest_box_to_refine,
   hier::IntVector& largest_patch,
   hier::IntVector& extend_ghosts,
   const int level_number,
   const bool for_building_finer) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(d_dim,
      smallest_patch,
      smallest_box_to_refine,
      largest_patch,
      extend_ghosts);

   TBOX_ASSERT((level_number >= 0) && (level_number < d_hierarchy->getMaxNumberOfLevels()));

   /*
    * Determine maximum ghost cell width needed over all variables
    * currently known to the patch descriptor, and set the smallest
    * patch size.  The maximum number of ghosts is multiplied by the
    * error coarsen ratio (which should always be 1 unless regridding
    * uses error estimation).  This assures that when levels are
    * coarsened during error estimation, the coarser level patches
    * will meet the ghost cell constraint.
    */
   bool allow_patches_smaller_than_ghost_width =
      d_hierarchy->allowPatchesSmallerThanGhostWidth();

   hier::IntVector max_ghosts(
      d_hierarchy->getPatchDescriptor()->getMaxGhostWidth(d_dim));
   max_ghosts = max_ghosts * d_tag_init_strategy->getErrorCoarsenRatio();
   smallest_patch = d_hierarchy->getSmallestPatchSize(level_number);
   if (!allow_patches_smaller_than_ghost_width) {
      smallest_patch.max(max_ghosts);
   } else {
      const hier::IntVector periodic_dirs(
         d_hierarchy->getGridGeometry()->getPeriodicShift(hier::IntVector::getOne(
               d_dim)));

      for (int i = 0; i < d_dim.getValue(); i++) {
         if (periodic_dirs(i)) {
            smallest_patch(i) =
               tbox::MathUtilities<int>::Max(smallest_patch(i), max_ghosts(i));
         }
      }
   }

   /*
    * Set largest patch size.
    */
   largest_patch = d_hierarchy->getLargestPatchSize(level_number);

   /*
    * Following if-check prevents changing a negative largest_patch
    * bacause TreeLoadBalancer interprets the non-negative value as
    * dissabling the upper limit on patch size.
    */
   if (largest_patch > hier::IntVector::getZero(d_dim)) {
      largest_patch.max(smallest_patch);
   }

   /*
    * Set the smallest box to refine based on the number of cells that
    * coarsened patches must accomodate to meet ghost cell needs of variables.
    * On the finest level, the smallest box to refine is the smallest patch.
    * On coarser levels, it is a function of the error coarsen ratio and
    * the ratio to the next finer level.
    *
    * If we are accessing gridding parameters for a level that is being
    * reconstructed, the smallest box to refine is not applicable so we
    * set it to -1 to indicate an invalid entry in case it is used.
    */
   if (for_building_finer) {

      smallest_box_to_refine = smallest_patch;

      /*
       * Shouldn't this division be rounded up because it
       * represents a coarsening?  BTNG.
       */
      smallest_box_to_refine /=
         d_hierarchy->getRatioToCoarserLevel(level_number);
      /*
       * den = ratio from level_number to Richardson-coarsened
       * version of level_number+1.
       */
      const hier::IntVector den(
         d_hierarchy->getRatioToCoarserLevel(level_number)
         / d_tag_init_strategy->getErrorCoarsenRatio());
      /*
       * sz = max ghosts on Richardson-coarsened level_number+1, as
       * seen on level_number.
       */
      const hier::IntVector sz(hier::IntVector::ceilingDivide(max_ghosts, den));
      smallest_box_to_refine.max(sz);

   } else {

      smallest_box_to_refine = hier::IntVector(d_dim, -1);

   }

   /*
    * Determine number of cells box may be extended to physical
    * domain boundary to accomodate ghost cells.
    */
   extend_ghosts = max_ghosts;

}

/*
 *************************************************************************
 *************************************************************************
 */

void
GriddingAlgorithm::warnIfDomainTooSmallInPeriodicDir() const
{
   const hier::PeriodicShiftCatalog* shift_catalog =
      hier::PeriodicShiftCatalog::getCatalog(d_dim);

   if (shift_catalog->isPeriodic()) {

      hier::IntVector periodic_shift(
         d_hierarchy->getGridGeometry()->getPeriodicShift(
            hier::IntVector::getOne(d_dim)));

      hier::IntVector domain_bounding_box_size(
         d_hierarchy->getDomainBoxLevel().
         getGlobalBoundingBox(0).numberCells());

      for (int ln = 0; ln < d_hierarchy->getNumberOfLevels(); ++ln) {

         if (ln > 0) {
            periodic_shift *= d_hierarchy->getRatioToCoarserLevel(ln);
            domain_bounding_box_size *= d_hierarchy->getRatioToCoarserLevel(ln);
         }

         hier::IntVector smallest_patch_size(d_dim);
         hier::IntVector largest_patch_size(d_dim);
         hier::IntVector extend_ghosts(d_dim);
         hier::IntVector smallest_box_to_refine(d_dim);
         // "false" argument: for_building_finer level = false
         getGriddingParameters(
            smallest_patch_size,
            smallest_box_to_refine,
            largest_patch_size,
            extend_ghosts,
            ln,
            false);

         for (int d = 0; d < d_dim.getValue(); ++d) {
            if (periodic_shift(d) > 0 &&
                domain_bounding_box_size(d) < smallest_patch_size(d)) {
               TBOX_WARNING("GriddingAlgorithm: domain bounding box size\n"
                  << domain_bounding_box_size << " is smaller\n"
                  << "than the smallest patch size "
                  << smallest_patch_size << " on level "
                  << ln << " in dimension " << d << "\n");
               break;
            }
         }

      }

   }
}

/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::printClassData(
   std::ostream& os) const
{
   os << "\nGriddingAlgorithm::printClassData..." << std::endl;
   os << "   static data members:" << std::endl;
   for (int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; d++) {
      os << "      (*s_tag_indx)[" << d << "] = "
         << (*s_tag_indx)[d] << std::endl;
      os << "      (*s_buf_tag_indx)[" << d << "] = "
         << (*s_buf_tag_indx)[d] << std::endl;
   }
   os << "GriddingAlgorithm: this = "
      << (GriddingAlgorithm *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_tag_init_strategy = "
      << d_tag_init_strategy.get() << std::endl;
   os << "d_box_generator = "
      << d_box_generator.get() << std::endl;
   os << "d_load_balancer = "
      << d_load_balancer.get() << std::endl;
   os << "d_load_balancer0 = "
      << d_load_balancer0.get() << std::endl;
   os << "d_tag = " << d_tag.get() << std::endl;
   os << "d_tag_indx = " << d_tag_indx << std::endl;
   os << "d_buf_tag_indx = " << d_buf_tag_indx << std::endl;
   os << "d_true_tag = " << d_true_tag << std::endl;
   os << "d_false_tag = " << d_false_tag << std::endl;

   int ln;

   os << "d_efficiency_tolerance..." << std::endl;
   for (ln = 0; ln < d_efficiency_tolerance.getSize(); ln++) {
      os << "    d_efficiency_tolerance[" << ln << "] = "
         << d_efficiency_tolerance[ln] << std::endl;
   }
   os << "d_combine_efficiency..." << std::endl;
   for (ln = 0; ln < d_combine_efficiency.getSize(); ln++) {
      os << "    d_combine_efficiency[" << ln << "] = "
         << d_combine_efficiency[ln] << std::endl;
   }
}

/*
 *************************************************************************
 *
 * Write out class version number and data members to database.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::putToDatabase(
   const boost::shared_ptr<tbox::Database>& db) const
{
   TBOX_ASSERT(db);

   db->putInteger("ALGS_GRIDDING_ALGORITHM_VERSION",
      ALGS_GRIDDING_ALGORITHM_VERSION);

   db->putInteger("d_true_tag", d_true_tag);
   db->putInteger("d_false_tag", d_false_tag);

   db->putDoubleArray("d_efficiency_tolerance", d_efficiency_tolerance);
   db->putDoubleArray("d_combine_efficiency", d_combine_efficiency);

   db->putBool("d_sequentialize_patch_indices", d_sequentialize_patch_indices);
}

/*
 *************************************************************************
 *
 * If simulation is not from restart, read data from input database.
 * Otherwise, override data members initialized from restart with
 * values in the input database.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::getFromInput(
   const boost::shared_ptr<tbox::Database>& db,
   bool is_from_restart)
{
   NULL_USE(is_from_restart);

   TBOX_ASSERT(db);

   d_check_overflow_nesting =
      db->getBoolWithDefault("check_overflow_nesting", d_check_overflow_nesting);
   d_check_proper_nesting =
      db->getBoolWithDefault("check_proper_nesting", d_check_proper_nesting);
   d_check_connectors =
      db->getBoolWithDefault("check_connectors", d_check_connectors);
   d_print_steps =
      db->getBoolWithDefault("print_steps", d_print_steps);
   d_log_metadata_statistics =
      db->getBoolWithDefault("log_metadata_statistics", d_log_metadata_statistics);

   /*
    * Read input for efficiency tolerance.
    */

   if (db->keyExists("efficiency_tolerance")) {
      tbox::Array<double> efficiency_tolerance = db->getDoubleArray("efficiency_tolerance");

      int ln;
      for (ln = 0;
           ln < efficiency_tolerance.getSize() && ln < d_hierarchy->getMaxNumberOfLevels();
           ++ln) {
         if ((efficiency_tolerance[ln] <= 0.0e0) ||
             (efficiency_tolerance[ln] >= 1.0e0)) {
            TBOX_ERROR(d_object_name << ":  "
                                     << "Key data `efficiency_tolerance' has values"
                                     << " out of range 0.0 < tol < 1.0.");
         }
         d_efficiency_tolerance[ln] = efficiency_tolerance[ln];
      }
      for ( ; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
         d_efficiency_tolerance[ln] = efficiency_tolerance.back();
      }

   }

   /*
    * Read input for combine efficiency.
    */

   if (db->keyExists("combine_efficiency")) {
      tbox::Array<double> combine_efficiency = db->getDoubleArray("combine_efficiency");

      int ln;
      for (ln = 0;
           ln < combine_efficiency.getSize() && ln < d_hierarchy->getMaxNumberOfLevels();
           ++ln) {
         if ((combine_efficiency[ln] <= 0.0e0) ||
             (combine_efficiency[ln] >= 1.0e0)) {
            TBOX_ERROR(
               d_object_name << ":  "
                             << "Key data `combine_efficiency' has values"
                             << " out of range 0.0 < tol < 1.0.");
         }
         d_combine_efficiency[ln] = combine_efficiency[ln];
      }
      for ( ; ln < d_hierarchy->getMaxNumberOfLevels(); ++ln) {
         d_combine_efficiency[ln] = combine_efficiency.back();
      }

   }

   d_proper_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels(),
      hier::BoxLevel(d_dim));
   d_to_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels());
   d_from_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels());

   std::string tmp_str;

   tmp_str = db->getStringWithDefault("check_nonrefined_tags",
         std::string("WARN"));
   d_check_nonrefined_tags = char(tolower(*tmp_str.c_str()));
   if (d_check_nonrefined_tags != 'i' &&
       d_check_nonrefined_tags != 'w' &&
       d_check_nonrefined_tags != 'e') {
      TBOX_ERROR("GriddingAlgorithm: input parameter check_nonrefined_tags\n"
         << "can only be \"IGNORE\", \"WARN\" or \"ERROR\"");
   }

   tmp_str = db->getStringWithDefault("check_overlapping_patches",
         std::string("IGNORE"));
   d_check_overlapping_patches = char(tolower(*tmp_str.c_str()));
   if (d_check_overlapping_patches != 'i' &&
       d_check_overlapping_patches != 'w' &&
       d_check_overlapping_patches != 'e') {
      TBOX_ERROR(
         "GriddingAlgorithm: input parameter check_overlapping_patches\n"
         << "can only be \"IGNORE\", \"WARN\" or \"ERROR\"");
   }

   tmp_str = db->getStringWithDefault("check_nonnesting_user_boxes",
         std::string("ERROR"));
   d_check_nonnesting_user_boxes = char(tolower(*tmp_str.c_str()));
   if (d_check_nonnesting_user_boxes != 'i' &&
       d_check_nonnesting_user_boxes != 'w' &&
       d_check_nonnesting_user_boxes != 'e') {
      TBOX_ERROR("GriddingAlgorithm: input parameter check_nonnesting_user_boxes\n"
         << "can only be \"IGNORE\", \"WARN\" or \"ERROR\"");
   }

   tmp_str = db->getStringWithDefault("check_boundary_proximity_violation",
         std::string("ERROR"));
   d_check_boundary_proximity_violation = char(tolower(*tmp_str.c_str()));
   if (d_check_boundary_proximity_violation != 'i' &&
       d_check_boundary_proximity_violation != 'w' &&
       d_check_boundary_proximity_violation != 'e') {
      TBOX_ERROR("GriddingAlgorithm: input parameter check_boundary_proximity_violation\n"
         << "can only be \"IGNORE\", \"WARN\" or \"ERROR\"");
   }

   d_sequentialize_patch_indices =
      db->getBoolWithDefault("sequentialize_patch_indices",
         d_sequentialize_patch_indices);

   d_enforce_proper_nesting =
      db->getBoolWithDefault("enforce_proper_nesting", d_enforce_proper_nesting);
   d_extend_to_domain_boundary =
      db->getBoolWithDefault("extend_to_domain_boundary",
         d_extend_to_domain_boundary);
   d_load_balance =
      db->getBoolWithDefault("load_balance", d_load_balance);

   d_barrier_and_time =
      db->getBoolWithDefault("barrier_and_time", d_barrier_and_time);
}

/*
 *************************************************************************
 *
 * Gets the database in the root database that corresponds to the object
 * name.  This method then checks to make sure that the version number
 * of the class is that same as the version number in the restart file.
 * If these values are equal, the data members are read in from the
 * restart database.
 *
 *************************************************************************
 */

void
GriddingAlgorithm::getFromRestart()
{
   boost::shared_ptr<tbox::Database> root_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in restart file.");
   }
   boost::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));

   int ver = db->getInteger("ALGS_GRIDDING_ALGORITHM_VERSION");
   if (ver != ALGS_GRIDDING_ALGORITHM_VERSION) {
      TBOX_ERROR(
         d_object_name << ":  "
                       << "Restart file version different than class version.");
   }

   d_proper_nesting_complement.resize(d_hierarchy->getMaxNumberOfLevels(),
      hier::BoxLevel(d_dim));

   d_efficiency_tolerance = db->getDoubleArray("d_efficiency_tolerance");
   d_combine_efficiency = db->getDoubleArray("d_combine_efficiency");

   d_sequentialize_patch_indices = db->getBool("d_sequentialize_patch_indices");

}

/*
 *************************************************************************
 * Log metadata statistics after generating a new level.
 *
 * Log the given level, its peer connector and if requested, the
 * connectors to the next finer and next coarser levels.  Connectors
 * logged will have unit width.
 *************************************************************************
 */
void
GriddingAlgorithm::logMetadataStatistics(
   const char *caller_name,
   int ln,
   bool log_fine_connector,
   bool log_coarse_connector) const
{
   const std::string name("L" + tbox::Utilities::levelToString(ln));
   const hier::BoxLevel &level = *d_hierarchy->getPatchLevel(ln)->getBoxLevel();
   hier::PersistentOverlapConnectors &poc = level.getPersistentOverlapConnectors();
   const hier::IntVector &one_vector = hier::IntVector::getOne(d_dim);

   tbox::plog << "GriddingAlgorithm::" << caller_name << " added " << name << ":\n"
              << level.format("\t",0)
              << name << " statistics:\n"
              << level.formatStatistics("\t");

   const hier::Connector &peer_conn = poc.findOrCreateConnector(level, one_vector, true);
   tbox::plog << "Peer connector:\n" << peer_conn.format("\t",0)
              << "Peer connector statistics:\n" << peer_conn.formatStatistics("\t");

   if ( log_fine_connector ) {
      const hier::BoxLevel &fine_level = *d_hierarchy->getPatchLevel(ln+1)->getBoxLevel();
      const hier::Connector &fine_conn = poc.findOrCreateConnector(fine_level, one_vector, true);
      tbox::plog << "Fine connector:\n" << fine_conn.format("\t",0)
                 << "Fine connector statistics:\n" << fine_conn.formatStatistics("\t");
   }

   if ( log_coarse_connector ) {
      const hier::BoxLevel &crse_level = *d_hierarchy->getPatchLevel(ln-1)->getBoxLevel();
      const hier::Connector &crse_conn = poc.findOrCreateConnector(crse_level, one_vector, true);
      tbox::plog << "Coarse connector:\n" << crse_conn.format("\t",0)
                 << "Coarse connector statistics:\n" << crse_conn.formatStatistics("\t");
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
GriddingAlgorithm::allocateTimers()
{
   /*
    * Timers:  for gathering performance information about box
    * calculus and other regridding operations.
    */
   t_load_balance = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::load_balance");
   t_load_balance0 = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::load_balance0");
   t_load_balance_setup = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::load_balance_setup");
   t_bdry_fill_tags_create = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bdry_fill_tags_create");
   t_make_coarsest = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeCoarsestLevel()");
   t_make_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeFinerLevel()");
   t_make_finer_setup = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeFinerLevel()_setup");
   t_make_finer_tagging = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeFinerLevel()_tagging");
   t_make_finer_create = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeFinerLevel()_create");
   t_regrid_all_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::regridAllFinerLevels()");
   t_regrid_finer_create = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::regridFinerLevel()_create");
   t_fill_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::fillTagsFromBoxLevel()");
   t_tag_cells_for_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::tag_cells_for_refinement");
   t_buffer_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bufferTagsOnLevel()");
   t_second_finer_tagging = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::second_finer_tagging");
   t_bdry_fill_tags_comm = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bdry_fill_tags_comm");
   t_find_refinement = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::findRefinementBoxes()");
   t_find_boxes_containing_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::find_boxes_containing_tags");
   t_enforce_nesting = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::enforce_nesting");
   t_make_nesting_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeProperNestingMap()");
   t_make_nesting_map_compute = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeProperNestingMap()_compute");
   t_make_nesting_map_convert = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeProperNestingMap()_convert");
   t_use_nesting_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::use_nesting_map");
   t_make_overflow_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeOverflowNestingMap()");
   t_make_overflow_map_compute = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeOverflowNestingMap()_compute");
   t_make_overflow_map_convert = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeOverflowNestingMap()_convert");
   t_use_overflow_map = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::use_overflow_map");
   t_compute_external_parts = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::compute_external_parts");
   t_compute_nesting_violator = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::computeNestingViolator()");
   t_extend_to_domain_boundary = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::extend_to_domain_boundary");
   t_extend_within_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::extend_within_domain");
   t_grow_boxes_within_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::grow_boxes_within_domain");
   t_sort_nodes = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::renumberBoxes()");
   t_find_new_to_new = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::find_new_to_new");
   t_bridge_new_to_new = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_new_to_new");
   t_bridge_new_to_coarser = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_new_to_coarser");
   t_bridge_new_to_finer = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_new_to_finer");
   t_bridge_new_to_old = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_new_to_old");
   t_bridge_links = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::bridge_links");
   t_modify_connector = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::modify_connector");
   t_make_domain = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeCoarsestLevel()_make_domain");
   t_get_balance = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::get_balance");
   t_use_balance = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::use_balance");
   t_make_new = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::makeCoarsestLevel()_make_new");
   t_process_error = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::process_error");
   t_reset_hier = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::reset_hierarchy_config");
   t_misc1 = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::misc1");
   t_misc2 = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::misc2");
   t_misc3 = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::misc3");
   t_misc4 = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::misc4");
   t_misc5 = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::misc5");
   t_limit_overflow = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::limit_overflow");
   t_box_massage = tbox::TimerManager::getManager()->
      getTimer("mesh::GriddingAlgorithm::box_massage");
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
