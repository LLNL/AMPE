/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   An AMR hierarchy of patch levels
 *
 ************************************************************************/

#ifndef included_hier_PatchHierarchy_C
#define included_hier_PatchHierarchy_C

#include "SAMRAI/hier/PatchHierarchy.h"

#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace hier {

const int PatchHierarchy::HIER_PATCH_HIERARCHY_VERSION = 3;

std::vector<const PatchHierarchy::ConnectorWidthRequestorStrategy *>
PatchHierarchy::s_class_cwrs;

tbox::StartupShutdownManager::Handler
PatchHierarchy::s_initialize_finalize_handler(
   PatchHierarchy::initializeCallback,
   0,
   0,
   PatchHierarchy::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *
 * Instantiate the patch hierarchy and set default values.
 * Initialize from restart if necessary.
 *
 *************************************************************************
 */

PatchHierarchy::PatchHierarchy(
   const std::string& object_name,
   const boost::shared_ptr<BaseGridGeometry>& geometry,
   const boost::shared_ptr<tbox::Database>& database,
   bool register_for_restart):

   d_dim(geometry->getDim()),

   d_max_levels(1),
   d_ratio_to_coarser(1, IntVector(d_dim, 1)),
   d_proper_nesting_buffer(d_max_levels - 1, 1),
   d_smallest_patch_size(1, IntVector(d_dim, 1)),
   d_largest_patch_size(1, IntVector(d_dim, tbox::MathUtilities<int>::getMax())),
   d_allow_patches_smaller_than_ghostwidth(false),
   d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps(false),
   d_self_connector_widths(),
   d_fine_connector_widths(),
   d_connector_widths_are_computed(false),
   d_individual_cwrs(),
   d_domain_mapped_box_level(d_dim)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(geometry);

   d_object_name = object_name;
   d_registered_for_restart = register_for_restart;
   d_number_levels = 0;
   d_grid_geometry = geometry;
   d_patch_descriptor = VariableDatabase::getDatabase()->
      getPatchDescriptor();
   d_patch_factory.reset(new PatchFactory);
   d_patch_level_factory.reset(new PatchLevelFactory);
   d_number_blocks = d_grid_geometry->getNumberBlocks();

   /*
    * Grab the physical domain (including periodic images) from the
    * grid geometry and set up domain data dependent on it.
    */
   d_domain_mapped_box_level.initialize(
      IntVector::getOne(d_dim),
      getGridGeometry(),
      tbox::SAMRAI_MPI::getSAMRAIWorld(),
      BoxLevel::GLOBALIZED);
   d_grid_geometry->computePhysicalDomain(d_domain_mapped_box_level,
      IntVector::getOne(d_dim));
   d_domain_mapped_box_level.finalize();

   d_individual_cwrs = s_class_cwrs;

   if (database) {
      getFromInput(database);
   } else {
      /*
       * Without input database, the default is single-level with no
       * patch size constraints.
       */
   }

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(d_object_name, this);
   }
}

/*
 **************************************************************************
 *
 * The destructor tells the tbox::RestartManager to remove this hierarchy
 * from the list of restart items and automatically deletes all
 * allocated resources through smart pointers and arrays.
 *
 **************************************************************************
 */

PatchHierarchy::~PatchHierarchy()
{
   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
   }
}

/*
 *************************************************************************
 * If simulation is not from restart, read data from input database.
 * Otherwise, override data members initialized from restart with
 * values in the input database.
 *************************************************************************
 */

void
PatchHierarchy::getFromInput(
   const boost::shared_ptr<tbox::Database>& db)
{
   TBOX_ASSERT(db);

   /*
    * Read input for maximum number of levels.
    */

   d_max_levels = db->getIntegerWithDefault("max_levels", d_max_levels);

   if (d_max_levels != int(d_ratio_to_coarser.size())) {
      d_ratio_to_coarser.resize(d_max_levels, d_ratio_to_coarser.back());
      d_smallest_patch_size.resize(d_max_levels, d_smallest_patch_size.back());
      d_largest_patch_size.resize(d_max_levels, d_largest_patch_size.back());
   }

   std::vector<std::string> level_names(d_max_levels, std::string("level_"));
   for (int ln = 0; ln < d_max_levels; ++ln) {
      level_names[ln] += tbox::Utilities::intToString(ln);
   }

   // Read in ratio_to_coarser.
   if (db->isDatabase("ratio_to_coarser")) {
      const boost::shared_ptr<tbox::Database> tmp_db(
         db->getDatabase("ratio_to_coarser"));
      for (int ln = 1; ln < d_max_levels; ++ln) {
         if (tmp_db->isInteger(level_names[ln])) {
            tmp_db->getIntegerArray(level_names[ln],
               &d_ratio_to_coarser[ln][0],
               d_dim.getValue());
         } else {
            d_ratio_to_coarser[ln] = d_ratio_to_coarser[ln - 1];
         }
      }
   }

   // Read in smallest_patch_size.
   if (db->isDatabase("smallest_patch_size")) {
      const boost::shared_ptr<tbox::Database> tmp_db(
         db->getDatabase("smallest_patch_size"));
      for (int ln = 0; ln < d_max_levels; ++ln) {
         if (tmp_db->isInteger(level_names[ln])) {
            tmp_db->getIntegerArray(level_names[ln],
               &d_smallest_patch_size[ln][0],
               d_dim.getValue());
         } else {
            d_smallest_patch_size[ln] = d_smallest_patch_size[ln - 1];
         }
      }
   }

   // Read in largest_patch_size.
   if (db->isDatabase("largest_patch_size")) {
      const boost::shared_ptr<tbox::Database> tmp_db(
         db->getDatabase("largest_patch_size"));
      for (int ln = 0; ln < d_max_levels; ++ln) {
         if (tmp_db->isInteger(level_names[ln])) {
            tmp_db->getIntegerArray(level_names[ln],
               &d_largest_patch_size[ln][0],
               d_dim.getValue());
         } else {
            d_largest_patch_size[ln] = d_largest_patch_size[ln - 1];
         }
      }
   }

   tbox::Array<int> proper_nesting_buffer(1, 1);
   if (db->isInteger("proper_nesting_buffer")) {
      proper_nesting_buffer = db->getIntegerArray("proper_nesting_buffer");
   }
   d_proper_nesting_buffer.clear();
   for (int ln = 0; ln < d_max_levels - 1; ++ln) {
      if (ln < proper_nesting_buffer.size()) {
         d_proper_nesting_buffer.push_back(int(proper_nesting_buffer[ln]));
      } else {
         d_proper_nesting_buffer.push_back(int(d_proper_nesting_buffer[ln - 1]));
      }
   }
   for (size_t ln = 0; ln < d_proper_nesting_buffer.size(); ln++) {
      if (d_proper_nesting_buffer[ln] < 0) {
         TBOX_ERROR(
            d_object_name << ":  "
                          << "Key data `proper_nesting_buffer' has values < 0.");
      }
      if (d_proper_nesting_buffer[ln] == 0) {
         TBOX_WARNING(
            d_object_name << ":  "
                          << "Using zero `proper_nesting_buffer' values.");
      }
   }

   d_allow_patches_smaller_than_ghostwidth =
      db->getBoolWithDefault("allow_patches_smaller_than_ghostwidth",
         d_allow_patches_smaller_than_ghostwidth);

   d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
      db->getBoolWithDefault(
         "allow_patches_smaller_than_minimum_size_to_prevent_overlaps",
         d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps);
   if (d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
      TBOX_WARNING(
         d_object_name << ":  "
                       << "Allowing patches smaller than the given "
                       << "smallest patch size.  Note:  If periodic "
                       << "boundary conditions are used, this flag is "
                       << "ignored in the periodic directions.");
   }

}

/*
 *************************************************************************
 * Adds a ConnectorWidthRequestorStrategy to be used when this
 * PatchHierarchy computes its required Connector width.
 *************************************************************************
 */
void
PatchHierarchy::registerConnectorWidthRequestor(
   const ConnectorWidthRequestorStrategy& cwrs)
{
   if (d_connector_widths_are_computed) {
      TBOX_ERROR("PatchHierarchy::registerConnectorWidthRequestor:\n"
         << "Registering a new ConnectorWidthRequestorStrategy after\n"
         << "calling getRequiredConnectorWidth() is disallowed because\n"
         << "it may cause getRequiredConnectorWidth() to return\n"
         << "conflicting requirements.");
   }

   size_t i;
   for (i = 0; i < d_individual_cwrs.size(); ++i) {
      if (d_individual_cwrs[i] == &cwrs) {
         break;
      }
   }
   if (i == d_individual_cwrs.size()) {
      d_individual_cwrs.push_back(&cwrs);
   }
}

/*
 *************************************************************************
 * Adds a ConnectorWidthRequestorStrategy to be automatically
 * registered with all PatchHierarchy objects during their
 * construction (if they are not constructed with the flag to bypass
 * the auto-registration mechanism).
 *************************************************************************
 */
void
PatchHierarchy::registerAutoConnectorWidthRequestorStrategy(
   const ConnectorWidthRequestorStrategy& cwrs)
{
   size_t i;
   for (i = 0; i < s_class_cwrs.size(); ++i) {
      if (s_class_cwrs[i] == &cwrs) {
         break;
      }
   }
   if (i == s_class_cwrs.size()) {
      s_class_cwrs.push_back(&cwrs);
   }
}

/*
 ***************************************************************************
 * Clear out static registry.
 ***************************************************************************
 */

void
PatchHierarchy::finalizeCallback()
{
   for (int i = 0; i < int(s_class_cwrs.size()); ++i) {
      s_class_cwrs[i] = NULL;
   }
   s_class_cwrs.clear();
   /*
    * Hopefully, reserving 0 will free memory, making memory checkers
    * happy.
    */
   s_class_cwrs.reserve(0);
}

/*
 *************************************************************************
 *************************************************************************
 */
IntVector
PatchHierarchy::getRequiredConnectorWidth(
   int base_ln,
   int head_ln) const
{
   TBOX_ASSERT(head_ln >= 0);
   TBOX_ASSERT(head_ln < d_max_levels);
   TBOX_ASSERT(base_ln >= 0);
   TBOX_ASSERT(base_ln < d_max_levels);

   if (!d_connector_widths_are_computed) {

      const IntVector& zero_vector(IntVector::getZero(d_dim));
      d_self_connector_widths.clear();
      d_self_connector_widths.insert(d_self_connector_widths.begin(),
         d_max_levels,
         zero_vector);
      d_fine_connector_widths.clear();
      if (d_max_levels > 1) {
         d_fine_connector_widths.insert(d_fine_connector_widths.begin(),
            d_max_levels - 1,
            zero_vector);
      }

      /*
       * Get the required widths satisfying all registered
       * ConnectorWidthRequestorStrategy objects.
       */

      std::vector<IntVector> self_connector_widths;
      std::vector<IntVector> fine_connector_widths;
      for (size_t i = 0; i < d_individual_cwrs.size(); ++i) {
         d_individual_cwrs[i]->computeRequiredConnectorWidths(
            self_connector_widths,
            fine_connector_widths,
            *this);
         TBOX_ASSERT(self_connector_widths.size() == static_cast<unsigned int>(d_max_levels));
         TBOX_ASSERT(fine_connector_widths.size() == static_cast<unsigned int>(d_max_levels - 1));
         for (int ln = 0; ln < d_max_levels; ++ln) {
            d_self_connector_widths[ln].max(self_connector_widths[ln]);
         }
         for (int ln = 0; ln < d_max_levels - 1; ++ln) {
            d_fine_connector_widths[ln].max(fine_connector_widths[ln]);
         }
      }

      /*
       * Make sure the self connector widths are at least as big as
       * the fine.  This is required because self Connectors at the
       * tag level is used to compute the fine Connectors.  This
       * requirement is due to the GriddingAlgorithm, so perhaps it
       * should be moved there!  On the other hand, GriddingAlgorithm
       * cannot be expected know about fine_connector_width
       * requirements of other width requestors.
       */
      for (int ln = 0; ln < d_max_levels - 1; ++ln) {
         d_self_connector_widths[ln].max(d_fine_connector_widths[ln]);
      }
      d_connector_widths_are_computed = true;
   }

   if (base_ln != head_ln) {
      if (head_ln == base_ln + 1) {
         // Width is for fine Connector.
         return d_fine_connector_widths[base_ln];
      } else if (base_ln == head_ln + 1) {
         // Width is for coarse Connector.
         return d_fine_connector_widths[head_ln] * d_ratio_to_coarser[base_ln];
      }
      TBOX_ERROR("PatchHierarchy::getRequiredConnectorWidth: base_ln and\n"
         << "head_ln should differ by at most 1.\n"
         << "base_ln=" << base_ln << "  head_ln=" << head_ln);
   }
   return d_self_connector_widths[base_ln];
}

/*
 *************************************************************************
 * Get the Connector between 2 levels in the hierarchy.
 * The levels must be the same or adjacent levels.
 * The Connector will have the width required by the hierarchy.
 *
 * This is primarily a convenience.  Its functionality can be duplicated
 * by getting the required Connector widths and getting the Connector
 * with that width from the PatchLevels' BoxLevels.
 *************************************************************************
 */
const Connector&
PatchHierarchy::getConnector(
   const int base_ln,
   const int head_ln) const
{
   TBOX_ASSERT(base_ln >= 0);
   TBOX_ASSERT(base_ln < d_number_levels);
   TBOX_ASSERT(head_ln >= 0);
   TBOX_ASSERT(head_ln < d_number_levels);
   const BoxLevel& base(*d_patch_levels[base_ln]->getBoxLevel());
   const BoxLevel& head(*d_patch_levels[head_ln]->getBoxLevel());
   const IntVector width(getRequiredConnectorWidth(base_ln, head_ln));
   const Connector& con = base.getPersistentOverlapConnectors().
      findConnector(head, width);
   return con;
}

/*
 *************************************************************************
 *
 * Create a copy of this patch hierarchy with each level refined by
 * the given ratio and return a pointer to it.
 *
 *************************************************************************
 */

boost::shared_ptr<PatchHierarchy>
PatchHierarchy::makeRefinedPatchHierarchy(
   const std::string& fine_hierarchy_name,
   const IntVector& refine_ratio,
   bool register_for_restart) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, refine_ratio);
   TBOX_ASSERT(!fine_hierarchy_name.empty());
   TBOX_ASSERT(fine_hierarchy_name != d_object_name);
   TBOX_ASSERT(refine_ratio > IntVector::getZero(refine_ratio.getDim()));

   boost::shared_ptr<BaseGridGeometry> fine_geometry(
      d_grid_geometry->makeRefinedGridGeometry(
         fine_hierarchy_name + "GridGeometry",
         refine_ratio,
         register_for_restart));

   PatchHierarchy* fine_hierarchy =
      new PatchHierarchy(fine_hierarchy_name,
         fine_geometry,
         boost::shared_ptr<tbox::Database>(),
         register_for_restart);

   // Set hierarchy parameters.

   fine_hierarchy->d_max_levels = d_max_levels;
   fine_hierarchy->d_ratio_to_coarser = d_ratio_to_coarser;
   fine_hierarchy->d_smallest_patch_size = d_smallest_patch_size;
   fine_hierarchy->d_largest_patch_size = d_largest_patch_size;
   fine_hierarchy->d_individual_cwrs = d_individual_cwrs;
   fine_hierarchy->d_proper_nesting_buffer = d_proper_nesting_buffer;
   fine_hierarchy->d_allow_patches_smaller_than_ghostwidth =
      d_allow_patches_smaller_than_ghostwidth;
   fine_hierarchy->d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
      d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps;

   for (int ln = 0; ln < d_number_levels; ln++) {

      // SGS TODOD why isn't this a ctor with more args and not all this setting of new_level?
      // What happened to ctor is initialization?
      boost::shared_ptr<PatchLevel> new_level(
         boost::make_shared<PatchLevel>(d_dim));
      new_level->setRefinedPatchLevel(d_patch_levels[ln],
         refine_ratio,
         fine_geometry);

      new_level->setLevelNumber(ln);
      new_level->setNextCoarserHierarchyLevelNumber(ln - 1);
      new_level->setLevelInHierarchy(true);
      new_level->setRatioToCoarserLevel(
         d_patch_levels[ln]->getRatioToCoarserLevel());
      if (ln >= fine_hierarchy->d_number_levels) {
         fine_hierarchy->d_number_levels = ln + 1;
         fine_hierarchy->d_patch_levels.resizeArray(d_number_levels);
      }
      fine_hierarchy->d_patch_levels[ln] = new_level;
   }

   return boost::shared_ptr<PatchHierarchy>(fine_hierarchy);

}

/*
 *************************************************************************
 *                                                                       *
 * Create a copy of this patch hierarchy with each level coarsened by    *
 * the given ratio and return a pointer to it.                           *
 *                                                                       *
 *************************************************************************
 */

boost::shared_ptr<PatchHierarchy>
PatchHierarchy::makeCoarsenedPatchHierarchy(
   const std::string& coarse_hierarchy_name,
   const IntVector& coarsen_ratio,
   bool register_for_restart) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, coarsen_ratio);
   TBOX_ASSERT(!coarse_hierarchy_name.empty());
   TBOX_ASSERT(coarse_hierarchy_name != d_object_name);
   TBOX_ASSERT(coarsen_ratio > IntVector::getZero(coarsen_ratio.getDim()));

   boost::shared_ptr<BaseGridGeometry> coarse_geometry(
      d_grid_geometry->makeCoarsenedGridGeometry(
         coarse_hierarchy_name + "GridGeometry",
         coarsen_ratio,
         register_for_restart));

   PatchHierarchy* coarse_hierarchy =
      new PatchHierarchy(coarse_hierarchy_name,
         coarse_geometry,
         boost::shared_ptr<tbox::Database>(),
         register_for_restart);

   // Set hierarchy parameters.

   coarse_hierarchy->d_max_levels = d_max_levels;
   coarse_hierarchy->d_ratio_to_coarser = d_ratio_to_coarser;
   coarse_hierarchy->d_smallest_patch_size = d_smallest_patch_size;
   coarse_hierarchy->d_largest_patch_size = d_largest_patch_size;
   coarse_hierarchy->d_individual_cwrs = d_individual_cwrs;
   coarse_hierarchy->d_proper_nesting_buffer = d_proper_nesting_buffer;

   for (int ln = 0; ln < d_number_levels; ln++) {
      boost::shared_ptr<PatchLevel> new_level(new PatchLevel(d_dim));
      new_level->setCoarsenedPatchLevel(d_patch_levels[ln],
         coarsen_ratio,
         coarse_geometry);
      new_level->setLevelNumber(ln);
      new_level->setNextCoarserHierarchyLevelNumber(ln - 1);
      new_level->setLevelInHierarchy(true);
      new_level->setRatioToCoarserLevel(
         d_patch_levels[ln]->getRatioToCoarserLevel());
      if (ln >= coarse_hierarchy->d_number_levels) {
         coarse_hierarchy->d_number_levels = ln + 1;
         coarse_hierarchy->d_patch_levels.resizeArray(d_number_levels);
      }
      coarse_hierarchy->d_patch_levels[ln] = new_level;
   }

   return boost::shared_ptr<PatchHierarchy>(coarse_hierarchy);

}

/*
 *************************************************************************
 *                                                                       *
 * Create a new patch level in the hierarchy.                            *
 *                                                                       *
 *************************************************************************
 */

void
PatchHierarchy::makeNewPatchLevel(
   const int ln,
   const BoxLevel& new_mapped_box_level)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, new_mapped_box_level);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ln >= 0);
   for (int i = 0; i < d_dim.getValue(); i++) {
      TBOX_ASSERT(
         new_mapped_box_level.getRefinementRatio() > IntVector::getZero(d_dim));
   }
#endif

   /*
    * Make sure the level conforms to certain parameters preset
    * for the hierarchy.  We are not (yet) checking everything we
    * should.
    */
   if (ln >= d_max_levels) {
      TBOX_ERROR("PatchHierarchy::makeNewPatchLevel: Cannot make\n"
         << "level " << ln << " in a PatchHierarchy with a\n"
         << "max of " << d_max_levels << ".\n"
         << "Use setMaxNumberOfLevels() to change the max.\n");
   }
   if (ln > 0) {
      const IntVector expected_ratio(
         d_ratio_to_coarser[ln] * (d_patch_levels[ln - 1]->getRatioToLevelZero()));
      if (ln > 0 &&
          (new_mapped_box_level.getRefinementRatio() != expected_ratio)) {
         TBOX_ERROR("PatchHierarchy::makeNewPatchLevel: patch level "
            << ln << " has refinement ratio "
            << new_mapped_box_level.getRefinementRatio()
            << ", it should be " << expected_ratio);
      }
   }

   if (ln >= d_number_levels) {
      d_number_levels = ln + 1;
      d_patch_levels.resizeArray(d_number_levels);
   }

   d_patch_levels[ln] = d_patch_level_factory->allocate(
         new_mapped_box_level,
         d_grid_geometry, d_patch_descriptor, d_patch_factory);
   d_patch_levels[ln]->getBoxLevel()->cacheGlobalReducedData();

   d_patch_levels[ln]->setLevelNumber(ln);
   d_patch_levels[ln]->setNextCoarserHierarchyLevelNumber(ln - 1);
   d_patch_levels[ln]->setLevelInHierarchy(true);

   if ((ln > 0) && d_patch_levels[ln - 1]) {
      d_patch_levels[ln]->setRatioToCoarserLevel(
         d_patch_levels[ln]->getRatioToLevelZero()
         / (d_patch_levels[ln - 1]->getRatioToLevelZero()));
   }

}

/*
 *************************************************************************
 *                                                                       *
 * Remove the specified patch level from the hierarchy.                  *
 *                                                                       *
 *************************************************************************
 */

void
PatchHierarchy::removePatchLevel(
   const int l)
{
   TBOX_ASSERT((l >= 0) && (l < d_number_levels));

   d_patch_levels[l].reset();
   if (d_number_levels == l + 1) {
      d_number_levels--;
   }
}

/*
 *************************************************************************
 *
 * Writes out the class version number and the number of levels in the
 * hierarchy and has each patch_level write itself out.
 * The database keys for the patch levels are given by
 * "level#" where # is the level number for the patch_level.
 * The patchdata that are written to the database are determined by
 * which those bits in the VariableDatabase restart table.
 *
 * Asserts that the database pointer passed in is not NULL.
 *
 *************************************************************************
 */

void
PatchHierarchy::putToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   putToDatabase(database,
      VariableDatabase::getDatabase()->getPatchDataRestartTable());
}

/*
 *************************************************************************
 *
 * Writes out the class version number and the number of levels in the
 * hierarchy and has each patch_level write itself out.
 * The database keys for the patch levels are given by
 * "level#" where # is the level number for the patch_level.
 * The patchdata that are written to the database are determined by
 * which those bits in the specified ComponentSelector that are
 * set.
 *
 * Asserts that the database pointer passed in is not NULL.
 *
 *************************************************************************
 */

void
PatchHierarchy::putToDatabase(
   const boost::shared_ptr<tbox::Database>& database,
   const ComponentSelector& patchdata_write_table) const
{
   TBOX_ASSERT(database);

   database->putInteger("HIER_PATCH_HIERARCHY_VERSION",
      HIER_PATCH_HIERARCHY_VERSION);
   database->putInteger("d_number_levels", d_number_levels);

   std::vector<std::string> level_names(d_max_levels);
   const std::string prefix("level_");
   for (int ln = 0; ln < d_max_levels; ++ln) {
      level_names[ln] = prefix + tbox::Utilities::levelToString(ln);
   }

   /*
    * Write hierarchy parameters.
    */
   database->putInteger("d_max_levels", d_max_levels);
   database->putBool("d_allow_patches_smaller_than_ghostwidth",
      d_allow_patches_smaller_than_ghostwidth);
   database->putBool("d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps",
      d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps);
   if (d_max_levels > 1) {
      database->putIntegerArray("d_proper_nesting_buffer",
         &d_proper_nesting_buffer[0],
         static_cast<int>(d_proper_nesting_buffer.size()));
   }

   boost::shared_ptr<tbox::Database> ratio_to_coarser_db(
      database->putDatabase("d_ratio_to_coarser"));
   for (unsigned int ln = 0; ln < d_ratio_to_coarser.size(); ln++) {
      const int* tmp_array = &d_ratio_to_coarser[ln][0];
      ratio_to_coarser_db->putIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> largest_patch_db(
      database->putDatabase("d_largest_patch_size"));
   for (unsigned int ln = 0; ln < d_largest_patch_size.size(); ln++) {
      const int* tmp_array = &d_largest_patch_size[ln][0];
      largest_patch_db->putIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> smallest_patch_db(
      database->putDatabase("d_smallest_patch_size"));
   for (unsigned int ln = 0; ln < d_smallest_patch_size.size(); ln++) {
      const int* tmp_array = &d_smallest_patch_size[ln][0];
      smallest_patch_db->putIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> self_connector_widths_db(
      database->putDatabase("d_self_connector_widths"));
   for (unsigned int ln = 0; ln < d_self_connector_widths.size(); ln++) {
      int* tmp_array = &d_self_connector_widths[ln][0];
      self_connector_widths_db->putIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> fine_connector_widths_db(
      database->putDatabase("d_fine_connector_widths"));
   for (unsigned int ln = 0; ln < d_fine_connector_widths.size(); ln++) {
      int* tmp_array = &d_fine_connector_widths[ln][0];
      fine_connector_widths_db->putIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   for (int i = 0; i < d_number_levels; i++) {

     boost::shared_ptr<tbox::Database> level_database(
         database->putDatabase(level_names[i]));

      d_patch_levels[i]->putUnregisteredToDatabase(
         level_database,
         patchdata_write_table);
   }
}

/*
 *************************************************************************
 *
 * Gets the database in the root database that corresponds to the object
 * name.  This method then checks the class version against restart
 * file version.  If they match, it creates each hierarchy level and
 * reads in the level data.   The number of levels read from restart is
 * the minimum of the argument max levels and the number of levels in
 * the restart file.
 *
 *************************************************************************
 */
void
PatchHierarchy::getFromRestart()
{

  boost::shared_ptr<tbox::Database> restart_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!restart_db->isDatabase(d_object_name)) {
      TBOX_ERROR("PatchHierarchy::getFromRestart() error...\n"
         << "   Restart database with name "
         << d_object_name << " not found in restart file" << std::endl);
   }
   boost::shared_ptr<tbox::Database> database(
      restart_db->getDatabase(d_object_name));

   int ver = database->getInteger("HIER_PATCH_HIERARCHY_VERSION");
   if (ver != HIER_PATCH_HIERARCHY_VERSION) {
      TBOX_ERROR("PatchHierarchy::getFromRestart error...\n"
         << "  object name = " << d_object_name
         << " : Restart file version different than class version" << std::endl);
   }

   getFromDatabase(
      database,
      VariableDatabase::getDatabase()->getPatchDataRestartTable());
}

void
PatchHierarchy::getFromDatabase(
   const boost::shared_ptr<tbox::Database>& database,
   const ComponentSelector& component_selector)
{
   TBOX_ASSERT(database);

   d_number_levels = database->getInteger("d_number_levels");
   if (d_number_levels <= 0) {
      TBOX_ERROR("PatchHierarchy::getFromDatabase error ...\n"
         << "  object name = " << d_object_name
         << " : `d_number_levels' is <= zero in restart file");
   }

   d_number_levels = tbox::MathUtilities<int>::Min(d_number_levels,
         d_max_levels);

   d_patch_levels.resizeArray(d_number_levels);

   std::vector<std::string> level_names(d_max_levels);
   const std::string prefix("level_");
   for (int ln = 0; ln < d_max_levels; ++ln) {
      level_names[ln] = prefix + tbox::Utilities::levelToString(ln);
   }

   /*
    * Read hierarchy paremeters.
    */
   if (d_max_levels > 0) {
      d_max_levels = tbox::MathUtilities<int>::Min(
            d_max_levels, database->getInteger("d_max_levels"));
   } else {
      d_max_levels = database->getInteger("d_max_levels");
   }

   d_allow_patches_smaller_than_ghostwidth = database->getBool(
         "d_allow_patches_smaller_than_ghostwidth");
   d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps = database->getBool(
         "d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps");
   d_proper_nesting_buffer.resize(d_max_levels - 1, 0);
   if (d_max_levels > 1) {
      database->getIntegerArray("d_proper_nesting_buffer",
         &d_proper_nesting_buffer[0],
         d_max_levels - 1);
   }

   boost::shared_ptr<tbox::Database> ratio_to_coarser_db(
      database->getDatabase("d_ratio_to_coarser"));
   for (unsigned int ln = 0; ln < d_ratio_to_coarser.size(); ln++) {
      int* tmp_array = &d_ratio_to_coarser[ln][0];
      ratio_to_coarser_db->getIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> largest_patch_db(
      database->getDatabase("d_largest_patch_size"));
   for (unsigned int ln = 0; ln < d_largest_patch_size.size(); ln++) {
      int* tmp_array = &d_largest_patch_size[ln][0];
      largest_patch_db->getIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> smallest_patch_db(
      database->getDatabase("d_smallest_patch_size"));
   for (unsigned int ln = 0; ln < d_smallest_patch_size.size(); ln++) {
      int* tmp_array = &d_smallest_patch_size[ln][0];
      smallest_patch_db->getIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> self_connector_widths_db(
      database->getDatabase("d_self_connector_widths"));
   for (unsigned int ln = 0; ln < d_self_connector_widths.size(); ln++) {
      int* tmp_array = &d_self_connector_widths[ln][0];
      self_connector_widths_db->getIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> fine_connector_widths_db(
      database->getDatabase("d_fine_connector_widths"));
   for (unsigned int ln = 0; ln < d_fine_connector_widths.size(); ln++) {
      int* tmp_array = &d_fine_connector_widths[ln][0];
      fine_connector_widths_db->getIntegerArray(level_names[ln], tmp_array,
         d_dim.getValue());
   }

   for (int i = 0; i < d_number_levels; i++) {

     boost::shared_ptr<tbox::Database> level_database(
         database->getDatabase(level_names[i]));

      d_patch_levels[i] = d_patch_level_factory->allocate(
            level_database,
            d_grid_geometry,
            d_patch_descriptor,
            component_selector,
            d_patch_factory,
            false);
   }
   /*
    * Compute Connectors.
    * BTNG FIXME: This should be replaced by writing edges to
    * restart and reading them back.
    */
   for (int i = 0; i < d_number_levels; ++i) {
      for (int j = i - 1; j <= i + 1; ++j) {
         if (j >= 0 && j < d_number_levels) {
            const IntVector connector_width = getRequiredConnectorWidth(i, j);
            d_patch_levels[i]->getBoxLevel()->
            getPersistentOverlapConnectors().findOrCreateConnector(
               *d_patch_levels[j]->getBoxLevel(),
               connector_width);
         }
      }
   }

}

int
PatchHierarchy::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int depth)
{
   int totl_npatches = 0;
   int totl_ncells = 0;
   int nlevels = getNumberOfLevels();
   os << border << "Domain of hierarchy:\n" << d_domain_mapped_box_level.format(border, 2) << '\n'
      << border << "Number of levels = " << nlevels << '\n';
   if (depth > 0) {
      int ln;
      for (ln = 0; ln < nlevels; ++ln) {
         os << border << "Level " << ln << '/' << nlevels << "\n";
         boost::shared_ptr<PatchLevel> level(getPatchLevel(ln));
         level->recursivePrint(os, border + "\t", depth - 1);
         totl_npatches += level->getGlobalNumberOfPatches();
         totl_ncells += level->getBoxLevel()->getGlobalNumberOfCells();
      }
      os << border << "Total number of patches = " << totl_npatches << "\n";
      os << border << "Total number of cells = " << totl_ncells << "\n";
   }
   return 0;
}

PatchHierarchy::ConnectorWidthRequestorStrategy::~ConnectorWidthRequestorStrategy()
{
}

}
}

#endif
