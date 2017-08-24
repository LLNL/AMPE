/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   An AMR hierarchy of patch levels
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchHierarchy.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace hier {

const int PatchHierarchy::HIER_PATCH_HIERARCHY_VERSION = 3;

std::vector<const PatchHierarchy::ConnectorWidthRequestorStrategy *>
PatchHierarchy::s_class_cwrs;

tbox::StartupShutdownManager::Handler
PatchHierarchy::s_finalize_handler(
   0,
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
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(geometry->getDim()),
   d_number_levels(0),
   d_patch_descriptor(VariableDatabase::getDatabase()->getPatchDescriptor()),
   d_patch_factory(new PatchFactory),
   d_patch_level_factory(new PatchLevelFactory),
   d_max_levels(1),
   d_ratio_to_coarser(1, IntVector(IntVector::getOne(d_dim), geometry->getNumberBlocks())),
   d_proper_nesting_buffer(d_max_levels - 1, 1),
   d_smallest_patch_size(1, IntVector(d_dim, 1)),
   d_largest_patch_size(1, IntVector(d_dim, tbox::MathUtilities<int>::getMax())),
   d_allow_patches_smaller_than_ghostwidth(false),
   d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps(false),
   d_self_connector_widths(1, IntVector(IntVector::getOne(d_dim), geometry->getNumberBlocks())),
   d_fine_connector_widths(1, IntVector(IntVector::getOne(d_dim), geometry->getNumberBlocks())),
   d_connector_widths_committed(false),
   d_individual_cwrs()
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(geometry);

   d_object_name = object_name;
   d_grid_geometry = geometry;
   d_number_blocks = d_grid_geometry->getNumberBlocks();

   /*
    * Grab the physical domain (including periodic images) from the
    * grid geometry and set up domain data dependent on it.
    */
   d_domain_box_level.reset(new BoxLevel(
      d_ratio_to_coarser[0],
      getGridGeometry(),
      tbox::SAMRAI_MPI::getSAMRAIWorld(),
      BoxLevel::GLOBALIZED));
   d_grid_geometry->computePhysicalDomain(*d_domain_box_level,
      d_ratio_to_coarser[0]);
   d_domain_box_level->finalize();

   d_individual_cwrs = s_class_cwrs;

   /*
    * Without input database, the default is single-level with no
    * patch size constraints.
    */
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   d_grid_geometry->setUpRatios(d_ratio_to_coarser);

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
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
   tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
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
   const boost::shared_ptr<tbox::Database>& input_db,
   bool is_from_restart)
{
   if (input_db) {
      if (!is_from_restart) {

         /*
          * Read input for maximum number of levels.
          */

         d_max_levels =
            input_db->getIntegerWithDefault("max_levels", 1);
         if (!(d_max_levels >= 1)) {
            INPUT_RANGE_ERROR("max_levels");
         }

         if (d_max_levels != int(d_ratio_to_coarser.size())) {
            d_ratio_to_coarser.resize(d_max_levels, d_ratio_to_coarser.back());
            d_smallest_patch_size.resize(d_max_levels,
               d_smallest_patch_size.back());
            d_largest_patch_size.resize(d_max_levels,
               d_largest_patch_size.back());
         }

         std::vector<std::string> level_names(d_max_levels,
                                              std::string("level_"));
         for (int ln = 0; ln < d_max_levels; ++ln) {
            level_names[ln] += tbox::Utilities::intToString(ln);
         }

         // Read in ratio_to_coarser.
         d_ratio_to_coarser[0].setAll(IntVector::getOne(d_dim));
         if (d_number_blocks == 1) {
            if (input_db->isDatabase("ratio_to_coarser")) {
               const boost::shared_ptr<tbox::Database> tmp_db(
                  input_db->getDatabase("ratio_to_coarser"));
               for (int ln = 1; ln < d_max_levels; ++ln) {
                  if (tmp_db->isInteger(level_names[ln])) {
                     IntVector read_vector(d_dim);
                     tmp_db->getIntegerArray(level_names[ln],
                        &read_vector[0], 
                        //&d_ratio_to_coarser[ln].getBlockVector(BlockId(0))[0],
                        d_dim.getValue());
                     for (int i = 0; i < d_dim.getValue(); ++i) {
                        if (!(read_vector[i] > 0)) {
                           INPUT_RANGE_ERROR("ratio_to_coarser");
                        }
                     }
                     d_ratio_to_coarser[ln].setAll(read_vector);
                  } else {
                     d_ratio_to_coarser[ln] = d_ratio_to_coarser[ln - 1];
                  }
               }
            }
         } else {
            std::vector< std::vector<IntVector> > block_ratio(d_max_levels);
            for (int ln = 1; ln < d_max_levels; ++ln) {
               block_ratio[ln].resize(d_number_blocks, IntVector::getZero(d_dim));
            }
            for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
               std::string ratio_name("ratio_to_coarser_");
               ratio_name += tbox::Utilities::intToString(static_cast<int>(b));
               if (input_db->isDatabase(ratio_name)) {
                  const boost::shared_ptr<tbox::Database> tmp_db(
                     input_db->getDatabase(ratio_name));
                  for (int ln = 1; ln < d_max_levels; ++ln) {
                     if (tmp_db->isInteger(level_names[ln])) {
                        tmp_db->getIntegerArray(level_names[ln],
                           &block_ratio[ln][b][0],
                           d_dim.getValue());
                        for (int i = 0; i < d_dim.getValue(); ++i) {
                           if (!(block_ratio[ln][b][i] > 0)) {
                              INPUT_RANGE_ERROR(ratio_name);
                           }
                        }
                     } else {
                        block_ratio[ln][b] = block_ratio[ln-1][b];
                     }
                  }
               } else if (input_db->isDatabase("ratio_to_coarser")) {
                  const boost::shared_ptr<tbox::Database> tmp_db(
                     input_db->getDatabase("ratio_to_coarser"));
                  for (int ln = 1; ln < d_max_levels; ++ln) {
                     if (tmp_db->isInteger(level_names[ln])) {
                        tmp_db->getIntegerArray(level_names[ln],
                           &block_ratio[ln][b][0],
                           d_dim.getValue());
                        for (int i = 0; i < d_dim.getValue(); ++i) {
                           if (!(block_ratio[ln][b][i] > 0)) {
                              INPUT_RANGE_ERROR(ratio_name);
                           }
                        }
                     } else {
                        block_ratio[ln][b] = block_ratio[ln-1][b];
                     }
                  }
               }
            }
            for (int ln = 1; ln < d_max_levels; ++ln) {
               for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
                  for (int i = 0; i < d_dim.getValue(); ++i) {
                     d_ratio_to_coarser[ln](b,i) = block_ratio[ln][b][i];
                  }
               }
            }
         } 

         // Read in smallest_patch_size.
         if (input_db->isDatabase("smallest_patch_size")) {
            const boost::shared_ptr<tbox::Database> tmp_db(
               input_db->getDatabase("smallest_patch_size"));
            for (int ln = 0; ln < d_max_levels; ++ln) {
               if (tmp_db->isInteger(level_names[ln])) {
                  tmp_db->getIntegerArray(level_names[ln],
                     &d_smallest_patch_size[ln][0],
                     d_dim.getValue());
                  for (int i = 0; i < d_dim.getValue(); ++i) {
                     if (d_smallest_patch_size[ln][i] <= 0) {
                        INPUT_RANGE_ERROR("smallest_patch_size");
                     }
                  }
               } else {
                  d_smallest_patch_size[ln] = d_smallest_patch_size[ln - 1];
               }
            }
         }

         // Read in largest_patch_size.
         if (input_db->isDatabase("largest_patch_size")) {
            const boost::shared_ptr<tbox::Database> tmp_db(
               input_db->getDatabase("largest_patch_size"));
            for (int ln = 0; ln < d_max_levels; ++ln) {
               if (tmp_db->isInteger(level_names[ln])) {
                  tmp_db->getIntegerArray(level_names[ln],
                     &d_largest_patch_size[ln][0],
                     d_dim.getValue());
                  for (int i = 0; i < d_dim.getValue(); ++i) {
                     if (!(d_largest_patch_size[ln][i] < 0 ||
                           d_largest_patch_size[ln][i] >= d_smallest_patch_size[ln][i])) {
                        INPUT_RANGE_ERROR("largest_patch_size");
                     }
                     /*
                      * If largest patch size is input as negative, that means
                      * no largest size restriction is desired. We store
                      * an INT_MAX value in this case.
                      */
                     if (d_largest_patch_size[ln][i] < 0) {
                        d_largest_patch_size[ln][i] =
                           tbox::MathUtilities<int>::getMax();
                     }
                  }
               } else {
                  d_largest_patch_size[ln] = d_largest_patch_size[ln - 1];
               }
            }
         }

         std::vector<int> proper_nesting_buffer(1, 1);
         if (input_db->isInteger("proper_nesting_buffer")) {
            proper_nesting_buffer = input_db->getIntegerVector(
                  "proper_nesting_buffer");
         }
         d_proper_nesting_buffer.clear();
         for (int ln = 0; ln < d_max_levels - 1; ++ln) {
            if (ln < static_cast<int>(proper_nesting_buffer.size())) {
               d_proper_nesting_buffer.push_back(proper_nesting_buffer[ln]);
            } else {
               d_proper_nesting_buffer.push_back(d_proper_nesting_buffer[ln - 1]);
            }
         }
         for (size_t ln = 0; ln < d_proper_nesting_buffer.size(); ++ln) {
            if (d_proper_nesting_buffer[ln] < 0) {
               TBOX_ERROR(
                  d_object_name << ":  "
                                << "Key data `proper_nesting_buffer' has values < 0."
                                << std::endl);
            }
            if (d_proper_nesting_buffer[ln] == 0) {
               TBOX_WARNING(
                  d_object_name << ":  "
                                << "Using zero `proper_nesting_buffer' values."
                                << std::endl);
            }
         }

         d_allow_patches_smaller_than_ghostwidth =
            input_db->getBoolWithDefault("allow_patches_smaller_than_ghostwidth", false);

         d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
            input_db->getBoolWithDefault(
               "allow_patches_smaller_than_minimum_size_to_prevent_overlaps",
               false);
         if (d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
            TBOX_WARNING(
               d_object_name << ":  "
                             << "Allowing patches smaller than the given "
                             << "smallest patch size to prevent overlaps.\n"
                             << "Note:  If periodic "
                             << "boundary conditions are used, this flag is "
                             << "ignored in the periodic directions."
                             << std::endl);
         }
      } else {
         bool read_on_restart =
            input_db->getBoolWithDefault("read_on_restart", false);
         if (!read_on_restart) {
            return;
         }

         /*
          * Read input for maximum number of levels.
          */

         int new_max_levels =
            input_db->getIntegerWithDefault("max_levels", d_max_levels);
         if (new_max_levels > d_max_levels) {
            TBOX_ERROR("PatchHierarchy::getFromInput error...\n"
               << "max_levels must not increase on restart." << std::endl);
         }
         d_max_levels = new_max_levels;

         if (d_max_levels != int(d_ratio_to_coarser.size())) {
            d_ratio_to_coarser.resize(d_max_levels, d_ratio_to_coarser.back());
            d_smallest_patch_size.resize(d_max_levels,
               d_smallest_patch_size.back());
            d_largest_patch_size.resize(d_max_levels,
               d_largest_patch_size.back());
         }

         std::vector<std::string> level_names(d_max_levels,
                                              std::string("level_"));
         for (int ln = 0; ln < d_max_levels; ++ln) {
            level_names[ln] += tbox::Utilities::intToString(ln);
         }

         // Read in smallest_patch_size.
         if (input_db->isDatabase("smallest_patch_size")) {
            const boost::shared_ptr<tbox::Database> tmp_db(
               input_db->getDatabase("smallest_patch_size"));
            for (int ln = 0; ln < d_max_levels; ++ln) {
               if (tmp_db->isInteger(level_names[ln])) {
                  tmp_db->getIntegerArray(level_names[ln],
                     &d_smallest_patch_size[ln][0],
                     d_dim.getValue());
                  for (int i = 0; i < d_dim.getValue(); ++i) {
                     if (d_smallest_patch_size[ln][i] < 1) {
                        TBOX_ERROR("PatchHierarchy::getFromInput error...\n"
                           << "smallest_patch_size must be > 0." << std::endl);
                     }
                  }
               } else {
                  d_smallest_patch_size[ln] = d_smallest_patch_size[ln - 1];
               }
            }
         }

         // Read in largest_patch_size.
         if (input_db->isDatabase("largest_patch_size")) {
            const boost::shared_ptr<tbox::Database> tmp_db(
               input_db->getDatabase("largest_patch_size"));
            for (int ln = 0; ln < d_max_levels; ++ln) {
               if (tmp_db->isInteger(level_names[ln])) {
                  tmp_db->getIntegerArray(level_names[ln],
                     &d_largest_patch_size[ln][0],
                     d_dim.getValue());
                  for (int i = 0; i < d_dim.getValue(); ++i) {
                     if (d_largest_patch_size[ln][i] >= 0 &&
                         d_largest_patch_size[ln][i] <
                         d_smallest_patch_size[ln][i]) {
                        TBOX_ERROR("PatchHierarchy::getFromInput error...\n"
                           << "largest_patch_size must be >= smallest_patch_size."
                           << std::endl);
                     }
                     /*
                      * If largest patch size is input as negative, that means
                      * no largest size restriction is desired. We store
                      * an INT_MAX value in this case.
                      */
                     if (d_largest_patch_size[ln][i] < 0) {
                        d_largest_patch_size[ln][i] =
                           tbox::MathUtilities<int>::getMax();
                     }

                  }
               } else {
                  d_largest_patch_size[ln] = d_largest_patch_size[ln - 1];
               }
            }
         }

         if (input_db->keyExists("proper_nesting_buffer")) {
            std::vector<int> proper_nesting_buffer(1, 1);
            proper_nesting_buffer = input_db->getIntegerVector(
                  "proper_nesting_buffer");
            for (int ln = 0; ln < d_max_levels - 1; ++ln) {
               int val;
               if (ln < static_cast<int>(proper_nesting_buffer.size())) {
                  val = proper_nesting_buffer[ln];
               } else {
                  val = proper_nesting_buffer[static_cast<int>(proper_nesting_buffer.size()) - 1];
               }
               if (val != d_proper_nesting_buffer[ln]) {
                  TBOX_WARNING("PatchHierarchy::getFromInput warning...\n"
                     << "proper_nesting_buffer may not be changed on restart."
                     << std::endl);
                  break;
               }
            }
         }

         if (input_db->keyExists("allow_patches_smaller_than_ghostwidth")) {
            bool tmp = input_db->getBool("allow_patches_smaller_than_ghostwidth");
            if (tmp != d_allow_patches_smaller_than_ghostwidth) {
               TBOX_WARNING("PatchHierarchy::getFromInput warning...\n"
                  << "allow_patches_smaller_than_ghostwidth\n"
                  << "may not be changed on restart." << std::endl);
            }
         }

         d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
            input_db->getBoolWithDefault(
               "allow_patches_smaller_than_minimum_size_to_prevent_overlaps",
               d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps);
      }
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
   if (d_connector_widths_committed) {
      TBOX_ERROR("PatchHierarchy::registerConnectorWidthRequestor:\n"
         << "Registering a new ConnectorWidthRequestorStrategy is not\n"
         << "allowed for this hierarchy because some code has called\n"
         << "getRequiredConnectorWidth() with commit = true.  See the\n"
         << "documentation for getRequiredConnectorWidth()."
         << std::endl);
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
      s_class_cwrs[i] = 0;
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
   int head_ln,
   bool commit) const
{
   TBOX_ASSERT(head_ln >= 0);
   TBOX_ASSERT(head_ln < d_max_levels);
   TBOX_ASSERT(base_ln >= 0);
   TBOX_ASSERT(base_ln < d_max_levels);
   TBOX_ASSERT(abs(base_ln - head_ln) <= 1);

   if (!d_connector_widths_committed) {
      computeRequiredConnectorWidths();
      if (commit) {
         d_connector_widths_committed = true;
      }
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
         << "base_ln=" << base_ln << "  head_ln=" << head_ln << std::endl);
   }
   return d_self_connector_widths[base_ln];
}

/*
 *************************************************************************
 * Compute required Connector widths using all the registered
 * ConnectorWidthRequestorStrategy objects.
 *************************************************************************
 */
void
PatchHierarchy::computeRequiredConnectorWidths() const
{
   IntVector zero_vector(d_dim, 0, d_number_blocks);
   d_self_connector_widths.clear();
   d_self_connector_widths.resize(d_max_levels, zero_vector);
   for (int ln = 0; ln < d_max_levels; ++ln) {
      d_self_connector_widths[ln].setAll(zero_vector);
   }
   d_fine_connector_widths.clear();
   if (d_max_levels > 1) {
      d_fine_connector_widths.resize(d_max_levels - 1, zero_vector);
      for (int ln = 0; ln < d_max_levels-1; ++ln) {
         d_fine_connector_widths[ln].setAll(zero_vector);
      }
   }

   /*
    * Get the required widths satisfying all registered
    * ConnectorWidthRequestorStrategy objects.
    */

   
   std::vector<IntVector> self_connector_widths(0, d_self_connector_widths[0]);
   std::vector<IntVector> fine_connector_widths(0, d_fine_connector_widths[0]);
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
         for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
            d_fine_connector_widths[ln].max(fine_connector_widths[ln]);
         }
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
   const IntVector& refine_ratio) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, refine_ratio);
   TBOX_ASSERT(!fine_hierarchy_name.empty());
   TBOX_ASSERT(fine_hierarchy_name != d_object_name);
   TBOX_ASSERT(refine_ratio > IntVector::getZero(refine_ratio.getDim()));

   boost::shared_ptr<BaseGridGeometry> fine_geometry(
      d_grid_geometry->makeRefinedGridGeometry(
         fine_hierarchy_name + "GridGeometry",
         refine_ratio));

   PatchHierarchy* fine_hierarchy =
      new PatchHierarchy(fine_hierarchy_name,
         fine_geometry,
         boost::shared_ptr<tbox::Database>());

   // Set hierarchy parameters.

   fine_hierarchy->d_number_blocks = d_number_blocks;
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
   fine_hierarchy->d_grid_geometry->setUpRatios(d_ratio_to_coarser);

   for (int ln = 0; ln < d_number_levels; ++ln) {
      BoxContainer refined_boxes(d_patch_levels[ln]->getBoxLevel()->getBoxes());
      refined_boxes.refine(refine_ratio);
      boost::shared_ptr<BoxLevel> refined_box_level(
         boost::make_shared<BoxLevel>(
            d_patch_levels[ln]->getBoxLevel()->getRefinementRatio(),
            fine_geometry,
            d_patch_levels[ln]->getBoxLevel()->getMPI()));
      refined_box_level->swapInitialize(
         refined_boxes,
         d_patch_levels[ln]->getBoxLevel()->getRefinementRatio(),
         fine_geometry,
         d_patch_levels[ln]->getBoxLevel()->getMPI());
      fine_hierarchy->makeNewPatchLevel(ln, refined_box_level);
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
   const IntVector& coarsen_ratio) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, coarsen_ratio);
   TBOX_ASSERT(!coarse_hierarchy_name.empty());
   TBOX_ASSERT(coarse_hierarchy_name != d_object_name);
   TBOX_ASSERT(coarsen_ratio > IntVector::getZero(coarsen_ratio.getDim()));

   boost::shared_ptr<BaseGridGeometry> coarse_geometry(
      d_grid_geometry->makeCoarsenedGridGeometry(
         coarse_hierarchy_name + "GridGeometry",
         coarsen_ratio));

   PatchHierarchy* coarse_hierarchy =
      new PatchHierarchy(coarse_hierarchy_name,
         coarse_geometry,
         boost::shared_ptr<tbox::Database>());

   // Set hierarchy parameters.

   coarse_hierarchy->d_number_blocks = d_number_blocks;
   coarse_hierarchy->d_max_levels = d_max_levels;
   coarse_hierarchy->d_ratio_to_coarser = d_ratio_to_coarser;
   coarse_hierarchy->d_smallest_patch_size = d_smallest_patch_size;
   coarse_hierarchy->d_largest_patch_size = d_largest_patch_size;
   coarse_hierarchy->d_individual_cwrs = d_individual_cwrs;
   coarse_hierarchy->d_proper_nesting_buffer = d_proper_nesting_buffer;
   coarse_hierarchy->d_grid_geometry->setUpRatios(d_ratio_to_coarser);

   for (int ln = 0; ln < d_number_levels; ++ln) {
      BoxContainer coarsened_boxes(d_patch_levels[ln]->getBoxLevel()->getBoxes());
      coarsened_boxes.coarsen(coarsen_ratio);
      boost::shared_ptr<BoxLevel> coarsened_box_level(
         boost::make_shared<BoxLevel>(
            d_patch_levels[ln]->getBoxLevel()->getRefinementRatio(),
            coarse_geometry,
            d_patch_levels[ln]->getBoxLevel()->getMPI()));
      coarsened_box_level->swapInitialize(
         coarsened_boxes,
         d_patch_levels[ln]->getBoxLevel()->getRefinementRatio(),
         coarse_geometry,
         d_patch_levels[ln]->getBoxLevel()->getMPI());
      coarse_hierarchy->makeNewPatchLevel(ln, coarsened_box_level);
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
   const BoxLevel& new_box_level)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, new_box_level);
   TBOX_ASSERT(ln >= 0);
   TBOX_ASSERT(new_box_level.getRefinementRatio() > IntVector::getZero(d_dim));
   TBOX_ASSERT(new_box_level.getGridGeometry() == d_grid_geometry);
   TBOX_ASSERT(d_domain_box_level->getGridGeometry() == d_grid_geometry);

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
      if (new_box_level.getRefinementRatio() != expected_ratio) {
         TBOX_ERROR("PatchHierarchy::makeNewPatchLevel: patch level "
            << ln << " has refinement ratio "
            << new_box_level.getRefinementRatio()
            << ", it should be " << expected_ratio << std::endl);
      }
   }
   if (static_cast<int>(d_patch_levels.size()) > ln &&
       d_patch_levels[ln].get() != 0) {
      TBOX_ERROR("PatchHierarchy::makeNewPatchLevel: patch level "
         << ln << " already exists. "
         << "Remove old level from the hierarchy before making "
         << "a new level in its place." << std::endl);
   }

   if (ln >= d_number_levels) {
      d_number_levels = ln + 1;
      d_patch_levels.resize(d_number_levels);
   }

   d_patch_levels[ln] = d_patch_level_factory->allocate(
         new_box_level,
         d_grid_geometry,
         d_patch_descriptor,
         d_patch_factory);
   d_patch_levels[ln]->getBoxLevel()->cacheGlobalReducedData();

   d_patch_levels[ln]->setLevelNumber(ln);
   d_patch_levels[ln]->setNextCoarserHierarchyLevelNumber(ln - 1);
   d_patch_levels[ln]->setLevelInHierarchy(true);

   if ((ln > 0) && d_patch_levels[ln - 1]) {
      IntVector ratio = d_patch_levels[ln]->getRatioToLevelZero() /
         d_patch_levels[ln-1]->getRatioToLevelZero();
      d_patch_levels[ln]->setRatioToCoarserLevel(ratio);
   }

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
   const boost::shared_ptr<BoxLevel> new_box_level)
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(d_dim, *new_box_level);
   TBOX_ASSERT(ln >= 0);
   TBOX_ASSERT(new_box_level->getRefinementRatio() > IntVector::getZero(d_dim));
   TBOX_ASSERT(new_box_level->getGridGeometry() == d_grid_geometry);
   TBOX_ASSERT(d_domain_box_level->getGridGeometry() == d_grid_geometry);

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
      if (new_box_level->getRefinementRatio() != expected_ratio) {
         TBOX_ERROR("PatchHierarchy::makeNewPatchLevel: patch level "
            << ln << " has refinement ratio "
            << new_box_level->getRefinementRatio()
            << ", it should be " << expected_ratio << std::endl);
      }
   }
   if (static_cast<int>(d_patch_levels.size()) > ln &&
       d_patch_levels[ln].get() != 0) {
      TBOX_ERROR("PatchHierarchy::makeNewPatchLevel: patch level "
         << ln << " already exists. "
         << "Remove old level from the hierarchy before making "
         << "a new level in its place." << std::endl);
   }

   if (ln >= d_number_levels) {
      d_number_levels = ln + 1;
      d_patch_levels.resize(d_number_levels);
   }

   d_patch_levels[ln] = d_patch_level_factory->allocate(
         new_box_level,
         d_grid_geometry,
         d_patch_descriptor,
         d_patch_factory);
   d_patch_levels[ln]->getBoxLevel()->cacheGlobalReducedData();

   d_patch_levels[ln]->setLevelNumber(ln);
   d_patch_levels[ln]->setNextCoarserHierarchyLevelNumber(ln - 1);
   d_patch_levels[ln]->setLevelInHierarchy(true);

   if ((ln > 0) && d_patch_levels[ln - 1]) {
      IntVector ratio = d_patch_levels[ln]->getRatioToLevelZero() /
         d_patch_levels[ln-1]->getRatioToLevelZero();
      d_patch_levels[ln]->setRatioToCoarserLevel(ratio);
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
      --d_number_levels;
   }
}

/*
 *************************************************************************
 * Log the given level, its peer connector and if requested, the
 * connectors to the next finer and next coarser levels.  Connectors
 * logged will have width required by the hierarchy.
 *************************************************************************
 */
void
PatchHierarchy::logMetadataStatistics(
   const char* note,
   int ln,
   int cycle,
   double level_time,
   bool log_fine_connector,
   bool log_coarse_connector) const
{
   const std::string name("L" + tbox::Utilities::levelToString(ln));
   const boost::shared_ptr<PatchLevel> level =
      getPatchLevel(ln);
   const BoxLevel& box_level = *level->getBoxLevel();

   tbox::plog << "PatchHierarchy metadata statistics '"
              << note << "', at cycle " << cycle
              << ", time " << level_time << ", added "
              << name << ":\n"
              << box_level.format("\t", 0)
              << '\t' << name << " statistics:\n"
              << box_level.formatStatistics("\t\t");

   const Connector& peer_conn =
      level->findConnector(*level,
         getRequiredConnectorWidth(ln, ln),
         CONNECTOR_CREATE,
         true);
   tbox::plog << "\tL" << ln
              << " Peer connector:\n" << peer_conn.format("\t\t", 0)
              << "\tL"
              << ln << " peer Connector statistics:\n" << peer_conn.formatStatistics("\t\t");

   if (log_fine_connector) {
      const Connector& to_fine =
         level->findConnector(*getPatchLevel(ln + 1),
            getRequiredConnectorWidth(ln, ln + 1),
            CONNECTOR_CREATE,
            true);
      tbox::plog << "\tL" << ln << "->L" << ln + 1
                 << " Connector:\n" << to_fine.format("\t\t", 0)
                 << "\tL" << ln << "->L" << ln + 1
                 << " Connector statistics:\n" << to_fine.formatStatistics("\t\t");
      const Connector& from_fine =
         getPatchLevel(ln + 1)->findConnector(*level,
            getRequiredConnectorWidth(ln + 1, ln),
            CONNECTOR_CREATE,
            true);
      tbox::plog << "\tL" << ln + 1 << "->L" << ln
                 << " Connector:\n" << from_fine.format("\t\t", 0)
                 << "\tL" << ln + 1 << "->L" << ln
                 << " Connector statistics:\n" << from_fine.formatStatistics("\t\t");
   }

   if (log_coarse_connector) {
      const Connector& to_crse =
         level->findConnector(*getPatchLevel(ln - 1),
            getRequiredConnectorWidth(ln, ln - 1),
            CONNECTOR_CREATE,
            true);
      tbox::plog << "\tL" << ln << "->L" << ln - 1
                 << " Connector:\n" << to_crse.format("\t\t", 0)
                 << "\tL" << ln << "->L" << ln - 1
                 << " Connector statistics:\n" << to_crse.formatStatistics("\t\t");
      const Connector& from_crse =
         getPatchLevel(ln - 1)->findConnector(*level,
            getRequiredConnectorWidth(ln - 1, ln),
            CONNECTOR_CREATE,
            true);
      tbox::plog << "\tL" << ln - 1 << "->L" << ln
                 << " Connector:\n" << from_crse.format("\t\t", 0)
                 << "\tL" << ln - 1 << "->L" << ln
                 << " Connector statistics:\n" << from_crse.formatStatistics("\t\t");
   }
}

/*
 *************************************************************************
 *
 * Writes the class version number and the number of levels in the
 * hierarchy to the restart database.  Each patch_level write itself out
 * to the restart database.  The database keys for the patch levels are
 * given by "level#" where # is the level number for the patch_level.
 * The patchdata that are written to the database are determined by
 * which those bits in the VariableDatabase restart table.
 *
 * Asserts that the restart_db pointer passed in is not NULL.
 *
 *************************************************************************
 */

void
PatchHierarchy::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("HIER_PATCH_HIERARCHY_VERSION",
      HIER_PATCH_HIERARCHY_VERSION);

   restart_db->putInteger("d_number_levels", d_number_levels);

   std::vector<std::string> level_names(d_max_levels);
   const std::string prefix("level_");
   for (int ln = 0; ln < d_max_levels; ++ln) {
      level_names[ln] = prefix + tbox::Utilities::levelToString(ln);
   }

   /*
    * Write hierarchy parameters.
    */
   restart_db->putInteger("max_levels", d_max_levels);

   boost::shared_ptr<tbox::Database> ratio_to_coarser_db(
      restart_db->putDatabase("ratio_to_coarser"));
   for (int ln = 0; ln < d_max_levels; ++ln) {
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         std::vector<int> ratio_vec(d_dim.getValue());
         for (int d = 0; d < d_dim.getValue(); ++d) {
            ratio_vec[d] = (d_ratio_to_coarser[ln](b,d));
         }
         std::string level_block_name(level_names[ln]);
         level_block_name += "_";
         level_block_name += tbox::Utilities::intToString(static_cast<int>(b));
         ratio_to_coarser_db->putIntegerArray(level_block_name,
            &ratio_vec[0],
            d_dim.getValue());
      }
   }

   boost::shared_ptr<tbox::Database> smallest_patch_db(
      restart_db->putDatabase("smallest_patch_size"));
   for (int ln = 0; ln < d_max_levels; ++ln) {
      smallest_patch_db->putIntegerArray(level_names[ln],
         &d_smallest_patch_size[ln][0],
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> largest_patch_db(
      restart_db->putDatabase("largest_patch_size"));
   for (int ln = 0; ln < d_max_levels; ++ln) {
      largest_patch_db->putIntegerArray(level_names[ln],
         &d_largest_patch_size[ln][0],
         d_dim.getValue());
   }

   if (d_max_levels > 1) {
      restart_db->putIntegerVector("proper_nesting_buffer",
         d_proper_nesting_buffer);
   }

   restart_db->putBool("allow_patches_smaller_than_ghostwidth",
      d_allow_patches_smaller_than_ghostwidth);

   restart_db->putBool("allow_patches_smaller_than_minimum_size_to_prevent_overlaps",
      d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps);

   boost::shared_ptr<tbox::Database> self_connector_widths_db(
      restart_db->putDatabase("d_self_connector_widths"));
   d_self_connector_widths.resize(d_max_levels, IntVector(d_dim, 1, d_number_blocks));
   for (int ln = 0; ln < d_max_levels; ++ln) {

      std::vector<int> put_self_widths(d_number_blocks * d_dim.getValue());
      int ic = 0;
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         for (int d = 0; d < d_dim.getValue(); ++d) {
            put_self_widths[ic] = d_self_connector_widths[ln](b,d);
            ++ic;
         }
      }
      self_connector_widths_db->putIntegerVector(level_names[ln],
         put_self_widths);
   }

   boost::shared_ptr<tbox::Database> fine_connector_widths_db(
      restart_db->putDatabase("d_fine_connector_widths"));
   d_fine_connector_widths.resize(d_max_levels-1, IntVector(d_dim, 1, d_number_blocks));
   for (int ln = 0; ln < d_max_levels - 1; ++ln) {

      std::vector<int> put_fine_widths(d_number_blocks * d_dim.getValue());
      int ic = 0;
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         for (int d = 0; d < d_dim.getValue(); ++d) {
            put_fine_widths[ic] = d_self_connector_widths[ln](b,d);
            ++ic;
         }
      }
      fine_connector_widths_db->putIntegerVector(level_names[ln],
         put_fine_widths);
   }

   for (int i = 0; i < d_number_levels; ++i) {

      boost::shared_ptr<tbox::Database> level_database(
         restart_db->putDatabase(level_names[i]));

      d_patch_levels[i]->putToRestart(level_database);
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

   /*
    * Read hierarchy paremeters.
    */

   int ver = database->getInteger("HIER_PATCH_HIERARCHY_VERSION");
   if (ver != HIER_PATCH_HIERARCHY_VERSION) {
      TBOX_ERROR("PatchHierarchy::getFromRestart error...\n"
         << "  object name = " << d_object_name
         << " : Restart file version different than class version" << std::endl);
   }

   d_number_levels = database->getInteger("d_number_levels");
   if (d_number_levels <= 0) {
      TBOX_ERROR("PatchHierarchy::getFromRestart error ...\n"
         << "  object name = " << d_object_name
         << " : `d_number_levels' is <= zero in restart file" << std::endl);
   }

   d_max_levels = database->getInteger("max_levels");

   std::vector<std::string> level_names(d_max_levels);
   const std::string prefix("level_");
   for (int ln = 0; ln < d_max_levels; ++ln) {
      level_names[ln] = prefix + tbox::Utilities::levelToString(ln);
   }

   boost::shared_ptr<tbox::Database> ratio_to_coarser_db(
      database->getDatabase("ratio_to_coarser"));
   d_ratio_to_coarser.resize(d_max_levels, d_ratio_to_coarser.back());
   for (int ln = 0; ln < d_max_levels; ++ln) {
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         std::string level_block_name(level_names[ln]);
         level_block_name += "_";
         level_block_name += tbox::Utilities::intToString(static_cast<int>(b));
         std::vector<int> ratio_vec(d_dim.getValue());
         ratio_to_coarser_db->getIntegerArray(level_block_name,
            &ratio_vec[0],
            d_dim.getValue());
         for (int d = 0; d < d_dim.getValue(); ++d) {
            d_ratio_to_coarser[ln](b,d) = ratio_vec[d];
         }
      }
   }

   boost::shared_ptr<tbox::Database> smallest_patch_db(
      database->getDatabase("smallest_patch_size"));
   d_smallest_patch_size.resize(d_max_levels, d_smallest_patch_size.back());
   for (int ln = 0; ln < d_max_levels; ++ln) {
      smallest_patch_db->getIntegerArray(level_names[ln],
         &d_smallest_patch_size[ln][0],
         d_dim.getValue());
   }

   boost::shared_ptr<tbox::Database> largest_patch_db(
      database->getDatabase("largest_patch_size"));
   d_largest_patch_size.resize(d_max_levels, d_largest_patch_size.back());
   for (int ln = 0; ln < d_max_levels; ++ln) {
      largest_patch_db->getIntegerArray(level_names[ln],
         &d_largest_patch_size[ln][0],
         d_dim.getValue());
   }

   d_proper_nesting_buffer.resize(d_max_levels - 1, 0);
   if (d_max_levels > 1) {
      d_proper_nesting_buffer =
         database->getIntegerVector("proper_nesting_buffer");
   }

   d_allow_patches_smaller_than_ghostwidth = database->getBool(
         "allow_patches_smaller_than_ghostwidth");

   d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps =
      database->getBool(
         "allow_patches_smaller_than_minimum_size_to_prevent_overlaps");
   if (d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps) {
      TBOX_WARNING(
         d_object_name << ":  "
                       << "Allowing patches smaller than the given "
                       << "smallest patch size.  Note:  If periodic "
                       << "boundary conditions are used, this flag is "
                       << "ignored in the periodic directions." << std::endl);
   }

   boost::shared_ptr<tbox::Database> self_connector_widths_db(
      database->getDatabase("d_self_connector_widths"));
   d_self_connector_widths.resize(d_max_levels, IntVector(d_dim, 1, d_number_blocks));
   for (int ln = 0; ln < d_max_levels; ++ln) {
      std::vector<int> get_self_widths =
         self_connector_widths_db->getIntegerVector(level_names[ln]);

      d_self_connector_widths[ln] = IntVector(d_dim, 0, d_number_blocks);
      int ic = 0;
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         for (int d = 0; d < d_dim.getValue(); ++d) {
            d_self_connector_widths[ln](b,d) = get_self_widths[ic];
            ++ic;
         }
      }
   }

   boost::shared_ptr<tbox::Database> fine_connector_widths_db(
      database->getDatabase("d_fine_connector_widths"));
   d_fine_connector_widths.resize(d_max_levels - 1, IntVector(d_dim, 1, d_number_blocks));
   for (int ln = 0; ln < d_max_levels - 1; ++ln) {
      std::vector<int> get_fine_widths =
         fine_connector_widths_db->getIntegerVector(level_names[ln]);

      d_fine_connector_widths[ln] = IntVector(d_dim, 0, d_number_blocks);
      int ic = 0;
      for (BlockId::block_t b = 0; b < d_number_blocks; ++b) {
         for (int d = 0; d < d_dim.getValue(); ++d) {
            d_fine_connector_widths[ln](b,d) = get_fine_widths[ic];
            ++ic;
         }
      }
   }
}

void
PatchHierarchy::initializeHierarchy()
{
   boost::shared_ptr<tbox::Database> restart_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!restart_db->isDatabase(d_object_name)) {
      TBOX_ERROR("PatchHierarchy::initializeHierarchy() error...\n"
         << "   Restart database with name "
         << d_object_name << " not found in restart file" << std::endl);
   }
   boost::shared_ptr<tbox::Database> database(
      restart_db->getDatabase(d_object_name));

   d_patch_levels.resize(d_number_levels);
   for (int i = 0; i < d_number_levels; ++i) {
      std::string level_name = "level_" + tbox::Utilities::levelToString(i);

      boost::shared_ptr<tbox::Database> level_database(
         database->getDatabase(level_name));

      d_patch_levels[i] = d_patch_level_factory->allocate(
            level_database,
            d_grid_geometry,
            d_patch_descriptor,
            d_patch_factory,
            false);
   }
   /*
    * Compute Connectors.
    * BTNG TODO: This should be replaced by writing edges to
    * restart and reading them back.
    */
   for (int i = 0; i < d_number_levels; ++i) {
      d_patch_levels[i]->findConnector(*d_patch_levels[i],
         getRequiredConnectorWidth(i, i),
         CONNECTOR_CREATE);
      if (i < d_number_levels - 1) {
         d_patch_levels[i]->findConnectorWithTranspose(*d_patch_levels[i + 1],
            getRequiredConnectorWidth(i, i + 1),
            getRequiredConnectorWidth(i + 1, i),
            CONNECTOR_CREATE);
      }
   }

}

int
PatchHierarchy::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int depth)
{
   size_t totl_npatches = 0;
   size_t totl_ncells = 0;
   int nlevels = getNumberOfLevels();
   os << border << "Domain of hierarchy:\n" << d_domain_box_level->format(border, 2) << '\n'
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
