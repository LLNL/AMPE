/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A collection of patches at one level of the AMR hierarchy
 *
 ************************************************************************/

#ifndef included_hier_PatchLevel_C
#define included_hier_PatchLevel_C

#include "SAMRAI/hier/PatchLevel.h"

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

#include <boost/make_shared.hpp>
#include <cstdio>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int PatchLevel::HIER_PATCH_LEVEL_VERSION = 3;

static boost::shared_ptr<tbox::Timer> t_level_constructor;
static boost::shared_ptr<tbox::Timer> t_constructor_setup;
static boost::shared_ptr<tbox::Timer> t_constructor_phys_domain;
static boost::shared_ptr<tbox::Timer> t_constructor_touch_boundaries;
static boost::shared_ptr<tbox::Timer> t_constructor_set_geometry;
static boost::shared_ptr<tbox::Timer> t_set_patch_touches;
static boost::shared_ptr<tbox::Timer> t_constructor_compute_shifts;

tbox::StartupShutdownManager::Handler
PatchLevel::s_initialize_finalize_handler(
   PatchLevel::initializeCallback,
   0,
   0,
   PatchLevel::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 *************************************************************************
 *
 * Default patch level constructor sets default (non-usable) state.
 *
 *************************************************************************
 */

PatchLevel::PatchLevel(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_has_globalized_data(false),
   d_ratio_to_level_zero(IntVector::getZero(dim)),
   d_factory(boost::make_shared<PatchFactory>()),
   d_local_number_patches(0),
   d_physical_domain(0),
   d_ratio_to_coarser_level(IntVector::getZero(dim)),
   d_level_number(-1),
   d_next_coarser_level_number(-1),
   d_in_hierarchy(false)
{
   t_level_constructor->start();

   t_level_constructor->stop();
}

/*
 *************************************************************************
 *
 * Create a new patch level using the specified boxes and processor
 * mapping.  Only those patches that are local to the processor are
 * allocated.  Allocate patches using the specified patch factory or
 * the standard patch factory if none is explicitly specified.
 *
 *************************************************************************
 */

PatchLevel::PatchLevel(
   const BoxLevel& mapped_box_level,
   const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
   const boost::shared_ptr<PatchDescriptor>& descriptor,
   const boost::shared_ptr<PatchFactory>& factory,
   bool defer_boundary_box_creation):
   d_dim(grid_geometry->getDim()),
   d_mapped_box_level(boost::make_shared<BoxLevel>(mapped_box_level)),
   d_has_globalized_data(false),
   d_ratio_to_level_zero(d_mapped_box_level->getRefinementRatio()),
   d_factory(factory ? factory : boost::make_shared<PatchFactory>()),
   d_physical_domain(grid_geometry->getNumberBlocks()),
   d_ratio_to_coarser_level(grid_geometry->getDim(), 0)

{
   d_number_blocks = grid_geometry->getNumberBlocks();

   TBOX_DIM_ASSERT_CHECK_ARGS2(mapped_box_level, *grid_geometry);

   t_level_constructor->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(grid_geometry);
   TBOX_ASSERT(descriptor);
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(mapped_box_level.getRefinementRatio() !=
      IntVector::getZero(getDim()));

   if (getDim().getValue() > 1) {
      for (int i = 0; i < getDim().getValue(); i++) {
         TBOX_ASSERT((mapped_box_level.getRefinementRatio() (i)
                      * mapped_box_level.getRefinementRatio() ((i
                                                                + 1)
                         % getDim().getValue()) > 0)
            || (mapped_box_level.getRefinementRatio() (i) == 1)
            || (mapped_box_level.getRefinementRatio() ((i + 1) % getDim().getValue()) ==
                1));
      }
   }
#endif

   t_constructor_setup->start();

   d_local_number_patches =
      static_cast<int>(d_mapped_box_level->getLocalNumberOfBoxes());
   d_descriptor = descriptor;

   d_geometry = grid_geometry;

   d_level_number = -1;
   d_next_coarser_level_number = -1;
   d_in_hierarchy = false;

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   for (RealBoxConstIterator ni(mapped_boxes.realBegin());
        ni != mapped_boxes.realEnd(); ++ni) {
      const Box& mapped_box = *ni;
      const BoxId& ip = mapped_box.getId();
      boost::shared_ptr<Patch>& patch(d_patches[ip]);
      patch = d_factory->allocate(mapped_box, d_descriptor);
      patch->setPatchLevelNumber(d_level_number);
      patch->setPatchInHierarchy(d_in_hierarchy);
   }

   d_boundary_boxes_created = false;
   t_constructor_setup->stop();

   t_constructor_phys_domain->start();
   for (int nb = 0; nb < d_number_blocks; nb++) {
      grid_geometry->computePhysicalDomain(d_physical_domain[nb],
         d_ratio_to_level_zero, BlockId(nb));
   }
   t_constructor_phys_domain->stop();

   t_constructor_touch_boundaries->start();
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;
   grid_geometry->findPatchesTouchingBoundaries(
      touches_regular_bdry,
      touches_periodic_bdry,
      *this);
   t_constructor_touch_boundaries->stop();

   t_constructor_set_geometry->start();
   grid_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);
   t_constructor_set_geometry->stop();

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   t_level_constructor->stop();
}

/*
 *************************************************************************
 *
 * Create a new patch level from information in the given database.
 *
 *************************************************************************
 */

PatchLevel::PatchLevel(
   const boost::shared_ptr<tbox::Database>& level_database,
   const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
   const boost::shared_ptr<PatchDescriptor>& descriptor,
   const boost::shared_ptr<PatchFactory>& factory,
   const ComponentSelector& component_selector,
   bool defer_boundary_box_creation):
   d_dim(grid_geometry->getDim()),
   d_has_globalized_data(false),
   d_ratio_to_level_zero(IntVector(grid_geometry->getDim(),
                                   tbox::MathUtilities<int>::getMax())),
   d_factory(factory ? factory : boost::make_shared<PatchFactory>()),
   d_physical_domain(grid_geometry->getNumberBlocks()),
   d_ratio_to_coarser_level(IntVector(grid_geometry->getDim(),
                                      tbox::MathUtilities<int>::getMax()))
{
   d_number_blocks = grid_geometry->getNumberBlocks();

   TBOX_ASSERT(level_database);
   TBOX_ASSERT(grid_geometry);
   TBOX_ASSERT(descriptor);

   t_level_constructor->start();

   d_geometry = grid_geometry;
   d_descriptor = descriptor;

   getFromDatabase(level_database, component_selector);

   d_boundary_boxes_created = false;

   t_constructor_touch_boundaries->start();
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;
   grid_geometry->findPatchesTouchingBoundaries(
      touches_regular_bdry,
      touches_periodic_bdry,
      *this);
   t_constructor_touch_boundaries->stop();

   t_constructor_set_geometry->start();
   grid_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);
   t_constructor_set_geometry->stop();

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

   t_level_constructor->stop();
}

PatchLevel::~PatchLevel()
{
}

/*
 * ************************************************************************
 *
 * Set data members of this patch level by refining information on
 * the argument level by the given ratio.
 *
 * ************************************************************************
 */

void
PatchLevel::setRefinedPatchLevel(
   const boost::shared_ptr<PatchLevel>& coarse_level,
   const IntVector& refine_ratio,
   const boost::shared_ptr<BaseGridGeometry>& fine_grid_geometry,
   bool defer_boundary_box_creation)
{
   TBOX_ASSERT(coarse_level);
   TBOX_ASSERT(refine_ratio > IntVector::getZero(getDim()));
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, *coarse_level, refine_ratio);
   if (fine_grid_geometry) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *fine_grid_geometry);
   }
#endif

   /*
    * The basic state of the new patch level is initialized from the state of
    * the given existing patch level.
    */

   // d_global_number_patches = coarse_level->d_global_number_patches;
   d_descriptor = coarse_level->d_descriptor;
   d_factory = coarse_level->d_factory;
   // d_mapping.setProcessorMapping( (coarse_level->d_mapping).getProcessorMapping() );

   /*
    * Compute the ratio to coarsest level (reference level in hierarchy --
    * usually level  zero) and set grid geometry for this (fine) level.  If
    * pointer to given fine grid geometry is null, then it is assumed that
    * this level is to use the same grid geometry as the given coarse level
    * and the ratio to level zero is set relative to the give coarse level.
    * Otherwise, use given grid geometry and copy ratio to level zero from
    * given coarse level.
    */

   if (!fine_grid_geometry) {

      d_geometry = coarse_level->d_geometry;

      const IntVector& coarse_ratio = coarse_level->getRatioToLevelZero();
      for (int i = 0; i < getDim().getValue(); i++) {
         int coarse_rat = coarse_ratio(i);
         int refine_rat = refine_ratio(i);
         if (coarse_rat < 0) {
            if (tbox::MathUtilities<int>::Abs(coarse_rat) >= refine_rat) {
               d_ratio_to_level_zero(i) =
                  -(tbox::MathUtilities<int>::Abs(coarse_rat / refine_rat));
            } else {
               d_ratio_to_level_zero(i) =
                  tbox::MathUtilities<int>::Abs(refine_rat / coarse_rat);
            }
         } else {
            d_ratio_to_level_zero(i) = coarse_rat * refine_rat;
         }

      }

   } else {

      d_geometry = fine_grid_geometry;

      d_ratio_to_level_zero = coarse_level->d_ratio_to_level_zero;
   }

   /*
    * Set global box array and index space for level based on refining
    * coarse level information.
    */

   d_boxes = coarse_level->d_boxes;
   d_boxes.refine(refine_ratio);

   {
      d_mapped_box_level.reset(
         new BoxLevel(
            d_ratio_to_level_zero,
            d_geometry,
            coarse_level->d_mapped_box_level->getMPI()));
      coarse_level->d_mapped_box_level->refineBoxes(
         *d_mapped_box_level,
         refine_ratio,
         d_ratio_to_level_zero);
      d_mapped_box_level->finalize();
   }
   d_local_number_patches = coarse_level->getLocalNumberOfPatches();
   d_number_blocks = coarse_level->d_number_blocks;

   d_physical_domain.resizeArray(d_number_blocks);
   for (int nb = 0; nb < d_number_blocks; nb++) {
      d_physical_domain[nb] = coarse_level->d_physical_domain[nb];
      d_physical_domain[nb].refine(refine_ratio);
   }

   /*
    * Allocate arrays of patches and patch information.  Then, allocate and
    * initialize patch objects.  Finally, set patch geometry and remaining
    * domain information.
    */

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   for (RealBoxConstIterator ni(mapped_boxes.realBegin());
        ni != mapped_boxes.realEnd(); ++ni) {
      const Box& mapped_box = *ni;
      const BoxId& mapped_box_id = mapped_box.getId();
      d_patches[mapped_box_id] = d_factory->allocate(mapped_box, d_descriptor);
      d_patches[mapped_box_id]->setPatchLevelNumber(d_level_number);
      d_patches[mapped_box_id]->setPatchInHierarchy(d_in_hierarchy);
   }

   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;

   for (iterator ip(coarse_level->begin()); ip != coarse_level->end(); ip++) {
      boost::shared_ptr<PatchGeometry> coarse_pgeom((*ip)->getPatchGeometry());

      /* If map does not contain values create them */
      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_regular_bdry(
         touches_regular_bdry.find(ip->getBox().getId()));
      if (iter_touches_regular_bdry == touches_regular_bdry.end()) {
         iter_touches_regular_bdry = touches_regular_bdry.insert(
               iter_touches_regular_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_periodic_bdry(
         touches_periodic_bdry.find(ip->getBox().getId()));
      if (iter_touches_periodic_bdry == touches_periodic_bdry.end()) {
         iter_touches_periodic_bdry = touches_periodic_bdry.insert(
               iter_touches_periodic_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      PatchGeometry::TwoDimBool&
      touches_regular_bdry_ip((*iter_touches_regular_bdry).second);

      for (int axis = 0; axis < getDim().getValue(); axis++) {
         for (int side = 0; side < 2; side++) {
            touches_regular_bdry_ip(axis, side) =
               coarse_pgeom->getTouchesRegularBoundary(axis, side);
         }
      }
   }

   d_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

}

/*
 * ************************************************************************
 *
 * Set data members of this patch level by coarsening information on
 * the argument level by the given ratio.
 *
 * ************************************************************************
 */

void
PatchLevel::setCoarsenedPatchLevel(
   const boost::shared_ptr<PatchLevel>& fine_level,
   const IntVector& coarsen_ratio,
   const boost::shared_ptr<BaseGridGeometry>& coarse_grid_geom,
   bool defer_boundary_box_creation)
{
   TBOX_ASSERT(fine_level);
   TBOX_ASSERT(coarsen_ratio > IntVector::getZero(getDim()));

#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, *fine_level, coarsen_ratio);
   if (coarse_grid_geom) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *coarse_grid_geom);
   }
#endif

   /*
    * The basic state of the new patch level is initialized from the state of
    * the given existing patch level.
    */

   // d_global_number_patches = fine_level->d_global_number_patches;
   d_descriptor = fine_level->d_descriptor;
   d_factory = fine_level->d_factory;
   // d_mapping.setProcessorMapping( (fine_level->d_mapping).getProcessorMapping() );

   /*
    * Compute the ratio to coarsest level (reference level in hierarchy --
    * usually level zero) and set grid geometry for this (coarse) level.  If
    * pointer to a given coarse grid geometry is null, then it is assumed
    * that this level is to use the same grid geometry as the given fine
    * level and the ratio to level zero is set relative to the given fine
    * level.  Otherwise, use given grid geometry and copy ratio to level zero
    * from given fine level.
    */

   if (!coarse_grid_geom) {

      d_geometry = fine_level->d_geometry;

      const IntVector& fine_ratio =
         fine_level->d_ratio_to_level_zero;

      for (int i = 0; i < getDim().getValue(); i++) {
         int fine_rat = fine_ratio(i);
         int coarsen_rat = coarsen_ratio(i);
         if (fine_rat > 0) {
            if (fine_rat >= coarsen_rat) {
               d_ratio_to_level_zero(i) = fine_rat / coarsen_rat;
            } else {
               d_ratio_to_level_zero(i) =
                  -(tbox::MathUtilities<int>::Abs(coarsen_rat / fine_rat));
            }
         } else {
            d_ratio_to_level_zero(i) =
               -(tbox::MathUtilities<int>::Abs(fine_rat * coarsen_rat));
         }
      }

   } else {

      d_geometry = coarse_grid_geom;

      d_ratio_to_level_zero = fine_level->d_ratio_to_level_zero;
   }

   /*
    * Set global box array and index space for level based on coarsening
    * of fine level information.
    */

   d_boxes = fine_level->d_boxes;
   d_boxes.coarsen(coarsen_ratio);

   /*
    * Set coarse mapped_box_level to be the coarsened version of fine mapped_box_level.
    *
    * NOTE: Some parts of SAMRAI (CoarsenSchedule in particular)
    * assumes that the mapped_box identities are the same between the
    * fine and coarsened levels.
    */
   const BoxLevel& fine_mapped_box_level =
      *fine_level->d_mapped_box_level;
   d_mapped_box_level.reset(
      new BoxLevel(
         d_ratio_to_level_zero,
         d_geometry,
         fine_mapped_box_level.getMPI()));
   fine_level->d_mapped_box_level->coarsenBoxes(
      *d_mapped_box_level,
      coarsen_ratio,
      d_ratio_to_level_zero);
   d_mapped_box_level->finalize();
   d_local_number_patches = fine_level->getNumberOfPatches();
   d_number_blocks = fine_level->d_number_blocks;

   d_physical_domain.resizeArray(d_number_blocks);
   for (int nb = 0; nb < d_number_blocks; nb++) {
      d_physical_domain[nb] = fine_level->d_physical_domain[nb];
      d_physical_domain[nb].coarsen(coarsen_ratio);
   }

   /*
    * Allocate arrays of patches and patch information.  Then, allocate and
    * initialize patch objects.  Finally, set patch geometry and remaining
    * domain information.
    */

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   for (RealBoxConstIterator ni(mapped_boxes.realBegin());
        ni != mapped_boxes.realEnd(); ++ni) {
      const Box& mapped_box = *ni;
      const BoxId& mapped_box_id = mapped_box.getId();
      d_patches[mapped_box_id] = d_factory->allocate(mapped_box, d_descriptor);
      d_patches[mapped_box_id]->setPatchLevelNumber(d_level_number);
      d_patches[mapped_box_id]->setPatchInHierarchy(d_in_hierarchy);
   }

   d_boundary_boxes_created = false;

   std::map<BoxId, PatchGeometry::TwoDimBool> touches_regular_bdry;
   std::map<BoxId, PatchGeometry::TwoDimBool> touches_periodic_bdry;

   for (iterator ip(fine_level->begin()); ip != fine_level->end(); ip++) {
      boost::shared_ptr<PatchGeometry> fine_pgeom((*ip)->getPatchGeometry());

      /* If map does not contain values create them */
      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_regular_bdry(
         touches_regular_bdry.find(ip->getBox().getId()));
      if (iter_touches_regular_bdry == touches_regular_bdry.end()) {
         iter_touches_regular_bdry = touches_regular_bdry.insert(
               iter_touches_regular_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      std::map<BoxId,
               PatchGeometry::TwoDimBool>::iterator iter_touches_periodic_bdry(
         touches_periodic_bdry.find(ip->getBox().getId()));
      if (iter_touches_periodic_bdry == touches_periodic_bdry.end()) {
         iter_touches_periodic_bdry = touches_periodic_bdry.insert(
               iter_touches_periodic_bdry,
               std::pair<BoxId, PatchGeometry::TwoDimBool>(ip->getBox().getId(),
                  PatchGeometry::TwoDimBool(getDim())));
      }

      PatchGeometry::TwoDimBool&
      touches_regular_bdry_ip((*iter_touches_regular_bdry).second);

      for (int axis = 0; axis < getDim().getValue(); axis++) {
         for (int side = 0; side < 2; side++) {
            touches_regular_bdry_ip(axis, side) =
               fine_pgeom->getTouchesRegularBoundary(axis, side);
         }
      }
   }

   d_geometry->setGeometryOnPatches(
      *this,
      d_ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry,
      defer_boundary_box_creation);

   if (!defer_boundary_box_creation) {
      d_boundary_boxes_created = true;
   }

}

/*
 * ************************************************************************
 * ************************************************************************
 */
void
PatchLevel::getBoxes(
   BoxContainer& boxes,
   const BlockId& block_id) const
{
   if (!d_has_globalized_data) {
      initializeGlobalizedBoxLevel();
   }

   boxes.clear();
   const BoxContainer& global_mapped_boxes =
      d_mapped_box_level->getGlobalizedVersion().getGlobalBoxes();

   for (BoxContainerSingleBlockIterator gi(global_mapped_boxes.begin(block_id));
        gi != global_mapped_boxes.end(block_id); ++gi) {
      boxes.pushBack(*gi);
   }
}

/*
 * ************************************************************************
 *
 *  Check that class version and restart file number are the same.  If
 *  so, read in data from database and build patch level from data.
 *
 * ************************************************************************
 */

void
PatchLevel::getFromDatabase(
   const boost::shared_ptr<tbox::Database>& database,
   const ComponentSelector& component_selector)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("HIER_PATCH_LEVEL_VERSION");
   if (ver != HIER_PATCH_LEVEL_VERSION) {
      TBOX_ERROR("PatchLevel::getFromDatabase() error...\n"
         << "   Restart file version different than class version.");
   }

   if (database->keyExists("d_boxes")) {
      d_boxes = database->getDatabaseBoxArray("d_boxes");
   }

   int* temp_ratio = &d_ratio_to_level_zero[0];
   database->getIntegerArray(
      "d_ratio_to_level_zero",
      temp_ratio,
      getDim().getValue());

   d_number_blocks = database->getInteger("d_number_blocks");

   d_physical_domain.resizeArray(d_number_blocks);
   for (int nb = 0; nb < d_number_blocks; nb++) {
      std::string domain_name = "d_physical_domain_"
         + tbox::Utilities::blockToString(nb);
      d_physical_domain[nb] = database->getDatabaseBoxArray(domain_name);
      for (BoxContainer::iterator bi = d_physical_domain[nb].begin();
           bi != d_physical_domain[nb].end(); ++bi) {
         bi->setBlockId(BlockId(nb));
      }
   }

   d_level_number = database->getInteger("d_level_number");
   d_next_coarser_level_number =
      database->getInteger("d_next_coarser_level_number");
   d_in_hierarchy = database->getBool("d_in_hierarchy");

   temp_ratio = &d_ratio_to_coarser_level[0];
   database->getIntegerArray(
      "d_ratio_to_coarser_level",
      temp_ratio,
      getDim().getValue());

   /*
    * Put local patches in database.
    */

   boost::shared_ptr<tbox::Database> mbl_database(
      database->getDatabase("mapped_box_level"));
   boost::shared_ptr<BoxLevel> mapped_box_level(
      boost::make_shared<BoxLevel>(getDim()));
   boost::shared_ptr<const BaseGridGeometry> grid_geometry(getGridGeometry());
   mapped_box_level->getFromDatabase(*mbl_database, grid_geometry);
   d_mapped_box_level = mapped_box_level;

   d_patches.clear();

   const BoxContainer& mapped_boxes = d_mapped_box_level->getBoxes();
   for (RealBoxConstIterator ni(mapped_boxes.realBegin());
        ni != mapped_boxes.realEnd(); ++ni) {
      const Box& mapped_box = *ni;
      const LocalId& local_id = mapped_box.getLocalId();
      const BoxId& mapped_box_id = mapped_box.getId();

      std::string patch_name = "level_" + tbox::Utilities::levelToString(
            d_level_number)
         + "-patch_" + tbox::Utilities::patchToString(local_id.getValue())
         + "-block_"
         + tbox::Utilities::blockToString(
            mapped_box.getBlockId().getBlockValue());
      if (!(database->isDatabase(patch_name))) {
         TBOX_ERROR("PatchLevel::getFromDatabase() error...\n"
            << "   patch name " << patch_name
            << " not found in database" << std::endl);
      }

      boost::shared_ptr<Patch>& patch(d_patches[mapped_box_id]);
      patch = d_factory->allocate(mapped_box, d_descriptor);
      patch->setPatchLevelNumber(d_level_number);
      patch->setPatchInHierarchy(d_in_hierarchy);
      patch->getFromDatabase(
         database->getDatabase(patch_name),
         component_selector);
   }

}

/*
 * ************************************************************************
 *
 *  Write out class version number and patch_level data members to the
 *  database, then has each patch on the local processor write itself
 *  to the database.   The following are written out to the database:
 *  d_physical_domain, d_ratio_to_level_zero, d_boxes, d_mapping,
 *  d_global_number_patches, d_level_number, d_next_coarser_level_number,
 *  d_in_hierarchy, d_patches[].
 *  The database key for all data members except for d_patches is
 *  the same as the variable name.  For the patches, the database keys
 *  are "level_Xpatch_Y" where X is the level number and Y is the index
 *  position of the patch in the patch in d_patches.
 *
 * ************************************************************************
 */
void
PatchLevel::putUnregisteredToDatabase(
   const boost::shared_ptr<tbox::Database>& database,
   const ComponentSelector& patchdata_write_table) const
{
   TBOX_ASSERT(database);

   database->putInteger("HIER_PATCH_LEVEL_VERSION", HIER_PATCH_LEVEL_VERSION);

   database->putBool("d_is_patch_level", true);

   tbox::Array<tbox::DatabaseBox> temp_boxes = d_boxes;
   if (temp_boxes.getSize() > 0) {
      database->putDatabaseBoxArray("d_boxes", temp_boxes);
   }

   // database->putInteger("d_global_number_patches",d_global_number_patches);

   // database->putIntegerArray("d_mapping", d_mapping.getProcessorMapping());

   const int* temp_ratio_to_level_zero = &d_ratio_to_level_zero[0];
   database->putIntegerArray("d_ratio_to_level_zero",
      temp_ratio_to_level_zero, getDim().getValue());

   database->putInteger("d_number_blocks", d_number_blocks);

   for (int nb = 0; nb < d_number_blocks; nb++) {
      tbox::Array<tbox::DatabaseBox> temp_domain = d_physical_domain[nb];
      std::string domain_name = "d_physical_domain_"
         + tbox::Utilities::blockToString(nb);
      database->putDatabaseBoxArray(domain_name, temp_domain);
   }
   database->putInteger("d_level_number", d_level_number);
   database->putInteger("d_next_coarser_level_number",
      d_next_coarser_level_number);
   database->putBool("d_in_hierarchy", d_in_hierarchy);

   const int* temp_ratio_to_coarser_level = &d_ratio_to_coarser_level[0];
   database->putIntegerArray("d_ratio_to_coarser_level",
      temp_ratio_to_coarser_level, getDim().getValue());

   /*
    * Put local patches in database.
    */

   boost::shared_ptr<tbox::Database> mbl_database(
      database->putDatabase("mapped_box_level"));
   d_mapped_box_level->putUnregisteredToDatabase(mbl_database);

   for (iterator ip(begin()); ip != end(); ip++) {

      std::string patch_name = "level_" + tbox::Utilities::levelToString(
            d_level_number)
         + "-patch_"
         + tbox::Utilities::patchToString(ip->getLocalId().getValue())
         + "-block_"
         + tbox::Utilities::blockToString(
            ip->getBox().getBlockId().getBlockValue());

      ip->putUnregisteredToDatabase(
         database->putDatabase(patch_name),
         patchdata_write_table);
   }

}

int
PatchLevel::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int depth)
{
   int npatch = getGlobalNumberOfPatches();

// Disable Intel warnings on conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#pragma warning (disable:857)
#endif

   os << border << "Local/Global number of patches and cells = "
      << getLocalNumberOfPatches() << "/" << getGlobalNumberOfPatches() << "  "
      << getLocalNumberOfCells() << "/" << getGlobalNumberOfCells() << "\n";
   os << getBoxLevel()->format(border, 2) << std::endl;

   if (depth > 0) {
      for (iterator pi(begin()); pi != end(); pi++) {
         const boost::shared_ptr<Patch>& patch = *pi;
         os << border << "Patch " << patch->getLocalId() << '/' << npatch << "\n";
         patch->recursivePrint(os, border + "\t", depth - 1);

      }
   }
   return 0;
}

/*
 *************************************************************************
 * Private utility function to gather and store globalized data, if needed.
 *************************************************************************
 */
void
PatchLevel::initializeGlobalizedBoxLevel() const
{
   if (!d_has_globalized_data) {

      const BoxLevel& globalized_mapped_box_level(
         d_mapped_box_level->getGlobalizedVersion());

      const int nboxes = globalized_mapped_box_level.getGlobalNumberOfBoxes();
      d_boxes.clear();
      d_mapping.setMappingSize(nboxes);

      /*
       * Backward compatibility with things requiring global sequential
       * indices (such as the VisIt writer) is provided by the implicit
       * ordering of the mapped_boxes in the nested loops below.
       *
       * Due to this necessary renumbering, the patch number obtained
       * by the PatchLevel::Iterator does not correspond to the
       * global sequential index.
       */
      int count = 0;
      const BoxContainer& mapped_boxes =
         globalized_mapped_box_level.getGlobalBoxes();
      for (RealBoxConstIterator ni(mapped_boxes.realBegin());
           ni != mapped_boxes.realEnd();
           ++ni) {
         d_mapping.setProcessorAssignment(count, ni->getOwnerRank());
         d_boxes.pushBack(*ni);
         ++count;
      }

      d_has_globalized_data = true;
   }
}

/*
 * ************************************************************************
 * ************************************************************************
 */

void
PatchLevel::initializeCallback()
{
   t_level_constructor = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::level_constructor");
   t_constructor_setup = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_setup");
   t_constructor_phys_domain = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_phys_domain");
   t_constructor_touch_boundaries = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_touch_boundaries");
   t_constructor_set_geometry = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::set_geometry");
   t_constructor_compute_shifts = tbox::TimerManager::getManager()->
      getTimer("hier::PatchLevel::constructor_compute_shifts");
}

/*
 ***************************************************************************
 *
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 *
 ***************************************************************************
 */

void
PatchLevel::finalizeCallback()
{
   t_level_constructor.reset();
   t_constructor_setup.reset();
   t_constructor_phys_domain.reset();
   t_constructor_touch_boundaries.reset();
   t_constructor_set_geometry.reset();
   t_set_patch_touches.reset();
   t_constructor_compute_shifts.reset();
}

/*
 *************************************************************************
 * Copy constructor.
 *************************************************************************
 */
PatchLevel::Iterator::Iterator(
   const PatchLevel::Iterator& r):
   d_iterator(r.d_iterator),
   d_patches(NULL /* Unused since not backward compatibility not needed */)
{
}

// Support for backward compatible interface by new PatchLevel::Iterator.
PatchLevel::Iterator::Iterator(
   const PatchLevel* patch_level,
   bool begin):
   d_iterator(begin ? patch_level->d_patches.begin() :
                      patch_level->d_patches.end()),
   d_patches(&patch_level->d_patches)
{
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
