/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   For describing coarse-fine boundary interfaces
 *
 ************************************************************************/

#ifndef included_hier_CoarseFineBoundary_C
#define included_hier_CoarseFineBoundary_C

#include "SAMRAI/hier/CoarseFineBoundary.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/PatchLevel.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

CoarseFineBoundary::CoarseFineBoundary(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_initialized(1, false)
{
}

CoarseFineBoundary::CoarseFineBoundary(
   const CoarseFineBoundary& rhs):
   d_dim(rhs.d_dim),
   d_initialized(1, false),
   d_boundary_boxes(rhs.d_boundary_boxes)
{
   /*
    * This needs to be written this way since STL vector for bools
    * causes uninitialized memory reads because it is poorly implemented.
    * So the vector is initialized and then copied.
    */
   d_initialized = rhs.d_initialized;
}

CoarseFineBoundary::CoarseFineBoundary(
   const PatchHierarchy& hierarchy,
   int level_num,
   const IntVector& max_ghost_width):
   d_dim(max_ghost_width.getDim()),
   d_initialized(hierarchy.getGridGeometry()->getNumberBlocks(), false)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, max_ghost_width);
   TBOX_ASSERT(max_ghost_width > IntVector(d_dim, -1));

   if (hierarchy.getGridGeometry()->getNumberBlocks() == 1) {
      const PatchLevel& level =
         dynamic_cast<const PatchLevel&>(*hierarchy.getPatchLevel(level_num));
      const Connector& level_to_level =
         level.getBoxLevel()->getPersistentOverlapConnectors().
         findOrCreateConnector(
            *level.getBoxLevel(),
            max_ghost_width);
      const Connector& level_to_domain =
         level.getBoxLevel()->getPersistentOverlapConnectors().
         findOrCreateConnector(
            hierarchy.getDomainBoxLevel(),
            max_ghost_width);
      computeFromLevel(level,
         level_to_domain,
         level_to_level,
         max_ghost_width);
   } else {
      const PatchLevel& level =
         dynamic_cast<const PatchLevel&>(*hierarchy.getPatchLevel(level_num));
      const PatchLevel& level0 =
         dynamic_cast<const PatchLevel&>(*hierarchy.getPatchLevel(0));
      computeFromLevel(
         level,
         level0,
         max_ghost_width);
   }

}

CoarseFineBoundary::CoarseFineBoundary(
   const PatchLevel& level,
   const Connector& mapped_box_level_to_domain,
   const Connector& mapped_box_level_to_self,
   const IntVector& max_ghost_width):
   d_dim(max_ghost_width.getDim()),
   d_initialized(1, false)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, max_ghost_width);
   TBOX_ASSERT(max_ghost_width > IntVector(d_dim, -1));

   computeFromLevel(level,
      mapped_box_level_to_domain,
      mapped_box_level_to_self,
      max_ghost_width);

}

CoarseFineBoundary::~CoarseFineBoundary()
{
}

/*
 ************************************************************************
 * Use grid_geometry.computeBoundaryGeometry function,
 * setting up the arguments in a way that will generate
 * the coarse-fine boundary (instead of the domain boundary).
 ************************************************************************
 */
void
CoarseFineBoundary::computeFromLevel(
   const PatchLevel& level,
   const Connector& mapped_box_level_to_domain,
   const Connector& mapped_box_level_to_self,
   const IntVector& max_ghost_width)
{
//this is single block version.
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, max_ghost_width);

   clear();

   const BoxLevel& mapped_box_level = *level.getBoxLevel();
   const IntVector& ratio = level.getRatioToLevelZero();

   boost::shared_ptr<BaseGridGeometry> grid_geometry (level.getGridGeometry());

   /*
    * Get the domain's periodic shift.
    */
   const IntVector periodic_shift(grid_geometry->getPeriodicShift(ratio));

   bool is_periodic = false;
   for (int i = 0; i < d_dim.getValue(); ++i) {
      is_periodic = is_periodic || periodic_shift(i);
   }

   /*
    * Here we add some boxes outside of non-periodic boundaries to the
    * adjusted level.  For each patch that touches a regular boundary,
    * grow the patch box (and any periodic images of the patch box) by
    * the max ghost width.  Remove intersections with the periodic
    * adjusted physical domain.  Add what remains to the adjusted level.
    *
    * This will ensure that ensuing call to create boundary boxes will not
    * create boundary boxes at the locations where the level touches a
    * non-periodic physical boundary, but only where there is a coarse-fine
    * interface in the domain interior (A periodic boundary is considered
    * part of the domain interior for this purpose).
    */

   /*
    * Build a fake domain in fake_domain_list.
    *
    * The fake domain should be such that when fed to computeBoundaryBoxesOnLevel,
    * the coarse-fine boundaries are computed rather than the physical boundary.
    * computeBoundaryBoxesOnLevel defines boundaries of a patch to be the
    * parts of the grown patch box that lie outside the "domain".  So se make
    * the fake domain be everywhere there is NOT a coarse-fine boundary--or
    * everywhere there IS a physical boundary or a fine-boundary.
    */
   tbox::Array<BoxContainer> fake_domain(1);
   BoxContainer& fake_domain_list = fake_domain[0];

   // Every mapped_box should connect to the domain mapped_box_level.
   TBOX_ASSERT(mapped_box_level_to_domain.getLocalNumberOfNeighborSets() ==
               static_cast<int>(mapped_box_level.getLocalNumberOfBoxes()));

   // Add physical boundaries to the fake domain.
   for (Connector::ConstNeighborhoodIterator ei = mapped_box_level_to_domain.begin();
        ei != mapped_box_level_to_domain.end(); ++ei) {
      const Box& mapped_box = *mapped_box_level.getBoxStrict(*ei);
      BoxContainer refined_domain_nabrs;
      for (Connector::ConstNeighborIterator ni = mapped_box_level_to_domain.begin(ei);
           ni != mapped_box_level_to_domain.end(ei); ++ni) {
         refined_domain_nabrs.insert(refined_domain_nabrs.end(), *ni);
      }
      refined_domain_nabrs.refine(ratio);
      Box box = mapped_box;
      box.grow(max_ghost_width);
      BoxContainer physical_boundary_portion(box);
      physical_boundary_portion.removeIntersections(refined_domain_nabrs);
      fake_domain_list.spliceBack(physical_boundary_portion);
   }

   // Add fine-fine boundaries to the fake domain.
   TBOX_ASSERT(mapped_box_level_to_self.getConnectorWidth() >=
      IntVector::getOne(d_dim));
   mapped_box_level_to_self.getLocalNeighbors(fake_domain_list);

   /*
    * Call BaseGridGeometry::computeBoundaryGeometry with arguments contrived
    * such that they give the coarse-fine boundaries instead of the domain
    * boundaries.  The basic algorithm used by
    * BaseGridGeometry::computeBoundaryGeometry is
    * 1. grow boxes by ghost width
    * 2. remove intersection with domain
    * 3. reorganize and classify resulting boxes
    *
    * This is how we get BaseGridGeometry::computeBoundaryGeometry to
    * compute the coarse-fine boundary instead of the physical boundary.
    *
    * Since we handle the periodic boundaries ourselves, do not treat
    * them differently from regular boundaries.  State that all boundaries
    * are non-periodic boundaries.
    *
    * Send the periodic-adjusted level boxes as the domain for the
    * remove-intersection-with-domain operation.  This causes that
    * operation to remove non-coarse-fine (that is, fine-fine) boxes
    * along the periodic boundaries, leaving the coarse-fine boundary
    * boxes.
    *
    * Send the periodic-adjusted domain for the limit-domain intersect
    * operation.  This removes the boundaries that are on the non-periodic
    * boundaries, which is what we want because there is no possibility
    * of a coarse-fine boundary there.
    */
   bool do_all_patches = true;
   const IntVector use_periodic_shift(d_dim, 0);
   grid_geometry->computeBoundaryBoxesOnLevel(
      d_boundary_boxes,
      level,
      use_periodic_shift,
      max_ghost_width,
      fake_domain,
      do_all_patches);

   d_initialized[0] = true;
}

void
CoarseFineBoundary::computeFromLevel(
   const PatchLevel& level,
   const PatchLevel& level0,
   const IntVector& max_ghost_width)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, max_ghost_width);

   /*
    * Get all the boxes on level and level0.  These will be used later.
    */
   const BoxContainer& all_boxes_on_level =
      level.getBoxLevel()->getGlobalizedVersion().getGlobalBoxes();
   const BoxContainer& all_boxes_on_level0 =
      level0.getBoxLevel()->getGlobalizedVersion().getGlobalBoxes();

   /*
    * Get the dimension and number of blocks from the grid geometry.
    */
   boost::shared_ptr<BaseGridGeometry> grid_geometry(level.getGridGeometry());
   int nblocks = grid_geometry->getNumberBlocks();

   tbox::Array<BoxContainer> adjusted_level_domain(nblocks);

   /*
    * Loop over each block.
    */
   for (int i = 0; i < nblocks; ++i) {

      clear();

      BlockId block_id(i);
      /*
       * Construct an iterator which filters only level's boxes in this block.
       */
      BoxContainerSingleBlockIterator itr(
         all_boxes_on_level.begin(block_id));

      /*
       * Only do work if there any boxes in this block.
       */
      if (itr != all_boxes_on_level.end(block_id)) {

         /*
          * Construct the array of boxes on level and level0 in this block.
          */
         BoxContainer level_domain(all_boxes_on_level, block_id);
         BoxContainer phys_domain(all_boxes_on_level0, block_id);

         const IntVector& ratio = level.getRatioToLevelZero();
         phys_domain.refine(ratio);

         /*
          * Create a pseudo-domain -- the union of the physical domain boxes
          * of the current block with the physical domain boxes of all of the
          * current block's neighbors.  These are all represented in the
          * current block's index space, and refined by level's refinement
          * ratio.
          */

         BoxContainer pseudo_domain(phys_domain);
         pseudo_domain.unorder();
         const std::list<BaseGridGeometry::Neighbor>& nbr_list =
            grid_geometry->getNeighbors(block_id);
         for (std::list<BaseGridGeometry::Neighbor>::const_iterator ni =
              nbr_list.begin();
              ni != nbr_list.end(); ni++) {

            BoxContainer neighbor_domain(ni->getTransformedDomain());
            neighbor_domain.refine(ratio);

            pseudo_domain.spliceFront(neighbor_domain);

         }

         /*
          * Make a list containing the level boxes for the current block,
          * then add more boxes as a buffer around physical domain boundaries.
          * This prevents physical boundaries from being identified as
          * coarse-fine boundaries.
          */

         adjusted_level_domain[i] = level_domain;
         adjusted_level_domain[i].unorder();

         for (PatchLevel::iterator p(level.begin()); p != level.end(); ++p) {
            if ((*p)->getBox().getBlockId() == i &&
                (*p)->getPatchGeometry()->getTouchesRegularBoundary()) {

               const Box& patch_box = (*p)->getBox();

               BoxContainer no_shift_boxes(patch_box);
               no_shift_boxes.grow(max_ghost_width);
               no_shift_boxes.removeIntersections(pseudo_domain);
               adjusted_level_domain[i].spliceFront(no_shift_boxes);
            }
         }

         /*
          * Add buffer of boxes that exist on the current level across
          * block boundaries from the current block.  This prevents block
          * boundaries from being identified as coarse-fine boundaries when
          * they are not.
          */
         for (std::list<BaseGridGeometry::Neighbor>::const_iterator ni =
              nbr_list.begin();
              ni != nbr_list.end(); ni++) {

            /*
             * Construct the array of boxes on level in this neighbor's block.
             */
            BlockId nbr_block_id(ni->getBlockId());
            BoxContainer neighbor_boxes(all_boxes_on_level, nbr_block_id);

            if (neighbor_boxes.size()) {
               neighbor_boxes.unorder();
               grid_geometry->transformBoxContainer(neighbor_boxes,
                  ratio,
                  block_id,
                  nbr_block_id);

               BoxContainer neighbor_boxes_to_add(phys_domain);
               neighbor_boxes_to_add.unorder();
               neighbor_boxes_to_add.grow(max_ghost_width);

               neighbor_boxes_to_add.intersectBoxes(neighbor_boxes);

               adjusted_level_domain[i].spliceFront(neighbor_boxes_to_add);
            }
         }

      }

      d_initialized[i] = true;
   }

   d_boundary_boxes.clear();

   /*
    * Call BaseGridGeometry::computeBoundaryGeometry with arguments contrived
    * such that they give the coarse-fine boundaries instead of the
    * domain boundaries.  The basic algorithm used by
    * BaseGridGeometry::computeBoundaryGeometry is
    * 1. grow boxes by ghost width
    * 2. remove intersection with domain
    * 3. reorganize and classify resulting boxes
    *
    * This is how we get BaseGridGeometry::computeBoundaryGeometry to
    * compute the coarse-fine boundary instead of the physical boundary.
    *
    * Send the adjusted level boxes as the domain for the
    * remove-intersection-with-domain operation.  This causes that
    * operation to remove non-coarse-fine (that is, fine-fine) boxes
    * along the periodic boundaries, leaving the coarse-fine boundary
    * boxes.
    *
    * Send the adjusted domain for the limit-domain intersect
    * operation.  This removes the boundaries that are on the physical
    * boundaries, which is what we want because there is no possibility
    * of a coarse-fine boundary there.
    */
   bool do_all_patches = true;
   IntVector use_periodic_shift(d_dim, 0);
   grid_geometry->computeBoundaryBoxesOnLevel(
      d_boundary_boxes,
      level,
      use_periodic_shift,
      max_ghost_width,
      adjusted_level_domain,
      do_all_patches);

}

const tbox::Array<BoundaryBox>&
CoarseFineBoundary::getBoundaries(
   const GlobalId& global_id,
   const int boundary_type,
   const BlockId& block_id) const
{
   const int& block_num = block_id.getBlockValue();
   if (!d_initialized[block_num]) {
      TBOX_ERROR("The boundary boxes have not been computed.");
   }

   BoxId mapped_box_id(global_id);
   std::map<BoxId, PatchBoundaries>::const_iterator
      mi = d_boundary_boxes.find(mapped_box_id);
   TBOX_ASSERT(mi != d_boundary_boxes.end());
   return (*mi).second[boundary_type - 1];
}

void
CoarseFineBoundary::printClassData(
   std::ostream& os) const {
   os << "\nCoarseFineBoundary::printClassData...";
   for (std::map<BoxId, PatchBoundaries>::const_iterator
        mi = d_boundary_boxes.begin(); mi != d_boundary_boxes.end(); ++mi) {
      os << "\n         patch " << (*mi).first;
      for (unsigned int btype = 0; btype < d_dim.getValue(); ++btype) {
         os << "\n                type " << btype;
         const tbox::Array<BoundaryBox>
         & array_of_boxes = (*mi).second[btype];
         int num_boxes = array_of_boxes.getSize();
         int bn;
         for (bn = 0; bn < num_boxes; ++bn) {
            os << "\n                           box "
               << bn << "/" << num_boxes << ":";
            os << array_of_boxes[bn].getBox();
         }
      }
   }
   os << "\n";
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
