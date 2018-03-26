/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelBorderAndInteriorFillPattern_C
#define included_xfer_PatchLevelBorderAndInteriorFillPattern_C

#include "SAMRAI/xfer/PatchLevelBorderAndInteriorFillPattern.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */

PatchLevelBorderAndInteriorFillPattern::PatchLevelBorderAndInteriorFillPattern():
   d_max_fill_boxes(0)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

PatchLevelBorderAndInteriorFillPattern::~PatchLevelBorderAndInteriorFillPattern()
{
}

/*
 *************************************************************************
 *
 * computeFillBoxesAndNeighborhoodSets
 *
 *************************************************************************
 */
void
PatchLevelBorderAndInteriorFillPattern::computeFillBoxesAndNeighborhoodSets(
   hier::BoxLevel& fill_mapped_boxes,
   hier::Connector& dst_to_fill,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& dst_to_dst,
   const hier::Connector& dst_to_src,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_to_src);
   NULL_USE(src_to_dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   const hier::BoxContainer& dst_mapped_boxes =
      dst_mapped_box_level.getBoxes();

   /*
    * Grow each patch box and remove the level from it, except the
    * patch box itself.  (Do not fill ghost cells that are normally
    * filled by same mapped_box_level.  Do fill ghost cells that are
    * normally filled by coarser mapped_box_level.)
    */
   hier::LocalId last_id = dst_mapped_box_level.getLastLocalId();
   for (hier::BoxContainer::const_iterator ni = dst_mapped_boxes.begin();
        ni != dst_mapped_boxes.end(); ++ni) {

      const hier::BoxId& gid(ni->getId());
      const hier::Box& dst_mapped_box =
         *dst_mapped_box_level.getBox(gid);
      hier::BoxContainer fill_boxes(
         hier::Box::grow(dst_mapped_box, fill_ghost_width));
      hier::Connector::ConstNeighborhoodIterator nabrs =
         dst_to_dst.find(dst_mapped_box.getId());

      for (hier::Connector::ConstNeighborIterator na = dst_to_dst.begin(nabrs);
           na != dst_to_dst.end(nabrs); ++na) {
         if (!ni->isSpatiallyEqual(*na)) {
            if (dst_mapped_box.getBlockId() == na->getBlockId()) {
               fill_boxes.removeIntersections(*na);
            } else {

               boost::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
                  dst_mapped_box_level.getGridGeometry());

               const hier::BlockId& dst_block_id = dst_mapped_box.getBlockId();
               const hier::BlockId& nbr_block_id = na->getBlockId();

               TBOX_ASSERT(grid_geometry->areNeighbors(dst_block_id,
                     nbr_block_id));

               hier::Transformation::RotationIdentifier rotation =
                  grid_geometry->getRotationIdentifier(dst_block_id,
                     nbr_block_id);
               hier::IntVector offset(
                  grid_geometry->getOffset(dst_block_id, nbr_block_id));

               offset *= (dst_mapped_box_level.getRefinementRatio());

               hier::Transformation transformation(rotation, offset,
                                                   nbr_block_id, dst_block_id);

               hier::Box nbr_box(*na);
               transformation.transform(nbr_box);

               fill_boxes.removeIntersections(nbr_box);

            }
         }
      }

      if (!fill_boxes.isEmpty()) {
         d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
               fill_boxes.size());

         hier::Connector::NeighborhoodIterator base_box_itr =
            dst_to_fill.makeEmptyLocalNeighborhood(gid);
         for (hier::BoxContainer::iterator li(fill_boxes);
              li != fill_boxes.end(); ++li) {
            hier::Box fill_mapped_box(*li,
                                      ++last_id,
                                      dst_mapped_box.getOwnerRank());
            TBOX_ASSERT(fill_mapped_box.getBlockId() ==
                        dst_mapped_box.getBlockId());
            fill_mapped_boxes.addBoxWithoutUpdate(fill_mapped_box);
            dst_to_fill.insertLocalNeighbor(fill_mapped_box, base_box_itr);
         }
      }
   }
   fill_mapped_boxes.finalize();
}

void
PatchLevelBorderAndInteriorFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_fill_boxes_on_src_proc);
   NULL_USE(dst_mapped_box_level);
   NULL_USE(src_to_dst);
   NULL_USE(fill_ghost_width);
   if (!needsToCommunicateDestinationFillBoxes()) {
      TBOX_ERROR(
         "PatchLevelBorderAndInteriorFillPattern cannot compute destination:\n"
         << "fill boxes on the source processor.\n");
   }
}

bool
PatchLevelBorderAndInteriorFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return true;
}

bool
PatchLevelBorderAndInteriorFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return true;
}

int
PatchLevelBorderAndInteriorFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

bool
PatchLevelBorderAndInteriorFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

bool
PatchLevelBorderAndInteriorFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

}
}
#endif
