/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelBorderFillPattern_C
#define included_xfer_PatchLevelBorderFillPattern_C

#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
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

PatchLevelBorderFillPattern::PatchLevelBorderFillPattern():
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

PatchLevelBorderFillPattern::~PatchLevelBorderFillPattern()
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
PatchLevelBorderFillPattern::computeFillBoxesAndNeighborhoodSets(
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
    * To get the level border, grow each patch box and remove
    * the level from it.
    */
   hier::LocalId last_id = dst_mapped_box_level.getLastLocalId();
   for (hier::RealBoxConstIterator ni(dst_mapped_boxes.realBegin());
        ni != dst_mapped_boxes.realEnd(); ++ni) {
      const hier::Box& dst_mapped_box = *ni;
      hier::BoxContainer fill_boxes(
         hier::Box::grow(dst_mapped_box, fill_ghost_width));
      hier::Connector::ConstNeighborhoodIterator nabrs =
         dst_to_dst.find(dst_mapped_box.getId());
      for (hier::Connector::ConstNeighborIterator na = dst_to_dst.begin(nabrs);
           na != dst_to_dst.end(nabrs); ++na) {
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

      if (!fill_boxes.isEmpty()) {
         d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
               fill_boxes.size());
         hier::Connector::NeighborhoodIterator base_box_itr =
            dst_to_fill.makeEmptyLocalNeighborhood(dst_mapped_box.getId());
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
PatchLevelBorderFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_mapped_box_level);
   NULL_USE(src_to_dst);
   NULL_USE(fill_ghost_width);
   NULL_USE(dst_fill_boxes_on_src_proc);
   if (!needsToCommunicateDestinationFillBoxes()) {
      TBOX_ERROR(
         "PatchLevelBorderFillPattern cannot compute destination:\n"
         << "fill boxes on the source processor.\n");
   }
}

bool
PatchLevelBorderFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return true;
}

bool
PatchLevelBorderFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return false;
}

bool
PatchLevelBorderFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

bool
PatchLevelBorderFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

int
PatchLevelBorderFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

}
}
#endif
