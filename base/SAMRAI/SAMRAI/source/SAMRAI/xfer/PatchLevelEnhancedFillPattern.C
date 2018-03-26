/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelEnhancedFillPattern_C
#define included_xfer_PatchLevelEnhancedFillPattern_C

#include "SAMRAI/xfer/PatchLevelEnhancedFillPattern.h"
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

PatchLevelEnhancedFillPattern::PatchLevelEnhancedFillPattern():
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

PatchLevelEnhancedFillPattern::~PatchLevelEnhancedFillPattern()
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
PatchLevelEnhancedFillPattern::computeFillBoxesAndNeighborhoodSets(
   hier::BoxLevel& fill_mapped_boxes,
   hier::Connector& dst_to_fill,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& dst_to_dst,
   const hier::Connector& dst_to_src,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(dst_to_dst);
   NULL_USE(dst_to_src);
   NULL_USE(src_to_dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   boost::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
      dst_mapped_box_level.getGridGeometry());

   const hier::BoxContainer& dst_mapped_boxes =
      dst_mapped_box_level.getBoxes();

   hier::LocalId last_id = dst_mapped_box_level.getLastLocalId();
   for (hier::RealBoxConstIterator ni(dst_mapped_boxes.realBegin());
        ni != dst_mapped_boxes.realEnd(); ++ni) {
      const hier::Box& dst_mapped_box = *ni;
      const hier::BoxId& dst_mapped_box_id = dst_mapped_box.getId();
      hier::BoxContainer fill_boxes(
         hier::Box::grow(dst_mapped_box, fill_ghost_width));

      const std::list<hier::BaseGridGeometry::Neighbor>& neighbors =
         grid_geometry->getNeighbors(dst_mapped_box.getBlockId());

      hier::BoxContainer constructed_fill_boxes;

      hier::Connector::NeighborhoodIterator base_box_itr =
         dst_to_fill.findLocal(dst_mapped_box_id);
      bool has_base_box = base_box_itr != dst_to_fill.end();

      for (std::list<hier::BaseGridGeometry::Neighbor>::const_iterator ni =
           neighbors.begin();
           ni != neighbors.end(); ni++) {

         if (ni->isSingularity()) {

            hier::BoxContainer encon_boxes(ni->getTransformedDomain());
            encon_boxes.refine(dst_mapped_box_level.getRefinementRatio());
            encon_boxes.intersectBoxes(fill_boxes);
            encon_boxes.removeIntersections(constructed_fill_boxes);

            if (encon_boxes.size()) {

               if (!has_base_box) {
                  base_box_itr = dst_to_fill.makeEmptyLocalNeighborhood(
                     dst_mapped_box_id);
                  has_base_box = true;
               }
               for (hier::BoxContainer::iterator ei(encon_boxes);
                    ei != encon_boxes.end(); ei++) {

                  hier::Box fill_mapped_box(
                     *ei,
                     ++last_id,
                     dst_mapped_box.getOwnerRank());

                  TBOX_ASSERT(fill_mapped_box.getBlockId() ==
                              dst_mapped_box.getBlockId());

                  fill_mapped_boxes.addBoxWithoutUpdate(fill_mapped_box);

                  dst_to_fill.insertLocalNeighbor(
                     fill_mapped_box,
                     base_box_itr);

                  constructed_fill_boxes.pushBack(*ei);
               }
            }
         }
      }

      d_max_fill_boxes = tbox::MathUtilities<int>::Max(
            d_max_fill_boxes,
            constructed_fill_boxes.size());
   }
   fill_mapped_boxes.finalize();
}

void
PatchLevelEnhancedFillPattern::computeDestinationFillBoxesOnSourceProc(
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
         "PatchLevelEnhancedFillPattern cannot compute destination:\n"
         << "fill boxes on the source processor.\n");
   }
}

bool
PatchLevelEnhancedFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return true;
}

bool
PatchLevelEnhancedFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return false;
}

bool
PatchLevelEnhancedFillPattern::fillingCoarseFineGhosts() const
{
   return true;
}

bool
PatchLevelEnhancedFillPattern::fillingEnhancedConnectivityOnly() const
{
   return true;
}

int
PatchLevelEnhancedFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
}

}
}
#endif
