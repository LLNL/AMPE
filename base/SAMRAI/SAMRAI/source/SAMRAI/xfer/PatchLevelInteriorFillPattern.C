/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelInteriorFillPattern_C
#define included_xfer_PatchLevelInteriorFillPattern_C

#include "SAMRAI/xfer/PatchLevelInteriorFillPattern.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */

PatchLevelInteriorFillPattern::PatchLevelInteriorFillPattern():
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

PatchLevelInteriorFillPattern::~PatchLevelInteriorFillPattern()
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
PatchLevelInteriorFillPattern::computeFillBoxesAndNeighborhoodSets(
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
   NULL_USE(fill_ghost_width);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   const hier::BoxContainer& dst_mapped_boxes =
      dst_mapped_box_level.getBoxes();

   /*
    * Fill just the interior.  Disregard gcw.
    */
   for (hier::BoxContainer::const_iterator ni = dst_mapped_boxes.begin();
        ni != dst_mapped_boxes.end(); ++ni) {
      const hier::BoxId& gid = ni->getId();
      const hier::Box& dst_mapped_box =
         *dst_mapped_box_level.getBox(gid);
      fill_mapped_boxes.addBoxWithoutUpdate(dst_mapped_box);
      dst_to_fill.insertLocalNeighbor(dst_mapped_box, gid);
   }
   fill_mapped_boxes.finalize();
}

/*
 *************************************************************************
 *
 * computeDestinationFillBoxesOnSourceProc
 *
 *************************************************************************
 */

void
PatchLevelInteriorFillPattern::computeDestinationFillBoxesOnSourceProc(
   FillSet& dst_fill_boxes_on_src_proc,
   const hier::BoxLevel& dst_mapped_box_level,
   const hier::Connector& src_to_dst,
   const hier::IntVector& fill_ghost_width)
{
   NULL_USE(fill_ghost_width);
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_mapped_box_level, fill_ghost_width);

   const tbox::Dimension& dim(fill_ghost_width.getDim());
   const hier::IntVector& ratio(dst_mapped_box_level.getRefinementRatio());

   bool is_periodic = false;
   if (dst_mapped_box_level.getGridGeometry()->getPeriodicShift(ratio) != 
       hier::IntVector::getZero(dim)) {
      is_periodic = true;
   }

   /*
    * src_to_dst initialized only when there is a src mapped_box_level.
    * Without the src mapped_box_level, we do not need to compute
    * dst_fill_boxes_on_src_proc.
    *
    * For PatchLevelInteriorFillPattern, the src owner can compute fill
    * boxes for all its dst neighbors using local data.  This info is
    * stored in dst_fill_boxes_on_src_proc.
    */
   bool ordered = true;
   hier::BoxContainer all_dst_nabrs(ordered);
   if (is_periodic) {
      hier::BoxContainer tmp_nabrs(ordered);
      src_to_dst.getLocalNeighbors(tmp_nabrs);
      tmp_nabrs.unshiftPeriodicImageBoxes(
         all_dst_nabrs,
         dst_mapped_box_level.getRefinementRatio());
   } else {
      src_to_dst.getLocalNeighbors(all_dst_nabrs);
   }
   for (hier::BoxContainer::const_iterator na = all_dst_nabrs.begin();
        na != all_dst_nabrs.end(); ++na) {
      FillSet::Iterator dst_fill_boxes_iter =
         dst_fill_boxes_on_src_proc.insert(na->getId()).first;
      dst_fill_boxes_on_src_proc.insert(dst_fill_boxes_iter, *na);
      d_max_fill_boxes = tbox::MathUtilities<int>::Max(d_max_fill_boxes,
            static_cast<int>(dst_fill_boxes_on_src_proc.numNeighbors(
               dst_fill_boxes_iter)));
   }
}

bool
PatchLevelInteriorFillPattern::needsToCommunicateDestinationFillBoxes() const
{
   return false;
}

bool
PatchLevelInteriorFillPattern::doesSourceLevelCommunicateToDestination() const
{
   return true;
}

bool
PatchLevelInteriorFillPattern::fillingCoarseFineGhosts() const
{
   return false;
}

bool
PatchLevelInteriorFillPattern::fillingEnhancedConnectivityOnly() const
{
   return false;
}

int
PatchLevelInteriorFillPattern::getMaxFillBoxes() const
{
   return d_max_fill_boxes;
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
