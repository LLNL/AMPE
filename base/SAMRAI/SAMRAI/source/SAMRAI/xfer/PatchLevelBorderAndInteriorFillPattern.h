/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelBorderAndInteriorFillPattern
#define included_xfer_PatchLevelBorderAndInteriorFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief PatchLevelFillPattern implementation for filling at PatchLevel
 * boundaries and interiors
 *
 * For documentation on this interface see @ref xfer::PatchLevelFillPattern
 *
 * The fill boxes will consist of the interior of the destination
 * level as well as ghost regions lying outside of the level interior
 * at physical and coarse-fine boundaries.  Ghost regions that overlap
 * the interiors of other boxes on the destination level will not be
 * included.
 */

class PatchLevelBorderAndInteriorFillPattern:public PatchLevelFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   PatchLevelBorderAndInteriorFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~PatchLevelBorderAndInteriorFillPattern();

   /*!
    * @brief Compute the mapped boxes to be filled and related communication
    * data.
    *
    * The computed fill_mapped_boxes will cover the ghost regions around the
    * boxes of dst_mapped_box_level at coarse-fine and physical boundaries,
    * as well as the interior of the boxes of dst_mapped_box_level.  The
    * width of the ghost regions will be determined by fill_ghost_width.
    *
    * @param[out] fill_mapped_boxes    Output BoxLevel to be filled
    * @param[out] dst_to_fill          Output Connector between
    *                                  dst_mapped_box_level and
    *                                  and fill_mapped_boxes
    * @param[in] dst_mapped_box_level  destination level
    * @param[in] dst_to_dst            Connector of destination to itself
    * @param[in] dst_to_src            Connector of destination to source
    * @param[in] src_to_dst            Connector of source to destination
    * @param[in] fill_ghost_width      Ghost width being filled by refine
    *                                  schedule
    */
   void
   computeFillBoxesAndNeighborhoodSets(
      hier::BoxLevel& fill_mapped_boxes,
      hier::Connector& dst_to_fill,
      const hier::BoxLevel& dst_mapped_box_level,
      const hier::Connector& dst_to_dst,
      const hier::Connector& dst_to_src,
      const hier::Connector& src_to_dst,
      const hier::IntVector& fill_ghost_width);

   /*!
    * @brief Return true to indicate source patch owners cannot compute
    * fill boxes without using communication.
    */
   bool
   needsToCommunicateDestinationFillBoxes() const;

   /*!
    * @brief Virtual method to compute the destination fill boxes.
    *
    * Unimplemented for this class as needsToCommunicateDestinationFillBoxes()
    * returns true.
    *
    * @param[out] dst_fill_boxes_on_src_proc FillSet storing the destination
    *                                        neighbors of the source mapped
    *                                        boxes
    * @param[in] dst_mapped_box_level  destination level
    * @param[in] src_to_dst                Connector of source to destination
    * @param[in] fill_ghost_width          Ghost width being filled by refine
    *                                      schedule
    */
   void
   computeDestinationFillBoxesOnSourceProc(
      FillSet& dst_fill_boxes_on_src_proc,
      const hier::BoxLevel& dst_mapped_box_level,
      const hier::Connector& src_to_dst,
      const hier::IntVector& fill_ghost_width);

   /*!
    * @brief Tell RefineSchedule to communicate data directly from source
    * level to destination level.
    *
    * RefineSchedule should attempt to fill as much of the
    * destination level as possible from the source level at the same
    * level of resolution for this fill pattern.
    *
    * @return true
    */
   bool
   doesSourceLevelCommunicateToDestination() const;

   /*!
    * @brief Return the maximum number of fill boxes.
    *
    * @return maximum number of fill boxes.
    */
   int
   getMaxFillBoxes() const;

   /*!
    * @brief Returns true because this fill pattern fills coarse fine ghost
    * data.
    */
   bool
   fillingCoarseFineGhosts() const;

   /*!
    * @brief Returns false because this fill pattern is specialized for
    * enhanced connectivity only.
    */
   bool
   fillingEnhancedConnectivityOnly() const;

private:
   PatchLevelBorderAndInteriorFillPattern(
      const PatchLevelBorderAndInteriorFillPattern&);    // not implemented
   void
   operator = (
      const PatchLevelBorderAndInteriorFillPattern&);    // not implemented

   /*!
    * @brief Maximum number of fill boxes across all destination patches.
    */
   int d_max_fill_boxes;

};

}
}

#endif
