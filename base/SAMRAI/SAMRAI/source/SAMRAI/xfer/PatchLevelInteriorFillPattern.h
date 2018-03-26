/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   fill pattern class for filling interiors only
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelInteriorFillPattern
#define included_xfer_PatchLevelInteriorFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief PatchLevelFillPattern implementation for patch interior filling.
 *
 * For documentation on this interface see @ref xfer::PatchLevelFillPattern
 *
 * Those fill boxes for this PatchLevelFillPattern will consist of the
 * patch interiors on the destination level only.
 *
 * @see xfer::RefineAlgorithm
 * @see xfer::RefineSchedule
 */

class PatchLevelInteriorFillPattern:public PatchLevelFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   PatchLevelInteriorFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~PatchLevelInteriorFillPattern();

   /*!
    * @brief Compute the mapped boxes to be filled and related communication
    * data.
    *
    * The computed fill_mapped_boxes will be the boxes of dst_mapped_box_level,
    * the "interior" of the destination level.
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
    * @brief Return false because source patch owner can compute destination
    * fill boxes.
    *
    * For this fill pattern, the source owner can compute fill boxes for
    * all of its destination neighbors using local data, so this method
    * returns false, allowing a communication step to be skipped.
    *
    * @return  false.
    */
   bool
   needsToCommunicateDestinationFillBoxes() const;

   /*!
    * @brief Compute the destination fill boxes.
    *
    * @param[in] dst_mapped_box_level  Destination level
    * @param[in] src_to_dst            Connector of source to destination
    * @param[in] fill_ghost_width      Ghost width filled by refine schedule
    * @param[out] dst_fill_boxes_on_src_proc FillSet storing the destination
    *                                        neighbors of the source mapped
    *                                        boxes
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
    * RefineSchedule should attempt to fill the destination level from
    * the source level on the same resolution to the extent possible.
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
    * @brief Returns false because this fill pattern fills only interior data.
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
   PatchLevelInteriorFillPattern(
      const PatchLevelInteriorFillPattern&);         // not implemented
   void
   operator = (
      const PatchLevelInteriorFillPattern&);        // not implemented

   /*!
    * @brief Maximum number of fill boxes across all destination patches.
    */
   int d_max_fill_boxes;

};

}
}

#endif
