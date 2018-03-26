/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelFullFillPattern
#define included_xfer_PatchLevelFullFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief PatchLevelFullFillPattern is a PatchLevelFillPattern that
 * fills the entire region the destination level, both interior and
 * ghost.
 *
 * For documentation on this interface see @ref xfer::PatchLevelFillPattern
 *
 * The fill boxes for this "All" PatchLevelFillPattern will consist of
 * the entire region of the destination level that can be filled, both
 * interior and ghost regions.
 *
 * If a RefineSchedule is created using an
 * RefineAlgorithm::createSchedule which takes no
 * PatchLevelFillPattern argument, this class will be used as the
 * default PatchLevelFillPattern.
 *
 * @see xfer::RefineAlgorithm
 * @see xfer::RefineSchedule
 */

class PatchLevelFullFillPattern:public PatchLevelFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   PatchLevelFullFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~PatchLevelFullFillPattern();

   /*!
    * @copydoc PatchLevelFillPattern::computeFillBoxesAndNeighborhoodSets()
    *
    * The computed fill_mapped_boxes for this fill pattern will be the
    * boxes of dst_mapped_box_level grown by the fill_ghost_width.
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
    * @copydoc PatchLevelFillPattern::needsToCommunicateDestinationFillBoxes()
    *
    * For this fill pattern, the source owner can compute fill boxes for
    * all of its destination neighbors using local data, so this method
    * returns false, allowing a communication step to be skipped.
    */
   bool
   needsToCommunicateDestinationFillBoxes() const;

   /*!
    * @copydoc PatchLevelFillPattern::computeDestinationFillBoxesOnSourceProc()
    */
   void
   computeDestinationFillBoxesOnSourceProc(
      FillSet& dst_fill_boxes_on_src_proc,
      const hier::BoxLevel& dst_mapped_box_level,
      const hier::Connector& src_to_dst,
      const hier::IntVector& fill_ghost_width);

   /*!
    * @copydoc PatchLevelFillPattern::doesSourceLevelCommunicateToDestination()
    *
    * RefineSchedule should attempt to fill the destination level from
    * the source level on the same resolution to the extent possible.
    */
   bool
   doesSourceLevelCommunicateToDestination() const;

   /*!
    * @copydoc PatchLevelFillPattern::getMaxFillBoxes()
    */
   int
   getMaxFillBoxes() const;

   /*!
    * @copydoc PatchLevelFillPattern::fillingCoarseFineGhosts()
    */
   bool
   fillingCoarseFineGhosts() const;

   /*!
    * @copydoc PatchLevelFillPattern::fillingEnhancedConnectivityOnly()
    */
   bool
   fillingEnhancedConnectivityOnly() const;

private:
   PatchLevelFullFillPattern(
      const PatchLevelFullFillPattern&);             // not implemented
   void
   operator = (
      const PatchLevelFullFillPattern&);             // not implemented

   /*!
    * @brief Maximum number of fill boxes across all destination patches.
    */
   int d_max_fill_boxes;
};

}
}

#endif
