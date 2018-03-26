/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeOverlap
#define included_pdat_NodeOverlap

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace pdat {

/**
 * Class NodeOverlap represents the intersection between two node
 * centered geometry boxes.  It is a subclass of hier::BoxOverlap and records
 * the portions of index space that needs to be copied between two objects
 * with node centered geometry.
 *
 * @see hier::BoxOverlap
 * @see pdat::NodeOverlap
 */

class NodeOverlap:public hier::BoxOverlap
{
public:
   /**
    * The constructor takes the list of boxes and the tranformation from
    * source to destination index spaces.  This information is used later
    * in the generation of communication schedules.
    */
   NodeOverlap(
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation);

   /**
    * The virtual destructor does nothing interesting except deallocate
    * box data.
    */
   virtual ~NodeOverlap();

   /**
    * Return whether there is an empty intersection between the two
    * node centered boxes.  This method over-rides the virtual function
    * in the hier::BoxOverlap base class.
    */
   virtual bool
   isOverlapEmpty() const;

   /**
    * Return the list of boxes (in node centered index space) that
    * constitute the intersection.  The boxes are given in the
    * destination coordinate space and must be shifted by
    * -(getSourceOffset()) to lie in the source index space.  This
    * method over-rides the virtual function in the
    * hier::BoxOverlap base class.
    */
   virtual const hier::BoxContainer&
   getDestinationBoxContainer() const;

   /**
    * Return the offset between the destination and source index spaces.
    * The destination index space is the source index spaced shifted
    * by this amount.
    */
   virtual const hier::IntVector&
   getSourceOffset() const;

   virtual const hier::Transformation&
   getTransformation() const;

private:
   bool d_is_overlap_empty;
   hier::Transformation d_transformation;
   hier::BoxContainer d_dst_boxes;

};

}
}
#endif
